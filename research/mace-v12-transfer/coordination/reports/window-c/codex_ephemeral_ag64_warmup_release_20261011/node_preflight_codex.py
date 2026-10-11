#!/usr/bin/env python3
"""Read-only owner-side resource/identity receipt; never creates roots or starts sampling."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

from codex_release_gate import MODE_EPHEMERAL, load_signed_config, machine_fingerprint, sha256_file


def nearest_existing(path: Path) -> Path:
    current = path
    while not current.exists():
        current = current.parent
    return current


def cgroup_value(paths):
    for path in paths:
        if path.is_file():
            raw = path.read_text().strip()
            if raw != "max":
                try:
                    return int(raw)
                except ValueError:
                    continue
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--owner", choices=["D", "E", "F"], required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--potential-path", type=Path)
    args = parser.parse_args()
    try:
        config, config_sha = load_signed_config(args.config, args.config_sha256)
        hosts = config.get("hosts", {})
        if args.owner not in hosts:
            raise ValueError("OWNER_NOT_REGISTERED")
        host = hosts[args.owner]
        if host.get("owner") != args.owner or args.task_id not in host.get("allowed_task_ids", []):
            raise ValueError("OWNER_OR_TASK_ID_MISMATCH")
        reasons = []
        checks = {}
        checks["storage_mode"] = config.get("storage_mode") == MODE_EPHEMERAL
        if not checks["storage_mode"]:
            reasons.append("CODEX_EPHEMERAL_MODE_REQUIRED")
        checks["persistence_verified_is_false"] = host.get("persistence_verified") is False
        if not checks["persistence_verified_is_false"]:
            reasons.append("PERSISTENCE_FLAG_MUST_REMAIN_FALSE")
        checks["C_release_signed"] = config.get("config_status") == "C_SIGNED_FOR_WARMUP" and host.get("C_host_release") is True
        if not checks["C_release_signed"]:
            reasons.append("C_CONFIG_OR_OWNER_RELEASE_HOLD")

        root = Path(host["sampler_work_root"])
        staging = Path(host["archive_staging_root"])
        script_path = Path(__file__).resolve()
        checkout = next((p for p in script_path.parents if (p / ".git").exists()), None)
        checks["workspace_paths_absolute"] = root.is_absolute() and staging.is_absolute()
        symlink_free = True
        for candidate in (root, staging):
            for component in (candidate, *candidate.parents):
                if component.is_symlink():
                    symlink_free = False
        checks["root_components_have_no_symlinks"] = symlink_free
        if not symlink_free:
            reasons.append("ROOT_PATH_CONTAINS_SYMLINK")
        checks["roots_outside_checkout"] = bool(checkout and
             not (root == checkout or checkout in root.parents) and
             not (staging == checkout or checkout in staging.parents))
        checks["sampler_root_absent"] = not root.exists()
        if not checks["sampler_root_absent"]:
            reasons.append("FROZEN_TASK_ROOT_ALREADY_EXISTS_REVIEW_REQUIRED")
        root_pairs = [(name, field, Path(item[field])) for name, item in hosts.items()
                      for field in ("sampler_work_root", "archive_staging_root")]
        distinct = True
        for i, (owner_a, field_a, path_a) in enumerate(root_pairs):
            for owner_b, field_b, path_b in root_pairs[i + 1:]:
                if path_a == path_b or path_a in path_b.parents or path_b in path_a.parents:
                    distinct = False
        checks["configured_roots_nonoverlap"] = distinct
        if not distinct:
            reasons.append("CONFIGURED_ROOTS_OVERLAP")

        anchor = nearest_existing(root.parent)
        archive_anchor = nearest_existing(staging.parent)
        usage = shutil.disk_usage(anchor)
        vfs = os.statvfs(anchor)
        min_free = int(config["sampling"]["min_free_bytes"])
        checkpoint_reserve = int(config["sampling"]["checkpoint_reserve_bytes"])
        free_inodes = vfs.f_favail
        rss_bytes = None
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                rss_bytes = int(line.split()[1]) * 1024
                break
        cpu_quota = None
        cpu_max = Path("/sys/fs/cgroup/cpu.max")
        if cpu_max.is_file():
            quota, period = cpu_max.read_text().split()
            if quota != "max":
                cpu_quota = int(quota) / int(period)
        else:
            quota_path = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
            period_path = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
            if quota_path.is_file() and period_path.is_file() and int(quota_path.read_text()) > 0:
                cpu_quota = int(quota_path.read_text()) / int(period_path.read_text())
        memory_limit = cgroup_value((Path("/sys/fs/cgroup/memory.max"),
                                     Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")))
        checks["workspace_min_free_space"] = usage.free >= min_free + checkpoint_reserve
        checks["workspace_min_free_inodes"] = free_inodes >= int(config["codex_archive_policy"].get("min_free_inodes", 10000))
        checks["rss_below_task_limit"] = rss_bytes is not None and rss_bytes < int(config["sampling"]["max_rss_bytes"])
        checks["cpu_quota_visible"] = cpu_quota is not None
        checks["memory_limit_sufficient"] = memory_limit is not None and memory_limit >= int(config["sampling"]["max_rss_bytes"])
        checks["archive_staging_same_filesystem"] = anchor.stat().st_dev == archive_anchor.stat().st_dev
        for name, ok in checks.items():
            if not ok and name not in {"C_release_signed", "storage_mode", "persistence_verified_is_false"}:
                reasons.append(name.upper() + "_FAIL")
        remote = host.get("remote_archive", {})
        checks["synthetic_remote_roundtrip"] = (remote.get("roundtrip_status") == "PASS"
                                                 and remote.get("synthetic_fixture_verified") is True)
        if not checks["synthetic_remote_roundtrip"]:
            reasons.append("REMOTE_SYNTHETIC_ROUNDTRIP_NOT_VERIFIED")
        fingerprint = machine_fingerprint(Path("/workspace"))
        checks["machine_fingerprint_matches_record"] = fingerprint == host.get("machine_fingerprint_sha256")
        if not checks["machine_fingerprint_matches_record"]:
            reasons.append("NODE_FINGERPRINT_CHANGED_REQUIRES_C_RESIGN")
        potential_result = {"provided": args.potential_path is not None, "path": None, "basename": None,
                            "sha256": None, "matches_exact_frozen_sha": False}
        if args.potential_path is not None:
            potential_path = args.potential_path.resolve(strict=False)
            potential_result.update({"path": str(potential_path), "basename": potential_path.name})
            if potential_path.is_file():
                potential_sha = sha256_file(potential_path)
                potential_result["sha256"] = potential_sha
                potential_result["matches_exact_frozen_sha"] = (
                    potential_path.name == "Ag.eam.alloy" and
                    potential_sha == config["potential"]["sha256"])
        checks["exact_node_local_potential"] = potential_result["matches_exact_frozen_sha"]
        if not checks["exact_node_local_potential"]:
            reasons.append("NODE_LOCAL_AG_EAM_FILE_SHA_NOT_VERIFIED_NO_SUBSTITUTE")
        result = {
            "owner": args.owner, "task_id": args.task_id, "config_sha256": config_sha,
            "machine_fingerprint_sha256": fingerprint,
            "configured_machine_fingerprint_sha256": host.get("machine_fingerprint_sha256"),
            "sampler_work_root": str(root), "archive_staging_root": str(staging),
            "filesystem_anchor": str(anchor), "workspace_available_bytes": usage.free,
            "workspace_free_inodes": free_inodes, "current_rss_bytes": rss_bytes,
            "cpu_quota_cores": cpu_quota, "cgroup_memory_limit_bytes": memory_limit,
            "potential": potential_result, "checks": checks,
            "physical_host_exclusivity": "UNVERIFIED_NOT_CLAIMED",
            "persistence_verified": False, "root_created": False,
            "science_execution_started": False,
            "reasons": reasons,
            "PREFLIGHT_STATUS": "PASS" if not reasons else "HOLD",
        }
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if not reasons else 2
    except Exception as exc:
        print(json.dumps({"owner": args.owner, "task_id": args.task_id,
                          "status": "PREFLIGHT_ERROR", "reason": str(exc),
                          "science_execution_started": False}, indent=2), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
