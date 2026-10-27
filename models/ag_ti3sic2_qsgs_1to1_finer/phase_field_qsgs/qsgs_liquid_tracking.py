"""Temporal identity and geometry tracking for liquid QSGS voxel clusters."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qsgs_liquid_connectivity import liquid_cluster_labels


def _minimal_periodic_width(indices: np.ndarray, count: int, spacing_A: float) -> float:
    """Smallest periodic interval covering voxel centers, reported in voxel widths."""
    unique = np.unique(indices)
    if unique.size == 0:
        return 0.0
    if unique.size == count:
        return float(count * spacing_A)
    gaps = np.diff(np.concatenate((unique, unique[:1] + count)))
    covered_voxels = count - int(gaps.max()) + 1
    return float(covered_voxels * spacing_A)


def _component_geometry(
    labels: np.ndarray,
    component_count: int,
    spacing_A: np.ndarray,
    domain_size_A: np.ndarray,
) -> list[dict[str, int | float | list[float]]]:
    records = []
    shape = np.asarray(labels.shape)
    for component_id in range(1, component_count + 1):
        coords = np.where(labels == component_id)
        voxel_count = int(coords[0].size)
        if not voxel_count:
            continue
        centers = [(coords[axis].astype(float) + 0.5) * spacing_A[axis] for axis in range(3)]
        centroids = []
        for axis in (0, 1):
            angles = 2.0 * np.pi * centers[axis] / domain_size_A[axis]
            angle = np.arctan2(np.sin(angles).mean(), np.cos(angles).mean()) % (2.0 * np.pi)
            centroids.append(float(angle * domain_size_A[axis] / (2.0 * np.pi)))
        centroids.append(float(centers[2].mean()))
        z_min_A = float(coords[2].min() * spacing_A[2])
        z_max_A = float((coords[2].max() + 1) * spacing_A[2])
        records.append(
            {
                "component_id": component_id,
                "voxel_count": voxel_count,
                "volume_A3": float(voxel_count * np.prod(spacing_A)),
                "centroid_A": centroids,
                "depth_from_top_A": float(domain_size_A[2] - centroids[2]),
                "z_min_A": z_min_A,
                "z_max_A": z_max_A,
                "lateral_size_x_A": _minimal_periodic_width(coords[0], int(shape[0]), float(spacing_A[0])),
                "lateral_size_y_A": _minimal_periodic_width(coords[1], int(shape[1]), float(spacing_A[1])),
            }
        )
    return records


@dataclass
class LiquidClusterTracker:
    """Assign persistent IDs by voxel overlap and retain merge/split lineage.

    When one component splits, the child with the largest overlap inherits its
    parent's ID; other children receive new IDs with that parent recorded. When
    components merge, the parent with the largest overlap (with ties resolved
    by component size and then ID) supplies the continuing ID. These are
    explicit bookkeeping rules, not a physical coalescence model.
    """

    spacing_A: np.ndarray | tuple[float, float, float]
    domain_size_A: np.ndarray | tuple[float, float, float]
    liquid_threshold: float = 0.5
    previous_track_labels: np.ndarray | None = None
    previous_time_ps: float | None = None
    next_track_id: int = 1
    track_start_ps: dict[int, float] = field(default_factory=dict)
    closed_tracks: list[dict[str, int | float | str]] = field(default_factory=list)
    history: list[dict[str, int | float | str | bool | list[int] | list[float]]] = field(default_factory=list)
    first_liquid_time_ps: float | None = None
    first_connected_liquid_time_ps: float | None = None
    first_multivoxel_pool_time_ps: float | None = None

    def __post_init__(self) -> None:
        self.spacing_A = np.asarray(self.spacing_A, dtype=float)
        self.domain_size_A = np.asarray(self.domain_size_A, dtype=float)
        if self.spacing_A.shape != (3,) or self.domain_size_A.shape != (3,):
            raise ValueError("spacing_A and domain_size_A must each have three entries")
        if np.any(self.spacing_A <= 0.0) or np.any(self.domain_size_A <= 0.0):
            raise ValueError("spacing_A and domain_size_A must be positive")
        if not 0.0 < self.liquid_threshold <= 1.0:
            raise ValueError("liquid_threshold must be in (0, 1]")
        if self.previous_track_labels is not None:
            self.previous_track_labels = np.asarray(self.previous_track_labels, dtype=np.int32)

    def update(
        self,
        time_ps: float,
        liquid_fraction: np.ndarray,
        ag_mask: np.ndarray,
    ) -> tuple[dict[str, int | float | list[float] | None], list[dict[str, object]]]:
        liquid_fraction = np.asarray(liquid_fraction, dtype=float)
        ag_mask = np.asarray(ag_mask, dtype=bool)
        if liquid_fraction.shape != ag_mask.shape:
            raise ValueError("liquid_fraction and ag_mask must have matching shapes")
        if self.previous_track_labels is not None and self.previous_track_labels.shape != ag_mask.shape:
            raise ValueError("Previous tracker labels do not match the current geometry")
        if self.first_liquid_time_ps is None and np.any(ag_mask & (liquid_fraction > 1.0e-12)):
            self.first_liquid_time_ps = float(time_ps)

        liquid_mask = ag_mask & (liquid_fraction >= self.liquid_threshold)
        component_labels, component_count = liquid_cluster_labels(liquid_mask)
        geometry = _component_geometry(
            component_labels, component_count, self.spacing_A, self.domain_size_A
        )
        if component_count and self.first_connected_liquid_time_ps is None:
            self.first_connected_liquid_time_ps = float(time_ps)
        if any(item["voxel_count"] >= 2 for item in geometry) and self.first_multivoxel_pool_time_ps is None:
            self.first_multivoxel_pool_time_ps = float(time_ps)

        overlaps: dict[tuple[int, int], int] = {}
        if self.previous_track_labels is not None and component_count:
            both = (self.previous_track_labels > 0) & (component_labels > 0)
            if np.any(both):
                previous_ids, previous_inverse = np.unique(
                    self.previous_track_labels[both], return_inverse=True
                )
                current_ids = component_labels[both].astype(np.int64)
                encoded = previous_inverse.astype(np.int64) * (component_count + 1) + current_ids
                values, counts = np.unique(encoded, return_counts=True)
                for value, count in zip(values, counts, strict=True):
                    previous_index, current_id = divmod(int(value), component_count + 1)
                    overlaps[(int(previous_ids[previous_index]), current_id)] = int(count)

        parents_by_component: dict[int, list[tuple[int, int]]] = {i: [] for i in range(1, component_count + 1)}
        children_by_parent: dict[int, list[tuple[int, int]]] = {}
        for (parent_id, component_id), count in overlaps.items():
            parents_by_component[component_id].append((parent_id, count))
            children_by_parent.setdefault(parent_id, []).append((component_id, count))
        preferred_child: dict[int, int] = {}
        size_by_component = {int(item["component_id"]): int(item["voxel_count"]) for item in geometry}
        for parent_id, children in children_by_parent.items():
            preferred_child[parent_id] = min(
                children,
                key=lambda pair: (-pair[1], -size_by_component[pair[0]], pair[0]),
            )[0]

        current_track_labels = np.zeros(component_labels.shape, dtype=np.int32)
        used_track_ids: set[int] = set()
        current_records: list[dict[str, object]] = []
        previous_track_ids = (
            set(int(value) for value in np.unique(self.previous_track_labels) if value > 0)
            if self.previous_track_labels is not None
            else set()
        )
        for component in geometry:
            component_id = int(component["component_id"])
            parents = sorted(
                parents_by_component[component_id],
                key=lambda pair: (-pair[1], -size_by_component[component_id], pair[0]),
            )
            inheritable = [
                (parent_id, count)
                for parent_id, count in parents
                if preferred_child.get(parent_id) == component_id and parent_id not in used_track_ids
            ]
            candidates = inheritable or [
                (parent_id, count) for parent_id, count in parents if parent_id not in used_track_ids
            ]
            if candidates:
                track_id = min(candidates, key=lambda pair: (-pair[1], pair[0]))[0]
                used_track_ids.add(track_id)
            else:
                track_id = self.next_track_id
                self.next_track_id += 1
                self.track_start_ps[track_id] = float(time_ps)
            self.track_start_ps.setdefault(track_id, float(time_ps))

            split = any(len(children_by_parent.get(parent_id, [])) > 1 for parent_id, _ in parents)
            merged = len(parents) > 1
            if not parents:
                event = "birth"
            elif split and merged:
                event = "merge+split"
            elif merged:
                event = "merge"
            elif split and preferred_child.get(parents[0][0]) != component_id:
                event = "split_child"
            elif split:
                event = "split_parent_continues"
            else:
                event = "continuation"

            component_mask = component_labels == component_id
            current_track_labels[component_mask] = track_id
            record = {
                "time_ps": float(time_ps),
                "track_id": int(track_id),
                "parent_track_ids": [int(parent_id) for parent_id, _ in parents],
                "event": event,
                "is_main_cluster": False,
                "continuous_duration_ps": float(time_ps - self.track_start_ps[track_id]),
                **component,
            }
            current_records.append(record)

        if current_records:
            main_index = max(
                range(len(current_records)),
                key=lambda idx: (int(current_records[idx]["voxel_count"]), -int(current_records[idx]["track_id"])),
            )
            current_records[main_index]["is_main_cluster"] = True

        current_track_ids = set(used_track_ids)
        for track_id in sorted(previous_track_ids - current_track_ids):
            start = self.track_start_ps.get(track_id, float(time_ps))
            self.closed_tracks.append(
                {
                    "track_id": track_id,
                    "start_time_ps": start,
                    "end_time_ps": float(self.previous_time_ps if self.previous_time_ps is not None else time_ps),
                    "sampled_continuous_duration_ps": float(
                        (self.previous_time_ps if self.previous_time_ps is not None else time_ps) - start
                    ),
                    "end_event": "no_overlap_in_current_frame",
                }
            )
        self.history.extend(current_records)
        self.previous_track_labels = current_track_labels
        self.previous_time_ps = float(time_ps)

        largest = max(geometry, key=lambda item: int(item["voxel_count"]), default=None)
        summary = {
            "connected_liquid_cluster_count": component_count,
            "liquid_voxel_count": int(liquid_mask.sum()),
            "largest_liquid_cluster_voxels": int(largest["voxel_count"]) if largest else 0,
            "largest_liquid_cluster_volume_A3": float(largest["volume_A3"]) if largest else 0.0,
            "largest_liquid_cluster_centroid_A": largest["centroid_A"] if largest else None,
            "largest_liquid_track_id": next(
                (int(item["track_id"]) for item in current_records if item["is_main_cluster"]), None
            ),
            "largest_liquid_track_duration_ps": next(
                (float(item["continuous_duration_ps"]) for item in current_records if item["is_main_cluster"]),
                0.0,
            ),
        }
        return summary, current_records

    def to_metadata(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "liquid_threshold": float(self.liquid_threshold),
            "previous_time_ps": self.previous_time_ps,
            "next_track_id": int(self.next_track_id),
            "track_start_ps": {str(key): value for key, value in self.track_start_ps.items()},
            "closed_tracks": self.closed_tracks,
            "history": self.history,
            "first_liquid_time_ps": self.first_liquid_time_ps,
            "first_connected_liquid_time_ps": self.first_connected_liquid_time_ps,
            "first_multivoxel_pool_time_ps": self.first_multivoxel_pool_time_ps,
        }

    def load_metadata(self, metadata: dict[str, object], previous_track_labels: np.ndarray | None) -> None:
        if int(metadata.get("schema_version", -1)) != 1:
            raise ValueError("Unsupported liquid tracker checkpoint schema")
        if not np.isclose(float(metadata["liquid_threshold"]), self.liquid_threshold, rtol=0.0, atol=0.0):
            raise ValueError("Liquid threshold differs from the checkpoint")
        self.previous_time_ps = (
            None if metadata.get("previous_time_ps") is None else float(metadata["previous_time_ps"])
        )
        self.next_track_id = int(metadata["next_track_id"])
        self.track_start_ps = {int(key): float(value) for key, value in metadata["track_start_ps"].items()}
        self.closed_tracks = list(metadata["closed_tracks"])
        self.history = list(metadata["history"])
        self.first_liquid_time_ps = metadata.get("first_liquid_time_ps")
        self.first_connected_liquid_time_ps = metadata.get("first_connected_liquid_time_ps")
        self.first_multivoxel_pool_time_ps = metadata.get("first_multivoxel_pool_time_ps")
        self.previous_track_labels = (
            None if previous_track_labels is None else np.asarray(previous_track_labels, dtype=np.int32)
        )
