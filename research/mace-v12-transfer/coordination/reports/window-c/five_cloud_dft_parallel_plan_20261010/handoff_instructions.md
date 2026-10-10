# D–H 云端可复制交接指令（只读预检版）

**这份交接仅授权环境/数据来源核查和任务划分，不是结构生成或 DFT 放行。** 不运行 MD/AIMD、GPAW 单点、MACE 训练、LAMMPS 或 TTM。不得改 A/B 队列、输入、进程、checkpoint 或成果。

每个云端只收到自己的段落。任务 ID 是 C 的规划编号；先确认 owner、现有任务状态和几何去重，不能按 task ID 推断任务已经派发。D/E/F 如有运行中/排队任务，原样保留；仅评估尚未启动构型。G/H 在六点 pilot 阶段预留。

## 通用只读资源预检

在各自独立 checkout/计算机执行，避免 `set -x`。不要把 shell 环境、凭证或原始私有路径贴回 GitHub。仅回报脱敏数值、版本、组件/许可状态和路径 basename（如确有需要）：

```bash
python3 - <<'PY'
import os, sys, shutil, importlib.metadata as md
from pathlib import Path
print('python', sys.version.split()[0])
print('cpu_logical', os.cpu_count())
try:
    rows = Path('/proc/cpuinfo').read_text().split('\n\n')
    pairs = set()
    for row in rows:
        d = dict(line.split(':', 1) for line in row.splitlines() if ':' in line)
        if 'physical id' in d and 'core id' in d:
            pairs.add((d['physical id'].strip(), d['core id'].strip()))
    print('cpu_physical', len(pairs) or 'unknown')
except OSError:
    print('cpu_physical', 'unknown')
try:
    mem = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        k, v = line.split(':', 1)
        if k in ('MemTotal', 'MemAvailable'):
            mem[k] = int(v.split()[0])
    print('memory_total_GiB', round(mem.get('MemTotal', 0)/1048576, 1))
    print('memory_available_GiB', round(mem.get('MemAvailable', 0)/1048576, 1))
except OSError:
    print('memory', 'unknown')
st = os.statvfs('.')
print('workdir_free_GB', round(st.f_bavail*st.f_frsize/1e9, 1))
for package in ('gpaw', 'gpaw-data', 'ase', 'numpy', 'scipy'):
    try: print(package, md.version(package))
    except md.PackageNotFoundError: print(package, 'not-installed')
print('mpiexec_available', bool(shutil.which('mpiexec') or shutil.which('mpirun')))
try:
    import gpaw_data
    p = Path(gpaw_data.datapath())
    print('gpaw_data_directory_present', p.is_dir())
except Exception:
    print('gpaw_data_import_or_directory', 'unavailable')
print('gpaw_setup_path_configured', bool(os.environ.get('GPAW_SETUP_PATH')))
PY
if command -v mpiexec >/dev/null 2>&1; then mpiexec --version | sed -n '1,2p'; fi
```

另外由当地 GPAW 环境 owner 核对并脱敏填写：GPAW build 是否启用 MPI、MPI implementation/runtime ABI、BLAS/LAPACK、ScaLAPACK、LibXC 能力；PAW 数据版本/来源/许可/可用性与 SHA256。不要用上述资源脚本启动 GPAW，不进行 k 点试算或性能基准。

每台给 C 一条状态行：`host | running IDs (count only if IDs sensitive) | queued | unstarted | cores | memory | free GB | GPAW/ASE/MPI/PAW versions and license status | comparable-job elapsed/rank count | blocker`。无需提供绝对私有路径。任务/数据角色在 C 核对后才能固定。

## 云端 D

拟议 full-30 分配：4 liquid Ag（含 2 个首批 1250 K pilot 候选）+ 1 Ti₃SiC₂ + 1 solid Ag。首批 pilot 候选：

- `C-D-AGL-1250K-rho097-seed01`
- `C-D-AGL-1250K-rho103-seed02`

请先只读盘点 D 本机状态和资源；如存在运行中/排队任务，保护原样并标识，不能重复分派。**当前不启动以上 pilot，也不生成结构。** 液体密度、采样器、许可和候选坐标待审。

## 云端 E

拟议 full-30 分配：4 liquid Ag（含 2 个首批 1250 K pilot 候选）+ 1 Ti₃SiC₂ + 1 solid Ag。首批 pilot 候选：

- `C-E-AGL-1250K-rho097-seed02`
- `C-E-AGL-1250K-rho100-seed01`

只提交脱敏环境和已有任务状态。活动任务/checkpoint 不动；未启动任务只有在后续明确放行后才接单。当前不生成结构或跑 DFT。

## 云端 F

拟议 full-30 分配：4 liquid Ag（含 2 个首批 1250 K pilot 候选）+ 1 Ti₃SiC₂ + 1 solid Ag。首批 pilot 候选：

- `C-F-AGL-1250K-rho100-seed02`
- `C-F-AGL-1250K-rho103-seed01`

只读检查主机资源、GPAW/MPI/PAW 来源与已运行任务；不安装/替换 GPAW 栈，不覆盖任何状态。当前不生成结构或跑 DFT。

## 云端 G

拟议 full-30 分配：3 liquid Ag + 2 Ti₃SiC₂ + 1 solid Ag。**六点 pilot 不派给 G**；G 在 pilot 阶段预留给固相和之后经批准的正式批次。仅提交资源、来源/许可及任务状态预检；不得提前开始 30 点矩阵中的任何计算。

## 云端 H

拟议 full-30 分配：3 liquid Ag + 1 Ti₃SiC₂ + 2 solid Ag。**六点 pilot 不派给 H**；H 在 pilot 阶段预留给固相和之后经批准的正式批次。仅提交资源、来源/许可及任务状态预检；不得提前开始 30 点矩阵中的任何计算。

## 后续 Git 和文件隔离约定

获得计算放行和几何输入后，各 owner 使用自己的 Git branch 和独立工作根目录。每个 task 的目录必须在本机、按 task ID 独立，不共用 checkpoint。源输入只读冻结并先核 SHA；多主机副本的输入 SHA 必须完全相同。不要 force push、不要推送 `main`、不要上传活跃 `.gpw` checkpoint。若主机 task inventory 与本表冲突，先向 C 报告，保留已运行结果，只重排 unstarted case。
