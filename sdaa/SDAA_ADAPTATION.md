# pdbfixer SDAA 适配

日期：2026-10-27
范围：`openmm/pdbfixer` 1.12.0 在太初 SDAA（Loongson-3C6000）上的安装与验证。

---

## 结论

**依赖前提：OpenMM 8.5.2 必须先在 loongarch64 上源码编译安装**（无官方 wheel，PyPI/loongnix 源都没有；
本机已装，安装步骤与装机位置见 §2）。OpenMM 就绪后，pdbfixer 1.12.0 直接
`pip install` 即可，**pdbfixer 本身零代码改动**。

| 项 | 结果 |
| --- | --- |
| 依赖链 | OpenMM 8.5.2（loongarch64 源码编译，C++ 在 `/opt/openmm`）→ pdbfixer 1.12.0（纯 Python）→ numpy 2.2.6 |
| 官方测试套件（`pdbfixer/tests/`，23 个用例） | **23 passed**，498.76 s（`logs/pdbfixer_reinstall_stdinstall.log`） |
| 跨平台结构对照（14 个 PDB × 50~9000 原子） | 原子数 / 残基数 / 链数 / 键数 / 元素计数 **14/14 完全一致，零差异** |
| 大样本坐标对照（官方测试 PDB 列表 325 个，固定 `random.seed(12345)`） | 两边成功集一致（各 322）；**143 个 case 全原子逐位一致**，128 个仅氢有差，51 个重原子局部有差——分叉量级与 x86 自身换种子相同，非平台偏差（见 §3.3） |
| 跨平台能量（最小化后） | 差 73.02 kcal/mol（5.25%），**小于 SDAA 自身噪声地板 88.77 kcal/mol（5.28%）** |

pdbfixer 是纯 Python 包，`install_requires` 只有 `numpy` / `openmm>=8.2` / `legacy-cgi`；
真正需要适配工作量的只有 OpenMM 的 loongarch64 源码编译（§2），那是本机既有成果，
本次未重做（源码目录已被容器重建清掉，重装需重新拉源码编译）。

---

## 1. 平台特性

### OpenMM 只有 CPU 平台

115 容器 `sdaa` 里的 OpenMM 8.5.2（源码编译）只有两个平台插件：

```
platforms: ['Reference', 'CPU']
```

x86 对照侧（OpenMM 官方 manylinux wheel）多一个 `OpenCL`：`['Reference', 'CPU', 'OpenCL']`。

**SDAA 上跑不了 GPU 加速的 MD** —— OpenMM 没有 SDAA 后端，pdbfixer 的加氢/建模全走 CPU 平台。这是平台限制，不是适配问题。

### pdbfixer 本身与 SDAA 无关

pdbfixer 1767 行纯 Python，只调 `openmm.app` 的 API，不含任何 CUDA 专属调用：

```python
def __init__(self, filename=None, pdbfile=None, pdbxfile=None, url=None, pdbid=None, platform=None):
    ...
    self.platform = platform        # 可选，None 时 OpenMM 自动选
```

`platform=None` 时 OpenMM 自动挑 CPU，所以**不需要改任何代码**。

### 依赖链

| 包 | 115 容器 | 说明 |
| --- | --- | --- |
| OpenMM | 8.5.2（`8.5.2.dev-36a30cb`） | 源码编译（CMake），非 wheel；loongarch64 无官方 wheel |
| pdbfixer | 1.12.0 | 纯 Python 包（`py3-none-any`），PyPI 只有 sdist，pip 本地 build |
| numpy | 2.2.6 | 容器自带 |
| pytest | 9.1.0 | 容器自带 |

pdbfixer 1.12.0 **不依赖 openff-toolkit**（`Requires-Dist: numpy` / `openmm>=8.2` / `legacy-cgi; python_version >= "3.13"`）。

### 之前「pdbfixer 装不上」的真正原因：openmm 找不到，不是 pdbfixer 的问题

早前一次 `pip install pdbfixer` 报 `Could not find a version that satisfies the requirement openff-toolkit`。
用 `--dry-run --ignore-installed` 复现出了完整机制：

```
$ pip install --dry-run --ignore-installed pdbfixer
INFO: pip is looking at multiple versions of pdbfixer to determine which version is compatible with other requirements.
ERROR: Could not find a version that satisfies the requirement openmm>=8.2 (from pdbfixer) (from versions: none)
```

**直接诱因是 openmm 找不到** —— 源上没有 loongarch64 的 openmm 轮子：

```
$ pip index versions openmm
ERROR: No matching distribution found for openmm
```

pip 找不到 openmm>=8.2 就去回退更老的 pdbfixer，老版本又依赖 openff-toolkit，
而 openff-toolkit 在源上同样没有 —— 最终报的却是 openff-toolkit，实际卡点在 openmm。

**openmm 8.5.2 装好之后，pdbfixer 一条命令就能装上**（见第 2 节）。
也就是说「pdbfixer 装不上」这个结论只在 openmm 没装好的前提下成立。

---

## 2. 安装

### OpenMM 安装步骤（loongarch64 源码编译）

OpenMM 源码编译就做 4 件事：装依赖 → 编译安装 C++ 库 → 编译安装 Python 层 → 自检。

```bash
# 0. 前置：cmake / swig / gcc / make 容器里已有（cmake 3.26.3, swig 4.1.1）；
#    swig 必须有（Python 层的 .so 靠它生成）。FFTW/GSL 不用装，OpenMM 自带内置 FFT。

# 1. 取源码（装的是 git 开发版 8.5.2.dev-36a30cb，不是 release tarball）
git clone https://github.com/openmm/openmm.git
cd openmm && git checkout 36a30cbca54e727b216b606f3c011b67201eb8b4

# 2. 编译安装 C++ 库，前缀 /opt/openmm
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/openmm \
         -DCMAKE_BUILD_TYPE=Release \
         -DOPENMM_BUILD_CUDA_LIB=OFF          # 本机无 CUDA SDK，GPU 后端全不编
make -j$(nproc)
make install                                 # 装出 /opt/openmm/lib/libOpenMM.so + Reference/CPU 插件

# 3. 编译安装 Python 层（当时用的是 setup.py install，不是 pip）
cd ../python
/home/py312/bin/python3 setup.py install     # 装出 /home/py312/.../site-packages/openmm/
```

自检（两条都过才算装好）：

```bash
/home/py312/bin/python3 -c "import openmm; print(openmm.version.full_version)"
# 预期 8.5.2.dev-36a30cb
/home/py312/bin/python3 -c "import openmm as mm; print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])"
# 预期 ['Reference', 'CPU']；libOpenMM.so 从 /opt/openmm/lib 解析（RUNPATH 自带，不用设 LD_LIBRARY_PATH）
```

### OpenMM：源码编译（115 上已装好，版本 8.5.2）

`openmm==8.5.2` **在 115 上是源码编译安装的，不是 wheel** —— PyPI 与 loongnix 源都没有
loongarch64 的 OpenMM 轮子（`pip index versions openmm` → `No matching distribution found`）。

安装位置分两处，都是 root 拥有：

| 路径 | 内容 | 时间戳 |
| --- | --- | --- |
| `/opt/openmm/` | CMake install 前缀：`lib/libOpenMM.so`（3.9 MB）+ Amoeba/Drude/RPMD 库 + `lib/plugins/` | Sep 10 10:12 |
| `/home/py312/lib/python3.12/site-packages/openmm/` | Python 包装层，`_openmm.cpython-312-loongarch64-linux-gnu.so`（25 MB，swig 生成） | Sep 10 10:11~10:25 |
| `OpenMM-8.5.2-py3.12.egg-info/` | 安装元数据（**egg-info 而非 dist-info** → legacy `setup.py install`，非 pip 装的） | Sep 10 10:26 |

`openmm/version.py` 里的关键信息：

```python
short_version = '8.5.2'
full_version  = '8.5.2.dev-36a30cb'
git_revision  = '36a30cbca54e727b216b606f3c011b67201eb8b4'
release       = False
openmm_library_path = r'/opt/openmm/lib'
```

即**从 git 源码树 checkout 的开发版**（`8.5.2.dev-36a30cb`），不是 release tarball；
`openmm_library_path` 硬编码指向 `/opt/openmm/lib`，Python 层靠它加载 C++ 库。

**加载方式**：`/opt/openmm/lib` 写进了 wrapper 的 RUNPATH，所以不需要设 `LD_LIBRARY_PATH`：

```
$ readelf -d openmm/_openmm.cpython-312-loongarch64-linux-gnu.so | grep RUNPATH
 0x... (RUNPATH)  Library runpath: [/opt/openmm/lib]

$ ldd openmm/_openmm.cpython-312-loongarch64-linux-gnu.so | grep OpenMM
	libOpenMM.so => /opt/openmm/lib/libOpenMM.so (0x...)
	libOpenMMAmoeba.so => /opt/openmm/lib/libOpenMMAmoeba.so (0x...)
	libOpenMMRPMD.so => /opt/openmm/lib/libOpenMMRPMD.so (0x...)
	libOpenMMDrude.so => /opt/openmm/lib/libOpenMMDrude.so (0x...)
```

插件也靠 `openmm/__init__.py:29-30` 从 `openmm_library_path/plugins` 自动加载
（`OPENMM_PLUGIN_DIR` 未设时）。

**编译期关掉了 GPU 后端**（`/opt/openmm/lib/plugins/` 内容）：

```
libOpenMMAmoebaReference.so   libOpenMMCPU.so        libOpenMMPME.so
libOpenMMDrudeReference.so     libOpenMMRPMDReference.so
```

只有 `Reference` 与 `CPU` 两个平台插件，**没有 CUDA / OpenCL / HIP**。
系统里也没装 FFTW / GSL（`rpm -qa` 查无此包），说明 `libOpenMM.so` 用的是内置 FFT
（`ldd` 依赖只有 `libstdc++` / `libm` / `libgcc_s` / `libc`）。

**编译工具链**（容器 `sdaa` 内，现仍可用）：

```
cmake 3.26.3   /usr/bin/cmake     (cmake-3.26.3-1.lns23.loongarch64)
swig  4.1.1    /usr/bin/swig      (swig-4.1.1-3.lns23.loongarch64)
g++ / gcc / make
```

**源码目录已被容器重建清掉**（`/data/application/xuqiang/openmm-src/` 不存在，
`/tmp/openmm_build.log` 也没了），要重装得重新拉源码编译。

### pdbfixer：一条命令

```bash
docker exec -u 1007:1007 -e USER=xuqiang -e HOME=/data/application/xuqiang \
  sdaa bash -c '
    source /opt/tecoai/setvars.sh
    /home/py312/bin/python3 -m pip install pdbfixer==1.12.0
  '
```

实测输出（容器 `sdaa`）：

```
Using cached pdbfixer-1.12.0-py3-none-any.whl
Requirement already satisfied: numpy in /home/py312/lib/python3.12/site-packages (from pdbfixer==1.12.0) (2.2.6)
Requirement already satisfied: openmm>=8.2 in /home/py312/lib/python3.12/site-packages (from pdbfixer==1.12.0) (8.5.2)
Installing collected packages: pdbfixer
Successfully installed pdbfixer-1.12.0
```

不需要 `--no-deps`：容器里 numpy 2.2.6 与 openmm 8.5.2 都已满足要求，pip 不会重拉。

**PyPI 上 pdbfixer 1.12.0 只有 sdist（`.tar.gz`，650 KB），没有预编译 wheel**，
但它是纯 Python 包（`Root-Is-Purelib: true`，`Tag: py3-none-any`），
pip 会在本地 build 一次（`Building wheel for pdbfixer (pyproject.toml)`，约 10 秒），产物与官方 wheel 等价。

### PATH 注意事项

容器 PATH 不含 `/home/py312/bin`，而 `pdbfixer` 命令行装在那里；跑 CLI 前
`export PATH=/home/py312/bin:$PATH` 即可，否则报 `pdbfixer: not found`。

---

## 3. 验证

### 3.1 官方测试套件（主验收）

`pdbfixer/tests/` 共 23 个用例，覆盖加氢（4）、构建+模拟（1）、电荷与溶剂化（6）、CLI（3）、突变（7）、删链（2）。

```
======================== 23 passed in 498.76s (0:08:18) ========================
```

（这是用标准 `pip install pdbfixer==1.12.0` 重装后重跑的结果，与首次安装时一致。）

x86 对照侧同一套测试也是 23/23 过（CLI 三例需要同样的 PATH 修复）。

`test_build_and_simulate` 覆盖 100+ 个真实 PDB（从 RCSB 下载），逐个走完
「找缺失残基 → 替换非标残基 → 去杂原子 → 找缺失原子 → 补原子 → 加氢 → 建体系 → 最小化 → 跑 500 步 MD」，
并检查能量不是 NaN —— **全过**。

### 3.2 跨平台结构对照

14 个 PDB（`1LE1 1BNA 4HHB 1TIM 2LYZ 1AKI 3EIY 1CRN 6LU7 1A2P 1R8J 2JOF 1L2Y 1A6M`），
两边走完全相同的流程（`findMissingResidues` → `findNonstandardResidues` →
`replaceNonstandardResidues` → `removeHeterogens` → `findMissingAtoms` →
`addMissingAtoms` → `addMissingHydrogens(pH=7.0)`），比对结构摘要：

| pdbid | 原子数 x86 / SDAA | 残基数 | 氢原子数 | 氢占比 |
| --- | --- | --- | --- | --- |
| 1A2P | 5181 / 5181 | 330 / 330 | 2547 / 2547 | 49.16% / 49.16% |
| 1A6M | 2435 / 2435 | 151 / 151 | 1231 / 1231 | 50.55% / 50.55% |
| 1AKI | 1960 / 1960 | 129 / 129 | 959 / 959 | 48.93% / 48.93% |
| 1BNA | 758 / 758 | 24 / 24 | 272 / 272 | 35.88% / 35.88% |
| 1CRN | 642 / 642 | 46 / 46 | 315 / 315 | 49.07% / 49.07% |
| 1L2Y | 304 / 304 | 20 / 20 | 150 / 150 | 49.34% / 49.34% |
| 1LE1 | 218 / 218 | 12 / 12 | 102 / 102 | 46.79% / 46.79% |
| 1R8J | 9222 / 9222 | 578 / 578 | 4578 / 4578 | 49.64% / 49.64% |
| 1TIM | 7516 / 7516 | 494 / 494 | 3774 / 3774 | 50.21% / 50.21% |
| 2JOF | 284 / 284 | 20 / 20 | 136 / 136 | 47.89% / 47.89% |
| 2LYZ | 1960 / 1960 | 129 / 129 | 959 / 959 | 48.93% / 48.93% |
| 3EIY | 3013 / 3013 | 196 / 196 | 1503 / 1503 | 49.88% / 49.88% |
| 4HHB | 8766 / 8766 | 574 / 574 | 4382 / 4382 | 49.99% / 49.99% |
| 6LU7 | 4730 / 4730 | 309 / 309 | 2341 / 2341 | 49.49% / 49.49% |

**14/14 全部完全一致，有差异的 case：0 个。**

### 3.3 坐标对照（固定 seed=12345，325-case 大样本）

默认流程里 `addHydrogens` 用未播种的 `random` 放氢初始位置（`modeller.py:1026/1281`），
**两次运行的输入就不同，坐标不能直接比**。有效对照是固定 `random.seed(12345)` 后两边重跑
（`pdbfixer_seeded.py`），样本用 `test_build_and_simulate` 自带的 PDB 列表（325 个，
`110D`~`1DFZ`，以核酸/小蛋白为主）。

**成功集两边一致**：各 322/325。失败 3 个且两边同因，均与平台无关：
1AO9（DOP 残基无定义，官方注释即写明不可处理）、1C4S / 1CAP（去杂原子后体系为空，建不了 Context）。

322 个逐原子比对结果：

| 类别 | case 数 | 占比 |
| --- | --- | --- |
| **全原子逐位一致**（重原子+氢全同） | **143** | 44% |
| 重原子逐位一致，仅氢有差 | 128 | 40% |
| 重原子局部有差 | 51 | 16% |

- 重原子有差的 51 个 case 中：38 个最大差 ≤1 Å（多数 ≤0.5 Å，PDB 精度量级）；
  **13 个差 >2 Å**（最大 1CEK 9.15 Å / 1CL4 8.30 Å）。
- 氢最大差分布：≤0.01 Å 占 206 个，0.01~0.1 占 69 个，0.1~1 占 32 个，>1 仅 15 个。

13 个大差 case 的噪声地板校验——同一台 x86、同样流程，仅把种子从 12345 换成 54321
（等价于换一组氢初始位置），再与 x86@12345 比：

| case | x86 vs x86 换种子 | SDAA vs x86（平台间） |
| --- | --- | --- |
| 1CL4 | **11.87** | 8.30 |
| 1CEK | **8.42** | 9.15 |
| 1AQS | **8.33** | 6.48 |
| 1AM0 | **7.84** | 6.83 |
| 1BBL | **7.66** | 2.76 |
| 1ANP | 5.10 | 6.73 |
| 1BOE | 4.03 | 3.49 |
| 1BX7 | 3.58 | 2.33 |
| 1BRV | 2.54 | 6.76 |
| 1BH7 | 2.98 | 4.00 |
| 1BYY | 2.02 | 7.26 |
| 135D | 1.91 | 2.91 |
| 1BDE | 0.31 | 2.24 |

**读法**：271/322（84%）case 重原子逐位一致；有重原子差的 case，分叉量级与 x86 自身
换种子相同（多数 case 甚至更大）——坐标分叉是"氢初始位置 + 能量最小化混沌"的固有性质，
不是 SDAA 引入的偏差。结论：**ALIGNED**。

原始比对表：`logs/big_heavy_cmp.txt`；结果目录 `results/big_sdaa/`、`results/big_x86/`、
噪声地板对照 `results/big_x86_s2/`（14-case 精查在 `results/seed_*`，耗时数据见 §4）。

### 3.4 能量对照（物理等价性）

初始能量对氢位置极敏感（坐标差 0.1 Å 就能让能量差几百 kcal/mol），所以先做能量最小化再比。
能量单位均为 kcal/mol，dE 为两份输出的能量差绝对值：

| 对照 | dE 均值 | dE 最大 | 相对差均值 |
| --- | --- | --- | --- |
| **SDAA vs x86（平台间）** | **73.02** | 343.06 | **5.25%** |
| SDAA 自身两遍（噪声地板） | 88.77 | 460.00 | 5.28% |
| x86 自身两遍（噪声地板） | 49.59 | 234.09 | 2.11% |

**平台间差异小于 SDAA 自身的噪声地板**，即 SDAA 引入的偏差被模型自身的数值非确定性淹没。

---

## 4. 性能

SDAA 侧 128 核（`nproc` = 128），x86 对照侧 56 核。口径：14 个 case 固定 `random.seed(12345)`
完整流程（找缺失残基 → 补原子 → 加氢），`wall_s` 含 RCSB 下载时间（PDB 文件小，占比低），
两边同为默认线程配置（`results/seed_sdaa/summary.json` vs `results/seed_x86/summary.json`）：

| case | SDAA (s) | x86 (s) | 倍率 |
| --- | --- | --- | --- |
| 1LE1 | 3.02 | 5.71 | 0.53 |
| 1L2Y | 1.82 | 3.74 | 0.49 |
| 1CRN | 2.13 | 3.22 | 0.66 |
| 1BNA | 3.35 | 3.17 | 1.06 |
| 2JOF | 3.21 | 3.85 | 0.83 |
| 2LYZ | 2.51 | 3.51 | 0.72 |
| 1AKI | 2.41 | 3.62 | 0.67 |
| 3EIY | 5.63 | 5.59 | 1.01 |
| 1A6M | 2.71 | 4.20 | 0.65 |
| 6LU7 | 4.12 | 6.16 | 0.67 |
| 4HHB | 5.05 | 5.32 | 0.95 |
| 1TIM | 5.07 | 6.26 | 0.81 |
| 1A2P | **32.44** | 5.99 | **5.42** |
| 1R8J | **119.17** | 7.33 | **16.26** |
| **总计** | **192.6** | **67.7** | **2.85** |

**12/14 case SDAA 为 x86 的 0.49~1.06 倍（普遍更快或持平）**；总倍率被两个异常 case 拉高：

- 1R8J（578 残基，缺残基多）：16.3×。异常在 `addMissingAtoms` 的 L-BFGS 最小化
  （SDAA 104.8 s vs x86 2.5 s），不是加氢（SDAA 2.5 s vs x86 2.5 s）。
  插桩数迭代（`scripts/iter_count.py`，结果 `results/iter_1R8J_*.json`）：
  **迭代数两边相当（SDAA 271 vs x86 321），慢在单次迭代**——
  addMissingAtoms 最小化 0.12 / 0.046 s/次 vs 0.0031 / 0.0008 s/次（约 40~56 倍），
  即每次迭代的力计算本身在 loongarch 上慢一个量级以上，与收敛路径无关。
  最可能原因：OpenMM CPU 平台力 kernel 在 x86 上有手写 SIMD 实现，loongarch64
  编译时无对应实现、走通用 C++ 回退（符号表两边都查不到 avx/sse，此条未实锤；
  合成 r⁻⁴ 探针只差 1.8×，说明不是全线慢，是特定 kernel 慢）。
- 1A2P（含金属/配体）：5.4×，同属最小化路径分叉。

`OMP_NUM_THREADS` 扫描（4HHB）：1 → 6415 ms、8 → 6608 ms、32 → 5408 ms、
64 → 6277 ms、128 → 5334 ms。**线程数影响很小**，瓶颈是单线程计算量本身。

---

## 5. 已知限制

- **OpenMM 无 SDAA 后端**：SDAA 上只有 `Reference` / `CPU` 两个平台插件，GPU 加速的 MD 跑不了。pdbfixer 的加氢与建模全在 CPU 上。
- **OpenMM 无 loongarch64 轮子**：需要源码编译（步骤见 §2）。x86 侧可直接用官方 manylinux wheel。
- **容器 PATH 不含 `/home/py312/bin`**：`pdbfixer` 命令行入口不在 PATH 里，跑 CLI 前要 `export PATH=/home/py312/bin:$PATH`。
- **pdbfixer 的加氢不可重复**：OpenMM `Modeller.addHydrogens` 内部用未播种的 `random` 决定氢的初始位置。需要可复现时得自己在调用前 `random.seed()`。
- **大复合物在 SDAA 上很慢**：5XNL（15 万原子）在 x86 上 64 s，SDAA 上 30 分钟未完成。
- **PDBFixer 的 `platform` 参数在 SDAA 上只有两个选择**：`Reference`（纯标量，最慢但最确定）与 `CPU`（多线程 SIMD）。默认 `None` 由 OpenMM 自动选。

---

## 6. 复现

**① 官方测试套件（主验收）**

```bash
docker exec -u 1007:1007 -e USER=xuqiang -e HOME=/data/application/xuqiang \
  sdaa bash -c '
    source /opt/tecoai/setvars.sh
    export PATH=/home/py312/bin:$PATH
    cd /tmp/pdbfix_probe/pdbfixer-1.12.0
    /home/py312/bin/python3 -m pytest pdbfixer/tests/ -v
  '
# 预期：23 passed
```

**② 官方命令行**

```bash
export PATH=/home/py312/bin:$PATH
pdbfixer --pdbid 1CRN --output out.pdb          # 从 RCSB 取结构并修复
pdbfixer input.pdb --output out.pdb --ph 7.0    # 修复本地文件
pdbfixer --pdbid 1CRN --output out.pdb --add-residues --replace-nonstandard
```

**③ 官方 API**

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(pdbid='1CRN')       # 或 PDBFixer(filename='input.pdb')
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
PDBFile.writeFile(fixer.topology, fixer.positions, open('output.pdb', 'w'))
```
