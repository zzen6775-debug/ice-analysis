# Ice Growth Analysis

从 GROMACS 轨迹中定量分析**冰/水两相体系的冰前沿位置与生长速率**，并导出可直接在 Origin 中作图的表格。

适用于板状/薄膜冰体系（一侧或两侧为过冷液态水），也可用于研究聚合物、离子、表面等对冰生长的影响。

> 全部脚本平铺在根目录，克隆即可用。

---

## 目录

- [快速开始](#快速开始)
- [完整流程](#完整流程)
- [两种生长量：thickness 与 ice amount](#两种生长量thickness-与-ice-amount)
- [文件说明](#文件说明)
- [方法学](#方法学)
- [Origin 作图](#origin-作图)
- [踩过的坑](#踩过的坑)
- [引用](#引用)

---

## 快速开始

仓库自带示例数据（4 个体系 × 2 ns 间隔），不用跑模拟就能验证环境：

```bash
pip install -r requirements.txt

# 用示例数据直接拟合
python3 growth_fit.py fronts_track.csv --windows 0 20 40 60
```

---

## 完整流程

### 1. 在集群上抽帧（每 2 ns，只取水氧）

`extract_frames2.sbatch` 从 `npt_260K.xtc` 里每 2000 ps 取一帧，只导出水氧原子（q(z) 只需要 O 的位置，体积缩小到 1/4）。

抽帧需要一个只含水氧的索引文件，用 `make_ow_index.py` 生成（拓扑顺序是 SOL 在前，OW 为第 1、5、9… 个原子）：

```bash
python3 make_ow_index.py                    # 生成 ow_<tag>.ndx
scp ow_*.ndx extract_frames2.sbatch USER@HOST:/scratch/USER/zsmx/
ssh USER@HOST "cd /scratch/USER/zsmx && sbatch extract_frames2.sbatch"
```

**注意**：必须提交到计算节点（`--partition=pub`），不要在登录节点跑 `trjconv`。

### 2. 取回并分析

```bash
scp USER@HOST:/scratch/USER/zsmx/_frames2/frames2.tar.gz .
bash run_pipeline.sh frames2.tar.gz
```

`run_pipeline.sh` 一路做完：解包 → 拆分 → 前沿定位 → 拟合 → 导出 Origin 表格。

---

## 两种生长量：thickness 与 ice amount

| 量 | 定义 | 优点 | 缺点 |
|---|---|---|---|
| `thickness` | q 首次超过阈值处的 z 上界 − 下界 | 直观，与"冰层厚度"对应 | 界面层被当成纯冰，且受分箱量化影响 |
| `ice_amount` | ∫(冰分子占比) dz | **完全连续，无量化误差** | 数值比 thickness 小约 25%（界面按比例计入） |

**两者都是合法的速率定义，但必须统一使用，不能混着比较。**

推荐把 `ice_amount` 作为主指标：它不需要定位前沿，因此对噪声和分箱完全不敏感，R² 通常 > 0.96。

示例结果（4 个体系，0–66 ns，每 2 ns）：

| 体系 | v_thickness (nm/ns) | v_ice_amount (nm/ns) |
|---|---|---|
| HEMA | 0.2050 ± 0.002 | 0.1521 |
| MMA | 0.1885 ± 0.001 | 0.1386 |
| PMMA10 | 0.2075 ± 0.003 | 0.1528 |
| PMMA20 | **0.1531 ± 0.006** | **0.1142** |

两个指标给出同一排序：PMMA20 慢 20–25%。

---

## 文件说明

| 文件 | 作用 |
|---|---|
| `ice_front.py` | 基础模块：读 `.gro`、算四面体序参量 q |
| `fronts_track.py` | **主分析**：0.2 nm 网格 + 锚点窗口搜索 + 插值，给出亚纳米前沿位置和 ice amount |
| `growth_fit.py` | 线性/分段拟合，`--windows` 出分窗口斜率 |
| `export_for_origin.py` | 生成 Origin 宽表 |
| `centre_fronts.py` | 把前沿位置平移到各自的冰板中心 |
| `split_frames.py` | 多帧 `.gro` 拆成单帧 |
| `make_ow_index.py` | 生成只含水氧的 GROMACS 索引文件 |
| `run_pipeline.sh` | 一键流程 |
| `extract_frames2.sbatch` | 集群端抽帧（每 2 ns，只取 OW） |
| `plot_thickness.ogs` | Origin：冰层厚度 |
| `plot_ice_amount.ogs` | Origin：冰量积分（推荐） |
| `plot_fronts_centered.ogs` | Origin：8 条前沿，**按冰板中心居中**（推荐） |
| `plot_fronts_absolute.ogs` | Origin：8 条前沿，绝对 z（仅限单体系内部比较） |
| `fronts_track.csv` | 示例输出 |
| `origin_fronts_centered.csv` | 示例 Origin 表（居中前沿） |
| `origin_ice_amount.csv` | 示例 Origin 表（冰量） |

---

## 方法学

### 四面体序参量

对每个水氧原子 $i$，取截断半径（0.35 nm）内的邻居，对所有邻居对 $(j,k)$ 的夹角求和：

$$
q_i = 1 - \frac{3}{8}\sum_{j<k}\left(\cos\theta_{jik} + \frac{1}{3}\right)^2
$$

| 状态 | q |
|---|---|
| 理想四面体 | 1.0 |
| 冰（TIP4P/Ice 实测） | 0.94 – 0.95 |
| 过冷液态水 | 0.75 – 0.80 |

冰相判据：q > 0.85，且该分箱内至少 8 个水分子。

### 前沿定位：锚点跟踪

1. 第一帧：取**最长的一段连续冰分箱**作为初始冰区
2. 后续帧：以上一帧的前沿位置为锚点，只在 **±1.2 nm 窗口**内搜索冰/水边界
3. 边界位置由 q 跨越 0.85 处的**线性插值**给出（亚 0.2 nm 精度）

这样既避免了全局搜索导致的"跳带"，又没有整纳米量化。

### 生长速率

对前沿位置或冰量做线性拟合，斜率即速率。下前沿斜率为负，取绝对值报告。

---

## Origin 作图

LabTalk 脚本**必须用完整字面路径**——X-Function 参数不支持 `a$ + b$` 拼接。

每个脚本顶部有两行路径需要改成你自己的数据目录：

```labtalk
impasc fname:="你的目录\origin_data\origin_thickness.csv";
expGraph type:=png path:="你的输出目录" filename:="ice_thickness" overwrite:=replace;
```

然后 Origin 里按 `Shift+F3` 打开命令窗口：

```
run.file("完整路径\plot_ice_amount.ogs")
```

### 手动作图

1. `文件 > 导入 > 逗号分隔(CSV)`（**必须激活工作簿窗口**，图窗口激活时没有"导入"）
2. 右键第一列 → `设置为 > X`
3. 选中所有 Y 列 → `绘图 > 基础2D图 > 线+符号图`
4. 双击轴标题手动改文字（默认会用第一个 Y 列的长名称）
5. `文件 > 导出图形`

---

## 断点续跑（当 mdrun 因为数值不稳定中断时）

长轨迹里偶尔会有 `Too many LINCS warnings` 导致 mdrun 退出。用 checkpoint 续跑即可，
但要遵守两条：

1. **续跑必须带 `-cpi <deffnm>.cpt -append`**，否则会从头开始，丢掉已有进度。
2. **参数不能改**（步长、约束、压浴、参考坐标文件），否则 grompp 会拒绝，或者物理状态被破坏。

一个实用的重试策略（见 `extract_frames2.sbatch` 同目录的思路，实际写在运行脚本里）：

```
最多 40 次尝试，每次从最新 checkpoint 续跑；
若连续 3 次各自推进 < 0.2 ns，则停止并报告——这个模式说明是真卡住，
不是运气不好，继续磨只会白烧机时。
```

### 续跑会不会破坏轨迹的连续性？

不会。GROMACS 的 checkpoint 按**步数**记录，而 `nstxout-compressed` 也是按步数，
续跑时帧网格天然对齐。用 `gmx check` 可以核实：

```bash
gmx check -f traj.xtc
# Item        #frames   Timestep (ps)
# Time          1799          10        <- 1799 帧 x 10 ps = 17980 ps
```

若 `#frames × 间隔` 恰好等于 `末时刻 − 首时刻`，说明**既没有缺失也没有重复**。
一本 17.98 ns 的轨迹经历了 6 次 checkpoint 续跑，帧数依然精确吻合。

## 踩过的坑

1. **整纳米台阶**：1 nm 分箱把前沿量化成整数，2 ns 间隔内只推进约 0.2 nm，导致多个点落在同一台阶上。解决办法是锚点跟踪 + 插值，而不是单纯细化分箱。

2. **细网格 + 全局"最长段"会跳带**：0.2 nm 网格下噪声放大，最长连续段会在候选之间切换，前沿位置突跳几纳米（实测出现过 30 nm 的跳变）。必须用锚点窗口约束搜索范围。

3. **索引基准错误最阴险**：从窗口切片取出的相对索引若当成绝对分箱号使用，前沿会永远落在盒子边缘（z ≈ 1.2 nm），而且表面上看还"有数据"。判断标准：每个体系的不同位置数应接近帧数，最大单步应远小于窗口宽度。

4. **绝对 z 不能跨体系比较**：降温后盒长不同（67.8–73.6 nm），同一个冰板落在不同绝对 z 上。必须先减去各自的冰板中心（见 `centre_fronts.py`）。

5. **两点/三点拟合的 R² 是假象**：`R² = 1.000` 往往只是点太少恰好共线。`growth_fit.py` 现在要求每段至少 3 个点。

6. **分段拟合总能找到"拐点"**：拐点位置是拟合参数，在量化数据上尤其不稳定。改用固定窗口斜率（`--windows`）判断趋势更可靠。

7. **`.gro` 是定宽格式**：原子数超过 99999 时原子名与编号粘连，必须按列切片，不能用 `split()`。

8. **`-pbc whole` 会把原子推出盒子**：KD 树会报 `Some input data are greater than the size of the periodic box`，脚本里已做 `% box` 折回。

9. **LabTalk 不支持字符串拼接**：`impasc fname:=root$ + "\a.csv"` 会报 `value + is missing a variable name`，必须写完整字面路径。

10. **不要从日志的 `Step/Time` 行推断轨迹是否重叠**：日志里最后一次 `Step/Time` 是**能量输出行**，不是崩溃点；续跑点可能比它晚一个帧间隔，看起来像"缺口"，其实完全连续。判断依据只能是 `gmx check -f traj.xtc` 的帧数与时间跨度是否自洽。

---

## 引用

- 四面体序参量：Chau P.-L., Hardwick A. J. *Mol. Phys.* **1998**, 93, 511.
- TIP4P/Ice 水模型：Abascal J. L. F. et al. *J. Chem. Phys.* **2005**, 122, 234511.

## 许可

MIT，见 `LICENSE`。
