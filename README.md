# Ice Growth Analysis

从 GROMACS 轨迹中定量分析**冰/水两相体系的冰前沿位置和生长速率**，并导出可直接在 Origin 中作图的表格。

适用于薄膜/板状冰体系（一侧或两侧为过冷液态水），也适用于研究聚合物、离子、表面等对冰生长的影响。

---

## 目录

- [功能](#功能)
- [环境要求](#环境要求)
- [完整流程](#完整流程)
- [文件说明](#文件说明)
- [方法学](#方法学)
- [输出文件](#输出文件)
- [Origin 作图](#origin-作图)
- [常见问题](#常见问题)
- [已知的坑](#已知的坑)
- [引用](#引用)

---

## 功能

| 能力 | 说明 |
|---|---|
| 抽帧 | 从 `.xtc` 轨迹中按指定时刻导出 `.gro` 构型 |
| 序参量 | 计算每个水氧原子的四面体序参量 q |
| 前沿定位 | 沿 z 方向扫描，确定下/上两个冰前沿的位置 |
| 生长速率 | 分别拟合下前沿、上前沿和总增厚速率 |
| 数据导出 | 生成 Origin 可直接导入的宽表 CSV |

---

## 环境要求

**本地分析（必需）**

- Python ≥ 3.9
- `numpy`
- `scipy`（用于周期边界 KD 树）

```bash
pip install -r requirements.txt
```

**集群抽帧（可选）**

- 任意带 GROMACS 和 Slurm 的集群
- 脚本默认适配 `gmx_mpi` + `pub` 分区

**作图（可选）**

- Origin 2018 或更新版本（脚本在 Origin 2024 中文版上验证）

---

## 完整流程

### 步骤 1：在集群上抽帧

编辑 `extract_frames.sbatch`，修改：

```bash
# 体系名:标签 的对照，按你的目录结构改
for pair in "HEMA_on_ice_slab:HEMA" "MMA_on_ice_slab:MMA" ...
# 抽取时刻（单位 ps）
for t in 0 5000 10000 20000
```

提交：

```bash
cd /scratch/YOUR_NAME/zsmx
sbatch extract_frames.sbatch
```

脚本会把所有帧写到 `zsmx/_frames/`，并打包成 `frames.tar.gz`。

**注意**：必须用计算节点（`--partition=pub`），不要在登录节点跑 `trjconv`。

### 步骤 2：取回本地

```bash
scp USER@HOST:/scratch/USER/zsmx/_frames/frames.tar.gz .
tar xzf frames.tar.gz -C frames/
```

### 步骤 3：计算前沿位置

```bash
python3 ice_front.py --frames frames/ --out fronts.csv
```

输出示例：

```
  HEMA       t=  0.00 ns  ice 32.0–35.0 nm (3.0 nm)  chains 26.84/40.26
  HEMA       t=  5.00 ns  ice 31.0–35.0 nm (4.0 nm)  chains 26.86/40.31
  HEMA       t= 10.00 ns  ice 31.0–36.0 nm (5.0 nm)  chains 26.81/40.17
  HEMA       t= 20.00 ns  ice 30.0–37.0 nm (7.0 nm)  chains 26.76/40.31
```

参数：

| 参数 | 默认 | 含义 |
|---|---|---|
| `--frames` | 必需 | 含 `<标签>_<时刻>ps.gro` 的目录 |
| `--out` | `fronts.csv` | 输出 CSV |
| `--bin` | `1.0` | z 方向分箱宽度（nm） |
| `--cut` | `0.35` | 第一壳层截断半径（nm） |
| `--qcut` | `0.85` | 判定为冰的 q 阈值 |

### 步骤 4：拟合生长速率

```bash
python3 growth_fit.py fronts.csv
```

输出：

```
system       v_total (nm/ns)   v_per_front     R^2   v_lower   v_upper   n
HEMA              0.2000±0.000      0.1000   1.000    0.0714    0.1286   3
MMA               0.1429±0.049      0.0714   0.893    0.0714    0.0714   3
```

- `v_total`：冰层总增厚速率
- `v_lower` / `v_upper`：下前沿、上前沿各自的推进速率（**推荐分别报告**）

其他选项：

```bash
python3 growth_fit.py fronts.csv --piecewise      # 尝试一个折点（前沿接触障碍物后）
python3 growth_fit.py fronts.csv --plot out.png   # 出图（需 matplotlib）
python3 growth_fit.py fronts.csv --tube 10        # 只拟合 t >= 10 ns
```

### 步骤 5：导出 Origin 表格

```bash
python3 export_for_origin.py fronts.csv --outdir origin_data
```

---

## 文件说明

| 文件 | 作用 |
|---|---|
| `ice_front.py` | 逐帧计算 q(z)、冰层范围、两个前沿位置、聚合物位置 |
| `growth_fit.py` | 从 CSV 拟合生长速率（线性 / 分段），可出图 |
| `export_for_origin.py` | 生成 Origin 可直接导入的宽表 |
| `run_analysis.sh` | 一键串起"拉帧 → 分析 → 拟合"（WSL/Linux） |
| `extract_frames.sbatch` | 集群端抽帧作业 |
| `origin/plot_ice.ogs` | Origin LabTalk：厚度曲线作图 + 导出 PNG |
| `origin/plot_fronts_both.ogs` | Origin LabTalk：上下两个前沿同时作图 |
| `sample_data/` | 示例数据（4 体系 × 4 时刻），可用来试跑 |

---

## 方法学

### 四面体序参量

对每个水氧原子 $i$，取截断半径内的所有邻居，计算所有邻居对 $(j,k)$ 的夹角：

$$
q_i = 1 - \frac{3}{8}\sum_{j<k}\left(\cos\theta_{jik} + \frac{1}{3}\right)^2
$$

| 状态 | q 值 |
|---|---|
| 理想四面体 | 1.0 |
| 冰（TIP4P/Ice，本体系实测） | 0.94 – 0.95 |
| 过冷液态水 | 0.75 – 0.80 |

### 冰相判据

某个 z 分箱内满足以下两条即视为冰：

1. 平均 q > `qcut`（默认 0.85）
2. 该分箱内至少 20 个水分子

### 前沿定义

沿 z 扫描所有分箱，取**最长的一段连续冰分箱**，其下边界为下前沿、上边界为上前沿。这样定义可以避免把远处孤立的六元环误判成冰。

### 生长速率

对前沿位置 $z(t)$ 做线性拟合，斜率即推进速率。下前沿斜率为负（向下推进），取绝对值报告。

---

## 输出文件

### `fronts.csv`

| 列 | 含义 |
|---|---|
| `system` | 体系标签 |
| `t_ns` | 时间（ns） |
| `n_waters` | 水分子数 |
| `box_z` | 盒子 z 长度（nm） |
| `ice_lo` / `ice_hi` | 下 / 上冰前沿位置（nm） |
| `thickness` | 冰层厚度 = `ice_hi - ice_lo` |
| `chainA_z` / `chainB_z` | 两条聚合物链的 z 中位数 |
| `gapA` / `gapB` | 链 A 到下前沿、链 B 到上前沿的距离 |

### `fronts.csv.qprofile`

| 列 | 含义 |
|---|---|
| `system`, `t_ns` | 体系与时刻 |
| `z_nm` | z 坐标 |
| `q` | 该分箱平均 q |
| `ice_fraction` | 该分箱 q > qcut 的分子占比 |

### `origin_data/`

| 文件 | 内容 |
|---|---|
| `origin_thickness.csv` | 冰层厚度 vs 时间 |
| `origin_front_lower.csv` | 下前沿位置 vs 时间 |
| `origin_front_upper.csv` | 上前沿位置 vs 时间 |
| `origin_fronts_both.csv` | 8 列：4 体系 × 上下前沿 |
| `origin_gap_chainA.csv` / `origin_gap_chainB.csv` | 聚合物到前沿的距离 |
| `origin_qprofile_<T>ns.csv` | 各时刻的 q(z) 剖面 |

---

## Origin 作图

### 方式 A：跑现成脚本（推荐）

1. 用记事本打开 `origin/plot_ice.ogs`，把文件开头的路径改成你自己的：

   ```labtalk
   string outdir$ = "你的数据目录\origin_output";
   ...
   impasc fname:="你的数据目录\origin_data\origin_thickness.csv";
   ```

2. Origin 中按 `Shift+F3` 打开**命令窗口**，输入：

   ```
   run.file("完整路径\plot_ice.ogs")
   ```

3. 脚本自动完成：新建工作簿 → 导入 → 设 X 列 → 线+符号图 → 设轴标题 → 导出 PNG。

### 方式 B：手动

1. `文件 > 导入 > 逗号分隔(CSV)`（**必须激活工作簿窗口**，图窗口激活时菜单里没有"导入"）
2. 右键第一列 → `设置为 > X`
3. 选中所有 Y 列 → `绘图 > 基础2D图 > 线+符号图`
4. 双击轴标题改成你要的文字（默认会用第一个 Y 列的长名称，容易误导）
5. `文件 > 导出图形`

---

## 常见问题

**Q：纵轴标题显示成了第一个体系的名字？**

Origin 默认用第一列的长名称作纵轴标题。用 LabTalk 覆盖：

```
xb.text$ = "t (ns)";
yl.text$ = "Ice thickness (nm)";
```

或直接双击轴标题手动改。

**Q：`plotxy iy:=2:5` 报错？**

语法要写成 `(1,2:5)`：X 取第 1 列，Y 取第 2–5 列。

**Q：`impasc` 报 `Failed to resolve tree node ... names.FNameToLabels`？**

新版 Origin 的导入选项树路径变了，去掉 `options.*` 参数只保留 `fname:=` 即可。

**Q：`yl.from = 0` 报 `illegal name`？**

2024 版不支持该属性名，删掉这行，在图上双击轴手动设范围。

**Q：一条 LabTalk 命令报错后，后面的命令都不执行？**

LabTalk 是顺序执行、遇错即停。建议先用单行命令逐个验证，再合成脚本。

---

## 已知的坑

1. **分辨率陷阱**：前沿推进速率约 0.1 nm/ns，1 nm 分箱下 **2 ns 内只有 0.2 nm 位移，看不出来**。至少需要 10 ns 的数据点才能分辨趋势。曾用 1.5 ns 数据误判为"前沿不动"。

2. **`-pbc whole` 让原子跑出盒子**：KD 树会报 `Some input data are greater than the size of the periodic box`。脚本已做 `% box` 折回处理。

3. **`.gro` 是定宽格式**：原子数超过 99999 时原子名会和编号粘连，必须按列切片（`[10:15]`、`[20:28]`…），不能用 `split()`。

4. **分箱边缘效应**：前沿恰在分箱边界时厚度会跳变整 1 nm，判断趋势要看多个时间点的单调性，不要看单点绝对值。

5. **上下前沿必须分开看**：两个前沿到聚合物的距离往往不同，合并成厚度曲线会把"一侧被封住、另一侧继续长"的效应平均掉。

6. **两点拟合没有意义**：`R² = 1.000` 常常只是点数太少（3–4 个点恰好共线），不要当成拟合良好的证据。

---

## 引用

如果这个方法对你的工作有帮助，请引用：

- 四面体序参量：Chau P.-L., Hardwick A. J. *Mol. Phys.* **1998**, 93, 511.
- TIP4P/Ice 水模型：Abascal J. L. F. et al. *J. Chem. Phys.* **2005**, 122, 234511.

---

## 许可

MIT License，见 `LICENSE`。
