# GNSS-IMU Fusion Baseline

## 项目简介

这是一个纯 `CMake + C++17` 的 GNSS-IMU 融合基线工程，用统一的数据接口比较三种后端：

- `ESKF`：误差状态卡尔曼滤波
- `UKF`：无迹卡尔曼滤波
- `FGO`：固定滑窗因子图优化

当前版本的重点不是“把三个方案写在一个大文件里做分支判断”，而是把它们整理成独立对象、共享底层数据通路、共享状态定义与 IMU 传播模型，便于后续单独优化某一个方案而不影响其他模块。

详细设计和公式说明见 [docs/design.md](docs/design.md)。

## 方案介绍

### 1. ESKF

`ESKF` 使用统一的名义状态：

```text
x = [p, v, q, b_a, b_g]
```

以及 15 维误差状态：

```text
dx = [dp, dv, dtheta, dba, dbg]
```

特点：

- 用 IMU 中值积分做名义状态传播
- 用误差状态雅可比传播协方差
- GNSS 位置和可选速度做量测更新
- 姿态修正通过四元数指数映射注入

适用性：

- 结构清晰
- 计算量最低
- 是当前仓库里最适合作为工程基线的方案

### 2. UKF

`UKF` 与 `ESKF` 使用相同状态定义，但更新方式不同：

- 在 15 维切空间生成 sigma points
- 每个 sigma state 通过同一套 IMU 传播模型前向传播
- 用加权均值重建状态
- 用加权协方差完成量测更新

特点：

- 不依赖一阶线性化量测更新
- 对更强非线性场景更有表达力
- 计算量高于 `ESKF`

### 3. FGO

`FGO` 后端采用固定滑窗图优化思路：

- 每个 GNSS 时刻生成一个图节点
- 相邻 GNSS 节点之间建立 IMU 边
- 窗口内同时加入 GNSS 因子、IMU 因子和先验因子
- 窗口超长后，通过 Schur 补把最老节点边缘化为新的先验

当前实现特点：

- 采用 `Eigen` 自带线性代数实现 Gauss-Newton
- 不依赖 `Ceres` 或 `GTSAM`
- 是“低依赖、可编译、可跑”的固定滑窗 FGO 框架
- 目前仍是 GNSS 位置/速度观测图，不是严格的原始伪距/多普勒紧耦合图

## 代码框架

当前版本的代码结构如下：

```text
include/gnss_imu_fusion/
  fusion.h

src/
  fusion.cpp
  main.cpp
  benchmark.cpp
  data_manager.cpp
  dataset_io.cpp

  internal/
    estimator_common.h
    estimator_common.cpp

  estimators/
    estimator_base.h
    eskf_estimator.h
    eskf_estimator.cpp
    ukf_estimator.h
    ukf_estimator.cpp
    fgo_estimator.h
    fgo_estimator.cpp

scripts/
  run_fixposition_backtest.py

docs/
  design.md
  fixposition_backtest.md
  refactor_2026-03-20.md
  rigor_review_2026-03-16.md
```

模块职责：

- `include/gnss_imu_fusion/fusion.h`
  - 对外公共接口
  - 数据结构定义
  - benchmark API
- `src/fusion.cpp`
  - 后端工厂
  - 根据字符串创建 `ESKF / UKF / FGO` 对象
- `src/estimators/estimator_base.h`
  - 公共估计器生命周期
  - 通用状态保存与轨迹记录
- `src/estimators/*.h/.cpp`
  - 各自算法的独立实现
- `src/internal/estimator_common.*`
  - 公共状态流形运算
  - IMU 传播
  - 协方差传播
  - 数值稳定性辅助函数
- `src/data_manager.cpp`
  - IMU 缓冲
  - GNSS 到来时构造 `MeasureGroup`
- `src/dataset_io.cpp`
  - 合成数据生成
  - CSV 数据读写
- `src/benchmark.cpp`
  - 跑三种算法
  - 计算 RMSE 与耗时
  - 导出轨迹和指标

运行时数据流：

```text
IMU/GNSS input
  -> DataManager
  -> MeasureGroup
  -> FusionEstimator(eskf | ukf | fgo)
  -> trajectory / metrics csv
```

## 使用方法

### 1. 编译

```bash
cmake -S . -B build
cmake --build build -j
```

### 2. 合成数据测试

```bash
./build/gnss_imu_bench \
  --mode synthetic \
  --duration 10 \
  --algorithm all \
  --output_dir results/smoke
```

### 3. CSV 数据测试

```bash
./build/gnss_imu_bench \
  --mode csv \
  --imu data/imu.csv \
  --gnss data/gnss.csv \
  --gt data/ground_truth.csv \
  --algorithm eskf,ukf,fgo \
  --output_dir results
```

支持的 CSV 格式：

```text
imu.csv  : timestamp,ax,ay,az,gx,gy,gz
gnss.csv : timestamp,px,py,pz[,vx,vy,vz][,pos_std][,vel_std]
gt.csv   : timestamp,px,py,pz,vx,vy,vz
```

### 4. Fixposition rosbag 回放

```bash
python3 scripts/run_fixposition_backtest.py \
  --bag_dir /media/auto/新加卷/downloads/rosbag2_2026_01_05-11_14_35_0 \
  --output_dir results/fixposition_full
```

该脚本会：

- 读取 `/fixposition/rawimu`
- 读取 `/fixposition/gnss1` 与 `/fixposition/gnss2`
- 读取 `/fixposition/inspvax` 作为参考轨迹
- 将双天线 GNSS 融合成中点位置观测
- 用有限差分估计 GNSS 速度
- 调用本仓库的 `ESKF / UKF / FGO`
- 生成轨迹、误差图、KML 和总结报告

## 输出结果

典型输出包括：

- `results/metrics.csv`
- `results/eskf_trajectory.csv`
- `results/ukf_trajectory.csv`
- `results/fgo_trajectory.csv`
- `results/ground_truth.csv`

Fixposition 回放还会额外输出：

- `benchmark/summary.md`
- `benchmark/trajectory_xy.png`
- `benchmark/position_error.png`
- `benchmark/kml/*.kml`

## 核心函数

下面这些函数/类是当前框架的主干：

### 公共入口

- `CreateEstimator(const std::string&, const FusionConfig&)`
  - 按算法名称创建后端对象
- `RunBenchmark(const Dataset&, const FusionConfig&, const std::vector<std::string>&, const std::string&)`
  - 批量运行算法并导出评估结果

### 数据组织

- `DataManager::PushImu(const ImuData&)`
  - 推入高频 IMU
- `DataManager::PushGnss(const GnssData&)`
  - 推入 GNSS，并在时机成熟时生成一个 `MeasureGroup`
- `LoadDatasetFromCsv(...)`
  - 读取离线 CSV 数据集
- `GenerateSyntheticDataset(...)`
  - 生成用于烟测或回归的合成数据

### 公共数学与传播

- `InitializeStateFromGroup(...)`
  - 用首个 GNSS 和 IMU 均值初始化导航状态
- `PropagateNominal(...)`
  - 公共 IMU 名义状态传播
- `ContinuousTimeJacobian(...)`
  - 构造误差状态连续时间雅可比
- `UpdateCovarianceEskf(...)`
  - ESKF 协方差传播
- `PropagateEdgeLinearized(...)`
  - 对一段 IMU 边做传播并输出状态转移与噪声累积，供 UKF/FGO 使用
- `ApplyError(...)`
  - 将 15 维误差状态注入到名义状态
- `StateMinus(...)`
  - 计算两个状态在切空间下的差值

### 各后端核心逻辑

- `EskfEstimator::Process(const MeasureGroup&)`
  - ESKF 的一次预测-更新
- `UkfEstimator::Process(const MeasureGroup&)`
  - UKF 的一次 sigma 点传播与量测更新
- `FgoEstimator::Process(const MeasureGroup&)`
  - FGO 的节点扩展、优化、边缘化与回退控制


### 1. ESKF

- 参考其误差状态与四元数误差注入思想
- 但实现上做了工程化裁剪与数值稳健性处理
- 不是逐公式逐细节复现某一篇论文的实验配置

### 2. UKF

- 参考了 sigma point、均值/协方差重建的经典 UKF 思路
- 结合导航状态四元数表示做了切空间实现
- 加入了面向本仓库数据流的过程噪声累积方式

### 3. FGO

`FGO` 模块有更明确的直接参考背景，但仍不是严格论文复现。

当前仓库中明确提到的参考论文：

- *Real-time tightly coupled GNSS and IMU integration via Factor Graph Optimization*
- `arXiv:2603.03556v1`
- 发布时间：`2026-03-03`

本仓库对应关系：

- 参考了实时固定滑窗 FGO 的整体架构思路
- 参考了“GNSS/IMU 因子 + 边缘化先验 + 因子图实时优化”的框架方向
- 但当前输入只有 `NavSatFix`，没有原始伪距和多普勒
- 因此当前实现是 GNSS 位置/速度观测图，不是严格的紧耦合原始观测图

## 当前限制

- 假定工作坐标系是本地 ENU，重力沿 `-Z`
- `FGO` 当前是低依赖固定滑窗图优化，不是完整紧耦合原始 GNSS 因子图
- CPU 利用率统计方式为：

```text
algorithm_compute_time / dataset_duration
```

- 如果后续要继续往论文级紧耦合方向推进，至少还需要：
  - 原始伪距/多普勒观测建模
  - 接收机钟差与钟漂状态
  - 更完整的 GNSS 因子设计
  - 视需求切换到 `Ceres` 或 `GTSAM`

## 当前版本配套文档

- [docs/design.md](docs/design.md)：结构与公式总览
