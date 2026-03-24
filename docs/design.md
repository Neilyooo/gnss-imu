# GNSS-IMU Fusion Design

## 1. Current code structure

The project has been split into explicit modules so each algorithm is an independent object instead of being buried inside one monolithic source file.

```text
include/gnss_imu_fusion/
  fusion.h                         # public data types + benchmark API

src/
  fusion.cpp                       # algorithm factory
  data_manager.cpp                 # IMU buffering and GNSS grouping
  dataset_io.cpp                   # synthetic data + CSV loading/export
  benchmark.cpp                    # runtime benchmark and metrics export

  internal/
    estimator_common.h/.cpp        # shared state-manifold math and IMU propagation

  estimators/
    estimator_base.h               # common estimator lifecycle and trajectory storage
    eskf_estimator.h/.cpp          # ESKF backend
    ukf_estimator.h/.cpp           # UKF backend
    fgo_estimator.h/.cpp           # fixed-lag FGO backend
```

This separation is intentional:

- public API stays stable in `fusion.h`
- shared math stays in one internal layer and is not copy-pasted across algorithms
- each backend owns its own state, update rule, and internal buffers
- the factory only decides which estimator object to instantiate

## 2. Runtime data flow

The runtime path is shared up to the estimator boundary:

1. IMU samples are pushed into `DataManager`.
2. Each GNSS epoch closes one `MeasureGroup`.
3. `MeasureGroup` is passed to one `FusionEstimator`.
4. The estimator updates its internal state and appends one trajectory sample.
5. `benchmark.cpp` exports trajectories and metrics.

In short:

```text
IMU/GNSS input
  -> DataManager
  -> MeasureGroup
  -> FusionEstimator(eskf | ukf | fgo)
  -> Trajectory / metrics / csv
```

This keeps the hot path simple: the high-rate data path does not branch on algorithm internals until the final estimator call.

## 3. Common state convention

All three backends use the same nominal navigation state:

$$
\mathbf{x} =
\begin{bmatrix}
\mathbf{p} &
\mathbf{v} &
\mathbf{q} &
\mathbf{b}_a &
\mathbf{b}_g
\end{bmatrix}
$$

where:

- $\mathbf{p} \in \mathbb{R}^3$: position
- $\mathbf{v} \in \mathbb{R}^3$: velocity
- $\mathbf{q} \in SO(3)$: orientation quaternion
- $\mathbf{b}_a \in \mathbb{R}^3$: accelerometer bias
- $\mathbf{b}_g \in \mathbb{R}^3$: gyroscope bias

The shared 15D tangent error state is:

$$
\delta \mathbf{x} =
\begin{bmatrix}
\delta \mathbf{p} &
\delta \mathbf{v} &
\delta \boldsymbol{\theta} &
\delta \mathbf{b}_a &
\delta \mathbf{b}_g
\end{bmatrix}^{\top}
\in \mathbb{R}^{15}
$$

Using one state convention gives two direct benefits:

- all three algorithms can be benchmarked fairly on the same dataset
- propagation, error injection, and covariance utilities can be reused across modules

## 4. Shared IMU propagation model

All backends use the same midpoint IMU propagation in `src/internal/estimator_common.cpp`.

Given body-frame acceleration $\mathbf{f}$ and angular velocity $\boldsymbol{\omega}$:

$$
\mathbf{f}_m = \frac{1}{2}
\left[
\left(\mathbf{a}_{k-1} - \mathbf{b}_a\right) +
\left(\mathbf{a}_{k} - \mathbf{b}_a\right)
\right]
$$

$$
\boldsymbol{\omega}_m = \frac{1}{2}
\left[
\left(\boldsymbol{\omega}_{k-1} - \mathbf{b}_g\right) +
\left(\boldsymbol{\omega}_{k} - \mathbf{b}_g\right)
\right]
$$

$$
\mathbf{q}_{k+\frac{1}{2}} =
\mathbf{q}_k \otimes \exp\left(\frac{1}{2}\boldsymbol{\omega}_m \Delta t\right)
$$

$$
\mathbf{a}^w_m = \mathbf{R}(\mathbf{q}_{k+\frac{1}{2}})\mathbf{f}_m + \mathbf{g}
$$

$$
\mathbf{p}_{k+1} = \mathbf{p}_k + \mathbf{v}_k \Delta t + \frac{1}{2}\mathbf{a}^w_m \Delta t^2
$$

$$
\mathbf{v}_{k+1} = \mathbf{v}_k + \mathbf{a}^w_m \Delta t
$$

$$
\mathbf{q}_{k+1} = \mathbf{q}_k \otimes \exp\left(\boldsymbol{\omega}_m \Delta t\right)
$$

This common propagation layer is why later backend-specific optimization does not require changing the rest of the data pipeline.

## 5. ESKF backend

Source mapping:

- `src/estimators/eskf_estimator.h`
- `src/estimators/eskf_estimator.cpp`

### 5.1 Prediction

The ESKF keeps one nominal state and one covariance:

$$
\mathbf{P}_{k+1|k} = \mathbf{\Phi}_k \mathbf{P}_{k|k} \mathbf{\Phi}_k^\top + \mathbf{Q}_k
$$

with continuous-time linearization:

$$
\delta \dot{\mathbf{x}} = \mathbf{F}\delta \mathbf{x} + \mathbf{n}
$$

and second-order discretization:

$$
\mathbf{\Phi}_k \approx \mathbf{I} + \mathbf{F}\Delta t + \frac{1}{2}\left(\mathbf{F}\Delta t\right)^2
$$

### 5.2 GNSS update

The measurement model is:

$$
\mathbf{z}_k =
\begin{bmatrix}
\mathbf{p}^{gnss}_k \\
\mathbf{v}^{gnss}_k
\end{bmatrix}
,\quad
\mathbf{h}(\mathbf{x}_k) =
\begin{bmatrix}
\mathbf{p}_k \\
\mathbf{v}_k
\end{bmatrix}
$$

The Kalman gain is:

$$
\mathbf{K}_k = \mathbf{P}_{k|k-1}\mathbf{H}^\top
\left(\mathbf{H}\mathbf{P}_{k|k-1}\mathbf{H}^\top + \mathbf{R}\right)^{-1}
$$

The correction is injected on the manifold:

$$
\delta \hat{\mathbf{x}} = \mathbf{K}_k
\left(\mathbf{z}_k - \mathbf{h}(\mathbf{x}_{k|k-1})\right)
$$

$$
\mathbf{x}_{k|k} = \mathbf{x}_{k|k-1} \boxplus \delta \hat{\mathbf{x}}
$$

where $\boxplus$ means additive updates for Euclidean parts and quaternion exponential-map update for attitude.

## 6. UKF backend

Source mapping:

- `src/estimators/ukf_estimator.h`
- `src/estimators/ukf_estimator.cpp`

### 6.1 Sigma-point generation

The UKF uses the same 15D tangent space:

$$
\lambda = \alpha^2 (n + \kappa) - n
$$

$$
\mathbf{\chi}_0 = \bar{\mathbf{x}}
$$

$$
\mathbf{\chi}_i = \bar{\mathbf{x}} \boxplus
\left(+\sqrt{n+\lambda}\,\mathbf{L}_i\right)
$$

$$
\mathbf{\chi}_{i+n} = \bar{\mathbf{x}} \boxplus
\left(-\sqrt{n+\lambda}\,\mathbf{L}_i\right)
$$

where $\mathbf{L}$ is the Cholesky factor of the covariance.

### 6.2 Propagation and mean reconstruction

Each sigma state is propagated through the same IMU model:

$$
\mathbf{\chi}^{-}_i = f(\mathbf{\chi}_i, \mathbf{u}_{imu})
$$

The Euclidean components use weighted averaging; orientation uses iterative quaternion averaging:

$$
\bar{\mathbf{p}} = \sum_i w_i^{(m)} \mathbf{p}_i,\quad
\bar{\mathbf{v}} = \sum_i w_i^{(m)} \mathbf{v}_i
$$

The covariance is reconstructed in tangent space:

$$
\mathbf{P}^{-} =
\sum_i w_i^{(c)}
\left(\mathbf{\chi}^{-}_i \boxminus \bar{\mathbf{x}}^{-}\right)
\left(\mathbf{\chi}^{-}_i \boxminus \bar{\mathbf{x}}^{-}\right)^\top
+ \mathbf{Q}
$$

In this implementation, $\mathbf{Q}$ is accumulated from the actual IMU segment rather than approximated by one coarse GNSS-interval lump.

### 6.3 GNSS update

Measurement sigma points are:

$$
\mathbf{z}_i =
\begin{bmatrix}
\mathbf{p}_i \\
\mathbf{v}_i
\end{bmatrix}
$$

Then:

$$
\bar{\mathbf{z}} = \sum_i w_i^{(m)} \mathbf{z}_i
$$

$$
\mathbf{S} = \sum_i w_i^{(c)} (\mathbf{z}_i - \bar{\mathbf{z}})(\mathbf{z}_i - \bar{\mathbf{z}})^\top + \mathbf{R}
$$

$$
\mathbf{P}_{xz} = \sum_i w_i^{(c)}
\left(\mathbf{\chi}^{-}_i \boxminus \bar{\mathbf{x}}^{-}\right)
(\mathbf{z}_i - \bar{\mathbf{z}})^\top
$$

$$
\mathbf{K} = \mathbf{P}_{xz}\mathbf{S}^{-1}
$$

$$
\mathbf{x}^{+} = \bar{\mathbf{x}}^{-} \boxplus \mathbf{K}(\mathbf{z} - \bar{\mathbf{z}})
$$

## 7. FGO backend

Source mapping:

- `src/estimators/fgo_estimator.h`
- `src/estimators/fgo_estimator.cpp`

### 7.1 Graph structure

The graph contains:

- one node per GNSS epoch
- one IMU edge between adjacent nodes
- one prior factor for the sliding-window anchor

For a window of states $\mathcal{X} = \{\mathbf{x}_0, \dots, \mathbf{x}_N\}$, the cost is:

$$
J(\mathcal{X}) =
\|\mathbf{r}_{prior}\|^2_{\mathbf{\Lambda}_{prior}}
+ \sum_{k=0}^{N}\|\mathbf{r}^{gnss}_k\|^2_{\mathbf{\Lambda}^{gnss}_k}
+ \sum_{k=0}^{N-1}\|\mathbf{r}^{imu}_{k,k+1}\|^2_{\mathbf{\Lambda}^{imu}_{k,k+1}}
$$

with:

$$
\mathbf{r}^{gnss}_k =
\begin{bmatrix}
\mathbf{p}_k - \mathbf{p}^{gnss}_k \\
\mathbf{v}_k - \mathbf{v}^{gnss}_k
\end{bmatrix}
$$

$$
\mathbf{r}^{imu}_{k,k+1} =
\mathbf{x}_{k+1} \boxminus \hat{\mathbf{x}}_{k+1|k}
$$

where $\hat{\mathbf{x}}_{k+1|k}$ is the propagated state produced by the shared IMU model.

### 7.2 Gauss-Newton step

At each iteration the backend builds:

$$
\mathbf{H} = \sum_j \mathbf{J}_j^\top \mathbf{\Lambda}_j \mathbf{J}_j
$$

$$
\mathbf{b} = \sum_j \mathbf{J}_j^\top \mathbf{\Lambda}_j \mathbf{r}_j
$$

and solves:

$$
\mathbf{H}\Delta \mathbf{x} = -\mathbf{b}
$$

Then each node state is updated by manifold error injection.

The implementation also applies:

- Huber-style robust weighting on GNSS and IMU factors
- step scaling to prevent large unstable updates
- a post-optimization consistency gate against the current GNSS epoch

### 7.3 Marginalization

When the window exceeds `fgo_window_size`, the oldest node is marginalized by Schur complement.

Partition the normal equation as:

$$
\begin{bmatrix}
\mathbf{H}_{mm} & \mathbf{H}_{mr} \\
\mathbf{H}_{rm} & \mathbf{H}_{rr}
\end{bmatrix}
\begin{bmatrix}
\Delta \mathbf{x}_m \\
\Delta \mathbf{x}_r
\end{bmatrix}
\begin{bmatrix}
\mathbf{b}_m \\
\mathbf{b}_r
\end{bmatrix}
$$

Then the reduced prior becomes:

$$
\mathbf{H}_{sc} = \mathbf{H}_{rr} - \mathbf{H}_{rm}\mathbf{H}_{mm}^{-1}\mathbf{H}_{mr}
$$

$$
\mathbf{b}_{sc} = \mathbf{b}_{r} - \mathbf{H}_{rm}\mathbf{H}_{mm}^{-1}\mathbf{b}_{m}
$$

This reduced system is stored as the new prior on the remaining oldest node.

## 8. Why this framework is more maintainable

Compared with a monolithic fusion file, the current layout improves maintainability in four ways:

- each backend is now an explicit class and can be extended independently
- shared state-manifold math is centralized, so bug fixes land once instead of three times
- benchmark/data IO code no longer couples to algorithm-specific implementation details
- the factory is the only selection point, which keeps the rest of the code closed to branching

This is the structure needed for later work such as:

- swapping the FGO solver to `Ceres` or `GTSAM`
- adding a tightly coupled raw-GNSS backend without touching ESKF/UKF
- experimenting with different process-noise or initialization strategies per backend

## 9. Real-bag limitation

The Fixposition bag used in this repository contains `NavSatFix`, not raw pseudorange and Doppler.

So the current FGO backend is:

- a real-time fixed-lag GNSS-position graph
- not yet a strict tightly-coupled raw-GNSS factor graph

To reach that target, the next required steps are:

1. extend `GnssData` to raw satellite observations
2. add clock bias and drift to the graph state
3. replace the current GNSS residual with pseudorange and Doppler factors
4. evaluate whether Eigen Gauss-Newton is still sufficient or whether a dedicated solver is needed

## 10. Verification status

Current local verification on this workspace:

- `cmake --build build -j4`
- synthetic smoke benchmark through `gnss_imu_bench`

This means the refactored module structure remains buildable after separating the estimator classes.
