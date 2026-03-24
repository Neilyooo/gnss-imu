#include "gnss_imu_fusion/fusion.h"

#include "internal/estimator_common.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <stdexcept>

namespace gnss_imu {

namespace {

double VectorRmse(double squared_sum, std::size_t count) {
  if (count == 0) {
    return 0.0;
  }
  return std::sqrt(squared_sum / static_cast<double>(count));
}

std::optional<GroundTruthData> InterpolateGroundTruth(
    const AlignedDeque<GroundTruthData>& ground_truth,
    double timestamp) {
  if (ground_truth.empty()) {
    return std::nullopt;
  }

  auto it = std::lower_bound(
      ground_truth.begin(),
      ground_truth.end(),
      timestamp,
      [](const GroundTruthData& sample, double value) { return sample.timestamp < value; });

  if (it == ground_truth.begin()) {
    return *it;
  }
  if (it == ground_truth.end()) {
    return ground_truth.back();
  }

  const auto& upper = *it;
  const auto& lower = *(it - 1);
  const double span = upper.timestamp - lower.timestamp;
  if (span < 1e-9) {
    return lower;
  }

  const double alpha = (timestamp - lower.timestamp) / span;
  GroundTruthData interpolated;
  interpolated.timestamp = timestamp;
  interpolated.position = (1.0 - alpha) * lower.position + alpha * upper.position;
  interpolated.velocity = (1.0 - alpha) * lower.velocity + alpha * upper.velocity;
  interpolated.orientation = lower.orientation.slerp(alpha, upper.orientation);
  return interpolated;
}

BenchmarkMetrics ComputeMetrics(const std::string& algorithm,
                                const AlignedVector<TrajectorySample>& trajectory,
                                const AlignedDeque<GroundTruthData>& ground_truth,
                                double total_compute_ms,
                                double avg_update_ms,
                                double max_update_ms,
                                std::size_t updates,
                                double dataset_duration_seconds) {
  BenchmarkMetrics metrics;
  metrics.algorithm = algorithm;
  metrics.total_compute_ms = total_compute_ms;
  metrics.avg_update_ms = avg_update_ms;
  metrics.max_update_ms = max_update_ms;
  metrics.updates = updates;
  metrics.cpu_utilization_percent =
      dataset_duration_seconds > 0.0
          ? (total_compute_ms / 1000.0) / dataset_duration_seconds * 100.0
          : 0.0;

  if (ground_truth.empty()) {
    return metrics;
  }

  double position_squared_error = 0.0;
  double velocity_squared_error = 0.0;
  std::size_t count = 0;
  for (const auto& sample : trajectory) {
    const auto truth = InterpolateGroundTruth(ground_truth, sample.timestamp);
    if (!truth.has_value()) {
      continue;
    }
    position_squared_error +=
        (sample.state.position - truth->position).squaredNorm();
    velocity_squared_error +=
        (sample.state.velocity - truth->velocity).squaredNorm();
    ++count;
  }

  metrics.rmse_position = VectorRmse(position_squared_error, count);
  metrics.rmse_velocity = VectorRmse(velocity_squared_error, count);
  return metrics;
}

}  // namespace

std::vector<BenchmarkMetrics> RunBenchmark(const Dataset& dataset,
                                           const FusionConfig& config,
                                           const std::vector<std::string>& algorithms,
                                           const std::string& output_dir) {
  if (dataset.gnss.empty()) {
    throw std::invalid_argument("Benchmark requires at least one GNSS sample");
  }

  std::filesystem::create_directories(output_dir);

  if (!dataset.ground_truth.empty()) {
    AlignedVector<TrajectorySample> gt_trajectory;
    gt_trajectory.reserve(dataset.ground_truth.size());
    for (const auto& sample : dataset.ground_truth) {
      NavState state;
      state.timestamp = sample.timestamp;
      state.position = sample.position;
      state.velocity = sample.velocity;
      state.orientation = sample.orientation;
      gt_trajectory.push_back(TrajectorySample{sample.timestamp, state});
    }
    WriteTrajectoryCsv(output_dir + "/ground_truth.csv", gt_trajectory);
  }

  const double dataset_duration =
      dataset.gnss.back().timestamp - dataset.gnss.front().timestamp;

  std::vector<BenchmarkMetrics> metrics;
  metrics.reserve(algorithms.size());
  for (const auto& algorithm : algorithms) {
    auto estimator = CreateEstimator(algorithm, config);
    DataManager manager;

    std::size_t imu_index = 0;
    std::size_t updates = 0;
    double total_compute_ms = 0.0;
    double max_update_ms = 0.0;

    for (const auto& gnss : dataset.gnss) {
      while (imu_index < dataset.imu.size() &&
             dataset.imu[imu_index].timestamp <= gnss.timestamp + internal::kEps) {
        manager.PushImu(dataset.imu[imu_index]);
        ++imu_index;
      }

      const auto group = manager.PushGnss(gnss);
      if (!group.has_value()) {
        continue;
      }

      const auto start = std::chrono::steady_clock::now();
      estimator->Process(*group);
      const auto end = std::chrono::steady_clock::now();
      const double elapsed_ms =
          std::chrono::duration<double, std::milli>(end - start).count();
      total_compute_ms += elapsed_ms;
      max_update_ms = std::max(max_update_ms, elapsed_ms);
      ++updates;
    }

    const double avg_update_ms =
        updates > 0 ? total_compute_ms / static_cast<double>(updates) : 0.0;
    const std::string normalized = internal::ToLower(algorithm);
    WriteTrajectoryCsv(output_dir + "/" + normalized + "_trajectory.csv",
                       estimator->Trajectory());
    metrics.push_back(ComputeMetrics(normalized,
                                     estimator->Trajectory(),
                                     dataset.ground_truth,
                                     total_compute_ms,
                                     avg_update_ms,
                                     max_update_ms,
                                     updates,
                                     dataset_duration));
  }

  WriteMetricsCsv(output_dir + "/metrics.csv", metrics);
  return metrics;
}

}  // namespace gnss_imu
