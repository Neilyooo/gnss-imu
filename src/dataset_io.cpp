#include "gnss_imu_fusion/fusion.h"

#include "internal/estimator_common.h"

#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>

namespace gnss_imu {

namespace {

std::vector<double> ParseNumericRow(const std::string& line) {
  std::vector<double> values;
  std::stringstream stream(line);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (token.empty()) {
      continue;
    }
    values.push_back(std::stod(token));
  }
  return values;
}

}  // namespace

Dataset GenerateSyntheticDataset(const FusionConfig& config,
                                 double duration_seconds,
                                 double imu_rate_hz,
                                 double gnss_rate_hz,
                                 std::uint32_t seed) {
  Dataset dataset;
  if (duration_seconds <= 0.0 || imu_rate_hz <= 0.0 || gnss_rate_hz <= 0.0) {
    throw std::invalid_argument("Synthetic dataset parameters must be positive");
  }

  std::mt19937 generator(seed);
  std::normal_distribution<double> accel_noise(0.0, config.accel_noise_density);
  std::normal_distribution<double> gyro_noise(0.0, config.gyro_noise_density);
  std::normal_distribution<double> gnss_pos_noise(0.0, config.gnss_position_std);
  std::normal_distribution<double> gnss_vel_noise(0.0, config.gnss_velocity_std);
  std::normal_distribution<double> random_walk(0.0, 1.0);

  const double imu_dt = 1.0 / imu_rate_hz;
  const double gnss_dt = 1.0 / gnss_rate_hz;

  NavState truth;
  truth.orientation = Eigen::Quaterniond::Identity();
  truth.position = Vec3::Zero();
  truth.velocity = Vec3::Zero();

  Vec3 accel_bias = Vec3::Zero();
  Vec3 gyro_bias = Vec3::Zero();
  double next_gnss_time = 0.0;

  for (double t = 0.0; t <= duration_seconds + internal::kEps; t += imu_dt) {
    GroundTruthData truth_sample;
    truth_sample.timestamp = t;
    truth_sample.position = truth.position;
    truth_sample.velocity = truth.velocity;
    truth_sample.orientation = truth.orientation;
    dataset.ground_truth.push_back(truth_sample);

    const Vec3 true_gyro(
        0.03 * std::sin(0.41 * t),
        0.04 * std::cos(0.27 * t),
        0.18 + 0.05 * std::sin(0.17 * t));
    const Vec3 true_specific_force(
        0.6 + 0.25 * std::sin(0.13 * t),
        0.18 * std::cos(0.09 * t),
        0.10 * std::sin(0.31 * t));

    ImuData imu;
    imu.timestamp = t;
    imu.accel =
        true_specific_force + accel_bias +
        Vec3(accel_noise(generator), accel_noise(generator), accel_noise(generator));
    imu.gyro =
        true_gyro + gyro_bias +
        Vec3(gyro_noise(generator), gyro_noise(generator), gyro_noise(generator));
    dataset.imu.push_back(imu);

    if (t + internal::kEps >= next_gnss_time) {
      GnssData gnss;
      gnss.timestamp = t;
      gnss.position =
          truth.position +
          Vec3(gnss_pos_noise(generator), gnss_pos_noise(generator), gnss_pos_noise(generator));
      gnss.velocity =
          truth.velocity +
          Vec3(gnss_vel_noise(generator), gnss_vel_noise(generator), gnss_vel_noise(generator));
      gnss.position_covariance =
          Mat3::Identity() * config.gnss_position_std * config.gnss_position_std;
      gnss.velocity_covariance =
          Mat3::Identity() * config.gnss_velocity_std * config.gnss_velocity_std;
      gnss.has_velocity = true;
      dataset.gnss.push_back(gnss);
      next_gnss_time += gnss_dt;
    }

    accel_bias += Vec3(random_walk(generator), random_walk(generator), random_walk(generator)) *
                  config.accel_bias_random_walk * std::sqrt(imu_dt);
    gyro_bias += Vec3(random_walk(generator), random_walk(generator), random_walk(generator)) *
                 config.gyro_bias_random_walk * std::sqrt(imu_dt);

    truth.orientation =
        internal::NormalizeQuaternion(truth.orientation * internal::ExpQuaternion(true_gyro * imu_dt));
    const Vec3 accel_world = truth.orientation * true_specific_force + internal::Gravity(config);
    truth.position += truth.velocity * imu_dt + 0.5 * accel_world * imu_dt * imu_dt;
    truth.velocity += accel_world * imu_dt;
    truth.timestamp = t + imu_dt;
  }

  return dataset;
}

Dataset LoadDatasetFromCsv(const std::string& imu_path,
                           const std::string& gnss_path,
                           const std::optional<std::string>& gt_path) {
  Dataset dataset;

  {
    std::ifstream imu_file(imu_path);
    if (!imu_file.is_open()) {
      throw std::runtime_error("Failed to open IMU csv: " + imu_path);
    }
    std::string line;
    while (std::getline(imu_file, line)) {
      if (line.empty()) {
        continue;
      }
      try {
        const auto values = ParseNumericRow(line);
        if (values.size() < 7) {
          continue;
        }
        ImuData imu;
        imu.timestamp = values[0];
        imu.accel = Vec3(values[1], values[2], values[3]);
        imu.gyro = Vec3(values[4], values[5], values[6]);
        dataset.imu.push_back(imu);
      } catch (const std::exception&) {
        continue;
      }
    }
  }

  {
    std::ifstream gnss_file(gnss_path);
    if (!gnss_file.is_open()) {
      throw std::runtime_error("Failed to open GNSS csv: " + gnss_path);
    }
    std::string line;
    while (std::getline(gnss_file, line)) {
      if (line.empty()) {
        continue;
      }
      try {
        const auto values = ParseNumericRow(line);
        if (values.size() < 4) {
          continue;
        }
        GnssData gnss;
        gnss.timestamp = values[0];
        gnss.position = Vec3(values[1], values[2], values[3]);
        if (values.size() >= 7) {
          gnss.velocity = Vec3(values[4], values[5], values[6]);
          gnss.has_velocity = true;
        }
        const double pos_std = values.size() >= 8 ? values[7] : 1.5;
        const double vel_std = values.size() >= 9 ? values[8] : 0.3;
        gnss.position_covariance = Mat3::Identity() * pos_std * pos_std;
        gnss.velocity_covariance = Mat3::Identity() * vel_std * vel_std;
        dataset.gnss.push_back(gnss);
      } catch (const std::exception&) {
        continue;
      }
    }
  }

  if (gt_path.has_value()) {
    std::ifstream gt_file(*gt_path);
    if (!gt_file.is_open()) {
      throw std::runtime_error("Failed to open ground-truth csv: " + *gt_path);
    }
    std::string line;
    while (std::getline(gt_file, line)) {
      if (line.empty()) {
        continue;
      }
      try {
        const auto values = ParseNumericRow(line);
        if (values.size() < 7) {
          continue;
        }
        GroundTruthData sample;
        sample.timestamp = values[0];
        sample.position = Vec3(values[1], values[2], values[3]);
        sample.velocity = Vec3(values[4], values[5], values[6]);
        dataset.ground_truth.push_back(sample);
      } catch (const std::exception&) {
        continue;
      }
    }
  }

  return dataset;
}

void WriteTrajectoryCsv(const std::string& path,
                        const AlignedVector<TrajectorySample>& trajectory) {
  std::ofstream output(path);
  if (!output.is_open()) {
    throw std::runtime_error("Failed to open output trajectory: " + path);
  }
  output << std::setprecision(10);
  output << "timestamp,px,py,pz,vx,vy,vz,qw,qx,qy,qz,bax,bay,baz,bgx,bgy,bgz\n";
  for (const auto& sample : trajectory) {
    output << sample.timestamp << ','
           << sample.state.position.x() << ','
           << sample.state.position.y() << ','
           << sample.state.position.z() << ','
           << sample.state.velocity.x() << ','
           << sample.state.velocity.y() << ','
           << sample.state.velocity.z() << ','
           << sample.state.orientation.w() << ','
           << sample.state.orientation.x() << ','
           << sample.state.orientation.y() << ','
           << sample.state.orientation.z() << ','
           << sample.state.accel_bias.x() << ','
           << sample.state.accel_bias.y() << ','
           << sample.state.accel_bias.z() << ','
           << sample.state.gyro_bias.x() << ','
           << sample.state.gyro_bias.y() << ','
           << sample.state.gyro_bias.z() << '\n';
  }
}

void WriteMetricsCsv(const std::string& path,
                     const std::vector<BenchmarkMetrics>& metrics) {
  std::ofstream output(path);
  if (!output.is_open()) {
    throw std::runtime_error("Failed to open output metrics: " + path);
  }
  output << std::setprecision(10);
  output << "algorithm,rmse_position,rmse_velocity,cpu_utilization_percent,total_compute_ms,"
            "avg_update_ms,max_update_ms,updates\n";
  for (const auto& metric : metrics) {
    output << metric.algorithm << ','
           << metric.rmse_position << ','
           << metric.rmse_velocity << ','
           << metric.cpu_utilization_percent << ','
           << metric.total_compute_ms << ','
           << metric.avg_update_ms << ','
           << metric.max_update_ms << ','
           << metric.updates << '\n';
  }
}

std::string Usage() {
  return
      "Usage:\n"
      "  gnss_imu_bench [--mode synthetic|csv] [--algorithm all|eskf,ukf,fgo]\n"
      "                 [--output_dir path] [--duration seconds]\n"
      "                 [--imu_rate hz] [--gnss_rate hz] [--seed value]\n"
      "                 [--imu path/to/imu.csv] [--gnss path/to/gnss.csv]\n"
      "                 [--gt path/to/ground_truth.csv]\n"
      "\n"
      "CSV formats:\n"
      "  imu.csv  : timestamp,ax,ay,az,gx,gy,gz\n"
      "  gnss.csv : timestamp,px,py,pz[,vx,vy,vz][,pos_std][,vel_std]\n"
      "  gt.csv   : timestamp,px,py,pz,vx,vy,vz\n";
}

}  // namespace gnss_imu
