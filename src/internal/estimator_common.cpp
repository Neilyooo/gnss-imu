#include "internal/estimator_common.h"

#include <algorithm>
#include <stdexcept>

namespace gnss_imu::internal {

namespace {

Mat3 Skew(const Vec3& vector) {
  Mat3 result = Mat3::Zero();
  result(0, 1) = -vector.z();
  result(0, 2) = vector.y();
  result(1, 0) = vector.z();
  result(1, 2) = -vector.x();
  result(2, 0) = -vector.y();
  result(2, 1) = vector.x();
  return result;
}

Vec3 LogQuaternion(const Eigen::Quaterniond& q_in) {
  const Eigen::Quaterniond q = NormalizeQuaternion(q_in);
  const Vec3 imag = q.vec();
  const double imag_norm = imag.norm();
  if (imag_norm < 1e-12) {
    return 2.0 * imag;
  }

  const double angle = 2.0 * std::atan2(imag_norm, q.w());
  return imag / imag_norm * angle;
}

NavState InitializeStateFromGnss(const GnssData& gnss) {
  NavState state;
  state.timestamp = gnss.timestamp;
  state.position = gnss.position;
  state.velocity = gnss.has_velocity ? gnss.velocity : Vec3::Zero();
  if (gnss.has_velocity) {
    const Eigen::Vector2d horizontal_velocity(gnss.velocity.x(), gnss.velocity.y());
    if (horizontal_velocity.norm() > 1e-3) {
      const double yaw = std::atan2(horizontal_velocity.y(), horizontal_velocity.x());
      state.orientation =
          NormalizeQuaternion(Eigen::Quaterniond(Eigen::AngleAxisd(yaw, Vec3::UnitZ())));
    }
  }
  return state;
}

std::optional<std::pair<Vec3, Vec3>> AverageImuSamples(
    const AlignedDeque<ImuData>& imu_samples) {
  if (imu_samples.empty()) {
    return std::nullopt;
  }

  Vec3 accel_sum = Vec3::Zero();
  Vec3 gyro_sum = Vec3::Zero();
  std::size_t count = 0;
  for (const auto& imu : imu_samples) {
    if (!imu.accel.allFinite() || !imu.gyro.allFinite()) {
      continue;
    }
    accel_sum += imu.accel;
    gyro_sum += imu.gyro;
    ++count;
  }

  if (count == 0) {
    return std::nullopt;
  }
  return std::make_pair(accel_sum / static_cast<double>(count),
                        gyro_sum / static_cast<double>(count));
}

}  // namespace

Vec3 Gravity(const FusionConfig& config) {
  return Vec3(0.0, 0.0, -config.gravity_mps2);
}

std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

Eigen::Quaterniond NormalizeQuaternion(const Eigen::Quaterniond& q_in) {
  Eigen::Quaterniond q = q_in.normalized();
  if (q.w() < 0.0) {
    q.coeffs() *= -1.0;
  }
  return q;
}

Eigen::Quaterniond ExpQuaternion(const Vec3& delta_theta) {
  const double theta = delta_theta.norm();
  if (theta < 1e-12) {
    Eigen::Quaterniond q(1.0,
                         0.5 * delta_theta.x(),
                         0.5 * delta_theta.y(),
                         0.5 * delta_theta.z());
    return NormalizeQuaternion(q);
  }

  const Vec3 axis = delta_theta / theta;
  const double half_theta = 0.5 * theta;
  return NormalizeQuaternion(Eigen::Quaterniond(
      std::cos(half_theta),
      axis.x() * std::sin(half_theta),
      axis.y() * std::sin(half_theta),
      axis.z() * std::sin(half_theta)));
}

Vec15 StateMinus(const NavState& lhs, const NavState& rhs) {
  Vec15 delta = Vec15::Zero();
  delta.segment<3>(kPos) = lhs.position - rhs.position;
  delta.segment<3>(kVel) = lhs.velocity - rhs.velocity;
  delta.segment<3>(kTheta) =
      LogQuaternion(rhs.orientation.conjugate() * lhs.orientation);
  delta.segment<3>(kAccelBias) = lhs.accel_bias - rhs.accel_bias;
  delta.segment<3>(kGyroBias) = lhs.gyro_bias - rhs.gyro_bias;
  return delta;
}

NavState ApplyError(const NavState& state, const Vec15& delta) {
  NavState updated = state;
  updated.position += delta.segment<3>(kPos);
  updated.velocity += delta.segment<3>(kVel);
  updated.orientation =
      NormalizeQuaternion(state.orientation * ExpQuaternion(delta.segment<3>(kTheta)));
  updated.accel_bias += delta.segment<3>(kAccelBias);
  updated.gyro_bias += delta.segment<3>(kGyroBias);
  return updated;
}

Mat15 Symmetrize(const Mat15& matrix) {
  return 0.5 * (matrix + matrix.transpose());
}

Eigen::MatrixXd SymmetrizeDynamic(const Eigen::MatrixXd& matrix) {
  return 0.5 * (matrix + matrix.transpose());
}

Mat15 StateTransitionMatrix(const Mat15& f, double dt) {
  const Mat15 fdt = f * dt;
  return Mat15::Identity() + fdt + 0.5 * fdt * fdt;
}

Mat15 ProcessNoiseCovariance(double dt, const FusionConfig& config) {
  Mat15 q = Mat15::Zero();

  const double dt2 = dt * dt;
  const double dt3 = dt2 * dt;
  const double dt4 = dt2 * dt2;

  // Discrete-time zero-order-hold IMU noise model:
  //   delta_v = -R * n_a * dt
  //   delta_p = -0.5 * R * n_a * dt^2
  //   delta_theta = -n_g * dt
  // The public configuration fields keep legacy "density" names, but the
  // current implementation interprets them as per-sample standard deviations.
  const double accel_var = config.accel_noise_density * config.accel_noise_density;
  const double gyro_var = config.gyro_noise_density * config.gyro_noise_density;
  const double accel_bias_var =
      config.accel_bias_random_walk * config.accel_bias_random_walk;
  const double gyro_bias_var =
      config.gyro_bias_random_walk * config.gyro_bias_random_walk;

  q.block<3, 3>(kPos, kPos) =
      Mat3::Identity() * std::max(0.25 * accel_var * dt4, 1e-10);
  q.block<3, 3>(kPos, kVel) =
      Mat3::Identity() * std::max(0.5 * accel_var * dt3, 1e-10);
  q.block<3, 3>(kVel, kPos) = q.block<3, 3>(kPos, kVel).transpose();
  q.block<3, 3>(kVel, kVel) =
      Mat3::Identity() * std::max(accel_var * dt2, 1e-10);
  q.block<3, 3>(kTheta, kTheta) =
      Mat3::Identity() * std::max(gyro_var * dt2, 1e-10);
  q.block<3, 3>(kAccelBias, kAccelBias) =
      Mat3::Identity() * std::max(accel_bias_var * dt, 1e-12);
  q.block<3, 3>(kGyroBias, kGyroBias) =
      Mat3::Identity() * std::max(gyro_bias_var * dt, 1e-12);
  return q;
}

bool IsFinite(const NavState& state) {
  return std::isfinite(state.timestamp) &&
         state.position.allFinite() &&
         state.velocity.allFinite() &&
         state.orientation.coeffs().allFinite() &&
         state.accel_bias.allFinite() &&
         state.gyro_bias.allFinite();
}

NavState InitializeStateFromGroup(const MeasureGroup& group,
                                  const FusionConfig& config) {
  NavState state = InitializeStateFromGnss(group.gnss);
  const auto imu_average = AverageImuSamples(group.imu_samples);
  if (!imu_average.has_value()) {
    return state;
  }

  const Vec3 avg_accel = imu_average->first;
  const Vec3 avg_gyro = imu_average->second;
  const double accel_norm = avg_accel.norm();
  if (std::isfinite(accel_norm) &&
      accel_norm > 1.0 &&
      std::abs(accel_norm - config.gravity_mps2) < 4.0) {
    const Eigen::Quaterniond roll_pitch =
        NormalizeQuaternion(Eigen::Quaterniond::FromTwoVectors(
            avg_accel.normalized(), Vec3::UnitZ()));
    state.orientation = roll_pitch;

    if (group.gnss.has_velocity) {
      const Eigen::Vector2d horizontal_velocity(
          group.gnss.velocity.x(), group.gnss.velocity.y());
      if (horizontal_velocity.norm() > 0.5) {
        const double yaw =
            std::atan2(horizontal_velocity.y(), horizontal_velocity.x());
        state.orientation = NormalizeQuaternion(
            Eigen::Quaterniond(Eigen::AngleAxisd(yaw, Vec3::UnitZ())) *
            roll_pitch);
      }
    }
  }

  if (avg_gyro.allFinite() &&
      avg_gyro.norm() <= config.initial_gyro_bias_max_norm) {
    state.gyro_bias = avg_gyro;
  }
  return state;
}

void PropagateNominal(NavState* state,
                      const ImuData& imu_prev,
                      const ImuData& imu_curr,
                      const FusionConfig& config) {
  const double dt = imu_curr.timestamp - imu_prev.timestamp;
  if (dt <= 0.0 || !std::isfinite(dt)) {
    return;
  }

  const Vec3 gyro_prev = imu_prev.gyro - state->gyro_bias;
  const Vec3 gyro_curr = imu_curr.gyro - state->gyro_bias;
  const Vec3 omega_mid = 0.5 * (gyro_prev + gyro_curr);

  const Vec3 accel_prev = imu_prev.accel - state->accel_bias;
  const Vec3 accel_curr = imu_curr.accel - state->accel_bias;
  const Vec3 accel_body_mid = 0.5 * (accel_prev + accel_curr);

  const Eigen::Quaterniond q_mid =
      NormalizeQuaternion(state->orientation * ExpQuaternion(0.5 * omega_mid * dt));
  const Vec3 accel_world = q_mid * accel_body_mid + Gravity(config);

  state->position += state->velocity * dt + 0.5 * accel_world * dt * dt;
  state->velocity += accel_world * dt;
  state->orientation =
      NormalizeQuaternion(state->orientation * ExpQuaternion(omega_mid * dt));
  state->timestamp = imu_curr.timestamp;
}

Mat15 ContinuousTimeJacobian(const NavState& state,
                             const ImuData& imu_prev,
                             const ImuData& imu_curr) {
  const double dt = imu_curr.timestamp - imu_prev.timestamp;
  const Vec3 accel_body =
      0.5 * ((imu_prev.accel - state.accel_bias) + (imu_curr.accel - state.accel_bias));
  const Vec3 gyro_body =
      0.5 * ((imu_prev.gyro - state.gyro_bias) + (imu_curr.gyro - state.gyro_bias));
  const Eigen::Quaterniond q_mid =
      dt > 0.0
          ? NormalizeQuaternion(state.orientation * ExpQuaternion(0.5 * gyro_body * dt))
          : state.orientation;
  const Mat3 rotation = q_mid.toRotationMatrix();

  Mat15 f = Mat15::Zero();
  f.block<3, 3>(kPos, kVel) = Mat3::Identity();
  f.block<3, 3>(kVel, kTheta) = -rotation * Skew(accel_body);
  f.block<3, 3>(kVel, kAccelBias) = -rotation;
  f.block<3, 3>(kTheta, kTheta) = -Skew(gyro_body);
  f.block<3, 3>(kTheta, kGyroBias) = -Mat3::Identity();
  return f;
}

Mat15 UpdateCovarianceEskf(const NavState& state,
                           const ImuData& imu_prev,
                           const ImuData& imu_curr,
                           const FusionConfig& config,
                           const Mat15& covariance) {
  const double dt = imu_curr.timestamp - imu_prev.timestamp;
  if (dt <= 0.0 || !std::isfinite(dt)) {
    return covariance;
  }

  const Mat15 f = ContinuousTimeJacobian(state, imu_prev, imu_curr);
  const Mat15 phi = StateTransitionMatrix(f, dt);
  Mat15 predicted = phi * covariance * phi.transpose() + ProcessNoiseCovariance(dt, config);
  predicted = Symmetrize(predicted);
  RegularizeInPlace(&predicted, 1e-9);
  return predicted;
}

UnscentedWeights MakeUnscentedWeights(const FusionConfig& config) {
  const int n = kErrorStateDim;
  const double lambda =
      config.ukf_alpha * config.ukf_alpha * (n + config.ukf_kappa) - n;
  if (n + lambda <= 1e-9) {
    throw std::invalid_argument("UKF scaling parameters produced a non-positive n + lambda");
  }
  const double scale = std::sqrt(n + lambda);

  UnscentedWeights weights;
  weights.mean = Eigen::VectorXd::Constant(2 * n + 1, 0.5 / (n + lambda));
  weights.covariance = weights.mean;
  weights.mean(0) = lambda / (n + lambda);
  weights.covariance(0) =
      lambda / (n + lambda) + (1.0 - config.ukf_alpha * config.ukf_alpha + config.ukf_beta);
  weights.scale = scale;
  return weights;
}

Mat15 RobustCholeskyInput(const Mat15& covariance) {
  Mat15 adjusted = Symmetrize(covariance);
  for (int attempt = 0; attempt < 6; ++attempt) {
    Eigen::LLT<Mat15> llt(adjusted);
    if (llt.info() == Eigen::Success) {
      return adjusted;
    }
    RegularizeInPlace(&adjusted, std::pow(10.0, static_cast<double>(attempt - 9)));
  }
  RegularizeInPlace(&adjusted, 1e-3);
  return adjusted;
}

Eigen::Quaterniond WeightedQuaternionAverage(
    const AlignedVector<NavState>& sigma_states,
    const Eigen::VectorXd& weights) {
  Eigen::Quaterniond mean = sigma_states.front().orientation;
  for (int iteration = 0; iteration < 8; ++iteration) {
    Vec3 accumulated = Vec3::Zero();
    for (std::size_t i = 0; i < sigma_states.size(); ++i) {
      accumulated +=
          weights(static_cast<Eigen::Index>(i)) *
          LogQuaternion(mean.conjugate() * sigma_states[i].orientation);
    }
    if (accumulated.norm() < 1e-10) {
      break;
    }
    mean = NormalizeQuaternion(mean * ExpQuaternion(accumulated));
  }
  return mean;
}

std::optional<ImuData> PropagateOverEdge(NavState* state,
                                         const std::optional<ImuData>& last_imu,
                                         const ImuEdge& edge,
                                         const FusionConfig& config) {
  if (!last_imu.has_value()) {
    return std::nullopt;
  }

  ImuData previous = edge.start_imu;
  if (previous.timestamp + kEps < last_imu->timestamp) {
    previous = *last_imu;
  }

  for (const auto& imu : edge.imu_samples) {
    if (imu.timestamp <= previous.timestamp + kEps) {
      previous = imu;
      continue;
    }
    PropagateNominal(state, previous, imu, config);
    previous = imu;
  }

  if (previous.timestamp + kEps < edge.end_time) {
    ImuData hold = previous;
    hold.timestamp = edge.end_time;
    PropagateNominal(state, previous, hold, config);
    previous = hold;
  }

  return previous;
}

std::optional<ImuData> PropagateOverGroup(NavState* state,
                                          const std::optional<ImuData>& last_imu,
                                          const MeasureGroup& group,
                                          const FusionConfig& config) {
  if (!last_imu.has_value()) {
    return std::nullopt;
  }

  ImuEdge edge;
  edge.start_imu = *last_imu;
  edge.imu_samples = group.imu_samples;
  edge.end_time = group.gnss.timestamp;
  return PropagateOverEdge(state, last_imu, edge, config);
}

EdgeLinearization PropagateEdgeLinearized(const NavState& start_state,
                                          const ImuEdge& edge,
                                          const FusionConfig& config) {
  EdgeLinearization result;
  result.predicted = start_state;
  result.predicted.timestamp = start_state.timestamp;

  ImuData previous = edge.start_imu;
  auto integrate_step = [&](const ImuData& imu) {
    const double dt = imu.timestamp - previous.timestamp;
    if (dt <= 0.0 || !std::isfinite(dt)) {
      previous = imu;
      return;
    }

    const Mat15 f = ContinuousTimeJacobian(result.predicted, previous, imu);
    const Mat15 phi = StateTransitionMatrix(f, dt);
    result.transition = phi * result.transition;
    result.covariance =
        phi * result.covariance * phi.transpose() + ProcessNoiseCovariance(dt, config);
    result.covariance = Symmetrize(result.covariance);
    PropagateNominal(&result.predicted, previous, imu, config);
    previous = imu;
  };

  for (const auto& imu : edge.imu_samples) {
    if (imu.timestamp <= previous.timestamp + kEps) {
      previous = imu;
      continue;
    }
    integrate_step(imu);
  }

  if (previous.timestamp + kEps < edge.end_time) {
    ImuData hold = previous;
    hold.timestamp = edge.end_time;
    integrate_step(hold);
  }

  result.predicted.timestamp = edge.end_time;
  result.final_imu = previous;
  RegularizeInPlace(&result.covariance, 1e-9);
  return result;
}

Mat3 SafeCovarianceInverse(const Mat3& covariance, double default_std) {
  Mat3 adjusted = covariance;
  if (!adjusted.allFinite()) {
    adjusted = Mat3::Identity() * default_std * default_std;
  }
  for (int i = 0; i < 3; ++i) {
    adjusted(i, i) = std::max(adjusted(i, i), 1e-6);
  }
  return adjusted.inverse();
}

Mat15 SafeInformationInverse(const Mat15& covariance, double epsilon) {
  return SafeInverseSpd(Symmetrize(covariance), epsilon);
}

}  // namespace gnss_imu::internal
