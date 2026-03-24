#include "estimators/ukf_estimator.h"

#include <stdexcept>

namespace gnss_imu::internal {

UkfEstimator::UkfEstimator(const FusionConfig& config)
    : EstimatorBase(config), covariance_(Mat15::Zero()), weights_(MakeUnscentedWeights(config)) {
  covariance_.block<3, 3>(kPos, kPos) =
      Mat3::Identity() * config_.initial_position_std * config_.initial_position_std;
  covariance_.block<3, 3>(kVel, kVel) =
      Mat3::Identity() * config_.initial_velocity_std * config_.initial_velocity_std;
  covariance_.block<3, 3>(kTheta, kTheta) =
      Mat3::Identity() * config_.initial_attitude_std * config_.initial_attitude_std;
  covariance_.block<3, 3>(kAccelBias, kAccelBias) =
      Mat3::Identity() * config_.initial_bias_std * config_.initial_bias_std;
  covariance_.block<3, 3>(kGyroBias, kGyroBias) =
      Mat3::Identity() * config_.initial_bias_std * config_.initial_bias_std;
}

std::string UkfEstimator::Name() const { return "ukf"; }

void UkfEstimator::Process(const MeasureGroup& group) {
  if (!initialized_) {
    InitializeFromGroup(group);
    last_imu_ = group.imu_samples.empty() ? std::nullopt
                                          : std::optional<ImuData>(group.imu_samples.back());
    return;
  }

  const AlignedVector<NavState> predicted_sigma = PredictSigmaStates(group);
  NavState predicted_mean = current_state_;
  predicted_mean.position.setZero();
  predicted_mean.velocity.setZero();
  predicted_mean.accel_bias.setZero();
  predicted_mean.gyro_bias.setZero();
  for (std::size_t i = 0; i < predicted_sigma.size(); ++i) {
    const double w = weights_.mean(static_cast<Eigen::Index>(i));
    predicted_mean.position += w * predicted_sigma[i].position;
    predicted_mean.velocity += w * predicted_sigma[i].velocity;
    predicted_mean.accel_bias += w * predicted_sigma[i].accel_bias;
    predicted_mean.gyro_bias += w * predicted_sigma[i].gyro_bias;
  }
  predicted_mean.orientation = WeightedQuaternionAverage(predicted_sigma, weights_.mean);
  predicted_mean.timestamp = group.gnss.timestamp;

  Mat15 predicted_covariance = Mat15::Zero();
  for (std::size_t i = 0; i < predicted_sigma.size(); ++i) {
    const Vec15 error = StateMinus(predicted_sigma[i], predicted_mean);
    predicted_covariance +=
        weights_.covariance(static_cast<Eigen::Index>(i)) * error * error.transpose();
  }
  const double total_dt = std::max(group.gnss.timestamp - current_state_.timestamp, 0.0);
  Mat15 additive_process_noise = ProcessNoiseCovariance(total_dt, config_);
  if (last_imu_.has_value()) {
    ImuEdge edge;
    edge.start_imu = *last_imu_;
    edge.imu_samples = group.imu_samples;
    edge.end_time = group.gnss.timestamp;
    additive_process_noise =
        PropagateEdgeLinearized(current_state_, edge, config_).covariance;
  }
  predicted_covariance += additive_process_noise;
  predicted_covariance = Symmetrize(predicted_covariance);
  RegularizeInPlace(&predicted_covariance, 1e-9);

  Vec6 z_mean = Vec6::Zero();
  std::vector<Vec6> measurement_sigma;
  measurement_sigma.reserve(predicted_sigma.size());
  for (std::size_t i = 0; i < predicted_sigma.size(); ++i) {
    Vec6 z = Vec6::Zero();
    z.segment<3>(0) = predicted_sigma[i].position;
    z.segment<3>(3) = predicted_sigma[i].velocity;
    measurement_sigma.push_back(z);
    z_mean += weights_.mean(static_cast<Eigen::Index>(i)) * z;
  }

  Mat6 s = Mat6::Zero();
  Eigen::Matrix<double, 15, 6> pxz = Eigen::Matrix<double, 15, 6>::Zero();
  for (std::size_t i = 0; i < predicted_sigma.size(); ++i) {
    const Vec15 state_error = StateMinus(predicted_sigma[i], predicted_mean);
    const Vec6 measurement_error = measurement_sigma[i] - z_mean;
    s += weights_.covariance(static_cast<Eigen::Index>(i)) *
         measurement_error * measurement_error.transpose();
    pxz += weights_.covariance(static_cast<Eigen::Index>(i)) *
           state_error * measurement_error.transpose();
  }

  Mat6 measurement_covariance = Mat6::Zero();
  measurement_covariance.block<3, 3>(0, 0) =
      group.gnss.position_covariance.allFinite()
          ? group.gnss.position_covariance
          : Mat3::Identity() * config_.gnss_position_std * config_.gnss_position_std;
  measurement_covariance.block<3, 3>(3, 3) =
      group.gnss.has_velocity
          ? (group.gnss.velocity_covariance.allFinite()
                 ? group.gnss.velocity_covariance
                 : Mat3::Identity() * config_.gnss_velocity_std * config_.gnss_velocity_std)
          : Mat3::Identity() * 1e6;
  s += measurement_covariance;

  const Mat6 s_inverse = SafeInverseSpd(s, 1e-9);
  const Eigen::Matrix<double, 15, 6> gain = pxz * s_inverse;

  Vec6 observation = Vec6::Zero();
  observation.segment<3>(0) = group.gnss.position;
  observation.segment<3>(3) =
      group.gnss.has_velocity ? group.gnss.velocity : predicted_mean.velocity;

  const Vec15 correction = gain * (observation - z_mean);
  current_state_ = ApplyError(predicted_mean, correction);
  current_state_.timestamp = group.gnss.timestamp;
  covariance_ = Symmetrize(predicted_covariance - gain * s * gain.transpose());
  RegularizeInPlace(&covariance_, 1e-9);

  if (!group.imu_samples.empty()) {
    last_imu_ = group.imu_samples.back();
  } else if (last_imu_.has_value() && last_imu_->timestamp + kEps < group.gnss.timestamp) {
    last_imu_->timestamp = group.gnss.timestamp;
  }

  if (!IsFinite(current_state_)) {
    throw std::runtime_error("UKF produced a non-finite state");
  }
  AppendTrajectory();
}

AlignedVector<NavState> UkfEstimator::PredictSigmaStates(const MeasureGroup& group) {
  const Mat15 cholesky_input = RobustCholeskyInput(covariance_);
  const Eigen::LLT<Mat15> llt(cholesky_input);
  const Mat15 square_root = weights_.scale * llt.matrixL().toDenseMatrix();

  AlignedVector<NavState> sigma_states(2 * kErrorStateDim + 1, current_state_);
  for (int column = 0; column < kErrorStateDim; ++column) {
    const Vec15 delta = square_root.col(column);
    sigma_states[static_cast<std::size_t>(column + 1)] = ApplyError(current_state_, delta);
    sigma_states[static_cast<std::size_t>(column + 1 + kErrorStateDim)] =
        ApplyError(current_state_, -delta);
  }

  for (auto& sigma_state : sigma_states) {
    const auto final_imu = PropagateOverGroup(&sigma_state, last_imu_, group, config_);
    if (final_imu.has_value()) {
      sigma_state.timestamp = final_imu->timestamp;
    } else {
      sigma_state.timestamp = group.gnss.timestamp;
    }
  }

  return sigma_states;
}

}  // namespace gnss_imu::internal
