#include "estimators/eskf_estimator.h"

#include <stdexcept>

namespace gnss_imu::internal {

EskfEstimator::EskfEstimator(const FusionConfig& config)
    : EstimatorBase(config), covariance_(Mat15::Zero()) {}

std::string EskfEstimator::Name() const { return "eskf"; }

void EskfEstimator::Process(const MeasureGroup& group) {
  if (!initialized_) {
    InitializeFromGroup(group);
    last_imu_ = group.imu_samples.empty() ? std::nullopt
                                          : std::optional<ImuData>(group.imu_samples.back());
    covariance_.setZero();
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
    return;
  }

  if (last_imu_.has_value()) {
    ImuData previous = *last_imu_;
    for (const auto& imu : group.imu_samples) {
      if (imu.timestamp <= previous.timestamp + kEps) {
        previous = imu;
        continue;
      }
      covariance_ = UpdateCovarianceEskf(current_state_, previous, imu, config_, covariance_);
      PropagateNominal(&current_state_, previous, imu, config_);
      previous = imu;
    }
    if (previous.timestamp + kEps < group.gnss.timestamp) {
      ImuData hold = previous;
      hold.timestamp = group.gnss.timestamp;
      covariance_ = UpdateCovarianceEskf(current_state_, previous, hold, config_, covariance_);
      PropagateNominal(&current_state_, previous, hold, config_);
      previous = hold;
    }
    last_imu_ = previous;
  } else if (!group.imu_samples.empty()) {
    last_imu_ = group.imu_samples.back();
  }

  Eigen::Matrix<double, 6, 15> h = Eigen::Matrix<double, 6, 15>::Zero();
  h.block<3, 3>(0, kPos) = Mat3::Identity();
  h.block<3, 3>(3, kVel) = Mat3::Identity();

  Vec6 innovation = Vec6::Zero();
  innovation.segment<3>(0) = group.gnss.position - current_state_.position;
  if (group.gnss.has_velocity) {
    innovation.segment<3>(3) = group.gnss.velocity - current_state_.velocity;
  }

  Mat6 measurement_covariance = Mat6::Identity();
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

  const Mat6 s = h * covariance_ * h.transpose() + measurement_covariance;
  const Mat6 s_inverse = SafeInverseSpd(s, 1e-9);
  const Eigen::Matrix<double, 15, 6> k =
      covariance_ * h.transpose() * s_inverse;
  const Vec15 delta = k * innovation;

  current_state_ = ApplyError(current_state_, delta);
  current_state_.timestamp = group.gnss.timestamp;
  covariance_ = Symmetrize(
      (Mat15::Identity() - k * h) * covariance_ * (Mat15::Identity() - k * h).transpose() +
      k * measurement_covariance * k.transpose());
  RegularizeInPlace(&covariance_, 1e-9);

  if (!IsFinite(current_state_)) {
    throw std::runtime_error("ESKF produced a non-finite state");
  }
  AppendTrajectory();
}

}  // namespace gnss_imu::internal
