#pragma once

#include "internal/estimator_common.h"

#include <utility>

namespace gnss_imu::internal {

class EstimatorBase : public FusionEstimator {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  explicit EstimatorBase(FusionConfig config) : config_(std::move(config)) {}

  const NavState& CurrentState() const override { return current_state_; }
  const AlignedVector<TrajectorySample>& Trajectory() const override {
    return trajectory_;
  }
  bool IsInitialized() const override { return initialized_; }

 protected:
  void InitializeFromGroup(const MeasureGroup& group) {
    current_state_ = InitializeStateFromGroup(group, config_);
    initialized_ = true;
    AppendTrajectory();
  }

  void AppendTrajectory() {
    trajectory_.push_back(TrajectorySample{current_state_.timestamp, current_state_});
  }

  FusionConfig config_;
  NavState current_state_;
  AlignedVector<TrajectorySample> trajectory_;
  bool initialized_ = false;
};

}  // namespace gnss_imu::internal
