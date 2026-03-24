#pragma once

#include "estimators/estimator_base.h"

#include <optional>
#include <string>

namespace gnss_imu::internal {

class UkfEstimator final : public EstimatorBase {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  explicit UkfEstimator(const FusionConfig& config);

  std::string Name() const override;
  void Process(const MeasureGroup& group) override;

 private:
  AlignedVector<NavState> PredictSigmaStates(const MeasureGroup& group);

  Mat15 covariance_ = Mat15::Zero();
  UnscentedWeights weights_;
  std::optional<ImuData> last_imu_;
};

}  // namespace gnss_imu::internal
