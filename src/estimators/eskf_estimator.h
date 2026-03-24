#pragma once

#include "estimators/estimator_base.h"

#include <optional>
#include <string>

namespace gnss_imu::internal {

class EskfEstimator final : public EstimatorBase {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  explicit EskfEstimator(const FusionConfig& config);

  std::string Name() const override;
  void Process(const MeasureGroup& group) override;

 private:
  Mat15 covariance_ = Mat15::Zero();
  std::optional<ImuData> last_imu_;
};

}  // namespace gnss_imu::internal
