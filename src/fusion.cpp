#include "gnss_imu_fusion/fusion.h"

#include "estimators/eskf_estimator.h"
#include "estimators/fgo_estimator.h"
#include "estimators/ukf_estimator.h"
#include "internal/estimator_common.h"

#include <stdexcept>

namespace gnss_imu {

std::unique_ptr<FusionEstimator> CreateEstimator(const std::string& algorithm,
                                                 const FusionConfig& config) {
  const std::string normalized = internal::ToLower(algorithm);
  if (normalized == "eskf") {
    return std::make_unique<internal::EskfEstimator>(config);
  }
  if (normalized == "ukf") {
    return std::make_unique<internal::UkfEstimator>(config);
  }
  if (normalized == "fgo") {
    return std::make_unique<internal::FgoEstimator>(config);
  }
  throw std::invalid_argument("Unsupported algorithm: " + algorithm);
}

}  // namespace gnss_imu
