#pragma once

#include "estimators/estimator_base.h"

#include <optional>
#include <string>

namespace gnss_imu::internal {

struct GraphNode {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  NavState state;
  GnssData gnss;
};

struct PriorConstraint {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  bool valid = false;
  NavState state;
  Mat15 information = Mat15::Identity();
};

struct Linearization {
  Eigen::MatrixXd hessian;
  Eigen::VectorXd gradient;
  double cost = 0.0;
};

class FgoEstimator final : public EstimatorBase {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  explicit FgoEstimator(const FusionConfig& config);

  std::string Name() const override;
  void Process(const MeasureGroup& group) override;

 private:
  static double HuberScale(double squared_norm, double delta);

  void AddPriorFactor(Linearization* linearization) const;
  void AddGnssFactors(Linearization* linearization) const;
  void AddImuFactors(Linearization* linearization) const;
  Linearization BuildLinearization() const;
  void OptimizeWindow();
  void MarginalizeOldestNode();

  AlignedDeque<GraphNode> nodes_;
  AlignedDeque<ImuEdge> imu_edges_;
  PriorConstraint prior_;
  std::optional<ImuData> last_imu_;
};

}  // namespace gnss_imu::internal
