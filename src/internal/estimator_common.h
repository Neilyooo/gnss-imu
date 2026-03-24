#pragma once

#include "gnss_imu_fusion/fusion.h"

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <cmath>
#include <deque>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace gnss_imu::internal {

constexpr int kErrorStateDim = 15;
constexpr double kEps = 1e-9;

using Mat6 = Eigen::Matrix<double, 6, 6>;
using Vec6 = Eigen::Matrix<double, 6, 1>;

enum ErrorStateIndex {
  kPos = 0,
  kVel = 3,
  kTheta = 6,
  kAccelBias = 9,
  kGyroBias = 12
};

template <typename MatrixType>
void RegularizeInPlace(MatrixType* matrix, double epsilon = 1e-9) {
  matrix->diagonal().array() += epsilon;
}

template <typename MatrixType>
MatrixType SafeInverseSpd(MatrixType matrix, double epsilon = 1e-9) {
  matrix = 0.5 * (matrix + matrix.transpose());
  for (int attempt = 0; attempt < 6; ++attempt) {
    RegularizeInPlace(&matrix, epsilon * std::pow(10.0, static_cast<double>(attempt)));
    Eigen::LDLT<MatrixType> ldlt(matrix);
    if (ldlt.info() != Eigen::Success) {
      continue;
    }
    MatrixType inverse = ldlt.solve(MatrixType::Identity());
    if (inverse.allFinite()) {
      return inverse;
    }
  }
  return MatrixType::Identity();
}

template <typename MatrixType>
MatrixType SafeSqrtInformation(MatrixType information, double epsilon = 1e-9) {
  information = 0.5 * (information + information.transpose());
  for (int attempt = 0; attempt < 6; ++attempt) {
    RegularizeInPlace(&information, epsilon * std::pow(10.0, static_cast<double>(attempt)));
    Eigen::LLT<MatrixType> llt(information);
    if (llt.info() == Eigen::Success) {
      return llt.matrixL().transpose();
    }
  }
  return MatrixType::Identity();
}

Vec3 Gravity(const FusionConfig& config);
std::string ToLower(std::string value);
Eigen::Quaterniond NormalizeQuaternion(const Eigen::Quaterniond& q_in);
Eigen::Quaterniond ExpQuaternion(const Vec3& delta_theta);
Vec15 StateMinus(const NavState& lhs, const NavState& rhs);
NavState ApplyError(const NavState& state, const Vec15& delta);
Mat15 Symmetrize(const Mat15& matrix);
Eigen::MatrixXd SymmetrizeDynamic(const Eigen::MatrixXd& matrix);
Mat15 StateTransitionMatrix(const Mat15& f, double dt);
Mat15 ProcessNoiseCovariance(double dt, const FusionConfig& config);
bool IsFinite(const NavState& state);
NavState InitializeStateFromGroup(const MeasureGroup& group,
                                  const FusionConfig& config);
void PropagateNominal(NavState* state,
                      const ImuData& imu_prev,
                      const ImuData& imu_curr,
                      const FusionConfig& config);
Mat15 ContinuousTimeJacobian(const NavState& state,
                             const ImuData& imu_prev,
                             const ImuData& imu_curr);
Mat15 UpdateCovarianceEskf(const NavState& state,
                           const ImuData& imu_prev,
                           const ImuData& imu_curr,
                           const FusionConfig& config,
                           const Mat15& covariance);

struct UnscentedWeights {
  Eigen::VectorXd mean;
  Eigen::VectorXd covariance;
  double scale = 1.0;
};

UnscentedWeights MakeUnscentedWeights(const FusionConfig& config);
Mat15 RobustCholeskyInput(const Mat15& covariance);
Eigen::Quaterniond WeightedQuaternionAverage(
    const AlignedVector<NavState>& sigma_states,
    const Eigen::VectorXd& weights);

struct ImuEdge {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ImuData start_imu;
  AlignedDeque<ImuData> imu_samples;
  double end_time = 0.0;
};

std::optional<ImuData> PropagateOverEdge(NavState* state,
                                         const std::optional<ImuData>& last_imu,
                                         const ImuEdge& edge,
                                         const FusionConfig& config);
std::optional<ImuData> PropagateOverGroup(NavState* state,
                                          const std::optional<ImuData>& last_imu,
                                          const MeasureGroup& group,
                                          const FusionConfig& config);

struct EdgeLinearization {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  NavState predicted;
  Mat15 transition = Mat15::Identity();
  Mat15 covariance = Mat15::Zero();
  std::optional<ImuData> final_imu;
};

EdgeLinearization PropagateEdgeLinearized(const NavState& start_state,
                                          const ImuEdge& edge,
                                          const FusionConfig& config);

Mat3 SafeCovarianceInverse(const Mat3& covariance, double default_std);
Mat15 SafeInformationInverse(const Mat15& covariance, double epsilon = 1e-9);

}  // namespace gnss_imu::internal
