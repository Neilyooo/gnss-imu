#include "estimators/fgo_estimator.h"

#include <stdexcept>

namespace gnss_imu::internal {

FgoEstimator::FgoEstimator(const FusionConfig& config) : EstimatorBase(config) {}

std::string FgoEstimator::Name() const { return "fgo"; }

void FgoEstimator::Process(const MeasureGroup& group) {
  if (!initialized_) {
    InitializeFromGroup(group);
    GraphNode first_node;
    first_node.state = current_state_;
    first_node.gnss = group.gnss;
    nodes_.push_back(first_node);
    prior_.valid = true;
    prior_.state = current_state_;
    prior_.information = Mat15::Zero();
    prior_.information.block<3, 3>(kPos, kPos) =
        Mat3::Identity() /
        (config_.gnss_position_std * config_.gnss_position_std);
    prior_.information.block<3, 3>(kVel, kVel) =
        Mat3::Identity() /
        (config_.gnss_velocity_std * config_.gnss_velocity_std);
    prior_.information.block<3, 3>(kTheta, kTheta) =
        Mat3::Identity() /
        (config_.initial_attitude_std * config_.initial_attitude_std);
    prior_.information.block<3, 3>(kAccelBias, kAccelBias) =
        Mat3::Identity() /
        (config_.initial_bias_std * config_.initial_bias_std);
    prior_.information.block<3, 3>(kGyroBias, kGyroBias) =
        Mat3::Identity() /
        (config_.initial_bias_std * config_.initial_bias_std);
    last_imu_ = group.imu_samples.empty() ? std::nullopt
                                          : std::optional<ImuData>(group.imu_samples.back());
    return;
  }

  GraphNode node;
  node.state = nodes_.back().state;
  node.gnss = group.gnss;
  node.state.timestamp = nodes_.back().state.timestamp;

  ImuEdge edge;
  if (last_imu_.has_value()) {
    edge.start_imu = *last_imu_;
  } else if (!group.imu_samples.empty()) {
    edge.start_imu = group.imu_samples.front();
  } else {
    edge.start_imu.timestamp = nodes_.back().state.timestamp;
  }
  edge.imu_samples = group.imu_samples;
  edge.end_time = group.gnss.timestamp;

  const auto final_imu = PropagateOverEdge(&node.state, last_imu_, edge, config_);
  node.state.timestamp = group.gnss.timestamp;
  nodes_.push_back(node);
  imu_edges_.push_back(edge);

  const auto nodes_before_opt = nodes_;
  const auto imu_edges_before_opt = imu_edges_;
  const PriorConstraint prior_before_opt = prior_;
  const NavState predicted_tail_before_opt = nodes_.back().state;

  OptimizeWindow();
  bool reject_optimization = false;
  if (nodes_.empty() || !IsFinite(nodes_.back().state)) {
    reject_optimization = true;
  } else {
    const double position_error =
        (nodes_.back().state.position - group.gnss.position).norm();
    const double motion_jump =
        (nodes_.back().state.position - predicted_tail_before_opt.position).norm();
    const double velocity_error =
        group.gnss.has_velocity
            ? (nodes_.back().state.velocity - group.gnss.velocity).norm()
            : 0.0;
    reject_optimization =
        position_error > 5.0 ||
        motion_jump > 5.0 ||
        (group.gnss.has_velocity && velocity_error > 5.0);
  }
  if (reject_optimization) {
    nodes_ = nodes_before_opt;
    imu_edges_ = imu_edges_before_opt;
    prior_ = prior_before_opt;
    nodes_.back().state = predicted_tail_before_opt;
    nodes_.back().state.position = group.gnss.position;
    if (group.gnss.has_velocity) {
      nodes_.back().state.velocity = group.gnss.velocity;
    }
    nodes_.back().state.timestamp = group.gnss.timestamp;
  }

  while (static_cast<int>(nodes_.size()) > config_.fgo_window_size) {
    MarginalizeOldestNode();
  }

  current_state_ = nodes_.back().state;
  current_state_.timestamp = group.gnss.timestamp;
  last_imu_ = final_imu.has_value()
                  ? final_imu
                  : (group.imu_samples.empty() ? last_imu_
                                               : std::optional<ImuData>(group.imu_samples.back()));

  if (!IsFinite(current_state_)) {
    throw std::runtime_error("FGO produced a non-finite state");
  }
  AppendTrajectory();
}

double FgoEstimator::HuberScale(double squared_norm, double delta) {
  if (!std::isfinite(squared_norm) || squared_norm <= 0.0) {
    return 1.0;
  }
  const double norm = std::sqrt(squared_norm);
  if (norm <= delta) {
    return 1.0;
  }
  return std::sqrt(delta / norm);
}

void FgoEstimator::AddPriorFactor(Linearization* linearization) const {
  if (!prior_.valid || nodes_.empty()) {
    return;
  }
  const Vec15 residual = StateMinus(nodes_.front().state, prior_.state);
  const Mat15 sqrt_information = SafeSqrtInformation(prior_.information, 1e-9);
  const Vec15 weighted_residual = sqrt_information * residual;

  linearization->cost += weighted_residual.squaredNorm();
  linearization->hessian.block<15, 15>(0, 0) +=
      sqrt_information.transpose() * sqrt_information;
  linearization->gradient.segment<15>(0) +=
      sqrt_information.transpose() * weighted_residual;
}

void FgoEstimator::AddGnssFactors(Linearization* linearization) const {
  for (std::size_t index = 0; index < nodes_.size(); ++index) {
    const auto& node = nodes_[index];
    Eigen::Matrix<double, 6, 15> jacobian = Eigen::Matrix<double, 6, 15>::Zero();
    jacobian.block<3, 3>(0, kPos) = Mat3::Identity();
    jacobian.block<3, 3>(3, kVel) = Mat3::Identity();

    Vec6 residual = Vec6::Zero();
    residual.segment<3>(0) = node.state.position - node.gnss.position;
    if (node.gnss.has_velocity) {
      residual.segment<3>(3) = node.state.velocity - node.gnss.velocity;
    }

    Mat6 information = Mat6::Zero();
    information.block<3, 3>(0, 0) =
        SafeCovarianceInverse(node.gnss.position_covariance, config_.gnss_position_std);
    information.block<3, 3>(3, 3) = node.gnss.has_velocity
                                        ? SafeCovarianceInverse(node.gnss.velocity_covariance,
                                                                config_.gnss_velocity_std)
                                        : Mat3::Identity() * 1e-6;
    const Mat6 sqrt_information = SafeSqrtInformation(information, 1e-9);
    Eigen::Matrix<double, 6, 15> weighted_jacobian =
        sqrt_information * jacobian;
    Vec6 weighted_residual = sqrt_information * residual;
    const double robust_scale =
        HuberScale(weighted_residual.squaredNorm(), 6.0);
    weighted_jacobian *= robust_scale;
    weighted_residual *= robust_scale;

    const int offset = static_cast<int>(index) * kErrorStateDim;
    linearization->cost += weighted_residual.squaredNorm();
    linearization->hessian.block(offset, offset, 15, 15) +=
        weighted_jacobian.transpose() * weighted_jacobian;
    linearization->gradient.segment(offset, 15) +=
        weighted_jacobian.transpose() * weighted_residual;
  }
}

void FgoEstimator::AddImuFactors(Linearization* linearization) const {
  for (std::size_t index = 1; index < nodes_.size(); ++index) {
    const auto& edge = imu_edges_[index - 1];
    const auto& first = nodes_[index - 1].state;
    const auto& second = nodes_[index].state;
    const EdgeLinearization edge_linearization =
        PropagateEdgeLinearized(first, edge, config_);
    const Vec15 residual = StateMinus(second, edge_linearization.predicted);
    const Mat15 jacobian_first = -edge_linearization.transition;
    const Mat15 jacobian_second = Mat15::Identity();
    const Mat15 information =
        0.1 * SafeInformationInverse(edge_linearization.covariance, 1e-8);
    const Mat15 sqrt_information = SafeSqrtInformation(information, 1e-9);

    Mat15 weighted_jacobian_first = sqrt_information * jacobian_first;
    Mat15 weighted_jacobian_second = sqrt_information * jacobian_second;
    Vec15 weighted_residual = sqrt_information * residual;
    const double robust_scale =
        HuberScale(weighted_residual.squaredNorm(), 8.0);
    weighted_jacobian_first *= robust_scale;
    weighted_jacobian_second *= robust_scale;
    weighted_residual *= robust_scale;

    const int first_offset = static_cast<int>(index - 1) * kErrorStateDim;
    const int second_offset = static_cast<int>(index) * kErrorStateDim;

    linearization->cost += weighted_residual.squaredNorm();
    linearization->hessian.block(first_offset, first_offset, 15, 15) +=
        weighted_jacobian_first.transpose() * weighted_jacobian_first;
    linearization->hessian.block(first_offset, second_offset, 15, 15) +=
        weighted_jacobian_first.transpose() * weighted_jacobian_second;
    linearization->hessian.block(second_offset, first_offset, 15, 15) +=
        weighted_jacobian_second.transpose() * weighted_jacobian_first;
    linearization->hessian.block(second_offset, second_offset, 15, 15) +=
        weighted_jacobian_second.transpose() * weighted_jacobian_second;
    linearization->gradient.segment(first_offset, 15) +=
        weighted_jacobian_first.transpose() * weighted_residual;
    linearization->gradient.segment(second_offset, 15) +=
        weighted_jacobian_second.transpose() * weighted_residual;
  }
}

Linearization FgoEstimator::BuildLinearization() const {
  const int total_dim = static_cast<int>(nodes_.size()) * kErrorStateDim;
  Linearization linearization;
  linearization.hessian = Eigen::MatrixXd::Zero(total_dim, total_dim);
  linearization.gradient = Eigen::VectorXd::Zero(total_dim);

  AddPriorFactor(&linearization);
  AddGnssFactors(&linearization);
  AddImuFactors(&linearization);

  linearization.hessian = SymmetrizeDynamic(linearization.hessian);
  return linearization;
}

void FgoEstimator::OptimizeWindow() {
  if (nodes_.size() < 2) {
    return;
  }

  for (int iteration = 0; iteration < config_.fgo_max_iterations; ++iteration) {
    Linearization linearization = BuildLinearization();
    Eigen::MatrixXd hessian = linearization.hessian;
    RegularizeInPlace(&hessian, 1e-6);
    const Eigen::VectorXd step =
        -hessian.ldlt().solve(linearization.gradient);
    if (!step.allFinite()) {
      break;
    }

    double max_step = 0.0;
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
      const Vec15 delta =
          step.segment<kErrorStateDim>(static_cast<Eigen::Index>(index) * kErrorStateDim);
      max_step = std::max(max_step, delta.norm());
    }
    const double step_scale =
        max_step > 1.0 ? (1.0 / max_step) : 1.0;
    double applied_max_step = 0.0;
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
      const Vec15 delta =
          step.segment<kErrorStateDim>(static_cast<Eigen::Index>(index) * kErrorStateDim) *
          step_scale;
      applied_max_step = std::max(applied_max_step, delta.norm());
      nodes_[index].state = ApplyError(nodes_[index].state, delta);
    }
    if (applied_max_step < 1e-5) {
      break;
    }
  }
}

void FgoEstimator::MarginalizeOldestNode() {
  if (nodes_.size() <= 1 || imu_edges_.empty()) {
    return;
  }

  Linearization linearization = BuildLinearization();
  const int total_dim = linearization.hessian.rows();
  if (total_dim <= kErrorStateDim) {
    return;
  }

  const Mat15 hmm =
      linearization.hessian.block<15, 15>(0, 0) + Mat15::Identity() * 1e-6;
  const Eigen::MatrixXd hmr =
      linearization.hessian.block(0, kErrorStateDim, kErrorStateDim, total_dim - kErrorStateDim);
  const Eigen::MatrixXd hrm =
      linearization.hessian.block(kErrorStateDim, 0, total_dim - kErrorStateDim, kErrorStateDim);
  const Eigen::MatrixXd hrr = linearization.hessian.block(
      kErrorStateDim, kErrorStateDim, total_dim - kErrorStateDim, total_dim - kErrorStateDim);
  const Vec15 bm = linearization.gradient.segment<15>(0);
  const Eigen::VectorXd br =
      linearization.gradient.segment(kErrorStateDim, total_dim - kErrorStateDim);

  const Eigen::MatrixXd solved_hmr = hmm.ldlt().solve(hmr);
  const Eigen::VectorXd solved_bm = hmm.ldlt().solve(bm);
  const Eigen::MatrixXd h_sc = hrr - hrm * solved_hmr;
  const Eigen::VectorXd b_sc = br - hrm * solved_bm;

  prior_.valid = true;
  prior_.state = nodes_[1].state;
  prior_.information = Symmetrize(h_sc.block<15, 15>(0, 0));
  RegularizeInPlace(&prior_.information, 1e-6);
  if (b_sc.size() >= kErrorStateDim) {
    const Vec15 delta =
        -prior_.information.ldlt().solve(b_sc.segment<15>(0));
    if (delta.allFinite()) {
      prior_.state = ApplyError(prior_.state, delta);
    }
  }

  nodes_.pop_front();
  imu_edges_.pop_front();
}

}  // namespace gnss_imu::internal
