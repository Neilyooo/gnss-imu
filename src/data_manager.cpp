#include "gnss_imu_fusion/fusion.h"

#include "internal/estimator_common.h"

namespace gnss_imu {

DataManager::DataManager(double max_buffer_seconds)
    : max_buffer_seconds_(max_buffer_seconds) {}

void DataManager::PushImu(const ImuData& imu) {
  imu_buffer_.push_back(imu);
  TrimBuffer(imu.timestamp);
}

std::optional<MeasureGroup> DataManager::PushGnss(const GnssData& gnss) {
  MeasureGroup group;
  group.is_initial = !last_gnss_time_.has_value();
  group.start_time = last_gnss_time_.value_or(gnss.timestamp);
  group.end_time = gnss.timestamp;
  group.gnss = gnss;

  while (!imu_buffer_.empty() &&
         imu_buffer_.front().timestamp + internal::kEps < group.start_time) {
    imu_buffer_.pop_front();
  }

  while (!imu_buffer_.empty() &&
         imu_buffer_.front().timestamp <= gnss.timestamp + internal::kEps) {
    if (group.is_initial ||
        imu_buffer_.front().timestamp > group.start_time + internal::kEps) {
      group.imu_samples.push_back(imu_buffer_.front());
    }
    imu_buffer_.pop_front();
  }

  last_gnss_time_ = gnss.timestamp;
  return group;
}

void DataManager::TrimBuffer(double current_time) {
  while (!imu_buffer_.empty() &&
         imu_buffer_.front().timestamp + max_buffer_seconds_ < current_time) {
    imu_buffer_.pop_front();
  }
}

}  // namespace gnss_imu
