#include "gnss_imu_fusion/fusion.h"

#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> ParseAlgorithms(const std::string& value) {
  if (value == "all") {
    return {"eskf", "ukf", "fgo"};
  }

  std::vector<std::string> algorithms;
  std::stringstream stream(value);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (!token.empty()) {
      algorithms.push_back(token);
    }
  }
  return algorithms;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace gnss_imu;

  std::string mode = "synthetic";
  std::string output_dir = "results";
  std::string imu_path;
  std::string gnss_path;
  std::optional<std::string> gt_path;
  std::vector<std::string> algorithms = {"eskf", "ukf", "fgo"};
  double duration = 60.0;
  double imu_rate = 200.0;
  double gnss_rate = 10.0;
  std::uint32_t seed = 1;
  FusionConfig config;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      std::cout << Usage();
      return 0;
    }
    if (i + 1 >= argc) {
      std::cerr << "Missing value for argument: " << arg << '\n';
      std::cerr << Usage();
      return 1;
    }

    const std::string value = argv[++i];
    if (arg == "--mode") {
      mode = value;
    } else if (arg == "--algorithm") {
      algorithms = ParseAlgorithms(value);
    } else if (arg == "--output_dir") {
      output_dir = value;
    } else if (arg == "--duration") {
      duration = std::stod(value);
    } else if (arg == "--imu_rate") {
      imu_rate = std::stod(value);
    } else if (arg == "--gnss_rate") {
      gnss_rate = std::stod(value);
    } else if (arg == "--seed") {
      seed = static_cast<std::uint32_t>(std::stoul(value));
    } else if (arg == "--imu") {
      imu_path = value;
    } else if (arg == "--gnss") {
      gnss_path = value;
    } else if (arg == "--gt") {
      gt_path = value;
    } else if (arg == "--window_size") {
      config.fgo_window_size = std::stoi(value);
    } else if (arg == "--fgo_iterations") {
      config.fgo_max_iterations = std::stoi(value);
    } else {
      std::cerr << "Unknown argument: " << arg << '\n';
      std::cerr << Usage();
      return 1;
    }
  }

  try {
    Dataset dataset;
    if (mode == "synthetic") {
      dataset = GenerateSyntheticDataset(config, duration, imu_rate, gnss_rate, seed);
    } else if (mode == "csv") {
      if (imu_path.empty() || gnss_path.empty()) {
        throw std::invalid_argument(
            "CSV mode requires both --imu and --gnss");
      }
      dataset = LoadDatasetFromCsv(imu_path, gnss_path, gt_path);
    } else {
      throw std::invalid_argument("Unsupported mode: " + mode);
    }

    const auto metrics = RunBenchmark(dataset, config, algorithms, output_dir);
    std::cout << "algorithm,rmse_position,rmse_velocity,cpu_utilization_percent,"
                 "avg_update_ms,max_update_ms,updates\n";
    for (const auto& metric : metrics) {
      std::cout << metric.algorithm << ','
                << metric.rmse_position << ','
                << metric.rmse_velocity << ','
                << metric.cpu_utilization_percent << ','
                << metric.avg_update_ms << ','
                << metric.max_update_ms << ','
                << metric.updates << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
