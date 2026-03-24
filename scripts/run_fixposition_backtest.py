#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rosbags.rosbag2 import Reader
from rosbags.typesys import Stores, get_types_from_msg, get_typestore


IMU_TOPIC = "/fixposition/rawimu"
GNSS1_TOPIC = "/fixposition/gnss1"
GNSS2_TOPIC = "/fixposition/gnss2"
REFERENCE_TOPIC = "/fixposition/inspvax"
TARGET_TOPICS = {IMU_TOPIC, GNSS1_TOPIC, GNSS2_TOPIC, REFERENCE_TOPIC}


@dataclass
class GnssSample:
    timestamp: float
    latitude: float
    longitude: float
    altitude: float
    covariance: np.ndarray


@dataclass
class ReferenceSample:
    timestamp: float
    latitude: float
    longitude: float
    altitude: float
    velocity_enu: np.ndarray
    roll_deg: float
    pitch_deg: float
    azimuth_deg: float


@dataclass
class GeoReference:
    origin_latitude_deg: float
    origin_longitude_deg: float
    origin_altitude_m: float
    origin_ecef: np.ndarray
    ecef_to_enu: np.ndarray


def build_cli() -> argparse.ArgumentParser:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Run ESKF/UKF/FGO backtest on a Fixposition rosbag."
    )
    parser.add_argument("--bag_dir", required=True, type=Path, help="ROS2 bag directory")
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=repo_root / "results" / "fixposition_backtest",
        help="Directory for extracted CSVs, benchmark outputs, plots, and reports",
    )
    parser.add_argument(
        "--build_dir",
        type=Path,
        default=repo_root / "build",
        help="CMake build directory containing gnss_imu_bench",
    )
    parser.add_argument(
        "--algorithm",
        default="eskf,ukf,fgo",
        help="Algorithms to run, e.g. eskf,ukf,fgo or all",
    )
    parser.add_argument(
        "--skip_build",
        action="store_true",
        help="Skip cmake configure/build before running benchmark",
    )
    parser.add_argument(
        "--gnss_pair_tolerance",
        type=float,
        default=0.03,
        help="Maximum timestamp delta in seconds for pairing gnss1 and gnss2",
    )
    parser.add_argument(
        "--velocity_smoothing_window",
        type=int,
        default=5,
        help="Moving-average window for GNSS finite-difference velocities",
    )
    parser.add_argument(
        "--max_duration_s",
        type=float,
        default=0.0,
        help="If positive, keep only the first N seconds after the common time origin",
    )
    return parser


def load_typestore() -> object:
    typestore = get_typestore(Stores.ROS2_FOXY)
    msg_root = Path("/home/auto/common/install/msfl_msgs/share/msfl_msgs/msg")
    for message_name in ("HeaderBinaryData", "INSPVAXData"):
        msg_path = msg_root / f"{message_name}.msg"
        if msg_path.exists():
            typestore.register(
                get_types_from_msg(msg_path.read_text(), f"msfl_msgs/msg/{message_name}")
            )
    return typestore


def stamp_to_sec(stamp: object) -> float:
    return float(stamp.sec) + float(stamp.nanosec) * 1e-9


def lla_to_ecef(latitude_deg: float, longitude_deg: float, altitude_m: float) -> np.ndarray:
    a = 6378137.0
    e2 = 6.69437999014e-3
    lat = math.radians(latitude_deg)
    lon = math.radians(longitude_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    return np.array(
        [
            (n + altitude_m) * cos_lat * cos_lon,
            (n + altitude_m) * cos_lat * sin_lon,
            (n * (1.0 - e2) + altitude_m) * sin_lat,
        ],
        dtype=np.float64,
    )


def ecef_to_enu_matrix(latitude_deg: float, longitude_deg: float) -> np.ndarray:
    lat = math.radians(latitude_deg)
    lon = math.radians(longitude_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    return np.array(
        [
            [-sin_lon, cos_lon, 0.0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
        ],
        dtype=np.float64,
    )


def lla_to_enu(
    latitude_deg: float,
    longitude_deg: float,
    altitude_m: float,
    origin_latitude_deg: float,
    origin_longitude_deg: float,
    origin_ecef: np.ndarray,
    ecef_to_enu: np.ndarray,
) -> np.ndarray:
    ecef = lla_to_ecef(latitude_deg, longitude_deg, altitude_m)
    return ecef_to_enu @ (ecef - origin_ecef)


def enu_to_ecef(enu: np.ndarray, origin_ecef: np.ndarray, ecef_to_enu: np.ndarray) -> np.ndarray:
    return ecef_to_enu.T @ enu + origin_ecef


def ecef_to_lla(ecef: np.ndarray) -> tuple[float, float, float]:
    x = float(ecef[0])
    y = float(ecef[1])
    z = float(ecef[2])
    a = 6378137.0
    e2 = 6.69437999014e-3
    b = a * math.sqrt(1.0 - e2)
    ep2 = (a * a - b * b) / (b * b)
    p = math.sqrt(x * x + y * y)
    theta = math.atan2(a * z, b * p)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    longitude = math.atan2(y, x)
    latitude = math.atan2(
        z + ep2 * b * sin_theta * sin_theta * sin_theta,
        p - e2 * a * cos_theta * cos_theta * cos_theta,
    )
    sin_latitude = math.sin(latitude)
    n = a / math.sqrt(1.0 - e2 * sin_latitude * sin_latitude)
    altitude = p / max(math.cos(latitude), 1e-12) - n
    return math.degrees(latitude), math.degrees(longitude), altitude


def moving_average(samples: np.ndarray, window: int) -> np.ndarray:
    if window <= 1 or samples.shape[0] < window:
        return samples.copy()
    kernel = np.ones(window, dtype=np.float64) / float(window)
    result = np.empty_like(samples)
    for axis in range(samples.shape[1]):
        result[:, axis] = np.convolve(samples[:, axis], kernel, mode="same")
    return result


def finite_difference_velocity(timestamps: np.ndarray, positions: np.ndarray) -> np.ndarray:
    velocities = np.zeros_like(positions)
    if positions.shape[0] < 2:
        return velocities
    for index in range(positions.shape[0]):
        if index == 0:
            delta_t = max(timestamps[1] - timestamps[0], 1e-3)
            velocities[index] = (positions[1] - positions[0]) / delta_t
        elif index == positions.shape[0] - 1:
            delta_t = max(timestamps[-1] - timestamps[-2], 1e-3)
            velocities[index] = (positions[-1] - positions[-2]) / delta_t
        else:
            delta_t = max(timestamps[index + 1] - timestamps[index - 1], 1e-3)
            velocities[index] = (positions[index + 1] - positions[index - 1]) / delta_t
    return velocities


def pair_dual_gnss(
    gnss1: list[GnssSample], gnss2: list[GnssSample], tolerance_s: float
) -> list[tuple[float, GnssSample | None, GnssSample | None]]:
    pairs: list[tuple[float, GnssSample | None, GnssSample | None]] = []
    index1 = 0
    index2 = 0
    while index1 < len(gnss1) and index2 < len(gnss2):
        sample1 = gnss1[index1]
        sample2 = gnss2[index2]
        delta = sample1.timestamp - sample2.timestamp
        if abs(delta) <= tolerance_s:
            pairs.append((0.5 * (sample1.timestamp + sample2.timestamp), sample1, sample2))
            index1 += 1
            index2 += 1
            continue
        if delta < 0.0:
            pairs.append((sample1.timestamp, sample1, None))
            index1 += 1
        else:
            pairs.append((sample2.timestamp, None, sample2))
            index2 += 1

    for index in range(index1, len(gnss1)):
        sample = gnss1[index]
        pairs.append((sample.timestamp, sample, None))
    for index in range(index2, len(gnss2)):
        sample = gnss2[index]
        pairs.append((sample.timestamp, None, sample))
    pairs.sort(key=lambda item: item[0])
    return pairs


def read_fixposition_bag(
    bag_dir: Path,
) -> tuple[list[tuple[float, float, float, float, float, float, float]], list[GnssSample], list[GnssSample], list[ReferenceSample]]:
    typestore = load_typestore()
    imu_samples: list[tuple[float, float, float, float, float, float, float]] = []
    gnss1_samples: list[GnssSample] = []
    gnss2_samples: list[GnssSample] = []
    reference_samples: list[ReferenceSample] = []

    with Reader(str(bag_dir)) as reader:
        connections = [connection for connection in reader.connections if connection.topic in TARGET_TOPICS]
        for connection, _, rawdata in reader.messages(connections=connections):
            message = typestore.deserialize_cdr(rawdata, connection.msgtype)
            if connection.topic == IMU_TOPIC:
                imu_samples.append(
                    (
                        stamp_to_sec(message.header.stamp),
                        float(message.linear_acceleration.x),
                        float(message.linear_acceleration.y),
                        float(message.linear_acceleration.z),
                        float(message.angular_velocity.x),
                        float(message.angular_velocity.y),
                        float(message.angular_velocity.z),
                    )
                )
                continue

            if connection.topic in (GNSS1_TOPIC, GNSS2_TOPIC):
                covariance = np.asarray(message.position_covariance, dtype=np.float64).reshape(3, 3)
                sample = GnssSample(
                    timestamp=stamp_to_sec(message.header.stamp),
                    latitude=float(message.latitude),
                    longitude=float(message.longitude),
                    altitude=float(message.altitude),
                    covariance=covariance,
                )
                if connection.topic == GNSS1_TOPIC:
                    gnss1_samples.append(sample)
                else:
                    gnss2_samples.append(sample)
                continue

            reference_samples.append(
                ReferenceSample(
                    timestamp=stamp_to_sec(message.header.stamp),
                    latitude=float(message.latitude),
                    longitude=float(message.longitude),
                    altitude=float(message.height + message.undulation),
                    velocity_enu=np.array(
                        [
                            float(message.east_velocity),
                            float(message.north_velocity),
                            float(message.up_velocity),
                        ],
                        dtype=np.float64,
                    ),
                    roll_deg=float(message.roll),
                    pitch_deg=float(message.pitch),
                    azimuth_deg=float(message.azimuth),
                )
            )

    imu_samples.sort(key=lambda row: row[0])
    gnss1_samples.sort(key=lambda sample: sample.timestamp)
    gnss2_samples.sort(key=lambda sample: sample.timestamp)
    reference_samples.sort(key=lambda sample: sample.timestamp)
    return imu_samples, gnss1_samples, gnss2_samples, reference_samples


def build_dataset_frames(
    imu_samples: list[tuple[float, float, float, float, float, float, float]],
    gnss1_samples: list[GnssSample],
    gnss2_samples: list[GnssSample],
    reference_samples: list[ReferenceSample],
    tolerance_s: float,
    velocity_smoothing_window: int,
) -> tuple[GeoReference, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not imu_samples or not reference_samples or (not gnss1_samples and not gnss2_samples):
        raise RuntimeError("Bag does not contain enough IMU/GNSS/reference data to build a dataset.")

    gnss_pairs = pair_dual_gnss(gnss1_samples, gnss2_samples, tolerance_s)
    first_pair = gnss_pairs[0]
    first_latitude = (
        first_pair[1].latitude if first_pair[1] is not None else first_pair[2].latitude
    )
    first_longitude = (
        first_pair[1].longitude if first_pair[1] is not None else first_pair[2].longitude
    )
    first_altitude = (
        0.5 * (first_pair[1].altitude + first_pair[2].altitude)
        if first_pair[1] is not None and first_pair[2] is not None
        else (first_pair[1].altitude if first_pair[1] is not None else first_pair[2].altitude)
    )
    origin_ecef = lla_to_ecef(first_latitude, first_longitude, first_altitude)
    rotation = ecef_to_enu_matrix(first_latitude, first_longitude)
    georeference = GeoReference(
        origin_latitude_deg=first_latitude,
        origin_longitude_deg=first_longitude,
        origin_altitude_m=first_altitude,
        origin_ecef=origin_ecef,
        ecef_to_enu=rotation,
    )

    gnss_measurements_rows: list[dict[str, float]] = []
    for timestamp, sample1, sample2 in gnss_pairs:
        positions = []
        covariances = []
        source_count = 0
        for sample in (sample1, sample2):
            if sample is None:
                continue
            positions.append(
                lla_to_enu(
                    sample.latitude,
                    sample.longitude,
                    sample.altitude,
                    first_latitude,
                    first_longitude,
                    origin_ecef,
                    rotation,
                )
            )
            covariances.append(sample.covariance)
            source_count += 1
        position = np.mean(np.vstack(positions), axis=0)
        covariance = (
            sum(covariances) / float(source_count * source_count)
            if source_count > 1
            else covariances[0]
        )
        gnss_measurements_rows.append(
            {
                "timestamp": timestamp,
                "px": position[0],
                "py": position[1],
                "pz": position[2],
                "pos_std": math.sqrt(max(float(np.mean(np.diag(covariance))), 1e-6)),
                "source_count": float(source_count),
            }
        )

    gnss_frame = pd.DataFrame(gnss_measurements_rows)
    gnss_positions = gnss_frame[["px", "py", "pz"]].to_numpy(dtype=np.float64)
    gnss_timestamps = gnss_frame["timestamp"].to_numpy(dtype=np.float64)
    gnss_velocities = finite_difference_velocity(gnss_timestamps, gnss_positions)
    gnss_velocities = moving_average(gnss_velocities, velocity_smoothing_window)
    median_dt = float(np.median(np.diff(gnss_timestamps))) if len(gnss_timestamps) > 1 else 0.2
    velocity_std = np.maximum(
        0.8, 2.0 * gnss_frame["pos_std"].to_numpy(dtype=np.float64) / max(median_dt, 0.1)
    )
    gnss_frame["vx"] = gnss_velocities[:, 0]
    gnss_frame["vy"] = gnss_velocities[:, 1]
    gnss_frame["vz"] = gnss_velocities[:, 2]
    gnss_frame["vel_std"] = velocity_std

    reference_rows: list[dict[str, float]] = []
    for sample in reference_samples:
        position = lla_to_enu(
            sample.latitude,
            sample.longitude,
            sample.altitude,
            first_latitude,
            first_longitude,
            origin_ecef,
            rotation,
        )
        reference_rows.append(
            {
                "timestamp": sample.timestamp,
                "px": position[0],
                "py": position[1],
                "pz": position[2],
                "latitude": sample.latitude,
                "longitude": sample.longitude,
                "altitude": sample.altitude,
                "vx": sample.velocity_enu[0],
                "vy": sample.velocity_enu[1],
                "vz": sample.velocity_enu[2],
                "roll_deg": sample.roll_deg,
                "pitch_deg": sample.pitch_deg,
                "azimuth_deg": sample.azimuth_deg,
            }
        )

    imu_frame = pd.DataFrame(
        imu_samples, columns=["timestamp", "ax", "ay", "az", "gx", "gy", "gz"]
    )
    reference_frame = pd.DataFrame(reference_rows)

    time_origin = min(
        float(imu_frame["timestamp"].iloc[0]),
        float(gnss_frame["timestamp"].iloc[0]),
        float(reference_frame["timestamp"].iloc[0]),
    )
    for frame in (imu_frame, gnss_frame, reference_frame):
        frame["timestamp"] = frame["timestamp"] - time_origin

    return georeference, imu_frame, gnss_frame, reference_frame, pd.DataFrame(
        gnss_measurements_rows
    ).assign(timestamp=lambda frame: frame["timestamp"] - time_origin)


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.10f")


def write_origin_metadata(georeference: GeoReference, path: Path) -> None:
    payload = {
        "origin_latitude_deg": georeference.origin_latitude_deg,
        "origin_longitude_deg": georeference.origin_longitude_deg,
        "origin_altitude_m": georeference.origin_altitude_m,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def configure_and_build(repo_root: Path, build_dir: Path) -> None:
    subprocess.run(["cmake", "-S", str(repo_root), "-B", str(build_dir)], check=True, cwd=repo_root)
    subprocess.run(["cmake", "--build", str(build_dir), "-j4"], check=True, cwd=repo_root)


def run_benchmark(
    repo_root: Path,
    build_dir: Path,
    output_dir: Path,
    algorithm: str,
) -> None:
    bench = build_dir / "gnss_imu_bench"
    if not bench.exists():
        raise FileNotFoundError(f"Missing benchmark executable: {bench}")

    data_dir = output_dir / "dataset"
    benchmark_dir = output_dir / "benchmark"
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(bench),
        "--mode",
        "csv",
        "--imu",
        str(data_dir / "imu.csv"),
        "--gnss",
        str(data_dir / "gnss.csv"),
        "--gt",
        str(data_dir / "gt.csv"),
        "--algorithm",
        algorithm,
        "--output_dir",
        str(benchmark_dir),
    ]
    subprocess.run(command, check=True, cwd=repo_root)


def interpolate_columns(
    reference_frame: pd.DataFrame, query_timestamps: np.ndarray, columns: Iterable[str]
) -> np.ndarray:
    valid = (query_timestamps >= reference_frame["timestamp"].iloc[0]) & (
        query_timestamps <= reference_frame["timestamp"].iloc[-1]
    )
    timestamps = reference_frame["timestamp"].to_numpy(dtype=np.float64)
    stacked = []
    for column in columns:
        values = reference_frame[column].to_numpy(dtype=np.float64)
        interpolated = np.interp(query_timestamps[valid], timestamps, values)
        stacked.append(interpolated)
    return valid, np.vstack(stacked).T


def trajectory_frame_to_lla(frame: pd.DataFrame, georeference: GeoReference) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for row in frame.itertuples(index=False):
        enu = np.array([row.px, row.py, row.pz], dtype=np.float64)
        ecef = enu_to_ecef(enu, georeference.origin_ecef, georeference.ecef_to_enu)
        latitude, longitude, altitude = ecef_to_lla(ecef)
        rows.append(
            {
                "timestamp": float(row.timestamp),
                "latitude": latitude,
                "longitude": longitude,
                "altitude": altitude,
            }
        )
    return pd.DataFrame(rows)


def kml_color_for_name(name: str) -> str:
    return {
        "eskf": "ff1b9e77",
        "ukf": "ff7570b3",
        "fgo": "ffe7298a",
        "ground_truth_inspvax": "ff000000",
    }.get(name, "ff0055aa")


def write_kml_line(lla_frame: pd.DataFrame, path: Path, name: str) -> None:
    coordinates = "\n".join(
        f"{row.longitude:.10f},{row.latitude:.10f},{row.altitude:.3f}"
        for row in lla_frame.itertuples(index=False)
    )
    contents = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>{name}</name>
    <Style id="lineStyle">
      <LineStyle>
        <color>{kml_color_for_name(name)}</color>
        <width>4</width>
      </LineStyle>
    </Style>
    <Placemark>
      <name>{name}</name>
      <styleUrl>#lineStyle</styleUrl>
      <LineString>
        <tessellate>1</tessellate>
        <altitudeMode>absolute</altitudeMode>
        <coordinates>
{coordinates}
        </coordinates>
      </LineString>
    </Placemark>
  </Document>
</kml>
"""
    path.write_text(contents, encoding="utf-8")


def export_kml_outputs(
    benchmark_dir: Path, reference_frame: pd.DataFrame, georeference: GeoReference
) -> None:
    kml_dir = benchmark_dir / "kml"
    kml_dir.mkdir(parents=True, exist_ok=True)

    ground_truth_lla = reference_frame[["timestamp", "latitude", "longitude", "altitude"]].copy()
    ground_truth_lla.to_csv(
        kml_dir / "ground_truth_inspvax_lla.csv", index=False, float_format="%.10f"
    )
    write_kml_line(
        ground_truth_lla, kml_dir / "ground_truth_inspvax.kml", "ground_truth_inspvax"
    )

    for algorithm in ("eskf", "ukf", "fgo"):
        trajectory_path = benchmark_dir / f"{algorithm}_trajectory.csv"
        if not trajectory_path.exists():
            continue
        trajectory_frame = pd.read_csv(trajectory_path)
        lla_frame = trajectory_frame_to_lla(trajectory_frame, georeference)
        lla_frame.to_csv(
            kml_dir / f"{algorithm}_trajectory_lla.csv", index=False, float_format="%.10f"
        )
        write_kml_line(lla_frame, kml_dir / f"{algorithm}.kml", algorithm)


def compute_algorithm_summary(
    benchmark_dir: Path, reference_frame: pd.DataFrame
) -> pd.DataFrame:
    metrics_path = benchmark_dir / "metrics.csv"
    metrics_frame = pd.read_csv(metrics_path)
    rows: list[dict[str, float | str]] = []
    for _, metric in metrics_frame.iterrows():
        algorithm = str(metric["algorithm"])
        trajectory_path = benchmark_dir / f"{algorithm}_trajectory.csv"
        trajectory_frame = pd.read_csv(trajectory_path)
        timestamps = trajectory_frame["timestamp"].to_numpy(dtype=np.float64)
        valid, reference_position = interpolate_columns(
            reference_frame, timestamps, ("px", "py", "pz")
        )
        reference_velocity_valid, reference_velocity = interpolate_columns(
            reference_frame, timestamps, ("vx", "vy", "vz")
        )
        if not np.array_equal(valid, reference_velocity_valid):
            raise RuntimeError("Position and velocity interpolation masks do not match.")
        estimate_position = trajectory_frame.loc[valid, ["px", "py", "pz"]].to_numpy(dtype=np.float64)
        estimate_velocity = trajectory_frame.loc[valid, ["vx", "vy", "vz"]].to_numpy(dtype=np.float64)
        position_error = estimate_position - reference_position
        velocity_error = estimate_velocity - reference_velocity
        aligned_position_error = position_error - position_error[0]
        rows.append(
            {
                "algorithm": algorithm,
                "rmse_position_raw_m": float(np.sqrt(np.mean(np.sum(position_error * position_error, axis=1)))),
                "rmse_position_aligned_m": float(
                    np.sqrt(np.mean(np.sum(aligned_position_error * aligned_position_error, axis=1)))
                ),
                "rmse_velocity_mps": float(np.sqrt(np.mean(np.sum(velocity_error * velocity_error, axis=1)))),
                "final_position_error_m": float(np.linalg.norm(position_error[-1])),
                "cpu_utilization_percent": float(metric["cpu_utilization_percent"]),
                "avg_update_ms": float(metric["avg_update_ms"]),
                "max_update_ms": float(metric["max_update_ms"]),
                "updates": int(metric["updates"]),
            }
        )
    summary = pd.DataFrame(rows).sort_values("rmse_position_aligned_m").reset_index(drop=True)
    summary.to_csv(benchmark_dir / "metrics_summary.csv", index=False, float_format="%.6f")
    return summary


def plot_trajectories(benchmark_dir: Path, dataset_dir: Path) -> None:
    gt_frame = pd.read_csv(benchmark_dir / "ground_truth.csv")
    gnss_frame = pd.read_csv(dataset_dir / "gnss_measurements.csv")
    plt.figure(figsize=(11, 8))
    plt.plot(gt_frame["px"], gt_frame["py"], label="inspvax", linewidth=2.0, color="black")
    plt.scatter(
        gnss_frame["px"], gnss_frame["py"], label="dual-gnss midpoint", s=6, alpha=0.35, color="#d95f02"
    )
    for algorithm, color in (("eskf", "#1b9e77"), ("ukf", "#7570b3"), ("fgo", "#e7298a")):
        trajectory_path = benchmark_dir / f"{algorithm}_trajectory.csv"
        if trajectory_path.exists():
            frame = pd.read_csv(trajectory_path)
            plt.plot(frame["px"], frame["py"], label=algorithm, linewidth=1.5, color=color)
    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.title("Fixposition Backtest Trajectory Comparison")
    plt.legend()
    plt.grid(True, alpha=0.25)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(benchmark_dir / "trajectory_xy.png", dpi=180)
    plt.close()


def plot_position_errors(benchmark_dir: Path, reference_frame: pd.DataFrame) -> None:
    plt.figure(figsize=(12, 7))
    for algorithm, color in (("eskf", "#1b9e77"), ("ukf", "#7570b3"), ("fgo", "#e7298a")):
        trajectory_path = benchmark_dir / f"{algorithm}_trajectory.csv"
        if not trajectory_path.exists():
            continue
        trajectory_frame = pd.read_csv(trajectory_path)
        timestamps = trajectory_frame["timestamp"].to_numpy(dtype=np.float64)
        valid, reference_position = interpolate_columns(reference_frame, timestamps, ("px", "py", "pz"))
        estimate_position = trajectory_frame.loc[valid, ["px", "py", "pz"]].to_numpy(dtype=np.float64)
        position_error = estimate_position - reference_position
        error_norm = np.linalg.norm(position_error, axis=1)
        plt.plot(timestamps[valid], error_norm, label=algorithm, linewidth=1.3, color=color)
    plt.xlabel("Time [s]")
    plt.ylabel("3D Position Error [m]")
    plt.title("Position Error Relative to /fixposition/inspvax")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(benchmark_dir / "position_error.png", dpi=180)
    plt.close()


def write_summary_report(
    benchmark_dir: Path,
    summary_frame: pd.DataFrame,
    imu_frame: pd.DataFrame,
    gnss_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> None:
    duration_s = float(reference_frame["timestamp"].iloc[-1] - reference_frame["timestamp"].iloc[0])
    report_path = benchmark_dir / "summary.md"
    with report_path.open("w", encoding="utf-8") as report:
        report.write("# Fixposition GNSS-IMU Backtest Summary\n\n")
        report.write(
            "This run uses `/fixposition/rawimu` as the propagation source, fuses "
            "`/fixposition/gnss1` and `/fixposition/gnss2` into a midpoint GNSS observation, "
            "and evaluates against `/fixposition/inspvax`.\n\n"
        )
        report.write(
            "Important note: the bag contains `NavSatFix`, not raw pseudorange/Doppler, so the FGO backend "
            "is a real-time fixed-lag GNSS-position graph rather than the full tightly coupled paper formulation.\n\n"
        )
        report.write(
            f"- Duration: {duration_s:.2f} s\n"
            f"- IMU samples: {len(imu_frame)}\n"
            f"- GNSS epochs: {len(gnss_frame)}\n"
            f"- Reference samples: {len(reference_frame)}\n\n"
        )
        report.write(
            "KML exports:\n"
            "- `benchmark/kml/eskf.kml`\n"
            "- `benchmark/kml/ukf.kml`\n"
            "- `benchmark/kml/fgo.kml`\n"
            "- `benchmark/kml/ground_truth_inspvax.kml`\n\n"
        )
        report.write("## Metrics\n\n")
        header = (
            "| algorithm | raw pos rmse [m] | aligned pos rmse [m] | vel rmse [m/s] "
            "| final pos err [m] | cpu [%] | avg update [ms] | max update [ms] |\n"
        )
        separator = (
            "|---|---:|---:|---:|---:|---:|---:|---:|\n"
        )
        report.write(header)
        report.write(separator)
        for _, row in summary_frame.iterrows():
            report.write(
                "| "
                f"{row['algorithm']} | "
                f"{row['rmse_position_raw_m']:.3f} | "
                f"{row['rmse_position_aligned_m']:.3f} | "
                f"{row['rmse_velocity_mps']:.3f} | "
                f"{row['final_position_error_m']:.3f} | "
                f"{row['cpu_utilization_percent']:.3f} | "
                f"{row['avg_update_ms']:.3f} | "
                f"{row['max_update_ms']:.3f} |\n"
            )


def main() -> None:
    parser = build_cli()
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir.resolve()
    dataset_dir = output_dir / "dataset"
    benchmark_dir = output_dir / "benchmark"
    dataset_dir.mkdir(parents=True, exist_ok=True)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    imu_samples, gnss1_samples, gnss2_samples, reference_samples = read_fixposition_bag(
        args.bag_dir.resolve()
    )
    (
        georeference,
        imu_frame,
        gnss_frame,
        reference_frame,
        gnss_measurements_frame,
    ) = build_dataset_frames(
        imu_samples,
        gnss1_samples,
        gnss2_samples,
        reference_samples,
        args.gnss_pair_tolerance,
        args.velocity_smoothing_window,
    )

    if args.max_duration_s > 0.0:
        limit = float(args.max_duration_s)
        imu_frame = imu_frame.loc[imu_frame["timestamp"] <= limit].reset_index(drop=True)
        gnss_frame = gnss_frame.loc[gnss_frame["timestamp"] <= limit].reset_index(drop=True)
        reference_frame = reference_frame.loc[reference_frame["timestamp"] <= limit].reset_index(drop=True)
        gnss_measurements_frame = (
            gnss_measurements_frame.loc[gnss_measurements_frame["timestamp"] <= limit]
            .reset_index(drop=True)
        )

    write_frame(imu_frame, dataset_dir / "imu.csv")
    write_frame(gnss_frame[["timestamp", "px", "py", "pz", "vx", "vy", "vz", "pos_std", "vel_std"]], dataset_dir / "gnss.csv")
    write_frame(reference_frame[["timestamp", "px", "py", "pz", "vx", "vy", "vz"]], dataset_dir / "gt.csv")
    write_frame(reference_frame, dataset_dir / "inspvax_reference.csv")
    write_frame(gnss_measurements_frame, dataset_dir / "gnss_measurements.csv")
    write_origin_metadata(georeference, dataset_dir / "geodetic_origin.json")

    if not args.skip_build:
        configure_and_build(repo_root, args.build_dir.resolve())
    run_benchmark(repo_root, args.build_dir.resolve(), output_dir, args.algorithm)
    summary_frame = compute_algorithm_summary(benchmark_dir, reference_frame)
    plot_trajectories(benchmark_dir, dataset_dir)
    plot_position_errors(benchmark_dir, reference_frame)
    export_kml_outputs(benchmark_dir, reference_frame, georeference)
    write_summary_report(benchmark_dir, summary_frame, imu_frame, gnss_frame, reference_frame)


if __name__ == "__main__":
    main()
