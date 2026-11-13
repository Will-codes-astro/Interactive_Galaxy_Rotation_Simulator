"""Visualise the baryonic orbit CSV produced by analysis/simulate_baryonic_orbits.py."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent.parent
CSV_PATH = BASE_DIR / "results" / "star_trajectories_baryonic.csv"
PLOTS_DIR = BASE_DIR / "plots"


def load_trajectories(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Trajectory file not found: {path}\n"
            "Run analysis/simulate_baryonic_orbits.py first."
        )

    data = pd.read_csv(path)
    required = {"frame", "time_myr", "star_id", "x_kpc", "y_kpc", "z_kpc"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Trajectory CSV is missing columns: {sorted(missing)}")
    if data.empty:
        raise ValueError("Trajectory CSV contains no rows.")
    return data.sort_values(["frame", "star_id"]).reset_index(drop=True)


def plot_snapshot(data: pd.DataFrame, output_path: Path, frame: int | None) -> None:
    selected_frame = int(data["frame"].max() if frame is None else frame)
    frame_data = data[data["frame"] == selected_frame]
    if frame_data.empty:
        raise ValueError(f"No trajectory data exists for frame {selected_frame}.")

    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    axis.scatter(
        frame_data["x_kpc"],
        frame_data["y_kpc"],
        frame_data["z_kpc"],
        s=12,
        alpha=0.85,
        label=f"Particles at {frame_data['time_myr'].iloc[0]:.1f} Myr",
    )
    axis.scatter([0], [0], [0], color="gold", edgecolor="black", s=90, label="Galaxy centre")

    limit = float(
        np.nanmax(
            np.abs(
                data[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)
            )
        )
    )
    limit = max(limit * 1.05, 1.0)
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_zlim(-limit, limit)
    axis.set_xlabel("x (kpc)")
    axis.set_ylabel("y (kpc)")
    axis.set_zlabel("z (kpc)")
    axis.set_title("NGC 3198 baryonic orbit model")
    axis.legend(loc="upper left")
    axis.text2D(
        0.02,
        0.02,
        "This model is currently a 2D mid-plane simulation: z = 0 for every particle.",
        transform=axis.transAxes,
        fontsize=9,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def animate_trajectories(
    data: pd.DataFrame,
    output_path: Path,
    frame_step: int,
    fps: int,
) -> None:
    source_frames = sorted(data["frame"].unique())[::frame_step]
    if source_frames[-1] != data["frame"].max():
        source_frames.append(int(data["frame"].max()))

    frame_data_by_frame = {
        frame: data[data["frame"] == frame].sort_values("star_id")
        for frame in source_frames
    }
    animation_frames = np.arange(len(source_frames) * 2 - 1)

    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    limit = float(
        np.nanmax(np.abs(data[["x_kpc", "y_kpc", "z_kpc"]].to_numpy(dtype=float)))
    )
    limit = max(limit * 1.05, 1.0)
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_zlim(-limit, limit)
    axis.set_xlabel("x (kpc)")
    axis.set_ylabel("y (kpc)")
    axis.set_zlabel("z (kpc)")
    axis.scatter([0], [0], [0], color="gold", edgecolor="black", s=90)
    scatter = axis.scatter([], [], [], s=12)
    title = axis.set_title("")

    def update(animation_frame: int):
        source_index = animation_frame // 2
        blend = 0.5 if animation_frame % 2 else 0.0
        frame_data = frame_data_by_frame[source_frames[source_index]]
        if blend:
            next_frame_data = frame_data_by_frame[source_frames[source_index + 1]]
            positions = (
                frame_data[["x_kpc", "y_kpc", "z_kpc"]].to_numpy()
                * (1.0 - blend)
                + next_frame_data[["x_kpc", "y_kpc", "z_kpc"]].to_numpy() * blend
            )
        else:
            positions = frame_data[["x_kpc", "y_kpc", "z_kpc"]].to_numpy()
        scatter._offsets3d = (
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
        )
        time_myr = frame_data["time_myr"].iloc[0]
        if blend:
            time_myr = (time_myr + next_frame_data["time_myr"].iloc[0]) / 2.0
        title.set_text(
            f"NGC 3198 baryonic orbit model — "
            f"{time_myr:.1f} Myr"
        )
        return scatter, title

    animation = FuncAnimation(
        figure,
        update,
        frames=animation_frames,
        interval=1000 / fps,
        blit=False,
    )
    animation.save(output_path, writer=PillowWriter(fps=fps))
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frame",
        type=int,
        default=None,
        help="Frame to plot; defaults to the final frame.",
    )
    parser.add_argument(
        "--animate",
        action="store_true",
        help="Also create an animated GIF of the particle positions.",
    )
    parser.add_argument(
        "--frame-step",
        type=int,
        default=3,
        help="Use every Nth frame in the GIF (default: 3).",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=20,
        help="GIF playback rate; 20 FPS with interpolation preserves the slower orbit timing (default: 20).",
    )
    args = parser.parse_args()
    if args.frame_step < 1:
        parser.error("--frame-step must be at least 1")
    if args.fps < 1:
        parser.error("--fps must be at least 1")

    data = load_trajectories(CSV_PATH)
    PLOTS_DIR.mkdir(exist_ok=True)
    snapshot_path = PLOTS_DIR / "ngc3198_baryonic_orbits_3d.png"
    plot_snapshot(data, snapshot_path, args.frame)
    print(f"[SUCCESS] Saved 3D snapshot to: {snapshot_path}")

    if args.animate:
        animation_path = PLOTS_DIR / "ngc3198_baryonic_orbits_3d.gif"
        animate_trajectories(data, animation_path, args.frame_step, fps=args.fps)
        print(f"[SUCCESS] Saved 3D animation to: {animation_path}")


if __name__ == "__main__":
    main()
