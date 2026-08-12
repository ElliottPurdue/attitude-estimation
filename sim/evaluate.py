"""Score a filter's estimates against the truth the simulator recorded.

Reports the geodesic attitude error -- the single rotation angle separating the
estimate from truth -- rather than per-axis Euler differences, which exaggerate
error near vertical and are not a distance.

Tilt error is reported separately from total error because a 6-DOF IMU cannot
observe yaw. Total error will grow without bound as yaw drifts; tilt error is
the part the sensor set can actually be held responsible for, and it is the
number to compare filters on.
"""

import argparse
import pathlib
import subprocess
import sys

import numpy as np


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def quat_array(rows, prefix=""):
    return np.column_stack([rows[f"{prefix}qw"], rows[f"{prefix}qx"],
                            rows[f"{prefix}qy"], rows[f"{prefix}qz"]])


def geodesic_error(a, b):
    """Angle in radians between two arrays of quaternions.

    Uses atan2 of the error quaternion's parts rather than acos of its scalar,
    which loses roughly half its significant digits as the angle approaches
    zero -- exactly the regime being measured.
    """
    aw, ax, ay, az = a.T
    bw, bx, by, bz = b.T
    # conj(a) * b
    ew = aw * bw + ax * bx + ay * by + az * bz
    ex = aw * bx - ax * bw - ay * bz + az * by
    ey = aw * by + ax * bz - ay * bw - az * bx
    ez = aw * bz - ax * by + ay * bx - az * bw
    vector = np.sqrt(ex * ex + ey * ey + ez * ez)
    return 2.0 * np.arctan2(vector, np.abs(ew))


def gravity_in_body(quats):
    """World +z expressed in the body frame, for each quaternion."""
    w, x, y, z = quats.T
    return np.column_stack([
        2.0 * (x * z - w * y),
        2.0 * (y * z + w * x),
        1.0 - 2.0 * (x * x + y * y),
    ])


def tilt_error(estimate, truth):
    up_estimate = gravity_in_body(estimate)
    up_truth = gravity_in_body(truth)
    dot = np.clip(np.sum(up_estimate * up_truth, axis=1), -1.0, 1.0)
    return np.arccos(dot)


def summarize(name, truth_rows, estimate_rows, settle_fraction=0.2):
    truth = quat_array(truth_rows)
    estimate = quat_array(estimate_rows)
    count = min(len(truth), len(estimate))
    truth, estimate = truth[:count], estimate[:count]

    # The opening transient is excluded from the steady-state figures. It is
    # reported separately rather than dropped silently, since a filter that
    # converges slowly is a different problem from one that converges wrongly.
    start = int(count * settle_fraction)

    total = np.degrees(geodesic_error(estimate, truth))
    tilt = np.degrees(tilt_error(estimate, truth))

    return {
        "name": name,
        "samples": count,
        "tilt_rms": float(np.sqrt(np.mean(tilt[start:] ** 2))),
        "tilt_max": float(np.max(tilt[start:])),
        "tilt_final": float(tilt[-1]),
        "total_final": float(total[-1]),
        "settle_tilt": float(np.sqrt(np.mean(tilt[:start] ** 2))) if start else 0.0,
        "bias": (float(estimate_rows["bx"][-1]),
                 float(estimate_rows["by"][-1]),
                 float(estimate_rows["bz"][-1])),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", default=None,
                        help="IMU log from generate.py")
    parser.add_argument("--runner", default=None,
                        help="compiled bench/run_filter executable")
    parser.add_argument("--kp", type=float, default=1.0)
    parser.add_argument("--ki", type=float, default=0.05)
    args = parser.parse_args()

    here = pathlib.Path(__file__).parent
    data = pathlib.Path(args.data) if args.data else here / "data" / "gentle.csv"
    runner = pathlib.Path(args.runner) if args.runner else (
        here.parent / "build" / "run_filter.exe")

    if not data.exists():
        sys.exit(f"no IMU log at {data}; run sim/generate.py first")
    if not runner.exists():
        sys.exit(f"no runner at {runner}; run make bench first")

    result = subprocess.run([str(runner), str(data), str(args.kp), str(args.ki)],
                            capture_output=True, text=True, check=True)
    estimates = here / "data" / f"{data.stem}_estimates.csv"
    estimates.write_text(result.stdout, encoding="utf-8")

    stats = summarize(f"complementary kp={args.kp} ki={args.ki}",
                      load_csv(data), load_csv(estimates))

    print(f"\n  {data.name}  ({stats['samples']:,} samples)")
    print(f"  {'-' * 58}")
    print(f"  {stats['name']}")
    print(f"    tilt error, RMS steady-state   {stats['tilt_rms']:8.3f} deg")
    print(f"    tilt error, worst              {stats['tilt_max']:8.3f} deg")
    print(f"    tilt error, final              {stats['tilt_final']:8.3f} deg")
    print(f"    total error, final             {stats['total_final']:8.3f} deg"
          f"   (includes unobservable yaw)")
    print(f"    convergence-window tilt RMS    {stats['settle_tilt']:8.3f} deg")
    print(f"    estimated gyro bias            "
          f"{stats['bias'][0]:+.4f} {stats['bias'][1]:+.4f} "
          f"{stats['bias'][2]:+.4f} rad/s")
    print()


if __name__ == "__main__":
    main()
