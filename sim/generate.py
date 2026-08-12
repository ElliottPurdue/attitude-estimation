"""Generate synthetic IMU data from a known attitude trajectory.

The point of a simulator here is ground truth. On real hardware you can watch an
estimate look plausible; you cannot say it was accurate to 0.4 degrees, because
nothing tells you the true attitude. Here the trajectory is prescribed, so every
error is measurable.

Truth is produced by integrating a prescribed angular velocity at a much finer
step than the sensors are sampled at. Defining omega analytically and integrating
to get attitude guarantees the two are consistent; deriving omega by
differentiating a prescribed attitude introduces numerical error into the very
signal the filter is supposed to trust.
"""

import argparse
import csv
import pathlib

import numpy as np

# Integration substeps per sensor sample. The truth trajectory is integrated at
# this multiple of the output rate so that its own discretisation error stays far
# below the sensor noise being studied.
SUBSTEPS = 20

GRAVITY = np.array([0.0, 0.0, 1.0])   # world frame, in g

# Local magnetic field, unit length. World +x is its horizontal direction by
# definition of the frame, and +z is up, so a northern-hemisphere dip angle gives
# a negative vertical component.
DEFAULT_DIP_DEGREES = 60.0


def quat_multiply(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return np.array([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ])


def quat_from_rotation_vector(rotation):
    angle = np.linalg.norm(rotation)
    if angle < 1e-12:
        return np.array([1.0, 0.0, 0.0, 0.0])
    half = 0.5 * angle
    axis = rotation / angle * np.sin(half)
    return np.array([np.cos(half), axis[0], axis[1], axis[2]])


def quat_rotate_inverse(q, v):
    """World-frame vector into the body frame."""
    w, x, y, z = q
    conj = np.array([w, -x, -y, -z])
    u = conj[1:]
    t = np.cross(u, v) + conj[0] * v
    return v + 2.0 * np.cross(u, t)


# --------------------------------------------------------------------------
# Motion profiles
# --------------------------------------------------------------------------

def omega_static(t):
    """Perfectly still. Isolates bias and noise from any tracking error."""
    return np.zeros(3)


def omega_gentle(t):
    """Slow multi-axis weaving, of the kind a hand-held device produces."""
    return np.array([
        0.6 * np.sin(0.7 * t),
        0.4 * np.sin(0.5 * t + 1.0),
        0.3 * np.sin(0.3 * t + 2.0),
    ])


def omega_aggressive(t):
    """Fast, large-amplitude motion. Exposes integrator error and any filter
    whose correction gain is too slow to keep up."""
    return np.array([
        3.0 * np.sin(2.3 * t),
        2.5 * np.sin(1.9 * t + 0.5),
        2.0 * np.sin(3.1 * t + 1.5),
    ])


PROFILES = {
    "static": omega_static,
    "gentle": omega_gentle,
    "aggressive": omega_aggressive,
}


def magnetic_field(dip_degrees=DEFAULT_DIP_DEGREES):
    dip = np.radians(dip_degrees)
    return np.array([np.cos(dip), 0.0, -np.sin(dip)])


def generate(profile, duration, rate, gyro_noise, accel_noise, gyro_bias,
             linear_accel, seed, mag_noise=0.02, dip_degrees=DEFAULT_DIP_DEGREES):
    """Return (times, gyro, accel, mag, truth_quats) sampled at `rate` Hz."""
    generator = np.random.default_rng(seed)
    omega_of = PROFILES[profile]

    dt = 1.0 / rate
    fine_dt = dt / SUBSTEPS
    samples = int(duration * rate)

    q = np.array([1.0, 0.0, 0.0, 0.0])
    times = np.zeros(samples)
    truth = np.zeros((samples, 4))
    gyro = np.zeros((samples, 3))
    accel = np.zeros((samples, 3))
    mag = np.zeros((samples, 3))
    field = magnetic_field(dip_degrees)

    bias = np.asarray(gyro_bias, dtype=float)

    for i in range(samples):
        t = i * dt
        times[i] = t
        truth[i] = q

        true_rate = omega_of(t)

        # Gyro: truth plus a constant bias plus white noise. The bias is what
        # makes an unaided integration diverge, and what the filters must find.
        gyro[i] = true_rate + bias + generator.normal(0.0, gyro_noise, 3)

        # Accelerometer: gravity in the body frame, plus optional body-frame
        # linear acceleration, plus noise. With linear_accel above zero the
        # measurement is no longer purely gravity, which is what the filters'
        # magnitude gate exists to detect.
        specific_force = GRAVITY.copy()
        if linear_accel > 0.0:
            specific_force = specific_force + generator.normal(
                0.0, linear_accel, 3)
        accel[i] = (quat_rotate_inverse(q, specific_force)
                    + generator.normal(0.0, accel_noise, 3))

        # Magnetometer: the local field in the body frame. Noise only -- hard and
        # soft iron distortion would be a separate, systematic effect.
        mag[i] = (quat_rotate_inverse(q, field)
                  + generator.normal(0.0, mag_noise, 3))

        # Advance truth with substeps, using the noiseless rate.
        for step in range(SUBSTEPS):
            sub_t = t + step * fine_dt
            q = quat_multiply(q, quat_from_rotation_vector(
                omega_of(sub_t) * fine_dt))
            q = q / np.linalg.norm(q)

    return times, gyro, accel, mag, truth


def write_csv(path, times, gyro, accel, mag, truth):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t", "gx", "gy", "gz", "ax", "ay", "az",
                         "mx", "my", "mz", "qw", "qx", "qy", "qz"])
        for i in range(len(times)):
            writer.writerow([f"{times[i]:.6f}"]
                            + [f"{v:.9g}" for v in gyro[i]]
                            + [f"{v:.9g}" for v in accel[i]]
                            + [f"{v:.9g}" for v in mag[i]]
                            + [f"{v:.9g}" for v in truth[i]])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--profile", default="gentle", choices=sorted(PROFILES))
    parser.add_argument("--duration", type=float, default=60.0, help="seconds")
    parser.add_argument("--rate", type=float, default=200.0, help="Hz")
    parser.add_argument("--gyro-noise", type=float, default=0.01,
                        help="rad/s, one sigma")
    parser.add_argument("--accel-noise", type=float, default=0.02,
                        help="g, one sigma")
    parser.add_argument("--gyro-bias", type=float, nargs=3,
                        default=[0.02, -0.015, 0.01], metavar=("X", "Y", "Z"),
                        help="rad/s, constant")
    parser.add_argument("--linear-accel", type=float, default=0.0,
                        help="g, one sigma of body linear acceleration")
    parser.add_argument("--mag-noise", type=float, default=0.02,
                        help="one sigma, as a fraction of field strength")
    parser.add_argument("--dip", type=float, default=DEFAULT_DIP_DEGREES,
                        help="magnetic dip angle in degrees")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    times, gyro, accel, mag, truth = generate(
        args.profile, args.duration, args.rate, args.gyro_noise,
        args.accel_noise, args.gyro_bias, args.linear_accel, args.seed,
        args.mag_noise, args.dip)

    out = pathlib.Path(args.out) if args.out else (
        pathlib.Path(__file__).parent / "data" / f"{args.profile}.csv")
    write_csv(out, times, gyro, accel, mag, truth)

    print(f"wrote {out}")
    print(f"  {len(times):,} samples at {args.rate:g} Hz "
          f"({args.duration:g} s, profile '{args.profile}')")
    print(f"  gyro bias {args.gyro_bias} rad/s, "
          f"gyro noise {args.gyro_noise:g} rad/s, "
          f"accel noise {args.accel_noise:g} g")


if __name__ == "__main__":
    main()
