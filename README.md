# Attitude Estimation in C

Quaternion attitude estimation from a 6-DOF IMU, written in dependency-free C99
for microcontrollers, and validated against synthetic trajectories with known
ground truth.

The point of the simulator is that ground truth exists. On hardware you can watch
an estimate look plausible; you cannot say it was accurate to 0.4 degrees, because
nothing tells you the true attitude. Here the trajectory is prescribed, so every
error is measurable.

---

## Results

60 seconds at 200 Hz, gyro bias `(0.020, -0.015, 0.010)` rad/s, gyro noise
0.01 rad/s, accelerometer noise 0.02 g.

| Motion profile | Tilt RMS | Tilt worst | Yaw drift @60s | z-bias estimated |
|---|---|---|---|---|
| static | **0.407°** | 0.967° | 34.6° | +0.0000 |
| gentle | 0.568° | 1.344° | 15.2° | +0.0093 |
| aggressive | 0.789° | 1.611° | **8.2°** | +0.0040 |

*(true z-bias is 0.010 rad/s)*

**Tilt error** is the angle between the true and estimated gravity directions.
**Yaw drift** is the unobservable remainder — see below.

### Gyro bias observability depends on motion

The static case has the best tilt accuracy and the *worst* yaw drift, and the
z-axis bias is estimated as exactly zero. That is not a bug.

The accelerometer correction is a cross product between measured and expected
gravity. When the device is still, gravity lies along body-z, and the cross
product of two nearly-parallel z-vectors has no z-component — so the integral
term never sees evidence about z-bias and never learns it. Yaw then integrates
the full uncorrected bias.

Motion rotates gravity through the body frame, giving the correction a
z-component and making the bias partially observable. Under aggressive motion
yaw drift falls to 8.2°, a quarter of the static case, purely because the filter
could finally see the bias it needed to cancel.

The tradeoff runs the other way for tilt: aggressive motion nearly doubles tilt
RMS, because the correction spends its authority tracking the trajectory instead
of settling.

### What this sensor set cannot do

Gravity is invariant under rotation about the vertical, so **a 6-DOF IMU carries
no yaw information at all.** Yaw is integrated open-loop and drifts without
bound. Correcting it requires a magnetometer or an external heading reference.

This is a property of the sensors, not a shortcoming of the algorithm, and the
test suite asserts it (`test_yaw_is_not_observable_from_a_six_dof_imu`) rather
than hiding it. Reported metrics separate the observable part (tilt) from the
total, so the yaw drift is visible rather than buried in an average.

### Accelerometer rejection

Adding 0.4 g of body linear acceleration degrades tilt RMS from 0.568° to 3.958°
rather than diverging. The filter gates the accelerometer on magnitude: a reading
far from 1 g is not gravity, so it is discarded and the bias integrator is frozen
for that sample. Integrating during a long acceleration would poison the bias
estimate, and bias is the one state that persists after the disturbance ends.

---

## Design notes

**Hamilton convention, scalar-first, body-to-world.** Stated in the header
because the aerospace literature is split between Hamilton and JPL conventions,
and mixing them silently inverts rotations while everything still looks
self-consistent. Pinned by `test_ninety_degree_yaw_maps_x_to_y`.

**Exact exponential map, not `q + 0.5·ω·q`.** The first-order approximation
drifts off the unit sphere at a rate proportional to angular velocity — worst
exactly when an estimator can least afford it. `test_integrating_many_small_steps_does_not_drift`
integrates 4,000 steps into one full revolution and checks the residual.

**The error metric uses `atan2`, not `acos`.** For float32, `acos(1-ε) ≈ √(2ε)`,
so rounding alone puts a 0.04° floor under any small-angle measurement. The
initial `acos` implementation reported **exactly zero** for rotations of 1e-6 and
1e-4 rad — it could not resolve small angles at all, which is the entire job of
an attitude error metric. Caught by the test suite before it reached the filters.

**Float, not double.** On a Cortex-M or Xtensa core with a single-precision FPU,
doubles are emulated in software at roughly an order of magnitude more cost per
operation — the difference between closing a 1 kHz loop and not.

**No dynamic allocation, no I/O in `src/`.** The library compiles unchanged for a
microcontroller. File handling lives in `bench/`, which exists only to connect
the C filters to the Python harness.

---

## Layout

```
src/        vec3.h  quaternion.{h,c}  complementary.{h,c}     the library
tests/      23 tests, no framework                            host only
sim/        generate.py  evaluate.py                          truth + scoring
bench/      run_filter.c                                       CSV driver
```

## Build and run

Requires a C99 compiler and `make`. Python 3 with NumPy for the simulator.

```
make test                     # 23 tests
make bench                    # build the CSV driver

python sim/generate.py --profile gentle --duration 60 --rate 200
python sim/evaluate.py
```

`generate.py` takes `--gyro-bias`, `--gyro-noise`, `--accel-noise` and
`--linear-accel` to sweep sensor quality, and three motion profiles.

---

## Status

Working: quaternion algebra, complementary (Mahony-style) filter with gyro bias
estimation, simulation harness, error scoring, 23 tests.

Next: an Extended Kalman Filter to compare against on the same harness, and an
Xtensa build measuring per-update cost on an ESP32-S3.

## License

MIT.
