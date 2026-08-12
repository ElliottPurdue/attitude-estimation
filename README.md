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
0.01 rad/s, accelerometer noise 0.02 g. Tilt error is the angle between the true
and estimated gravity directions.

| Motion | Complementary | EKF | |
|---|---|---|---|
| static | 0.407° | **0.140°** | 2.9x |
| gentle | 0.568° | **0.148°** | 3.8x |
| aggressive | 0.789° | **0.264°** | 3.0x |

Worst-case tilt error follows the same pattern: 0.97°/1.34°/1.61° against
0.44°/0.40°/0.65°.

Gyro bias, gentle profile, against a true `(0.020, -0.015, 0.010)`:

```
complementary   +0.0172  -0.0120  +0.0093
EKF             +0.0196  -0.0155  +0.0101
```

The complementary filter's fixed gain has to compromise between converging
quickly and rejecting noise. The EKF re-derives that tradeoff every step from its
own covariance, leaning on the accelerometer while uncertain and largely ignoring
it once confident — which is worth roughly a factor of three, and near-exact bias
recovery rather than a systematic 15% shortfall.

### Gyro bias observability depends on motion

The static case has the best tilt accuracy and the *worst* yaw drift, and the
complementary filter estimates its z-axis bias as exactly zero.

The accelerometer correction is a cross product between measured and expected
gravity. When the device is still, gravity lies along body-z, and the cross
product of two nearly-parallel z-vectors has no z-component — so nothing ever
informs the z-bias. Yaw then integrates it uncorrected.

Motion rotates gravity through the body frame and makes the bias partially
observable. Yaw drift over 60 seconds falls from 34.6° when static to 8.2° under
aggressive motion, purely because the filter could finally see what it needed to
cancel.

The EKF behaves differently and worse here, in an instructive way. Sensor noise
jitters the attitude slightly off level, giving the measurement Jacobian a small
non-zero yaw column, and the filter treats that noise-driven observability as
real information: its z-bias estimate converges to a value that varies with the
noise realization rather than toward the truth, while its stated variance shrinks.
Adding bias process noise does not fix it — the sweep from 0 to 5e-2 rad/s²/√Hz
changes the answer without improving it, because the problem is that the
information is absent, not that the filter is over-weighting it.

**A covariance filter can become confident about a quantity it cannot observe.**
The complementary filter has no covariance to be confident with, so it simply
leaves the state alone, which here is the better failure mode. Neither filter is
wrong; the sensor set is incomplete, and the fix is a magnetometer, not a better
estimator.

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
src/        vec3.h  quaternion.{h,c}  matrix.{h,c}            the library
            complementary.{h,c}  ekf.{h,c}
tests/      37 tests, no framework                            host only
sim/        generate.py  evaluate.py                          truth + scoring
bench/      run_filter.c                                       CSV driver
```

## Build and run

Requires a C99 compiler and `make`. Python 3 with NumPy for the simulator.

```
make test                     # 37 tests
make bench                    # build the CSV driver

python sim/generate.py --profile gentle --duration 60 --rate 200
python sim/evaluate.py --filter complementary
python sim/evaluate.py --filter ekf
```

`generate.py` takes `--gyro-bias`, `--gyro-noise`, `--accel-noise` and
`--linear-accel` to sweep sensor quality, and three motion profiles.

---

## Status

Working: quaternion algebra, complementary (Mahony-style) filter, a
multiplicative EKF over a 6-element error state, fixed-size matrix arithmetic,
simulation harness, error scoring, 37 tests.

Next: an Xtensa build measuring per-update cost on an ESP32-S3.

## License

MIT.
