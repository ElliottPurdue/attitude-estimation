# Attitude Estimation in C

Quaternion attitude estimation from a 6- or 9-DOF IMU, written in dependency-free
C99 for microcontrollers, and validated against synthetic trajectories with known
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
| static | 0.402° | **0.130°** | 3.1x |
| gentle | 0.576° | **0.135°** | 4.3x |
| aggressive | 0.789° | **0.250°** | 3.2x |

Worst-case tilt error follows the same pattern: 0.98°/1.27°/1.58° against
0.32°/0.33°/0.51°.

Gyro bias, gentle profile, against a true `(0.020, -0.015, 0.010)`:

```
complementary   +0.0173  -0.0115  +0.0092
EKF             +0.0199  -0.0147  +0.0101
```

The complementary filter's fixed gain has to compromise between converging
quickly and rejecting noise. The EKF re-derives that tradeoff every step from its
own covariance, leaning on the accelerometer while uncertain and largely ignoring
it once confident. That is worth roughly a factor of three, plus near-exact bias
recovery rather than a systematic 15% shortfall.

### Heading, with and without a magnetometer

Yaw RMS error over the same runs. The 6-DOF columns are open-loop integration
with nothing to correct them, so their values are not accuracy figures so much as
a record of how far each drifted before the run ended.

| Motion | Complementary | EKF | EKF + magnetometer |
|---|---|---|---|
| static | 22.1° | 63.1° | **0.17°** |
| gentle | 12.0° | 19.9° | **0.32°** |
| aggressive | 5.8° | 11.2° | **0.46°** |

**Tilt accuracy is unchanged:** 0.130°/0.134°/0.250° against the 6-DOF EKF's
0.130°/0.135°/0.250°, a difference of one thousandth of a degree on one profile
and nothing at all on the other two. That is the design goal, not a
coincidence: the magnetometer update is scalar and its Jacobian points along the
world vertical, so every correction it can produce is a rotation about that
vertical. It is arithmetically incapable of moving roll or pitch, and a
disturbance test asserts it.

The z-axis gyro bias is the clearest evidence that the heading is now genuinely
observed rather than merely constrained. On the static profile, against a true
`+0.0100` rad/s:

```
EKF                  +0.0343     3.4x the true value, and still moving
EKF + magnetometer   +0.0099
```

### Gyro bias observability depends on motion

The static case has the best tilt accuracy and the *worst* yaw drift, and the
complementary filter estimates its z-axis bias as exactly zero.

The accelerometer correction is a cross product between measured and expected
gravity. When the device is still, gravity lies along body-z, and the cross
product of two nearly-parallel z-vectors has no z-component, so nothing ever
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
Adding bias process noise does not fix it. The sweep from 0 to 5e-2 rad/s²/√Hz
changes the answer without improving it, because the problem is that the
information is absent, not that the filter is over-weighting it.

**A covariance filter can become confident about a quantity it cannot observe.**
The complementary filter has no covariance to be confident with, so it simply
leaves the state alone, which here is the better failure mode. Neither filter is
wrong; the sensor set is incomplete, and the fix is a magnetometer, not a better
estimator.

Adding one confirms the diagnosis. The same EKF, same tuning, same log, with a
heading measurement supplied: the z-bias estimate goes from `+0.0343` to
`+0.0099` against a true `+0.0100`. Nothing about the estimator changed. The
missing information was the whole problem.

### Cost, and the tradeoff it forces

Measured on the host, and cross-compiled for the ESP32-S3 to size it:

| | us/update | Max rate | State | Flash (Xtensa) |
|---|---|---|---|---|
| complementary | 0.20 | ~5 MHz | 48 B | 1,160 B |
| EKF | 4.02 | 249 kHz | 204 B | 3,762 B |
| EKF + magnetometer | 5.51 | 182 kHz | 204 B | (same object) |

**The EKF costs 20x the compute and 4x the RAM for roughly 3x the accuracy.**
Whether that is worth paying depends entirely on the loop rate and the processor,
which is the point of measuring rather than assuming.

Heading correction adds 1.45 us per cycle, about 36% of the accelerometer update,
and no state at all, since the two extra floats in the struct are tuning parameters.
It is cheaper than the accelerometer update because the measurement is scalar:
a 1x1 innovation and a rank-one gain, so there is no matrix inverse. Most
magnetometers also run far slower than the gyro, so in practice the update is
called at a fraction of the loop rate and the amortised cost is lower still.

The whole library is **7.2 KB of flash and zero static RAM** on the ESP32-S3 --
no `.data`, no `.bss`, because there is no global state and nothing is
dynamically allocated. Every filter's state lives in a struct the caller owns.
Quaternion and matrix code accounts for the remaining 2.4 KB and is shared.

Flash figures are `.text` per object file, built `-Os` with
`xtensa-esp32s3-elf-gcc`, before the linker discards unused functions -- so they
are an upper bound on what a sketch using only one filter would pay:

```
xtensa-esp32s3-elf-gcc -std=c99 -Os -mlongcalls -c src/*.c
xtensa-esp32s3-elf-size -t *.o
```

`bench/esp32s3_bench/` is an Arduino sketch that runs the same measurement on
device, for the number that actually decides whether either filter fits a given
control loop.

The host benchmark binary is named `bench_cost` rather than anything containing
"update". Windows UAC applies installer-detection heuristics by filename, and
auto-elevates executables matching *update*, *setup*, *install* or *patch* -- the
original name failed to run with "requires elevation" on an ordinary shell.

### What each sensor set cannot do

Gravity is invariant under rotation about the vertical, so **a 6-DOF IMU carries
no yaw information at all.** Yaw is integrated open-loop and drifts without
bound. This is a property of the sensors, not a shortcoming of the algorithm, and
the test suite asserts it (`test_yaw_is_not_observable_from_a_six_dof_imu`)
rather than hiding it. Reported metrics separate the observable part (tilt) from
the total, so the drift is visible rather than buried in an average.

A magnetometer removes that limitation and introduces different ones:

- **It measures magnetic north, not true north.** Everything here is relative to
  the local field. Converting to true heading needs the local declination, which
  depends on where on Earth you are and is the caller's business.
- **Hard and soft iron distortion is not handled.** The simulator models sensor
  noise only. A real magnetometer needs calibration against the fixed offset from
  nearby ferrous material and the ellipsoidal scaling from field warping, and an
  uncalibrated one has a heading error that varies with orientation, the exact
  signature that looks like a filter bug.
- **A sustained disturbance degrades heading**, and the filter cannot tell it
  from a real rotation. It gates on the horizontal projection being long enough
  to carry a heading, which catches the field being cancelled or the sensor
  sitting near a magnetic pole, but a motor that merely bends the field will be
  believed. Roll and pitch stay correct throughout, which is the point of the
  scalar formulation.

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
drifts off the unit sphere at a rate proportional to angular velocity, worst
exactly when an estimator can least afford it. `test_integrating_many_small_steps_does_not_drift`
integrates 4,000 steps into one full revolution and checks the residual.

**The error metric uses `atan2`, not `acos`.** For float32, `acos(1-ε) ≈ √(2ε)`,
so rounding alone puts a 0.04° floor under any small-angle measurement. The
initial `acos` implementation reported **exactly zero** for rotations of 1e-6 and
1e-4 rad. It could not resolve small angles at all, which is the entire job of
an attitude error metric. Caught by the test suite before it reached the filters.

**The measurement Jacobian lives in the error state's frame.** The heading
measurement is a rotation about the *world* vertical, but the error state is a
rotation in the *body* frame, so H over the attitude block is the world vertical
expressed in body coordinates, not `[0 0 1]`. Getting this wrong is invisible
while level, where the two are the same vector, and the first nine magnetometer
tests all passed with the wrong version because they held the attitude near
level. The simulator caught it immediately: tilt error went from 0.14° to 10.6°
and the y-bias estimate ran into its clamp. `test_tracks_heading_through_a_tumble`
now pins it, and it has to be a *moving* test. At a fixed attitude even a wrong
gain direction eventually nulls the innovation, so a static test converges and
reports success.

**Float, not double.** On a Cortex-M or Xtensa core with a single-precision FPU,
doubles are emulated in software at roughly an order of magnitude more cost per
operation. That is the difference between closing a 1 kHz loop and not.

**No dynamic allocation, no I/O in `src/`.** The library compiles unchanged for a
microcontroller. File handling lives in `bench/`, which exists only to connect
the C filters to the Python harness.

---

## Layout

```
src/        vec3.h  quaternion.{h,c}  matrix.{h,c}            the library
            complementary.{h,c}  ekf.{h,c}
tests/      47 tests, no framework                            host only
sim/        generate.py  evaluate.py                          truth + scoring
bench/      run_filter.c  bench_cost.c                       CSV driver, timing
            esp32s3_bench/                                     on-device sketch
```

## Build and run

Requires a C99 compiler and `make`. Python 3 with NumPy for the simulator.

```
make test                     # 47 tests
make bench                    # build the CSV driver
make time                     # per-update cost and memory footprint

python sim/generate.py --profile gentle --duration 60 --rate 200
python sim/evaluate.py --filter complementary
python sim/evaluate.py --filter ekf
python sim/evaluate.py --filter ekf9      # EKF with the magnetometer
```

`generate.py` takes `--gyro-bias`, `--gyro-noise`, `--accel-noise`,
`--linear-accel`, `--mag-noise` and `--dip` to sweep sensor quality and magnetic
environment, and three motion profiles.

---

## Status

Working: quaternion algebra, complementary (Mahony-style) filter, a
multiplicative EKF over a 6-element error state, scalar magnetometer heading
correction, fixed-size matrix arithmetic, simulation harness, error scoring,
47 tests.

Cross-compiles clean for the ESP32-S3 at 5.9 KB flash and zero static RAM.
`bench/esp32s3_bench/` measures per-update cost on device; the figures in the
cost table above are from the host.

Possible next steps: hard and soft iron calibration, and a control law on top of
the estimator to make it a complete GNC stack.

## License

MIT.
