/* Complementary (Mahony-style) attitude filter for a 6-DOF IMU.
 *
 * Fuses two sensors with opposite failure modes. A gyroscope is accurate over
 * short intervals but its bias integrates into unbounded drift. An
 * accelerometer has no drift -- gravity does not wander -- but is corrupted by
 * any linear acceleration of the vehicle. Trusting the gyro at high frequency
 * and the accelerometer at low frequency gives an estimate better than either.
 *
 * The correction is applied as a rate added to the measured angular velocity,
 * rather than by blending two orientations. That keeps the filter on the
 * quaternion manifold throughout and gives the integral term a physical
 * meaning: it converges to the gyroscope's bias.
 *
 * WHAT THIS FILTER CANNOT DO: gravity is invariant under rotation about the
 * vertical, so a 6-DOF IMU carries no yaw information whatsoever. Yaw is
 * integrated open-loop and drifts without bound. Correcting it requires a
 * magnetometer or an external heading reference. This is a property of the
 * sensor set, not a shortcoming of the algorithm, and the test suite asserts it
 * rather than hiding it.
 */

#ifndef COMPLEMENTARY_H
#define COMPLEMENTARY_H

#include "quaternion.h"

typedef struct {
    quat orientation;       /* body-to-world rotation estimate */
    vec3 gyro_bias;         /* rad/s, estimated by the integral term */

    float kp;               /* proportional gain on the gravity error */
    float ki;               /* integral gain; drives bias estimation */

    /* Accelerometer readings whose magnitude departs from 1 g by more than this
     * fraction are discarded. Under linear acceleration the measured vector is
     * gravity plus that acceleration, and using it would drag the estimate
     * toward a direction that is not down. A vehicle in sustained coordinated
     * flight can fool this check, which is why the gate is loose enough to keep
     * ordinary motion and tight enough to reject impacts and free fall. */
    float accel_tolerance;

    /* Magnitude the accelerometer reports at rest, in the caller's units. The
     * filter uses only the direction of the vector, so this exists purely to
     * scale the trust gate above: set it to 1 for g's, or 9.80665 for m/s^2. */
    float gravity;

    float max_bias;         /* rad/s, symmetric clamp on the bias estimate */
} complementary_filter;

/* Sensible defaults for a consumer MEMS IMU at a few hundred Hz. */
#define COMPLEMENTARY_DEFAULT_KP     1.0f
#define COMPLEMENTARY_DEFAULT_KI     0.05f
#define COMPLEMENTARY_ACCEL_TOL      0.15f     /* +/- 15% of gravity */
#define COMPLEMENTARY_MAX_BIAS       0.35f     /* rad/s, about 20 deg/s */
#define COMPLEMENTARY_GRAVITY_G      1.0f      /* default: accel supplied in g */
#define COMPLEMENTARY_GRAVITY_MS2    9.80665f

void complementary_init(complementary_filter *filter, float kp, float ki);

/* Sets orientation from a single accelerometer sample, assuming the device is
 * at rest. Removes the initial convergence transient. Yaw is set to zero
 * because the sensor set cannot observe it. */
void complementary_set_from_accel(complementary_filter *filter, vec3 accel);

/* Advances the estimate by dt seconds.
 *
 * gyro is body-frame angular velocity in rad/s. accel is body-frame specific
 * force in any consistent unit; only its direction is used, and a zero or
 * out-of-tolerance vector skips the correction. */
void complementary_update(complementary_filter *filter,
                          vec3 gyro, vec3 accel, float dt);

#endif /* COMPLEMENTARY_H */
