/* Multiplicative Extended Kalman Filter for 6-DOF IMU attitude.
 *
 * Where the complementary filter applies a fixed correction gain, this one
 * derives the gain from how uncertain it currently is. When the estimate is
 * fresh or the gyro has been integrating unaided, it leans on the
 * accelerometer; once confident, it largely ignores it. That is worth the extra
 * arithmetic mainly during convergence and after disturbances.
 *
 * MULTIPLICATIVE, not additive. The state carried is a unit quaternion, but the
 * covariance is kept over a three-element *error* rotation rather than over the
 * quaternion's four components. A quaternion has four numbers describing three
 * degrees of freedom, so a 4x4 attitude covariance is necessarily rank
 * deficient, and enforcing unit norm afterwards is inconsistent with the update
 * that produced it. Tracking a small-angle error instead keeps the covariance
 * full rank and minimal, and the error is folded back into the reference
 * quaternion by multiplication at the end of each update -- hence the name.
 *
 * Error state, six elements:
 *   [0..2]  attitude error as a rotation vector, radians
 *   [3..5]  gyro bias error, rad/s
 *
 * The same caveat as every 6-DOF filter applies: gravity is invariant under
 * rotation about the vertical, so yaw is unobservable and its variance grows
 * without bound. The filter reports that honestly in P rather than pretending
 * to a heading it cannot have.
 */

#ifndef EKF_H
#define EKF_H

#include "quaternion.h"

#define EKF_STATES 6

typedef struct {
    quat orientation;       /* reference attitude, body to world */
    vec3 gyro_bias;         /* rad/s */

    float P[EKF_STATES * EKF_STATES];   /* error covariance */

    /* Process noise, as spectral densities rather than per-step variances, so
     * tuning survives a change of sample rate. */
    float gyro_noise;       /* rad/s/sqrt(Hz), random walk driving attitude */
    float bias_noise;       /* rad/s^2/sqrt(Hz), how fast bias is allowed to move */

    /* Measurement noise on the normalized gravity direction, per axis. */
    float accel_noise;

    float accel_tolerance;  /* fractional gate on |accel| vs gravity */
    float gravity;          /* magnitude at rest, in the caller's units */
    float max_bias;         /* rad/s, symmetric clamp */
} ekf_filter;

#define EKF_DEFAULT_GYRO_NOISE   0.01f
#define EKF_DEFAULT_BIAS_NOISE   0.0005f
#define EKF_DEFAULT_ACCEL_NOISE  0.05f
#define EKF_ACCEL_TOL            0.15f
#define EKF_GRAVITY_G            1.0f
#define EKF_MAX_BIAS             0.35f

void ekf_init(ekf_filter *filter);

/* Sets the reference attitude from one accelerometer sample, assuming rest, and
 * shrinks the tilt covariance to match that confidence. Yaw variance is left
 * large because the sample says nothing about heading. */
void ekf_set_from_accel(ekf_filter *filter, vec3 accel);

void ekf_update(ekf_filter *filter, vec3 gyro, vec3 accel, float dt);

/* One-sigma tilt uncertainty in radians, from the attitude block of P.
 * Useful for asserting that the filter's own confidence is calibrated. */
float ekf_tilt_sigma(const ekf_filter *filter);

/* One-sigma yaw uncertainty in radians. Expected to grow without bound on a
 * 6-DOF sensor set. */
float ekf_yaw_sigma(const ekf_filter *filter);

#endif /* EKF_H */
