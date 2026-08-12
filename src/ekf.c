#include "ekf.h"
#include "matrix.h"

#include <math.h>

#define N EKF_STATES

static float clampf(float value, float limit)
{
    if (value > limit) {
        return limit;
    }
    if (value < -limit) {
        return -limit;
    }
    return value;
}

/* Skew-symmetric matrix of v, so that skew(v) * u == cross(v, u). */
static void skew(vec3 v, float *out)
{
    MAT_AT(out, 3, 0, 0) =  0.0f;  MAT_AT(out, 3, 0, 1) = -v.z;  MAT_AT(out, 3, 0, 2) =  v.y;
    MAT_AT(out, 3, 1, 0) =  v.z;   MAT_AT(out, 3, 1, 1) =  0.0f; MAT_AT(out, 3, 1, 2) = -v.x;
    MAT_AT(out, 3, 2, 0) = -v.y;   MAT_AT(out, 3, 2, 1) =  v.x;  MAT_AT(out, 3, 2, 2) =  0.0f;
}

void ekf_init(ekf_filter *filter)
{
    filter->orientation = quat_identity();
    filter->gyro_bias = vec3_make(0.0f, 0.0f, 0.0f);

    filter->gyro_noise = EKF_DEFAULT_GYRO_NOISE;
    filter->bias_noise = EKF_DEFAULT_BIAS_NOISE;
    filter->accel_noise = EKF_DEFAULT_ACCEL_NOISE;
    filter->accel_tolerance = EKF_ACCEL_TOL;
    filter->gravity = EKF_GRAVITY_G;
    filter->max_bias = EKF_MAX_BIAS;

    /* Initial uncertainty: attitude unknown to about a radian, bias to a few
     * degrees per second. Starting too confident would make the filter reject
     * the very measurements it needs to converge. */
    mat_zero(filter->P, N, N);
    for (int i = 0; i < 3; ++i) {
        MAT_AT(filter->P, N, i, i) = 1.0f;
        MAT_AT(filter->P, N, i + 3, i + 3) = 0.01f;
    }
}

void ekf_set_from_accel(ekf_filter *filter, vec3 accel)
{
    vec3 measured = vec3_normalize(accel);
    if (vec3_norm(measured) < 0.5f) {
        return;
    }

    vec3 reference = vec3_make(0.0f, 0.0f, 1.0f);
    vec3 axis = vec3_cross(measured, reference);
    float sine = vec3_norm(axis);
    float cosine = vec3_dot(measured, reference);

    if (sine < 1e-6f) {
        filter->orientation = (cosine > 0.0f)
            ? quat_identity()
            : quat_make(0.0f, 1.0f, 0.0f, 0.0f);
    } else {
        float angle = atan2f(sine, cosine);
        filter->orientation =
            quat_from_rotation_vector(vec3_scale(vec3_normalize(axis), angle));
    }

    /* Tilt is now known to roughly the accelerometer's own noise; heading is
     * not known at all. Encoding that asymmetry is the point -- a filter told it
     * knows its yaw would never accept a heading correction later. */
    float tilt_variance = filter->accel_noise * filter->accel_noise;
    MAT_AT(filter->P, N, 0, 0) = tilt_variance;
    MAT_AT(filter->P, N, 1, 1) = tilt_variance;
    MAT_AT(filter->P, N, 2, 2) = 1.0f;
}

void ekf_update(ekf_filter *filter, vec3 gyro, vec3 accel, float dt)
{
    if (dt <= 0.0f) {
        return;
    }

    /* ---------------------------------------------------------------- */
    /* Predict                                                           */
    /* ---------------------------------------------------------------- */

    vec3 rate = vec3_sub(gyro, filter->gyro_bias);

    filter->orientation = quat_normalize(quat_multiply(
        filter->orientation, quat_from_rotation_vector(vec3_scale(rate, dt))));

    /* Error-state transition. Over one step the attitude error rotates with the
     * body and is corrupted by whatever bias error exists:
     *
     *   F = [ I - skew(rate)*dt   -I*dt ]
     *       [        0              I   ]
     *
     * The lower-left block is zero because bias error does not depend on
     * attitude error; the lower-right is identity because bias is modelled as a
     * random walk with no deterministic drift. */
    float F[N * N];
    mat_identity(F, N);

    float S3[9];
    skew(rate, S3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            MAT_AT(F, N, i, j) -= MAT_AT(S3, 3, i, j) * dt;
        }
        MAT_AT(F, N, i, i + 3) = -dt;
    }

    float FP[N * N], P_new[N * N];
    mat_mul(F, filter->P, FP, N, N, N);
    mat_mul_bt(FP, F, P_new, N, N, N);

    /* Process noise, scaled by dt so that the tuning constants are spectral
     * densities and stay meaningful if the sample rate changes. */
    float gyro_variance = filter->gyro_noise * filter->gyro_noise * dt;
    float bias_variance = filter->bias_noise * filter->bias_noise * dt;
    for (int i = 0; i < 3; ++i) {
        MAT_AT(P_new, N, i, i) += gyro_variance;
        MAT_AT(P_new, N, i + 3, i + 3) += bias_variance;
    }

    mat_copy(filter->P, P_new, N, N);
    mat_symmetrize(filter->P, N);

    /* ---------------------------------------------------------------- */
    /* Update, from the accelerometer                                     */
    /* ---------------------------------------------------------------- */

    float magnitude = vec3_norm(accel);
    if (magnitude < 1e-6f || filter->gravity < 1e-6f) {
        return;
    }
    if (fabsf(magnitude / filter->gravity - 1.0f) >= filter->accel_tolerance) {
        return;     /* not gravity; predict-only step */
    }

    vec3 measured = vec3_scale(accel, 1.0f / magnitude);
    vec3 expected = quat_rotate_inverse(filter->orientation,
                                        vec3_make(0.0f, 0.0f, 1.0f));

    /* Measurement Jacobian. A small body-frame rotation error d rotates the
     * predicted gravity direction by skew(expected) * d, and the measurement
     * says nothing directly about bias, so the bias block is zero. */
    float H[3 * N];
    mat_zero(H, 3, N);
    float Hs[9];
    skew(expected, Hs);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            MAT_AT(H, N, i, j) = MAT_AT(Hs, 3, i, j);
        }
    }

    /* S = H P H^T + R */
    float PHt[N * 3], S[9], HP[3 * N];
    mat_mul(H, filter->P, HP, 3, N, N);
    mat_mul_bt(HP, H, S, 3, N, 3);
    float accel_variance = filter->accel_noise * filter->accel_noise;
    for (int i = 0; i < 3; ++i) {
        MAT_AT(S, 3, i, i) += accel_variance;
    }

    float S_inv[9];
    if (!mat3_inverse(S, S_inv)) {
        return;     /* measurement carries nothing usable this step */
    }

    /* K = P H^T S^-1 */
    float K[N * 3];
    mat_mul_bt(filter->P, H, PHt, N, N, 3);
    mat_mul(PHt, S_inv, K, N, 3, 3);

    vec3 innovation = vec3_sub(measured, expected);
    float y[3] = { innovation.x, innovation.y, innovation.z };

    float correction[N];
    for (int i = 0; i < N; ++i) {
        correction[i] = 0.0f;
        for (int j = 0; j < 3; ++j) {
            correction[i] += MAT_AT(K, 3, i, j) * y[j];
        }
    }

    /* Inject the attitude error multiplicatively and the bias error additively.
     * This is what keeps the reference quaternion on the unit sphere: the error
     * is a rotation, and rotations compose by multiplication. */
    filter->orientation = quat_normalize(quat_multiply(
        filter->orientation,
        quat_from_rotation_vector(vec3_make(correction[0],
                                            correction[1],
                                            correction[2]))));

    filter->gyro_bias = vec3_make(
        clampf(filter->gyro_bias.x + correction[3], filter->max_bias),
        clampf(filter->gyro_bias.y + correction[4], filter->max_bias),
        clampf(filter->gyro_bias.z + correction[5], filter->max_bias));

    /* Joseph form: P = (I-KH) P (I-KH)^T + K R K^T.
     *
     * The shorter P = (I-KH)P is algebraically equal but subtracts two nearly
     * equal quantities, and in float that loses symmetry and eventually
     * positive-definiteness. Joseph is a sum of two quadratic forms, so it
     * cannot go indefinite through rounding alone. The extra cost buys a filter
     * that does not silently diverge after a few million updates. */
    float KH[N * N], IKH[N * N], temp[N * N], JP[N * N];
    mat_mul(K, H, KH, N, 3, N);
    mat_identity(IKH, N);
    mat_sub(IKH, KH, IKH, N, N);

    mat_mul(IKH, filter->P, temp, N, N, N);
    mat_mul_bt(temp, IKH, JP, N, N, N);

    float KR[N * 3], KRKt[N * N];
    mat_scale(K, accel_variance, KR, N, 3);
    mat_mul_bt(KR, K, KRKt, N, 3, N);

    mat_add(JP, KRKt, filter->P, N, N);
    mat_symmetrize(filter->P, N);
}

float ekf_tilt_sigma(const ekf_filter *filter)
{
    float variance = MAT_AT(filter->P, N, 0, 0) + MAT_AT(filter->P, N, 1, 1);
    return sqrtf(variance > 0.0f ? variance : 0.0f);
}

float ekf_yaw_sigma(const ekf_filter *filter)
{
    float variance = MAT_AT(filter->P, N, 2, 2);
    return sqrtf(variance > 0.0f ? variance : 0.0f);
}
