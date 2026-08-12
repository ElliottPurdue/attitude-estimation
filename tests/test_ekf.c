/* EKF behaviour, and the matrix algebra underneath it.
 *
 * Beyond the accuracy checks the complementary filter also gets, these assert
 * things only a covariance filter can be held to: that its stated uncertainty
 * shrinks where the sensors inform it and grows where they do not, and that the
 * covariance stays a valid covariance over long runs.
 */

#include "test.h"
#include "../src/ekf.h"
#include "../src/matrix.h"

#include <math.h>

#define PI 3.14159265358979323846f
#define DEG (PI / 180.0f)
#define N EKF_STATES

static vec3 gravity_in_body(quat truth)
{
    return quat_rotate_inverse(truth, vec3_make(0.0f, 0.0f, 1.0f));
}

static float tilt_error(quat estimate, quat truth)
{
    vec3 a = gravity_in_body(truth), b = gravity_in_body(estimate);
    float cosine = vec3_dot(a, b);
    if (cosine > 1.0f) {
        cosine = 1.0f;
    } else if (cosine < -1.0f) {
        cosine = -1.0f;
    }
    return acosf(cosine);
}

static void settle(ekf_filter *filter, quat truth, vec3 bias,
                   float seconds, float dt)
{
    int steps = (int)(seconds / dt);
    for (int i = 0; i < steps; ++i) {
        ekf_update(filter, bias, gravity_in_body(truth), dt);
    }
}

/* ---------------------------------------------------------------- */
/* Matrix algebra                                                     */
/* ---------------------------------------------------------------- */

static void test_multiply_matches_hand_computation(void)
{
    float a[6] = { 1, 2, 3, 4, 5, 6 };      /* 2x3 */
    float b[6] = { 7, 8, 9, 10, 11, 12 };   /* 3x2 */
    float out[4];
    mat_mul(a, b, out, 2, 3, 2);

    CHECK_NEAR(out[0], 58.0f, 1e-5f);
    CHECK_NEAR(out[1], 64.0f, 1e-5f);
    CHECK_NEAR(out[2], 139.0f, 1e-5f);
    CHECK_NEAR(out[3], 154.0f, 1e-5f);
}

static void test_mul_bt_matches_transposing_first(void)
{
    float a[6] = { 1, 2, 3, 4, 5, 6 };      /* 2x3 */
    float b[6] = { 2, 0, 1, 3, 5, 7 };      /* 2x3, used as b^T (3x2) */
    float bt[6], expected[4], actual[4];

    mat_transpose(b, bt, 2, 3);
    mat_mul(a, bt, expected, 2, 3, 2);
    mat_mul_bt(a, b, actual, 2, 3, 2);

    for (int i = 0; i < 4; ++i) {
        CHECK_NEAR(actual[i], expected[i], 1e-5f);
    }
}

static void test_3x3_inverse_round_trips(void)
{
    float a[9] = { 4, 7, 2, 3, 6, 1, 2, 5, 3 };
    float inv[9], product[9];

    CHECK(mat3_inverse(a, inv));
    mat_mul(a, inv, product, 3, 3, 3);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            CHECK_NEAR(MAT_AT(product, 3, i, j), (i == j) ? 1.0f : 0.0f, 1e-4f);
        }
    }
}

static void test_singular_inverse_is_refused(void)
{
    /* Third row is the sum of the first two, so the matrix has no inverse.
     * Returning failure lets the filter skip the update instead of applying an
     * enormous gain built from garbage. */
    float singular[9] = { 1, 2, 3, 4, 5, 6, 5, 7, 9 };
    float out[9];
    CHECK(!mat3_inverse(singular, out));
}

static void test_symmetrize_removes_asymmetry(void)
{
    float a[9] = { 1, 2, 3, 0, 4, 5, 0, 0, 6 };
    mat_symmetrize(a, 3);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            CHECK_NEAR(MAT_AT(a, 3, i, j), MAT_AT(a, 3, j, i), 1e-6f);
        }
    }
    CHECK_NEAR(MAT_AT(a, 3, 0, 1), 1.0f, 1e-6f);    /* mean of 2 and 0 */
}

/* ---------------------------------------------------------------- */
/* Filter behaviour                                                   */
/* ---------------------------------------------------------------- */

static void test_converges_to_a_static_attitude(void)
{
    quat truth = quat_from_euler(-15.0f * DEG, 40.0f * DEG, 0.0f);

    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, truth, vec3_make(0, 0, 0), 20.0f, 0.005f);

    CHECK_BELOW(tilt_error(filter.orientation, truth), 1.0f * DEG);

    vec3 euler = quat_to_euler(filter.orientation);
    CHECK_NEAR(euler.x / DEG, -15.0f, 1.0f);
    CHECK_NEAR(euler.y / DEG, 40.0f, 1.0f);
}

static void test_estimates_a_constant_gyro_bias(void)
{
    vec3 bias = vec3_make(3.0f * DEG, -2.0f * DEG, 0.0f);

    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, quat_identity(), bias, 200.0f, 0.005f);

    CHECK_BELOW(fabsf(filter.gyro_bias.x - bias.x), 0.5f * DEG);
    CHECK_BELOW(fabsf(filter.gyro_bias.y - bias.y), 0.5f * DEG);
}

static void test_tilt_confidence_grows_as_it_converges(void)
{
    /* The distinguishing property of a covariance filter: it knows how well it
     * knows. Tilt uncertainty must fall once the accelerometer has been
     * informing it for a while. */
    ekf_filter filter;
    ekf_init(&filter);
    float initial = ekf_tilt_sigma(&filter);

    settle(&filter, quat_identity(), vec3_make(0, 0, 0), 10.0f, 0.005f);
    float settled = ekf_tilt_sigma(&filter);

    CHECK_BELOW(settled, initial);
    CHECK_BELOW(settled, 5.0f * DEG);
}

static void test_yaw_uncertainty_grows_without_bound(void)
{
    /* And the converse: nothing in a 6-DOF sensor set informs heading, so its
     * variance must increase. A filter reporting shrinking yaw confidence would
     * be fabricating information, and would later reject a real heading
     * correction as an outlier. */
    ekf_filter filter;
    ekf_init(&filter);
    ekf_set_from_accel(&filter, vec3_make(0.0f, 0.0f, 1.0f));

    float early = ekf_yaw_sigma(&filter);
    settle(&filter, quat_identity(), vec3_make(0, 0, 0), 60.0f, 0.005f);
    float late = ekf_yaw_sigma(&filter);

    CHECK(late >= early);
}

static void test_covariance_stays_symmetric_and_positive(void)
{
    /* Run long enough that a naive covariance update would have lost symmetry
     * and gone indefinite. Diagonal entries are variances and cannot be
     * negative. */
    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, quat_from_euler(0.2f, -0.3f, 0.0f),
           vec3_make(0.01f, -0.01f, 0.02f), 400.0f, 0.002f);

    for (int i = 0; i < N; ++i) {
        CHECK(MAT_AT(filter.P, N, i, i) >= 0.0f);
        CHECK(!isnan(MAT_AT(filter.P, N, i, i)));
        for (int j = 0; j < N; ++j) {
            CHECK_NEAR(MAT_AT(filter.P, N, i, j), MAT_AT(filter.P, N, j, i), 1e-5f);
        }
    }
    CHECK_NEAR(quat_norm(filter.orientation), 1.0f, 1e-4f);
}

static void test_rejects_heavy_acceleration(void)
{
    ekf_filter filter;
    ekf_init(&filter);
    ekf_set_from_accel(&filter, vec3_make(0.0f, 0.0f, 1.0f));

    for (int i = 0; i < 2000; ++i) {
        ekf_update(&filter, vec3_make(0, 0, 0), vec3_make(3.0f, 0, 0), 0.005f);
    }

    CHECK_BELOW(tilt_error(filter.orientation, quat_identity()), 1.0f * DEG);
}

static void test_survives_free_fall(void)
{
    ekf_filter filter;
    ekf_init(&filter);

    for (int i = 0; i < 2000; ++i) {
        ekf_update(&filter, vec3_make(0, 0, 0), vec3_make(0, 0, 0), 0.005f);
    }

    CHECK(!isnan(filter.orientation.w));
    CHECK_NEAR(quat_norm(filter.orientation), 1.0f, 1e-4f);
}

static void test_tracks_a_rotating_body(void)
{
    float rate = 45.0f * DEG, dt = 0.002f;

    ekf_filter filter;
    ekf_init(&filter);
    ekf_set_from_accel(&filter, vec3_make(0.0f, 0.0f, 1.0f));

    quat truth = quat_identity();
    for (int i = 0; i < 1000; ++i) {
        truth = quat_normalize(quat_multiply(
            truth, quat_from_rotation_vector(vec3_make(rate * dt, 0, 0))));
        ekf_update(&filter, vec3_make(rate, 0, 0), gravity_in_body(truth), dt);
    }

    CHECK_BELOW(tilt_error(filter.orientation, truth), 2.0f * DEG);
}

static void test_a_zero_timestep_is_a_no_op(void)
{
    ekf_filter filter;
    ekf_init(&filter);
    quat before = filter.orientation;

    ekf_update(&filter, vec3_make(1, 1, 1), vec3_make(0, 0, 1), 0.0f);
    CHECK_NEAR(quat_angle_between(filter.orientation, before), 0.0f, 1e-9f);
}

void register_ekf_tests(void)
{
    RUN(test_multiply_matches_hand_computation);
    RUN(test_mul_bt_matches_transposing_first);
    RUN(test_3x3_inverse_round_trips);
    RUN(test_singular_inverse_is_refused);
    RUN(test_symmetrize_removes_asymmetry);
    RUN(test_converges_to_a_static_attitude);
    RUN(test_estimates_a_constant_gyro_bias);
    RUN(test_tilt_confidence_grows_as_it_converges);
    RUN(test_yaw_uncertainty_grows_without_bound);
    RUN(test_covariance_stays_symmetric_and_positive);
    RUN(test_rejects_heavy_acceleration);
    RUN(test_survives_free_fall);
    RUN(test_tracks_a_rotating_body);
    RUN(test_a_zero_timestep_is_a_no_op);
}
