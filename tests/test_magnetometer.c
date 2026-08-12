/* Heading correction from a magnetometer.
 *
 * Several of these are the deliberate inverses of assertions in test_ekf.c. That
 * suite pins what a 6-DOF sensor set cannot do -- yaw unobservable, its variance
 * growing without bound. Adding a magnetometer should turn each of those into its
 * opposite, and nothing else should change.
 *
 * The most important test here is the one asserting what the magnetometer must
 * NOT do: a magnetic disturbance must leave roll and pitch alone.
 */

#include "test.h"
#include "../src/ekf.h"
#include "../src/matrix.h"

#include <math.h>

#define PI 3.14159265358979323846f
#define DEG (PI / 180.0f)
#define N EKF_STATES

/* Local field: north and downward at a mid-latitude dip angle. World +x is the
 * horizontal direction of the field, +z is up, so the vertical component is
 * negative. */
#define DIP (60.0f * DEG)

static vec3 world_field(void)
{
    return vec3_make(cosf(DIP), 0.0f, -sinf(DIP));
}

static vec3 gravity_in_body(quat truth)
{
    return quat_rotate_inverse(truth, vec3_make(0.0f, 0.0f, 1.0f));
}

static vec3 field_in_body(quat truth)
{
    return quat_rotate_inverse(truth, world_field());
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

/* Runs both updates at a fixed true attitude. */
static void settle(ekf_filter *filter, quat truth, vec3 bias,
                   float seconds, float dt)
{
    int steps = (int)(seconds / dt);
    for (int i = 0; i < steps; ++i) {
        ekf_update(filter, bias, gravity_in_body(truth), dt);
        ekf_update_magnetometer(filter, field_in_body(truth));
    }
}

static void test_yaw_becomes_observable(void)
{
    /* test_ekf.c asserts that a 90 degree yaw error survives indefinitely on a
     * 6-DOF sensor set. With a magnetometer the same error must be corrected. */
    quat truth = quat_from_euler(0.0f, 0.0f, 90.0f * DEG);

    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, truth, vec3_make(0, 0, 0), 30.0f, 0.005f);

    CHECK_BELOW(quat_angle_between(filter.orientation, truth), 2.0f * DEG);
}

static void test_recovers_yaw_from_any_quadrant(void)
{
    /* Wrapping matters here: an estimate near +180 and a measurement near -180
     * are close together, and an unwrapped innovation would drive the correction
     * the long way round. */
    const float headings[] = { -170.0f, -95.0f, -20.0f, 45.0f, 130.0f, 179.0f };

    for (int i = 0; i < 6; ++i) {
        quat truth = quat_from_euler(10.0f * DEG, -15.0f * DEG,
                                     headings[i] * DEG);
        ekf_filter filter;
        ekf_init(&filter);
        settle(&filter, truth, vec3_make(0, 0, 0), 30.0f, 0.005f);

        CHECK_BELOW(quat_angle_between(filter.orientation, truth), 3.0f * DEG);
    }
}

/* Runs the filter along a constant-rate rotation, returning the true attitude
 * reached. The gyro is fed truth plus the bias it is supposed to discover. */
static quat run_motion(ekf_filter *filter, vec3 rate, vec3 bias,
                       float seconds, float dt)
{
    quat truth = quat_identity();
    int steps = (int)(seconds / dt);

    for (int i = 0; i < steps; ++i) {
        ekf_update(filter, vec3_add(rate, bias), gravity_in_body(truth), dt);
        ekf_update_magnetometer(filter, field_in_body(truth));
        truth = quat_normalize(quat_multiply(
            truth, quat_from_rotation_vector(vec3_scale(rate, dt))));
    }
    return truth;
}

static void test_tracks_heading_through_a_tumble(void)
{
    /* Regression test for a wrong measurement Jacobian.
     *
     * The error state is a body-frame rotation, so a heading measurement -- which
     * is about the WORLD vertical -- has H equal to that vertical expressed in
     * body coordinates. Using [0 0 1] instead is correct only while level.
     *
     * A static attitude will not reveal the difference: any gain with a positive
     * component along the true error still drives the innovation to zero given
     * enough time, so the filter converges anyway and the test passes. It takes
     * sustained motion through steep attitudes, where the wrong direction keeps
     * mis-attributing heading error into roll, pitch and bias faster than it can
     * settle. This tumble puts the sensor through every orientation.
     */
    vec3 rate = vec3_make(0.8f, 0.3f, 0.2f);
    vec3 bias = vec3_make(0.02f, -0.015f, 0.01f);

    ekf_filter filter;
    ekf_init(&filter);
    quat truth = run_motion(&filter, rate, bias, 60.0f, 0.005f);

    CHECK_BELOW(quat_angle_between(filter.orientation, truth), 3.0f * DEG);
    CHECK_BELOW(vec3_norm(vec3_sub(filter.gyro_bias, bias)), 1.0f * DEG);
}

static void test_yaw_uncertainty_now_shrinks(void)
{
    /* The converse of test_yaw_uncertainty_grows_without_bound. Heading is now
     * measured, so the filter's own confidence in it must improve. */
    ekf_filter filter;
    ekf_init(&filter);
    float initial = ekf_yaw_sigma(&filter);

    settle(&filter, quat_identity(), vec3_make(0, 0, 0), 20.0f, 0.005f);
    float settled = ekf_yaw_sigma(&filter);

    CHECK_BELOW(settled, initial);
    CHECK_BELOW(settled, 10.0f * DEG);
}

static void test_z_axis_bias_becomes_estimable(void)
{
    /* The static z-bias that neither filter could see without a magnetometer.
     * With heading observable, the bias driving heading is observable too. */
    vec3 bias = vec3_make(0.0f, 0.0f, 3.0f * DEG);

    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, quat_identity(), bias, 200.0f, 0.005f);

    CHECK_BELOW(fabsf(filter.gyro_bias.z - bias.z), 1.0f * DEG);
}

static void test_a_magnetic_disturbance_leaves_tilt_alone(void)
{
    /* The reason for the scalar, heading-only formulation.
     *
     * A large disturbance is applied to the magnetometer only. Roll and pitch
     * come from the accelerometer and are accurate; a three-axis magnetometer
     * update would let this pull them off, trading a working measurement for a
     * corrupted one. Heading may degrade -- that is the sensor being wrong --
     * but tilt must not move.
     */
    quat truth = quat_from_euler(20.0f * DEG, -30.0f * DEG, 0.0f);

    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, truth, vec3_make(0, 0, 0), 20.0f, 0.005f);
    float tilt_before = tilt_error(filter.orientation, truth);

    /* Something ferrous nearby: the field direction is badly wrong. */
    vec3 disturbed = quat_rotate_inverse(truth,
        vec3_make(cosf(DIP) * cosf(70.0f * DEG),
                  cosf(DIP) * sinf(70.0f * DEG),
                  -sinf(DIP)));

    for (int i = 0; i < 4000; ++i) {
        ekf_update(&filter, vec3_make(0, 0, 0), gravity_in_body(truth), 0.005f);
        ekf_update_magnetometer(&filter, disturbed);
    }

    float tilt_after = tilt_error(filter.orientation, truth);
    CHECK_BELOW(tilt_after, 1.0f * DEG);
    CHECK_BELOW(tilt_after - tilt_before, 0.5f * DEG);
}

static void test_ignores_a_reading_with_no_horizontal_component(void)
{
    /* A field pointing straight down carries no heading. Its horizontal
     * projection is noise, and using it would inject a random heading. */
    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, quat_identity(), vec3_make(0, 0, 0), 10.0f, 0.005f);
    quat before = filter.orientation;

    for (int i = 0; i < 1000; ++i) {
        ekf_update_magnetometer(&filter, vec3_make(0.0f, 0.0f, -1.0f));
    }

    CHECK_NEAR(quat_angle_between(filter.orientation, before), 0.0f, 1e-4f);
}

static void test_ignores_a_zero_reading(void)
{
    ekf_filter filter;
    ekf_init(&filter);
    quat before = filter.orientation;

    ekf_update_magnetometer(&filter, vec3_make(0.0f, 0.0f, 0.0f));

    CHECK_NEAR(quat_angle_between(filter.orientation, before), 0.0f, 1e-9f);
    CHECK(!isnan(filter.orientation.w));
}

static void test_covariance_stays_valid_with_both_updates(void)
{
    ekf_filter filter;
    ekf_init(&filter);
    settle(&filter, quat_from_euler(0.2f, -0.3f, 1.0f),
           vec3_make(0.01f, -0.01f, 0.02f), 300.0f, 0.002f);

    for (int i = 0; i < N; ++i) {
        CHECK(MAT_AT(filter.P, N, i, i) >= 0.0f);
        CHECK(!isnan(MAT_AT(filter.P, N, i, i)));
        for (int j = 0; j < N; ++j) {
            CHECK_NEAR(MAT_AT(filter.P, N, i, j), MAT_AT(filter.P, N, j, i), 1e-5f);
        }
    }
    CHECK_NEAR(quat_norm(filter.orientation), 1.0f, 1e-4f);
}

static void test_works_at_steep_dip_angles(void)
{
    /* Near the magnetic poles the field is almost vertical and its horizontal
     * component is small. Heading should still be recoverable while the
     * projection stays above the gate. */
    ekf_filter filter;
    ekf_init(&filter);
    quat truth = quat_from_euler(0.0f, 0.0f, 60.0f * DEG);

    vec3 steep = vec3_make(cosf(80.0f * DEG), 0.0f, -sinf(80.0f * DEG));
    for (int i = 0; i < 6000; ++i) {
        ekf_update(&filter, vec3_make(0, 0, 0), gravity_in_body(truth), 0.005f);
        ekf_update_magnetometer(&filter, quat_rotate_inverse(truth, steep));
    }

    CHECK_BELOW(quat_angle_between(filter.orientation, truth), 3.0f * DEG);
}

void register_magnetometer_tests(void)
{
    RUN(test_yaw_becomes_observable);
    RUN(test_recovers_yaw_from_any_quadrant);
    RUN(test_tracks_heading_through_a_tumble);
    RUN(test_yaw_uncertainty_now_shrinks);
    RUN(test_z_axis_bias_becomes_estimable);
    RUN(test_a_magnetic_disturbance_leaves_tilt_alone);
    RUN(test_ignores_a_reading_with_no_horizontal_component);
    RUN(test_ignores_a_zero_reading);
    RUN(test_covariance_stays_valid_with_both_updates);
    RUN(test_works_at_steep_dip_angles);
}
