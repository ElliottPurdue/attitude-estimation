/* Complementary filter behaviour.
 *
 * These assert what the filter must do, what it must refuse to do, and the one
 * thing it fundamentally cannot do. That last group matters most: a filter that
 * quietly reports a confident yaw from a 6-DOF IMU is worse than one that
 * admits it has none.
 */

#include "test.h"
#include "../src/complementary.h"

#include <math.h>

#define PI 3.14159265358979323846f
#define DEG (PI / 180.0f)

/* Body-frame accelerometer reading for a given true attitude, at rest and with
 * no noise: world +z rotated into the body frame. */
static vec3 gravity_in_body(quat truth)
{
    return quat_rotate_inverse(truth, vec3_make(0.0f, 0.0f, 1.0f));
}

/* Angle between the true and estimated gravity directions.
 *
 * This, and not the full quaternion angle, is what a 6-DOF IMU can be held to.
 * All orientations sharing a gravity direction produce identical accelerometer
 * readings, so they form an equivalence class the sensor cannot distinguish
 * within, and any two members differ by a rotation about vertical.
 *
 * The distinction is easy to miss because it is invisible for a pure roll or a
 * pure pitch, where the minimal-tilt rotation and the Z-Y-X Euler composition
 * happen to coincide. Combine the two and they diverge -- at roll -15 and pitch
 * 40 the filter recovers both angles exactly while sitting 5.5 degrees away in
 * yaw. Asserting on the full quaternion there measures the filter's failure to
 * guess information it was never given.
 */
static float tilt_error(quat estimate, quat truth)
{
    vec3 up_true = gravity_in_body(truth);
    vec3 up_estimate = gravity_in_body(estimate);

    float cosine = vec3_dot(up_true, up_estimate);
    if (cosine > 1.0f) {
        cosine = 1.0f;
    } else if (cosine < -1.0f) {
        cosine = -1.0f;
    }
    return acosf(cosine);
}

/* Runs the filter for `seconds` holding a fixed true attitude, feeding a
 * gyroscope that reports only `bias`. */
static void settle(complementary_filter *filter, quat truth, vec3 bias,
                   float seconds, float dt)
{
    int steps = (int)(seconds / dt);
    for (int i = 0; i < steps; ++i) {
        complementary_update(filter, bias, gravity_in_body(truth), dt);
    }
}

static void test_converges_to_a_static_attitude(void)
{
    /* Started at identity, given an accelerometer that says the device is
     * rolled 30 degrees. The estimate must find it. */
    quat truth = quat_from_euler(30.0f * DEG, 0.0f, 0.0f);

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);
    settle(&filter, truth, vec3_make(0, 0, 0), 20.0f, 0.005f);

    CHECK_BELOW(tilt_error(filter.orientation, truth), 1.0f * DEG);
}

static void test_converges_in_pitch_as_well_as_roll(void)
{
    quat truth = quat_from_euler(-15.0f * DEG, 40.0f * DEG, 0.0f);

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);
    settle(&filter, truth, vec3_make(0, 0, 0), 20.0f, 0.005f);

    CHECK_BELOW(tilt_error(filter.orientation, truth), 1.0f * DEG);

    /* Roll and pitch specifically, since tilt alone would not catch the two
     * being swapped. Yaw is deliberately not checked. */
    vec3 euler = quat_to_euler(filter.orientation);
    CHECK_NEAR(euler.x / DEG, -15.0f, 1.0f);
    CHECK_NEAR(euler.y / DEG, 40.0f, 1.0f);
}

static void test_initialising_from_accel_skips_the_transient(void)
{
    quat truth = quat_from_euler(25.0f * DEG, -35.0f * DEG, 0.0f);

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);
    complementary_set_from_accel(&filter, gravity_in_body(truth));

    /* Correct immediately, before a single update call. */
    CHECK_BELOW(tilt_error(filter.orientation, truth), 0.5f * DEG);
}

static void test_estimates_a_constant_gyro_bias(void)
{
    /* The gyroscope reports 3 deg/s about x while the device is motionless.
     * Without the integral term the estimate would tilt until the proportional
     * term balanced the false rate; with it, the bias is identified. */
    vec3 bias = vec3_make(3.0f * DEG, -2.0f * DEG, 0.0f);
    quat truth = quat_identity();

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP, 0.5f);
    settle(&filter, truth, bias, 300.0f, 0.005f);

    CHECK_BELOW(fabsf(filter.gyro_bias.x - bias.x), 0.5f * DEG);
    CHECK_BELOW(fabsf(filter.gyro_bias.y - bias.y), 0.5f * DEG);
    CHECK_BELOW(tilt_error(filter.orientation, truth), 1.0f * DEG);
}

static void test_bias_estimate_is_clamped(void)
{
    /* A wildly wrong gyro must not drive the bias state to infinity. */
    vec3 absurd = vec3_make(50.0f, 0.0f, 0.0f);

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP, 1.0f);
    settle(&filter, quat_identity(), absurd, 60.0f, 0.005f);

    CHECK(fabsf(filter.gyro_bias.x) <= filter.max_bias + 1e-6f);
    CHECK(!isnan(filter.orientation.w));
}

static void test_yaw_is_not_observable_from_a_six_dof_imu(void)
{
    /* The honest test. Gravity is unchanged by rotation about the vertical, so
     * no amount of accelerometer data can correct a yaw error. The filter must
     * converge in roll and pitch while leaving yaw exactly where it started.
     *
     * If this ever passes, the filter is inventing information it does not
     * have. */
    quat truth = quat_from_euler(0.0f, 0.0f, 90.0f * DEG);

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);
    settle(&filter, truth, vec3_make(0, 0, 0), 30.0f, 0.005f);

    vec3 euler = quat_to_euler(filter.orientation);
    CHECK_BELOW(fabsf(euler.x), 1.0f * DEG);      /* roll corrected */
    CHECK_BELOW(fabsf(euler.y), 1.0f * DEG);      /* pitch corrected */
    CHECK_BELOW(fabsf(euler.z), 1.0f * DEG);      /* yaw never moved */

    /* And so the total attitude error remains the full 90 degrees. */
    CHECK_NEAR(quat_angle_between(filter.orientation, truth) / DEG, 90.0f, 1.0f);
}

static void test_rejects_an_accelerometer_under_heavy_acceleration(void)
{
    /* A 3 g reading is not gravity. Accepting it would tilt the estimate
     * toward the resultant of gravity and the vehicle's own acceleration. */
    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);

    vec3 sideways = vec3_make(3.0f, 0.0f, 0.0f);
    for (int i = 0; i < 2000; ++i) {
        complementary_update(&filter, vec3_make(0, 0, 0), sideways, 0.005f);
    }

    CHECK_BELOW(quat_angle_between(filter.orientation, quat_identity()),
                0.5f * DEG);
}

static void test_rejects_free_fall(void)
{
    /* In free fall the accelerometer reads zero and carries no direction. */
    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);

    for (int i = 0; i < 2000; ++i) {
        complementary_update(&filter, vec3_make(0, 0, 0),
                             vec3_make(0, 0, 0), 0.005f);
    }

    CHECK(!isnan(filter.orientation.w));
    CHECK_NEAR(quat_norm(filter.orientation), 1.0f, 1e-4f);
}

static void test_tracks_a_rotating_body(void)
{
    /* Rotating steadily about x at 45 deg/s for two seconds, with a perfect
     * gyro. The filter should follow it, and the accelerometer correction must
     * not fight the motion. */
    float rate = 45.0f * DEG;
    float dt = 0.002f;
    int steps = 1000;

    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);

    quat truth = quat_identity();
    for (int i = 0; i < steps; ++i) {
        truth = quat_normalize(quat_multiply(
            truth, quat_from_rotation_vector(vec3_make(rate * dt, 0, 0))));
        complementary_update(&filter, vec3_make(rate, 0, 0),
                             gravity_in_body(truth), dt);
    }

    CHECK_BELOW(tilt_error(filter.orientation, truth), 2.0f * DEG);
}

static void test_orientation_stays_normalised_over_a_long_run(void)
{
    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);

    quat truth = quat_from_euler(10.0f * DEG, 20.0f * DEG, 0.0f);
    settle(&filter, truth, vec3_make(0.01f, -0.01f, 0.02f), 600.0f, 0.002f);

    CHECK_NEAR(quat_norm(filter.orientation), 1.0f, 1e-4f);
    CHECK(!isnan(filter.orientation.w));
}

static void test_a_zero_timestep_is_a_no_op(void)
{
    complementary_filter filter;
    complementary_init(&filter, COMPLEMENTARY_DEFAULT_KP,
                       COMPLEMENTARY_DEFAULT_KI);
    quat before = filter.orientation;

    complementary_update(&filter, vec3_make(1, 1, 1),
                         vec3_make(0, 0, 1), 0.0f);

    CHECK_NEAR(quat_angle_between(filter.orientation, before), 0.0f, 1e-9f);
}

void register_complementary_tests(void)
{
    RUN(test_converges_to_a_static_attitude);
    RUN(test_converges_in_pitch_as_well_as_roll);
    RUN(test_initialising_from_accel_skips_the_transient);
    RUN(test_estimates_a_constant_gyro_bias);
    RUN(test_bias_estimate_is_clamped);
    RUN(test_yaw_is_not_observable_from_a_six_dof_imu);
    RUN(test_rejects_an_accelerometer_under_heavy_acceleration);
    RUN(test_rejects_free_fall);
    RUN(test_tracks_a_rotating_body);
    RUN(test_orientation_stays_normalised_over_a_long_run);
    RUN(test_a_zero_timestep_is_a_no_op);
}
