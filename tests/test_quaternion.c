/* Quaternion algebra.
 *
 * These are the invariants everything above depends on. If the exponential map
 * drifts off the unit sphere, or the rotation convention is inverted, the
 * filters still run and still produce plausible-looking output -- which is
 * precisely why the algebra is pinned here rather than checked by eye.
 */

#include "test.h"
#include "../src/quaternion.h"

#include <math.h>

#define PI 3.14159265358979323846f

static void test_identity_is_unit_and_neutral(void)
{
    quat identity = quat_identity();
    CHECK_NEAR(quat_norm(identity), 1.0f, 1e-6f);

    quat q = quat_normalize(quat_make(0.3f, -0.5f, 0.7f, 0.2f));
    quat left = quat_multiply(identity, q);
    CHECK_NEAR(quat_angle_between(left, q), 0.0f, 1e-6f);
}

static void test_conjugate_undoes_a_rotation(void)
{
    quat q = quat_from_euler(0.3f, -0.7f, 1.1f);
    quat round_trip = quat_multiply(q, quat_conjugate(q));
    CHECK_NEAR(quat_angle_between(round_trip, quat_identity()), 0.0f, 1e-5f);
}

static void test_rotation_preserves_vector_length(void)
{
    quat q = quat_from_euler(0.9f, 0.4f, -1.2f);
    vec3 v = vec3_make(1.0f, -2.0f, 3.0f);
    vec3 rotated = quat_rotate(q, v);
    CHECK_NEAR(vec3_norm(rotated), vec3_norm(v), 1e-5f);
}

static void test_inverse_rotation_returns_the_original_vector(void)
{
    quat q = quat_from_euler(-0.4f, 1.0f, 0.25f);
    vec3 v = vec3_make(0.2f, 5.0f, -1.5f);
    vec3 there_and_back = quat_rotate_inverse(q, quat_rotate(q, v));

    CHECK_NEAR(there_and_back.x, v.x, 1e-4f);
    CHECK_NEAR(there_and_back.y, v.y, 1e-4f);
    CHECK_NEAR(there_and_back.z, v.z, 1e-4f);
}

static void test_ninety_degree_yaw_maps_x_to_y(void)
{
    /* Pins the rotation convention. A quaternion built from a +90 degree yaw
     * must carry the body x-axis onto the world y-axis. If this reads -y, the
     * library is using the opposite handedness and every downstream result is
     * mirrored while still looking self-consistent. */
    quat yaw90 = quat_from_euler(0.0f, 0.0f, PI / 2.0f);
    vec3 mapped = quat_rotate(yaw90, vec3_make(1.0f, 0.0f, 0.0f));

    CHECK_NEAR(mapped.x, 0.0f, 1e-5f);
    CHECK_NEAR(mapped.y, 1.0f, 1e-5f);
    CHECK_NEAR(mapped.z, 0.0f, 1e-5f);
}

static void test_euler_round_trip(void)
{
    float roll = 0.35f, pitch = -0.62f, yaw = 2.10f;
    vec3 recovered = quat_to_euler(quat_from_euler(roll, pitch, yaw));

    CHECK_NEAR(recovered.x, roll, 1e-4f);
    CHECK_NEAR(recovered.y, pitch, 1e-4f);
    CHECK_NEAR(recovered.z, yaw, 1e-4f);
}

static void test_euler_survives_near_vertical_pitch(void)
{
    /* Just short of gimbal lock, where asinf's argument reaches 1 and an
     * unclamped implementation returns NaN. */
    quat q = quat_from_euler(0.0f, PI / 2.0f - 1e-4f, 0.0f);
    vec3 euler = quat_to_euler(q);

    CHECK(!isnan(euler.x) && !isnan(euler.y) && !isnan(euler.z));
    CHECK_NEAR(euler.y, PI / 2.0f - 1e-4f, 1e-3f);
}

static void test_rotation_vector_matches_a_known_rotation(void)
{
    /* A rotation vector of pi/2 about z is the same rotation as a 90 degree
     * yaw, which is what ties the exponential map to the Euler constructor. */
    quat exponential = quat_from_rotation_vector(vec3_make(0.0f, 0.0f, PI / 2.0f));
    quat euler = quat_from_euler(0.0f, 0.0f, PI / 2.0f);
    CHECK_NEAR(quat_angle_between(exponential, euler), 0.0f, 1e-5f);
}

static void test_rotation_vector_stays_on_the_unit_sphere(void)
{
    /* The first-order approximation q + 0.5*omega*q leaves the unit sphere by
     * an amount that grows with rate. The exponential map does not, and this
     * checks it at a rate high enough to expose the difference. */
    quat q = quat_from_rotation_vector(vec3_make(0.0f, 0.0f, 1.0f));
    CHECK_NEAR(quat_norm(q), 1.0f, 1e-6f);

    quat tiny = quat_from_rotation_vector(vec3_make(1e-9f, 0.0f, 0.0f));
    CHECK_NEAR(quat_norm(tiny), 1.0f, 1e-6f);
}

static void test_integrating_many_small_steps_does_not_drift(void)
{
    /* 4000 steps of 0.09 degrees about z should land on exactly one full turn.
     * Drift in the integrator shows up here as a nonzero residual angle. */
    quat q = quat_identity();
    vec3 step = vec3_make(0.0f, 0.0f, (2.0f * PI) / 4000.0f);
    for (int i = 0; i < 4000; ++i) {
        q = quat_multiply(q, quat_from_rotation_vector(step));
    }

    CHECK_NEAR(quat_norm(q), 1.0f, 1e-4f);
    CHECK_NEAR(quat_angle_between(q, quat_identity()), 0.0f, 1e-3f);
}

static void test_normalize_rejects_a_degenerate_quaternion(void)
{
    quat zero = quat_make(0.0f, 0.0f, 0.0f, 0.0f);
    quat result = quat_normalize(zero);
    CHECK_NEAR(quat_norm(result), 1.0f, 1e-6f);
}

static void test_angle_between_ignores_quaternion_sign(void)
{
    /* q and -q are the same orientation. Reporting ~180 degrees for that pair
     * is a classic bug, and it corrupts any error metric built on top. */
    quat q = quat_from_euler(0.2f, 0.3f, 0.4f);
    quat negated = quat_make(-q.w, -q.x, -q.y, -q.z);
    CHECK_NEAR(quat_angle_between(q, negated), 0.0f, 1e-5f);
}

void register_quaternion_tests(void)
{
    RUN(test_identity_is_unit_and_neutral);
    RUN(test_conjugate_undoes_a_rotation);
    RUN(test_rotation_preserves_vector_length);
    RUN(test_inverse_rotation_returns_the_original_vector);
    RUN(test_ninety_degree_yaw_maps_x_to_y);
    RUN(test_euler_round_trip);
    RUN(test_euler_survives_near_vertical_pitch);
    RUN(test_rotation_vector_matches_a_known_rotation);
    RUN(test_rotation_vector_stays_on_the_unit_sphere);
    RUN(test_integrating_many_small_steps_does_not_drift);
    RUN(test_normalize_rejects_a_degenerate_quaternion);
    RUN(test_angle_between_ignores_quaternion_sign);
}
