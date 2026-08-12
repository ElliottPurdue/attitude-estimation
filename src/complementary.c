#include "complementary.h"

#include <math.h>

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

void complementary_init(complementary_filter *filter, float kp, float ki)
{
    filter->orientation = quat_identity();
    filter->gyro_bias = vec3_make(0.0f, 0.0f, 0.0f);
    filter->kp = kp;
    filter->ki = ki;
    filter->accel_tolerance = COMPLEMENTARY_ACCEL_TOL;
    filter->max_bias = COMPLEMENTARY_MAX_BIAS;
    filter->gravity = COMPLEMENTARY_GRAVITY_G;
}

void complementary_set_from_accel(complementary_filter *filter, vec3 accel)
{
    vec3 measured = vec3_normalize(accel);
    if (vec3_norm(measured) < 0.5f) {
        return;                     /* unusable sample; leave the estimate alone */
    }

    /* At rest the accelerometer measures specific force, which points along +z
     * in the world frame -- opposite to gravity. The rotation carrying the
     * body-frame reading onto world +z is the attitude, up to yaw.
     *
     * Built from the axis-angle pair between the two vectors: the axis is their
     * cross product, the angle is the angle they subtend. Yaw is left at zero
     * because no rotation about vertical changes the reading. */
    vec3 reference = vec3_make(0.0f, 0.0f, 1.0f);
    vec3 axis = vec3_cross(measured, reference);
    float sine = vec3_norm(axis);
    float cosine = vec3_dot(measured, reference);

    if (sine < 1e-6f) {
        /* Already aligned, or antiparallel. Antiparallel means the device is
         * inverted, and any horizontal axis serves. */
        filter->orientation = (cosine > 0.0f)
            ? quat_identity()
            : quat_make(0.0f, 1.0f, 0.0f, 0.0f);
        return;
    }

    float angle = atan2f(sine, cosine);
    filter->orientation =
        quat_from_rotation_vector(vec3_scale(vec3_normalize(axis), angle));
}

void complementary_update(complementary_filter *filter,
                          vec3 gyro, vec3 accel, float dt)
{
    if (dt <= 0.0f) {
        return;
    }

    /* The accelerometer is informative only when measuring gravity alone. Its
     * magnitude is the test: near one g when the vehicle is unaccelerated, far
     * from it in free fall or under impact. Comparing the ratio against 1
     * rather than the raw difference keeps the gate independent of whether the
     * caller works in g's or m/s^2, provided filter->gravity matches. */
    float magnitude = vec3_norm(accel);
    int accel_trusted = 0;
    vec3 correction = vec3_make(0.0f, 0.0f, 0.0f);

    if (magnitude > 1e-6f && filter->gravity > 1e-6f) {
        float ratio = magnitude / filter->gravity;
        accel_trusted = fabsf(ratio - 1.0f) < filter->accel_tolerance;
    }

    if (accel_trusted) {
        vec3 measured = vec3_scale(accel, 1.0f / magnitude);

        /* Where the filter currently believes "up" lies, in the body frame. A
         * perfect estimate would make this equal the measurement. */
        vec3 expected = quat_rotate_inverse(filter->orientation,
                                            vec3_make(0.0f, 0.0f, 1.0f));

        /* Cross product of two unit vectors: its direction is the axis that
         * rotates the estimate onto the measurement, its magnitude the sine of
         * the angle between them. For small errors that is the error angle
         * itself, giving a usable proportional term with no trig call. */
        correction = vec3_cross(measured, expected);
    }

    vec3 rate = vec3_sub(gyro, filter->gyro_bias);

    /* The integral term accumulates persistent gravity error. A constant gyro
     * bias produces exactly that -- a steady tilt the proportional term keeps
     * fighting -- so this converges on the bias and cancels it at its source.
     *
     * Frozen whenever the accelerometer is untrusted. Integrating a zero
     * correction would be harmless, but integrating during a long acceleration
     * with a stale correction would poison the bias estimate, and the bias is
     * the one state that persists after the disturbance ends. */
    if (filter->ki > 0.0f && accel_trusted) {
        filter->gyro_bias = vec3_sub(
            filter->gyro_bias,
            vec3_scale(correction, filter->ki * dt));

        filter->gyro_bias = vec3_make(
            clampf(filter->gyro_bias.x, filter->max_bias),
            clampf(filter->gyro_bias.y, filter->max_bias),
            clampf(filter->gyro_bias.z, filter->max_bias));
    }

    rate = vec3_add(rate, vec3_scale(correction, filter->kp));

    filter->orientation = quat_multiply(
        filter->orientation,
        quat_from_rotation_vector(vec3_scale(rate, dt)));

    /* Renormalized every step. The exponential map preserves norm in exact
     * arithmetic, but float rounding accumulates over millions of updates. */
    filter->orientation = quat_normalize(filter->orientation);
}
