#include "quaternion.h"

#include <math.h>

quat quat_identity(void)
{
    return quat_make(1.0f, 0.0f, 0.0f, 0.0f);
}

quat quat_make(float w, float x, float y, float z)
{
    quat q;
    q.w = w; q.x = x; q.y = y; q.z = z;
    return q;
}

quat quat_multiply(quat a, quat b)
{
    return quat_make(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w);
}

quat quat_conjugate(quat q)
{
    return quat_make(q.w, -q.x, -q.y, -q.z);
}

float quat_norm(quat q)
{
    return sqrtf(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
}

quat quat_normalize(quat q)
{
    float n = quat_norm(q);
    if (n < 1e-9f) {
        return quat_identity();
    }
    float inv = 1.0f / n;
    return quat_make(q.w * inv, q.x * inv, q.y * inv, q.z * inv);
}

vec3 quat_rotate(quat q, vec3 v)
{
    /* v' = v + 2 * cross(u, cross(u, v) + w*v), with u the vector part.
     *
     * Algebraically identical to q * v * q^-1 but avoids building the full
     * quaternion product: roughly 15 multiplies against 28.
     */
    vec3 u = vec3_make(q.x, q.y, q.z);
    vec3 t = vec3_add(vec3_cross(u, v), vec3_scale(v, q.w));
    return vec3_add(v, vec3_scale(vec3_cross(u, t), 2.0f));
}

vec3 quat_rotate_inverse(quat q, vec3 v)
{
    return quat_rotate(quat_conjugate(q), v);
}

quat quat_from_rotation_vector(vec3 rotation)
{
    float angle = vec3_norm(rotation);

    /* Below this the small-angle series is more accurate than the trigonometric
     * form, because sin(angle)/angle loses precision as angle approaches zero.
     * At 1 kHz this branch is the common case: a 1 deg/s rate gives an angle of
     * about 1.7e-5 rad per step. */
    if (angle < 1e-7f) {
        return quat_normalize(quat_make(1.0f,
                                        0.5f * rotation.x,
                                        0.5f * rotation.y,
                                        0.5f * rotation.z));
    }

    float half = 0.5f * angle;
    float scale = sinf(half) / angle;
    return quat_make(cosf(half),
                     rotation.x * scale,
                     rotation.y * scale,
                     rotation.z * scale);
}

vec3 quat_to_euler(quat q)
{
    vec3 euler;

    float sinr_cosp = 2.0f * (q.w * q.x + q.y * q.z);
    float cosr_cosp = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);
    euler.x = atan2f(sinr_cosp, cosr_cosp);

    /* Clamped because accumulated rounding can push this just past +/-1, and
     * asinf would then return NaN at exactly the attitudes -- straight up or
     * straight down -- where a vehicle most needs an estimate. */
    float sinp = 2.0f * (q.w * q.y - q.z * q.x);
    if (sinp > 1.0f) {
        sinp = 1.0f;
    } else if (sinp < -1.0f) {
        sinp = -1.0f;
    }
    euler.y = asinf(sinp);

    float siny_cosp = 2.0f * (q.w * q.z + q.x * q.y);
    float cosy_cosp = 1.0f - 2.0f * (q.y * q.y + q.z * q.z);
    euler.z = atan2f(siny_cosp, cosy_cosp);

    return euler;
}

quat quat_from_euler(float roll, float pitch, float yaw)
{
    float cr = cosf(roll * 0.5f),  sr = sinf(roll * 0.5f);
    float cp = cosf(pitch * 0.5f), sp = sinf(pitch * 0.5f);
    float cy = cosf(yaw * 0.5f),   sy = sinf(yaw * 0.5f);

    return quat_make(cr * cp * cy + sr * sp * sy,
                     sr * cp * cy - cr * sp * sy,
                     cr * sp * cy + sr * cp * sy,
                     cr * cp * sy - sr * sp * cy);
}

float quat_angle_between(quat a, quat b)
{
    quat error = quat_multiply(quat_conjugate(a), b);

    /* Computed as 2*atan2(|vector part|, |scalar part|) rather than the more
     * obvious 2*acos(w).
     *
     * acos is ill-conditioned exactly where this function matters most. Near
     * zero rotation w approaches 1, and acos(1-eps) ~ sqrt(2*eps), so a float's
     * 6e-8 of rounding becomes 7e-4 rad -- about 0.04 degrees of error floor on
     * a metric whose whole purpose is measuring small errors. atan2 stays
     * well-conditioned across the range and reports zero for identical inputs.
     *
     * q and -q are the same rotation, so the absolute value on the scalar part
     * picks the shorter arc and keeps the result in [0, pi]. Without it a small
     * rotation carrying a negative scalar part reads as nearly 180 degrees.
     */
    vec3 axis = vec3_make(error.x, error.y, error.z);
    return 2.0f * atan2f(vec3_norm(axis), fabsf(error.w));
}
