/* Three-element vectors and the small amount of linear algebra the filters need.
 *
 * Everything is float and by value. On a Cortex-M or Xtensa core with a
 * single-precision FPU, doubles are emulated in software and cost an order of
 * magnitude more per operation, which is the difference between an estimator
 * that closes a 1 kHz loop and one that does not.
 */

#ifndef VEC3_H
#define VEC3_H

#include <math.h>

typedef struct {
    float x, y, z;
} vec3;

static inline vec3 vec3_make(float x, float y, float z)
{
    vec3 v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

static inline vec3 vec3_add(vec3 a, vec3 b)
{
    return vec3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline vec3 vec3_sub(vec3 a, vec3 b)
{
    return vec3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline vec3 vec3_scale(vec3 a, float k)
{
    return vec3_make(a.x * k, a.y * k, a.z * k);
}

static inline float vec3_dot(vec3 a, vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline vec3 vec3_cross(vec3 a, vec3 b)
{
    return vec3_make(a.y * b.z - a.z * b.y,
                     a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x);
}

static inline float vec3_norm(vec3 a)
{
    return sqrtf(vec3_dot(a, a));
}

/* Returns the zero vector rather than dividing by a near-zero norm.
 *
 * A free-falling or violently shaken IMU reports an accelerometer magnitude at
 * or near zero. Callers must treat the zero vector as "no usable measurement"
 * and skip the correction step; normalizing it anyway would inject NaN into the
 * state and every subsequent estimate.
 */
static inline vec3 vec3_normalize(vec3 a)
{
    float n = vec3_norm(a);
    if (n < 1e-9f) {
        return vec3_make(0.0f, 0.0f, 0.0f);
    }
    return vec3_scale(a, 1.0f / n);
}

#endif /* VEC3_H */
