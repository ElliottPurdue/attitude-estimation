/* Unit quaternions for attitude.
 *
 * Convention: Hamilton, scalar-first (w, x, y, z), representing the rotation
 * that takes a vector from the body frame to the world frame. Every function
 * here assumes that convention; mixing it with the JPL convention used in parts
 * of the aerospace literature silently inverts rotations, so it is stated
 * rather than implied.
 *
 * Quaternions are used instead of Euler angles because they have no singularity
 * at 90 degrees of pitch, and instead of rotation matrices because they carry
 * four numbers rather than nine and renormalize cheaply.
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#include "vec3.h"

typedef struct {
    float w, x, y, z;
} quat;

quat quat_identity(void);
quat quat_make(float w, float x, float y, float z);

/* Hamilton product: the rotation b followed by the rotation a. */
quat quat_multiply(quat a, quat b);

/* Conjugate, which for a unit quaternion is also the inverse. */
quat quat_conjugate(quat q);

float quat_norm(quat q);

/* Scales to unit length. Returns identity for a degenerate quaternion rather
 * than dividing by zero. */
quat quat_normalize(quat q);

/* Rotates a body-frame vector into the world frame. */
vec3 quat_rotate(quat q, vec3 v);

/* Rotates a world-frame vector into the body frame. */
vec3 quat_rotate_inverse(quat q, vec3 v);

/* Small-angle update from a rotation vector, in radians.
 *
 * Uses the exact exponential map rather than the first-order approximation
 * q + 0.5*omega*q. The approximation drifts off the unit sphere and needs
 * renormalizing every step, and its error grows with rotation rate, which is
 * exactly when an estimator can least afford it.
 */
quat quat_from_rotation_vector(vec3 rotation);

/* Roll, pitch, yaw in radians, in the aerospace Z-Y-X sequence. Provided for
 * reporting and tests; the filters never convert to Euler internally. */
vec3 quat_to_euler(quat q);
quat quat_from_euler(float roll, float pitch, float yaw);

/* Geodesic angle between two orientations, in radians.
 *
 * This is the correct error metric for attitude: it is the single rotation
 * angle separating the two frames, invariant to how they are parameterized.
 * Comparing Euler angles component-wise instead exaggerates error near
 * gimbal lock and is not a distance at all.
 */
float quat_angle_between(quat a, quat b);

#endif /* QUATERNION_H */
