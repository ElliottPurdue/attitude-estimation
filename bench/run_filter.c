/* Drives a filter over a recorded IMU log and writes its estimates.
 *
 * Kept separate from the library so that src/ has no file I/O and no stdio
 * dependency at all -- the same objects compile for a microcontroller
 * unchanged. This program exists only to connect the C filters to the Python
 * harness that generates truth and scores the result.
 *
 * Reads the CSV written by sim/generate.py on stdin or from a path, and writes
 * t plus the estimated quaternion and bias to stdout.
 *
 * Filters: "complementary", "ekf" (6-DOF), "ekf9" (adds the magnetometer).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/complementary.h"
#include "../src/ekf.h"

#define MAX_LINE 512

typedef struct {
    float t;
    vec3 gyro;
    vec3 accel;
    vec3 mag;
    int has_mag;
    quat truth;
} sample;

/* Accepts both log formats: the 14-column one with magnetometer readings, and
 * the 11-column one without. Logs recorded before the magnetometer existed stay
 * readable, and a 6-DOF run over a 9-DOF log is still a fair comparison because
 * the extra columns are simply ignored. */
static int parse_line(const char *line, sample *out)
{
    int fields = sscanf(line, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
                        &out->t,
                        &out->gyro.x, &out->gyro.y, &out->gyro.z,
                        &out->accel.x, &out->accel.y, &out->accel.z,
                        &out->mag.x, &out->mag.y, &out->mag.z,
                        &out->truth.w, &out->truth.x, &out->truth.y,
                        &out->truth.z);

    if (fields == 14) {
        out->has_mag = 1;
        return 1;
    }

    if (fields == 11) {
        /* Without magnetometer columns the quaternion landed one slot early:
         * its four components are in mag.xyz and truth.w. Shift them back,
         * reading truth.w before it is overwritten. */
        float qz = out->truth.w;
        out->truth = quat_make(out->mag.x, out->mag.y, out->mag.z, qz);
        out->mag = vec3_make(0.0f, 0.0f, 0.0f);
        out->has_mag = 0;
        return 1;
    }

    return 0;
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : NULL;
    const char *which = (argc > 2) ? argv[2] : "complementary";
    int use_ekf = (which[0] == 'e');
    int use_mag = (strcmp(which, "ekf9") == 0);

    FILE *input = path ? fopen(path, "r") : stdin;
    if (!input) {
        fprintf(stderr, "cannot open %s\n", path);
        return 1;
    }

    char line[MAX_LINE];
    if (!fgets(line, sizeof(line), input)) {   /* discard the header row */
        fprintf(stderr, "empty input\n");
        return 1;
    }

    complementary_filter comp;
    ekf_filter ekf;
    complementary_init(&comp, COMPLEMENTARY_DEFAULT_KP, COMPLEMENTARY_DEFAULT_KI);
    ekf_init(&ekf);

    printf("t,qw,qx,qy,qz,bx,by,bz\n");

    sample previous;
    int have_previous = 0;
    int initialised = 0;
    int warned_no_mag = 0;
    long count = 0;

    while (fgets(line, sizeof(line), input)) {
        sample current;
        if (!parse_line(line, &current)) {
            continue;
        }

        if (use_mag && !current.has_mag && !warned_no_mag) {
            fprintf(stderr, "warning: log has no magnetometer columns; "
                            "ekf9 will behave as ekf\n");
            warned_no_mag = 1;
        }

        /* Seeded from the first accelerometer reading rather than started at
         * identity, so the run measures steady-state accuracy instead of a
         * convergence transient that would dominate a short log. */
        if (!initialised) {
            if (use_ekf) {
                ekf_set_from_accel(&ekf, current.accel);
            } else {
                complementary_set_from_accel(&comp, current.accel);
            }
            initialised = 1;
        }

        if (have_previous) {
            float dt = current.t - previous.t;
            if (use_ekf) {
                ekf_update(&ekf, current.gyro, current.accel, dt);
                if (use_mag && current.has_mag) {
                    ekf_update_magnetometer(&ekf, current.mag);
                }
            } else {
                complementary_update(&comp, current.gyro, current.accel, dt);
            }
        }

        quat q = use_ekf ? ekf.orientation : comp.orientation;
        vec3 b = use_ekf ? ekf.gyro_bias : comp.gyro_bias;

        printf("%.6f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
               current.t, q.w, q.x, q.y, q.z, b.x, b.y, b.z);

        previous = current;
        have_previous = 1;
        count++;
    }

    if (path) {
        fclose(input);
    }

    fprintf(stderr, "processed %ld samples with the %s filter\n", count, which);
    return 0;
}
