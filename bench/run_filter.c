/* Drives a filter over a recorded IMU log and writes its estimates.
 *
 * Kept separate from the library so that src/ has no file I/O and no stdio
 * dependency at all -- the same objects compile for a microcontroller
 * unchanged. This program exists only to connect the C filters to the Python
 * harness that generates truth and scores the result.
 *
 * Reads the CSV written by sim/generate.py on stdin or from a path, and writes
 * t plus the estimated quaternion and bias to stdout.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/complementary.h"

#define MAX_LINE 512

typedef struct {
    float t;
    vec3 gyro;
    vec3 accel;
    quat truth;
} sample;

static int parse_line(const char *line, sample *out)
{
    return sscanf(line, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
                  &out->t,
                  &out->gyro.x, &out->gyro.y, &out->gyro.z,
                  &out->accel.x, &out->accel.y, &out->accel.z,
                  &out->truth.w, &out->truth.x, &out->truth.y,
                  &out->truth.z) == 11;
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : NULL;
    float kp = (argc > 2) ? (float)atof(argv[2]) : COMPLEMENTARY_DEFAULT_KP;
    float ki = (argc > 3) ? (float)atof(argv[3]) : COMPLEMENTARY_DEFAULT_KI;

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

    complementary_filter filter;
    complementary_init(&filter, kp, ki);

    printf("t,qw,qx,qy,qz,bx,by,bz\n");

    sample previous;
    int have_previous = 0;
    int initialised = 0;
    long count = 0;

    while (fgets(line, sizeof(line), input)) {
        sample current;
        if (!parse_line(line, &current)) {
            continue;
        }

        /* Seeded from the first accelerometer reading rather than started at
         * identity, so the run measures steady-state accuracy instead of a
         * convergence transient that would dominate a short log. */
        if (!initialised) {
            complementary_set_from_accel(&filter, current.accel);
            initialised = 1;
        }

        if (have_previous) {
            float dt = current.t - previous.t;
            complementary_update(&filter, current.gyro, current.accel, dt);
        }

        printf("%.6f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
               current.t,
               filter.orientation.w, filter.orientation.x,
               filter.orientation.y, filter.orientation.z,
               filter.gyro_bias.x, filter.gyro_bias.y, filter.gyro_bias.z);

        previous = current;
        have_previous = 1;
        count++;
    }

    if (path) {
        fclose(input);
    }

    fprintf(stderr, "processed %ld samples (kp=%g ki=%g)\n", count, kp, ki);
    return 0;
}
