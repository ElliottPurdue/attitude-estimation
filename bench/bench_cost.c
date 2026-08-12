/* Measures the cost of one filter update.
 *
 * The number that decides whether an estimator is usable on a given
 * microcontroller is microseconds per update, because it sets the fastest loop
 * the processor can close. An estimator that is 3x more accurate and does not
 * fit in the control period is not the better choice.
 *
 * Builds for the host or, cross-compiled, for a target. Reports cost per update
 * and the maximum rate that leaves the CPU otherwise idle.
 *
 * Inputs vary between iterations. Feeding the same sample repeatedly lets the
 * branch predictor and the accelerometer trust gate settle into one path, which
 * flatters whichever filter branches more.
 */

#include <stdio.h>
#include <time.h>

#include "../src/complementary.h"
#include "../src/ekf.h"

#define ITERATIONS 200000
#define DT 0.005f

/* Deterministic pseudo-random input, so host and target runs are comparable and
 * neither depends on a library the target may lack. */
static unsigned rng_state = 1u;

static float noise(float scale)
{
    rng_state = rng_state * 1103515245u + 12345u;
    return scale * ((float)((rng_state >> 16) & 0x7fff) / 16383.5f - 1.0f);
}

static double seconds_now(void)
{
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

int main(void)
{
    printf("attitude filter update cost, %d iterations at dt=%.3f s\n\n",
           ITERATIONS, DT);

    complementary_filter comp;
    complementary_init(&comp, COMPLEMENTARY_DEFAULT_KP, COMPLEMENTARY_DEFAULT_KI);

    ekf_filter ekf;
    ekf_init(&ekf);

    /* Touch both so neither pays a cold-cache penalty the other avoids. */
    complementary_update(&comp, vec3_make(0, 0, 0), vec3_make(0, 0, 1), DT);
    ekf_update(&ekf, vec3_make(0, 0, 0), vec3_make(0, 0, 1), DT);

    rng_state = 1u;
    double start = seconds_now();
    for (int i = 0; i < ITERATIONS; ++i) {
        vec3 gyro = vec3_make(noise(0.5f), noise(0.5f), noise(0.5f));
        vec3 accel = vec3_make(noise(0.05f), noise(0.05f), 1.0f + noise(0.05f));
        complementary_update(&comp, gyro, accel, DT);
    }
    double comp_seconds = seconds_now() - start;

    rng_state = 1u;
    start = seconds_now();
    for (int i = 0; i < ITERATIONS; ++i) {
        vec3 gyro = vec3_make(noise(0.5f), noise(0.5f), noise(0.5f));
        vec3 accel = vec3_make(noise(0.05f), noise(0.05f), 1.0f + noise(0.05f));
        ekf_update(&ekf, gyro, accel, DT);
    }
    double ekf_seconds = seconds_now() - start;

    /* The same loop plus a magnetometer update, so the difference between the
     * two rows is the cost of heading correction rather than a separate
     * measurement of it. */
    rng_state = 1u;
    start = seconds_now();
    for (int i = 0; i < ITERATIONS; ++i) {
        vec3 gyro = vec3_make(noise(0.5f), noise(0.5f), noise(0.5f));
        vec3 accel = vec3_make(noise(0.05f), noise(0.05f), 1.0f + noise(0.05f));
        vec3 mag = vec3_make(0.5f + noise(0.02f), noise(0.02f),
                             -0.87f + noise(0.02f));
        ekf_update(&ekf, gyro, accel, DT);
        ekf_update_magnetometer(&ekf, mag);
    }
    double ekf9_seconds = seconds_now() - start;

    double comp_us = comp_seconds * 1e6 / ITERATIONS;
    double ekf_us = ekf_seconds * 1e6 / ITERATIONS;
    double ekf9_us = ekf9_seconds * 1e6 / ITERATIONS;

    printf("  %-16s %10s %14s %12s\n", "filter", "us/update", "max rate", "state bytes");
    printf("  %s\n", "------------------------------------------------------");
    printf("  %-16s %10.3f %11.0f Hz %12u\n", "complementary", comp_us,
           comp_us > 0.0 ? 1e6 / comp_us : 0.0,
           (unsigned)sizeof(complementary_filter));
    printf("  %-16s %10.3f %11.0f Hz %12u\n", "ekf", ekf_us,
           ekf_us > 0.0 ? 1e6 / ekf_us : 0.0,
           (unsigned)sizeof(ekf_filter));
    printf("  %-16s %10.3f %11.0f Hz %12u\n", "ekf + mag", ekf9_us,
           ekf9_us > 0.0 ? 1e6 / ekf9_us : 0.0,
           (unsigned)sizeof(ekf_filter));

    if (comp_us > 0.0) {
        printf("\n  the EKF costs %.1fx the complementary filter per update\n",
               ekf_us / comp_us);
    }
    if (ekf_us > 0.0) {
        printf("  heading correction adds %.3f us, %.0f%% of the accel update\n",
               ekf9_us - ekf_us, 100.0 * (ekf9_us - ekf_us) / ekf_us);
    }

    /* Prevents the optimiser discarding the loops as dead code. */
    printf("  (final states: %.6f %.6f)\n",
           (double)comp.orientation.w, (double)ekf.orientation.w);
    return 0;
}
