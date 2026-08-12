/* Measures filter update cost on an ESP32-S3.
 *
 * The host benchmark reports what the algorithms cost on a desktop CPU, which
 * settles nothing about whether they fit in a control loop on a microcontroller.
 * A 240 MHz Xtensa LX7 with a single-precision FPU is a different machine, and
 * the only honest way to get its number is to run it there.
 *
 * SETUP
 *   Copy the library into this folder -- the Arduino build compiles every
 *   source file beside the sketch, and will not reach outside it:
 *
 *     copy ..\..\src\*.c  .
 *     copy ..\..\src\*.h  .
 *
 *   Then select "ESP32S3 Dev Module", flash, and open the serial monitor at
 *   115200. Results print once and then idle.
 *
 * WHAT TO EXPECT
 *   Both filters are allocation-free and fixed-cost per update: no loops over
 *   data, no branches that depend on magnitude beyond the accelerometer gate.
 *   The measurement should be stable across runs, and a large variance would
 *   itself be worth investigating.
 */

extern "C" {
#include "complementary.h"
#include "ekf.h"
}

static const int ITERATIONS = 20000;
static const float DT = 0.005f;

/* Deterministic inputs, matching the host benchmark so the two are comparable.
 * Varying them matters: a repeated sample lets the branch predictor and the
 * accelerometer trust gate settle into one path, flattering whichever filter
 * branches more. */
static unsigned rng_state = 1u;

static float noise(float scale)
{
    rng_state = rng_state * 1103515245u + 12345u;
    return scale * ((float)((rng_state >> 16) & 0x7fff) / 16383.5f - 1.0f);
}

static void report(const char *name, unsigned long micros_elapsed,
                   unsigned state_bytes)
{
    float per_update = (float)micros_elapsed / (float)ITERATIONS;
    Serial.printf("  %-16s %8.3f us   %8.0f Hz max   %4u B state\n",
                  name, per_update,
                  per_update > 0.0f ? 1e6f / per_update : 0.0f,
                  state_bytes);
}

void setup()
{
    Serial.begin(115200);
    delay(2000);

    Serial.println();
    Serial.printf("ESP32-S3 attitude filter benchmark\n");
    Serial.printf("  CPU %u MHz, %d iterations, dt=%.3f s\n\n",
                  (unsigned)getCpuFrequencyMhz(), ITERATIONS, DT);

    complementary_filter comp;
    complementary_init(&comp, COMPLEMENTARY_DEFAULT_KP, COMPLEMENTARY_DEFAULT_KI);

    ekf_filter ekf;
    ekf_init(&ekf);

    /* One update each first, so neither pays a cold-cache cost the other avoids. */
    complementary_update(&comp, vec3_make(0, 0, 0), vec3_make(0, 0, 1), DT);
    ekf_update(&ekf, vec3_make(0, 0, 0), vec3_make(0, 0, 1), DT);

    Serial.printf("  %-16s %11s %14s %12s\n",
                  "filter", "us/update", "max rate", "state");
    Serial.println("  ------------------------------------------------------");

    rng_state = 1u;
    unsigned long start = micros();
    for (int i = 0; i < ITERATIONS; ++i) {
        vec3 gyro = vec3_make(noise(0.5f), noise(0.5f), noise(0.5f));
        vec3 accel = vec3_make(noise(0.05f), noise(0.05f), 1.0f + noise(0.05f));
        complementary_update(&comp, gyro, accel, DT);
    }
    report("complementary", micros() - start, sizeof(complementary_filter));

    rng_state = 1u;
    start = micros();
    for (int i = 0; i < ITERATIONS; ++i) {
        vec3 gyro = vec3_make(noise(0.5f), noise(0.5f), noise(0.5f));
        vec3 accel = vec3_make(noise(0.05f), noise(0.05f), 1.0f + noise(0.05f));
        ekf_update(&ekf, gyro, accel, DT);
    }
    report("ekf", micros() - start, sizeof(ekf_filter));

    /* Printed so the optimiser cannot discard the loops as dead code. */
    Serial.printf("\n  (final states: %.6f  %.6f)\n",
                  comp.orientation.w, ekf.orientation.w);
    Serial.printf("  free heap: %u bytes\n", (unsigned)ESP.getFreeHeap());
}

void loop()
{
    delay(10000);
}
