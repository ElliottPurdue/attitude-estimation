/* A minimal test harness.
 *
 * No framework, because the library has no dependencies and a test suite that
 * needs one would undercut that. It reports every failure with a file and line
 * rather than stopping at the first, so a single run shows the full damage.
 */

#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <math.h>

extern int tests_run;
extern int checks_failed;
extern const char *current_test;

#define RUN(fn)                     \
    do {                            \
        current_test = #fn;         \
        tests_run++;                \
        fn();                       \
    } while (0)

#define CHECK(condition)                                              \
    do {                                                              \
        if (!(condition)) {                                           \
            checks_failed++;                                          \
            printf("  FAIL  %s\n        %s:%d  %s\n",                 \
                   current_test, __FILE__, __LINE__, #condition);     \
        }                                                             \
    } while (0)

#define CHECK_NEAR(actual, expected, tolerance)                            \
    do {                                                                   \
        double a_ = (double)(actual), e_ = (double)(expected);             \
        if (!(fabs(a_ - e_) <= (double)(tolerance))) {                     \
            checks_failed++;                                               \
            printf("  FAIL  %s\n        %s:%d  %s\n"                       \
                   "        expected %.9g, got %.9g (tolerance %.3g)\n",   \
                   current_test, __FILE__, __LINE__, #actual,              \
                   e_, a_, (double)(tolerance));                           \
        }                                                                  \
    } while (0)

#define CHECK_BELOW(actual, limit)                                         \
    do {                                                                   \
        double a_ = (double)(actual), l_ = (double)(limit);                \
        if (!(a_ < l_)) {                                                  \
            checks_failed++;                                               \
            printf("  FAIL  %s\n        %s:%d  %s\n"                       \
                   "        expected below %.9g, got %.9g\n",              \
                   current_test, __FILE__, __LINE__, #actual, l_, a_);     \
        }                                                                  \
    } while (0)

#endif /* TEST_H */
