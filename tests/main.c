#include "test.h"

int tests_run = 0;
int checks_failed = 0;
const char *current_test = "(none)";

void register_quaternion_tests(void);
void register_complementary_tests(void);
void register_ekf_tests(void);

int main(void)
{
    printf("attitude-estimation test suite\n\n");

    register_quaternion_tests();
    register_complementary_tests();
    register_ekf_tests();

    printf("\n%d tests, %d failed check%s\n",
           tests_run, checks_failed, checks_failed == 1 ? "" : "s");
    return checks_failed == 0 ? 0 : 1;
}
