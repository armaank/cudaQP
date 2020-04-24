#include <stdio.h>
#include "../../../include/unittest.h"
#include "test_linsys.h"

int tests_run = 0;


static const char* all_tests() {

    ut_run_test(test_solve_linsys);

    return 0;
}

int main(void) {
    const char *result = all_tests();

    if (result != 0)
    {
        printf("%s\n", result);
    }
    else
    {
        printf("tests passed\n");
    }

    printf("tests run: %d\n", tests_run);

    return result != 0;
}
