/* OSQP TESTER MODULE */

/* THE CODE FOR MINIMAL UNIT TESTING HAS BEEN TAKEN FROM
   http://www.jera.com/techinfo/jtns/jtn002.html */

#include <stdio.h>

#include "minunit.h"
#include "osqp.h"
#include "glob_opts.h"
#include "osqp_tester.h"

// Include tests
#include "utils/c_test_utils.h" //helper functions
#include "svm/test_svm.h"

int tests_run = 0;

static const char* all_tests() {
  mu_run_test(test_svm);
  return 0;
}

int main(void) {
  const char *result = all_tests();

  if (result != 0) {
    printf("%s\n", result);
  }
  else {
    printf("i'm done\n");
  }
  return result != 0;
}
