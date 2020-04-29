#ifndef UNITTEST_H
#define UNITTEST_H
/* unit tests  */

#define ut_assert(message, test) \
  do { if (!(test)) return message; } while (0)
#define ut_run_test(test)                   \
  do { const char *message = test(); tests_run++; \
       if (message) return message; } while (0)
extern int tests_run;

#endif //UNITTEST_H
