

#include "ctrl.h"
# include <signal.h>

static int int_detected;
struct sigaction oact;
static void handle_ctrlc(int dummy) {
    int_detected = dummy ? dummy : -1;
}

void qp_start_interrupt_listener(void) {
    struct sigaction act;

    int_detected = 0;
    act.sa_flags = 0;
    sigemptyset(&act.sa_mask);
    act.sa_handler = handle_ctrlc;
    sigaction(SIGINT, &act, &oact);
}

void qp_end_interrupt_listener(void) {
    struct sigaction act;

    sigaction(SIGINT, &oact, &act);
}

int qp_is_interrupted(void) {
    return int_detected;
}

