#ifndef UTIL_H
# define UTIL_H

# include "types.h"
# include "constants.h"


/**********************
* Utility Functions  *
**********************/


/**
 * Copy settings creating a new settings structure (uses MALLOC)
 * @param  settings Settings to be copied
 * @return          New settings structure
 */
OSQPSettings* copy_settings(const OSQPSettings *settings);


/**
 * Custom string copy to avoid string.h library
 * @param dest   destination string
 * @param source source string
 */
void c_strcpy(char       dest[],
              const char source[]);



/**
 * Print Header before running the algorithm
 * @param work     osqp workspace
 */
void print_setup_header(const OSQPWorkspace *work);

/**
 * Print header with data to be displayed per iteration
 */
void print_header(void);

/**
 * Print iteration summary
 * @param work current workspace
 */
void print_summary(OSQPWorkspace *work);

/**
 * Print information after polish
 * @param work current workspace
 */
void print_polish(OSQPWorkspace *work);

/**
 * Print footer when algorithm terminates
 * @param info   info structure
 * @param polish is polish enabled?
 */
void print_footer(OSQPInfo *info,
                  int     polish);




/*********************************
* Timer Structs and Functions * *
*********************************/

/*! \cond PRIVATE */




/* Use POSIX clock_gettime() for timing on non-Windows machines */
#   include <time.h>
#   include <sys/time.h>


struct OSQP_TIMER {
  struct timespec tic;
  struct timespec toc;
};


/*! \endcond */

/**
 * Timer Methods
 */

/**
 * Start timer
 * @param t Timer object
 */
void    osqp_tic(OSQPTimer *t);

/**
 * Report time
 * @param  t Timer object
 * @return   Reported time
 */
float osqp_toc(OSQPTimer *t);



/* ================================= DEBUG FUNCTIONS ======================= */

/*! \cond PRIVATE */



/* Compare CSC matrices */
int is_eq_csc(csc    *A,
                csc    *B,
                float tol);

/* Convert sparse CSC to dense */
float* csc_to_dns(csc *M);



#  include <stdio.h>


/* Print a csc sparse matrix */
void print_csc_matrix(csc        *M,
                      const char *name);

/* Dump csc sparse matrix to file */
void dump_csc_matrix(csc        *M,
                     const char *file_name);

/* Print a triplet format sparse matrix */
void print_trip_matrix(csc        *M,
                       const char *name);

/* Print a dense matrix */
void print_dns_matrix(float    *M,
                      int       m,
                      int       n,
                      const char *name);

/* Print vector  */
void print_vec(float    *v,
               int       n,
               const char *name);

/* Dump vector to file */
void dump_vec(float    *v,
              int       len,
              const char *file_name);

// Print int array
void print_veint(int      *x,
                   int       n,
                   const char *name);


/*! \endcond */


#endif // ifndef UTIL_H
