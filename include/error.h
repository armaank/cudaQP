#ifndef ERROR_H
# define ERROR_H

/* QP error handling */




#  define osqp_error(error_code) _qp_error(error_code, __func__);
/**
 * Internal function to print error description and return error code.
 * @param  Error code
 * @param  Function name
 * @return Error code
 */
  int _qp_error(enum qp_error_type error_code,
		    const char * function_name);





#endif // ifndef ERROR_H