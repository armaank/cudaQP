import numpy as np
from scipy import sparse
import utils.codegen_utils as cu

"""
P = sparse.triu([[11., 0.], [0., 0.]], format='csc')
q = np.array([3., 4.])

A = sparse.csc_matrix(np.array([[-1., 0.], [0., -1.], [-1., 3.],
                                [2., 5.], [3., 4]]))
l = -np.inf * np.ones(A.shape[0])
u = np.array([0., 0., -15., 100., 80.])

n = P.shape[0]
m = A.shape[0]

# New data
q_new = np.array([1., 1.])
u_new = np.array([-2., 0., -20., 100., 80.])

# Generate problem solutions
sols_data = {'x_test': np.array([15., -0.]),
             'y_test': np.array([0., 508., 168., 0., 0.]),
             'obj_value_test': 1282.5,
             'status_test': 'optimal',
             'q_new': q_new,
             'u_new': u_new,
             'x_test_new': np.array([20., -0.]),
             'y_test_new': np.array([0., 664., 221., 0., 0.]),
             'obj_value_test_new': 2220.0,
             'status_test_new': 'optimal'}


# Generate problem data
cu.generate_problem_data(P, q, A, l, u, 'basic_svm', sols_data)
"""
n = 10
m = 1000
N = int(m / 2)
gamma = 1.0
b = np.hstack([np.ones(N), -np.ones(N)])
A_upp = sparse.random(N, n, density=0.5)
A_low = sparse.random(N, n, density=0.5)
Ad = sparse.vstack([
        A_upp / np.sqrt(n) + (A_upp != 0.).astype(float) / n,
        A_low / np.sqrt(n) - (A_low != 0.).astype(float) / n
     ], format='csc')

# OSQP data
Im = sparse.eye(m)
P = sparse.block_diag([sparse.eye(n), sparse.csc_matrix((m, m))], format='csc')
q = np.hstack([np.zeros(n), gamma*np.ones(m)])
A = sparse.vstack([
        sparse.hstack([sparse.diags(b).dot(Ad), -Im]),
        sparse.hstack([sparse.csc_matrix((m, n)), Im])
    ], format='csc')
l = np.hstack([-np.inf*np.ones(m), np.zeros(m)])
u = np.hstack([-np.ones(m), np.inf*np.ones(m)])

cu.generate_problem_data(P, q, A, l, u, 'svm')

