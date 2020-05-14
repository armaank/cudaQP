"""
generate random svm problems
"""

import numpy as np
import scipy as sp

from scipy import sparse
import utils.codegen_utils as cu

sp.random.seed(1234)
n = 50
m = 100
N = int(m / 2)
gamma = 1.0
b = np.hstack([np.ones(N), -np.ones(N)])
A_upp = sparse.random(N, n, density=0.5)
A_low = sparse.random(N, n, density=0.5)
Ad = sparse.vstack(
    [
        A_upp / np.sqrt(n) + (A_upp != 0.0).astype(float) / n,
        A_low / np.sqrt(n) - (A_low != 0.0).astype(float) / n,
    ],
    format="csc",
)

# osqp data
Im = sparse.eye(m)
P = sparse.block_diag([sparse.eye(n), sparse.csc_matrix((m, m))], format="csc")
q = np.hstack([np.zeros(n), gamma * np.ones(m)])
A = sparse.vstack(
    [
        sparse.hstack([sparse.diags(b).dot(Ad), -Im]),
        sparse.hstack([sparse.csc_matrix((m, n)), Im]),
    ],
    format="csc",
)
l = np.hstack([-np.inf * np.ones(m), np.zeros(m)])
u = np.hstack([-np.ones(m), np.inf * np.ones(m)])

cu.generate_problem_data(P, q, A, l, u, "svm")
