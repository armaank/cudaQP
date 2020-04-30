"""generate_svm

generate random svm problems
"""
import os

import numpy as np

from scipy import sparse

def svm_example():
   
    # Generate problem data
    sp.random.seed(1)
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

    Im = sparse.eye(m)
    P = sparse.block_diag([sparse.eye(n), sparse.csc_matrix((m, m))], format='csc')
    q = np.hstack([np.zeros(n), gamma*np.ones(m)])
    A = sparse.vstack([
            sparse.hstack([sparse.diags(b).dot(Ad), -Im]),
            sparse.hstack([sparse.csc_matrix((m, n)), Im])
        ], format='csc')
    l = np.hstack([-np.inf*np.ones(m), np.zeros(m)])
    u = np.hstack([-np.ones(m), np.inf*np.ones(m)])

    # this is passed to codegen to generate a data.h file for svm.c program 


def main():


    pass


if __name__=="__main__":

    main()

    pass
