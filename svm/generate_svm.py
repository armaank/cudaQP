"""generate_svm

generate random svm problems
"""
import os

import numpy as np
import scipy as sp

from scipy import sparse

import translate 

def make_svm(n, m, name):
    """
    generate random svm problem
    """
    sp.random.seed(1)
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
    translate.construct_svm_qp(P, q, A, l, u, name)   

    pass

def svm_example():
    """
    toy svm problem
    """
    sp.random.seed(1)
    n = 4 # 10
    m = 4 # 1000
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
    translate.construct_svm_qp(P, q, A, l, u, "svm_example")   

def main():
    n = 10
    m = 1000
    name = "svm_example"
    make_svm(n, m, name)
    pass


if __name__=="__main__":

    # add argparsing to pass to main
    main() 
    
    pass
