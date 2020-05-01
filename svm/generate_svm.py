"""generate_svm

generate random svm problems
"""
import argparse
import os

import numpy as np
import scipy as sp

from scipy import sparse

import translate 

def make_svm(n, m, seed, name):
    """
    generate random svm problem
    """

    sp.random.seed(seed)
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

def main(args):
    """
    generate svm problem
    """

    make_svm(args.n, args.m, args.rand_seed, args.name)
    
    pass


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="generate svm problem")
    
    parser.add_argument("--n", default=10, type=int, help="problem dimension")
    parser.add_argument("--m", default=1000, type=int, help="number of data points")
    parser.add_argument("--rand_seed", default=1, type=int, help="random seed")
    parser.add_argument("--name", default="m1000", type=str, help="experiment name")

    args = parser.parse_args()

    main(args)

    pass
