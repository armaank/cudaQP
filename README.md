# cudasvm
gpu accelerated solver for convex quadratic programs with benchmarks for svm classificaiton problems 

This is essentially a soft fork of ![osqp](https://github.com/oxfordcontrol/osqp) that has additional support for a CUDA `lin_sys` solver and optimizations for a Nvidia Jetson Nano target device. 

You can find my full report ![here](docs/report/report.pdf).

## Installation Requirements
To run the project on a Jetson Nano you must have access to the following:
```
# for problem generation
numpy
scipy
# CUDA libs
cuBLAS
cuSparse
# other
CMake

```
`cuBLAS` and `cuSPARSE` should be there, but you may need to install `cmake`, `numpy` and `scipy`, and possibly `gFortran` for `scipy`

## Basic Usage
To run an SVM classifier
```
# generate Make files
cmake CMakeLists.txt
# run Make
make
# run svm program
./out/svm_benchmark
```
To edit the size of the svm problem, you can edit it by modifying:
```
./benchmarks/svm/generate_problem.py
```
and then re-run `make` at the project root directory

In order to switch between using the LDL(cpu) and PCG(gpu) solver, you can directly edit `CMakeLists.txt` at the root directory to toggle the `Enable CUDA support` option. Then, re-generate the make files with `cmake` and run your problem. 

## References
[1] B. Stellato, G. Banjac, P. Goulart, A. Bemporad, and S. Boyd, “OSQP: An operator splitting
solver for quadratic programs,” Mathematical Programming Computation, 2020. [Online]. Available:
https://doi.org/10.1007/s12532-020-00179-2

[2] T. Hastie, R. Tibshirani, and J. Friedman, The Elements of Statistical Learning, ser. Springer Series in Statistics.
New York, NY, USA: Springer New York Inc., 2001.

[3] M. Schubiger, G. Banjac, and J. Lygeros, “GPU acceleration of ADMM for large-scale quadratic program-
ming,” arXiv:1912.04263, 2019.

