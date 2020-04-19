CC = gcc
CFLAGS  = -Wall

# make: test.c ldl.o
# 	$(CC) $(CFLAGS) test.c ldl_unittest.h ldl.o -o tests.out


# osqp


# lin sys

# ldl_interface
ldl_interface.o: ./src/scaling.c ldl.o kkt.o  ./include/glob_opts.h
	$(CC) $(CFLAGS) -c ./src/ldl_interface.c ldl.o kkt.o

# scaling - done?
scaling.o: ./src/scaling.c ./include/scaling.h ./include/lin_alg.h
	$(CC) $(CFLAGS) -c ./src/scaling.c 

# proj - done?
proj.o: ./src/proj.c ./include/proj.h
	$(CC) $(CFLAGS) -c ./src/proj.c

# lin alg - done?
lin_alg.o: ./src/lin_alg.c ./include/lin_alg.h ./include/glob_opts.h
	$(CC) $(CFLAGS) -c ./src/lin_alg.c 

# kkt - done?
kkt.o: ./src/kkt.c ./include/kkt.h
	$(CC) $(CFLAGS) -c ./src/kkt.c 

# error handling - done?
error.o: ./src/error.c ./include/error.h ./include/constants.h
	$(CC) $(CFLAGS) -c ./src/error.c 

# ctrl - done?
ctrlc.o: ./src/ctrlc.c ./include/ctrlc.h
	$(CC) $(CFLAGS) -c ./src/ctrlc.c 

# csc - done?
cs.o: ./src/cs.c ./include/cs.h
	$(CC) $(CFLAGS) -c ./src/cs.c 

# ldl lib - done?
ldl.o: ./ldl/ldl.c ./ldl/ldl.h
	$(CC) $(CFLAGS) -c ./ldl/ldl.c 

# util - done?
util.o: ./src/util.c ./include/util.h ./include/glob_opts.h
	$(CC) $(CFLAGS) -c ./src/util.c 

clean:
	rm *.o *.out
