CC = gcc
CFLAGS  = -Wall

# make: test.c ldl.o
# 	$(CC) $(CFLAGS) test.c ldl_unittest.h ldl.o -o tests.out


# osqp


# ldl lib - done?
ldl.o: ./ldl/ldl.c ./ldl/ldl.h
	$(CC) $(CFLAGS) -c ./ldl/ldl.c 

# util - done?
util.o: ./src/util.c ./include/util.h ./include/glob_opts.h
	$(CC) $(CFLAGS) -c ./src/util.c 

clean:
	rm *.o *.out
