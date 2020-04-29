CC = gcc
CFLAGS  = -Wall

test.out: test_basic_qp.c test_basic_qp.h
	$(CC) $(CFLAGS) test_basic_qp.c ldl_interface.o kkt.o csc.o linalg.o ldl.o linsys.o scaling.o timer.o ctrl.o proj.o error.o auxil.o polish.o qp.o -o tests.out -lm -lamd

qp.o: ../../qp.c ../../../include/qp.h
	$(CC) $(CFLAGS) -c ../../qp.c ../../polish.c ../../proj.c ../../linsys.c ../../auxil.c ../../scaling.c

polish.o: ../../polish.c ../../../include/polish.h
	$(CC) $(CFLAGS) -c ../../polish.c ../../proj.c ../../linsys.c ../../linalg.c ../../kkt.c ../../auxil.c 

auxil.o: ../../auxil.c ../../../include/auxil.h
	$(CC) $(CFLAGS) -c ../../auxil.c ../../scaling.c ../../proj.c

# might be able to remove linalg.c from scaling
scaling.o: ../../scaling.c ../../../include/scaling.h
	$(CC) $(CFLAGS) -c ../../scaling.c ../../linalg.c 

timer.o: ../../timer.c  ../../../include/timer.h
	$(CC) $(CFLAGS) -c ../../timer.c 

error.o: ../../error.c ../../../include/error.h
	$(CC) $(CFLAGS) -c ../../error.c 

ctrl.o: ../../ctrl.c ../../../include/ctrl.h
	$(CC) $(CFLAGS) -c ../../ctrl.c 

proj.o: ../../proj.c ../../../include/proj.h
	$(CC) $(CFLAGS) -c ../../proj.c 

linsys.o: ../../linsys.c ../../../include/linsys.h
	$(CC) $(CFLAGS) -c ../../linsys.c ../../../ldl/ldl_interface.c

ldl_interface.o: ../../../ldl/ldl_interface.c ../../../ldl/ldl_interface.h
	$(CC) $(CFLAGS) -c ../../../ldl/ldl_interface.c ../../../ldl/ldl_src/ldl.c ../../kkt.c -lamd

kkt.o: ../../kkt.c ../../../include/kkt.h
	$(CC) $(CFLAGS) -c ../../kkt.c ../../csc.c

csc.o: ../../csc.c ../../../include/csc.h
	$(CC) $(CFLAGS) -c ../../csc.c ../../linalg.c

linalg.o: ../../linalg.c ../../../include/linalg.h
	$(CC) $(CFLAGS) -c ../../linalg.c -lm

ldl.o: ../../../ldl/ldl_src/ldl.c ../../../ldl/ldl_src/ldl.h
	$(CC) $(CFLAGS) -c ../../../ldl/ldl_src/ldl.c

all: ldl.o timer.o linalg.o csc.o kkt.o ldl_interface.o linsys.o proj.o ctrl.o error.o scaling.o auxil.o polish.o qp.o

clean:
	rm *.o *.out
