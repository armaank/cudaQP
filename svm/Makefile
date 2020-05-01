CC = gcc
CFLAGS  = -Wall

svm.out: svm.c 
	$(CC) $(CFLAGS) svm.c ldl_interface.o kkt.o csc.o linalg.o ldl.o linsys.o scaling.o timer.o ctrl.o proj.o error.o auxil.o polish.o qp.o -o svm.out -lm -lamd

qp.o: ../src/qp.c ../include/qp.h
	$(CC) $(CFLAGS) -c ../src/qp.c ../src/polish.c ../src/proj.c ../src/linsys.c ../src/auxil.c ../src/scaling.c

polish.o: ../src/polish.c ../include/polish.h
	$(CC) $(CFLAGS) -c ../src/polish.c ../src/proj.c ../src/linsys.c ../src/linalg.c ../src/kkt.c ../src/auxil.c 

auxil.o: ../src/auxil.c ../include/auxil.h
	$(CC) $(CFLAGS) -c ../src/auxil.c ../src/scaling.c ../src/proj.c

# might be able to remove linalg.c from scaling
scaling.o: ../src/scaling.c ../include/scaling.h
	$(CC) $(CFLAGS) -c ../src/scaling.c ../src/linalg.c 

timer.o: ../src/timer.c  ../include/timer.h
	$(CC) $(CFLAGS) -c ../src/timer.c 

error.o: ../src/error.c ../include/error.h
	$(CC) $(CFLAGS) -c ../src/error.c 

ctrl.o: ../src/ctrl.c ../include/ctrl.h
	$(CC) $(CFLAGS) -c ../src/ctrl.c 

proj.o: ../src/proj.c ../include/proj.h
	$(CC) $(CFLAGS) -c ../src/proj.c 

linsys.o: ../src/linsys.c ../include/linsys.h
	$(CC) $(CFLAGS) -c ../src/linsys.c ../ldl/ldl_interface.c

ldl_interface.o: ../ldl/ldl_interface.c ../ldl/ldl_interface.h
	$(CC) $(CFLAGS) -c ../ldl/ldl_interface.c ../ldl/ldl_src/ldl.c ../src/kkt.c -lamd

kkt.o: ../src/kkt.c ../include/kkt.h
	$(CC) $(CFLAGS) -c ../src/kkt.c ../src/csc.c

csc.o: ../src/csc.c ../include/csc.h
	$(CC) $(CFLAGS) -c ../src/csc.c ../src/linalg.c

linalg.o: ../src/linalg.c ../include/linalg.h
	$(CC) $(CFLAGS) -c ../src/linalg.c -lm

ldl.o: ../ldl/ldl_src/ldl.c ../ldl/ldl_src/ldl.h
	$(CC) $(CFLAGS) -c ../ldl/ldl_src/ldl.c

all: ldl.o linalg.o csc.o kkt.o ldl_interface.o linsys.o proj.o ctrl.o error.o timer.o scaling.o auxil.o polish.o qp.o

clean:
	rm *.o *.out

clear: 
	rm *.o
