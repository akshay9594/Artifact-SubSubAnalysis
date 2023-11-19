CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

fdtd-2d:
	$(CC) -o fdtd-2d $(CCFLAGS) ../../../utilities/polybench.c fdtd-2d.c -DPOLYBENCH_TIME -lm 

clean:
	rm -f fdtd-2d $(OUTPUT)