CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

syrk:
	$(CC) -o syrk $(CCFLAGS) ../../../utilities/polybench.c syrk.c -DPOLYBENCH_TIME -lm 

clean:
	rm -f syrk $(OUTPUT)