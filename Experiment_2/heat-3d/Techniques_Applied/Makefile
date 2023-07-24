CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

heat-3d:
	$(CC) -o heat-3d $(CCFLAGS) ../../utilities/polybench.c heat-3d.c -DPOLYBENCH_TIME -lm 

clean:
	rm -f heat-3d $(OUTPUT)