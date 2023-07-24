CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

grshm:
	$(CC) -o grshm $(CCFLAGS) ../../utilities/polybench.c gramschmidt.c -DPOLYBENCH_TIME -lm 

clean:
	rm -f grshm $(OUTPUT)