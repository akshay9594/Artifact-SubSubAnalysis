CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

icm7_cetus:
	$(CC) -o icm7_cetus $(CCFLAGS) IncompleteCholesky_m7.c -lm 

clean:
	rm -f icm7_cetus $(OUTPUT)