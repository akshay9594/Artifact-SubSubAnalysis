CC = gcc
CCFLAGS = -g -fopenmp -O3
OUTPUT = *.out

gramschmidt:
	$(CC) -o gramschmidt $(CCFLAGS) ../../../utilities/polybench.c gramschmidt.c -DPOLYBENCH_TIME -lm 

clean:
	rm -f gramschmidt $(OUTPUT)