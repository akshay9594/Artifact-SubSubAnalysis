#BHEADER****************************************************************
# (c) 2007   The Regents of the University of California               *
#                                                                      *
# See the file COPYRIGHT_and_DISCLAIMER for a complete copyright       *
# notice and disclaimer.                                               *
#                                                                      *
#EHEADER****************************************************************
CC       = gcc

CFLAGS   = -O3 -fopenmp 

LDFLAGS  = -lm -mcmodel=small -fopenmp 

PROG     = ua.A.x

OBJS     = ua.o \
           convect.o adapt.o diffuse.o \
           move.o mason.o \
           precond.o utils.o \
           transfer.o verify.o \
           setup.o print_results.o \
           c_timers.o wtime.o

all : $(PROG)

$(PROG) : $(OBJS)
	${CC} -O3 -o $(PROG) $(OBJS) $(LDFLAGS)


clean :
	rm -f *.o $(PROG) core job.out *~ 


.SUFFIXES : .o  .c

#*****************************************************************
#* Rules for C object file dependencies
#*****************************************************************
.c.o :
	$(CC) -o $@ -c $(CFLAGS) $*.c 






