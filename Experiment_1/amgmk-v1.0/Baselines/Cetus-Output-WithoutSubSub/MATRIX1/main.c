/*
Copyright (C) 1991-2021 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it andor
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https:www.gnu.org/licenses/>. 
*/
/*
This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it. 
*/
/*
glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default. 
*/
/*
wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISOIEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters
*/
/*
BHEADER
 * (c) 2007   The Regents of the University of California               *
 *                                                                      *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright       *
 * notice and disclaimer.                                               *
 *                                                                      *
 *EHEADER
*/
/* -------------- */
/*  A micro kernel  */
/* -------------- */
#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "headers.h"
/*  */
const int testIter = 500;
double totalWallTime = 0.0;
double total_target_loop_time = 0.0;
/*  */
void test_Matvec();
void test_Relax();
void test_Axpy();
/*  */
int main(int argc, char * argv[])
{
	double t0 = 0.0, t1 = 0.0, del_wtime = 0.0;
	int max_num_threads;
	printf("\n");
	printf("//------------ \n");
	printf("// \n");
	printf("//  CORAL  AMGmk Benchmark Version 1.0 \n");
	printf("// \n");
	printf("//------------ \n");
	#pragma omp parallel 
	 #pragma omp master 
	max_num_threads=omp_get_num_threads();
	printf("\nmax_num_threads = %d \n\n", max_num_threads);
	printf("\n testIter   = %d \n\n", testIter);
	t0=omp_get_wtime();
	/* Matvec */
	totalWallTime=0.0;
	total_target_loop_time=0.0;
	test_Matvec();
	printf("\n");
	printf("//------------ \n");
	printf("// \n");
	printf("//   MATVEC\n");
	printf("// \n");
	printf("//------------ \n");
	printf("\nWall time = %f seconds. \n", totalWallTime);
	printf("Target loop time=%f seconds\n", total_target_loop_time);
	/* Relax */
	/* totalWallTime = 0.0; */
	/* test_Relax(); */
	/* printf("\n"); */
	/* printf("------------ \n"); */
	/* printf("\n"); */
	/* printf("  Relax\n"); */
	/* printf("\n"); */
	/* printf("------------ \n"); */
	/* printf("\nWall time = %f seconds. \n", totalWallTime); */
	/* Axpy */
	/* totalWallTime = 0.0; */
	/* test_Axpy(); */
	/* printf("\n"); */
	/* printf("------------ \n"); */
	/* printf("\n"); */
	/* printf("  Axpy\n"); */
	/* printf("\n"); */
	/* printf("------------ \n"); */
	/* printf("\nWall time = %f seconds. \n", totalWallTime); */
	/* t1 = omp_get_wtime();; */
	/*  del_wtime = t1 - t0; */
	/* printf("\nTotal Wall time = %f seconds. \n", del_wtime); */
	return 0;
}

void test_Matvec()
{
	double t0 = 0.0, t1 = 0.0;
	hypre_CSRMatrix * A;
	hypre_Vector * x, * y, * sol;
	int nx, ny, nz, i;
	double * values;
	double * y_data, * sol_data;
	double error, diff;
	double time_array[500];
	nx=60;
	/* size per proc nxny*nz */
	ny=60;
	nz=60;
	values=((double * )hypre_CAlloc((unsigned int)4, (unsigned int)sizeof (double)));
	values[0]=6;
	values[1]=( - 1);
	values[2]=( - 1);
	values[3]=( - 1);
	A=GenerateSeqLaplacian(nx, ny, nz, values,  & y,  & x,  & sol);
	hypre_CSRMatrixSetRownnz(A);
	hypre_SeqVectorSetConstantValues(x, 1);
	hypre_SeqVectorSetConstantValues(y, 0);
	t0=omp_get_wtime();
	#pragma cetus private(i) 
	#pragma loop name test_Matvec#0 
	for (i=0; i<testIter;  ++ i)
	{
		hypre_CSRMatrixMatvec(1, A, x, 0, y, time_array, i);
	}
	t1=omp_get_wtime();
	totalWallTime+=(t1-t0);
	#pragma cetus private(i) 
	#pragma loop name test_Matvec#1 
	#pragma cetus reduction(+: total_target_loop_time) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*testIter)))) private(i) reduction(+: total_target_loop_time)
	for (i=0; i<testIter; i ++ )
	{
		total_target_loop_time+=time_array[i];
	}
	y_data=y->data;
	sol_data=sol->data;
	error=0;
	#pragma cetus private(diff, i) 
	#pragma loop name test_Matvec#2 
	for (i=0; i<(((nx*ny)*nz)*5); i ++ )
	{
		diff=fabs(y_data[i]-sol_data[i]);
		if (diff>error)
		{
			error=diff;
		}
	}
	if (error>0)
	{
		printf(" \n Matvec: error: %e\n", error);
	}
	(hypre_Free((char * )values), (values=((void * )0)));
	hypre_CSRMatrixDestroy(A);
	hypre_SeqVectorDestroy(x);
	hypre_SeqVectorDestroy(y);
	hypre_SeqVectorDestroy(sol);
}

void test_Relax()
{
	double t0 = 0.0, t1 = 0.0;
	hypre_CSRMatrix * A;
	hypre_Vector * x, * y, * sol;
	int nx, ny, nz, i;
	double * values;
	double * x_data;
	double diff, error;
	nx=50;
	/* size per proc nxny*nz */
	ny=50;
	nz=50;
	values=((double * )hypre_CAlloc((unsigned int)4, (unsigned int)sizeof (double)));
	values[0]=6;
	values[1]=( - 1);
	values[2]=( - 1);
	values[3]=( - 1);
	A=GenerateSeqLaplacian(nx, ny, nz, values,  & y,  & x,  & sol);
	hypre_SeqVectorSetConstantValues(x, 1);
	t0=omp_get_wtime();
	#pragma cetus private(i) 
	#pragma loop name test_Relax#0 
	for (i=0; i<testIter;  ++ i)
	{
		hypre_BoomerAMGSeqRelax(A, sol, x);
	}
	t1=omp_get_wtime();
	totalWallTime+=(t1-t0);
	x_data=x->data;
	error=0;
	#pragma cetus private(diff, i) 
	#pragma loop name test_Relax#1 
	for (i=0; i<((nx*ny)*nz); i ++ )
	{
		diff=fabs(x_data[i]-1);
		if (diff>error)
		{
			error=diff;
		}
	}
	if (error>0)
	{
		printf(" \n Relax: error: %e\n", error);
	}
	(hypre_Free((char * )values), (values=((void * )0)));
	hypre_CSRMatrixDestroy(A);
	hypre_SeqVectorDestroy(x);
	hypre_SeqVectorDestroy(y);
	hypre_SeqVectorDestroy(sol);
}

void test_Axpy()
{
	double t0 = 0.0, t1 = 0.0;
	hypre_Vector * x, * y;
	int nx, i;
	double alpha = 0.5;
	double diff, error;
	double * y_data;
	nx=125000;
	/* size per proc  */
	x=hypre_SeqVectorCreate(nx);
	y=hypre_SeqVectorCreate(nx);
	hypre_SeqVectorInitialize(x);
	hypre_SeqVectorInitialize(y);
	hypre_SeqVectorSetConstantValues(x, 1);
	hypre_SeqVectorSetConstantValues(y, 1);
	t0=omp_get_wtime();
	#pragma cetus private(i) 
	#pragma loop name test_Axpy#0 
	for (i=0; i<testIter;  ++ i)
	{
		hypre_SeqVectorAxpy(alpha, x, y);
	}
	t1=omp_get_wtime();
	y_data=y->data;
	error=0;
	#pragma cetus private(diff, i) 
	#pragma loop name test_Axpy#1 
	for (i=0; i<nx; i ++ )
	{
		diff=fabs((y_data[i]-1)-(0.5*((double)testIter)));
		if (diff>error)
		{
			error=diff;
		}
	}
	if (error>0)
	{
		printf(" \n Axpy: error: %e\n", error);
	}
	totalWallTime+=(t1-t0);
	hypre_SeqVectorDestroy(x);
	hypre_SeqVectorDestroy(y);
}
