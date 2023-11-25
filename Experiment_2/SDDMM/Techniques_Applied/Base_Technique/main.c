/*
Copyright (C) 1991-2022 Free Software Foundation, Inc.
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

 The SDDMM kernel

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include "util.h"
void sddmm_CPU_CSR(int * col_ptr, int * col_ind, double * nnz_val, double * W, double * H, double * p, int n_rows, int k, int nonzeros);
/* void Par_sddmm_CPU_CSR(int par_row_ptr, int* col_ind, double* nnz_val, double* W, */
/*                      double H, double* p,int n_rows ,int k, int nonzeros); */
int * row_val;
int * col_val;
double * nnz_val;
int * col;
/* int rowP; */
int convertStrtoArr(char * str)
{
	/* get length of string str */
	int str_length = strlen(str);
	/* create an array with size as string */
	/* length and initialize with 0 */
	int arr = 0;
	int j = 0, i, sum = 0;
	/* Traverse the string */
	#pragma loop name convertStrtoArr#0 
	for (i=0; str[i]!='\0'; i ++ )
	{
		/* if str[i] is ', ' then split */
		if (str[i]==',')
		{
			continue;
		}
		if (str[i]==' ')
		{
			/* Increment j to point to next */
			/* array location */
			break;
		}
		else
		{
			/* subtract str[i] by 48 to convert it to int */
			/* Generate number by multiplying 10 and adding */
			/* (int)(str[i]) */
			if ((str[i]-48)>=0)
			{
				arr=((arr*10)+(str[i]-48));
			}
		}
	}
	return arr;
}

void sddmm_CPU_CSR(int * col_ptr, int * row_ind, double * nnz_val, double * W, double * H, double * p, int n_cols, int k, int nonzeros)
{
	/* reduction(+:rmse) */
	int i, r, ind, t, holder;
	double sm;
	holder=1;
	col_ptr[0]=0;
	r=col_val[0];
	#pragma cetus private(i) 
	#pragma loop name sddmm_CPU_CSR#0 
	for (i=0; i<nonzeros; i ++ )
	{
		if (col_val[i]!=r)
		{
			col_ptr[holder ++ ]=i;
			/* rowP[holder] = i; */
			r=col_val[i];
		}
	}
	col_ptr[holder]=nonzeros;
	/* #pragma omp parallel for private(sm,r,ind,t) */
	#pragma cetus private(ind, r, sm, t) 
	#pragma loop name sddmm_CPU_CSR#1 
	for (r=0; r<n_cols;  ++ r)
	{
		#pragma cetus private(ind, sm, t) 
		#pragma loop name sddmm_CPU_CSR#1#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((((1L+(5L*col_ptr[(1L+r)]))+(-5L*col_ptr[r]))+((3L*k)*col_ptr[(1L+r)]))+((-3L*k)*col_ptr[r])))) private(ind, sm, t)
		for (ind=col_ptr[r]; ind<col_ptr[r+1];  ++ ind)
		{
			sm=0;
			#pragma cetus private(t) 
			#pragma loop name sddmm_CPU_CSR#1#0#0 
			/* #pragma cetus reduction(+: sm)  */
			for (t=0; t<k;  ++ t)
			{
				sm+=(W[(r*k)+t]*H[(row_ind[ind]*k)+t]);
			}
			p[ind]=(sm*nnz_val[ind]);
			/* Scaling of non-zero elements of the sparse matrix */
		}
	}
}

/* void Par_sddmm_CPU_CSR(int par_row_ptr, int* col_ind, double* nnz_val, double* W, */
/*                      double H, double* p,int n_rows ,int k, int nonzeros){ */
	/*        reduction(+:rmse) */
	/*     int i,r, ind,t,holder; */
	/*     double sm; */
	/*      holder =0; */
	/*     par_row_ptr[0]=0; */
	/*      r = row_val[0]; */
	/*       for(i =0; i < nonzeros; i++){ */
		/*         if(row_val[i] != r){ */
			/*             par_row_ptr[holder++] = i; */
			/*            rowP[holder] = i; */
			/*             r = row_val[i]; */
		/*         } */
	/*     } */
	/*     par_row_ptr[holder] = nonzeros; */
	/*     #pragma omp parallel for private(sm,r,ind,t) */
	/*     for (r = 0; r < n_rows; ++r){ */
		/*         for (ind = par_row_ptr[r]; ind < par_row_ptr[r+1]; ++ind){ */
			/*             sm=0; */
			/*             for (t = 0; t < k; ++t){ */
				/*                 sm += W[r k + t] * H[col_ind[ind] * k + t]; */
			/*             } */
			/*             p[ind] = sm nnz_val[ind];     Scaling of non-zero elements of the sparse matrix */
		/*         }                 */
	/*     }  */
/* } */
int main(int argc, char * argv[])
{
	char * file_path = argv[1];
	int i, j, k, s_factor, holder, count;
	int num_rows, num_cols, nonzeros, num_runs, failed;
	size_t len, read;
	struct timeval start, end;
	FILE * fp;
	char * line = (void * )0;
	char * ptr;
	char * rowstr = (void * )0;
	char * temp_str;
	double seconds, total_time, init_time;
	double * W;
	double * H;
	double * P;
	double * ParallelP;
	char delim[] = " ";
	if (argc<2)
	{
		printf("Input missing arguments, you need to specify input list file\n");
	}
	/* std::string inputMatrix; */
	/* ,startTT, endTT; */
	s_factor=100;
	count=0;
	num_runs=1;
	fp=fopen(file_path, "r");
	if (fp==((void * )0))
	{
		exit(1);
	}
	while ((read=getline( & line,  & len, fp))!=( - 1))
	{
		if (line[0]=='%')
		{
			continue;
		}
		else
		{
			break;
		}
	}
	ptr=strtok(line, delim);
	rowstr=((char * )malloc(sizeof line));
	while (ptr!=((void * )0))
	{
		if (count==0)
		{
			memcpy(rowstr, ptr, sizeof ptr);
			num_rows=convertStrtoArr(ptr);
		}
		else
		{
			if (count==1)
			{
				num_cols=convertStrtoArr(ptr);
			}
			else
			{
				if (count==2)
				{
					nonzeros=convertStrtoArr(ptr);
				}
			}
		}
		count ++ ;
		ptr=strtok((void * )0, delim);
	}
	row_val=((int * )malloc(sizeof (int)*(nonzeros+1)));
	col_val=((int * )malloc(sizeof (int)*(nonzeros+1)));
	nnz_val=((double * )malloc(sizeof (double)*(nonzeros+1)));
	i=0;
	while ((read=getline( & line,  & len, fp))!=( - 1))
	{
		char * ptr = strtok(line, delim);
		count=0;
		while (ptr!=((void * )0))
		{
			if (count==0)
			{
				row_val[i]=convertStrtoArr(ptr);
			}
			else
			{
				if (count==1)
				{
					col_val[i]=convertStrtoArr(ptr);
				}
				else
				{
					if (count==2)
					{
						nnz_val[i]=strtod(ptr,  & temp_str);
					}
				}
			}
			count ++ ;
			ptr=strtok((void * )0, delim);
		}
		i ++ ;
	}
	col=malloc(sizeof (int)*(nonzeros+1));
	/* rowP = malloc(sizeof(int)(nonzeros+1)); */
	total_time=0.0;
	W=((double * )malloc(sizeof (double)*((num_rows*s_factor)+s_factor)));
	H=((double * )malloc(sizeof (double)*((num_cols*s_factor)+s_factor)));
	P=((double * )malloc(sizeof (double)*(nonzeros+1)));
	/* ParallelP =  (double)malloc(sizeof(double)*(nonzeros+1)); */
	initialize(W, num_rows, s_factor);
	initialize(H, num_cols, s_factor);
	#pragma cetus private(k, seconds) 
	#pragma loop name main#0 
	/* #pragma cetus reduction(+: total_time)  */
	for (k=0; k<num_runs; k ++ )
	{
		gettimeofday( & start, (void * )0);
		sddmm_CPU_CSR(col, row_val, nnz_val, W, H, P, num_cols, s_factor, nonzeros);
		gettimeofday( & end, (void * )0);
		seconds=((end.tv_sec+(((double)end.tv_usec)/1000000))-(start.tv_sec+(((double)start.tv_usec)/1000000)));
		total_time+=seconds;
		/* gettimeofday(&start,NULL); */
		/*    Parallel Run */
		/*   Par_sddmm_CPU_CSR(rowP,col_val,nnz_val,W,H,ParallelP,num_rows,s_factor,nonzeros); */
		/*   gettimeofday(&end, NULL); */
		/*   seconds = (end.tv_sec + (double)end.tv_usec1000000) - (start.tv_sec + (double)start.tv_usec/1000000);  */
		/*   total_timeP+= seconds; */
	}
	/* failed = 0; */
	/* for (i = 0; i < num_rows; ++i) { */
		/*     printf("P[%d]=%f,\n",i ,P[i]); */
		/*   if(P[i]-ParallelP[i] > 10e-4) failed=1; */
	/* } */
	/* if(failed == 1){ */
		/*   printf("Verification failed!!!"); */
		/*   exit(0); */
	/* } */
	printf("Input File Read successfully\n");
	printf("-->Avg time taken by the serial kernel= %f s\n",total_time/num_runs);
	/* printf("-->Avg time taken by the Parallel kernel for %d runs = %f s\n", num_runs,total_timePnum_runs); */
	fclose(fp);
	if (line)
	{
		free(line);
	}
	exit(0);
	return 0;
}
