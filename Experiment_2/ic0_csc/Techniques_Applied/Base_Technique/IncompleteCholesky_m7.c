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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 #include <omp.h> 
#include <sys/time.h>
#include <math.h>
#include <assert.h>
//#include "ic0_csc_m7.h"
int * row_val;
int * col_val;
double * nnz_val;
int * row;
int * col;
double * val;
double * val_par;
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
			arr=((arr*10)+(str[i]-48));
		}
	}
	return arr;
}

int main(int argc, char * argv[])
{
	char * file_path = argv[1];
	int i, j, k, m, l;
	struct timeval start, end;
	FILE * fp;
	char * line = (void * )0;
	int num_rows, num_cols, nonzeros, colCnt = 0, nnzCnt = 0, count = 0;
	size_t len = 0;
	size_t read;
	char matrix_details[100];
	
	int x, y, col_num = 0;
	int num_runs = 1;
	double seconds, total_time = 0.0;
	if (argc<2)
	{
		printf("Input missing arguments, you need to specify input list file\n");
	}
	/* std::string inputMatrix; */
	/* ,startTT, endTT; */
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

	char delim[] = " ";
	char * ptr = strtok(line, delim);
	char * rowstr = (void * )0;
	char * temp_str;
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
	row_val=((int * )malloc(sizeof (int)*nonzeros));
	col_val=((int * )malloc(sizeof (int)*nonzeros));
	nnz_val=((double * )malloc(sizeof (double)*nonzeros));
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
	row=malloc(sizeof (int)*nonzeros);
	col=malloc(sizeof (int)*nonzeros);
	val=malloc(sizeof (double)*nonzeros);
	col[0]=0;
	#pragma cetus private(j, y) 
	#pragma loop name main#0 
	for (j=0; j<nonzeros; j ++ )
	{
		y=col_val[j];
		y=(y-1);
		/* x = row_val[j]; */
		/* x=x-1; */
		if (y==col_num)
		{
			val[nnzCnt]=nnz_val[nnzCnt];
			/* row[nnzCnt]=x; */
			colCnt ++ ;
			nnzCnt ++ ;
		}
		else
		{
			/* New col */
			col[col_num+1]=(col[col_num]+colCnt);
			col_num ++ ;
			colCnt=1;
			val[nnzCnt]=nnz_val[nnzCnt];
			/* row[nnzCnt]=x; */
			nnzCnt ++ ;
		}
	}
	#pragma cetus private(i, j, x) 
	#pragma loop name main#1 
	for (i=0; i<nonzeros; i ++ )
	{
		#pragma cetus private(j, x) 
		#pragma loop name main#1#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(5L*col[(1L+i)]))+(-5L*col[i])))) private(j, x)
		for (j=col[i]; j<col[i+1]; j ++ )
		{
			x=row_val[j];
			x=(x-1);
			row[j]=x;
			/* printf("row[%d]=%d\n", j, row[j]); */
		}
	}
	col[num_cols]=(col[num_cols-1]+colCnt);
	/* last col */
	#pragma cetus private(i, k, l, m) 
	#pragma loop name main#2 
	for (k=0; k<num_runs; k ++ )
	{
		 gettimeofday(&start,NULL); 
		/* ic0_csr(num_cols,  val, col, row); */
		#pragma cetus private(i, l, m) 
		#pragma cetus lastprivate(k) 
		#pragma loop name main#2#0 
		for (i=0; i<(num_cols-1); i ++ )
		{
			val[col[i]]=sqrt(val[col[i]]);
			/* S1 */
			#pragma cetus private(m) 
			#pragma loop name main#2#0#0 
			for (m=(col[i]+1); m<col[i+1]; m ++ )
			{
				val[m]=(val[m]/val[col[i]]);
				/* S2 */
			}
			#pragma cetus private(l, m) 
			#pragma cetus lastprivate(k) 
			#pragma loop name main#2#0#1 
			for (m=(col[i]+1); m<col[i+1]; m ++ )
			{
				#pragma cetus private(l) 
				#pragma cetus lastprivate(k) 
				#pragma loop name main#2#0#1#0 
				for (k=col[row[m]]; k<col[row[m]+1]; k ++ )
				{
					#pragma cetus private(l) 
					#pragma loop name main#2#0#1#0#0 
					for (l=m; l<col[i+1]; l ++ )
					{
						if ((row[l]==row[k])&&(row[l+1]<=row[k]))
						{
							val[k]-=(val[m]*val[l]);
							/* S3 */
						}
					}
				}
			}
		}
		 gettimeofday(&end, NULL); 
		 double seconds = (end.tv_sec + (double)end.tv_usec/1000000) - (start.tv_sec + (double)start.tv_usec/1000000);  
		 total_time+= seconds; 
	}
	printf("Time taken by the Cetus Parallel kernel= %f s\n",total_time/num_runs);
	fclose(fp);
	if (line)
	{
		free(line);
	}
	exit(0);
	return 0;
}
