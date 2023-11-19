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
 * Copyright (c) 2006   The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by the HYPRE team. UCRL-CODE-222953.
 * All rights reserved.
 *
 * This file is part of HYPRE (see http:www.llnl.gov/CASC/hypre/).
 * Please see the COPYRIGHT_and_LICENSE file for the copyright notice, 
 * disclaimer, contact information and the GNU Lesser General Public License.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 *
 * HYPRE is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * $Revision: 2.10 $
 ***********************************************************************EHEADER
*/
/*

 *
 * Matvec functions for hypre_CSRMatrix class.
 *

*/
#include "headers.h"
#include <assert.h>
#include <stdlib.h>
 #include "omp.h"
/*
--------------------------------------------------------------------------
 hypre_CSRMatrixMatvec
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixMatvec(double alpha, hypre_CSRMatrix * A, hypre_Vector * x, double beta, hypre_Vector * y, double time_array[], int id)
{
	double * A_data = A->data;
	int * A_i = A->i;
	int * A_j = A->j;
	int num_rows = A->num_rows;
	int num_cols = A->num_cols;
	int * A_rownnz = A->rownnz;
	int num_rownnz = A->num_rownnz;
	double * x_data = x->data;
	double * y_data = y->data;
	int x_size = x->size;
	int y_size = y->size;
	int num_vectors = x->num_vectors;
	int idxstride_y = y->idxstride;
	int vecstride_y = y->vecstride;
	int idxstride_x = x->idxstride;
	int vecstride_x = x->vecstride;
	double temp, tempx;
	int i, j, jj;
	int m;
	double xpar = 0.7;
	int ierr = 0;
	double t0, t1;
	/*
	---------------------------------------------------------------------
	  Check for size compatibility.  Matvec returns ierr = 1 if
	    *  length of X doesn't equal the number of columns of A,
	    *  ierr = 2 if the length of Y doesn't equal the number of rows
	    *  of A, and ierr = 3 if both are true.
	    *
	    *  Because temporary vectors are often used in Matvec, none of 
	    *  these conditions terminates processing, and the ierr flag
	    *  is informational only.
	    *--------------------------------------------------------------------
	*/
	int adiag;
	int irownnz;
	int jj_0;
	int jj_1;
	int jj_2;
	if ( ! (num_vectors==y->num_vectors))
	{
		fprintf(stderr, "hypre_assert failed: %s\n", "num_vectors == hypre_VectorNumVectors(y)");
		hypre_error_handler("csr_matvec.c", 95, 1);
	}
	;
	if (num_cols!=x_size)
	{
		ierr=1;
	}
	if (num_rows!=y_size)
	{
		ierr=2;
	}
	if ((num_cols!=x_size)&&(num_rows!=y_size))
	{
		ierr=3;
	}
	/*
	-----------------------------------------------------------------------
	 Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
	    *-----------------------------------------------------------------------
	*/
	if (alpha==0.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvec#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+((3L*num_rows)*num_vectors)))) private(i)
		for (i=0; i<(num_rows*num_vectors); i ++ )
		{
			y_data[i]*=beta;
		}
		return ierr;
	}
	/*
	-----------------------------------------------------------------------
	 y = (beta/alpha)*y
	    *-----------------------------------------------------------------------
	*/
	temp=(beta/alpha);
	if (temp!=1.0)
	{
		if (temp==0.0)
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvec#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+((3L*num_rows)*num_vectors)))) private(i)
			for (i=0; i<(num_rows*num_vectors); i ++ )
			{
				y_data[i]=0.0;
			}
		}
		else
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvec#2 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+((3L*num_rows)*num_vectors)))) private(i)
			for (i=0; i<(num_rows*num_vectors); i ++ )
			{
				y_data[i]*=temp;
			}
		}
	}
	irownnz=0;
	#pragma cetus private(adiag, i) 
	#pragma loop name hypre_CSRMatrixMatvec#3 
	for (i=0; i<num_rows; i ++ )
	{
		adiag=(A_i[i+1]-A_i[i]);
		if (adiag>0)
		{
			A_rownnz[irownnz ++ ]=i;
		}
	}
	/* hypre_CSRMatrixRownnz(A) = A_rownnz; */
	/*
	-----------------------------------------------------------------
	 y += A*x
	    *-----------------------------------------------------------------
	*/
	/* use rownnz pointer to do the Ax multiplication  when num_rownnz is smaller than num_rows */
	t0=omp_get_wtime();
	if (num_rownnz<(xpar*num_rows))
	{
		#pragma cetus private(i, jj, jj_0, tempx) 
		#pragma loop name hypre_CSRMatrixMatvec#4 
		#pragma cetus parallel 
		/* #pragma omp parallel for private(i, jj, tempx) schedule(static) */
		#pragma omp parallel for if(((-1+num_rownnz)<=irownnz)) private(i, jj, jj_0, tempx)
		for (i=0; i<num_rownnz; i ++ )
		{
			/*
			
			 for (jj = A_i[m]; jj < A_i[m+1]; jj++)
			          * {
				          *         j = A_j[jj];   
				          *  y_data[m] += A_data[jj] * x_data[j];
			          * }
			*/
			if (num_vectors==1)
			{
				tempx=y_data[A_rownnz[i]];
				/* Normalized Loop */
				#pragma cetus lastprivate(jj_0) 
				#pragma loop name hypre_CSRMatrixMatvec#4#0 
				/* #pragma cetus reduction(+: tempx)  */
				for (jj_0=0; jj_0<=((-1+A_i[1+A_rownnz[i]])+(-1*A_i[A_rownnz[i]])); jj_0 ++ )
				{
					tempx+=(A_data[jj_0+A_i[A_rownnz[i]]]*x_data[A_j[jj_0+A_i[A_rownnz[i]]]]);
				}
				jj=(jj_0+A_i[A_rownnz[i]]);
				y_data[A_rownnz[i]]=tempx;
			}
			/* else */
			/*    for ( j=0; j<num_vectors; ++j ) */
			/*    { */
				/*       tempx = y_data[ jvecstride_y + A_rownnz[i]*idxstride_y ]; */
				/*       for (jj = A_i[A_rownnz[i]]; jj < A_i[A_rownnz[i]+1]; jj++)  */
				/*          tempx +=  A_data[jj] x_data[ j*vecstride_x + A_j[jj]*idxstride_x ]; */
				/*       y_data[ jvecstride_y + A_rownnz[i]*idxstride_y] = tempx; */
			/*    } */
		}
	}
	else
	{
		/* #pragma omp parallel for private(i,jj,temp) schedule(static) */
		#pragma cetus private(i, j, jj, jj_1, jj_2, temp) 
		#pragma loop name hypre_CSRMatrixMatvec#5 
		for (i=0; i<num_rows; i ++ )
		{
			if (num_vectors==1)
			{
				temp=y_data[i];
				/* Normalized Loop */
				#pragma cetus lastprivate(jj_1) 
				#pragma loop name hypre_CSRMatrixMatvec#5#0 
				#pragma cetus reduction(+: temp) 
				#pragma cetus parallel 
				#pragma omp parallel for if((10000<((1L+(3L*A_i[(1L+i)]))+(-3L*A_i[i])))) lastprivate(jj_1) reduction(+: temp)
				for (jj_1=0; jj_1<=((-1+A_i[1+i])+(-1*A_i[i])); jj_1 ++ )
				{
					temp+=(A_data[jj_1+A_i[i]]*x_data[A_j[jj_1+A_i[i]]]);
				}
				jj=(jj_1+A_i[i]);
				y_data[i]=temp;
			}
			else
			{
				#pragma cetus private(j, jj, jj_2, temp) 
				#pragma loop name hypre_CSRMatrixMatvec#5#1 
				for (j=0; j<num_vectors;  ++ j)
				{
					temp=y_data[(j*vecstride_y)+(i*idxstride_y)];
					/* Normalized Loop */
					#pragma cetus lastprivate(jj_2) 
					#pragma loop name hypre_CSRMatrixMatvec#5#1#0 
					#pragma cetus reduction(+: temp) 
					#pragma cetus parallel 
					#pragma omp parallel for if((10000<((1L+(3L*A_i[(1L+i)]))+(-3L*A_i[i])))) lastprivate(jj_2) reduction(+: temp)
					for (jj_2=0; jj_2<=((-1+A_i[1+i])+(-1*A_i[i])); jj_2 ++ )
					{
						temp+=(A_data[jj_2+A_i[i]]*x_data[(j*vecstride_x)+(A_j[jj_2+A_i[i]]*idxstride_x)]);
					}
					jj=(jj_2+A_i[i]);
					y_data[(j*vecstride_y)+(i*idxstride_y)]=temp;
				}
			}
		}
	}
	t1=omp_get_wtime();
	time_array[id]=(t1-t0);
	/*
	-----------------------------------------------------------------
	 y = alpha*y
	    *-----------------------------------------------------------------
	*/
	if (alpha!=1.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvec#6 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+((3L*num_rows)*num_vectors)))) private(i)
		for (i=0; i<(num_rows*num_vectors); i ++ )
		{
			y_data[i]*=alpha;
		}
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixMatvecT
 *
 *   Performs y <- alpha * A^T * x + beta * y
 *
 *   From Van Henson's modification of hypre_CSRMatrixMatvec.
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixMatvecT(double alpha, hypre_CSRMatrix * A, hypre_Vector * x, double beta, hypre_Vector * y)
{
	double * A_data = A->data;
	int * A_i = A->i;
	int * A_j = A->j;
	int num_rows = A->num_rows;
	int num_cols = A->num_cols;
	double * x_data = x->data;
	double * y_data = y->data;
	int x_size = x->size;
	int y_size = y->size;
	int num_vectors = x->num_vectors;
	int idxstride_y = y->idxstride;
	int vecstride_y = y->vecstride;
	int idxstride_x = x->idxstride;
	int vecstride_x = x->vecstride;
	double temp;
	int i, i1, j, jv, jj, ns, ne, size, rest;
	int num_threads;
	int ierr = 0;
	/*
	---------------------------------------------------------------------
	  Check for size compatibility.  MatvecT returns ierr = 1 if
	    *  length of X doesn't equal the number of rows of A,
	    *  ierr = 2 if the length of Y doesn't equal the number of 
	    *  columns of A, and ierr = 3 if both are true.
	    *
	    *  Because temporary vectors are often used in MatvecT, none of 
	    *  these conditions terminates processing, and the ierr flag
	    *  is informational only.
	    *--------------------------------------------------------------------
	*/
	int jj_0;
	int jj_1;
	if ( ! (num_vectors==y->num_vectors))
	{
		fprintf(stderr, "hypre_assert failed: %s\n", "num_vectors == hypre_VectorNumVectors(y)");
		hypre_error_handler("csr_matvec.c", 278, 1);
	}
	;
	if (num_rows!=x_size)
	{
		ierr=1;
	}
	if (num_cols!=y_size)
	{
		ierr=2;
	}
	if ((num_rows!=x_size)&&(num_cols!=y_size))
	{
		ierr=3;
	}
	/*
	-----------------------------------------------------------------------
	 Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
	    *-----------------------------------------------------------------------
	*/
	if (alpha==0.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvecT#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+((3L*num_cols)*num_vectors)))) private(i)
		for (i=0; i<(num_cols*num_vectors); i ++ )
		{
			y_data[i]*=beta;
		}
		return ierr;
	}
	/*
	-----------------------------------------------------------------------
	 y = (beta/alpha)*y
	    *-----------------------------------------------------------------------
	*/
	temp=(beta/alpha);
	if (temp!=1.0)
	{
		if (temp==0.0)
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvecT#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+((3L*num_cols)*num_vectors)))) private(i)
			for (i=0; i<(num_cols*num_vectors); i ++ )
			{
				y_data[i]=0.0;
			}
		}
		else
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvecT#2 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+((3L*num_cols)*num_vectors)))) private(i)
			for (i=0; i<(num_cols*num_vectors); i ++ )
			{
				y_data[i]*=temp;
			}
		}
	}
	/*
	-----------------------------------------------------------------
	 y += A^T*x
	    *-----------------------------------------------------------------
	*/
	num_threads=1;
	#pragma cetus private(i, j, jj, jj_0, jj_1, jv) 
	#pragma loop name hypre_CSRMatrixMatvecT#3 
	/* #pragma cetus reduction(+: y_data[((jidxstride_y)+(jv*vecstride_y))], y_data[j])  */
	for (i=0; i<num_rows; i ++ )
	{
		if (num_vectors==1)
		{
			/* Normalized Loop */
			#pragma cetus private(j) 
			#pragma cetus lastprivate(jj_0) 
			#pragma loop name hypre_CSRMatrixMatvecT#3#0 
			/* #pragma cetus reduction(+: y_data[j])  */
			for (jj_0=0; jj_0<=((-1+A_i[1+i])+(-1*A_i[i])); jj_0 ++ )
			{
				j=A_j[jj_0+A_i[i]];
				y_data[j]+=(A_data[jj_0+A_i[i]]*x_data[i]);
			}
			jj=(jj_0+A_i[i]);
		}
		else
		{
			#pragma cetus private(j, jj, jj_1, jv) 
			#pragma loop name hypre_CSRMatrixMatvecT#3#1 
			/* #pragma cetus reduction(+: y_data[((jidxstride_y)+(jv*vecstride_y))])  */
			for (jv=0; jv<num_vectors;  ++ jv)
			{
				/* Normalized Loop */
				#pragma cetus private(j) 
				#pragma cetus lastprivate(jj_1) 
				#pragma loop name hypre_CSRMatrixMatvecT#3#1#0 
				/* #pragma cetus reduction(+: y_data[((jidxstride_y)+(jv*vecstride_y))])  */
				for (jj_1=0; jj_1<=((-1+A_i[1+i])+(-1*A_i[i])); jj_1 ++ )
				{
					j=A_j[jj_1+A_i[i]];
					y_data[(j*idxstride_y)+(jv*vecstride_y)]+=(A_data[jj_1+A_i[i]]*x_data[(i*idxstride_x)+(jv*vecstride_x)]);
				}
				jj=(jj_1+A_i[i]);
			}
		}
	}
	/*
	-----------------------------------------------------------------
	 y = alpha*y
	    *-----------------------------------------------------------------
	*/
	if (alpha!=1.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvecT#4 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+((3L*num_cols)*num_vectors)))) private(i)
		for (i=0; i<(num_cols*num_vectors); i ++ )
		{
			y_data[i]*=alpha;
		}
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixMatvec_FF
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixMatvec_FF(double alpha, hypre_CSRMatrix * A, hypre_Vector * x, double beta, hypre_Vector * y, int * CF_marker_x, int * CF_marker_y, int fpt)
{
	double * A_data = A->data;
	int * A_i = A->i;
	int * A_j = A->j;
	int num_rows = A->num_rows;
	int num_cols = A->num_cols;
	double * x_data = x->data;
	double * y_data = y->data;
	int x_size = x->size;
	int y_size = y->size;
	double temp;
	int i, jj;
	int ierr = 0;
	/*
	---------------------------------------------------------------------
	  Check for size compatibility.  Matvec returns ierr = 1 if
	    *  length of X doesn't equal the number of columns of A,
	    *  ierr = 2 if the length of Y doesn't equal the number of rows
	    *  of A, and ierr = 3 if both are true.
	    *
	    *  Because temporary vectors are often used in Matvec, none of
	    *  these conditions terminates processing, and the ierr flag
	    *  is informational only.
	    *--------------------------------------------------------------------
	*/
	int jj_0;
	if (num_cols!=x_size)
	{
		ierr=1;
	}
	if (num_rows!=y_size)
	{
		ierr=2;
	}
	if ((num_cols!=x_size)&&(num_rows!=y_size))
	{
		ierr=3;
	}
	/*
	-----------------------------------------------------------------------
	 Do (alpha == 0.0) computation - RDF: USE MACHINE EPS
	    *-----------------------------------------------------------------------
	*/
	if (alpha==0.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvec_FF#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(4L*num_rows)))) private(i)
		for (i=0; i<num_rows; i ++ )
		{
			if (CF_marker_x[i]==fpt)
			{
				y_data[i]*=beta;
			}
		}
		return ierr;
	}
	/*
	-----------------------------------------------------------------------
	 y = (beta/alpha)*y
	    *-----------------------------------------------------------------------
	*/
	temp=(beta/alpha);
	if (temp!=1.0)
	{
		if (temp==0.0)
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvec_FF#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(4L*num_rows)))) private(i)
			for (i=0; i<num_rows; i ++ )
			{
				if (CF_marker_x[i]==fpt)
				{
					y_data[i]=0.0;
				}
			}
		}
		else
		{
			#pragma cetus private(i) 
			#pragma loop name hypre_CSRMatrixMatvec_FF#2 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(4L*num_rows)))) private(i)
			for (i=0; i<num_rows; i ++ )
			{
				if (CF_marker_x[i]==fpt)
				{
					y_data[i]*=temp;
				}
			}
		}
	}
	/*
	-----------------------------------------------------------------
	 y += A*x
	    *-----------------------------------------------------------------
	*/
	#pragma cetus private(i, jj, jj_0, temp) 
	#pragma loop name hypre_CSRMatrixMatvec_FF#3 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, jj, jj_0, temp)
	for (i=0; i<num_rows; i ++ )
	{
		if (CF_marker_x[i]==fpt)
		{
			temp=y_data[i];
			/* Normalized Loop */
			#pragma cetus lastprivate(jj_0) 
			#pragma loop name hypre_CSRMatrixMatvec_FF#3#0 
			/* #pragma cetus reduction(+: temp)  */
			for (jj_0=0; jj_0<=((-1+A_i[1+i])+(-1*A_i[i])); jj_0 ++ )
			{
				if (CF_marker_y[A_j[jj_0+A_i[i]]]==fpt)
				{
					temp+=(A_data[jj_0+A_i[i]]*x_data[A_j[jj_0+A_i[i]]]);
				}
			}
			jj=(jj_0+A_i[i]);
			y_data[i]=temp;
		}
	}
	/*
	-----------------------------------------------------------------
	 y = alpha*y
	    *-----------------------------------------------------------------
	*/
	if (alpha!=1.0)
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_CSRMatrixMatvec_FF#4 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(4L*num_rows)))) private(i)
		for (i=0; i<num_rows; i ++ )
		{
			if (CF_marker_x[i]==fpt)
			{
				y_data[i]*=alpha;
			}
		}
	}
	return ierr;
}
