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
 * $Revision: 2.12 $
 ***********************************************************************EHEADER
*/
/*

 *
 * Member functions for hypre_CSRMatrix class.
 *

*/
#include "headers.h"
/*
--------------------------------------------------------------------------
 hypre_CSRMatrixCreate
 *--------------------------------------------------------------------------
*/
hypre_CSRMatrix *hypre_CSRMatrixCreate(int num_rows, int num_cols, int num_nonzeros)
{
	hypre_CSRMatrix * matrix;
	matrix=((hypre_CSRMatrix * )hypre_CAlloc((unsigned int)1, (unsigned int)sizeof (hypre_CSRMatrix)));
	matrix->data=((void * )0);
	matrix->i=((void * )0);
	matrix->j=((void * )0);
	matrix->rownnz=((void * )0);
	matrix->num_rows=num_rows;
	matrix->num_cols=num_cols;
	matrix->num_nonzeros=num_nonzeros;
	/* set defaults */
	matrix->owns_data=1;
	printf("NNZs=%d, num_rows=%d, num_cols=%d\n", num_nonzeros, num_rows, num_cols);
	/* hypre_CSRMatrixNumRownnz(matrix) = num_rows; */
	return matrix;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixDestroy
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixDestroy(hypre_CSRMatrix * matrix)
{
	int ierr = 0;
	if (matrix)
	{
		(hypre_Free((char * )matrix->i), (matrix->i=((void * )0)));
		if (matrix->rownnz)
		{
			(hypre_Free((char * )matrix->rownnz), (matrix->rownnz=((void * )0)));
		}
		if (matrix->owns_data)
		{
			(hypre_Free((char * )matrix->data), (matrix->data=((void * )0)));
			(hypre_Free((char * )matrix->j), (matrix->j=((void * )0)));
		}
		(hypre_Free((char * )matrix), (matrix=((void * )0)));
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixInitialize
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixInitialize(hypre_CSRMatrix * matrix)
{
	int num_rows = matrix->num_rows;
	int num_nonzeros = matrix->num_nonzeros;
	/*   int  num_rownnz = hypre_CSRMatrixNumRownnz(matrix); */
	int ierr = 0;
	if (( ! matrix->data)&&num_nonzeros)
	{
		matrix->data=((double * )hypre_CAlloc((unsigned int)num_nonzeros, (unsigned int)sizeof (double)));
	}
	if ( ! matrix->i)
	{
		matrix->i=((int * )hypre_CAlloc((unsigned int)(num_rows+1), (unsigned int)sizeof (int)));
	}
	/*
	  if ( ! hypre_CSRMatrixRownnz(matrix) )
	      hypre_CSRMatrixRownnz(matrix)    = hypre_CTAlloc(int, num_rownnz);
	*/
	if (( ! matrix->j)&&num_nonzeros)
	{
		matrix->j=((int * )hypre_CAlloc((unsigned int)num_nonzeros, (unsigned int)sizeof (int)));
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixSetDataOwner
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixSetDataOwner(hypre_CSRMatrix * matrix, int owns_data)
{
	int ierr = 0;
	matrix->owns_data=owns_data;
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixSetRownnz
 *
 * function to set the substructure rownnz and num_rowsnnz inside the CSRMatrix
 * it needs the A_i substructure of CSRMatrix to find the nonzero rows.
 * It runs after the create CSR and when A_i is known..It does not check for
 * the existence of A_i or of the CSR matrix.
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixSetRownnz(hypre_CSRMatrix * matrix)
{
	int ierr = 0;
	int num_rows = matrix->num_rows;
	int * A_i = matrix->i;
	int * Arownnz;
	int i, adiag;
	int irownnz;
	irownnz=0;
	#pragma cetus private(adiag, i) 
	#pragma loop name hypre_CSRMatrixSetRownnz#0 
	#pragma cetus reduction(+: irownnz) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(5L*num_rows)))) private(adiag, i) reduction(+: irownnz)
	for (i=0; i<num_rows; i ++ )
	{
		adiag=(A_i[i+1]-A_i[i]);
		if (adiag>0)
		{
			irownnz ++ ;
		}
	}
	matrix->num_rownnz=irownnz;
	if ((irownnz==0)||(irownnz==num_rows))
	{
		matrix->rownnz=((void * )0);
	}
	else
	{
		Arownnz=((int * )hypre_CAlloc((unsigned int)irownnz, (unsigned int)sizeof (int)));
		/*  irownnz = 0; */
		/*  for (i=0; i < num_rows; i++) */
		/*  { */
			/*     adiag = A_i[i+1]-A_i[i]; */
			/*     if(adiag > 0) Arownnz[irownnz++] = i; */
		/*  } */
		matrix->rownnz=Arownnz;
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixRead
 *--------------------------------------------------------------------------
*/
hypre_CSRMatrix *hypre_CSRMatrixRead(char * file_name)
{
	hypre_CSRMatrix * matrix;
	FILE * fp;
	double * matrix_data;
	int * matrix_i;
	int * matrix_j;
	int num_rows;
	int num_nonzeros;
	int max_col = 0;
	int file_base = 1;
	int j;
	/*
	----------------------------------------------------------
	 Read in the data
	    *----------------------------------------------------------
	*/
	fp=fopen(file_name, "r");
	fscanf(fp, "%d",  & num_rows);
	matrix_i=((int * )hypre_CAlloc((unsigned int)(num_rows+1), (unsigned int)sizeof (int)));
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixRead#0 
	for (j=0; j<(num_rows+1); j ++ )
	{
		fscanf(fp, "%d",  & matrix_i[j]);
		matrix_i[j]-=file_base;
	}
	num_nonzeros=matrix_i[num_rows];
	matrix=hypre_CSRMatrixCreate(num_rows, num_rows, matrix_i[num_rows]);
	matrix->i=matrix_i;
	hypre_CSRMatrixInitialize(matrix);
	matrix_j=matrix->j;
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixRead#1 
	for (j=0; j<num_nonzeros; j ++ )
	{
		fscanf(fp, "%d",  & matrix_j[j]);
		matrix_j[j]-=file_base;
		if (matrix_j[j]>max_col)
		{
			max_col=matrix_j[j];
		}
	}
	matrix_data=matrix->data;
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixRead#2 
	for (j=0; j<matrix_i[num_rows]; j ++ )
	{
		fscanf(fp, "%le",  & matrix_data[j]);
	}
	fclose(fp);
	matrix->num_nonzeros=num_nonzeros;
	matrix->num_cols=( ++ max_col);
	return matrix;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixPrint
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixPrint(hypre_CSRMatrix * matrix, char * file_name)
{
	FILE * fp;
	double * matrix_data;
	int * matrix_i;
	int * matrix_j;
	int num_rows;
	int file_base = 1;
	int j;
	int ierr = 0;
	/*
	----------------------------------------------------------
	 Print the matrix data
	    *----------------------------------------------------------
	*/
	matrix_data=matrix->data;
	matrix_i=matrix->i;
	matrix_j=matrix->j;
	num_rows=matrix->num_rows;
	fp=fopen(file_name, "w");
	fprintf(fp, "%d\n", num_rows);
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixPrint#0 
	for (j=0; j<=num_rows; j ++ )
	{
		fprintf(fp, "%d\n", matrix_i[j]+file_base);
	}
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixPrint#1 
	for (j=0; j<matrix_i[num_rows]; j ++ )
	{
		fprintf(fp, "%d\n", matrix_j[j]+file_base);
	}
	if (matrix_data)
	{
		#pragma cetus private(j) 
		#pragma loop name hypre_CSRMatrixPrint#2 
		for (j=0; j<matrix_i[num_rows]; j ++ )
		{
			fprintf(fp, "%.14e\n", matrix_data[j]);
		}
	}
	else
	{
		fprintf(fp, "Warning: No matrix data!\n");
	}
	fclose(fp);
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixCopy:
 * copys A to B, 
 * if copy_data = 0 only the structure of A is copied to B.
 * the routine does not check if the dimensions of A and B match !!! 
 *--------------------------------------------------------------------------
*/
int hypre_CSRMatrixCopy(hypre_CSRMatrix * A, hypre_CSRMatrix * B, int copy_data)
{
	int ierr = 0;
	int num_rows = A->num_rows;
	int * A_i = A->i;
	int * A_j = A->j;
	double * A_data;
	int * B_i = B->i;
	int * B_j = B->j;
	double * B_data;
	int i, j;
	int j_0;
	int j_1;
	#pragma cetus private(i, j, j_0) 
	#pragma loop name hypre_CSRMatrixCopy#0 
	for (i=0; i<num_rows; i ++ )
	{
		B_i[i]=A_i[i];
		/* Normalized Loop */
		#pragma cetus lastprivate(j_0) 
		#pragma loop name hypre_CSRMatrixCopy#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*A_i[(1L+i)]))+(-3L*A_i[i])))) lastprivate(j_0)
		for (j_0=0; j_0<=((-1+A_i[1+i])+(-1*A_i[i])); j_0 ++ )
		{
			B_j[j_0+A_i[i]]=A_j[j_0+A_i[i]];
		}
		j=(j_0+A_i[i]);
	}
	B_i[num_rows]=A_i[num_rows];
	if (copy_data)
	{
		A_data=A->data;
		B_data=B->data;
		#pragma cetus private(i, j, j_1) 
		#pragma loop name hypre_CSRMatrixCopy#1 
		for (i=0; i<num_rows; i ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(j_1) 
			#pragma loop name hypre_CSRMatrixCopy#1#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<((1L+(3L*A_i[(1L+i)]))+(-3L*A_i[i])))) lastprivate(j_1)
			for (j_1=0; j_1<=((-1+A_i[1+i])+(-1*A_i[i])); j_1 ++ )
			{
				B_data[j_1+A_i[i]]=A_data[j_1+A_i[i]];
			}
			j=(j_1+A_i[i]);
		}
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixClone
 * Creates and returns a new copy of the argument, A.
 * Data is not copied, only structural information is reproduced.
 * Copying is a deep copy in that no pointers are copied; new arrays are
 * created where necessary.
 *--------------------------------------------------------------------------
*/
hypre_CSRMatrix *hypre_CSRMatrixClone(hypre_CSRMatrix * A)
{
	int num_rows = A->num_rows;
	int num_cols = A->num_cols;
	int num_nonzeros = A->num_nonzeros;
	hypre_CSRMatrix * B = hypre_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
	int * A_i;
	int * A_j;
	int * B_i;
	int * B_j;
	int i, j;
	hypre_CSRMatrixInitialize(B);
	A_i=A->i;
	A_j=A->j;
	B_i=B->i;
	B_j=B->j;
	#pragma cetus private(i) 
	#pragma loop name hypre_CSRMatrixClone#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(4L+(3L*num_rows)))) private(i)
	for (i=0; i<(num_rows+1);  ++ i)
	{
		B_i[i]=A_i[i];
	}
	#pragma cetus private(j) 
	#pragma loop name hypre_CSRMatrixClone#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*num_nonzeros)))) private(j)
	for (j=0; j<num_nonzeros;  ++ j)
	{
		B_j[j]=A_j[j];
	}
	B->num_rownnz=A->num_rownnz;
	if (A->rownnz)
	{
		hypre_CSRMatrixSetRownnz(B);
	}
	return B;
}

/*
--------------------------------------------------------------------------
 hypre_CSRMatrixUnion
 * Creates and returns a matrix whose elements are the union of those of A and B.
 * Data is not computed, only structural information is created.
 * A and B must have the same numbers of rows.
 * Nothing is done about Rownnz.
 *
 * If col_map_offd_A and col_map_offd_B are zero, A and B are expected to have
 * the same column indexing.  Otherwise, col_map_offd_A, col_map_offd_B should
 * be the arrays of that name from two ParCSRMatrices of which A and B are the
 * offd blocks.
 *
 * The algorithm can be expected to have reasonable efficiency only for very
 * sparse matrices (many rows, few nonzeros per row).
 * The nonzeros of a computed row are NOT necessarily in any particular order.
 *--------------------------------------------------------------------------
*/
hypre_CSRMatrix *hypre_CSRMatrixUnion(hypre_CSRMatrix * A, hypre_CSRMatrix * B, int * col_map_offd_A, int * col_map_offd_B, int * * col_map_offd_C)
{
	int num_rows = A->num_rows;
	int num_cols_A = A->num_cols;
	int num_cols_B = B->num_cols;
	int num_cols;
	int num_nonzeros;
	int * A_i = A->i;
	int * A_j = A->j;
	int * B_i = B->i;
	int * B_j = B->j;
	int * C_i;
	int * C_j;
	int * jC = (void * )0;
	int i, jA, jB, jBg;
	int ma, mb, mc, ma_min, ma_max, match;
	hypre_CSRMatrix * C;
	int mb_0;
	int ma_0;
	int ma_1;
	int mb_1;
	int ma_2;
	if ( ! (num_rows==B->num_rows))
	{
		fprintf(stderr, "hypre_assert failed: %s\n", "num_rows == hypre_CSRMatrixNumRows(B)");
		hypre_error_handler("csr_matrix.c", 426, 1);
	}
	;
	if (col_map_offd_B)
	{
		if ( ! col_map_offd_A)
		{
			fprintf(stderr, "hypre_assert failed: %s\n", "col_map_offd_A");
			hypre_error_handler("csr_matrix.c", 427, 1);
		}
	}
	;
	if (col_map_offd_A)
	{
		if ( ! col_map_offd_B)
		{
			fprintf(stderr, "hypre_assert failed: %s\n", "col_map_offd_B");
			hypre_error_handler("csr_matrix.c", 428, 1);
		}
	}
	;
	/* ==== First, go through the columns of A and B to count the columns of C. */
	if (col_map_offd_A==0)
	{
		/*
		The matrices are diagonal blocks.
		         Normally num_cols_A==num_cols_B, col_starts is the same, etc.
		     
		*/
		num_cols=((num_cols_A<num_cols_B) ? num_cols_B : num_cols_A);
	}
	else
	{
		/* The matrices are offdiagonal blocks. */
		jC=((int * )hypre_CAlloc((unsigned int)num_cols_B, (unsigned int)sizeof (int)));
		num_cols=num_cols_A;
		/* initialization; we'll compute the actual value */
		#pragma cetus private(jB, jBg, ma, match) 
		#pragma loop name hypre_CSRMatrixUnion#0 
		for (jB=0; jB<num_cols_B;  ++ jB)
		{
			match=0;
			jBg=col_map_offd_B[jB];
			#pragma cetus private(ma) 
			#pragma cetus lastprivate(match) 
			#pragma loop name hypre_CSRMatrixUnion#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(4L*num_cols_A)))) private(ma) lastprivate(match)
			for (ma=0; ma<num_cols_A;  ++ ma)
			{
				if (col_map_offd_A[ma]==jBg)
				{
					match=1;
				}
			}
			if (match==0)
			{
				jC[jB]=num_cols;
				 ++ num_cols;
			}
		}
	}
	/*
	==== If we're working on a ParCSRMatrix's offd block,
	      make and load col_map_offd_C
	*/
	if (col_map_offd_A)
	{
		( * col_map_offd_C)=((int * )hypre_CAlloc((unsigned int)num_cols, (unsigned int)sizeof (int)));
		#pragma cetus private(jA) 
		#pragma loop name hypre_CSRMatrixUnion#1 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(3L*num_cols_A)))) private(jA)
		for (jA=0; jA<num_cols_A;  ++ jA)
		{
			( * col_map_offd_C)[jA]=col_map_offd_A[jA];
		}
		#pragma cetus private(jB, jBg, ma, match) 
		#pragma loop name hypre_CSRMatrixUnion#2 
		for (jB=0; jB<num_cols_B;  ++ jB)
		{
			match=0;
			jBg=col_map_offd_B[jB];
			#pragma cetus private(ma) 
			#pragma cetus lastprivate(match) 
			#pragma loop name hypre_CSRMatrixUnion#2#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(4L*num_cols_A)))) private(ma) lastprivate(match)
			for (ma=0; ma<num_cols_A;  ++ ma)
			{
				if (col_map_offd_A[ma]==jBg)
				{
					match=1;
				}
			}
			if (match==0)
			{
				( * col_map_offd_C)[jC[jB]]=jBg;
			}
		}
	}
	/*
	==== The first run through A and B is to count the number of nonzero elements,
	      without double-counting duplicates.  Then we can create C.
	*/
	num_nonzeros=A->num_nonzeros;
	#pragma cetus private(i, jA, jB, ma, ma_0, ma_max, ma_min, match, mb, mb_0) 
	#pragma loop name hypre_CSRMatrixUnion#3 
	#pragma cetus reduction(+: num_nonzeros) 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, jA, jB, ma, ma_0, ma_max, ma_min, match, mb, mb_0) reduction(+: num_nonzeros)
	for (i=0; i<num_rows;  ++ i)
	{
		ma_min=A_i[i];
		ma_max=A_i[i+1];
		/* Normalized Loop */
		#pragma cetus private(jA, jB, ma, ma_0, match) 
		#pragma cetus lastprivate(mb_0) 
		#pragma loop name hypre_CSRMatrixUnion#3#0 
		/* #pragma cetus reduction(+: num_nonzeros)  */
		for (mb_0=0; mb_0<=((-1+B_i[1+i])+(-1*B_i[i])); mb_0 ++ )
		{
			jB=B_j[mb_0+B_i[i]];
			if (col_map_offd_B)
			{
				jB=col_map_offd_B[jB];
			}
			match=0;
			/* Normalized Loop */
			#pragma loop name hypre_CSRMatrixUnion#3#0#0 
			for (ma_0=0; ma_0<=((-1+ma_max)+(-1*ma_min)); ma_0 ++ )
			{
				jA=A_j[ma_0+ma_min];
				if (col_map_offd_A)
				{
					jA=col_map_offd_A[jA];
				}
				if (jB==jA)
				{
					match=1;
					if ((ma_0+ma_min)==ma_min)
					{
						 ++ ma_min;
					}
					break;
				}
			}
			ma=(ma_0+ma_min);
			if (match==0)
			{
				 ++ num_nonzeros;
			}
		}
		mb=(mb_0+B_i[i]);
	}
	C=hypre_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
	hypre_CSRMatrixInitialize(C);
	/*
	==== The second run through A and B is to pick out the column numbers
	      for each row, and put them in C.
	*/
	C_i=C->i;
	C_i[0]=0;
	C_j=C->j;
	mc=0;
	#pragma cetus private(i, jA, jB, ma, ma_1, ma_2, ma_max, ma_min, match, mb, mb_1) 
	#pragma loop name hypre_CSRMatrixUnion#4 
	for (i=0; i<num_rows;  ++ i)
	{
		ma_min=A_i[i];
		ma_max=A_i[i+1];
		/* Normalized Loop */
		#pragma cetus lastprivate(ma_1) 
		#pragma loop name hypre_CSRMatrixUnion#4#0 
		for (ma_1=0; ma_1<=((-1+ma_max)+(-1*ma_min)); ma_1 ++ )
		{
			C_j[mc]=A_j[ma_1+ma_min];
			 ++ mc;
		}
		ma=(ma_1+ma_min);
		/* Normalized Loop */
		#pragma cetus private(jA, jB, ma, ma_2, match) 
		#pragma cetus lastprivate(mb_1) 
		#pragma loop name hypre_CSRMatrixUnion#4#1 
		for (mb_1=0; mb_1<=((-1+B_i[1+i])+(-1*B_i[i])); mb_1 ++ )
		{
			jB=B_j[mb_1+B_i[i]];
			if (col_map_offd_B)
			{
				jB=col_map_offd_B[jB];
			}
			match=0;
			/* Normalized Loop */
			#pragma loop name hypre_CSRMatrixUnion#4#1#0 
			for (ma_2=0; ma_2<=((-1+ma_max)+(-1*ma_min)); ma_2 ++ )
			{
				jA=A_j[ma_2+ma_min];
				if (col_map_offd_A)
				{
					jA=col_map_offd_A[jA];
				}
				if (jB==jA)
				{
					match=1;
					if ((ma_2+ma_min)==ma_min)
					{
						 ++ ma_min;
					}
					break;
				}
			}
			ma=(ma_2+ma_min);
			if (match==0)
			{
				C_j[mc]=jC[B_j[mb_1+B_i[i]]];
				 ++ mc;
			}
		}
		mb=(mb_1+B_i[i]);
		C_i[i+1]=mc;
	}
	if ( ! (mc==num_nonzeros))
	{
		fprintf(stderr, "hypre_assert failed: %s\n", "mc == num_nonzeros");
		hypre_error_handler("csr_matrix.c", 550, 1);
	}
	;
	if (jC)
	{
		(hypre_Free((char * )jC), (jC=((void * )0)));
	}
	return C;
}
