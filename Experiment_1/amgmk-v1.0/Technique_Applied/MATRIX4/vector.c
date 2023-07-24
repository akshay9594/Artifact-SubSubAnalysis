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
 * $Revision: 2.8 $
 ***********************************************************************EHEADER
*/
/*

 *
 * Member functions for hypre_Vector class.
 *

*/
#include "headers.h"
#include <assert.h>
/*
--------------------------------------------------------------------------
 hypre_SeqVectorCreate
 *--------------------------------------------------------------------------
*/
hypre_Vector *hypre_SeqVectorCreate(int size)
{
	hypre_Vector * vector;
	vector=((hypre_Vector * )hypre_CAlloc((unsigned int)1, (unsigned int)sizeof (hypre_Vector)));
	vector->data=((void * )0);
	vector->size=size;
	vector->num_vectors=1;
	vector->multivec_storage_method=0;
	/* set defaults */
	vector->owns_data=1;
	return vector;
}

/*
--------------------------------------------------------------------------
 hypre_SeqMultiVectorCreate
 *--------------------------------------------------------------------------
*/
hypre_Vector *hypre_SeqMultiVectorCreate(int size, int num_vectors)
{
	hypre_Vector * vector = hypre_SeqVectorCreate(size);
	vector->num_vectors=num_vectors;
	return vector;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorDestroy
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorDestroy(hypre_Vector * vector)
{
	int ierr = 0;
	if (vector)
	{
		if (vector->owns_data)
		{
			(hypre_Free((char * )vector->data), (vector->data=((void * )0)));
		}
		(hypre_Free((char * )vector), (vector=((void * )0)));
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorInitialize
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorInitialize(hypre_Vector * vector)
{
	int size = vector->size;
	int ierr = 0;
	int num_vectors = vector->num_vectors;
	int multivec_storage_method = vector->multivec_storage_method;
	if ( ! vector->data)
	{
		vector->data=((double * )hypre_CAlloc((unsigned int)(num_vectors*size), (unsigned int)sizeof (double)));
	}
	if (multivec_storage_method==0)
	{
		vector->vecstride=size;
		vector->idxstride=1;
	}
	else
	{
		if (multivec_storage_method==1)
		{
			vector->vecstride=1;
			vector->idxstride=num_vectors;
		}
		else
		{
			 ++ ierr;
		}
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorSetDataOwner
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorSetDataOwner(hypre_Vector * vector, int owns_data)
{
	int ierr = 0;
	vector->owns_data=owns_data;
	return ierr;
}

/*
--------------------------------------------------------------------------
 ReadVector
 *--------------------------------------------------------------------------
*/
hypre_Vector *hypre_SeqVectorRead(char * file_name)
{
	hypre_Vector * vector;
	FILE * fp;
	double * data;
	int size;
	int j;
	/*
	----------------------------------------------------------
	 Read in the data
	    *----------------------------------------------------------
	*/
	fp=fopen(file_name, "r");
	fscanf(fp, "%d",  & size);
	vector=hypre_SeqVectorCreate(size);
	hypre_SeqVectorInitialize(vector);
	data=vector->data;
	#pragma cetus private(j) 
	#pragma loop name hypre_SeqVectorRead#0 
	for (j=0; j<size; j ++ )
	{
		fscanf(fp, "%le",  & data[j]);
	}
	fclose(fp);
	/* multivector code not written yet >>> */
	if ( ! (vector->num_vectors==1))
	{
		fprintf(stderr, "hypre_assert failed: %s\n", "hypre_VectorNumVectors(vector) == 1");
		hypre_error_handler("vector.c", 177, 1);
	}
	;
	return vector;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorPrint
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorPrint(hypre_Vector * vector, char * file_name)
{
	FILE * fp;
	double * data;
	int size, num_vectors, vecstride, idxstride;
	int i, j;
	int ierr = 0;
	num_vectors=vector->num_vectors;
	vecstride=vector->vecstride;
	idxstride=vector->idxstride;
	/*
	----------------------------------------------------------
	 Print in the data
	    *----------------------------------------------------------
	*/
	data=vector->data;
	size=vector->size;
	fp=fopen(file_name, "w");
	if (vector->num_vectors==1)
	{
		fprintf(fp, "%d\n", size);
	}
	else
	{
		fprintf(fp, "%d vectors of size %d\n", num_vectors, size);
	}
	if (num_vectors>1)
	{
		#pragma cetus private(i, j) 
		#pragma loop name hypre_SeqVectorPrint#0 
		for (j=0; j<num_vectors;  ++ j)
		{
			fprintf(fp, "vector %d\n", j);
			#pragma cetus private(i) 
			#pragma loop name hypre_SeqVectorPrint#0#0 
			for (i=0; i<size; i ++ )
			{
				fprintf(fp, "%.14e\n", data[(j*vecstride)+(i*idxstride)]);
			}
		}
	}
	else
	{
		#pragma cetus private(i) 
		#pragma loop name hypre_SeqVectorPrint#1 
		for (i=0; i<size; i ++ )
		{
			fprintf(fp, "%.14e\n", data[i]);
		}
	}
	fclose(fp);
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorSetConstantValues
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorSetConstantValues(hypre_Vector * v, double value)
{
	double * vector_data = v->data;
	int size = v->size;
	int i;
	int ierr = 0;
	size*=v->num_vectors;
	#pragma cetus private(i) 
	#pragma loop name hypre_SeqVectorSetConstantValues#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i)
	for (i=0; i<size; i ++ )
	{
		vector_data[i]=value;
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorCopy
 * copies data from x to y
 * y should have already been initialized at the same size as x
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorCopy(hypre_Vector * x, hypre_Vector * y)
{
	double * x_data = x->data;
	double * y_data = y->data;
	int size = x->size;
	int i;
	int ierr = 0;
	size*=x->num_vectors;
	#pragma cetus private(i) 
	#pragma loop name hypre_SeqVectorCopy#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i)
	for (i=0; i<size; i ++ )
	{
		y_data[i]=x_data[i];
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorCloneDeep
 * Returns a complete copy of x - a deep copy, with its own copy of the data.
 *--------------------------------------------------------------------------
*/
hypre_Vector *hypre_SeqVectorCloneDeep(hypre_Vector * x)
{
	int size = x->size;
	int num_vectors = x->num_vectors;
	hypre_Vector * y = hypre_SeqMultiVectorCreate(size, num_vectors);
	y->multivec_storage_method=x->multivec_storage_method;
	y->vecstride=x->vecstride;
	y->idxstride=x->idxstride;
	hypre_SeqVectorInitialize(y);
	hypre_SeqVectorCopy(x, y);
	return y;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorCloneShallow
 * Returns a complete copy of x - a shallow copy, pointing the data of x
 *--------------------------------------------------------------------------
*/
hypre_Vector *hypre_SeqVectorCloneShallow(hypre_Vector * x)
{
	int size = x->size;
	int num_vectors = x->num_vectors;
	hypre_Vector * y = hypre_SeqMultiVectorCreate(size, num_vectors);
	y->multivec_storage_method=x->multivec_storage_method;
	y->vecstride=x->vecstride;
	y->idxstride=x->idxstride;
	y->data=x->data;
	hypre_SeqVectorSetDataOwner(y, 0);
	hypre_SeqVectorInitialize(y);
	return y;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorScale
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorScale(double alpha, hypre_Vector * y)
{
	double * y_data = y->data;
	int size = y->size;
	int i;
	int ierr = 0;
	size*=y->num_vectors;
	#pragma cetus private(i) 
	#pragma loop name hypre_SeqVectorScale#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i)
	for (i=0; i<size; i ++ )
	{
		y_data[i]*=alpha;
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorAxpy
 *--------------------------------------------------------------------------
*/
int hypre_SeqVectorAxpy(double alpha, hypre_Vector * x, hypre_Vector * y)
{
	double * x_data = x->data;
	double * y_data = y->data;
	int size = x->size;
	int i;
	int ierr = 0;
	size*=x->num_vectors;
	#pragma cetus private(i) 
	#pragma loop name hypre_SeqVectorAxpy#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i)
	for (i=0; i<size; i ++ )
	{
		y_data[i]+=(alpha*x_data[i]);
	}
	return ierr;
}

/*
--------------------------------------------------------------------------
 hypre_SeqVectorInnerProd
 *--------------------------------------------------------------------------
*/
double hypre_SeqVectorInnerProd(hypre_Vector * x, hypre_Vector * y)
{
	double * x_data = x->data;
	double * y_data = y->data;
	int size = x->size;
	int i;
	double result = 0.0;
	size*=x->num_vectors;
	#pragma cetus private(i) 
	#pragma loop name hypre_SeqVectorInnerProd#0 
	#pragma cetus reduction(+: result) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i) reduction(+: result)
	for (i=0; i<size; i ++ )
	{
		result+=(y_data[i]*x_data[i]);
	}
	return result;
}

/*
--------------------------------------------------------------------------
 hypre_VectorSumElts:
 * Returns the sum of all vector elements.
 *--------------------------------------------------------------------------
*/
double hypre_VectorSumElts(hypre_Vector * vector)
{
	double sum = 0;
	double * data = vector->data;
	int size = vector->size;
	int i;
	#pragma cetus private(i) 
	#pragma loop name hypre_VectorSumElts#0 
	#pragma cetus reduction(+: sum) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*size)))) private(i) reduction(+: sum)
	for (i=0; i<size;  ++ i)
	{
		sum+=data[i];
	}
	return sum;
}
