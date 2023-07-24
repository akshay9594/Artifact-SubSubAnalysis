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
 * Relaxation scheme
 *

*/
#include "headers.h"
/* #include "omp.h" */
/*
--------------------------------------------------------------------------
 hypre_BoomerAMGSeqRelax
 *--------------------------------------------------------------------------
*/
int hypre_BoomerAMGSeqRelax(hypre_CSRMatrix * A, hypre_Vector * f, hypre_Vector * u)
{
	double * A_diag_data = A->data;
	int * A_diag_i = A->i;
	int * A_diag_j = A->j;
	int n = A->num_rows;
	double * u_data = u->data;
	double * f_data = f->data;
	double * tmp_data;
	double res;
	int i, j;
	int ii, jj;
	int ns, ne, size, rest;
	int relax_error = 0;
	/*   int		   index, start; */
	int num_threads;
	int i_0;
	int jj_0;
	int jj_1;
	num_threads=1;
	/*
	-----------------------------------------------------------------------
	 Switch statement to direct control based on relax_type:
	    *     relax_type = 3 -> hybrid: SOR-J mix off-processor, SOR on-processor
	    *     		    with outer relaxation parameters (forward solve)
	    *-----------------------------------------------------------------------
	*/
	/*
	-----------------------------------------------------------------
	 Relax all points.
	         *-----------------------------------------------------------------
	*/
	if (1)
	{
		tmp_data=((double * )hypre_CAlloc((unsigned int)n, (unsigned int)sizeof (double)));
		num_threads=omp_get_num_threads();
		#pragma cetus private(i) 
		#pragma loop name hypre_BoomerAMGSeqRelax#0 
		#pragma cetus parallel 
		/* #pragma omp for private(i) */
		#pragma omp parallel for if((10000<(1L+(3L*n)))) private(i)
		for (i=0; i<n; i ++ )
		{
			tmp_data[i]=u_data[i];
		}
		#pragma cetus private(i, i_0, ii, j, jj, jj_0, ne, ns, res, rest, size) 
		#pragma loop name hypre_BoomerAMGSeqRelax#1 
		/* #pragma omp for private(i, ii, j, jj, ne, ns, res, rest, size) */
		for (j=0; j<num_threads; j ++ )
		{
			size=(n/num_threads);
			rest=(n-(size*num_threads));
			if (j<rest)
			{
				ns=((j*size)+j);
				ne=((((j+1)*size)+j)+1);
			}
			else
			{
				ns=((j*size)+rest);
				ne=(((j+1)*size)+rest);
			}
			/* interior points first */
			/* Normalized Loop */
			#pragma cetus private(ii, jj, jj_0, res) 
			#pragma cetus lastprivate(i_0) 
			#pragma loop name hypre_BoomerAMGSeqRelax#1#0 
			for (i_0=0; i_0<=((-1+ne)+(-1*ns)); i_0 ++ )
			{
				/*
				-----------------------------------------------------------
				 If diagonal is nonzero, relax point i; otherwise, skip it.
				                *-----------------------------------------------------------
				*/
				if (A_diag_data[A_diag_i[i_0+ns]]!=0.0)
				{
					res=f_data[i_0+ns];
					/* Normalized Loop */
					#pragma cetus private(ii) 
					#pragma cetus lastprivate(jj_0) 
					#pragma loop name hypre_BoomerAMGSeqRelax#1#0#0 
					#pragma cetus reduction(+: res) 
					#pragma cetus parallel 
					#pragma omp parallel for if((10000<((-5L+(6L*A_diag_i[((1L+i_0)+ns)]))+(-6L*A_diag_i[(i_0+ns)])))) private(ii) lastprivate(jj_0) reduction(+: res)
					for (jj_0=0; jj_0<=((-2+A_diag_i[(1+i_0)+ns])+(-1*A_diag_i[i_0+ns])); jj_0 ++ )
					{
						ii=A_diag_j[(1+jj_0)+A_diag_i[i_0+ns]];
						if ((ii>=ns)&&(ii<ne))
						{
							res-=(A_diag_data[(1+jj_0)+A_diag_i[i_0+ns]]*u_data[ii]);
						}
						else
						{
							res-=(A_diag_data[(1+jj_0)+A_diag_i[i_0+ns]]*tmp_data[ii]);
						}
					}
					jj=((1+jj_0)+A_diag_i[i_0+ns]);
					u_data[i_0+ns]=(res/A_diag_data[A_diag_i[i_0+ns]]);
				}
			}
			i=(i_0+ns);
		}
		(hypre_Free((char * )tmp_data), (tmp_data=((void * )0)));
	}
	else
	{
		/* interior points first */
		#pragma cetus private(i, ii, jj, jj_1, res) 
		#pragma loop name hypre_BoomerAMGSeqRelax#2 
		for (i=0; i<n; i ++ )
		{
			/*
			-----------------------------------------------------------
			 If diagonal is nonzero, relax point i; otherwise, skip it.
			                *-----------------------------------------------------------
			*/
			if (A_diag_data[A_diag_i[i]]!=0.0)
			{
				res=f_data[i];
				/* Normalized Loop */
				#pragma cetus private(ii) 
				#pragma cetus lastprivate(jj_1) 
				#pragma loop name hypre_BoomerAMGSeqRelax#2#0 
				#pragma cetus reduction(+: res) 
				#pragma cetus parallel 
				#pragma omp parallel for if((10000<((-3L+(4L*A_diag_i[(1L+i)]))+(-4L*A_diag_i[i])))) private(ii) lastprivate(jj_1) reduction(+: res)
				for (jj_1=0; jj_1<=((-2+A_diag_i[1+i])+(-1*A_diag_i[i])); jj_1 ++ )
				{
					ii=A_diag_j[(1+jj_1)+A_diag_i[i]];
					res-=(A_diag_data[(1+jj_1)+A_diag_i[i]]*u_data[ii]);
				}
				jj=((1+jj_1)+A_diag_i[i]);
				u_data[i]=(res/A_diag_data[A_diag_i[i]]);
			}
		}
	}
	return relax_error;
}
