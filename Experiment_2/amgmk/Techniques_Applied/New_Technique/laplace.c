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
#include "headers.h"
/*
--------------------------------------------------------------------------
 hypre_GenerateLaplacian
 *--------------------------------------------------------------------------
*/
hypre_CSRMatrix *GenerateSeqLaplacian(int nx, int ny, int nz, double * value, hypre_Vector * * rhs_ptr, hypre_Vector * * x_ptr, hypre_Vector * * sol_ptr)
{
	hypre_CSRMatrix * A;
	hypre_Vector * rhs;
	hypre_Vector * x;
	hypre_Vector * sol;
	double * rhs_data;
	double * x_data;
	double * sol_data;
	int * A_i;
	int * A_j;
	double * A_data;
	int ix, iy, iz, im;
	int cnt;
	int row_index;
	int i, j;
	int grid_size;
	int nm = 5;
	int upper_bound = (nx*ny)*nz;
	int j_0;
	grid_size=(((nx*ny)*nz)*nm);
	A_i=((int * )hypre_CAlloc((unsigned int)(grid_size+1), (unsigned int)sizeof (int)));
	rhs_data=((double * )hypre_CAlloc((unsigned int)grid_size, (unsigned int)sizeof (double)));
	x_data=((double * )hypre_CAlloc((unsigned int)grid_size, (unsigned int)sizeof (double)));
	sol_data=((double * )hypre_CAlloc((unsigned int)grid_size, (unsigned int)sizeof (double)));
	#pragma cetus private(i) 
	#pragma loop name GenerateSeqLaplacian#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(5L*grid_size)))) private(i)
	for (i=0; i<grid_size; i ++ )
	{
		x_data[i]=0.0;
		sol_data[i]=0.0;
		rhs_data[i]=1.0;
	}
	cnt=0;
	#pragma cetus private(im, ix, iy, iz) 
	#pragma loop name GenerateSeqLaplacian#1 
	for (iz=0; iz<nz; iz ++ )
	{
		#pragma cetus private(im, ix, iy) 
		#pragma loop name GenerateSeqLaplacian#1#0 
		for (iy=0; iy<ny; iy ++ )
		{
			#pragma cetus private(im, ix) 
			#pragma loop name GenerateSeqLaplacian#1#0#0 
			for (ix=0; ix<nx; ix ++ )
			{
				#pragma cetus private(im) 
				#pragma loop name GenerateSeqLaplacian#1#0#0#0 
				for (im=0; im<nm; im ++ )
				{
					A_i[cnt ++ ]=0;
				}
			}
		}
	}
	cnt=1;
	A_i[0]=0;
	if ((-1+nx)>=0)
	{
		#pragma cetus private(ix, iy, iz) 
		#pragma loop name GenerateSeqLaplacian#2 
		for (iz=0; iz<nz; iz ++ )
		{
			#pragma cetus private(ix, iy) 
			#pragma loop name GenerateSeqLaplacian#2#0 
			for (iy=0; iy<ny; iy ++ )
			{
				#pragma cetus private(ix) 
				#pragma loop name GenerateSeqLaplacian#2#0#0 
				for (ix=0; ix<nx; ix ++ )
				{
					A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)]=A_i[(((1+ix)+(iy*nx))+((iz*nx)*ny))-1];
					A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					if (iz)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
					if (iy)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
					if (ix)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
					if ((ix+1)<nx)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
					if ((iy+1)<ny)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
					if ((iz+1)<nz)
					{
						A_i[((1+ix)+(iy*nx))+((iz*nx)*ny)] ++ ;
					}
				}
			}
		}
		cnt+=((nx*ny)*nz);
	}
	else
	{
		#pragma cetus private(ix, iy, iz) 
		#pragma loop name GenerateSeqLaplacian#3 
		for (iz=0; iz<nz; iz ++ )
		{
			#pragma cetus private(ix, iy) 
			#pragma loop name GenerateSeqLaplacian#3#0 
			for (iy=0; iy<ny; iy ++ )
			{
				#pragma cetus private(ix) 
				#pragma loop name GenerateSeqLaplacian#3#0#0 
				for (ix=0; ix<nx; ix ++ )
				{
					A_i[cnt]=A_i[cnt-1];
					A_i[cnt] ++ ;
					if (iz)
					{
						A_i[cnt] ++ ;
					}
					if (iy)
					{
						A_i[cnt] ++ ;
					}
					if (ix)
					{
						A_i[cnt] ++ ;
					}
					if ((ix+1)<nx)
					{
						A_i[cnt] ++ ;
					}
					if ((iy+1)<ny)
					{
						A_i[cnt] ++ ;
					}
					if ((iz+1)<nz)
					{
						A_i[cnt] ++ ;
					}
					cnt ++ ;
				}
			}
		}
	}
	/* printf("A_i generated\n"); */
	A_j=((int * )hypre_CAlloc((unsigned int)A_i[upper_bound], (unsigned int)sizeof (int)));
	A_data=((double * )hypre_CAlloc((unsigned int)A_i[upper_bound], (unsigned int)sizeof (double)));
	/* printf("bound=%d\n", A_i[upper_bound]); */
	row_index=0;
	cnt=0;
	if ((-1+nx)>=0)
	{
		#pragma cetus private(ix, iy, iz) 
		#pragma loop name GenerateSeqLaplacian#4 
		for (iz=0; iz<nz; iz ++ )
		{
			#pragma cetus private(ix, iy) 
			#pragma loop name GenerateSeqLaplacian#4#0 
			for (iy=0; iy<ny; iy ++ )
			{
				#pragma cetus private(ix) 
				#pragma loop name GenerateSeqLaplacian#4#0#0 
				for (ix=0; ix<nx; ix ++ )
				{
					A_j[cnt]=((ix+(iy*nx))+((iz*nx)*ny));
					A_data[cnt ++ ]=value[0];
					if (iz)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))-(nx*ny));
						A_data[cnt ++ ]=value[3];
					}
					if (iy)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))-nx);
						A_data[cnt ++ ]=value[2];
					}
					if (ix)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))-1);
						A_data[cnt ++ ]=value[1];
					}
					if ((ix+1)<nx)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))+1);
						A_data[cnt ++ ]=value[1];
					}
					if ((iy+1)<ny)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))+nx);
						A_data[cnt ++ ]=value[2];
					}
					if ((iz+1)<nz)
					{
						A_j[cnt]=(((ix+(iy*nx))+((iz*nx)*ny))+(nx*ny));
						A_data[cnt ++ ]=value[3];
					}
				}
			}
		}
		row_index+=((nx*ny)*nz);
	}
	else
	{
		#pragma cetus private(ix, iy, iz) 
		#pragma loop name GenerateSeqLaplacian#5 
		for (iz=0; iz<nz; iz ++ )
		{
			#pragma cetus private(ix, iy) 
			#pragma loop name GenerateSeqLaplacian#5#0 
			for (iy=0; iy<ny; iy ++ )
			{
				#pragma cetus private(ix) 
				#pragma loop name GenerateSeqLaplacian#5#0#0 
				for (ix=0; ix<nx; ix ++ )
				{
					A_j[cnt]=row_index;
					A_data[cnt ++ ]=value[0];
					if (iz)
					{
						A_j[cnt]=(row_index-(nx*ny));
						A_data[cnt ++ ]=value[3];
					}
					if (iy)
					{
						A_j[cnt]=(row_index-nx);
						A_data[cnt ++ ]=value[2];
					}
					if (ix)
					{
						A_j[cnt]=(row_index-1);
						A_data[cnt ++ ]=value[1];
					}
					if ((ix+1)<nx)
					{
						A_j[cnt]=(row_index+1);
						A_data[cnt ++ ]=value[1];
					}
					if ((iy+1)<ny)
					{
						A_j[cnt]=(row_index+nx);
						A_data[cnt ++ ]=value[2];
					}
					if ((iz+1)<nz)
					{
						A_j[cnt]=(row_index+(nx*ny));
						A_data[cnt ++ ]=value[3];
					}
					row_index ++ ;
				}
			}
		}
	}
	/* printf("A_j generated\n"); */
	A=hypre_CSRMatrixCreate(grid_size, grid_size, A_i[upper_bound]);
	rhs=hypre_SeqVectorCreate(grid_size);
	rhs->data=rhs_data;
	x=hypre_SeqVectorCreate(grid_size);
	x->data=x_data;
	#pragma cetus private(i, j, j_0) 
	#pragma loop name GenerateSeqLaplacian#6 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j, j_0)
	for (i=0; i<grid_size; i ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(j_0) 
		#pragma loop name GenerateSeqLaplacian#6#0 
		/* #pragma cetus reduction(+: sol_data[i])  */
		for (j_0=0; j_0<=((-1+A_i[1+i])+(-1*A_i[i])); j_0 ++ )
		{
			sol_data[i]+=A_data[j_0+A_i[i]];
		}
		j=(j_0+A_i[i]);
	}
	sol=hypre_SeqVectorCreate(grid_size);
	sol->data=sol_data;
	A->i=A_i;
	A->j=A_j;
	A->data=A_data;
	( * rhs_ptr)=rhs;
	( * x_ptr)=x;
	( * sol_ptr)=sol;
	return A;
}
