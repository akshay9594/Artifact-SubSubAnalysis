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
/* ------------------------------------------------------------------------- */
/*                                                                          */
/*  This benchmark is a serial C version of the NPB UA code. This C         */
/*  version is developed by the Center for Manycore Programming at Seoul    */
/*  National University and derived from the serial Fortran versions in     */
/*  "NPB3.3-SER" developed by NAS.                                          */
/*                                                                          */
/*  Permission to use, copy, distribute and modify this software for any    */
/*  purpose with or without fee is hereby granted. This software is         */
/*  provided "as is" without express or implied warranty.                   */
/*                                                                          */
/*  Information on NPB 3.3, including the technical report, the original    */
/*  specifications, source code, results and information on how to submit   */
/*  new results, is available at:                                           */
/*                                                                          */
/*           http:www.nas.nasa.govSoftware/NPB/                          */
/*                                                                          */
/*  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr   */
/*                                                                          */
/*          Center for Manycore Programming                                 */
/*          School of Computer Science and Engineering                      */
/*          Seoul National University                                       */
/*          Seoul 151-744, Korea                                            */
/*                                                                          */
/*          E-mail:  cmp@aces.snu.ac.kr                                     */
/*                                                                          */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,     */
/*          and Jaejin Lee                                                  */
/* ------------------------------------------------------------------------- */
#include <assert.h>
#include "header.h"
static void pc_corner(int imor);
static void com_dpc(int iside, int iel, int enumber, int n, int isize);
/* ------------------------------------------------------------------ */
/* Generate diagonal preconditioner for CG. */
/* Preconditioner computed in this subroutine is correct only */
/* for collocation point in element interior, on conforming face */
/* interior and conforming edge. */
/* ------------------------------------------------------------------ */
void setuppc()
{
	double dxtm1_2[5][5], rdtime;
	int ie, k, i, j, q, isize;
	#pragma cetus private(i, j) 
	#pragma loop name setuppc#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, j)
	*/
	for (j=0; j<5; j ++ )
	{
		#pragma cetus private(i) 
		#pragma loop name setuppc#0#0 
		for (i=0; i<5; i ++ )
		{
			dxtm1_2[j][i]=(dxtm1[j][i]*dxtm1[j][i]);
		}
	}
	rdtime=(1.0/dtime);
	#pragma cetus private(i, ie, isize, j, k, q) 
	#pragma loop name setuppc#1 
	for (ie=0; ie<nelt; ie ++ )
	{
		r_init(dpcelm[ie][0][0], (5*5)*5, 0.0);
		isize=size_e[ie];
		#pragma cetus private(i, j, k, q) 
		#pragma loop name setuppc#1#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, j, k, q)
		*/
		for (k=0; k<5; k ++ )
		{
			#pragma cetus private(i, j, q) 
			#pragma loop name setuppc#1#0#0 
			for (j=0; j<5; j ++ )
			{
				#pragma cetus private(i, q) 
				#pragma loop name setuppc#1#0#0#0 
				for (i=0; i<5; i ++ )
				{
					#pragma cetus private(q) 
					#pragma loop name setuppc#1#0#0#0#0 
					for (q=0; q<5; q ++ )
					{
						dpcelm[ie][k][j][i]=(((dpcelm[ie][k][j][i]+(g1m1_s[isize][k][j][q]*dxtm1_2[q][i]))+(g1m1_s[isize][k][q][i]*dxtm1_2[q][j]))+(g1m1_s[isize][q][j][i]*dxtm1_2[q][k]));
					}
					dpcelm[ie][k][j][i]=((0.005*dpcelm[ie][k][j][i])+(rdtime*bm1_s[isize][k][j][i]));
				}
			}
		}
	}
	/* do the stiffness summation */
	dssum();
	/* take inverse. */
	reciprocal((double * )dpcelm, ntot);
	/* compute preconditioner on mortar points. NOTE:  dpcmor for  */
	/* nonconforming cases will be corrected in subroutine setpcmo  */
	#pragma cetus private(i) 
	#pragma loop name setuppc#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*nmor)))) private(i)
	for (i=0; i<nmor; i ++ )
	{
		dpcmor[i]=(1.0/dpcmor[i]);
	}
}

/* -------------------------------------------------------------- */
/* pre-compute elemental contribution to preconditioner   */
/* for all situations */
/* -------------------------------------------------------------- */
void setpcmo_pre()
{
	int element_size, i, j, ii, jj, col;
	double p[5][5][5], p0[5][5][5], mtemp[5][5];
	double temp[5][5][5], temp1[5][5], tmp[5][5], tig[5];
	/* corners on face of type 3  */
	int i_0;
	int j_0;
	int i_1;
	int j_1;
	int i_2;
	int i_3;
	int j_2;
	int j_3;
	int i_4;
	int i_5;
	int j_4;
	int i_6;
	int i_7;
	r_init((double * )tcpre, 5*5, 0.0);
	r_init((double * )tmp, 5*5, 0.0);
	r_init(tig, 5, 0.0);
	tig[0]=1.0;
	tmp[0][0]=1.0;
	/* tcpre results from mapping a unit spike field (unity at  */
	/* collocation point (0,0), zero elsewhere) on an entire element */
	/* face to the (0,0) segment of a nonconforming face */
	/* Normalized Loop */
	#pragma cetus private(j) 
	#pragma cetus lastprivate(i_0) 
	#pragma loop name setpcmo_pre#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(j) lastprivate(i_0)
	*/
	for (i_0=0; i_0<=2; i_0 ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name setpcmo_pre#0#0 
		for (j=0; j<5; j ++ )
		{
			tmp[0][1+i_0]=(tmp[0][1+i_0]+(qbnew[0][j][(1+i_0)-1]*tig[j]));
		}
	}
	i=(1+i_0);
	#pragma cetus private(col, i, j, j_0) 
	#pragma loop name setpcmo_pre#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, j, j_0)
	*/
	for (col=0; col<5; col ++ )
	{
		tcpre[0][col]=tmp[0][col];
		/* Normalized Loop */
		#pragma cetus private(i) 
		#pragma cetus lastprivate(j_0) 
		#pragma loop name setpcmo_pre#1#0 
		for (j_0=0; j_0<=2; j_0 ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#1#0#0 
			for (i=0; i<5; i ++ )
			{
				tcpre[1+j_0][col]=(tcpre[1+j_0][col]+(qbnew[0][i][(1+j_0)-1]*tmp[i][col]));
			}
		}
		j=(1+j_0);
	}
	#pragma cetus private(element_size, i, i_1, i_2, i_3, i_4, i_5, i_6, i_7, ii, j, j_1, j_2, j_3, j_4, jj) 
	#pragma loop name setpcmo_pre#2 
	for (element_size=0; element_size<10; element_size ++ )
	{
		/* for conforming cases */
		/* pcmor_c[element_size][j][i] records the intermediate value  */
		/* (preconditioner=1pcmor_c) of the preconditor on collocation  */
		/* point (i,j) on a conforming face of an element of size  */
		/* element_size. */
		#pragma cetus private(i, i_1, j) 
		#pragma loop name setpcmo_pre#2#0 
		for (j=0; j<((5/2)+1); j ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_1) 
			#pragma loop name setpcmo_pre#2#0#0 
			for (i_1=0; i_1<=(2+(-1*j)); i_1 ++ )
			{
				r_init((double * )p, (5*5)*5, 0.0);
				p[0][j][i_1+j]=1.0;
				laplacian(temp, p, element_size);
				pcmor_c[element_size][j][i_1+j]=temp[0][j][i_1+j];
				pcmor_c[element_size][j][(5-1)-(i_1+j)]=temp[0][j][i_1+j];
				pcmor_c[element_size][i_1+j][j]=temp[0][j][i_1+j];
				pcmor_c[element_size][i_1+j][(5-1)-j]=temp[0][j][i_1+j];
				pcmor_c[element_size][(5-1)-(i_1+j)][j]=temp[0][j][i_1+j];
				pcmor_c[element_size][(5-1)-(i_1+j)][(5-1)-j]=temp[0][j][i_1+j];
				pcmor_c[element_size][(5-1)-j][i_1+j]=temp[0][j][i_1+j];
				pcmor_c[element_size][(5-1)-j][(5-1)-(i_1+j)]=temp[0][j][i_1+j];
			}
			i=(i_1+j);
		}
		/* for nonconforming cases  */
		/* nonconforming face interior */
		/* pcmor_nc1[element_size][jj][ii][j][i] records the intermediate  */
		/* preconditioner value on collocation point (i,j) on mortar  */
		/* (ii,jj)  on a nonconforming face of an element of size element_ */
		/* size */
		/* Normalized Loop */
		#pragma cetus private(i, i_2) 
		#pragma cetus lastprivate(j_1) 
		#pragma loop name setpcmo_pre#2#1 
		for (j_1=0; j_1<=3; j_1 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_2) 
			#pragma loop name setpcmo_pre#2#1#0 
			for (i_2=0; i_2<=(3+(-1*j_1)); i_2 ++ )
			{
				r_init((double * )mtemp, 5*5, 0.0);
				r_init((double * )p, (5*5)*5, 0.0);
				mtemp[1+j_1][(1+i_2)+j_1]=1.0;
				/* when i, j=LX1-1, mortar points are duplicated, so mtemp needs */
				/* to be doubled. */
				if (((1+i_2)+j_1)==(5-1))
				{
					mtemp[1+j_1][(1+i_2)+j_1]=(mtemp[1+j_1][(1+i_2)+j_1]*2.0);
				}
				if ((1+j_1)==(5-1))
				{
					mtemp[1+j_1][(1+i_2)+j_1]=(mtemp[1+j_1][(1+i_2)+j_1]*2.0);
				}
				transf_nc(mtemp, (double (* )[5])p);
				laplacian(temp, p, element_size);
				transfb_nc1(temp1, (double (* )[5])temp);
				/* values at points (i,j) and (j,i) are the same */
				pcmor_nc1[element_size][0][0][1+j_1][(1+i_2)+j_1]=temp1[1+j_1][(1+i_2)+j_1];
				pcmor_nc1[element_size][0][0][(1+i_2)+j_1][1+j_1]=temp1[1+j_1][(1+i_2)+j_1];
			}
			i=((1+i_2)+j_1);
			/* when i, j=LX1-1, mortar points are duplicated. so pcmor_nc1 needs */
			/* to be doubled on those points */
			pcmor_nc1[element_size][0][0][1+j_1][5-1]=(pcmor_nc1[element_size][0][0][1+j_1][5-1]*2.0);
			pcmor_nc1[element_size][0][0][5-1][1+j_1]=pcmor_nc1[element_size][0][0][1+j_1][5-1];
		}
		j=(1+j_1);
		pcmor_nc1[element_size][0][0][5-1][5-1]=(pcmor_nc1[element_size][0][0][5-1][5-1]*2.0);
		/* nonconforming edges */
		j=0;
		/* Normalized Loop */
		#pragma cetus private(ii, jj) 
		#pragma cetus lastprivate(i_3) 
		#pragma loop name setpcmo_pre#2#2 
		for (i_3=0; i_3<=3; i_3 ++ )
		{
			r_init((double * )mtemp, 5*5, 0.0);
			r_init((double * )p, (5*5)*5, 0.0);
			r_init((double * )p0, (5*5)*5, 0.0);
			mtemp[j][1+i_3]=1.0;
			if ((1+i_3)==(5-1))
			{
				mtemp[j][1+i_3]=2.0;
			}
			transf_nc(mtemp, (double (* )[5])p);
			laplacian(temp, p, element_size);
			transfb_nc1(temp1, (double (* )[5])temp);
			pcmor_nc1[element_size][0][0][j][1+i_3]=temp1[j][1+i_3];
			pcmor_nc1[element_size][0][0][1+i_3][j]=temp1[j][1+i_3];
			#pragma cetus private(ii, jj) 
			#pragma loop name setpcmo_pre#2#2#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(ii, jj)
			*/
			for (ii=0; ii<5; ii ++ )
			{
				/* p0 is for the case that a nonconforming edge is shared by */
				/* two conforming faces */
				p0[0][0][ii]=p[0][0][ii];
				#pragma cetus private(jj) 
				#pragma loop name setpcmo_pre#2#2#0#0 
				for (jj=0; jj<5; jj ++ )
				{
					/* now p is for the case that a nonconforming edge is shared */
					/* by nonconforming faces */
					p[jj][0][ii]=p[0][jj][ii];
				}
			}
			laplacian(temp, p, element_size);
			transfb_nc2(temp1, (double (* )[5])temp);
			/* pcmor_nc2[element_size][jj][ii][j][i] gives the intermediate */
			/* preconditioner value on collocation point (i,j) on a  */
			/* nonconforming face of an element with size size_element */
			pcmor_nc2[element_size][0][0][j][1+i_3]=(temp1[j][1+i_3]*2.0);
			pcmor_nc2[element_size][0][0][1+i_3][j]=pcmor_nc2[element_size][0][0][j][1+i_3];
			laplacian(temp, p0, element_size);
			transfb_nc0(temp1, temp);
			/* pcmor_nc0[element_size][jj][ii][j][i] gives the intermediate */
			/* preconditioner value on collocation point (i,j) on a  */
			/* conforming face of an element, which shares a nonconforming  */
			/* edge with another conforming face */
			pcmor_nc0[element_size][0][0][j][1+i_3]=temp1[j][1+i_3];
			pcmor_nc0[element_size][0][0][1+i_3][j]=temp1[j][1+i_3];
		}
		i=(1+i_3);
		pcmor_nc1[element_size][0][0][j][5-1]=(pcmor_nc1[element_size][0][0][j][5-1]*2.0);
		pcmor_nc1[element_size][0][0][5-1][j]=pcmor_nc1[element_size][0][0][j][5-1];
		pcmor_nc2[element_size][0][0][j][5-1]=(pcmor_nc2[element_size][0][0][j][5-1]*2.0);
		pcmor_nc2[element_size][0][0][5-1][j]=pcmor_nc2[element_size][0][0][j][5-1];
		pcmor_nc0[element_size][0][0][j][5-1]=(pcmor_nc0[element_size][0][0][j][5-1]*2.0);
		pcmor_nc0[element_size][0][0][5-1][j]=pcmor_nc0[element_size][0][0][j][5-1];
		/* symmetrical copy */
		#pragma cetus private(i) 
		#pragma loop name setpcmo_pre#2#3 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i)
		*/
		for (i=0; i<(5-1); i ++ )
		{
			pcmor_nc1[element_size][1][0][j][i]=pcmor_nc1[element_size][0][0][j][(5-1)-i];
			pcmor_nc0[element_size][1][0][j][i]=pcmor_nc0[element_size][0][0][j][(5-1)-i];
			pcmor_nc2[element_size][1][0][j][i]=pcmor_nc2[element_size][0][0][j][(5-1)-i];
		}
		/* Normalized Loop */
		#pragma cetus private(i) 
		#pragma cetus lastprivate(j_2) 
		#pragma loop name setpcmo_pre#2#4 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i) lastprivate(j_2)
		*/
		for (j_2=0; j_2<=3; j_2 ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#2#4#0 
			for (i=0; i<(5-1); i ++ )
			{
				pcmor_nc1[element_size][1][0][1+j_2][i]=pcmor_nc1[element_size][0][0][1+j_2][(5-1)-i];
			}
			i=(5-1);
			pcmor_nc1[element_size][1][0][1+j_2][i]=pcmor_nc1[element_size][0][0][1+j_2][(5-1)-i];
			pcmor_nc0[element_size][1][0][1+j_2][i]=pcmor_nc0[element_size][0][0][1+j_2][(5-1)-i];
			pcmor_nc2[element_size][1][0][1+j_2][i]=pcmor_nc2[element_size][0][0][1+j_2][(5-1)-i];
		}
		j=(1+j_2);
		j=0;
		i=0;
		pcmor_nc1[element_size][0][1][j][i]=pcmor_nc1[element_size][0][0][(5-1)-j][i];
		pcmor_nc0[element_size][0][1][j][i]=pcmor_nc0[element_size][0][0][(5-1)-j][i];
		pcmor_nc2[element_size][0][1][j][i]=pcmor_nc2[element_size][0][0][(5-1)-j][i];
		/* Normalized Loop */
		#pragma cetus private(i, i_4) 
		#pragma cetus lastprivate(j_3) 
		#pragma loop name setpcmo_pre#2#5 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, i_4) lastprivate(j_3)
		*/
		for (j_3=0; j_3<=2; j_3 ++ )
		{
			i=0;
			pcmor_nc1[element_size][0][1][1+j_3][i]=pcmor_nc1[element_size][0][0][(5-1)-(1+j_3)][i];
			pcmor_nc0[element_size][0][1][1+j_3][i]=pcmor_nc0[element_size][0][0][(5-1)-(1+j_3)][i];
			pcmor_nc2[element_size][0][1][1+j_3][i]=pcmor_nc2[element_size][0][0][(5-1)-(1+j_3)][i];
			/* Normalized Loop */
			#pragma cetus lastprivate(i_4) 
			#pragma loop name setpcmo_pre#2#5#0 
			for (i_4=0; i_4<=3; i_4 ++ )
			{
				pcmor_nc1[element_size][0][1][1+j_3][1+i_4]=pcmor_nc1[element_size][0][0][(5-1)-(1+j_3)][1+i_4];
			}
			i=(1+i_4);
		}
		j=(1+j_3);
		j=(5-1);
		/* Normalized Loop */
		#pragma cetus lastprivate(i_5) 
		#pragma loop name setpcmo_pre#2#6 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_5)
		*/
		for (i_5=0; i_5<=3; i_5 ++ )
		{
			pcmor_nc1[element_size][0][1][j][1+i_5]=pcmor_nc1[element_size][0][0][(5-1)-j][1+i_5];
			pcmor_nc0[element_size][0][1][j][1+i_5]=pcmor_nc0[element_size][0][0][(5-1)-j][1+i_5];
			pcmor_nc2[element_size][0][1][j][1+i_5]=pcmor_nc2[element_size][0][0][(5-1)-j][1+i_5];
		}
		i=(1+i_5);
		j=0;
		i=(5-1);
		pcmor_nc1[element_size][1][1][j][i]=pcmor_nc1[element_size][0][0][(5-1)-j][(5-1)-i];
		pcmor_nc0[element_size][1][1][j][i]=pcmor_nc0[element_size][0][0][(5-1)-j][(5-1)-i];
		pcmor_nc2[element_size][1][1][j][i]=pcmor_nc2[element_size][0][0][(5-1)-j][(5-1)-i];
		/* Normalized Loop */
		#pragma cetus private(i, i_6) 
		#pragma cetus lastprivate(j_4) 
		#pragma loop name setpcmo_pre#2#7 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, i_6) lastprivate(j_4)
		*/
		for (j_4=0; j_4<=2; j_4 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_6) 
			#pragma loop name setpcmo_pre#2#7#0 
			for (i_6=0; i_6<=2; i_6 ++ )
			{
				pcmor_nc1[element_size][1][1][1+j_4][1+i_6]=pcmor_nc1[element_size][0][0][(5-1)-(1+j_4)][(5-1)-(1+i_6)];
			}
			i=(1+i_6);
			i=(5-1);
			pcmor_nc1[element_size][1][1][1+j_4][i]=pcmor_nc1[element_size][0][0][(5-1)-(1+j_4)][(5-1)-i];
			pcmor_nc0[element_size][1][1][1+j_4][i]=pcmor_nc0[element_size][0][0][(5-1)-(1+j_4)][(5-1)-i];
			pcmor_nc2[element_size][1][1][1+j_4][i]=pcmor_nc2[element_size][0][0][(5-1)-(1+j_4)][(5-1)-i];
		}
		j=(1+j_4);
		j=(5-1);
		/* Normalized Loop */
		#pragma cetus lastprivate(i_7) 
		#pragma loop name setpcmo_pre#2#8 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_7)
		*/
		for (i_7=0; i_7<=2; i_7 ++ )
		{
			pcmor_nc1[element_size][1][1][j][1+i_7]=pcmor_nc1[element_size][0][0][(5-1)-j][(5-1)-(1+i_7)];
			pcmor_nc0[element_size][1][1][j][1+i_7]=pcmor_nc0[element_size][0][0][(5-1)-j][(5-1)-(1+i_7)];
			pcmor_nc2[element_size][1][1][j][1+i_7]=pcmor_nc2[element_size][0][0][(5-1)-j][(5-1)-(1+i_7)];
		}
		i=(1+i_7);
		/* vertices shared by at least one nonconforming face or edge */
		/* Among three edges and three faces sharing a vertex on an element */
		/* situation 1: only one edge is nonconforming */
		/* situation 2: two edges are nonconforming */
		/* situation 3: three edges are nonconforming */
		/* situation 4: one face is nonconforming  */
		/* situation 5: one face and one edge are nonconforming  */
		/* situation 6: two faces are nonconforming */
		/* situation 7: three faces are nonconforming */
		r_init((double * )p0, (5*5)*5, 0.0);
		p0[0][0][0]=1.0;
		laplacian(temp, p0, element_size);
		pcmor_cor[element_size][7]=temp[0][0][0];
		/* situation 1 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i) 
		#pragma loop name setpcmo_pre#2#9 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i)
		*/
		for (i=0; i<5; i ++ )
		{
			p0[0][0][i]=tcpre[0][i];
		}
		laplacian(temp, p0, element_size);
		transfb_cor_e(1,  & pcmor_cor[element_size][0], temp);
		/* situation 2 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i) 
		#pragma loop name setpcmo_pre#2#10 
		for (i=0; i<5; i ++ )
		{
			p0[0][0][i]=tcpre[0][i];
			p0[0][i][0]=tcpre[0][i];
		}
		laplacian(temp, p0, element_size);
		transfb_cor_e(2,  & pcmor_cor[element_size][1], temp);
		/* situation 3 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i) 
		#pragma loop name setpcmo_pre#2#11 
		for (i=0; i<5; i ++ )
		{
			p0[0][0][i]=tcpre[0][i];
			p0[0][i][0]=tcpre[0][i];
			p0[i][0][0]=tcpre[0][i];
		}
		laplacian(temp, p0, element_size);
		transfb_cor_e(3,  & pcmor_cor[element_size][2], temp);
		/* situation 4 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i, j) 
		#pragma loop name setpcmo_pre#2#12 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, j)
		*/
		for (j=0; j<5; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#2#12#0 
			for (i=0; i<5; i ++ )
			{
				p0[0][j][i]=tcpre[j][i];
			}
		}
		laplacian(temp, p0, element_size);
		transfb_cor_f(4,  & pcmor_cor[element_size][3], temp);
		/* situation 5 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i, j) 
		#pragma loop name setpcmo_pre#2#13 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, j)
		*/
		for (j=0; j<5; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#2#13#0 
			for (i=0; i<5; i ++ )
			{
				p0[0][j][i]=tcpre[j][i];
			}
		}
		#pragma cetus private(i) 
		#pragma loop name setpcmo_pre#2#14 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i)
		*/
		for (i=0; i<5; i ++ )
		{
			p0[i][0][0]=tcpre[0][i];
		}
		laplacian(temp, p0, element_size);
		transfb_cor_f(5,  & pcmor_cor[element_size][4], temp);
		/* situation 6 */
		r_init((double * )p0, (5*5)*5, 0.0);
		#pragma cetus private(i, j) 
		#pragma loop name setpcmo_pre#2#15 
		for (j=0; j<5; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#2#15#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(i)
			*/
			for (i=0; i<5; i ++ )
			{
				p0[0][j][i]=tcpre[j][i];
				p0[j][0][i]=tcpre[j][i];
			}
		}
		laplacian(temp, p0, element_size);
		transfb_cor_f(6,  & pcmor_cor[element_size][5], temp);
		/* situation 7 */
		#pragma cetus private(i, j) 
		#pragma loop name setpcmo_pre#2#16 
		for (j=0; j<5; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name setpcmo_pre#2#16#0 
			for (i=0; i<5; i ++ )
			{
				p0[0][j][i]=tcpre[j][i];
				p0[j][0][i]=tcpre[j][i];
				p0[j][i][0]=tcpre[j][i];
			}
		}
		laplacian(temp, p0, element_size);
		transfb_cor_f(7,  & pcmor_cor[element_size][6], temp);
	}
}

/* ------------------------------------------------------------------------ */
/* compute the preconditioner by identifying its geometry configuration */
/* and sum the values from the precomputed elemental contributions */
/* ------------------------------------------------------------------------ */
void setpcmo()
{
	int face2, nb1, nb2, sizei, imor, _enum, i, j, iel, iside, nn1, nn2;
	int j_0;
	int i_0;
	int j_1;
	int i_1;
	l_init(ifpcmor, nvertex, false);
	l_init((logical * )edgevis, 24*nelt, false);
	#pragma cetus private(_enum, face2, i, i_0, i_1, iel, imor, iside, j, j_0, j_1, nb1, nb2, nn1, nn2, sizei) 
	#pragma loop name setpcmo#0 
	for (iel=0; iel<nelt; iel ++ )
	{
		#pragma cetus private(_enum, face2, i, i_0, i_1, imor, iside, j, j_0, j_1, nb1, nb2, nn1, nn2, sizei) 
		#pragma loop name setpcmo#0#0 
		for (iside=0; iside<6; iside ++ )
		{
			/* for nonconforming faces */
			if (cbc[iel][iside]==3)
			{
				sizei=size_e[iel];
				/* vertices */
				/* ifpcmor[imor] = true indicates that mortar point imor has  */
				/* been visited */
				imor=idmo[iel][iside][0][0][0][0];
				if ( ! ifpcmor[imor])
				{
					/* compute the preconditioner on mortar point imor */
					pc_corner(imor);
					ifpcmor[imor]=true;
				}
				imor=idmo[iel][iside][1][0][0][5-1];
				if ( ! ifpcmor[imor])
				{
					pc_corner(imor);
					ifpcmor[imor]=true;
				}
				imor=idmo[iel][iside][0][1][5-1][0];
				if ( ! ifpcmor[imor])
				{
					pc_corner(imor);
					ifpcmor[imor]=true;
				}
				imor=idmo[iel][iside][1][1][5-1][5-1];
				if ( ! ifpcmor[imor])
				{
					pc_corner(imor);
					ifpcmor[imor]=true;
				}
				/* edges on nonconforming faces, _enum is local edge number */
				#pragma cetus private(_enum, face2, nb1, nb2) 
				#pragma loop name setpcmo#0#0#0 
				for (_enum=0; _enum<4; _enum ++ )
				{
					/* edgevis[iel][iside][_enum]=true indicates that local edge  */
					/* _enum of face iside of iel has been visited */
					if ( ! edgevis[iel][iside][_enum])
					{
						edgevis[iel][iside][_enum]=true;
						/* Examing neighbor element information, */
						/* calculateing the preconditioner value. */
						face2=f_e_ef[iside][_enum];
						if (cbc[iel][face2]==2)
						{
							nb1=sje[iel][face2][0][0];
							if (cbc[nb1][iside]==2)
							{
								/* Compute the preconditioner on local edge _enum on face */
								/* iside of element iel, 1 is neighborhood information got */
								/* by examing neighbors(nb1). For detailed meaning of 1,  */
								/* see subroutine com_dpc. */
								com_dpc(iside, iel, _enum, 1, sizei);
								nb2=sje[nb1][iside][0][0];
								edgevis[nb2][jjface[face2]][op[e_face2[iside][_enum]]]=true;
							}
							else
							{
								if (cbc[nb1][iside]==3)
								{
									com_dpc(iside, iel, _enum, 2, sizei);
									edgevis[nb1][iside][op[_enum]]=true;
								}
							}
						}
						else
						{
							if (cbc[iel][face2]==3)
							{
								edgevis[iel][face2][e_face2[iside][_enum]]=true;
								nb1=sje[iel][face2][1][0];
								if (cbc[nb1][iside]==1)
								{
									com_dpc(iside, iel, _enum, 3, sizei);
									nb2=sje[nb1][iside][0][0];
									edgevis[nb2][jjface[iside]][op[_enum]]=true;
									edgevis[nb2][jjface[face2]][op[e_face2[iside][_enum]]]=true;
								}
								else
								{
									if (cbc[nb1][iside]==2)
									{
										com_dpc(iside, iel, _enum, 4, sizei);
									}
								}
							}
							else
							{
								if (cbc[iel][face2]==0)
								{
									com_dpc(iside, iel, _enum, 0, sizei);
								}
							}
						}
					}
				}
				/* mortar element interior (not edge of mortar)  */
				#pragma cetus private(i, i_0, imor, j, j_0, nn1, nn2) 
				#pragma loop name setpcmo#0#0#1 
				for (nn1=0; nn1<2; nn1 ++ )
				{
					#pragma cetus private(i, i_0, imor, j, j_0, nn2) 
					#pragma loop name setpcmo#0#0#1#0 
					for (nn2=0; nn2<2; nn2 ++ )
					{
						/* Normalized Loop */
						#pragma cetus private(i, i_0, imor) 
						#pragma cetus lastprivate(j_0) 
						#pragma loop name setpcmo#0#0#1#0#0 
						for (j_0=0; j_0<=2; j_0 ++ )
						{
							/* Normalized Loop */
							#pragma cetus private(imor) 
							#pragma cetus lastprivate(i_0) 
							#pragma loop name setpcmo#0#0#1#0#0#0 
							for (i_0=0; i_0<=2; i_0 ++ )
							{
								imor=idmo[iel][iside][nn2][nn1][1+j_0][1+i_0];
								dpcmor[imor]=(1.0/(pcmor_nc1[sizei][nn2][nn1][1+j_0][1+i_0]+pcmor_c[sizei+1][1+j_0][1+i_0]));
							}
							i=(1+i_0);
						}
						j=(1+j_0);
					}
				}
				/* for i,j=LX1-1 there are duplicated mortar points, so  */
				/* pcmor_c needs to be doubled or quadrupled */
				i=(5-1);
				/* Normalized Loop */
				#pragma cetus private(imor) 
				#pragma cetus lastprivate(j_1) 
				#pragma loop name setpcmo#0#0#2 
				for (j_1=0; j_1<=2; j_1 ++ )
				{
					imor=idmo[iel][iside][0][0][1+j_1][i];
					dpcmor[imor]=(1.0/(pcmor_nc1[sizei][0][0][1+j_1][i]+(pcmor_c[sizei+1][1+j_1][i]*2.0)));
					imor=idmo[iel][iside][0][1][1+j_1][i];
					dpcmor[imor]=(1.0/(pcmor_nc1[sizei][0][1][1+j_1][i]+(pcmor_c[sizei+1][1+j_1][i]*2.0)));
				}
				j=(1+j_1);
				j=(5-1);
				imor=idmo[iel][iside][0][0][j][i];
				dpcmor[imor]=(1.0/(pcmor_nc1[sizei][0][0][j][i]+(pcmor_c[sizei+1][j][i]*4.0)));
				/* Normalized Loop */
				#pragma cetus private(imor) 
				#pragma cetus lastprivate(i_1) 
				#pragma loop name setpcmo#0#0#3 
				for (i_1=0; i_1<=2; i_1 ++ )
				{
					imor=idmo[iel][iside][0][0][j][1+i_1];
					dpcmor[imor]=(1.0/(pcmor_nc1[sizei][0][0][j][1+i_1]+(pcmor_c[sizei+1][j][1+i_1]*2.0)));
					imor=idmo[iel][iside][1][0][j][1+i_1];
					dpcmor[imor]=(1.0/(pcmor_nc1[sizei][1][0][j][1+i_1]+(pcmor_c[sizei+1][j][1+i_1]*2.0)));
				}
				i=(1+i_1);
			}
		}
	}
}

/* ------------------------------------------------------------------------ */
/* calculate preconditioner value for vertex with mortar index imor */
/* ------------------------------------------------------------------------ */
static void pc_corner(int imor)
{
	double tmortemp;
	int inemo, ie, sizei, cornernumber;
	int sface, sedge, iiface, iface, iiedge, iedge, n = 0;
	tmortemp=0.0;
	/* loop over all elements sharing this vertex */
	#pragma cetus private(cornernumber, ie, iedge, iface, iiedge, iiface, inemo, sedge, sface, sizei) 
	#pragma loop name pc_corner#0 
	/* #pragma cetus reduction(+: tmortemp)  */
	for (inemo=0; inemo<=nemo[imor]; inemo ++ )
	{
		ie=emo[imor][inemo][0];
		sizei=size_e[ie];
		cornernumber=emo[imor][inemo][1];
		sface=0;
		sedge=0;
		#pragma cetus private(iface, iiface) 
		#pragma loop name pc_corner#0#0 
		#pragma cetus reduction(+: sface) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(iface, iiface) reduction(+: sface)
		*/
		for (iiface=0; iiface<3; iiface ++ )
		{
			iface=f_c[cornernumber][iiface];
			/* sface sums the number of nonconforming faces sharing this vertex on */
			/* one element */
			if (cbc[ie][iface]==3)
			{
				sface=(sface+1);
			}
		}
		/* sedge sums the number of nonconforming edges sharing this vertex on */
		/* one element */
		#pragma cetus private(iedge, iiedge) 
		#pragma loop name pc_corner#0#1 
		#pragma cetus reduction(+: sedge) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(iedge, iiedge) reduction(+: sedge)
		*/
		for (iiedge=0; iiedge<3; iiedge ++ )
		{
			iedge=e_c[cornernumber][iiedge];
			if (ncon_edge[ie][iedge])
			{
				sedge=(sedge+1);
			}
		}
		/* each n indicates how many nonconforming faces and nonconforming */
		/* edges share this vertex on an element,  */
		if (sface==0)
		{
			if (sedge==0)
			{
				n=7;
			}
			else
			{
				if (sedge==1)
				{
					n=0;
				}
				else
				{
					if (sedge==2)
					{
						n=1;
					}
					else
					{
						if (sedge==3)
						{
							n=2;
						}
					}
				}
			}
		}
		else
		{
			if (sface==1)
			{
				if (sedge==1)
				{
					n=4;
				}
				else
				{
					n=3;
				}
			}
			else
			{
				if (sface==2)
				{
					n=5;
				}
				else
				{
					if (sface==3)
					{
						n=6;
					}
				}
			}
		}
		/* sum the intermediate pre-computed preconditioner values for  */
		/* all elements */
		tmortemp=(tmortemp+pcmor_cor[sizei][n]);
	}
	/* dpcmor[imor] is the value of the preconditioner on mortar point imor */
	dpcmor[imor]=(1.0/tmortemp);
}

/* ------------------------------------------------------------------------ */
/* Compute preconditioner for local edge enumber of face iside  */
/* on element iel. */
/* isize is element size, */
/* n is one of five different configurations */
/* anc1, ac, anc2, anc0 are coefficients for different edges.  */
/* nc0 refers to nonconforming edge shared by two conforming faces */
/* nc1 refers to nonconforming edge shared by one nonconforming face */
/* nc2 refers to nonconforming edges shared by two nonconforming faces */
/* c refers to conforming edge */
/* ------------------------------------------------------------------------ */
static void com_dpc(int iside, int iel, int enumber, int n, int isize)
{
	int nn1start, nn1end, nn2start;
	int nn2end, jstart, jend, istart, iend, i, j, nn1, nn2, imor = 0;
	double anc1, ac, anc2, anc0, temp = 0.0;
	/* different local edges have different loop ranges  */
	int nn2_0;
	int nn1_0;
	int j_0;
	int i_0;
	if (enumber==0)
	{
		nn1start=1;
		nn1end=1;
		nn2start=1;
		nn2end=2;
		jstart=1;
		jend=1;
		istart=2;
		iend=(5-1);
	}
	else
	{
		if (enumber==1)
		{
			nn1start=1;
			nn1end=2;
			nn2start=2;
			nn2end=2;
			jstart=2;
			jend=(5-1);
			istart=5;
			iend=5;
		}
		else
		{
			if (enumber==2)
			{
				nn1start=2;
				nn1end=2;
				nn2start=1;
				nn2end=2;
				jstart=5;
				jend=5;
				istart=2;
				iend=(5-1);
			}
			else
			{
				if (enumber==3)
				{
					nn1start=1;
					nn1end=2;
					nn2start=1;
					nn2end=1;
					jstart=2;
					jend=(5-1);
					istart=1;
					iend=1;
				}
				else
				{
					/* MUST NOT reachable!! */
					(((void)sizeof (0 ? 1 : 0)), ({
						if (0)
						{
							;
						}
						else
						{
							__assert_fail("0", "precond.c", 685, __PRETTY_FUNCTION__);
						}
					}));
					nn1start=0;
					nn1end=0;
					nn2start=0;
					nn2end=0;
					jstart=0;
					jend=(5-1);
					istart=0;
					iend=0;
				}
			}
		}
	}
	/* among the four elements sharing this edge */
	/* one has a smaller size */
	if (n==1)
	{
		anc1=2.0;
		ac=1.0;
		anc0=1.0;
		anc2=0.0;
		/* two (neighbored by a face) are of  smaller size */
	}
	else
	{
		if (n==2)
		{
			anc1=2.0;
			ac=2.0;
			anc0=0.0;
			anc2=0.0;
			/* two (neighbored by an edge) are of smaller size */
		}
		else
		{
			if (n==3)
			{
				anc2=2.0;
				ac=2.0;
				anc1=0.0;
				anc0=0.0;
				/* three are of smaller size */
			}
			else
			{
				if (n==4)
				{
					anc1=0.0;
					ac=3.0;
					anc2=1.0;
					anc0=0.0;
					/* on the boundary */
				}
				else
				{
					if (n==0)
					{
						anc1=1.0;
						ac=1.0;
						anc2=0.0;
						anc0=0.0;
					}
					else
					{
						/* MUST NOT reachable!! */
						(((void)sizeof (0 ? 1 : 0)), ({
							if (0)
							{
								;
							}
							else
							{
								__assert_fail("0", "precond.c", 734, __PRETTY_FUNCTION__);
							}
						}));
						anc1=0.0;
						ac=0.0;
						anc2=0.0;
						anc0=0.0;
					}
				}
			}
		}
	}
	/* edge interior */
	/* Normalized Loop */
	#pragma cetus private(i, i_0, j, j_0, nn1, nn1_0) 
	#pragma cetus lastprivate(imor, nn2_0, temp) 
	#pragma loop name com_dpc#0 
	for (nn2_0=0; nn2_0<=(nn2end+(-1*nn2start)); nn2_0 ++ )
	{
		/* Normalized Loop */
		#pragma cetus private(i, i_0, j, j_0) 
		#pragma cetus lastprivate(imor, nn1_0, temp) 
		#pragma loop name com_dpc#0#0 
		for (nn1_0=0; nn1_0<=(nn1end+(-1*nn1start)); nn1_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(i, i_0) 
			#pragma cetus lastprivate(imor, j_0, temp) 
			#pragma loop name com_dpc#0#0#0 
			for (j_0=0; j_0<=(jend+(-1*jstart)); j_0 ++ )
			{
				/* Normalized Loop */
				#pragma cetus lastprivate(i_0, imor, temp) 
				#pragma loop name com_dpc#0#0#0#0 
				for (i_0=0; i_0<=(iend+(-1*istart)); i_0 ++ )
				{
					imor=idmo[iel][iside][(-1+nn2_0)+nn2start][(-1+nn1_0)+nn1start][(-1+j_0)+jstart][(-1+i_0)+istart];
					temp=((((anc1*pcmor_nc1[isize][(-1+nn2_0)+nn2start][(-1+nn1_0)+nn1start][(-1+j_0)+jstart][(-1+i_0)+istart])+(ac*pcmor_c[isize+1][(-1+j_0)+jstart][(-1+i_0)+istart]))+(anc0*pcmor_nc0[isize][(-1+nn2_0)+nn2start][(-1+nn1_0)+nn1start][(-1+j_0)+jstart][(-1+i_0)+istart]))+(anc2*pcmor_nc2[isize][(-1+nn2_0)+nn2start][(-1+nn1_0)+nn1start][(-1+j_0)+jstart][(-1+i_0)+istart]));
					dpcmor[imor]=(1.0/temp);
				}
				i=((-1+i_0)+istart);
			}
			j=((-1+j_0)+jstart);
		}
		nn1=((-1+nn1_0)+nn1start);
	}
	nn2=((-1+nn2_0)+nn2start);
	/* local edge 0 */
	if (enumber==0)
	{
		imor=idmo[iel][iside][0][0][0][5-1];
		temp=((((anc1*pcmor_nc1[isize][0][0][0][5-1])+((ac*pcmor_c[isize+1][0][5-1])*2.0))+(anc0*pcmor_nc0[isize][0][0][0][5-1]))+(anc2*pcmor_nc2[isize][0][0][0][5-1]));
		/* local edge 1 */
	}
	else
	{
		if (enumber==1)
		{
			imor=idmo[iel][iside][1][0][5-1][5-1];
			temp=((((anc1*pcmor_nc1[isize][1][0][5-1][5-1])+((ac*pcmor_c[isize+1][5-1][5-1])*2.0))+(anc0*pcmor_nc0[isize][1][0][5-1][5-1]))+(anc2*pcmor_nc2[isize][1][0][5-1][5-1]));
			/* local edge 2 */
		}
		else
		{
			if (enumber==2)
			{
				imor=idmo[iel][iside][0][1][5-1][5-1];
				temp=((((anc1*pcmor_nc1[isize][0][1][5-1][5-1])+((ac*pcmor_c[isize+1][5-1][5-1])*2.0))+(anc0*pcmor_nc0[isize][0][1][5-1][5-1]))+(anc2*pcmor_nc2[isize][0][1][5-1][5-1]));
				/* local edge 3 */
			}
			else
			{
				if (enumber==3)
				{
					imor=idmo[iel][iside][0][0][5-1][0];
					temp=((((anc1*pcmor_nc1[isize][0][0][5-1][0])+((ac*pcmor_c[isize+1][5-1][0])*2.0))+(anc0*pcmor_nc0[isize][0][0][5-1][0]))+(anc2*pcmor_nc2[isize][0][0][5-1][0]));
				}
			}
		}
	}
	dpcmor[imor]=(1.0/temp);
}
