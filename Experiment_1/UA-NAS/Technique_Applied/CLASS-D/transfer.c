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
#include "header.h"
#include <stdio.h>
int count_init_transf = 0, count_init_transfb = 0;
/* ------------------------------------------------------------------ */
/* Map values from mortar(tmor) to element(tx) */
/* ------------------------------------------------------------------ */
void transf(double tmor[], double tx[])
{
	double tmp[2][5][5];
	int ig1, ig2, ig3, ig4, ie, iface, il1, il2, il3, il4;
	int nnje, ije1, ije2, col, i, j, ig, il;
	int iel, ntemp;
	/* zero out tx on element boundaries */
	int x;
	int y;
	int z;
	int i_0;
	int i_1;
	int col_0;
	int i_2;
	int i_3;
	int col_1;
	int i_4;
	int i_5;
	int i_6;
	int i_7;
	int i_8;
	int i_9;
	int i_10;
	int i_11;
	int i_12;
	col2(tx, (double * )tmult, ntot);
	v_end[0]=0;
	v_end[1]=(5-1);
	if (count_init_transf==0)
	{
		#pragma cetus private(i, iel, j, ntemp) 
		#pragma loop name transf#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(i, iel, j, ntemp)
		for (iel=0; iel<515000; iel ++ )
		{
			ntemp=(((5*5)*5)*iel);
			#pragma cetus private(i, j) 
			#pragma loop name transf#0#0 
			for (j=0; j<5; j ++ )
			{
				#pragma cetus private(i) 
				#pragma loop name transf#0#0#0 
				for (i=0; i<5; i ++ )
				{
					idel[iel][0][j][i]=((((ntemp+(i*5))+((j*5)*5))+5)-1);
					idel[iel][1][j][i]=((ntemp+(i*5))+((j*5)*5));
					idel[iel][2][j][i]=(((ntemp+(i*1))+((j*5)*5))+(5*(5-1)));
					idel[iel][3][j][i]=((ntemp+(i*1))+((j*5)*5));
					idel[iel][4][j][i]=(((ntemp+(i*1))+(j*5))+((5*5)*(5-1)));
					idel[iel][5][j][i]=((ntemp+(i*1))+(j*5));
				}
			}
		}
		count_init_transf ++ ;
	}
	#pragma cetus private(col, col_0, col_1, i, i_0, i_1, i_10, i_11, i_12, i_2, i_3, i_4, i_5, i_6, i_7, i_8, i_9, ie, iface, ig, ig1, ig2, ig3, ig4, ije1, ije2, il, il1, il2, il3, il4, j, nnje, tmp, x, y, z) 
	#pragma loop name transf#1 
	#pragma cetus parallel 
	#pragma omp parallel for if(((-1+nelt)<=514999)) private(col, col_0, col_1, i, i_0, i_1, i_10, i_11, i_12, i_2, i_3, i_4, i_5, i_6, i_7, i_8, i_9, ie, iface, ig, ig1, ig2, ig3, ig4, ije1, ije2, il, il1, il2, il3, il4, j, nnje, tmp, x, y, z)
	for (ie=0; ie<nelt; ie ++ )
	{
		#pragma cetus private(col, col_0, col_1, i, i_0, i_1, i_10, i_11, i_12, i_2, i_3, i_4, i_5, i_6, i_7, i_8, i_9, iface, ig, ig1, ig2, ig3, ig4, ije1, ije2, il, il1, il2, il3, il4, j, nnje, x, y, z) 
		#pragma cetus lastprivate(tmp) 
		#pragma loop name transf#1#0 
		for (iface=0; iface<6; iface ++ )
		{
			/* get the collocation point index of the four local corners on the */
			/* face iface of element ie */
			il1=idel[ie][iface][0][0];
			il2=idel[ie][iface][0][5-1];
			il3=idel[ie][iface][5-1][0];
			il4=idel[ie][iface][5-1][5-1];
			/* get the mortar indices of the four local corners */
			ig1=idmo[ie][iface][0][0][0][0];
			ig2=idmo[ie][iface][1][0][0][5-1];
			ig3=idmo[ie][iface][0][1][5-1][0];
			ig4=idmo[ie][iface][1][1][5-1][5-1];
			/* copy the value from tmor to tx for these four local corners */
			tx[il1]=tmor[ig1];
			tx[il2]=tmor[ig2];
			tx[il3]=tmor[ig3];
			tx[il4]=tmor[ig4];
			/* nnje=1 for conforming faces, nnje=2 for nonconforming faces */
			if (cbc[ie][iface]==3)
			{
				nnje=2;
			}
			else
			{
				nnje=1;
			}
			/* for nonconforming faces */
			if (nnje==2)
			{
				/* nonconforming faces have four pieces of mortar, first map them to */
				/* two intermediate mortars, stored in tmp */
				/* r_init((double)tmp, LX1*LX1*2, 0.0); */
				#pragma cetus private(x, y, z) 
				#pragma loop name transf#1#0#0 
				for (x=0; x<nnje; x ++ )
				{
					#pragma cetus private(y, z) 
					#pragma loop name transf#1#0#0#0 
					for (y=0; y<5; y ++ )
					{
						#pragma cetus private(z) 
						#pragma loop name transf#1#0#0#0#0 
						for (z=0; z<5; z ++ )
						{
							tmp[x][y][z]=0.0;
						}
					}
				}
				#pragma cetus firstprivate(tmp) 
				#pragma cetus private(col, i, i_0, ig, ije1, ije2, il, j) 
				#pragma cetus lastprivate(tmp) 
				#pragma loop name transf#1#0#1 
				for (ije1=0; ije1<nnje; ije1 ++ )
				{
					#pragma cetus firstprivate(tmp) 
					#pragma cetus private(col, i, i_0, ig, ije2, il, j) 
					#pragma cetus lastprivate(tmp) 
					#pragma loop name transf#1#0#1#0 
					for (ije2=0; ije2<nnje; ije2 ++ )
					{
						#pragma cetus private(col, i, i_0, ig, il, j) 
						#pragma loop name transf#1#0#1#0#0 
						for (col=0; col<5; col ++ )
						{
							/* in each row col, when coloumn i=1 or LX1, the value */
							/* in tmor is copied to tmp */
							i=v_end[ije2];
							ig=idmo[ie][iface][ije2][ije1][col][i];
							tmp[ije1][col][i]=tmor[ig];
							/* in each row col, value in the interior three collocation */
							/* points is computed by apply mapping matrix qbnew to tmor */
							/* Normalized Loop */
							#pragma cetus private(ig, il, j) 
							#pragma cetus lastprivate(i_0) 
							#pragma loop name transf#1#0#1#0#0#0 
							for (i_0=0; i_0<=2; i_0 ++ )
							{
								il=idel[ie][iface][col][1+i_0];
								#pragma cetus private(ig, j) 
								#pragma loop name transf#1#0#1#0#0#0#0 
								for (j=0; j<5; j ++ )
								{
									ig=idmo[ie][iface][ije2][ije1][col][j];
									tmp[ije1][col][1+i_0]=(tmp[ije1][col][1+i_0]+(qbnew[ije2][j][(1+i_0)-1]*tmor[ig]));
								}
							}
							i=(1+i_0);
						}
					}
				}
				/* mapping from two pieces of intermediate mortar tmp to element */
				/* face tx */
				#pragma cetus private(col, col_0, i, i_1, i_2, i_3, ije1, il, j) 
				#pragma loop name transf#1#0#2 
				for (ije1=0; ije1<nnje; ije1 ++ )
				{
					/* the first column, col=0, is an edge of face iface. */
					/* the value on the three interior collocation points, tx, is */
					/* computed by applying mapping matrices qbnew to tmp. */
					/* the mapping result is divided by 2, because there will be */
					/* duplicated contribution from another face sharing this edge. */
					col=0;
					/* Normalized Loop */
					#pragma cetus private(il, j) 
					#pragma cetus lastprivate(i_1) 
					#pragma loop name transf#1#0#2#0 
					for (i_1=0; i_1<=2; i_1 ++ )
					{
						il=idel[ie][iface][1+i_1][col];
						#pragma cetus private(j) 
						#pragma loop name transf#1#0#2#0#0 
						for (j=0; j<5; j ++ )
						{
							tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_1)-1]*tmp[ije1][j][col])*0.5));
						}
					}
					i=(1+i_1);
					/* for column 1 ~ lx-2 */
					/* Normalized Loop */
					#pragma cetus private(i, i_2, il, j) 
					#pragma cetus lastprivate(col_0) 
					#pragma loop name transf#1#0#2#1 
					for (col_0=0; col_0<=2; col_0 ++ )
					{
						/* when i=0 or LX1-1, the collocation points are also on an edge of */
						/* the face, so the mapping result also needs to be divided by 2 */
						i=v_end[ije1];
						il=idel[ie][iface][i][1+col_0];
						tx[il]=(tx[il]+(tmp[ije1][i][1+col_0]*0.5));
						/* compute the value at interior collocation points in */
						/* columns 1 ~ LX1-1 */
						/* Normalized Loop */
						#pragma cetus private(il, j) 
						#pragma cetus lastprivate(i_2) 
						#pragma loop name transf#1#0#2#1#0 
						for (i_2=0; i_2<=2; i_2 ++ )
						{
							il=idel[ie][iface][1+i_2][1+col_0];
							#pragma cetus private(j) 
							#pragma loop name transf#1#0#2#1#0#0 
							for (j=0; j<5; j ++ )
							{
								tx[il]=(tx[il]+(qbnew[ije1][j][(1+i_2)-1]*tmp[ije1][j][1+col_0]));
							}
						}
						i=(1+i_2);
					}
					col=(1+col_0);
					/* same as col=0 */
					col=(5-1);
					/* Normalized Loop */
					#pragma cetus private(il, j) 
					#pragma cetus lastprivate(i_3) 
					#pragma loop name transf#1#0#2#2 
					for (i_3=0; i_3<=2; i_3 ++ )
					{
						il=idel[ie][iface][1+i_3][col];
						#pragma cetus private(j) 
						#pragma loop name transf#1#0#2#2#0 
						for (j=0; j<5; j ++ )
						{
							tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_3)-1]*tmp[ije1][j][col])*0.5));
						}
					}
					i=(1+i_3);
				}
				/* for conforming faces */
			}
			else
			{
				/* face interior */
				/* Normalized Loop */
				#pragma cetus private(i, i_4, ig, il) 
				#pragma cetus lastprivate(col_1) 
				#pragma loop name transf#1#0#3 
				for (col_1=0; col_1<=2; col_1 ++ )
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(i_4) 
					#pragma loop name transf#1#0#3#0 
					for (i_4=0; i_4<=2; i_4 ++ )
					{
						il=idel[ie][iface][1+col_1][1+i_4];
						ig=idmo[ie][iface][0][0][1+col_1][1+i_4];
						tx[il]=tmor[ig];
					}
					i=(1+i_4);
				}
				col=(1+col_1);
				/* edges of conforming faces */
				/* if local edge 0 is a nonconforming edge */
				if (idmo[ie][iface][0][0][0][5-1]!=( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, ije1, il, j) 
					#pragma cetus lastprivate(i_5) 
					#pragma loop name transf#1#0#4 
					for (i_5=0; i_5<=2; i_5 ++ )
					{
						il=idel[ie][iface][0][1+i_5];
						#pragma cetus private(ig, ije1, j) 
						#pragma loop name transf#1#0#4#0 
						for (ije1=0; ije1<2; ije1 ++ )
						{
							#pragma cetus private(ig, j) 
							#pragma loop name transf#1#0#4#0#0 
							for (j=0; j<5; j ++ )
							{
								ig=idmo[ie][iface][ije1][0][0][j];
								tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_5)-1]*tmor[ig])*0.5));
							}
						}
					}
					i=(1+i_5);
					/* if local edge 0 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(i_6) 
					#pragma loop name transf#1#0#5 
					for (i_6=0; i_6<=2; i_6 ++ )
					{
						il=idel[ie][iface][0][1+i_6];
						ig=idmo[ie][iface][0][0][0][1+i_6];
						tx[il]=tmor[ig];
					}
					i=(1+i_6);
				}
				/* if local edge 1 is a nonconforming edge */
				if (idmo[ie][iface][1][0][1][5-1]!=( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, ije1, il, j) 
					#pragma cetus lastprivate(i_7) 
					#pragma loop name transf#1#0#6 
					for (i_7=0; i_7<=2; i_7 ++ )
					{
						il=idel[ie][iface][1+i_7][5-1];
						#pragma cetus private(ig, ije1, j) 
						#pragma loop name transf#1#0#6#0 
						for (ije1=0; ije1<2; ije1 ++ )
						{
							#pragma cetus private(ig, j) 
							#pragma loop name transf#1#0#6#0#0 
							for (j=0; j<5; j ++ )
							{
								ig=idmo[ie][iface][1][ije1][j][5-1];
								tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_7)-1]*tmor[ig])*0.5));
							}
						}
					}
					i=(1+i_7);
					/* if local edge 1 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(i_8) 
					#pragma loop name transf#1#0#7 
					for (i_8=0; i_8<=2; i_8 ++ )
					{
						il=idel[ie][iface][1+i_8][5-1];
						ig=idmo[ie][iface][0][0][1+i_8][5-1];
						tx[il]=tmor[ig];
					}
					i=(1+i_8);
				}
				/* if local edge 2 is a nonconforming edge */
				if (idmo[ie][iface][0][1][5-1][1]!=( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, ije1, il, j) 
					#pragma cetus lastprivate(i_9) 
					#pragma loop name transf#1#0#8 
					for (i_9=0; i_9<=2; i_9 ++ )
					{
						il=idel[ie][iface][5-1][1+i_9];
						#pragma cetus private(ig, ije1, j) 
						#pragma loop name transf#1#0#8#0 
						for (ije1=0; ije1<2; ije1 ++ )
						{
							#pragma cetus private(ig, j) 
							#pragma loop name transf#1#0#8#0#0 
							for (j=0; j<5; j ++ )
							{
								ig=idmo[ie][iface][ije1][1][5-1][j];
								tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_9)-1]*tmor[ig])*0.5));
							}
						}
					}
					i=(1+i_9);
					/* if local edge 2 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(i_10) 
					#pragma loop name transf#1#0#9 
					for (i_10=0; i_10<=2; i_10 ++ )
					{
						il=idel[ie][iface][5-1][1+i_10];
						ig=idmo[ie][iface][0][0][5-1][1+i_10];
						tx[il]=tmor[ig];
					}
					i=(1+i_10);
				}
				/* if local edge 3 is a nonconforming edge */
				if (idmo[ie][iface][0][0][5-1][0]!=( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, ije1, il, j) 
					#pragma cetus lastprivate(i_11) 
					#pragma loop name transf#1#0#10 
					for (i_11=0; i_11<=2; i_11 ++ )
					{
						il=idel[ie][iface][1+i_11][0];
						#pragma cetus private(ig, ije1, j) 
						#pragma loop name transf#1#0#10#0 
						for (ije1=0; ije1<2; ije1 ++ )
						{
							#pragma cetus private(ig, j) 
							#pragma loop name transf#1#0#10#0#0 
							for (j=0; j<5; j ++ )
							{
								ig=idmo[ie][iface][0][ije1][j][0];
								tx[il]=(tx[il]+((qbnew[ije1][j][(1+i_11)-1]*tmor[ig])*0.5));
							}
						}
					}
					i=(1+i_11);
					/* if local edge 3 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(i_12) 
					#pragma loop name transf#1#0#11 
					for (i_12=0; i_12<=2; i_12 ++ )
					{
						il=idel[ie][iface][1+i_12][0];
						ig=idmo[ie][iface][0][0][1+i_12][0];
						tx[il]=tmor[ig];
					}
					i=(1+i_12);
				}
			}
		}
	}
}

/* ------------------------------------------------------------------ */
/* Map from element(tx) to mortar(tmor). */
/* tmor sums contributions from all elements. */
/* ------------------------------------------------------------------ */
void transfb(double tmor[], double tx[])
{
	const double third = 1.0/3.0;
	int shift;
	double tmp, tmp1, temp[2][5][5], top[2][5];
	int il1, il2, il3, il4, ig1, ig2, ig3, ig4, ie, iface, nnje;
	int ije1, ije2, col, i, j, ije, ig, il;
	int iel, ntemp;
	int i_0;
	int j_0;
	int i_1;
	int col_0;
	int i_2;
	int i_3;
	int col_1;
	int j_1;
	int i_4;
	int j_2;
	int i_5;
	int j_3;
	int i_6;
	int j_4;
	int i_7;
	int j_5;
	r_init(tmor, nmor, 0.0);
	if (count_init_transfb==0)
	{
		#pragma cetus private(i, iel, j, ntemp) 
		#pragma loop name transfb#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(i, iel, j, ntemp)
		for (iel=0; iel<515000; iel ++ )
		{
			ntemp=(((5*5)*5)*iel);
			#pragma cetus private(i, j) 
			#pragma loop name transfb#0#0 
			for (j=0; j<5; j ++ )
			{
				#pragma cetus private(i) 
				#pragma loop name transfb#0#0#0 
				for (i=0; i<5; i ++ )
				{
					idel[iel][0][j][i]=((((ntemp+(i*5))+((j*5)*5))+5)-1);
					idel[iel][1][j][i]=((ntemp+(i*5))+((j*5)*5));
					idel[iel][2][j][i]=(((ntemp+(i*1))+((j*5)*5))+(5*(5-1)));
					idel[iel][3][j][i]=((ntemp+(i*1))+((j*5)*5));
					idel[iel][4][j][i]=(((ntemp+(i*1))+(j*5))+((5*5)*(5-1)));
					idel[iel][5][j][i]=((ntemp+(i*1))+(j*5));
				}
			}
		}
		count_init_transfb ++ ;
	}
	#pragma cetus private(col, col_0, col_1, i, i_0, i_1, i_2, i_3, i_4, i_5, i_6, i_7, ie, iface, ig, ig1, ig2, ig3, ig4, ije, ije1, ije2, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4, j_5, nnje, shift, tmp, tmp1) 
	#pragma loop name transfb#1 
	for (ie=0; ie<nelt; ie ++ )
	{
		#pragma cetus private(col, col_0, col_1, i, i_0, i_1, i_2, i_3, i_4, i_5, i_6, i_7, iface, ig, ig1, ig2, ig3, ig4, ije, ije1, ije2, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4, j_5, nnje, shift, tmp, tmp1) 
		#pragma loop name transfb#1#0 
		for (iface=0; iface<6; iface ++ )
		{
			/* nnje=1 for conforming faces, nnje=2 for nonconforming faces */
			if (cbc[ie][iface]==3)
			{
				nnje=2;
			}
			else
			{
				nnje=1;
			}
			/* get collocation point index of four local corners on the face */
			il1=idel[ie][iface][0][0];
			il2=idel[ie][iface][0][5-1];
			il3=idel[ie][iface][5-1][0];
			il4=idel[ie][iface][5-1][5-1];
			/* get the mortar indices of the four local corners */
			ig1=idmo[ie][iface][0][0][0][0];
			ig2=idmo[ie][iface][1][0][0][5-1];
			ig3=idmo[ie][iface][0][1][5-1][0];
			ig4=idmo[ie][iface][1][1][5-1][5-1];
			/* sum the values from tx to tmor for these four local corners */
			/* only 13 of the value is summed, since there will be two duplicated */
			/* contributions from the other two faces sharing this vertex */
			tmor[ig1]=(tmor[ig1]+(tx[il1]*third));
			tmor[ig2]=(tmor[ig2]+(tx[il2]*third));
			tmor[ig3]=(tmor[ig3]+(tx[il3]*third));
			tmor[ig4]=(tmor[ig4]+(tx[il4]*third));
			/* for nonconforming faces */
			if (nnje==2)
			{
				r_init((double * )temp, (5*5)*2, 0.0);
				/* nonconforming faces have four pieces of mortar, first map tx to */
				/* two intermediate mortars stored in temp */
				#pragma cetus private(col, i, i_0, i_1, ije2, il, j, j_0, shift, tmp) 
				#pragma loop name transfb#1#0#0 
				#pragma cetus parallel 
				#pragma omp parallel for if((10000<(1L+(499L*nnje)))) private(col, i, i_0, i_1, ije2, il, j, j_0, shift, tmp)
				for (ije2=0; ije2<nnje; ije2 ++ )
				{
					shift=ije2;
					#pragma cetus private(col, i, i_0, i_1, il, j, j_0, tmp) 
					#pragma loop name transfb#1#0#0#0 
					for (col=0; col<5; col ++ )
					{
						/* For mortar points on face edge (top and bottom), copy the */
						/* value from tx to temp */
						il=idel[ie][iface][v_end[ije2]][col];
						temp[ije2][v_end[ije2]][col]=tx[il];
						/* For mortar points on face edge (top and bottom), calculate */
						/* the interior points' contribution to them, i.e. top() */
						j=v_end[ije2];
						tmp=0.0;
						/* Normalized Loop */
						#pragma cetus private(il) 
						#pragma cetus lastprivate(i_0) 
						#pragma loop name transfb#1#0#0#0#0 
						/* #pragma cetus reduction(+: tmp)  */
						for (i_0=0; i_0<=2; i_0 ++ )
						{
							il=idel[ie][iface][1+i_0][col];
							tmp=(tmp+(qbnew[ije2][j][(1+i_0)-1]*tx[il]));
						}
						i=(1+i_0);
						top[ije2][col]=tmp;
						/* Use mapping matrices qbnew to map the value from tx to temp */
						/* for mortar points not on the top bottom face edge. */
						/* Normalized Loop */
						#pragma cetus private(i, i_1, il, tmp) 
						#pragma cetus lastprivate(j_0) 
						#pragma loop name transfb#1#0#0#0#1 
						for (j_0=0; j_0<=((3+ije2)+(-1*shift)); j_0 ++ )
						{
							tmp=0.0;
							/* Normalized Loop */
							#pragma cetus private(il) 
							#pragma cetus lastprivate(i_1) 
							#pragma loop name transfb#1#0#0#0#1#0 
							/* #pragma cetus reduction(+: tmp)  */
							for (i_1=0; i_1<=2; i_1 ++ )
							{
								il=idel[ie][iface][1+i_1][col];
								tmp=(tmp+(qbnew[ije2][(1+(-1*ije2))+j_0][(1+i_1)-1]*tx[il]));
							}
							i=(1+i_1);
							;
							temp[ije2][(1+(-1*ije2))+j_0][col]=(tmp+temp[ije2][(1+(-1*ije2))+j_0][col]);
						}
						j=((1+(-1*ije2))+j_0);
					}
				}
				/* mapping from temp to tmor */
				#pragma cetus private(col, col_0, i, i_2, i_3, ig, ije1, ije2, j, shift, tmp, tmp1) 
				#pragma loop name transfb#1#0#1 
				for (ije1=0; ije1<nnje; ije1 ++ )
				{
					shift=ije1;
					#pragma cetus private(col, col_0, i, i_2, i_3, ig, ije2, j, tmp, tmp1) 
					#pragma loop name transfb#1#0#1#0 
					for (ije2=0; ije2<nnje; ije2 ++ )
					{
						/* for each column of collocation points on a piece of mortar */
						/* Normalized Loop */
						#pragma cetus private(i, i_2, ig, j, tmp) 
						#pragma cetus lastprivate(col_0) 
						#pragma loop name transfb#1#0#1#0#0 
						for (col_0=0; col_0<=((3+ije1)+(-1*shift)); col_0 ++ )
						{
							/* For the end point, which is on an edge (local edge 1,3), */
							/* the contribution is halved since there will be duplicated */
							/* contribution from another face sharing this edge. */
							ig=idmo[ie][iface][ije2][ije1][(1+col_0)+(-1*ije1)][v_end[ije2]];
							tmor[ig]=(tmor[ig]+(temp[ije1][(1+col_0)+(-1*ije1)][v_end[ije2]]*0.5));
							/* In each row of collocation points on a piece of mortar, */
							/* sum the contributions from interior collocation points */
							/* (i=1,LX1-2) */
							#pragma cetus private(i, i_2, ig, j, tmp) 
							#pragma loop name transfb#1#0#1#0#0#0 
							for (j=0; j<5; j ++ )
							{
								tmp=0.0;
								/* Normalized Loop */
								#pragma cetus lastprivate(i_2) 
								#pragma loop name transfb#1#0#1#0#0#0#0 
								#pragma cetus reduction(+: tmp) 
								#pragma cetus parallel 
								/*
								Disabled due to low profitability: #pragma omp parallel for lastprivate(i_2) reduction(+: tmp)
								*/
								for (i_2=0; i_2<=2; i_2 ++ )
								{
									tmp=(tmp+(qbnew[ije2][j][(1+i_2)-1]*temp[ije1][(1+col_0)+(-1*ije1)][1+i_2]));
								}
								i=(1+i_2);
								ig=idmo[ie][iface][ije2][ije1][(1+col_0)+(-1*ije1)][j];
								tmor[ig]=(tmor[ig]+tmp);
							}
						}
						col=((1+col_0)+(-1*ije1));
						/* For tmor on local edge 0 and 2, tmp is the contribution from */
						/* an edge, so it is halved because of duplicated contribution */
						/* from another face sharing this edge. tmp1 is contribution */
						/* from face interior. */
						col=v_end[ije1];
						ig=idmo[ie][iface][ije2][ije1][col][v_end[ije2]];
						tmor[ig]=(tmor[ig]+(top[ije1][v_end[ije2]]*0.5));
						#pragma cetus private(i, i_3, ig, j, tmp, tmp1) 
						#pragma loop name transfb#1#0#1#0#1 
						for (j=0; j<5; j ++ )
						{
							tmp=0.0;
							tmp1=0.0;
							/* Normalized Loop */
							#pragma cetus lastprivate(i_3) 
							#pragma loop name transfb#1#0#1#0#1#0 
							#pragma cetus reduction(+: tmp, tmp1) 
							#pragma cetus parallel 
							/*
							Disabled due to low profitability: #pragma omp parallel for lastprivate(i_3) reduction(+: tmp, tmp1)
							*/
							for (i_3=0; i_3<=2; i_3 ++ )
							{
								tmp=(tmp+(qbnew[ije2][j][(1+i_3)-1]*temp[ije1][col][1+i_3]));
								tmp1=(tmp1+(qbnew[ije2][j][(1+i_3)-1]*top[ije1][1+i_3]));
							}
							i=(1+i_3);
							ig=idmo[ie][iface][ije2][ije1][col][j];
							tmor[ig]=((tmor[ig]+(tmp*0.5))+tmp1);
						}
					}
				}
				/* for conforming faces */
			}
			else
			{
				/* face interior */
				/* Normalized Loop */
				#pragma cetus private(ig, il, j, j_1) 
				#pragma cetus lastprivate(col_1) 
				#pragma loop name transfb#1#0#2 
				for (col_1=0; col_1<=2; col_1 ++ )
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_1) 
					#pragma loop name transfb#1#0#2#0 
					for (j_1=0; j_1<=2; j_1 ++ )
					{
						il=idel[ie][iface][1+col_1][1+j_1];
						ig=idmo[ie][iface][0][0][1+col_1][1+j_1];
						tmor[ig]=(tmor[ig]+tx[il]);
					}
					j=(1+j_1);
				}
				col=(1+col_1);
				/* edges of conforming faces */
				/* if local edge 0 is a nonconforming edge */
				if (idmo[ie][iface][0][0][0][5-1]!=( - 1))
				{
					#pragma cetus private(i, i_4, ig, ije, il, j, tmp) 
					#pragma loop name transfb#1#0#3 
					for (ije=0; ije<2; ije ++ )
					{
						#pragma cetus private(i, i_4, ig, il, j, tmp) 
						#pragma loop name transfb#1#0#3#0 
						for (j=0; j<5; j ++ )
						{
							tmp=0.0;
							/* Normalized Loop */
							#pragma cetus private(il) 
							#pragma cetus lastprivate(i_4) 
							#pragma loop name transfb#1#0#3#0#0 
							#pragma cetus reduction(+: tmp) 
							#pragma cetus parallel 
							/*
							Disabled due to low profitability: #pragma omp parallel for private(il) lastprivate(i_4) reduction(+: tmp)
							*/
							for (i_4=0; i_4<=2; i_4 ++ )
							{
								il=idel[ie][iface][0][1+i_4];
								tmp=(tmp+(qbnew[ije][j][(1+i_4)-1]*tx[il]));
							}
							i=(1+i_4);
							ig=idmo[ie][iface][ije][0][0][j];
							tmor[ig]=(tmor[ig]+(tmp*0.5));
						}
					}
					/* if local edge 0 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_2) 
					#pragma loop name transfb#1#0#4 
					for (j_2=0; j_2<=2; j_2 ++ )
					{
						il=idel[ie][iface][0][1+j_2];
						ig=idmo[ie][iface][0][0][0][1+j_2];
						tmor[ig]=(tmor[ig]+(tx[il]*0.5));
					}
					j=(1+j_2);
				}
				/* if local edge 1 is a nonconforming edge */
				if (idmo[ie][iface][1][0][1][5-1]!=( - 1))
				{
					#pragma cetus private(i, i_5, ig, ije, il, j, tmp) 
					#pragma loop name transfb#1#0#5 
					for (ije=0; ije<2; ije ++ )
					{
						#pragma cetus private(i, i_5, ig, il, j, tmp) 
						#pragma loop name transfb#1#0#5#0 
						for (j=0; j<5; j ++ )
						{
							tmp=0.0;
							/* Normalized Loop */
							#pragma cetus private(il) 
							#pragma cetus lastprivate(i_5) 
							#pragma loop name transfb#1#0#5#0#0 
							#pragma cetus reduction(+: tmp) 
							#pragma cetus parallel 
							/*
							Disabled due to low profitability: #pragma omp parallel for private(il) lastprivate(i_5) reduction(+: tmp)
							*/
							for (i_5=0; i_5<=2; i_5 ++ )
							{
								il=idel[ie][iface][1+i_5][5-1];
								tmp=(tmp+(qbnew[ije][j][(1+i_5)-1]*tx[il]));
							}
							i=(1+i_5);
							ig=idmo[ie][iface][1][ije][j][5-1];
							tmor[ig]=(tmor[ig]+(tmp*0.5));
						}
					}
					/* if local edge 1 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_3) 
					#pragma loop name transfb#1#0#6 
					for (j_3=0; j_3<=2; j_3 ++ )
					{
						il=idel[ie][iface][1+j_3][5-1];
						ig=idmo[ie][iface][0][0][1+j_3][5-1];
						tmor[ig]=(tmor[ig]+(tx[il]*0.5));
					}
					j=(1+j_3);
				}
				/* if local edge 2 is a nonconforming edge */
				if (idmo[ie][iface][0][1][5-1][1]!=( - 1))
				{
					#pragma cetus private(i, i_6, ig, ije, il, j, tmp) 
					#pragma loop name transfb#1#0#7 
					for (ije=0; ije<2; ije ++ )
					{
						#pragma cetus private(i, i_6, ig, il, j, tmp) 
						#pragma loop name transfb#1#0#7#0 
						for (j=0; j<5; j ++ )
						{
							tmp=0.0;
							/* Normalized Loop */
							#pragma cetus private(il) 
							#pragma cetus lastprivate(i_6) 
							#pragma loop name transfb#1#0#7#0#0 
							#pragma cetus reduction(+: tmp) 
							#pragma cetus parallel 
							/*
							Disabled due to low profitability: #pragma omp parallel for private(il) lastprivate(i_6) reduction(+: tmp)
							*/
							for (i_6=0; i_6<=2; i_6 ++ )
							{
								il=idel[ie][iface][5-1][1+i_6];
								tmp=(tmp+(qbnew[ije][j][(1+i_6)-1]*tx[il]));
							}
							i=(1+i_6);
							ig=idmo[ie][iface][ije][1][5-1][j];
							tmor[ig]=(tmor[ig]+(tmp*0.5));
						}
					}
					/* if local edge 2 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_4) 
					#pragma loop name transfb#1#0#8 
					for (j_4=0; j_4<=2; j_4 ++ )
					{
						il=idel[ie][iface][5-1][1+j_4];
						ig=idmo[ie][iface][0][0][5-1][1+j_4];
						tmor[ig]=(tmor[ig]+(tx[il]*0.5));
					}
					j=(1+j_4);
				}
				/* if local edge 3 is a nonconforming edge */
				if (idmo[ie][iface][0][0][5-1][0]!=( - 1))
				{
					#pragma cetus private(i, i_7, ig, ije, il, j, tmp) 
					#pragma loop name transfb#1#0#9 
					for (ije=0; ije<2; ije ++ )
					{
						#pragma cetus private(i, i_7, ig, il, j, tmp) 
						#pragma loop name transfb#1#0#9#0 
						for (j=0; j<5; j ++ )
						{
							tmp=0.0;
							/* Normalized Loop */
							#pragma cetus private(il) 
							#pragma cetus lastprivate(i_7) 
							#pragma loop name transfb#1#0#9#0#0 
							#pragma cetus reduction(+: tmp) 
							#pragma cetus parallel 
							/*
							Disabled due to low profitability: #pragma omp parallel for private(il) lastprivate(i_7) reduction(+: tmp)
							*/
							for (i_7=0; i_7<=2; i_7 ++ )
							{
								il=idel[ie][iface][1+i_7][0];
								tmp=(tmp+(qbnew[ije][j][(1+i_7)-1]*tx[il]));
							}
							i=(1+i_7);
							ig=idmo[ie][iface][0][ije][j][0];
							tmor[ig]=(tmor[ig]+(tmp*0.5));
						}
					}
					/* if local edge 3 is a conforming edge */
				}
				else
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_5) 
					#pragma loop name transfb#1#0#10 
					for (j_5=0; j_5<=2; j_5 ++ )
					{
						il=idel[ie][iface][1+j_5][0];
						ig=idmo[ie][iface][0][0][1+j_5][0];
						tmor[ig]=(tmor[ig]+(tx[il]*0.5));
					}
					j=(1+j_5);
				}
			}
			/* nnje=1 */
		}
	}
}

/* -------------------------------------------------------------- */
/* This subroutine performs the edge to mortar mapping and */
/* calculates the mapping result on the mortar point at a vertex */
/* under situation 1,2, or 3. */
/* n refers to the configuration of three edges sharing a vertex, */
/* n = 1: only one edge is nonconforming */
/* n = 2: two edges are nonconforming */
/* n = 3: three edges are nonconforming */
/* ------------------------------------------------------------------- */
void transfb_cor_e(int n, double * tmor, double tx[5][5][5])
{
	double tmp;
	int i;
	int i_0;
	int i_1;
	int i_2;
	tmp=tx[0][0][0];
	/* Normalized Loop */
	#pragma cetus lastprivate(i_0) 
	#pragma loop name transfb_cor_e#0 
	#pragma cetus reduction(+: tmp) 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for lastprivate(i_0) reduction(+: tmp)
	*/
	for (i_0=0; i_0<=2; i_0 ++ )
	{
		tmp=(tmp+(qbnew[0][0][(1+i_0)-1]*tx[0][0][1+i_0]));
	}
	i=(1+i_0);
	if (n>1)
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_1) 
		#pragma loop name transfb_cor_e#1 
		#pragma cetus reduction(+: tmp) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_1) reduction(+: tmp)
		*/
		for (i_1=0; i_1<=2; i_1 ++ )
		{
			tmp=(tmp+(qbnew[0][0][(1+i_1)-1]*tx[0][1+i_1][0]));
		}
		i=(1+i_1);
	}
	if (n==3)
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_2) 
		#pragma loop name transfb_cor_e#2 
		#pragma cetus reduction(+: tmp) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_2) reduction(+: tmp)
		*/
		for (i_2=0; i_2<=2; i_2 ++ )
		{
			tmp=(tmp+(qbnew[0][0][(1+i_2)-1]*tx[1+i_2][0][0]));
		}
		i=(1+i_2);
	}
	( * tmor)=tmp;
}

/* -------------------------------------------------------------- */
/* This subroutine performs the mapping from face to mortar. */
/* Output tmor is the mapping result on a mortar vertex */
/* of situations of three edges and three faces sharing a vertex: */
/* n=4: only one face is nonconforming */
/* n=5: one face and one edge are nonconforming */
/* n=6: two faces are nonconforming */
/* n=7: three faces are nonconforming */
/* -------------------------------------------------------------- */
void transfb_cor_f(int n, double * tmor, double tx[5][5][5])
{
	double temp[5], tmp;
	int col, i;
	int i_0;
	int i_1;
	int i_2;
	int i_3;
	int i_4;
	int col_0;
	int i_5;
	int i_6;
	r_init(temp, 5, 0.0);
	#pragma cetus private(col, i, i_0) 
	#pragma loop name transfb_cor_f#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, i_0)
	*/
	for (col=0; col<5; col ++ )
	{
		temp[col]=tx[0][0][col];
		/* Normalized Loop */
		#pragma cetus lastprivate(i_0) 
		#pragma loop name transfb_cor_f#0#0 
		for (i_0=0; i_0<=2; i_0 ++ )
		{
			temp[col]=(temp[col]+(qbnew[0][0][(1+i_0)-1]*tx[0][1+i_0][col]));
		}
		i=(1+i_0);
	}
	tmp=temp[0];
	/* Normalized Loop */
	#pragma cetus lastprivate(i_1) 
	#pragma loop name transfb_cor_f#1 
	#pragma cetus reduction(+: tmp) 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for lastprivate(i_1) reduction(+: tmp)
	*/
	for (i_1=0; i_1<=2; i_1 ++ )
	{
		tmp=(tmp+(qbnew[0][0][(1+i_1)-1]*temp[1+i_1]));
	}
	i=(1+i_1);
	if (n==5)
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_2) 
		#pragma loop name transfb_cor_f#2 
		#pragma cetus reduction(+: tmp) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_2) reduction(+: tmp)
		*/
		for (i_2=0; i_2<=2; i_2 ++ )
		{
			tmp=(tmp+(qbnew[0][0][(1+i_2)-1]*tx[1+i_2][0][0]));
		}
		i=(1+i_2);
	}
	if (n>=6)
	{
		r_init(temp, 5, 0.0);
		#pragma cetus private(col, i, i_3) 
		#pragma loop name transfb_cor_f#3 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(col, i, i_3)
		*/
		for (col=0; col<5; col ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_3) 
			#pragma loop name transfb_cor_f#3#0 
			for (i_3=0; i_3<=2; i_3 ++ )
			{
				temp[col]=(temp[col]+(qbnew[0][0][(1+i_3)-1]*tx[1+i_3][0][col]));
			}
			i=(1+i_3);
		}
		tmp=(tmp+temp[0]);
		/* Normalized Loop */
		#pragma cetus lastprivate(i_4) 
		#pragma loop name transfb_cor_f#4 
		#pragma cetus reduction(+: tmp) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_4) reduction(+: tmp)
		*/
		for (i_4=0; i_4<=2; i_4 ++ )
		{
			tmp=(tmp+(qbnew[0][0][(1+i_4)-1]*temp[1+i_4]));
		}
		i=(1+i_4);
	}
	if (n==7)
	{
		r_init(temp, 5, 0.0);
		/* Normalized Loop */
		#pragma cetus private(i, i_5) 
		#pragma cetus lastprivate(col_0) 
		#pragma loop name transfb_cor_f#5 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i, i_5) lastprivate(col_0)
		*/
		for (col_0=0; col_0<=2; col_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_5) 
			#pragma loop name transfb_cor_f#5#0 
			for (i_5=0; i_5<=2; i_5 ++ )
			{
				temp[1+col_0]=(temp[1+col_0]+(qbnew[0][0][(1+i_5)-1]*tx[1+i_5][1+col_0][0]));
			}
			i=(1+i_5);
		}
		col=(1+col_0);
		/* Normalized Loop */
		#pragma cetus lastprivate(i_6) 
		#pragma loop name transfb_cor_f#6 
		#pragma cetus reduction(+: tmp) 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for lastprivate(i_6) reduction(+: tmp)
		*/
		for (i_6=0; i_6<=2; i_6 ++ )
		{
			tmp=(tmp+(qbnew[0][0][(1+i_6)-1]*temp[1+i_6]));
		}
		i=(1+i_6);
	}
	( * tmor)=tmp;
}

/* ------------------------------------------------------------------------ */
/* Perform mortar to element mapping on a nonconforming face. */
/* This subroutin is used when all entries in tmor are zero except */
/* one tmor[j][i]=1. So this routine is simplified. Only one piece of */
/* mortar  (tmor only has two indices) and one piece of intermediate */
/* mortar (tmp) are involved. */
/* ------------------------------------------------------------------------ */
void transf_nc(double tmor[5][5], double tx[5][5])
{
	double tmp[5][5];
	int col, i, j;
	int i_0;
	int i_1;
	r_init((double * )tmp, 5*5, 0.0);
	#pragma cetus private(col, i, i_0, j) 
	#pragma loop name transf_nc#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, i_0, j)
	*/
	for (col=0; col<5; col ++ )
	{
		i=0;
		tmp[col][i]=tmor[col][i];
		/* Normalized Loop */
		#pragma cetus private(j) 
		#pragma cetus lastprivate(i_0) 
		#pragma loop name transf_nc#0#0 
		for (i_0=0; i_0<=2; i_0 ++ )
		{
			#pragma cetus private(j) 
			#pragma loop name transf_nc#0#0#0 
			for (j=0; j<5; j ++ )
			{
				tmp[col][1+i_0]=(tmp[col][1+i_0]+(qbnew[0][j][(1+i_0)-1]*tmor[col][j]));
			}
		}
		i=(1+i_0);
	}
	#pragma cetus private(col, i, i_1, j) 
	#pragma loop name transf_nc#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, i_1, j)
	*/
	for (col=0; col<5; col ++ )
	{
		i=0;
		tx[i][col]=(tx[i][col]+tmp[i][col]);
		/* Normalized Loop */
		#pragma cetus private(j) 
		#pragma cetus lastprivate(i_1) 
		#pragma loop name transf_nc#1#0 
		for (i_1=0; i_1<=2; i_1 ++ )
		{
			#pragma cetus private(j) 
			#pragma loop name transf_nc#1#0#0 
			for (j=0; j<5; j ++ )
			{
				tx[1+i_1][col]=(tx[1+i_1][col]+(qbnew[0][j][(1+i_1)-1]*tmp[j][col]));
			}
		}
		i=(1+i_1);
	}
}

/* ------------------------------------------------------------------------ */
/* Performs mapping from element to mortar when the nonconforming */
/* edges are shared by two conforming faces of an element. */
/* ------------------------------------------------------------------------ */
void transfb_nc0(double tmor[5][5], double tx[5][5][5])
{
	int i, j;
	int i_0;
	r_init((double * )tmor, 5*5, 0.0);
	#pragma cetus private(i, i_0, j) 
	#pragma loop name transfb_nc0#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, i_0, j)
	*/
	for (j=0; j<5; j ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_0) 
		#pragma loop name transfb_nc0#0#0 
		for (i_0=0; i_0<=2; i_0 ++ )
		{
			tmor[0][j]=(tmor[0][j]+(qbnew[0][j][(1+i_0)-1]*tx[0][0][1+i_0]));
		}
		i=(1+i_0);
	}
}

/* ------------------------------------------------------------------------ */
/* Maps values from element to mortar when the nonconforming edges are */
/* shared by two nonconforming faces of an element. */
/* Although each face shall have four pieces of mortar, only value in */
/* one piece (location (0,0)) is used in the calling routine so only */
/* the value in the first mortar is calculated in this subroutine. */
/* ------------------------------------------------------------------------ */
void transfb_nc2(double tmor[5][5], double tx[5][5])
{
	double bottom[5], temp[5][5];
	int col, j, i;
	int i_0;
	int j_0;
	int i_1;
	int i_2;
	int col_0;
	int i_3;
	r_init((double * )tmor, 5*5, 0.0);
	r_init((double * )temp, 5*5, 0.0);
	tmor[0][0]=tx[0][0];
	/* mapping from tx to intermediate mortar temp + bottom */
	#pragma cetus private(col, i, i_0, i_1, j, j_0) 
	#pragma loop name transfb_nc2#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, i_0, i_1, j, j_0)
	*/
	for (col=0; col<5; col ++ )
	{
		temp[0][col]=tx[0][col];
		j=0;
		bottom[col]=0.0;
		;
		/* Normalized Loop */
		#pragma cetus lastprivate(i_0) 
		#pragma loop name transfb_nc2#0#0 
		for (i_0=0; i_0<=2; i_0 ++ )
		{
			bottom[col]=(bottom[col]+(qbnew[0][j][(1+i_0)-1]*tx[1+i_0][col]));
		}
		i=(1+i_0);
		/* Normalized Loop */
		#pragma cetus private(i, i_1) 
		#pragma cetus lastprivate(j_0) 
		#pragma loop name transfb_nc2#0#1 
		for (j_0=0; j_0<=3; j_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_1) 
			#pragma loop name transfb_nc2#0#1#0 
			for (i_1=0; i_1<=2; i_1 ++ )
			{
				temp[1+j_0][col]=(temp[1+j_0][col]+(qbnew[0][1+j_0][(1+i_1)-1]*tx[1+i_1][col]));
			}
			i=(1+i_1);
		}
		j=(1+j_0);
	}
	/* from intermediate mortar to mortar */
	/* On the nonconforming edge, temp is divided by 2 as there will be */
	/* a duplicate contribution from another face sharing this edge */
	col=0;
	#pragma cetus private(i, i_2, j) 
	#pragma loop name transfb_nc2#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, i_2, j)
	*/
	for (j=0; j<5; j ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_2) 
		#pragma loop name transfb_nc2#1#0 
		for (i_2=0; i_2<=2; i_2 ++ )
		{
			tmor[col][j]=((tmor[col][j]+(qbnew[0][j][(1+i_2)-1]*bottom[1+i_2]))+((qbnew[0][j][(1+i_2)-1]*temp[col][1+i_2])*0.5));
		}
		i=(1+i_2);
	}
	/* Normalized Loop */
	#pragma cetus private(i, i_3, j) 
	#pragma cetus lastprivate(col_0) 
	#pragma loop name transfb_nc2#2 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, i_3, j) lastprivate(col_0)
	*/
	for (col_0=0; col_0<=3; col_0 ++ )
	{
		tmor[1+col_0][0]=(tmor[1+col_0][0]+temp[1+col_0][0]);
		#pragma cetus private(i, i_3, j) 
		#pragma loop name transfb_nc2#2#0 
		for (j=0; j<5; j ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_3) 
			#pragma loop name transfb_nc2#2#0#0 
			for (i_3=0; i_3<=2; i_3 ++ )
			{
				tmor[1+col_0][j]=(tmor[1+col_0][j]+(qbnew[0][j][(1+i_3)-1]*temp[1+col_0][1+i_3]));
			}
			i=(1+i_3);
		}
	}
	col=(1+col_0);
}

/* ------------------------------------------------------------------------ */
/* Maps values from element to mortar when the nonconforming edges are */
/* shared by a nonconforming face and a conforming face of an element */
/* ------------------------------------------------------------------------ */
void transfb_nc1(double tmor[5][5], double tx[5][5])
{
	double bottom[5], temp[5][5];
	int col, j, i;
	int i_0;
	int j_0;
	int i_1;
	int i_2;
	int col_0;
	int i_3;
	r_init((double * )tmor, 5*5, 0.0);
	r_init((double * )temp, 5*5, 0.0);
	tmor[0][0]=tx[0][0];
	/* Contribution from the nonconforming faces */
	/* Since the calling subroutine is only interested in the value on the */
	/* mortar (location (0,0)), only this piece of mortar is calculated. */
	#pragma cetus private(col, i, i_0, i_1, j, j_0) 
	#pragma loop name transfb_nc1#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(col, i, i_0, i_1, j, j_0)
	*/
	for (col=0; col<5; col ++ )
	{
		temp[0][col]=tx[0][col];
		j=0;
		bottom[col]=0.0;
		/* Normalized Loop */
		#pragma cetus lastprivate(i_0) 
		#pragma loop name transfb_nc1#0#0 
		for (i_0=0; i_0<=2; i_0 ++ )
		{
			bottom[col]=(bottom[col]+(qbnew[0][j][(1+i_0)-1]*tx[1+i_0][col]));
		}
		i=(1+i_0);
		/* Normalized Loop */
		#pragma cetus private(i, i_1) 
		#pragma cetus lastprivate(j_0) 
		#pragma loop name transfb_nc1#0#1 
		for (j_0=0; j_0<=3; j_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_1) 
			#pragma loop name transfb_nc1#0#1#0 
			for (i_1=0; i_1<=2; i_1 ++ )
			{
				temp[1+j_0][col]=(temp[1+j_0][col]+(qbnew[0][1+j_0][(1+i_1)-1]*tx[1+i_1][col]));
			}
			i=(1+i_1);
		}
		j=(1+j_0);
	}
	col=0;
	tmor[col][0]=(tmor[col][0]+bottom[0]);
	#pragma cetus private(i, i_2, j) 
	#pragma loop name transfb_nc1#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, i_2, j)
	*/
	for (j=0; j<5; j ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i_2) 
		#pragma loop name transfb_nc1#1#0 
		for (i_2=0; i_2<=2; i_2 ++ )
		{
			/* temp is not divided by 2 here. It includes the contribution */
			/* from the other conforming face. */
			tmor[col][j]=((tmor[col][j]+(qbnew[0][j][(1+i_2)-1]*bottom[1+i_2]))+(qbnew[0][j][(1+i_2)-1]*temp[col][1+i_2]));
		}
		i=(1+i_2);
	}
	/* Normalized Loop */
	#pragma cetus private(i, i_3, j) 
	#pragma cetus lastprivate(col_0) 
	#pragma loop name transfb_nc1#2 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(i, i_3, j) lastprivate(col_0)
	*/
	for (col_0=0; col_0<=3; col_0 ++ )
	{
		tmor[1+col_0][0]=(tmor[1+col_0][0]+temp[1+col_0][0]);
		#pragma cetus private(i, i_3, j) 
		#pragma loop name transfb_nc1#2#0 
		for (j=0; j<5; j ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i_3) 
			#pragma loop name transfb_nc1#2#0#0 
			for (i_3=0; i_3<=2; i_3 ++ )
			{
				tmor[1+col_0][j]=(tmor[1+col_0][j]+(qbnew[0][j][(1+i_3)-1]*temp[1+col_0][1+i_3]));
			}
			i=(1+i_3);
		}
	}
	col=(1+col_0);
}

/* ------------------------------------------------------------------- */
/* Prepare initial guess for cg. All values from conforming */
/* boundary are copied and summed on tmor. */
/* ------------------------------------------------------------------- */
void transfb_c(double tx[])
{
	const double third = 1.0/3.0;
	int il1, il2, il3, il4, ig1, ig2, ig3, ig4, ie, iface, col, j, ig, il;
	int col_0;
	int j_0;
	int j_1;
	int j_2;
	int j_3;
	int j_4;
	r_init(tmort, nmor, 0.0);
	#pragma cetus private(col, col_0, ie, iface, ig, ig1, ig2, ig3, ig4, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4) 
	#pragma loop name transfb_c#0 
	for (ie=0; ie<nelt; ie ++ )
	{
		#pragma cetus private(col, col_0, iface, ig, ig1, ig2, ig3, ig4, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4) 
		#pragma loop name transfb_c#0#0 
		for (iface=0; iface<6; iface ++ )
		{
			if (cbc[ie][iface]!=3)
			{
				il1=idel[ie][iface][0][0];
				il2=idel[ie][iface][0][5-1];
				il3=idel[ie][iface][5-1][0];
				il4=idel[ie][iface][5-1][5-1];
				ig1=idmo[ie][iface][0][0][0][0];
				ig2=idmo[ie][iface][1][0][0][5-1];
				ig3=idmo[ie][iface][0][1][5-1][0];
				ig4=idmo[ie][iface][1][1][5-1][5-1];
				tmort[ig1]=(tmort[ig1]+(tx[il1]*third));
				tmort[ig2]=(tmort[ig2]+(tx[il2]*third));
				tmort[ig3]=(tmort[ig3]+(tx[il3]*third));
				tmort[ig4]=(tmort[ig4]+(tx[il4]*third));
				/* Normalized Loop */
				#pragma cetus private(ig, il, j, j_0) 
				#pragma cetus lastprivate(col_0) 
				#pragma loop name transfb_c#0#0#0 
				for (col_0=0; col_0<=2; col_0 ++ )
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_0) 
					#pragma loop name transfb_c#0#0#0#0 
					for (j_0=0; j_0<=2; j_0 ++ )
					{
						il=idel[ie][iface][1+col_0][1+j_0];
						ig=idmo[ie][iface][0][0][1+col_0][1+j_0];
						tmort[ig]=(tmort[ig]+tx[il]);
					}
					j=(1+j_0);
				}
				col=(1+col_0);
				if (idmo[ie][iface][0][0][0][5-1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_1) 
					#pragma loop name transfb_c#0#0#1 
					for (j_1=0; j_1<=2; j_1 ++ )
					{
						il=idel[ie][iface][0][1+j_1];
						ig=idmo[ie][iface][0][0][0][1+j_1];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
					}
					j=(1+j_1);
				}
				if (idmo[ie][iface][1][0][1][5-1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_2) 
					#pragma loop name transfb_c#0#0#2 
					for (j_2=0; j_2<=2; j_2 ++ )
					{
						il=idel[ie][iface][1+j_2][5-1];
						ig=idmo[ie][iface][0][0][1+j_2][5-1];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
					}
					j=(1+j_2);
				}
				if (idmo[ie][iface][0][1][5-1][1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_3) 
					#pragma loop name transfb_c#0#0#3 
					for (j_3=0; j_3<=2; j_3 ++ )
					{
						il=idel[ie][iface][5-1][1+j_3];
						ig=idmo[ie][iface][0][0][5-1][1+j_3];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
					}
					j=(1+j_3);
				}
				if (idmo[ie][iface][0][0][5-1][0]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_4) 
					#pragma loop name transfb_c#0#0#4 
					for (j_4=0; j_4<=2; j_4 ++ )
					{
						il=idel[ie][iface][1+j_4][0];
						ig=idmo[ie][iface][0][0][1+j_4][0];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
					}
					j=(1+j_4);
				}
			}
		}
	}
}

/* ------------------------------------------------------------------- */
/* Prepare initial guess for CG. All values from conforming */
/* boundary are copied and summed in tmort. */
/* mormult is multiplicity, which is used to average tmort. */
/* ------------------------------------------------------------------- */
void transfb_c_2(double tx[])
{
	const double third = 1.0/3.0;
	int il1, il2, il3, il4, ig1, ig2, ig3, ig4, ie, iface, col, j, ig, il;
	int col_0;
	int j_0;
	int j_1;
	int j_2;
	int j_3;
	int j_4;
	r_init(tmort, nmor, 0.0);
	r_init(mormult, nmor, 0.0);
	#pragma cetus private(col, col_0, ie, iface, ig, ig1, ig2, ig3, ig4, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4) 
	#pragma loop name transfb_c_2#0 
	for (ie=0; ie<nelt; ie ++ )
	{
		#pragma cetus private(col, col_0, iface, ig, ig1, ig2, ig3, ig4, il, il1, il2, il3, il4, j, j_0, j_1, j_2, j_3, j_4) 
		#pragma loop name transfb_c_2#0#0 
		for (iface=0; iface<6; iface ++ )
		{
			if (cbc[ie][iface]!=3)
			{
				il1=idel[ie][iface][0][0];
				il2=idel[ie][iface][0][5-1];
				il3=idel[ie][iface][5-1][0];
				il4=idel[ie][iface][5-1][5-1];
				ig1=idmo[ie][iface][0][0][0][0];
				ig2=idmo[ie][iface][1][0][0][5-1];
				ig3=idmo[ie][iface][0][1][5-1][0];
				ig4=idmo[ie][iface][1][1][5-1][5-1];
				tmort[ig1]=(tmort[ig1]+(tx[il1]*third));
				tmort[ig2]=(tmort[ig2]+(tx[il2]*third));
				tmort[ig3]=(tmort[ig3]+(tx[il3]*third));
				tmort[ig4]=(tmort[ig4]+(tx[il4]*third));
				mormult[ig1]=(mormult[ig1]+third);
				mormult[ig2]=(mormult[ig2]+third);
				mormult[ig3]=(mormult[ig3]+third);
				mormult[ig4]=(mormult[ig4]+third);
				/* Normalized Loop */
				#pragma cetus private(ig, il, j, j_0) 
				#pragma cetus lastprivate(col_0) 
				#pragma loop name transfb_c_2#0#0#0 
				for (col_0=0; col_0<=2; col_0 ++ )
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_0) 
					#pragma loop name transfb_c_2#0#0#0#0 
					for (j_0=0; j_0<=2; j_0 ++ )
					{
						il=idel[ie][iface][1+col_0][1+j_0];
						ig=idmo[ie][iface][0][0][1+col_0][1+j_0];
						tmort[ig]=(tmort[ig]+tx[il]);
						mormult[ig]=(mormult[ig]+1.0);
					}
					j=(1+j_0);
				}
				col=(1+col_0);
				if (idmo[ie][iface][0][0][0][5-1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_1) 
					#pragma loop name transfb_c_2#0#0#1 
					for (j_1=0; j_1<=2; j_1 ++ )
					{
						il=idel[ie][iface][0][1+j_1];
						ig=idmo[ie][iface][0][0][0][1+j_1];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
						mormult[ig]=(mormult[ig]+0.5);
					}
					j=(1+j_1);
				}
				if (idmo[ie][iface][1][0][1][5-1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_2) 
					#pragma loop name transfb_c_2#0#0#2 
					for (j_2=0; j_2<=2; j_2 ++ )
					{
						il=idel[ie][iface][1+j_2][5-1];
						ig=idmo[ie][iface][0][0][1+j_2][5-1];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
						mormult[ig]=(mormult[ig]+0.5);
					}
					j=(1+j_2);
				}
				if (idmo[ie][iface][0][1][5-1][1]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_3) 
					#pragma loop name transfb_c_2#0#0#3 
					for (j_3=0; j_3<=2; j_3 ++ )
					{
						il=idel[ie][iface][5-1][1+j_3];
						ig=idmo[ie][iface][0][0][5-1][1+j_3];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
						mormult[ig]=(mormult[ig]+0.5);
					}
					j=(1+j_3);
				}
				if (idmo[ie][iface][0][0][5-1][0]==( - 1))
				{
					/* Normalized Loop */
					#pragma cetus private(ig, il) 
					#pragma cetus lastprivate(j_4) 
					#pragma loop name transfb_c_2#0#0#4 
					for (j_4=0; j_4<=2; j_4 ++ )
					{
						il=idel[ie][iface][1+j_4][0];
						ig=idmo[ie][iface][0][0][1+j_4][0];
						tmort[ig]=(tmort[ig]+(tx[il]*0.5));
						mormult[ig]=(mormult[ig]+0.5);
					}
					j=(1+j_4);
				}
			}
		}
	}
}
