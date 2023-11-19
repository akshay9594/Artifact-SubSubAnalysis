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
/* ------------------------------------------------------------------------- */
/*                                                                          */
/*  This benchmark is a serial C version of the NPB MG code. This C         */
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
/* --------------------------------------------------------------------- */
/*  program mg */
/* --------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "globals.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"
static void setup(int * n1, int * n2, int * n3);
static void mg3P(double u[], double v[], double r[], double a[4], double c[4], int n1, int n2, int n3);
static void psinv(void * or, void * ou, int n1, int n2, int n3, double c[4], int k);
static void resid(void * ou, void * ov, void * or, int n1, int n2, int n3, double a[4], int k);
static void rprj3(void * or, int m1k, int m2k, int m3k, void * os, int m1j, int m2j, int m3j, int k);
static void interp(void * oz, int mm1, int mm2, int mm3, void * ou, int n1, int n2, int n3, int k);
static void norm2u3(void * or, int n1, int n2, int n3, double * rnm2, double * rnmu, int nx, int ny, int nz);
static void rep_nrm(void * u, int n1, int n2, int n3, char * title, int kk);
static void comm3(void * ou, int n1, int n2, int n3, int kk);
static void zran3(void * oz, int n1, int n2, int n3, int nx, int ny, int k);
static void showall(void * oz, int n1, int n2, int n3);
static double power(double a, int n);
static void bubble(double ten[][2], int j1[][2], int j2[][2], int j3[][2], int m, int ind);
static void zero3(void * oz, int n1, int n2, int n3);
/* -------------------------------------------------------------------------c */
/* These arrays are in common because they are quite large */
/* and probably shouldn't be allocated on the stack. They */
/* are always passed as subroutine args.  */
/* -------------------------------------------------------------------------c */
/* commconnoautom */
static double u[(((((((((1*(2+(1<<8)))*(2+(1<<8)))*(2+(1<<8)))+((2+(1<<8))*(2+(1<<8))))+(5*(2+(1<<8))))+(7*8))+6)/7)*8)];
static double v[(((((((((1*(2+(1<<8)))*(2+(1<<8)))*(2+(1<<8)))+((2+(1<<8))*(2+(1<<8))))+(5*(2+(1<<8))))+(7*8))+6)/7)*8)];
static double r[(((((((((1*(2+(1<<8)))*(2+(1<<8)))*(2+(1<<8)))+((2+(1<<8))*(2+(1<<8))))+(5*(2+(1<<8))))+(7*8))+6)/7)*8)];
/* commongrid */
static int is1, is2, is3, ie1, ie2, ie3;
int main()
{
	/* -------------------------------------------------------------------------c */
	/* k is the current level. It is passed down through subroutine args */
	/* and is NOT global. it is the current iteration */
	/* -------------------------------------------------------------------------c */
	int k, it;
	double t, tinit, mflops;
	double a[4], c[4];
	double rnm2, rnmu, old2, oldu, epsilon;
	int n1, n2, n3, nit;
	double nn, verify_value, err;
	logical verified;
	int i;
	char * t_names[10];
	double tmax;
	FILE * fp;
	int i_0;
	int it_0;
	int i_1;
	#pragma cetus private(i) 
	#pragma loop name main#0 
	for (i=0; i<10; i ++ )
	{
		timer_clear(i);
	}
	timer_start(0);
	/* --------------------------------------------------------------------- */
	/* Read in and broadcast input data */
	/* --------------------------------------------------------------------- */
	if ((fp=fopen("timer.flag", "r"))!=((void * )0))
	{
		timeron=true;
		t_names[0]="init";
		t_names[1]="benchmk";
		t_names[2]="mg3P";
		t_names[3]="psinv";
		t_names[4]="resid";
		t_names[6]="rprj3";
		t_names[7]="interp";
		t_names[8]="norm2";
		t_names[9]="comm3";
		fclose(fp);
	}
	else
	{
		timeron=false;
	}
	printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - MG Benchmark\n\n");
	if ((fp=fopen("mg.input", "r"))!=((void * )0))
	{
		int result;
		printf(" Reading from input file mg.input\n");
		result=fscanf(fp, "%d\n",  & lt);
		while (fgetc(fp)!='\n')
		{
			;
		}
		result=fscanf(fp, "%d%d%d",  & nx[lt],  & ny[lt],  & nz[lt]);
		while (fgetc(fp)!='\n')
		{
			;
		}
		result=fscanf(fp, "%d",  & nit);
		while (fgetc(fp)!='\n')
		{
			;
		}
		#pragma cetus private(i, result) 
		#pragma loop name main#1 
		for (i=0; i<=7; i ++ )
		{
			result=fscanf(fp, "%d",  & debug_vec[i]);
		}
		fclose(fp);
	}
	else
	{
		printf(" No input file. Using compiled defaults \n");
		lt=8;
		nit=20;
		nx[lt]=256;
		ny[lt]=256;
		nz[lt]=256;
		#pragma cetus private(i) 
		#pragma loop name main#2 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(i)
		*/
		for (i=0; i<=7; i ++ )
		{
			debug_vec[i]=0;
		}
	}
	if ((nx[lt]!=ny[lt])||(nx[lt]!=nz[lt]))
	{
		Class='U';
	}
	else
	{
		if ((nx[lt]==32)&&(nit==4))
		{
			Class='S';
		}
		else
		{
			if ((nx[lt]==128)&&(nit==4))
			{
				Class='W';
			}
			else
			{
				if ((nx[lt]==256)&&(nit==4))
				{
					Class='A';
				}
				else
				{
					if ((nx[lt]==256)&&(nit==20))
					{
						Class='B';
					}
					else
					{
						if ((nx[lt]==512)&&(nit==20))
						{
							Class='C';
						}
						else
						{
							if ((nx[lt]==1024)&&(nit==50))
							{
								Class='D';
							}
							else
							{
								if ((nx[lt]==2048)&&(nit==50))
								{
									Class='E';
								}
								else
								{
									Class='U';
								}
							}
						}
					}
				}
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* Use these for debug info: */
	/* --------------------------------------------------------------------- */
	/*    debug_vec(0) = 1 !=> report all norms */
	/*    debug_vec(1) = 1 !=> some setup information */
	/*    debug_vec(1) = 2 !=> more setup information */
	/*    debug_vec(2) = k => at level k or below, show result of resid */
	/*    debug_vec(3) = k => at level k or below, show result of psinv */
	/*    debug_vec(4) = k => at level k or below, show result of rprj */
	/*    debug_vec(5) = k => at level k or below, show result of interp */
	/*    debug_vec(6) = 1 => (unused) */
	/*    debug_vec(7) = 1 => (unused) */
	/* --------------------------------------------------------------------- */
	a[0]=(( - 8.0)/3.0);
	a[1]=0.0;
	a[2]=(1.0/6.0);
	a[3]=(1.0/12.0);
	if (((Class=='A')||(Class=='S'))||(Class=='W'))
	{
		/* --------------------------------------------------------------------- */
		/* Coefficients for the S(a) smoother */
		/* --------------------------------------------------------------------- */
		c[0]=(( - 3.0)/8.0);
		c[1]=(( + 1.0)/32.0);
		c[2]=(( - 1.0)/64.0);
		c[3]=0.0;
	}
	else
	{
		/* --------------------------------------------------------------------- */
		/* Coefficients for the S(b) smoother */
		/* --------------------------------------------------------------------- */
		c[0]=(( - 3.0)/17.0);
		c[1]=(( + 1.0)/33.0);
		c[2]=(( - 1.0)/61.0);
		c[3]=0.0;
	}
	lb=1;
	k=lt;
	setup( & n1,  & n2,  & n3);
	zero3(u, n1, n2, n3);
	zran3(v, n1, n2, n3, nx[lt], ny[lt], k);
	norm2u3(v, n1, n2, n3,  & rnm2,  & rnmu, nx[lt], ny[lt], nz[lt]);
	printf(" Size: %4dx%4dx%4d  (class %c)\n", nx[lt], ny[lt], nz[lt], Class);
	printf(" Iterations: %3d\n", nit);
	printf("\n");
	resid(u, v, r, n1, n2, n3, a, k);
	norm2u3(r, n1, n2, n3,  & rnm2,  & rnmu, nx[lt], ny[lt], nz[lt]);
	old2=rnm2;
	oldu=rnmu;
	/* --------------------------------------------------------------------- */
	/* One iteration for startup */
	/* --------------------------------------------------------------------- */
	mg3P(u, v, r, a, c, n1, n2, n3);
	resid(u, v, r, n1, n2, n3, a, k);
	setup( & n1,  & n2,  & n3);
	zero3(u, n1, n2, n3);
	zran3(v, n1, n2, n3, nx[lt], ny[lt], k);
	timer_stop(0);
	tinit=timer_read(0);
	printf(" Initialization time: %15.3f seconds\n\n", tinit);
	/* Normalized Loop */
	#pragma cetus lastprivate(i_0) 
	#pragma loop name main#3 
	for (i_0=0; i_0<=8; i_0 ++ )
	{
		timer_clear(1+i_0);
	}
	i=(1+i_0);
	timer_start(1);
	if (timeron)
	{
		timer_start(5);
	}
	resid(u, v, r, n1, n2, n3, a, k);
	if (timeron)
	{
		timer_stop(5);
	}
	norm2u3(r, n1, n2, n3,  & rnm2,  & rnmu, nx[lt], ny[lt], nz[lt]);
	old2=rnm2;
	oldu=rnmu;
	/* Normalized Loop */
	#pragma cetus lastprivate(it_0) 
	#pragma loop name main#4 
	for (it_0=0; it_0<=(-1+nit); it_0 ++ )
	{
		if ((((1+it_0)==1)||((1+it_0)==nit))||(((1+it_0)%5)==0))
		{
			printf("  iter %3d\n", 1+it_0);
		}
		if (timeron)
		{
			timer_start(2);
		}
		mg3P(u, v, r, a, c, n1, n2, n3);
		if (timeron)
		{
			timer_stop(2);
		}
		if (timeron)
		{
			timer_start(5);
		}
		resid(u, v, r, n1, n2, n3, a, k);
		if (timeron)
		{
			timer_stop(5);
		}
	}
	it=(1+it_0);
	norm2u3(r, n1, n2, n3,  & rnm2,  & rnmu, nx[lt], ny[lt], nz[lt]);
	timer_stop(1);
	t=timer_read(1);
	verified=false;
	verify_value=0.0;
	printf("\n Benchmark completed\n");
	epsilon=1.0E-8;
	if (Class!='U')
	{
		if (Class=='S')
		{
			verify_value=5.307707005734E-5;
		}
		else
		{
			if (Class=='W')
			{
				verify_value=6.467329375339E-6;
			}
			else
			{
				if (Class=='A')
				{
					verify_value=2.433365309069E-6;
				}
				else
				{
					if (Class=='B')
					{
						verify_value=1.800564401355E-6;
					}
					else
					{
						if (Class=='C')
						{
							verify_value=5.70673228574E-7;
						}
						else
						{
							if (Class=='D')
							{
								verify_value=1.58327506044E-10;
							}
							else
							{
								if (Class=='E')
								{
									verify_value=8.157592357404E-11;
								}
							}
						}
					}
				}
			}
		}
		err=(fabs(rnm2-verify_value)/verify_value);
		/* err = fabs( rnm2 - verify_value ); */
		if (err<=epsilon)
		{
			verified=true;
			printf(" VERIFICATION SUCCESSFUL\n");
			printf(" L2 Norm is %20.13E\n", rnm2);
			printf(" Error is   %20.13E\n", err);
		}
		else
		{
			verified=false;
			printf(" VERIFICATION FAILED\n");
			printf(" L2 Norm is             %20.13E\n", rnm2);
			printf(" The correct L2 Norm is %20.13E\n", verify_value);
		}
	}
	else
	{
		verified=false;
		printf(" Problem size unknown\n");
		printf(" NO VERIFICATION PERFORMED\n");
		printf(" L2 Norm is %20.13E\n", rnm2);
	}
	nn=(((1.0*nx[lt])*ny[lt])*nz[lt]);
	if (t!=0.0)
	{
		mflops=((((58.0*nit)*nn)*1.0E-6)/t);
	}
	else
	{
		mflops=0.0;
	}
	print_results("MG", Class, nx[lt], ny[lt], nz[lt], nit, t, mflops, "          floating point", verified, "3.3.1", "19 Nov 2023", "gcc", "$(CC)", "-lm", "-I../common", "-g -Wall -O3 -fopenmp -mcmodel=large", "-g -O3 -fopenmp -mcmodel=large", "randdp");
	/* --------------------------------------------------------------------- */
	/* More timers */
	/* --------------------------------------------------------------------- */
	if (timeron)
	{
		tmax=timer_read(1);
		if (tmax==0.0)
		{
			tmax=1.0;
		}
		printf("  SECTION   Time (secs)\n");
		/* Normalized Loop */
		#pragma cetus private(t) 
		#pragma cetus lastprivate(i_1) 
		#pragma loop name main#5 
		for (i_1=0; i_1<=8; i_1 ++ )
		{
			t=timer_read(1+i_1);
			if ((1+i_1)==5)
			{
				t=(timer_read(4)-t);
				printf("    --> %8s:%9.3f  (%6.2f%%)\n", "mg-resid", t, (t*100.0)/tmax);
			}
			else
			{
				printf("  %-8s:%9.3f  (%6.2f%%)\n", t_names[1+i_1], t, (t*100.0)/tmax);
			}
		}
		i=(1+i_1);
	}
	return 0;
}

static void setup(int * n1, int * n2, int * n3)
{
	int k, j;
	int ax, mi[((8+1)+1)][3];
	int ng[((8+1)+1)][3];
	int k_0;
	int k_1;
	int k_2;
	int j_0;
	ng[lt][0]=nx[lt];
	ng[lt][1]=ny[lt];
	ng[lt][2]=nz[lt];
	/* Normalized Loop */
	#pragma cetus private(ax) 
	#pragma cetus lastprivate(k_0) 
	#pragma loop name setup#0 
	for (k_0=0; k_0<=(-2+lt); k_0 ++ )
	{
		#pragma cetus private(ax) 
		#pragma loop name setup#0#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(ax)
		*/
		for (ax=0; ax<3; ax ++ )
		{
			ng[(-1+(-1*k_0))+lt][ax]=(ng[((-1+(-1*k_0))+lt)+1][ax]/2);
		}
	}
	k=((-1+(-1*k_0))+lt);
	/* Normalized Loop */
	#pragma cetus lastprivate(k_1) 
	#pragma loop name setup#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(5L*lt)))) lastprivate(k_1)
	for (k_1=0; k_1<=(-1+lt); k_1 ++ )
	{
		nx[(-1*k_1)+lt]=ng[(-1*k_1)+lt][0];
		ny[(-1*k_1)+lt]=ng[(-1*k_1)+lt][1];
		nz[(-1*k_1)+lt]=ng[(-1*k_1)+lt][2];
	}
	k=((-1*k_1)+lt);
	/* Normalized Loop */
	#pragma cetus private(ax) 
	#pragma cetus lastprivate(k_2) 
	#pragma loop name setup#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(15L*lt)))) private(ax) lastprivate(k_2)
	for (k_2=0; k_2<=(-1+lt); k_2 ++ )
	{
		#pragma cetus private(ax) 
		#pragma loop name setup#2#0 
		for (ax=0; ax<3; ax ++ )
		{
			mi[(-1*k_2)+lt][ax]=(2+ng[(-1*k_2)+lt][ax]);
		}
		m1[(-1*k_2)+lt]=mi[(-1*k_2)+lt][0];
		m2[(-1*k_2)+lt]=mi[(-1*k_2)+lt][1];
		m3[(-1*k_2)+lt]=mi[(-1*k_2)+lt][2];
	}
	k=((-1*k_2)+lt);
	k=lt;
	is1=((2+ng[k][0])-ng[lt][0]);
	ie1=(1+ng[k][0]);
	( * n1)=((3+ie1)-is1);
	is2=((2+ng[k][1])-ng[lt][1]);
	ie2=(1+ng[k][1]);
	( * n2)=((3+ie2)-is2);
	is3=((2+ng[k][2])-ng[lt][2]);
	ie3=(1+ng[k][2]);
	( * n3)=((3+ie3)-is3);
	ir[lt]=0;
	/* Normalized Loop */
	#pragma cetus lastprivate(j_0) 
	#pragma loop name setup#3 
	for (j_0=0; j_0<=(-2+lt); j_0 ++ )
	{
		ir[(-1+(-1*j_0))+lt]=(ir[((-1+(-1*j_0))+lt)+1]+(((1*m1[((-1+(-1*j_0))+lt)+1])*m2[((-1+(-1*j_0))+lt)+1])*m3[((-1+(-1*j_0))+lt)+1]));
	}
	j=((-1+(-1*j_0))+lt);
	if (debug_vec[1]>=1)
	{
		printf(" in setup, \n");
		printf(" k  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
		printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n", k, lt, ng[k][0], ng[k][1], ng[k][2],  * n1,  * n2,  * n3, is1, is2, is3, ie1, ie2, ie3);
	}
}

/* --------------------------------------------------------------------- */
/* multigrid V-cycle routine */
/* --------------------------------------------------------------------- */
static void mg3P(double u[], double v[], double r[], double a[4], double c[4], int n1, int n2, int n3)
{
	int j, k;
	/* --------------------------------------------------------------------- */
	/* down cycle. */
	/* restrict the residual from the find grid to the coarse */
	/* --------------------------------------------------------------------- */
	int k_0;
	int k_1;
	/* Normalized Loop */
	#pragma cetus private(j) 
	#pragma cetus lastprivate(k_0) 
	#pragma loop name mg3P#0 
	for (k_0=0; k_0<=((-1+(-1*lb))+lt); k_0 ++ )
	{
		j=(((-1*k_0)+lt)-1);
		rprj3( & r[ir[(-1*k_0)+lt]], m1[(-1*k_0)+lt], m2[(-1*k_0)+lt], m3[(-1*k_0)+lt],  & r[ir[j]], m1[j], m2[j], m3[j], (-1*k_0)+lt);
	}
	k=((-1*k_0)+lt);
	k=lb;
	/* --------------------------------------------------------------------- */
	/* compute an approximate solution on the coarsest grid */
	/* --------------------------------------------------------------------- */
	zero3( & u[ir[k]], m1[k], m2[k], m3[k]);
	psinv( & r[ir[k]],  & u[ir[k]], m1[k], m2[k], m3[k], c, k);
	/* Normalized Loop */
	#pragma cetus private(j) 
	#pragma cetus lastprivate(k_1) 
	#pragma loop name mg3P#1 
	for (k_1=0; k_1<=((-2+(-1*lb))+lt); k_1 ++ )
	{
		j=(((1+k_1)+lb)-1);
		/* --------------------------------------------------------------------- */
		/* prolongate from level k-1  to k */
		/* --------------------------------------------------------------------- */
		zero3( & u[ir[(1+k_1)+lb]], m1[(1+k_1)+lb], m2[(1+k_1)+lb], m3[(1+k_1)+lb]);
		interp( & u[ir[j]], m1[j], m2[j], m3[j],  & u[ir[(1+k_1)+lb]], m1[(1+k_1)+lb], m2[(1+k_1)+lb], m3[(1+k_1)+lb], (1+k_1)+lb);
		/* --------------------------------------------------------------------- */
		/* compute residual for level k */
		/* --------------------------------------------------------------------- */
		resid( & u[ir[(1+k_1)+lb]],  & r[ir[(1+k_1)+lb]],  & r[ir[(1+k_1)+lb]], m1[(1+k_1)+lb], m2[(1+k_1)+lb], m3[(1+k_1)+lb], a, (1+k_1)+lb);
		/* --------------------------------------------------------------------- */
		/* apply smoother */
		/* --------------------------------------------------------------------- */
		psinv( & r[ir[(1+k_1)+lb]],  & u[ir[(1+k_1)+lb]], m1[(1+k_1)+lb], m2[(1+k_1)+lb], m3[(1+k_1)+lb], c, (1+k_1)+lb);
	}
	k=((1+k_1)+lb);
	j=(lt-1);
	k=lt;
	interp( & u[ir[j]], m1[j], m2[j], m3[j], u, n1, n2, n3, k);
	resid(u, v, r, n1, n2, n3, a, k);
	psinv(r, u, n1, n2, n3, c, k);
}

/* --------------------------------------------------------------------- */
/* psinv applies an approximate inverse as smoother:  u = u + Cr */
/*  */
/* This  implementation costs  15A + 4M per result, where */
/* A and M denote the costs of Addition and Multiplication.   */
/* Presuming coefficient c(3) is zero (the NPB assumes this, */
/* but it is thus not a general case), 2A + 1M may be eliminated, */
/* resulting in 13A + 3M. */
/* Note that this vectorizes, and is also fine for cache  */
/* based machines.   */
/* --------------------------------------------------------------------- */
static void psinv(void * or, void * ou, int n1, int n2, int n3, double c[4], int k)
{
	double (* r)[n2][n1] = (double (* )[n2][n1])or;
	double (* u)[n2][n1] = (double (* )[n2][n1])ou;
	int i3, i2, i1;
	double r1[((2+(1<<8))+1)], r2[((2+(1<<8))+1)];
	int i3_0;
	int i2_0;
	int i1_0;
	if (timeron)
	{
		timer_start(3);
	}
	#pragma cetus parallel 
	#pragma cetus private(i1, i1_0, i2, i2_0, r1, r2) 
	#pragma omp parallel if((10000<(((((((89L+(34L*n1))+(2L*n2))+(12L*n3))+((-20L*n1)*n2))+((-20L*n1)*n3))+((-1L*n2)*n3))+(((13L*n1)*n2)*n3)))) private(i1, i1_0, i2, i2_0, r1, r2)
	{
		double (* reduce0)[(-1+n2)][(-1+n1)] = (double (* )[(-1+n2)][(-1+n1)])malloc((((((((-1+n1)+n2)+n3)+((-1*n1)*n2))+((-1*n1)*n3))+((-1*n2)*n3))+((n1*n2)*n3))*sizeof (double));
		int reduce_span_0;
		int reduce_span_1;
		int reduce_span_2;
		for (reduce_span_0=0; reduce_span_0<(-1+n3); reduce_span_0 ++ )
		{
			for (reduce_span_1=0; reduce_span_1<(-1+n2); reduce_span_1 ++ )
			{
				for (reduce_span_2=0; reduce_span_2<(-1+n1); reduce_span_2 ++ )
				{
					reduce0[reduce_span_0][reduce_span_1][reduce_span_2]=0;
				}
			}
		}
		/* Normalized Loop */
		#pragma cetus lastprivate(i3_0) 
		#pragma loop name psinv#0 
		#pragma cetus for  
		#pragma omp for lastprivate(i3_0)
		for (i3_0=0; i3_0<=(-3+n3); i3_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(i1, i1_0) 
			#pragma cetus lastprivate(i2_0, r1, r2) 
			#pragma loop name psinv#0#0 
			/* #pragma cetus reduction(+: u[(1+i3_0)][(1+i2_0)][(1+i1_0)])  */
			for (i2_0=0; i2_0<=(-3+n2); i2_0 ++ )
			{
				#pragma cetus private(i1) 
				#pragma loop name psinv#0#0#0 
				for (i1=0; i1<n1; i1 ++ )
				{
					r1[i1]=(((r[1+i3_0][(1+i2_0)-1][i1]+r[1+i3_0][(1+i2_0)+1][i1])+r[(1+i3_0)-1][1+i2_0][i1])+r[(1+i3_0)+1][1+i2_0][i1]);
					r2[i1]=(((r[(1+i3_0)-1][(1+i2_0)-1][i1]+r[(1+i3_0)-1][(1+i2_0)+1][i1])+r[(1+i3_0)+1][(1+i2_0)-1][i1])+r[(1+i3_0)+1][(1+i2_0)+1][i1]);
				}
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_0) 
				#pragma loop name psinv#0#0#1 
				/* #pragma cetus reduction(+: u[(1+i3_0)][(1+i2_0)][(1+i1_0)])  */
				for (i1_0=0; i1_0<=(-3+n1); i1_0 ++ )
				{
					reduce0[1+i3_0][1+i2_0][1+i1_0]=(((reduce0[1+i3_0][1+i2_0][1+i1_0]+(c[0]*r[1+i3_0][1+i2_0][1+i1_0]))+(c[1]*((r[1+i3_0][1+i2_0][(1+i1_0)-1]+r[1+i3_0][1+i2_0][(1+i1_0)+1])+r1[1+i1_0])))+(c[2]*((r2[1+i1_0]+r1[(1+i1_0)-1])+r1[(1+i1_0)+1])));
					/* -------------------------------------------------------------------- */
					/* Assume c[3] = 0    (Enable line below if c[3] not= 0) */
					/* -------------------------------------------------------------------- */
					/*            + c[3] ( r2[i1-1] + r2[i1+1] ) */
					/* -------------------------------------------------------------------- */
				}
				i1=(1+i1_0);
			}
			i2=(1+i2_0);
		}
		#pragma cetus critical  
		#pragma omp critical
		{
			for (reduce_span_0=0; reduce_span_0<(-1+n3); reduce_span_0 ++ )
			{
				for (reduce_span_1=0; reduce_span_1<(-1+n2); reduce_span_1 ++ )
				{
					for (reduce_span_2=0; reduce_span_2<(-1+n1); reduce_span_2 ++ )
					{
						u[reduce_span_0][reduce_span_1][reduce_span_2]+=reduce0[reduce_span_0][reduce_span_1][reduce_span_2];
					}
				}
			}
		}
	}
	i3=(1+i3_0);
	if (timeron)
	{
		timer_stop(3);
	}
	/* --------------------------------------------------------------------- */
	/* exchange boundary points */
	/* --------------------------------------------------------------------- */
	comm3(u, n1, n2, n3, k);
	if (debug_vec[0]>=1)
	{
		rep_nrm(u, n1, n2, n3, "   psinv", k);
	}
	if (debug_vec[3]>=k)
	{
		showall(u, n1, n2, n3);
	}
}

/* --------------------------------------------------------------------- */
/* resid computes the residual:  r = v - Au */
/*  */
/* This  implementation costs  15A + 4M per result, where */
/* A and M denote the costs of Addition (or Subtraction) and  */
/* Multiplication, respectively.  */
/* Presuming coefficient a(1) is zero (the NPB assumes this, */
/* but it is thus not a general case), 3A + 1M may be eliminated, */
/* resulting in 12A + 3M. */
/* Note that this vectorizes, and is also fine for cache  */
/* based machines.   */
/* --------------------------------------------------------------------- */
static void resid(void * ou, void * ov, void * or, int n1, int n2, int n3, double a[4], int k)
{
	double (* u)[n2][n1] = (double (* )[n2][n1])ou;
	double (* v)[n2][n1] = (double (* )[n2][n1])ov;
	double (* r)[n2][n1] = (double (* )[n2][n1])or;
	int i3, i2, i1;
	double u1[((2+(1<<8))+1)], u2[((2+(1<<8))+1)];
	int i3_0;
	int i2_0;
	int i1_0;
	if (timeron)
	{
		timer_start(4);
	}
	/* Normalized Loop */
	#pragma cetus private(i1, i1_0, i2, i2_0, u1, u2) 
	#pragma cetus lastprivate(i3_0) 
	#pragma loop name resid#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((((((-11L+(28L*n1))+(2L*n2))+(6L*n3))+((-14L*n1)*n2))+((-14L*n1)*n3))+((-1L*n2)*n3))+(((7L*n1)*n2)*n3)))) private(i1, i1_0, i2, i2_0, u1, u2) lastprivate(i3_0)
	for (i3_0=0; i3_0<=(-3+n3); i3_0 ++ )
	{
		/* Normalized Loop */
		#pragma cetus private(i1, i1_0) 
		#pragma cetus lastprivate(i2_0, u1, u2) 
		#pragma loop name resid#0#0 
		for (i2_0=0; i2_0<=(-3+n2); i2_0 ++ )
		{
			#pragma cetus private(i1) 
			#pragma loop name resid#0#0#0 
			for (i1=0; i1<n1; i1 ++ )
			{
				u1[i1]=(((u[1+i3_0][(1+i2_0)-1][i1]+u[1+i3_0][(1+i2_0)+1][i1])+u[(1+i3_0)-1][1+i2_0][i1])+u[(1+i3_0)+1][1+i2_0][i1]);
				u2[i1]=(((u[(1+i3_0)-1][(1+i2_0)-1][i1]+u[(1+i3_0)-1][(1+i2_0)+1][i1])+u[(1+i3_0)+1][(1+i2_0)-1][i1])+u[(1+i3_0)+1][(1+i2_0)+1][i1]);
			}
			/* Normalized Loop */
			#pragma cetus lastprivate(i1_0) 
			#pragma loop name resid#0#0#1 
			for (i1_0=0; i1_0<=(-3+n1); i1_0 ++ )
			{
				/* ------------------------------------------------------------------- */
				/*  Assume a[1] = 0      (Enable 2 lines below if a[1] not= 0) */
				/* ------------------------------------------------------------------- */
				/*            - a[1] ( u[i3][i2][i1-1] + u[i3][i2][i1+1] */
				/*                     + u1[i1] ) */
				/* ------------------------------------------------------------------- */
				r[1+i3_0][1+i2_0][1+i1_0]=(((v[1+i3_0][1+i2_0][1+i1_0]-(a[0]*u[1+i3_0][1+i2_0][1+i1_0]))-(a[2]*((u2[1+i1_0]+u1[(1+i1_0)-1])+u1[(1+i1_0)+1])))-(a[3]*(u2[(1+i1_0)-1]+u2[(1+i1_0)+1])));
			}
			i1=(1+i1_0);
		}
		i2=(1+i2_0);
	}
	i3=(1+i3_0);
	if (timeron)
	{
		timer_stop(4);
	}
	/* --------------------------------------------------------------------- */
	/* exchange boundary data */
	/* --------------------------------------------------------------------- */
	comm3(r, n1, n2, n3, k);
	if (debug_vec[0]>=1)
	{
		rep_nrm(r, n1, n2, n3, "   resid", k);
	}
	if (debug_vec[2]>=k)
	{
		showall(r, n1, n2, n3);
	}
}

/* --------------------------------------------------------------------- */
/* rprj3 projects onto the next coarser grid,  */
/* using a trilinear Finite Element projection:  s = r' = P r */
/*      */
/* This  implementation costs  20A + 4M per result, where */
/* A and M denote the costs of Addition and Multiplication.   */
/* Note that this vectorizes, and is also fine for cache  */
/* based machines.   */
/* --------------------------------------------------------------------- */
static void rprj3(void * or, int m1k, int m2k, int m3k, void * os, int m1j, int m2j, int m3j, int k)
{
	double (* r)[m2k][m1k] = (double (* )[m2k][m1k])or;
	double (* s)[m2j][m1j] = (double (* )[m2j][m1j])os;
	int j3, j2, j1, i3, i2, i1, d1, d2, d3, j;
	double x1[((2+(1<<8))+1)], y1[((2+(1<<8))+1)], x2, y2;
	int j3_0;
	int j2_0;
	int j1_0;
	int j1_1;
	if (timeron)
	{
		timer_start(6);
	}
	if (m1k==3)
	{
		d1=2;
	}
	else
	{
		d1=1;
	}
	if (m2k==3)
	{
		d2=2;
	}
	else
	{
		d2=1;
	}
	if (m3k==3)
	{
		d3=2;
	}
	else
	{
		d3=1;
	}
	/* Normalized Loop */
	#pragma cetus private(i1, i2, i3, j1, j1_0, j1_1, j2, j2_0, x2, y2) 
	#pragma cetus lastprivate(j3_0) 
	#pragma loop name rprj3#0 
	for (j3_0=0; j3_0<=(-3+m3j); j3_0 ++ )
	{
		i3=((2*(1+j3_0))-d3);
		/* Normalized Loop */
		#pragma cetus private(i1, i2, j1, j1_0, j1_1, x2, y2) 
		#pragma cetus lastprivate(j2_0) 
		#pragma loop name rprj3#0#0 
		for (j2_0=0; j2_0<=(-3+m2j); j2_0 ++ )
		{
			i2=((2*(1+j2_0))-d2);
			/* Normalized Loop */
			#pragma cetus private(i1) 
			#pragma cetus lastprivate(j1_0) 
			#pragma loop name rprj3#0#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-4L+(5L*m1j)))) private(i1) lastprivate(j1_0)
			for (j1_0=0; j1_0<=(-2+m1j); j1_0 ++ )
			{
				i1=((2*(1+j1_0))-d1);
				x1[i1]=(((r[i3+1][i2][i1]+r[i3+1][i2+2][i1])+r[i3][i2+1][i1])+r[i3+2][i2+1][i1]);
				y1[i1]=(((r[i3][i2][i1]+r[i3+2][i2][i1])+r[i3][i2+2][i1])+r[i3+2][i2+2][i1]);
			}
			j1=(1+j1_0);
			/* Normalized Loop */
			#pragma cetus private(i1, x2, y2) 
			#pragma cetus lastprivate(j1_1) 
			#pragma loop name rprj3#0#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-11L+(6L*m1j)))) private(i1, x2, y2) lastprivate(j1_1)
			for (j1_1=0; j1_1<=(-3+m1j); j1_1 ++ )
			{
				i1=((2*(1+j1_1))-d1);
				y2=(((r[i3][i2][i1+1]+r[i3+2][i2][i1+1])+r[i3][i2+2][i1+1])+r[i3+2][i2+2][i1+1]);
				x2=(((r[i3+1][i2][i1+1]+r[i3+1][i2+2][i1+1])+r[i3][i2+1][i1+1])+r[i3+2][i2+1][i1+1]);
				s[1+j3_0][1+j2_0][1+j1_1]=((((0.5*r[i3+1][i2+1][i1+1])+(0.25*((r[i3+1][i2+1][i1]+r[i3+1][i2+1][i1+2])+x2)))+(0.125*((x1[i1]+x1[i1+2])+y2)))+(0.0625*(y1[i1]+y1[i1+2])));
			}
			j1=(1+j1_1);
		}
		j2=(1+j2_0);
	}
	j3=(1+j3_0);
	if (timeron)
	{
		timer_stop(6);
	}
	j=(k-1);
	comm3(s, m1j, m2j, m3j, j);
	if (debug_vec[0]>=1)
	{
		rep_nrm(s, m1j, m2j, m3j, "   rprj3", k-1);
	}
	if (debug_vec[4]>=k)
	{
		showall(s, m1j, m2j, m3j);
	}
}

/* --------------------------------------------------------------------- */
/* interp adds the trilinear interpolation of the correction */
/* from the coarser grid to the current approximation:  u = u + Qu' */
/*      */
/* Observe that this  implementation costs  16A + 4M, where */
/* A and M denote the costs of Addition and Multiplication.   */
/* Note that this vectorizes, and is also fine for cache  */
/* based machines.  Vector machines may get slightly better  */
/* performance however, with 8 separate "do i1" loops, rather than 4. */
/* --------------------------------------------------------------------- */
static void interp(void * oz, int mm1, int mm2, int mm3, void * ou, int n1, int n2, int n3, int k)
{
	double (* z)[mm2][mm1] = (double (* )[mm2][mm1])oz;
	double (* u)[n2][n1] = (double (* )[n2][n1])ou;
	int i3, i2, i1, d1, d2, d3, t1, t2, t3;
	/* note that m = 1037 in globals.h but for this only need to be */
	/* 535 to handle up to 1024^3 */
	/*      integer m */
	/*      parameter( m=535 ) */
	double z1[((2+(1<<8))+1)], z2[((2+(1<<8))+1)], z3[((2+(1<<8))+1)];
	int i3_0;
	int i2_0;
	int i1_0;
	int i1_1;
	int i2_1;
	int i1_2;
	int i1_3;
	int i3_1;
	int i2_2;
	int i1_4;
	int i1_5;
	int i2_3;
	int i1_6;
	int i1_7;
	if (timeron)
	{
		timer_start(7);
	}
	if (((n1!=3)&&(n2!=3))&&(n3!=3))
	{
		#pragma cetus private(i1, i2, i3, z1, z2, z3) 
		#pragma loop name interp#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((((((-11L+(21L*mm1))+(9L*mm2))+(12L*mm3))+((-21L*mm1)*mm2))+((-21L*mm1)*mm3))+((-9L*mm2)*mm3))+(((21L*mm1)*mm2)*mm3)))) private(i1, i2, i3, z1, z2, z3)
		for (i3=0; i3<(mm3-1); i3 ++ )
		{
			#pragma cetus private(i1, i2) 
			#pragma cetus lastprivate(z1, z2, z3) 
			#pragma loop name interp#0#0 
			for (i2=0; i2<(mm2-1); i2 ++ )
			{
				#pragma cetus private(i1) 
				#pragma loop name interp#0#0#0 
				for (i1=0; i1<mm1; i1 ++ )
				{
					z1[i1]=(z[i3][i2+1][i1]+z[i3][i2][i1]);
					z2[i1]=(z[i3+1][i2][i1]+z[i3][i2][i1]);
					z3[i1]=((z[i3+1][i2+1][i1]+z[i3+1][i2][i1])+z1[i1]);
				}
				#pragma cetus private(i1) 
				#pragma loop name interp#0#0#1 
				for (i1=0; i1<(mm1-1); i1 ++ )
				{
					u[2*i3][2*i2][2*i1]=(u[2*i3][2*i2][2*i1]+z[i3][i2][i1]);
					u[2*i3][2*i2][(2*i1)+1]=(u[2*i3][2*i2][(2*i1)+1]+(0.5*(z[i3][i2][i1+1]+z[i3][i2][i1])));
				}
				#pragma cetus private(i1) 
				#pragma loop name interp#0#0#2 
				for (i1=0; i1<(mm1-1); i1 ++ )
				{
					u[2*i3][(2*i2)+1][2*i1]=(u[2*i3][(2*i2)+1][2*i1]+(0.5*z1[i1]));
					u[2*i3][(2*i2)+1][(2*i1)+1]=(u[2*i3][(2*i2)+1][(2*i1)+1]+(0.25*(z1[i1]+z1[i1+1])));
				}
				#pragma cetus private(i1) 
				#pragma loop name interp#0#0#3 
				for (i1=0; i1<(mm1-1); i1 ++ )
				{
					u[(2*i3)+1][2*i2][2*i1]=(u[(2*i3)+1][2*i2][2*i1]+(0.5*z2[i1]));
					u[(2*i3)+1][2*i2][(2*i1)+1]=(u[(2*i3)+1][2*i2][(2*i1)+1]+(0.25*(z2[i1]+z2[i1+1])));
				}
				#pragma cetus private(i1) 
				#pragma loop name interp#0#0#4 
				for (i1=0; i1<(mm1-1); i1 ++ )
				{
					u[(2*i3)+1][(2*i2)+1][2*i1]=(u[(2*i3)+1][(2*i2)+1][2*i1]+(0.25*z3[i1]));
					u[(2*i3)+1][(2*i2)+1][(2*i1)+1]=(u[(2*i3)+1][(2*i2)+1][(2*i1)+1]+(0.125*(z3[i1]+z3[i1+1])));
				}
			}
		}
	}
	else
	{
		if (n1==3)
		{
			d1=2;
			t1=1;
		}
		else
		{
			d1=1;
			t1=0;
		}
		if (n2==3)
		{
			d2=2;
			t2=1;
		}
		else
		{
			d2=1;
			t2=0;
		}
		if (n3==3)
		{
			d3=2;
			t3=1;
		}
		else
		{
			d3=1;
			t3=0;
		}
		/* Normalized Loop */
		#pragma cetus private(i1, i1_0, i1_1, i1_2, i1_3, i2, i2_0, i2_1) 
		#pragma cetus lastprivate(i3_0) 
		#pragma loop name interp#1 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((((((((((((((((((1L+(-3L*d3))+(3L*mm3))+((-3L*d1)*d3))+((3L*d1)*mm3))+((3L*d2)*d3))+((-3L*d2)*mm3))+((6L*d3)*mm1))+((-6L*d3)*mm2))+((-6L*mm1)*mm3))+((6L*mm2)*mm3))+(((-3L*d1)*d2)*d3))+(((3L*d1)*d2)*mm3))+(((6L*d1)*d3)*mm2))+(((-6L*d1)*mm2)*mm3))+(((6L*d2)*d3)*mm1))+(((-6L*d2)*mm1)*mm3))+(((-12L*d3)*mm1)*mm2))+(((12L*mm1)*mm2)*mm3)))) private(i1, i1_0, i1_1, i1_2, i1_3, i2, i2_0, i2_1) lastprivate(i3_0)
		for (i3_0=0; i3_0<=((-1+(-1*d3))+mm3); i3_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(i1, i1_0, i1_1) 
			#pragma cetus lastprivate(i2_0) 
			#pragma loop name interp#1#0 
			for (i2_0=0; i2_0<=((-1+(-1*d2))+mm2); i2_0 ++ )
			{
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_0) 
				#pragma loop name interp#1#0#0 
				for (i1_0=0; i1_0<=((-1+(-1*d1))+mm1); i1_0 ++ )
				{
					u[((2*(d3+i3_0))-d3)-1][((2*(d2+i2_0))-d2)-1][((2*(d1+i1_0))-d1)-1]=(u[((2*(d3+i3_0))-d3)-1][((2*(d2+i2_0))-d2)-1][((2*(d1+i1_0))-d1)-1]+z[(d3+i3_0)-1][(d2+i2_0)-1][(d1+i1_0)-1]);
				}
				i1=(d1+i1_0);
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_1) 
				#pragma loop name interp#1#0#1 
				for (i1_1=0; i1_1<=(-2+mm1); i1_1 ++ )
				{
					u[((2*(d3+i3_0))-d3)-1][((2*(d2+i2_0))-d2)-1][((2*(1+i1_1))-t1)-1]=(u[((2*(d3+i3_0))-d3)-1][((2*(d2+i2_0))-d2)-1][((2*(1+i1_1))-t1)-1]+(0.5*(z[(d3+i3_0)-1][(d2+i2_0)-1][1+i1_1]+z[(d3+i3_0)-1][(d2+i2_0)-1][(1+i1_1)-1])));
				}
				i1=(1+i1_1);
			}
			i2=(d2+i2_0);
			/* Normalized Loop */
			#pragma cetus private(i1, i1_2, i1_3) 
			#pragma cetus lastprivate(i2_1) 
			#pragma loop name interp#1#1 
			for (i2_1=0; i2_1<=(-2+mm2); i2_1 ++ )
			{
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_2) 
				#pragma loop name interp#1#1#0 
				for (i1_2=0; i1_2<=((-1+(-1*d1))+mm1); i1_2 ++ )
				{
					u[((2*(d3+i3_0))-d3)-1][((2*(1+i2_1))-t2)-1][((2*(d1+i1_2))-d1)-1]=(u[((2*(d3+i3_0))-d3)-1][((2*(1+i2_1))-t2)-1][((2*(d1+i1_2))-d1)-1]+(0.5*(z[(d3+i3_0)-1][1+i2_1][(d1+i1_2)-1]+z[(d3+i3_0)-1][(1+i2_1)-1][(d1+i1_2)-1])));
				}
				i1=(d1+i1_2);
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_3) 
				#pragma loop name interp#1#1#1 
				for (i1_3=0; i1_3<=(-2+mm1); i1_3 ++ )
				{
					u[((2*(d3+i3_0))-d3)-1][((2*(1+i2_1))-t2)-1][((2*(1+i1_3))-t1)-1]=(u[((2*(d3+i3_0))-d3)-1][((2*(1+i2_1))-t2)-1][((2*(1+i1_3))-t1)-1]+(0.25*(((z[(d3+i3_0)-1][1+i2_1][1+i1_3]+z[(d3+i3_0)-1][(1+i2_1)-1][1+i1_3])+z[(d3+i3_0)-1][1+i2_1][(1+i1_3)-1])+z[(d3+i3_0)-1][(1+i2_1)-1][(1+i1_3)-1])));
				}
				i1=(1+i1_3);
			}
			i2=(1+i2_1);
		}
		i3=(d3+i3_0);
		/* Normalized Loop */
		#pragma cetus private(i1, i1_4, i1_5, i1_6, i1_7, i2, i2_2, i2_3) 
		#pragma cetus lastprivate(i3_1) 
		#pragma loop name interp#2 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((((((((((((((((-2L+(-3L*d1))+(3L*d2))+(6L*mm1))+(-6L*mm2))+(3L*mm3))+((-3L*d1)*d2))+((6L*d1)*mm2))+((3L*d1)*mm3))+((6L*d2)*mm1))+((-3L*d2)*mm3))+((-12L*mm1)*mm2))+((-6L*mm1)*mm3))+((6L*mm2)*mm3))+(((3L*d1)*d2)*mm3))+(((-6L*d1)*mm2)*mm3))+(((-6L*d2)*mm1)*mm3))+(((12L*mm1)*mm2)*mm3)))) private(i1, i1_4, i1_5, i1_6, i1_7, i2, i2_2, i2_3) lastprivate(i3_1)
		for (i3_1=0; i3_1<=(-2+mm3); i3_1 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(i1, i1_4, i1_5) 
			#pragma cetus lastprivate(i2_2) 
			#pragma loop name interp#2#0 
			for (i2_2=0; i2_2<=((-1+(-1*d2))+mm2); i2_2 ++ )
			{
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_4) 
				#pragma loop name interp#2#0#0 
				for (i1_4=0; i1_4<=((-1+(-1*d1))+mm1); i1_4 ++ )
				{
					u[((2*(1+i3_1))-t3)-1][((2*(d2+i2_2))-d2)-1][((2*(d1+i1_4))-d1)-1]=(u[((2*(1+i3_1))-t3)-1][((2*(d2+i2_2))-d2)-1][((2*(d1+i1_4))-d1)-1]+(0.5*(z[1+i3_1][(d2+i2_2)-1][(d1+i1_4)-1]+z[(1+i3_1)-1][(d2+i2_2)-1][(d1+i1_4)-1])));
				}
				i1=(d1+i1_4);
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_5) 
				#pragma loop name interp#2#0#1 
				for (i1_5=0; i1_5<=(-2+mm1); i1_5 ++ )
				{
					u[((2*(1+i3_1))-t3)-1][((2*(d2+i2_2))-d2)-1][((2*(1+i1_5))-t1)-1]=(u[((2*(1+i3_1))-t3)-1][((2*(d2+i2_2))-d2)-1][((2*(1+i1_5))-t1)-1]+(0.25*(((z[1+i3_1][(d2+i2_2)-1][1+i1_5]+z[1+i3_1][(d2+i2_2)-1][(1+i1_5)-1])+z[(1+i3_1)-1][(d2+i2_2)-1][1+i1_5])+z[(1+i3_1)-1][(d2+i2_2)-1][(1+i1_5)-1])));
				}
				i1=(1+i1_5);
			}
			i2=(d2+i2_2);
			/* Normalized Loop */
			#pragma cetus private(i1, i1_6, i1_7) 
			#pragma cetus lastprivate(i2_3) 
			#pragma loop name interp#2#1 
			for (i2_3=0; i2_3<=(-2+mm2); i2_3 ++ )
			{
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_6) 
				#pragma loop name interp#2#1#0 
				for (i1_6=0; i1_6<=((-1+(-1*d1))+mm1); i1_6 ++ )
				{
					u[((2*(1+i3_1))-t3)-1][((2*(1+i2_3))-t2)-1][((2*(d1+i1_6))-d1)-1]=(u[((2*(1+i3_1))-t3)-1][((2*(1+i2_3))-t2)-1][((2*(d1+i1_6))-d1)-1]+(0.25*(((z[1+i3_1][1+i2_3][(d1+i1_6)-1]+z[1+i3_1][(1+i2_3)-1][(d1+i1_6)-1])+z[(1+i3_1)-1][1+i2_3][(d1+i1_6)-1])+z[(1+i3_1)-1][(1+i2_3)-1][(d1+i1_6)-1])));
				}
				i1=(d1+i1_6);
				/* Normalized Loop */
				#pragma cetus lastprivate(i1_7) 
				#pragma loop name interp#2#1#1 
				for (i1_7=0; i1_7<=(-2+mm1); i1_7 ++ )
				{
					u[((2*(1+i3_1))-t3)-1][((2*(1+i2_3))-t2)-1][((2*(1+i1_7))-t1)-1]=(u[((2*(1+i3_1))-t3)-1][((2*(1+i2_3))-t2)-1][((2*(1+i1_7))-t1)-1]+(0.125*(((((((z[1+i3_1][1+i2_3][1+i1_7]+z[1+i3_1][(1+i2_3)-1][1+i1_7])+z[1+i3_1][1+i2_3][(1+i1_7)-1])+z[1+i3_1][(1+i2_3)-1][(1+i1_7)-1])+z[(1+i3_1)-1][1+i2_3][1+i1_7])+z[(1+i3_1)-1][(1+i2_3)-1][1+i1_7])+z[(1+i3_1)-1][1+i2_3][(1+i1_7)-1])+z[(1+i3_1)-1][(1+i2_3)-1][(1+i1_7)-1])));
				}
				i1=(1+i1_7);
			}
			i2=(1+i2_3);
		}
		i3=(1+i3_1);
	}
	if (timeron)
	{
		timer_stop(7);
	}
	if (debug_vec[0]>=1)
	{
		rep_nrm(z, mm1, mm2, mm3, "z: inter", k-1);
		rep_nrm(u, n1, n2, n3, "u: inter", k);
	}
	if (debug_vec[5]>=k)
	{
		showall(z, mm1, mm2, mm3);
		showall(u, n1, n2, n3);
	}
}

/* --------------------------------------------------------------------- */
/* norm2u3 evaluates approximations to the L2 norm and the */
/* uniform (or L-infinity or Chebyshev) norm, under the */
/* assumption that the boundaries are periodic or zero.  Add the */
/* boundaries in with half weight (quarter weight on the edges */
/* and eighth weight at the corners) for inhomogeneous boundaries. */
/* --------------------------------------------------------------------- */
static void norm2u3(void * or, int n1, int n2, int n3, double * rnm2, double * rnmu, int nx, int ny, int nz)
{
	double (* r)[n2][n1] = (double (* )[n2][n1])or;
	double s, a;
	int i3, i2, i1;
	double dn;
	int i3_0;
	int i2_0;
	int i1_0;
	if (timeron)
	{
		timer_start(8);
	}
	dn=(((1.0*nx)*ny)*nz);
	s=0.0;
	( * rnmu)=0.0;
	/* Normalized Loop */
	#pragma cetus private(a, i1, i1_0, i2, i2_0) 
	#pragma cetus lastprivate(i3_0) 
	#pragma loop name norm2u3#0 
	for (i3_0=0; i3_0<=(-3+n3); i3_0 ++ )
	{
		/* Normalized Loop */
		#pragma cetus private(a, i1, i1_0) 
		#pragma cetus lastprivate(i2_0) 
		#pragma loop name norm2u3#0#0 
		for (i2_0=0; i2_0<=(-3+n2); i2_0 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(a) 
			#pragma cetus lastprivate(i1_0) 
			#pragma loop name norm2u3#0#0#0 
			for (i1_0=0; i1_0<=(-3+n1); i1_0 ++ )
			{
				s=(s+pow(r[1+i3_0][1+i2_0][1+i1_0], 2.0));
				a=fabs(r[1+i3_0][1+i2_0][1+i1_0]);
				if (a>( * rnmu))
				{
					( * rnmu)=a;
				}
			}
			i1=(1+i1_0);
		}
		i2=(1+i2_0);
	}
	i3=(1+i3_0);
	( * rnm2)=sqrt(s/dn);
	if (timeron)
	{
		timer_stop(8);
	}
}

/* --------------------------------------------------------------------- */
/* report on norm */
/* --------------------------------------------------------------------- */
static void rep_nrm(void * u, int n1, int n2, int n3, char * title, int kk)
{
	double rnm2, rnmu;
	norm2u3(u, n1, n2, n3,  & rnm2,  & rnmu, nx[kk], ny[kk], nz[kk]);
	printf(" Level%2d in %8s: norms =%21.14E%21.14E\n", kk, title, rnm2, rnmu);
}

/* --------------------------------------------------------------------- */
/* comm3 organizes the communication on all borders  */
/* --------------------------------------------------------------------- */
static void comm3(void * ou, int n1, int n2, int n3, int kk)
{
	double (* u)[n2][n1] = (double (* )[n2][n1])ou;
	int i1, i2, i3;
	int i3_0;
	int i2_0;
	int i3_1;
	if (timeron)
	{
		timer_start(9);
	}
	/* Normalized Loop */
	#pragma cetus private(i2, i2_0) 
	#pragma cetus lastprivate(i3_0) 
	#pragma loop name comm3#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((9L+(-8L*n2))+(-4L*n3))+((4L*n2)*n3)))) private(i2, i2_0) lastprivate(i3_0)
	for (i3_0=0; i3_0<=(-3+n3); i3_0 ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(i2_0) 
		#pragma loop name comm3#0#0 
		for (i2_0=0; i2_0<=(-3+n2); i2_0 ++ )
		{
			u[1+i3_0][1+i2_0][0]=u[1+i3_0][1+i2_0][n1-2];
			u[1+i3_0][1+i2_0][n1-1]=u[1+i3_0][1+i2_0][1];
		}
		i2=(1+i2_0);
	}
	i3=(1+i3_0);
	/* Normalized Loop */
	#pragma cetus private(i1) 
	#pragma cetus lastprivate(i3_1) 
	#pragma loop name comm3#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((-5L+(-8L*n1))+(3L*n3))+((4L*n1)*n3)))) private(i1) lastprivate(i3_1)
	for (i3_1=0; i3_1<=(-3+n3); i3_1 ++ )
	{
		#pragma cetus private(i1) 
		#pragma loop name comm3#1#0 
		for (i1=0; i1<n1; i1 ++ )
		{
			u[1+i3_1][0][i1]=u[1+i3_1][n2-2][i1];
			u[1+i3_1][n2-1][i1]=u[1+i3_1][1][i1];
		}
	}
	i3=(1+i3_1);
	#pragma cetus private(i1, i2) 
	#pragma loop name comm3#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*n2))+((4L*n1)*n2)))) private(i1, i2)
	for (i2=0; i2<n2; i2 ++ )
	{
		#pragma cetus private(i1) 
		#pragma loop name comm3#2#0 
		for (i1=0; i1<n1; i1 ++ )
		{
			u[0][i2][i1]=u[n3-2][i2][i1];
			u[n3-1][i2][i1]=u[1][i2][i1];
		}
	}
	if (timeron)
	{
		timer_stop(9);
	}
}

/* --------------------------------------------------------------------- */
/* zran3  loads +1 at ten randomly chosen points, */
/* loads -1 at a different ten random points, */
/* and zero elsewhere. */
/* --------------------------------------------------------------------- */
static void zran3(void * oz, int n1, int n2, int n3, int nx, int ny, int k)
{
	double (* z)[n2][n1] = (double (* )[n2][n1])oz;
	int i0, m0, m1;
	int i1, i2, i3, d1, e1, e2, e3;
	double xx, x0, x1, a1, a2, ai;
	const int mm = 10;
	const double a = pow(5.0, 13.0);
	const double x = 3.14159265E8;
	double ten[mm][2], best;
	int i, j1[mm][2], j2[mm][2], j3[mm][2];
	int jg[4][mm][2];
	double rdummy;
	int i3_0;
	int i2_0;
	int i3_1;
	int i2_1;
	int i1_0;
	int i_0;
	int i_1;
	int i_2;
	a1=power(a, nx);
	a2=power(a, nx*ny);
	zero3(z, n1, n2, n3);
	i=((is1-2)+(nx*((is2-2)+(ny*(is3-2)))));
	ai=power(a, i);
	d1=((ie1-is1)+1);
	e1=((ie1-is1)+2);
	e2=((ie2-is2)+2);
	e3=((ie3-is3)+2);
	x0=x;
	rdummy=randlc( & x0, ai);
	/* Normalized Loop */
	#pragma cetus private(i2, i2_0, rdummy) 
	#pragma cetus lastprivate(i3_0) 
	#pragma loop name zran3#0 
	for (i3_0=0; i3_0<=(-2+e3); i3_0 ++ )
	{
		x1=x0;
		/* Normalized Loop */
		#pragma cetus private(rdummy) 
		#pragma cetus lastprivate(i2_0) 
		#pragma loop name zran3#0#0 
		for (i2_0=0; i2_0<=(-2+e2); i2_0 ++ )
		{
			xx=x1;
			vranlc(d1,  & xx, a,  & z[1+i3_0][1+i2_0][1]);
			rdummy=randlc( & x1, a1);
		}
		i2=(1+i2_0);
		rdummy=randlc( & x0, a2);
	}
	i3=(1+i3_0);
	/* --------------------------------------------------------------------- */
	/* comm3(z,n1,n2,n3); */
	/* showall(z,n1,n2,n3); */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* each processor looks for twenty candidates */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i) 
	#pragma loop name zran3#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(10L*mm)))) private(i)
	for (i=0; i<mm; i ++ )
	{
		ten[i][1]=0.0;
		j1[i][1]=0;
		j2[i][1]=0;
		j3[i][1]=0;
		ten[i][0]=1.0;
		j1[i][0]=0;
		j2[i][0]=0;
		j3[i][0]=0;
	}
	/* Normalized Loop */
	#pragma cetus private(i1, i1_0, i2, i2_1) 
	#pragma cetus lastprivate(i3_1) 
	#pragma loop name zran3#2 
	for (i3_1=0; i3_1<=(-3+n3); i3_1 ++ )
	{
		/* Normalized Loop */
		#pragma cetus private(i1, i1_0) 
		#pragma cetus lastprivate(i2_1) 
		#pragma loop name zran3#2#0 
		for (i2_1=0; i2_1<=(-3+n2); i2_1 ++ )
		{
			/* Normalized Loop */
			#pragma cetus lastprivate(i1_0) 
			#pragma loop name zran3#2#0#0 
			for (i1_0=0; i1_0<=(-3+n1); i1_0 ++ )
			{
				if (z[1+i3_1][1+i2_1][1+i1_0]>ten[0][1])
				{
					ten[0][1]=z[1+i3_1][1+i2_1][1+i1_0];
					j1[0][1]=(1+i1_0);
					j2[0][1]=(1+i2_1);
					j3[0][1]=(1+i3_1);
					bubble(ten, j1, j2, j3, mm, 1);
				}
				if (z[1+i3_1][1+i2_1][1+i1_0]<ten[0][0])
				{
					ten[0][0]=z[1+i3_1][1+i2_1][1+i1_0];
					j1[0][0]=(1+i1_0);
					j2[0][0]=(1+i2_1);
					j3[0][0]=(1+i3_1);
					bubble(ten, j1, j2, j3, mm, 0);
				}
			}
			i1=(1+i1_0);
		}
		i2=(1+i2_1);
	}
	i3=(1+i3_1);
	/* --------------------------------------------------------------------- */
	/* Now which of these are globally best? */
	/* --------------------------------------------------------------------- */
	i1=(mm-1);
	i0=(mm-1);
	/* Normalized Loop */
	#pragma cetus private(best) 
	#pragma cetus lastprivate(i_0) 
	#pragma loop name zran3#3 
	for (i_0=0; i_0<=(-1+mm); i_0 ++ )
	{
		best=0.0;
		if (best<ten[i1][1])
		{
			jg[0][(-1+(-1*i_0))+mm][1]=0;
			jg[1][(-1+(-1*i_0))+mm][1]=((is1-2)+j1[i1][1]);
			jg[2][(-1+(-1*i_0))+mm][1]=((is2-2)+j2[i1][1]);
			jg[3][(-1+(-1*i_0))+mm][1]=((is3-2)+j3[i1][1]);
			i1=(i1-1);
		}
		else
		{
			jg[0][(-1+(-1*i_0))+mm][1]=0;
			jg[1][(-1+(-1*i_0))+mm][1]=0;
			jg[2][(-1+(-1*i_0))+mm][1]=0;
			jg[3][(-1+(-1*i_0))+mm][1]=0;
		}
		best=1.0;
		if (best>ten[i0][0])
		{
			jg[0][(-1+(-1*i_0))+mm][0]=0;
			jg[1][(-1+(-1*i_0))+mm][0]=((is1-2)+j1[i0][0]);
			jg[2][(-1+(-1*i_0))+mm][0]=((is2-2)+j2[i0][0]);
			jg[3][(-1+(-1*i_0))+mm][0]=((is3-2)+j3[i0][0]);
			i0=(i0-1);
		}
		else
		{
			jg[0][(-1+(-1*i_0))+mm][0]=0;
			jg[1][(-1+(-1*i_0))+mm][0]=0;
			jg[2][(-1+(-1*i_0))+mm][0]=0;
			jg[3][(-1+(-1*i_0))+mm][0]=0;
		}
	}
	i=((-1+(-1*i_0))+mm);
	/*  m1 = i1+1; */
	/*  m0 = i0+1; */
	m1=0;
	m0=0;
	/*
	
	  int cnt = 0;
	  printf("  \n");
	  printf("  negative charges at\n");
	  for (i = 0; i < mm; i++) {
		    printf(" (%3d,%3d,%3d)", jg[1][i][0], jg[2][i][0], jg[3][i][0]);
		    if (++cnt % 5 == 0) printf("\n");
	  }
	
	  cnt = 0;
	  printf("  positive charges at\n");
	  for (i = 0; i < mm; i++) {
		    printf(" (%3d,%3d,%3d)", jg[1][i][1], jg[2][i][1], jg[3][i][1]);
		    if (++cnt % 5 == 0) printf("\n");
	  }
	
	  cnt = 0;
	  printf("  small random numbers were\n");
	  for (i = mm-1; i >= 0; i--) {
		    printf(" %15.8E", ten[i][0]);
		    if (++cnt % 5 == 0) printf("\n");
	  }
	
	  cnt = 0;
	  printf("  and they were found on processor number\n");
	  for (i = mm-1; i >= 0; i--) {
		    printf(" %4d", jg[0][i][0]);
		    if (++cnt % 10 == 0) printf("\n");
	  }
	
	  cnt = 0;
	  printf("  large random numbers were\n");
	  for (i = mm-1; i >= 0; i--) {
		    printf(" %15.8E", ten[i][1]);
		    if (++cnt % 5 == 0) printf("\n");
	  }
	
	  cnt = 0;
	  printf("  and they were found on processor number\n");
	  for (i = mm-1; i >= 0; i--) {
		    printf(" %4d", jg[0][i][1]);
		    if (++cnt % 10 == 0) printf("\n");
	  }
	 
	*/
	#pragma cetus private(i1, i2, i3) 
	#pragma loop name zran3#4 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*n3))+((3L*n2)*n3))+(((3L*n1)*n2)*n3)))) private(i1, i2, i3)
	for (i3=0; i3<n3; i3 ++ )
	{
		#pragma cetus private(i1, i2) 
		#pragma loop name zran3#4#0 
		for (i2=0; i2<n2; i2 ++ )
		{
			#pragma cetus private(i1) 
			#pragma loop name zran3#4#0#0 
			for (i1=0; i1<n1; i1 ++ )
			{
				z[i3][i2][i1]=0.0;
			}
		}
	}
	/* Normalized Loop */
	#pragma cetus lastprivate(i_1) 
	#pragma loop name zran3#5 
	for (i_1=0; i_1<=((-1+(-1*m0))+mm); i_1 ++ )
	{
		z[jg[3][(-1+(-1*i_1))+mm][0]][jg[2][(-1+(-1*i_1))+mm][0]][jg[1][(-1+(-1*i_1))+mm][0]]=( - 1.0);
	}
	i=((-1+(-1*i_1))+mm);
	/* Normalized Loop */
	#pragma cetus lastprivate(i_2) 
	#pragma loop name zran3#6 
	for (i_2=0; i_2<=((-1+(-1*m1))+mm); i_2 ++ )
	{
		z[jg[3][(-1+(-1*i_2))+mm][1]][jg[2][(-1+(-1*i_2))+mm][1]][jg[1][(-1+(-1*i_2))+mm][1]]=( + 1.0);
	}
	i=((-1+(-1*i_2))+mm);
	comm3(z, n1, n2, n3, k);
	/* --------------------------------------------------------------------- */
	/* showall(z,n1,n2,n3); */
	/* --------------------------------------------------------------------- */
}

static void showall(void * oz, int n1, int n2, int n3)
{
	double (* z)[n2][n1] = (double (* )[n2][n1])oz;
	int i1, i2, i3;
	int m1, m2, m3;
	m1=((n1<18) ? n1 : 18);
	m2=((n2<14) ? n2 : 14);
	m3=((n3<18) ? n3 : 18);
	printf("   \n");
	#pragma cetus private(i1, i2, i3) 
	#pragma loop name showall#0 
	for (i3=0; i3<m3; i3 ++ )
	{
		#pragma cetus private(i1, i2) 
		#pragma loop name showall#0#0 
		for (i1=0; i1<m1; i1 ++ )
		{
			#pragma cetus private(i2) 
			#pragma loop name showall#0#0#0 
			for (i2=0; i2<m2; i2 ++ )
			{
				printf("%6.3f", z[i3][i2][i1]);
			}
			printf("\n");
		}
		printf("  - - - - - - - \n");
	}
	printf("   \n");
}

/* --------------------------------------------------------------------- */
/* power  raises an integer, disguised as a double */
/* precision real, to an integer power */
/* --------------------------------------------------------------------- */
static double power(double a, int n)
{
	double aj;
	int nj;
	double rdummy;
	double power;
	power=1.0;
	nj=n;
	aj=a;
	while (nj!=0)
	{
		if ((nj%2)==1)
		{
			rdummy=randlc( & power, aj);
		}
		rdummy=randlc( & aj, aj);
		nj=(nj/2);
	}
	return power;
}

/* --------------------------------------------------------------------- */
/* bubble        does a bubble sort in direction dir */
/* --------------------------------------------------------------------- */
static void bubble(double ten[][2], int j1[][2], int j2[][2], int j3[][2], int m, int ind)
{
	double temp;
	int i, j_temp;
	if (ind==1)
	{
		#pragma loop name bubble#0 
		for (i=0; i<(m-1); i ++ )
		{
			if (ten[i][ind]>ten[i+1][ind])
			{
				temp=ten[i+1][ind];
				ten[i+1][ind]=ten[i][ind];
				ten[i][ind]=temp;
				j_temp=j1[i+1][ind];
				j1[i+1][ind]=j1[i][ind];
				j1[i][ind]=j_temp;
				j_temp=j2[i+1][ind];
				j2[i+1][ind]=j2[i][ind];
				j2[i][ind]=j_temp;
				j_temp=j3[i+1][ind];
				j3[i+1][ind]=j3[i][ind];
				j3[i][ind]=j_temp;
			}
			else
			{
				return ;
			}
		}
	}
	else
	{
		#pragma loop name bubble#1 
		for (i=0; i<(m-1); i ++ )
		{
			if (ten[i][ind]<ten[i+1][ind])
			{
				temp=ten[i+1][ind];
				ten[i+1][ind]=ten[i][ind];
				ten[i][ind]=temp;
				j_temp=j1[i+1][ind];
				j1[i+1][ind]=j1[i][ind];
				j1[i][ind]=j_temp;
				j_temp=j2[i+1][ind];
				j2[i+1][ind]=j2[i][ind];
				j2[i][ind]=j_temp;
				j_temp=j3[i+1][ind];
				j3[i+1][ind]=j3[i][ind];
				j3[i][ind]=j_temp;
			}
			else
			{
				return ;
			}
		}
	}
}

static void zero3(void * oz, int n1, int n2, int n3)
{
	double (* z)[n2][n1] = (double (* )[n2][n1])oz;
	int i1, i2, i3;
	#pragma cetus private(i1, i2, i3) 
	#pragma loop name zero3#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*n3))+((3L*n2)*n3))+(((3L*n1)*n2)*n3)))) private(i1, i2, i3)
	for (i3=0; i3<n3; i3 ++ )
	{
		#pragma cetus private(i1, i2) 
		#pragma loop name zero3#0#0 
		for (i2=0; i2<n2; i2 ++ )
		{
			#pragma cetus private(i1) 
			#pragma loop name zero3#0#0#0 
			for (i1=0; i1<n1; i1 ++ )
			{
				z[i3][i2][i1]=0.0;
			}
		}
	}
}
