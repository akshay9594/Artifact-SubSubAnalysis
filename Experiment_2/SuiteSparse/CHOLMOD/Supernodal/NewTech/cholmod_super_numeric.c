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
/* ========================================================================== */
/* === Supernodalcholmod_super_numeric ===================================== */
/* ========================================================================== */
/*
-----------------------------------------------------------------------------
 CHOLMOD/Supernodal Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * http:www.suitesparse.com
 * --------------------------------------------------------------------------
*/
/*
Computes the Cholesky factorization of A+betaI or A*F+beta*I.  Only the
 * the lower triangular part of A+beta*I or A*F+beta*I is accessed.  The
 * matrices A and F must already be permuted according to the fill-reduction
 * permutation L->Perm.  cholmod_factorize is an "easy" wrapper for this code
 * which applies that permutation.  beta is real.
 *
 * Symmetric case: A is a symmetric (lower) matrix.  F is not accessed.
 * With a fill-reducing permutation, A(p,p) should be passed instead, where is
 * p is L->Perm.
 *
 * Unsymmetric case: A is unsymmetric, and F must be present.  Normally, F=A'.
 * With a fill-reducing permutation, A(p,f) and A(p,f)' should be passed as A
 * and F, respectively, where f is a list of the subset of the columns of A.
 *
 * The input factorization L must be supernodal (L->is_super is TRUE).  It can
 * either be symbolic or numeric.  In the first case, L has been analyzed by
 * cholmod_analyze or cholmod_super_symbolic, but the matrix has not yet been
 * numerically factorized.  The numerical values are allocated here and the
 * factorization is computed.  In the second case, a prior matrix has been
 * analyzed and numerically factorized, and a new matrix is being factorized.
 * The numerical values of L are replaced with the new numerical factorization.
 *
 * L->is_ll is ignored, and set to TRUE.  This routine always computes an LL'
 * factorization.  Supernodal LDL' factorization is not (yet) supported.
 * FUTURE WORK: perform a supernodal LDL' factorization if L->is_ll is FALSE.
 *
 * Uses BLAS routines dsyrk, dgemm, dtrsm, and the LAPACK routine dpotrf.
 * The supernodal solver uses BLAS routines dtrsv, dgemv, dtrsm, and dgemm.
 *
 * If the matrix is not positive definite the routine returns TRUE, but sets
 * Common->status to CHOLMOD_NOT_POSDEF and L->minor is set to the column at
 * which the failure occurred.  The supernode containing the non-positive
 * diagonal entry is set to zero (this includes columns to the left of L->minor
 * in the same supernode), as are all subsequent supernodes.
 *
 * workspace: Flag (nrow), Head (nrow+1), Iwork (2*nrow + 5*nsuper).
 *	Allocates temporary space of size L->maxcsize * sizeof(double)
 *	(twice that for the complex/zomplex case).
 *
 * If L is supernodal symbolic on input, it is converted to a supernodal numeric
 * factor on output, with an xtype of real if A is real, or complex if A is
 * complex or zomplex.  If L is supernodal numeric on input, its xtype must
 * match A (except that L can be complex and A zomplex).  The xtype of A and F
 * must match.

*/
#include <sys/time.h>
#include "cholmod_internal.h"
#include "cholmod_supernodal.h"
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include "cholmod_template.h"
/* ========================================================================== */
/* === TEMPLATE codes for GPU and regular numeric factorization ============= */
/* ========================================================================== */
/* #define REAL */
/* #include "t_cholmod_super_numeric.c" */
/* #define COMPLEX */
/* #include "t_cholmod_super_numeric.c" */
/* #define ZOMPLEX */
/* #include "t_cholmod_super_numeric.c" */
/* -------------------------------------------------------------------------- */
/* A, F, and L are all real */
/* -------------------------------------------------------------------------- */
/* ========================================================================== */
/* === cholmod_super_numeric ================================================ */
/* ========================================================================== */
/*
Returns TRUE if successful, or if the matrix is not positive definite.
 Returns FALSE if out of memory, inputs are invalid, or other fatal error
 * occurs.

*/
/* ---- input ---- */
/* matrix to factorize */
/* F = A' or A(:,f)' */
/* betaI is added to diagonal of matrix to factorize */
/* ---- inout --- */
/* factorization */
/* --------------- */
int CHOLMOD(super_numeric)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    double beta [2],	/* beta*I is added to diagonal of matrix to factorize */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization */
    /* --------------- */
    cholmod_common *Common
)
{
	cholmod_dense * C;
	int * Super, * Map, * SuperMap;
	size_t maxcsize;
	int nsuper, n, i, k, s, stype, nrow;
	int ok = 1, symbolic;
	size_t t, w;
	/* ---------------------------------------------------------------------- */
	/* check inputs */
	/* ---------------------------------------------------------------------- */
	double one[2], zero[2], tstart, sum = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, thread_time[8], starttime, endtime, threadtime;
	double * Lx, * Ax, * Fx, * Az, * Fz, * Cw;
	int * Head, * Ls, * Lpi, * Lpx, * RelativeMap, * Next, * Lpos, * Fp, * Fi, * Fnz, * Ap, * Ai, * Anz, * Iwork, * Next_save, * Lpos_save, * Previous;
	int j, p, pend, k1, k2, nscol, psi, psx, psend, nsrow, pj, d, kd1, kd2, info, ndcol, ndrow, pdi, pdx, pdend, pdi1, pdi2, pdx1, ndrow1, ndrow2, px, dancestor, sparent, dnext, nsrow2, ndrow3, pk, pf, pfend, Apacked, Fpacked, q, imap, repeat_supernode, nscol2, ss, tail, nscol_new = 0;
	double fjk[2];
	struct timeval start, end, begin, terminate, init, stop, st, ed;
	int k_0;
	int p_0;
	int k_1;
	int pf_0;
	int p_1;
	int k_2;
	int i_0;
	int p_2;
	int ss_0;
	int p_3;
	if (Common==((void * )0))
	{
		return 0;
	}
	if ((Common->itype!=0)||(Common->dtype!=0))
	{
		Common->status=( - 4);
		return 0;
	}
	;
	if (L==((void * )0))
	{
		if (Common->status!=( - 2))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 150, "argument missing", Common);
		}
		return 0;
	}
	;
	if (A==((void * )0))
	{
		if (Common->status!=( - 2))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 151, "argument missing", Common);
		}
		return 0;
	}
	;
	if ((((A->xtype<1)||(A->xtype>3))||((A->xtype!=0)&&(A->x==((void * )0))))||((A->xtype==3)&&(A->z==((void * )0))))
	{
		if (Common->status!=( - 2))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 152, "invalid xtype", Common);
		}
		return 0;
	}
	;
	if ((((L->xtype<0)||(L->xtype>2))||((L->xtype!=0)&&(L->x==((void * )0))))||((L->xtype==3)&&(L->z==((void * )0))))
	{
		if (Common->status!=( - 2))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 153, "invalid xtype", Common);
		}
		return 0;
	}
	;
	stype=A->stype;
	if (stype<0)
	{
		if ((A->nrow!=A->ncol)||(A->nrow!=L->n))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 159, "invalid dimensions", Common);
			return 0;
		}
	}
	else
	{
		if (stype==0)
		{
			if (A->nrow!=L->n)
			{
				cholmod_error( - 4, "cholmod_super_numeric.c", 167, "invalid dimensions", Common);
				return 0;
			}
			if (F==((void * )0))
			{
				if (Common->status!=( - 2))
				{
					cholmod_error( - 4, "cholmod_super_numeric.c", 170, "argument missing", Common);
				}
				return 0;
			}
			;
			if ((((F->xtype<1)||(F->xtype>3))||((F->xtype!=0)&&(F->x==((void * )0))))||((F->xtype==3)&&(F->z==((void * )0))))
			{
				if (Common->status!=( - 2))
				{
					cholmod_error( - 4, "cholmod_super_numeric.c", 171, "invalid xtype", Common);
				}
				return 0;
			}
			;
			if (((A->nrow!=F->ncol)||(A->ncol!=F->nrow))||(F->stype!=0))
			{
				cholmod_error( - 4, "cholmod_super_numeric.c", 174, "F invalid", Common);
				return 0;
			}
			if (A->xtype!=F->xtype)
			{
				cholmod_error( - 4, "cholmod_super_numeric.c", 179, "A and F must have same xtype", Common);
				return 0;
			}
		}
		else
		{
			/* symmetric upper case not suppored */
			cholmod_error( - 4, "cholmod_super_numeric.c", 186, "symmetric upper case not supported", Common);
			return 0;
		}
	}
	if ( ! L->is_super)
	{
		cholmod_error( - 4, "cholmod_super_numeric.c", 191, "L not supernodal", Common);
		return 0;
	}
	if (L->xtype!=0)
	{
		if ( ! ((((A->xtype==1)&&(L->xtype==1))||((A->xtype==2)&&(L->xtype==2)))||((A->xtype==3)&&(L->xtype==2))))
		{
			cholmod_error( - 4, "cholmod_super_numeric.c", 200, "complex type mismatch", Common);
			return 0;
		}
	}
	Common->status=0;
	/* ---------------------------------------------------------------------- */
	/* allocate workspace in Common */
	/* ---------------------------------------------------------------------- */
	nsuper=L->nsuper;
	maxcsize=L->maxcsize;
	nrow=A->nrow;
	n=nrow;
	;
	;
	/* w = 2n + 5*nsuper */
	w=cholmod_mult_size_t(n, 2,  & ok);
	t=cholmod_mult_size_t(nsuper, 5,  & ok);
	w=cholmod_add_size_t(w, t,  & ok);
	if ( ! ok)
	{
		cholmod_error( - 3, "cholmod_super_numeric.c", 224, "problem too large", Common);
		return 0;
	}
	cholmod_allocate_work(n, w, 0, Common);
	if (Common->status<0)
	{
		return 0;
	}
	;
	/* ---------------------------------------------------------------------- */
	/* get the current factor L and allocate numerical part, if needed */
	/* ---------------------------------------------------------------------- */
	Super=L->super;
	symbolic=(L->xtype==0);
	if (symbolic)
	{
		/* convert to supernodal numeric by allocating L->x */
		cholmod_change_factor(((A->xtype==1) ? 1 : 2), 1, 1, 1, 1, L, Common);
		if (Common->status<0)
		{
			/* the factor L remains in symbolic supernodal form */
			return 0;
		}
	}
	;
	;
	/* supernodal LDL' is not supported */
	L->is_ll=1;
	/* ---------------------------------------------------------------------- */
	/* get more workspace */
	/* ---------------------------------------------------------------------- */
	C=cholmod_allocate_dense(maxcsize, 1, maxcsize, L->xtype, Common);
	if (Common->status<0)
	{
		int status = Common->status;
		if (symbolic)
		{
			/*
			Change L back to symbolic, since the numeric values are not
				 initialized.  This cannot fail.
			*/
			cholmod_change_factor(0, 1, 1, 1, 1, L, Common);
		}
		/* the factor L is now back to the form it had on input */
		Common->status=status;
		return 0;
	}
	if (stype!=0)
	{
		
	}
	
	/* ---------------------------------------------------------------------- */
	/* get workspace */
	/* ---------------------------------------------------------------------- */
	SuperMap=Common->Iwork;
	/* size n (ii/l) */
	Map=Common->Flag;
	/* size n, use Flag as workspace for Map array */
	#pragma cetus private(i) 
	#pragma loop name cholmod_super_numeric#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*n)))) private(i)
	for (i=0; i<n; i ++ )
	{
		Map[i]=( - 1);
	}
	/* ---------------------------------------------------------------------- */
	/* find the mapping of nodes to relaxed supernodes */
	/* ---------------------------------------------------------------------- */
	/* SuperMap [k] = s if column k is contained in supernode s */
	#pragma cetus private(k, k_0, s) 
	#pragma loop name cholmod_super_numeric#1 
	for (s=0; s<nsuper; s ++ )
	{
		/* Normalized Loop */
		#pragma cetus lastprivate(k_0) 
		#pragma loop name cholmod_super_numeric#1#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*Super[(1L+s)]))+(-3L*Super[s])))) lastprivate(k_0)
		for (k_0=0; k_0<=((-1+Super[1+s])+(-1*Super[s])); k_0 ++ )
		{
			SuperMap[k_0+Super[s]]=s;
		}
		k=(k_0+Super[s]);
	}
	/* ---------------------------------------------------------------------- */
	/* supernodal numerical factorization, using template routine */
	/* ---------------------------------------------------------------------- */
	/* ---------------------------------------------------------------------- */
	/* guard against integer overflow in the BLAS */
	/* ---------------------------------------------------------------------- */
	/*
	If integer overflow occurs in the BLAS, Common->status is set to
	 CHOLMOD_TOO_LARGE, and the contents of Lx are undefined.
	*/
	Common->blas_ok=1;
	/* ---------------------------------------------------------------------- */
	/* get inputs */
	/* ---------------------------------------------------------------------- */
	nsuper=L->nsuper;
	n=L->n;
	Cw=C->x;
	/* workspace of size L->maxcsize */
	one[0]=1.0;
	/* ALPHA forsyrk, *herk, *gemm, and *trsm */
	one[1]=0.0;
	zero[0]=0.0;
	/* BETA forsyrk, *herk, and *gemm */
	zero[1]=0.0;
	/*
	Iwork must be of size 2n + 5nsuper, allocated in the caller,
	            * cholmod_super_numeric.  The memory cannot be allocated here because the
	            * cholmod_super_numeric initializes SuperMap, and cholmod_allocate_work
	            * does not preserve existing workspace if the space needs to be increase
	            * in size.
	*/
	/* allocate integer workspace */
	Iwork=Common->Iwork;
	SuperMap=Iwork;
	/* size n (ii/l) */
	RelativeMap=(Iwork+n);
	/* size n (ii/l) */
	Next=(Iwork+(2*((size_t)n)));
	/* size nsuper */
	Lpos=((Iwork+(2*((size_t)n)))+nsuper);
	/* size nsuper */
	Next_save=((Iwork+(2*((size_t)n)))+(2*((size_t)nsuper)));
	/* size nsuper */
	Lpos_save=((Iwork+(2*((size_t)n)))+(3*((size_t)nsuper)));
	/* size nsuper */
	Previous=((Iwork+(2*((size_t)n)))+(4*((size_t)nsuper)));
	/* size nsuper */
	/* Map  = Common->Flag ;   size n, use Flag as workspace for Map array */
	Head=Common->Head;
	/* size n+1, only Head [0..nsuper-1] used */
	Ls=L->s;
	Lpi=L->pi;
	Lpx=L->px;
	Lx=L->x;
	/* clear GPU CPU statistics */
	Common->cholmod_cpu_gemm_calls=0;
	Common->cholmod_cpu_syrk_calls=0;
	Common->cholmod_cpu_trsm_calls=0;
	Common->cholmod_cpu_potrf_calls=0;
	Common->cholmod_gpu_gemm_calls=0;
	Common->cholmod_gpu_syrk_calls=0;
	Common->cholmod_gpu_trsm_calls=0;
	Common->cholmod_gpu_potrf_calls=0;
	Common->cholmod_cpu_gemm_time=0;
	Common->cholmod_cpu_syrk_time=0;
	Common->cholmod_cpu_trsm_time=0;
	Common->cholmod_cpu_potrf_time=0;
	Common->cholmod_gpu_gemm_time=0;
	Common->cholmod_gpu_syrk_time=0;
	Common->cholmod_gpu_trsm_time=0;
	Common->cholmod_gpu_potrf_time=0;
	Common->cholmod_assemble_time=0;
	Common->cholmod_assemble_time2=0;
	stype=A->stype;
	stype=0;
	Fp=F->p;
	Fi=F->i;
	Fx=F->x;
	Fz=F->z;
	Fnz=F->nz;
	Fpacked=F->packed;
	Ap=A->p;
	Ai=A->i;
	Ax=A->x;
	Az=A->z;
	Anz=A->nz;
	Apacked=A->packed;
	/* clear the Map so that changes in the pattern of A can be detected */
	#pragma cetus private(i) 
	#pragma loop name cholmod_super_numeric#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*n)))) private(i)
	for (i=0; i<n; i ++ )
	{
		Map[i]=( - 1);
	}
	/*
	If the matrix is not positive definite, the supernode s containing the
	 first zero or negative diagonal entry of L is repeated (but factorized
	            * only up to just before the problematic diagonal entry). The purpose is
	            * to provide MATLAB with [R,p]=chol(A); columns 1 to p-1 of L=R' are
	            * required, where L(p,p) is the problematic diagonal entry.  The
	            * repeat_supernode flag tells us whether this is the repeated supernode.
	            * Once supernode s is repeated, the factorization is terminated.
	*/
	repeat_supernode=0;
	/* ---------------------------------------------------------------------- */
	/* supernodal numerical factorization */
	/* ---------------------------------------------------------------------- */
	#pragma cetus private(d, fjk, i, i_0, j, k, k_1, k_2, kd1, kd2, p, p_0, p_1, p_2, p_3, pdend, pdi, pdi1, pdi2, pdx, pdx1, pf, pf_0, pk, px, q, ss_0, tstart) 
	#pragma loop name cholmod_super_numeric#3 
	/* #pragma cetus reduction(+: Common->cholmod_cpu_gemm_calls, Common->cholmod_cpu_gemm_time, Common->cholmod_cpu_potrf_calls, Common->cholmod_cpu_potrf_time, Common->cholmod_cpu_syrk_calls, Common->cholmod_cpu_syrk_time, Common->cholmod_cpu_trsm_calls, Common->cholmod_cpu_trsm_time)  */
	for (s=0; s<nsuper; s ++ )
	{
		/* ------------------------------------------------------------------ */
		/* get the size of supernode s */
		/* ------------------------------------------------------------------ */
		k1=Super[s];
		/* s contains columns k1 to k2-1 of L */
		k2=Super[s+1];
		nscol=(k2-k1);
		/* # of columns in all of s */
		psi=Lpi[s];
		/* pointer to first row of s in Ls */
		psx=Lpx[s];
		/* pointer to first row of s in Lx */
		psend=Lpi[s+1];
		/* pointer just past last row of s in Ls */
		nsrow=(psend-psi);
		/* # of rows in all of s */
		/* ------------------------------------------------------------------ */
		/* zero the supernode s */
		/* ------------------------------------------------------------------ */
		;
		pend=(psx+(nsrow*nscol));
		/* s is nsrow-by-nscol */
		/* Case of no GPU, zero individual supernodes */
		/* #pragma omp parallel for num_threads(CHOLMOD_OMP_NUM_THREADS)               schedule (static) if ( pend - psx > 1024 ) */
		/* Normalized Loop */
		#pragma cetus lastprivate(p_0) 
		#pragma loop name cholmod_super_numeric#3#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((((1L+((-3L*k1)*psend))+((3L*k1)*psi))+((3L*k2)*psend))+((-3L*k2)*psi)))) lastprivate(p_0)
		for (p_0=0; p_0<=((-1+pend)+(-1*psx)); p_0 ++ )
		{
			/* L_CLEAR (Lx,p); */
			Lx[p_0+psx]=0;
		}
		p=(p_0+psx);
		/* ------------------------------------------------------------------ */
		/* construct the scattered Map for supernode s */
		/* ------------------------------------------------------------------ */
		/*
		If row i is the kth row in s, then Map [i] = k.  Similarly, if
		 column j is the kth column in s, then  Map [j] = k.
		*/
		/* #pragma omp parallel for num_threads(CHOLMOD_OMP_NUM_THREADS)               if ( nsrow > 128 ) */
		#pragma cetus private(k) 
		#pragma loop name cholmod_super_numeric#3#1 
		for (k=0; k<nsrow; k ++ )
		{
			Map[Ls[psi+k]]=k;
		}
		/* ------------------------------------------------------------------ */
		/* when using GPU, reorder supernodes by levels. */
		/* (all supernodes in a level are independent) */
		/* ------------------------------------------------------------------ */
		pk=psx;
		/* stype = 0; */
		gettimeofday(&start,NULL);
		/* #pragma omp parallel for private ( p, pend, pfend, pf, i, j, imap, q )              num_threads(CHOLMOD_OMP_NUM_THREADS) if ( k2-k1 > 64 ) */
		/* Normalized Loop */
		#pragma cetus private(fjk, j, p, p_1, pf, pf_0) 
		#pragma cetus lastprivate(k_1) 
		#pragma loop name cholmod_super_numeric#3#2 
		/* #pragma cetus reduction(+: Lx[((Map[Ai[(p_1+Ap[j])]]+psx)+(((k1+k_1)-k1)nsrow))])  */
		#pragma cetus parallel 
		#pragma omp parallel for private(fjk, j, p, p_1, pf, pf_0) lastprivate(k_1)
		for (k_1=0; k_1<=((-1+(-1*k1))+k2); k_1 ++ )
		{
			/* Normalized Loop */
			#pragma cetus private(fjk, j, p, p_1) 
			#pragma cetus lastprivate(pf_0) 
			#pragma loop name cholmod_super_numeric#3#2#0 
			/* #pragma cetus reduction(+: Lx[((Map[Ai[(p_1+Ap[j])]]+psx)+(((k1+k_1)-k1)nsrow))])  */
			for (pf_0=0; pf_0<=((-1+Fp[(1+k1)+k_1])+(-1*Fp[k1+k_1])); pf_0 ++ )
			{
				j=Fi[pf_0+Fp[k1+k_1]];
				fjk[0]=Fx[pf_0+Fp[k1+k_1]];
				/* Normalized Loop */
				#pragma cetus lastprivate(p_1) 
				#pragma loop name cholmod_super_numeric#3#2#0#0 
				/* #pragma cetus reduction(+: Lx[((Map[Ai[(p_1+Ap[j])]]+psx)+(((k1+k_1)-k1)nsrow))])  */
				for (p_1=0; p_1<=((-1+Ap[1+j])+(-1*Ap[j])); p_1 ++ )
				{
					if (Ai[p_1+Ap[j]]>=(k1+k_1))
					{
						/* See the discussion of imap above. */
						if ((Map[Ai[p_1+Ap[j]]]>=0)&&(Map[Ai[p_1+Ap[j]]]<nsrow))
						{
							/* Lx [Map [i] + pk] += Ax [p] fjk ; */
							Lx[(Map[Ai[p_1+Ap[j]]]+psx)+(((k1+k_1)-k1)*nsrow)]+=(Ax[p_1+Ap[j]]*fjk[0]);
						}
					}
				}
				p=(p_1+Ap[j]);
			}
			pf=(pf_0+Fp[k1+k_1]);
		}
		k=(k1+k_1);
		gettimeofday(&end, NULL);
		sum += (end.tv_sec + (double)end.tv_usec/1000000) - (start.tv_sec + (double)start.tv_usec/1000000);
		/* add beta to the diagonal of the supernode, if nonzero */
		if (beta[0]!=0.0)
		{
			/* note that only the real part of beta is used */
			pk=psx;
			/* Normalized Loop */
			#pragma cetus lastprivate(k_2) 
			#pragma loop name cholmod_super_numeric#3#3 
			/* #pragma cetus reduction(+: Lx[pk])  */
			for (k_2=0; k_2<=((-1+(-1*k1))+k2); k_2 ++ )
			{
				/* Lx [pk] += beta [0] ; */
				Lx[pk]+=beta[0];
				pk+=(nsrow+1);
				/* advance to the next diagonal entry */
			}
			k=(k1+k_2);
		}
		/* ------------------------------------------------------------------ */
		/* saverestore the list of supernodes */
		/* ------------------------------------------------------------------ */
		if ( ! repeat_supernode)
		{
			/*
			Save the list of pending descendants in case s is not positive
			 definite.  Also save Lpos for each descendant d, so that we can
			                    * find which part of d is used to update s.
			*/
			#pragma cetus private(d) 
			#pragma loop name cholmod_super_numeric#3#4 
			for (d=Head[s]; d!=( - 1); d=Next[d])
			{
				Lpos_save[d]=Lpos[d];
				Next_save[d]=Next[d];
			}
		}
		else
		{
			#pragma cetus private(d) 
			#pragma loop name cholmod_super_numeric#3#5 
			for (d=Head[s]; d!=( - 1); d=Next[d])
			{
				Lpos[d]=Lpos_save[d];
				Next[d]=Next_save[d];
			}
		}
		/* ------------------------------------------------------------------ */
		/* update supernode s with each pending descendant d */
		/* ------------------------------------------------------------------ */
		/* GPU module not installed */
		dnext=Head[s];
		while (dnext!=( - 1))
		{
			/* GPU module not installed at compile time */
			d=dnext;
			/* -------------------------------------------------------------- */
			/* get the size of supernode d */
			/* -------------------------------------------------------------- */
			kd1=Super[d];
			/* d contains cols kd1 to kd2-1 of L */
			kd2=Super[d+1];
			ndcol=(kd2-kd1);
			/* # of columns in all of d */
			pdi=Lpi[d];
			/* pointer to first row of d in Ls */
			pdx=Lpx[d];
			/* pointer to first row of d in Lx */
			pdend=Lpi[d+1];
			/* pointer just past last row of d in Ls */
			ndrow=(pdend-pdi);
			/* # rows in all of d */
			;
			;
			/* -------------------------------------------------------------- */
			/* find the range of rows of d that affect rows k1 to k2-1 of s */
			/* -------------------------------------------------------------- */
			p=Lpos[d];
			/* offset of 1st row of d affecting s */
			pdi1=(pdi+p);
			/* ptr to 1st row of d affecting s in Ls */
			pdx1=(pdx+p);
			/* ptr to 1st row of d affecting s in Lx */
			/* there must be at least one row remaining in d to update s */
			;
			;
			;
			#pragma cetus lastprivate(pdi2) 
			#pragma loop name cholmod_super_numeric#3#6 
			for (pdi2=pdi1; (pdi2<pdend)&&(Ls[pdi2]<k2); pdi2 ++ )
			{
				;
			}
			ndrow1=(pdi2-pdi1);
			/* # rows in first part of d */
			ndrow2=(pdend-pdi1);
			/* # rows in remaining d */
			/*
			rows Ls [pdi1 ... pdi2-1] are in the range k1 to k2-1.  Since d
			 affects s, this set cannot be empty.
			*/
			;
			;
			;
			/* -------------------------------------------------------------- */
			/* construct the update matrix C for this supernode d */
			/* -------------------------------------------------------------- */
			/*
			C = L (k1:n-1, kd1:kd2-1) L (k1:k2-1, kd1:kd2-1)', except
			                    * that k1:n-1 refers to all of the rows in L, but many of the
			                    * rows are all zero.  Supernode d holds columns kd1 to kd2-1 of L.
			                    * Nonzero rows in the range k1:k2-1 are in the list
			                    * Ls [pdi1 ... pdi2-1], of size ndrow1.  Nonzero rows in the range
			                    * k2:n-1 are in the list Ls [pdi2 ... pdend], of size ndrow2.  Let
			                    * L1 = L (Ls [pdi1 ... pdi2-1], kd1:kd2-1), and let
			                    * L2 = L (Ls [pdi2 ... pdend],  kd1:kd2-1).  C is ndrow2-by-ndrow1.
			                    * Let C1 be the first ndrow1 rows of C and let C2 be the last
			                    * ndrow2-ndrow1 rows of C.  Only the lower triangular part of C1
			                    * needs to be computed since C1 is symmetric.
			                   
			*/
			/* maxcsize is the largest size of C for all pairs (d,s) */
			;
			/*
			compute leading ndrow1-by-ndrow1 lower triangular block of C,
			 C1 = L1*L1'
			*/
			ndrow3=(ndrow2-ndrow1);
			/* number of rows of C2 */
			;
			gettimeofday(&st, NULL);
			{
				/* GPU not installed, or not used */
				Common->cholmod_cpu_syrk_calls ++ ;
				tstart=SuiteSparse_time();
				{
					int N = ndrow1, K = ndcol, LDA = ndrow, LDC = ndrow2;
					;
					/* N, K: L1 is ndrow1-by-ndcol */
					/* ALPHA:  1 */
					/* A, LDA: L1, ndrow */
					/* BETA:   0 */
					dsyrk_("L", "N",  & N,  & K, one, Lx+(1*pdx1),  & LDA, zero, Cw,  & LDC);
				}
				;
				/* C, LDC: C1 */
				Common->cholmod_cpu_syrk_time+=(SuiteSparse_time()-tstart);
				/*
				compute remaining (ndrow2-ndrow1)-by-ndrow1 block of C,
				 C2 = L2*L1'
				*/
				if (ndrow3>0)
				{
					Common->cholmod_cpu_gemm_calls ++ ;
					tstart=SuiteSparse_time();
					{
						/* A, LDA: L2 */
						/* B, LDB: L1 */
						/* C, LDC: C2 */
						int M = ndrow3, N = ndrow1, K = ndcol, LDA = ndrow, LDB = ndrow, LDC = ndrow2;
						/* A, LDA: L2 */
						/* B, LDB: L1 */
						/* C, LDC: C2 */
						;
						/* M, N, K */
						/* ALPHA:  1 */
						/* ndrow */
						/* ndrow */
						/* BETA:   0 */
						dgemm_("N", "C",  & M,  & N,  & K, one, Lx+(1*(pdx1+ndrow1)),  & LDA, Lx+(1*pdx1),  & LDB, zero, Cw+(1*ndrow1),  & LDC);
					}
					;
					Common->cholmod_cpu_gemm_time+=(SuiteSparse_time()-tstart);
				}
				gettimeofday(&ed, NULL);
				sum3 += (ed.tv_sec + (double)ed.tv_usec/1000000) - (st.tv_sec + (double)st.tv_usec/1000000);
				/* ---------------------------------------------------------- */
				/* construct relative map to assemble d into s */
				/* ---------------------------------------------------------- */
				;
				if ((Ls[(pdi1+ndrow2)-1]>Ls[(psi+nsrow)-1])&&(Ls[pdi1]<Ls[0]))
				{
					break;
				}
				/* #pragma omp parallel for num_threads(CHOLMOD_OMP_NUM_THREADS)               if ( ndrow2 > 64 )schedule(dynamic) */
				#pragma cetus private(i) 
				#pragma loop name cholmod_super_numeric#3#7 
				#pragma cetus parallel 
				#pragma omp parallel for if((10000<(1L+(4L*ndrow2)))) private(i)
				for (i=0; i<ndrow2; i ++ )
				{
					RelativeMap[i]=Map[Ls[pdi1+i]];
					;
				}
				/* ---------------------------------------------------------- */
				/* assemble C into supernode s using the relative map */
				/* ---------------------------------------------------------- */
				 gettimeofday(&init, NULL);
				/* #pragma omp parallel for private ( j, i, px, q )                                        num_threads(CHOLMOD_OMP_NUM_THREADS) */
				/* cols k1:k2-1 */
				#pragma cetus private(i, i_0, j, px, q) 
				#pragma loop name cholmod_super_numeric#3#8 
				/* #pragma cetus reduction(+: )  */
				for (j=0; j<ndrow1; j ++ )
				{
					/* printf("executed\n"); */
					px=(psx+(RelativeMap[j]*nsrow));
					/* rows k1:n-1 */
					/* Normalized Loop */
					#pragma cetus private(q) 
					#pragma cetus lastprivate(i_0) 
					#pragma loop name cholmod_super_numeric#3#8#0 
					/* #pragma cetus reduction(+: Lx[q])  */
					for (i_0=0; i_0<=((-1+(-1*j))+ndrow2); i_0 ++ )
					{
						q=(px+RelativeMap[i_0+j]);
						Lx[q]-=Cw[(i_0+j)+(ndrow2*j)];
					}
					i=(i_0+j);
				}
				gettimeofday(&stop, NULL);
				sum1 += (stop.tv_sec + (double)stop.tv_usec/1000000) - (init.tv_sec + (double)init.tv_usec/1000000);
			}
			/* -------------------------------------------------------------- */
			/* prepare this supernode d for its next ancestor */
			/* -------------------------------------------------------------- */
			dnext=Next[d];
			if ( ! repeat_supernode)
			{
				/*
				If node s is being repeated, Head [dancestor] has already
				 been cleared (set to EMPTY).  It must remain EMPTY.  The
				                        * dancestor will not be factorized since the factorization
				                        * terminates at node s.
				*/
				Lpos[d]=(pdi2-pdi);
				if (Lpos[d]<ndrow)
				{
					dancestor=SuperMap[Ls[pdi2]];
					;
					/* place d in the link list of its next ancestor */
					Next[d]=Head[dancestor];
					Head[dancestor]=d;
				}
			}
		}
		/* end of descendant supernode loop */
		/* ------------------------------------------------------------------ */
		/* factorize diagonal block of supernode s in LL' */
		/* ------------------------------------------------------------------ */
		/*
		The current supernode s is ready to factorize.  It has been updated
		 by all descendant supernodes.  Let S = the current supernode, which
		                * holds rows k1:n-1 and columns k1:k2-1 of the updated matrix.   It
		                * splits into two parts:  the square diagonal block S1, and the
		                * rectangular part S2.  Here, S1 is factorized into L1*L1' and
		                * overwritten by L1.
		                *
		                * If supernode s is being repeated, only factorize it up to but not
		                * including the column containing the problematic entry.
		               
		*/
		nscol2=(repeat_supernode ? nscol_new : nscol);
		gettimeofday(&begin, NULL);
		{
			/* Note that the GPU will not be used for the triangular solve */
			Common->cholmod_cpu_potrf_calls ++ ;
			tstart=SuiteSparse_time();
			{
				int N = nscol2, LDA = nsrow, INFO = 1;
				;
				/* N: nscol2 */
				dpotrf_("L",  & N, Lx+(1*psx),  & LDA,  & INFO);
				/* A, LDA: S1, nsrow */
				info=INFO;
			}
			;
			/* INFO */
			Common->cholmod_cpu_potrf_time+=(SuiteSparse_time()-tstart);
		}
		gettimeofday(&terminate, NULL);
		sum2 += (terminate.tv_sec + (double)terminate.tv_usec/1000000) - (begin.tv_sec + (double)begin.tv_usec/1000000);
		/* ------------------------------------------------------------------ */
		/* check if the matrix is not positive definite */
		/* ------------------------------------------------------------------ */
		if (repeat_supernode)
		{
			/* the leading part has been refactorized; it must have succeeded */
			info=0;
			/* zero out the rest of this supernode */
			pend=(psx+(nsrow*nscol));
			/* s is nsrow-by-nscol */
			/* Normalized Loop */
			#pragma cetus lastprivate(p_2) 
			#pragma loop name cholmod_super_numeric#3#9 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<((1L+((3L*nscol)*nsrow))+((-3L*nscol_new)*nsrow)))) lastprivate(p_2)
			for (p_2=0; p_2<=(((-1+pend)+(-1*psx))+((-1*nscol_new)*nsrow)); p_2 ++ )
			{
				/* Lx [p] = 0 ; */
				Lx[(p_2+psx)+(nscol_new*nsrow)]=0;
			}
			p=((p_2+psx)+(nscol_new*nsrow));
		}
		/*
		info is set to one in LAPACK_potrf if blas_ok is FALSE.  It is
		                * set to zero in dpotrf/zpotrf if the factorization was successful.
		*/
		;
		if (info!=0)
		{
			/*
			Matrix is not positive definite.  dpotrfzpotrf do NOT report an
			                    * error if the diagonal of L has NaN's, only if it has a zero.
			*/
			if (Common->status==0)
			{
				cholmod_error(1, "cholmod_super_numeric.c", 1025, "matrix not positive definite", Common);
			}
			/*
			L->minor is the column of L that contains a zero or negative
			 diagonal term.
			*/
			L->minor=((k1+info)-1);
			/* clear the link lists of all subsequent supernodes */
			/* Normalized Loop */
			#pragma cetus lastprivate(ss_0) 
			#pragma loop name cholmod_super_numeric#3#10 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<((-2L+(3L*nsuper))+(-3L*s)))) lastprivate(ss_0)
			for (ss_0=0; ss_0<=((-2+nsuper)+(-1*s)); ss_0 ++ )
			{
				Head[(1+s)+ss_0]=( - 1);
			}
			ss=((1+s)+ss_0);
			/* zero this supernode, and all remaining supernodes */
			pend=L->xsize;
			/* Normalized Loop */
			#pragma cetus lastprivate(p_3) 
			#pragma loop name cholmod_super_numeric#3#11 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<((1L+(3L*pend))+(-3L*psx)))) lastprivate(p_3)
			for (p_3=0; p_3<=((-1+pend)+(-1*psx)); p_3 ++ )
			{
				/* Lx [p] = 0. ; */
				Lx[p_3+psx]=0;
			}
			p=(p_3+psx);
			/*
			If L is indefinite, it still contains useful information.
			 Supernodes 0 to s-1 are valid, similar to MATLAB [R,p]=chol(A),
			                    * where the 1-based p is identical to the 0-based L->minor.  Since
			                    * L->minor is in the current supernode s, it and any columns to the
			                    * left of it in supernode s are also all zero.  This differs from
			                    * [R,p]=chol(A), which contains nonzero rows 1 to p-1.  Fix this
			                    * by setting repeat_supernode to TRUE, and repeating supernode s.
			                    *
			                    * If Common->quick_return_if_not_posdef is true, then the entire
			                    * supernode s is not factorized; it is left as all zero.
			                   
			*/
			if ((info==1)||Common->quick_return_if_not_posdef)
			{
				/*
				If the first column of supernode s contains a zero or
				 negative diagonal entry, then it is already properly set to
				                        * zero.  Also, info will be 1 if integer overflow occured in
				                        * the BLAS.
				*/
				Head[s]=( - 1);
			}
			else
			{
				/*
				Repeat supernode s, but only factorize it up to but not
				 including the column containing the problematic diagonal
				                        * entry.
				*/
				repeat_supernode=1;
				s -- ;
				nscol_new=(info-1);
				continue;
			}
		}
		/* ------------------------------------------------------------------ */
		/* compute the subdiagonal block and prepare supernode for its parent */
		/* ------------------------------------------------------------------ */
		nsrow2=(nsrow-nscol2);
		if (nsrow2>0)
		{
			/*
			The current supernode is columns k1 to k2-1 of L.  Let L1 be the
			 diagonal block (factorized by dpotrf/zpotrf above; rows/cols
			                    * k1:k2-1), and L2 be rows k2:n-1 and columns k1:k2-1 of L.  The
			                    * triangular system to solve is L2*L1' = S2, where S2 is
			                    * overwritten with L2.  More precisely, L2 = S2 / L1' in MATLAB
			                    * notation.
			                   
			*/
			gettimeofday(&st, NULL);
			{
				Common->cholmod_cpu_trsm_calls ++ ;
				tstart=SuiteSparse_time();
				{
					/* B, LDB, L2, nsrow */
					int M = nsrow2, N = nscol2, LDA = nsrow, LDB = nsrow;
					/* B, LDB, L2, nsrow */
					;
					/* M, N */
					/* ALPHA: 1 */
					/* A, LDA: L1, nsrow */
					dtrsm_("R", "L", "C", "N",  & M,  & N, one, Lx+(1*psx),  & LDA, Lx+(1*(psx+nscol2)),  & LDB);
				}
				;
				Common->cholmod_cpu_trsm_time+=(SuiteSparse_time()-tstart);
			}
			gettimeofday(&ed, NULL);
			sum3 += (ed.tv_sec + (double)ed.tv_usec/1000000) - (st.tv_sec + (double)st.tv_usec/1000000); 
			;
			if ( ! repeat_supernode)
			{
				/* Lpos [s] is offset of first row of s affecting its parent */
				Lpos[s]=nscol;
				sparent=SuperMap[Ls[psi+nscol]];
				;
				;
				;
				;
				;
				/* place s in link list of its parent */
				Next[s]=Head[sparent];
				Head[sparent]=s;
			}
		}
		else
		{
			
		}
		Head[s]=( - 1);
		/* link list for supernode s no longer needed */
		/* clear the Map (debugging only, to detect changes in pattern of A) */
		;
		;
		if (repeat_supernode)
		{
			/*
			matrix is not positive definite; finished clean-up for supernode
			 containing negative diagonal
			*/
			return (Common->status>=0);
		}
	}
	/* success; matrix is positive definite */
	L->minor=n;
	printf("Time taken by 1st subsub loop = %lf s\n", sum);
	printf("Time taken by 2nd subsub loop = %lf s\n", sum1);
	printf("Time taken by lapack+BLAS routines = %lf s\n", sum2+sum3);
	printf("Time taken by kernel= %lf s\n", ((sum+sum1)+sum2)+sum3);
	ok=0;
	/* ---------------------------------------------------------------------- */
	/* clear Common workspace, free temp workspace C, and return */
	/* ---------------------------------------------------------------------- */
	/* Flag array was used as workspace, clear it */
	Common->mark=( - 1);
	/* CHOLMOD(clear_flag) (Common) ; */
	Common->mark ++ ;
	if (Common->mark<=0)
	{
		Common->mark=( - 1);
		cholmod_clear_flag(Common);
	}
	;
	;
	cholmod_free_dense( & C, Common);
	return ok;
}
