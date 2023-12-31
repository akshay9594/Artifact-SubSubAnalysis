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
/*  */
/*      C  _  P  R  I  N  T  _  R  E  S  U  L  T  S */
/*  */
#include <stdlib.h>
#include <stdio.h>
void c_print_results(char * name, char class, int n1, int n2, int n3, int niter, double t, double mops, char * optype, int passed_verification, char * npbversion, char * compiletime, char * cc, char * clink, char * c_lib, char * c_inc, char * cflags, char * clinkflags)
{
	printf("\n\n %s Benchmark Completed\n", name);
	printf(" Class           =                        %c\n", class);
	if (n3==0)
	{
		long nn = n1;
		if (n2!=0)
		{
			nn*=n2;
		}
		printf(" Size            =             %12ld\n", nn);
		/* as in IS */
	}
	else
	{
		printf(" Size            =             %4dx%4dx%4d\n", n1, n2, n3);
	}
	printf(" Iterations      =             %12d\n", niter);
	printf(" Time in seconds =             %12.2f\n", t);
	printf(" Mop/s total     =             %12.2f\n", mops);
	printf(" Operation type  = %24s\n", optype);
	if (passed_verification<0)
	{
		printf(" Verification    =            NOT PERFORMED\n");
	}
	else
	{
		if (passed_verification)
		{
			printf(" Verification    =               SUCCESSFUL\n");
		}
		else
		{
			printf(" Verification    =             UNSUCCESSFUL\n");
		}
	}
	printf(" Version         =             %12s\n", npbversion);
	printf(" Compile date    =             %12s\n", compiletime);
	printf("\n Compile options:\n");
	printf("    CC           = %s\n", cc);
	printf("    CLINK        = %s\n", clink);
	printf("    C_LIB        = %s\n", c_lib);
	printf("    C_INC        = %s\n", c_inc);
	printf("    CFLAGS       = %s\n", cflags);
	printf("    CLINKFLAGS   = %s\n", clinkflags);
	printf("\n--------------------------------------\n");
	printf(" Please send all errors/feedbacks to:\n");
	printf(" Center for Manycore Programming\n");
	printf(" cmp@aces.snu.ac.kr\n");
	printf(" http://aces.snu.ac.kr\n");
	printf("--------------------------------------\n");
}
