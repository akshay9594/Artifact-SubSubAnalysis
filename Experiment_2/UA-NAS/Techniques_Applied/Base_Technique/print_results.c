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
#include <stdio.h>
#include <math.h>
#include "type.h"
void print_results(char * name, char class, int n1, int n2, int n3, int niter, double t, double mops, char * optype, logical verified, char * npbversion, char * compiletime, char * cs1, char * cs2, char * cs3, char * cs4, char * cs5, char * cs6, char * cs7)
{
	char size[16];
	int j;
	printf("\n\n %s Benchmark Completed.\n", name);
	printf(" Class           =             %12c\n", class);
	/* If this is not a grid-based problem (EP, FT, CG), then */
	/* we only print n1, which contains some measure of the */
	/* problem size. In that case, n2 and n3 are both zero. */
	/* Otherwise, we print the grid size n1xn2xn3 */
	if ((n2==0)&&(n3==0))
	{
		if ((name[0]=='E')&&(name[1]=='P'))
		{
			sprintf(size, "%15.0lf", pow(2.0, n1));
			j=14;
			if (size[j]=='.')
			{
				size[j]=' ';
				j -- ;
			}
			size[j+1]='\0';
			printf(" Size            =          %15s\n", size);
		}
		else
		{
			printf(" Size            =             %12d\n", n1);
		}
	}
	else
	{
		printf(" Size            =           %4dx%4dx%4d\n", n1, n2, n3);
	}
	printf(" Iterations      =             %12d\n", niter);
	printf(" Time in seconds =             %12.2lf\n", t);
	printf(" Mop/s total     =          %15.2lf\n", mops);
	printf(" Operation type  = %24s\n", optype);
	if (verified)
	{
		printf(" Verification    =             %12s\n", "SUCCESSFUL");
	}
	else
	{
		printf(" Verification    =             %12s\n", "UNSUCCESSFUL");
	}
	printf(" Version         =             %12s\n", npbversion);
	printf(" Compile date    =             %12s\n", compiletime);
	printf("\n Compile options:\n""    CC           = %s\n", cs1);
	printf("    CLINK        = %s\n", cs2);
	printf("    C_LIB        = %s\n", cs3);
	printf("    C_INC        = %s\n", cs4);
	printf("    CFLAGS       = %s\n", cs5);
	printf("    CLINKFLAGS   = %s\n", cs6);
	printf("    RAND         = %s\n", cs7);
	printf("\n--------------------------------------\n"" Please send all errors/feedbacks to:\n"" Center for Manycore Programming\n"" cmp@aces.snu.ac.kr\n"" http://aces.snu.ac.kr\n""--------------------------------------\n\n");
}
