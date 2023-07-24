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
 * Written by the HYPRE team, UCRL-CODE-222953.
 * All rights reserved.
 *
 * This file is part of HYPRE (see http:www.llnl.gov/CASC/hypre/).
 * Please see the COPYRIGHT_and_LICENSE file for the copyright notice, 
 * disclaimer and the GNU Lesser General Public License.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the terms and conditions of the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * $Revision: 2.5 $
 ***********************************************************************EHEADER
*/
#include <stdio.h>
#include "utilities.h"
int hypre__global_error = 0;
/*
Process the error with code ierr raised in the given line of the
   given source file.
*/
void hypre_error_handler(char * filename, int line, int ierr)
{
	hypre__global_error|=ierr;
}

int HYPRE_GetError()
{
	return hypre__global_error;
}

int HYPRE_CheckError(int ierr, int hypre_error_code)
{
	return (ierr&hypre_error_code);
}

void HYPRE_DescribeError(int ierr, char * msg)
{
	if (ierr==0)
	{
		sprintf(msg, "[No error] ");
	}
	if (ierr&1)
	{
		sprintf(msg, "[Generic error] ");
	}
	if (ierr&2)
	{
		sprintf(msg, "[Memory error] ");
	}
	if (ierr&4)
	{
		sprintf(msg, "[Error in argument %d] ", HYPRE_GetErrorArg());
	}
	if (ierr&256)
	{
		sprintf(msg, "[Method did not converge] ");
	}
}

int HYPRE_GetErrorArg()
{
	return ((hypre__global_error>>3)&31);
}
