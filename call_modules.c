/*********************************************************
 * United States Geological Survey
 * used to call a Fortran version
 * $Id$
 *********************************************************/

#include <stdlib.h>
#include <string.h>
#include "mms.h"

extern long call_modules_ (char *, ftnlen);

int call_modules(char *arg) {
	 long retval;
	 ftnlen len;

	 len = (ftnlen)strlen(arg);
	 retval = call_modules_ (arg, len);
	 return((int)retval);
}
