/*+
 * United States Geological Survey
 * COMMENT  : saves the var data base to a file
 *             File name is passed in as arg
 * $Id$
-*/

/**1************************ INCLUDE FILES ****************************/
#define SAVE_VARS_C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : save_vars
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int save_vars (char *var_file_name) {
	FILE *var_file;
	PUBVAR *var;
	DIMEN *dim;
	long i,j;
	double *dvalptr;
	float *fvalptr;
	long *lvalptr;
    char *buf, *ptr;

/*
* get var file path, open file
*/
	if ((var_file = fopen (var_file_name, "w")) == NULL) {
		(void)fprintf(stderr, "ERROR - save_vars - creating file '%s'\n", var_file_name);
		return(1);
	}

/*
* write the run info string
*/
    buf = strdup (Mparaminfo);
    ptr = strchr (buf, '\n');

    if (ptr) *ptr = '\0';
   
	(void)fprintf (var_file, "%s\n", buf);

/*
* write nstep
*/
	(void)fprintf (var_file, "last nstep %ld\n", Mnsteps);

/*
* write Mnowtime->jt
*/
	(void)fprintf (var_file, "last julian time %f\n", Mnowtime->jt);

/*
* write delta time
*/
	(void)fprintf (var_file, "last delta time %f\n", Mdeltat);

/*
* write out dimensions
*/
	(void)fprintf(var_file, "** Dimensions **\n");
	for (i = 0; i < dim_db->count; i++) {
		dim = (DIMEN *)(dim_db->itm[i]);
		(void)fprintf(var_file, "####\n");
		(void)fprintf(var_file, "%s\n", dim->name);
		(void)fprintf(var_file, "%ld\n", dim->value);

	}

/*
* write out variable values
*/
	(void)fprintf(var_file, "** Variables **\n");

   for (i = 0; i < Mnvars; i++) {
      var = Mvarbase[i];
      (void)fprintf (var_file, "####\n");
      (void)fprintf (var_file, "%s\n", var->key);
      if (!(var->private)) {
         (void)fprintf (var_file, "%ld \n", var->ndimen);

         for (j = 0; j < var->ndimen; j++)
            (void)fprintf (var_file, "%s\n", var->dimen[j]->name);
      } else {
         (void)fprintf (var_file, "1\n");
         (void)fprintf (var_file, "PRIVATE\n");
      }

      (void)fprintf(var_file, "%ld \n", var->size);
      (void)fprintf(var_file, "%ld \n", var->type);

      switch (var->type) {
         case M_DOUBLE:
            dvalptr = (double *) var->value;
            for (j = 0; j < var->size; j++) {
               (void)fprintf(var_file, "%.20le \n", *dvalptr);
               dvalptr++;
            }
			break;

         case M_FLOAT:
            fvalptr = (float *) var->value;
            for (j = 0; j < var->size; j++) {
               (void)fprintf(var_file, "%.12e \n", *fvalptr);
               fvalptr++;
            }
            break;

         case M_LONG:
            lvalptr = (long *) var->value;
            for (j = 0; j < var->size; j++) {
               (void)fprintf(var_file, "%ld \n", *lvalptr);
               lvalptr++;
            }
            break;
      }

   }

   fclose(var_file);
   return(0);
}
