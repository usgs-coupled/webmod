/*
 * $Id: print_model_info.c 4098 2008-04-23 19:09:09Z markstro $
 */
#define PRINT_MODEL_INFO_C
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mms.h"

#define PRINTLEN 77

//typedef struct node_type {
//   char  *name;
//   int   x, y;
//   int   num_connectnions;
//} NODE;

//extern NODE nodes[];
//extern int con_index[];

/*--------------------------------------------------------------------*\
 | FUNCTION     : print_model_info
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int print_model_info (void) {

  char pathname[MAXPATHLEN];
  FILE *param_file;
  //NODE *np;
  char *module_names, str[21], *ptr, *end;
  int *xpos, *ypos, *con_index, *num_connectnions;
  int num_modules;
  long i, j, k;
  CONTROL *cp;

  /*
   * get param file path name, open file
   */


  (void)sprintf (pathname, "%s.mod_name", MAltContFile);

  if ((param_file = fopen (pathname, "w")) == NULL) {
    (void)fprintf(stderr,
	    "ERROR - print_model_info - creating file '%s'\n", pathname);
    perror("");
    return(1);
  }

  /*
   * write header
   */

  (void)fprintf(param_file, "%s\n", model_name);
  (void)fprintf(param_file, "============\n\n");

  (void)fprintf(param_file, "Printout of module call order, X, Y, and number of connections.\n\n");

  cp = control_addr ("module_names");
  num_modules = cp->size;
   module_names = (char *)(cp->start_ptr);

   cp = control_addr ("module_x");
   xpos = (int *)(cp->start_ptr);

   cp = control_addr ("module_y");
   ypos = (int *)(cp->start_ptr);

   cp = control_addr ("module_num_con");
   num_connectnions = (int *)(cp->start_ptr);

   cp = control_addr ("module_connections");
   con_index = (int *)(cp->start_ptr);

   ptr = module_names;
   str[20] = '\0';
   for (i = 0; i < num_modules; i++) {
		strncpy (str, "                     ", 20);
		strncpy (str, ptr, 20);
		end = strchr(str, ' ');
		if (end) *(end) = '\0';
	   (void)fprintf(param_file, "%s, %d, %d, %d\n", str, xpos[i], ypos[i], num_connectnions[i]);
	   ptr = ptr + 20;
   }

  (void)fprintf(param_file, "\n\nPrintout of connections.\n\n");

  j = 0;
   for (k = 0; k < num_modules; k++) {
       for (i = 0; i < num_connectnions[k]; i++) {
	      fprintf(param_file, "%d ", con_index[j++]);
       }
	}
	fprintf(param_file, "\n");


  fclose(param_file);

  return(0);

}
