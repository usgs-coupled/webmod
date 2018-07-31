/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * NAME     : read_control.c
 * AUTHOR   : CADSWES
 * DATE     : Mon 08 Apr 1996
 * FUNCTION :
 * COMMENT  :
 * read_control.c: reads the control data base from a file
 * File name is obtained from the environment variable "mms_control_file"
 * REF      :
 * REVIEW   :
 * PR NRS   :
 *
 * $Id: read_control.c 4627 2008-10-01 16:48:11Z markstro $
 *
   $Revision: 4627 $
        $Log: read_control.c,v $
        Revision 1.21  1996/08/28 15:24:10  markstro
        Unknown

 * Revision 1.20  1996/05/14  02:42:05  msbrewer
 * Cleaned up cvs conflicts. Bug fixes in dump_to_db.
 *
        Revision 1.19  1996/04/29 16:23:09  markstro
        Unknown

 * Revision 1.18  1996/04/09  21:04:09  markstro
 * (1) Work on control files
 * (2) Runtime graphs
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_CONTROL_C
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**2************************* LOCAL MACROS ****************************/

/**3************************ LOCAL TYPEDEFS ***************************/

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static char *rc (char *);
char *fgets_rc (char *, int , FILE *);

/**5*********************** LOCAL VARIABLES ***************************/

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : read_control
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_control (char *control_name) {
   static char *foo = NULL;
   char old[256], *cptr;

   if (foo) {
      strncpy (old, foo, 256);
      free (foo);
      foo = strdup (control_name);
   } else {
      strncpy (old, control_name, 256);
      foo = strdup (control_name);
   }

   cptr = rc (control_name);

   if (cptr) {
      rc (old);

      free (foo);
      foo = strdup (old);
   }

   return (cptr);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : rc
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *rc (char *control_name) {
   FILE   *control_file;
   CONTROL *cp;
   long   type, size, i;
   double   *dptr;
   float   *fptr;
   long   *lptr;
   char   **cptr;
   char   line[MAXLNLEN], key[MAXKEYLEN];
   static char      buf[256];

/*
* compute control path, open file
*/
   if ((control_file = fopen (control_name, "r")) == NULL) {
      (void)sprintf (buf, "read_control: Couldn't open %s", control_name);
      return (buf);
   }

   if (!fgets_rc(line, MAXLNLEN, control_file)) {
      fclose (control_file);
      (void)sprintf (buf, "read_control: Problems reading %s", control_name);
      return (buf);
   }

/*
**   Read in control variables.
*/
   while (!feof (control_file)) {

/*
**   Space fwd to #### header.
*/
      while (strncmp(line, "####", 4)) {
         if (fgets_rc(line, MAXLNLEN, control_file) == NULL) {
            fclose(control_file);
            return(NULL);
         }
      }
/*
**   get key
*/
      if (!fgets_rc (key, MAXKEYLEN, control_file)) {
         (void)sprintf (buf, "read_control: reading key; Early end-of-file");
         return (buf);
      }

      *(key + strlen (key) - 1) = '\0';
      cp = control_addr (key);
/*
      if (cp) {
         (void)sprintf (buf, "read_control: control variable %s already declared", key);
         return (buf);
      }
*/

/*
**   get size
*/
      if (!fgets_rc (line, MAXLNLEN, control_file)) {
         (void)sprintf (buf,"read_control: reading size; key = %s", key);
         return (buf);
      }

      if ((size = atol(line)) < 0) {
         (void)sprintf (buf, "read_control: negative size; key = %s, line = %s", key, line);
         return (buf);
      }
/*
**   get type
*/
      if (!fgets_rc (line, MAXLNLEN, control_file)) {
         (void)sprintf (buf, "read_control: reading type; key = %s", key);
         return (buf);
      }

      if (!(type = atol(line))) {
        (void)sprintf (buf, "read_control: invalid type; key = %s, line = %s", key, line);
        return (buf);
      }

      if (!cp) cp = add_control (key, type, size);

      switch (type) {
         case M_DOUBLE:
            dptr = (double *)(cp->start_ptr);
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXLNLEN, control_file) == NULL) {
                  (void)sprintf (buf, "read_control: key is %s.\n, file: %s", key, control_name);
                  return (buf);
               }
               dptr[i] = atof(line);
            }
            break;

         case M_FLOAT:
            fptr = (float *) cp->start_ptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXLNLEN, control_file) == NULL) {
                  (void)sprintf (buf, "read_control: key is %s.\n, file: %s", key, control_name);
                  return (buf);
               }
               fptr[i] = (float) atof(line);
            }
            break;

         case M_LONG:
            lptr = (long *) cp->start_ptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXLNLEN, control_file) == NULL) {
                  (void)sprintf (buf, "read_control: key is %s.\n, file: %s", key, control_name);
                  return (buf);
               }
               lptr[i] =  atol(line);
            }
            break;

         case M_STRING:
            cptr = (char **) cp->start_ptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXLNLEN, control_file) == NULL) {
                  (void)sprintf (buf, "read_control: key is %s.\n, file: %s", key, control_name);
                  return (buf);
               }
               line[strlen(line)-1] = '\0';
               *(cptr + i) = strdup (line);
            }
            break;
      }
   }

   fclose (control_file);
   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : fgets_rc
 | COMMENT      : replacement in read_control functions for fgets which
 |                stripps off comments.
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *fgets_rc (char *str, int num, FILE *stream) {
	char *ptr, *ptr2;
	// Four situations: (1) Line is blank, (2) line starts with //,
	//   (3) line has info and contains //, (4) line has info and
	//   does not contain //

	ptr = fgets(str, num, stream);
	if (!ptr) return ptr;

/*
**  A line that starts with "//" is a comment (as far as MMS is concerned).
*/
		if ((str[0] == '/') && (str[1] == '/')) {
			return fgets_rc (str, num, stream);

/*
**  Ignore blank lines
*/
		} else if (strlen (str) <= 1) {
			return fgets_rc (str, num, stream);
/*
** Assume anything else is a data line
*/

		} else if (strstr (str, "//")) {
/*
** comment in data line
*/
			ptr2 = strstr (str, "//");
			ptr2--;
			while (*ptr2 != ' ') ptr2--;
			*(ptr2 + 1) = '\0';
			return ptr;

		} else {
			return ptr;
		}
}

/**7****************** LOCAL FUNCTION DEFINITIONS *********************/

/**8************************** TEST DRIVER ****************************/

