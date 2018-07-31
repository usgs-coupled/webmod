/*********************************************************
 * call_modules.c: to replace the one created by 'mbuild',
 * used to call a Fortran version, such as for GSFLOW
 * Creation time: Wed Jan 18 15:52:21 2007
 * Creation time: Thu May 26 10:54:21 2005
 *********************************************************/

#include <stdlib.h>
#include <string.h>
#include "nodes.h"
#include "mms.h"

extern long call_modules_ (char *, ftnlen);

int call_modules(char *arg)

{

 long retval;
 long nmodules;
 long nconi;
 char *mname, *cptr;
 int *xpos, *ypos, *ncon, i;
 ftnlen len;

 if (strncmp (arg, "declare", 7) == 0) {
     nmodules = sizeof(nodes)/sizeof(NODE) - 1;  // subtract one becuase last element is NULL
     nconi = sizeof(con_index)/sizeof(int);

     mname = (char *)malloc (nmodules * sizeof (char) * 21);
     xpos = (int *)malloc (nmodules * sizeof (int));
     ypos = (int *)malloc (nmodules * sizeof (int));
     ncon = (int *)malloc (nmodules * sizeof (int));

     for (i = 0; i < nmodules; i++ ) {
         cptr = strrchr (nodes[i].name, '/');
         if (cptr == NULL) {
             cptr = nodes[i].name;
         } else {
             cptr++; // move to the next char after the /
         }

         strncpy (mname+(i * 20), cptr, 20);
         xpos[i] = nodes[i].x;
         ypos[i] = nodes[i].y;
         ncon[i] = nodes[i].num_connectnions;
     }


    decl_control("module_names", 4, nmodules, mname);
    decl_control("module_x", 1, nmodules, xpos);
    decl_control("module_y", 1, nmodules, ypos);
    decl_control("module_num_con", 1, nmodules, ncon);
    decl_control("module_connections", 1, nconi, con_index);
 }

 len = (ftnlen)strlen(arg);
 retval = call_modules_ (arg, len);
 return((int)retval);

}
