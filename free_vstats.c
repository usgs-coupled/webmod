/*+
 * United States Geological Survey
 * COMMENT  : free linked list for stats variables
 * $Id$
-*/

/**1************************ INCLUDE FILES ****************************/
#define FREE_VSTATS_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : free_vstats
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void free_vstats (void) {
  long nstatVars;
  STAT_LIST_TYPE *curr_stat_list, *prev_stat_list;

  nstatVars = *control_lvar("nstatVars");

  if (nstatVars > 0) {

    curr_stat_list  = Mfirst_stat_list;

    while (curr_stat_list->next != NULL) {
      	prev_stat_list = curr_stat_list;
		curr_stat_list = prev_stat_list->next;
    }

	Mfirst_stat_list = NULL;
  }
}
