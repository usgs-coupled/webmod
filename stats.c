/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : stats
 * COMMENT  : statistical analysis postprocessor
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define STATS_C
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "mms.h"

#define MAXCELLS 100

/**2************************* LOCAL MACROS ****************************/
#define SQR(A) (A) * (A)
#define CUBE(A) (A) * (A) * (A)

/**3************************ LOCAL TYPEDEFS ***************************/
  typedef struct {
    char varName[30];        /* variable name                  */
    char * elem_number;      /* element number                 */
    float Sx;                /* Sum of x                       */
    float Sx2;               /* Sum of x*x                     */
    float Sx3;               /* Sum of x*x*x                   */
    float mx;		     /* mean of x                      */
    float sdev;		     /* standard deviation             */
    float skew;		     /* skewness                       */
    float min;		     /* minimum                        */
    float max;		     /* maximum                        */
    float histmin;	     /* histogram minimum              */
    float histmax;	     /* histogram maximum              */
    int   ncells;	     /* number of cell in the histogram*/
    float width;             /* width of the histogram cell    */
    float histog[MAXCELLS];  /* histogram                      */
  }STATS;

/*--------------------------------------------------------------------*\
 | FUNCTION     : stats
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : int
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int stats (void) {
     
  int     nvars;                 /* number of variables in the file*/
  STATS   st[MAXSTATVARS];            /* array of statistics structures */
  FILE    *statvar_file;
  char    path[MAXDATALNLEN];
  char    line[MAXDATALNLEN];
  int     i,nvals;
  int     recNo;
  int     year;
  int     month;
  int     day;
  int     hour;
  int     minute;
  int     second;
  char    elem_number[MAXDATALNLEN];
  double  x[MAXSTATVARS];
  float   squared;

  /*
   * Open statvar file, and store number of variables and variable names 
   */

  if (*((char **) control_var("stat_var_file")) == NULL) {
    (void)fprintf(stderr, "ERROR - stats");
    (void)fprintf(stderr, "Check Control File, stat_var_file missing.\n");
    return(1);
  }

  (void)sprintf(path, "%s", *((char **) control_var("stat_var_file")));

  if ((statvar_file = fopen(path, "r")) == NULL) {
      (void)fprintf(stderr, "ERROR - stats");
      (void)fprintf(stderr, "Could not open statvar file '%s'\n",
	      path);
      return(1);
    }
  
  fscanf(statvar_file,"%d",&nvars);

	if (!nvars) {
		return(0);
	}


  for (i=0;i<nvars;i++) {
      memset (&st[i], 0, sizeof(STATS));
      st[i].min = 1e30;
      st[i].max = -1e30;
      fscanf(statvar_file,"%s %s", st[i].varName,
	     elem_number);
      st[i].elem_number = (char *)malloc(strlen(elem_number) + 1);
      (void)strcpy(st[i].elem_number, elem_number);

    }

  nvals = 0;

  while (EOF !=
	 fscanf(statvar_file, "%d %d %d %d %d %d %d",
		&recNo,&year,&month,&day,&hour,&minute,&second))
    {
      nvals++;
      for (i=0;i<nvars;i++) 
	{
	  fscanf(statvar_file, "%lf", &x[i]); 
	  st[i].Sx += x[i];
	  st[i].Sx2 += SQR(x[i]);
	  st[i].Sx3 += CUBE(x[i]);
	  st[i].min = MIN(st[i].min,x[i]);
	  st[i].max = MAX(st[i].max,x[i]);
	}
    }

  for (i=0;i<nvars;i++)
    {
      if (nvals > 1) {
	st[i].mx   = st[i].Sx/nvals;
	
	squared = (st[i].Sx2 - nvals*SQR(st[i].mx)) / (nvals-1);
	
	if (squared < 0.0)
	  squared = 0.0;
	
	st[i].sdev =  (float)sqrt(squared);
	
	if (st[i].sdev > 0.0) {
	  st[i].skew = ((st[i].Sx3/nvals -
			 3.0*st[i].mx*st[i].Sx2/nvals+2.0*CUBE(st[i].mx))
			/(CUBE(st[i].sdev)));
	}
	
	if(nvals) 
/* DANGER
	  st[i].ncells = MIN(1 + 3.3 * log10((double)nvals),MAXCELLS);
*/
	  st[i].ncells = MAXCELLS;
	else
	  {
	    /*
	     ** NOTE: No values exist for histogram
	     */
	    st[i].ncells = 0;
	  }
	
	st[i].histmin = st[i].min;
	st[i].histmax = st[i].max;
	
	if (st[i].ncells)
	  st[i].width = ((st[i].histmax - st[i].histmin)/st[i].ncells);
	else
	  st[i].width = 0.0;
	
      }
    }
  /*
   * rewind the statvar file
   */

  fseek(statvar_file, 0L, 0);

  /*
   * space fwd to data
   */

  for (i = 0; i < nvars+1; i++) {
    if (fgets(line, MAXDATALNLEN, statvar_file) == NULL) {
      (void)fprintf(stderr, "ERROR - stats.\n");
      (void)fprintf(stderr, "Reading statvar file for histogram comps.\n");
      perror(path);
      return(1);
    }
  }

  /*
   * re-read the data to load the histograms
   */

  while (EOF !=
	 fscanf(statvar_file,"%d %d %d %d %d %d %d",
		&recNo,&year,&month,&day,&hour,&minute,&second))
    {
      for (i=0;i<nvars;i++) 
	{
	  fscanf(statvar_file,"%lf",&x[i]); 
	  (st[i].histog[(int)((x[i]-st[i].histmin)/st[i].width)])++;
	}
    }

  /*
   * close statvar file
   */

  fclose(statvar_file);

  return(0);

}
