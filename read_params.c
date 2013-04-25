/*+
 * United States Geological Survey
 * reads the params data base from a file
 * File name is passed in as an argument
 *
 * $Id$
-*/ 

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
#define READ_PARAMS_C
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <stdlib.h>
#include "mms.h"

static char *READ_param_head (PARAM **, FILE **, char *, char[]);
static char *READ_param_values (PARAM *, FILE *, char []);
static char *rp (char *, int);

int nComments;
char **Comments;

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_params
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_params (char *param_file_name, int index) {
  	static char *foo = NULL;
  	char old[256], *cptr;

	if (foo) {
		strncpy (old, foo, 256);
		free (foo);
		foo = strdup (param_file_name);
	} else {
		strncpy (old, param_file_name, 256);
		foo = strdup (param_file_name);
	}

	cptr = rp (param_file_name, index);

	if (cptr) {
		rp (old, index);

		free (foo);
		foo = strdup (old);
	}

	return (cptr);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_dims
 | COMMENT	:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_dims (char *param_file_name) {
  FILE *param_file;
  DIMEN *dim;
  int dim_size, i, j;

  char line[MAXDATALNLEN], key[MAXDATALNLEN];
  static char buf[256];
  char *endptr;
  char *nch;
  int		done;


/*
* get param name, open file
*/
	if ((param_file = fopen (param_file_name, "r")) == NULL) {
		if (param_file_name) {
			(void)sprintf (buf, "ERROR: cannot open Parameter File: %s", param_file_name);
		} else {
			(void)sprintf (buf, "ERROR: cannot open Parameter File");
		}
		return (buf);
	}

/*
* read in run info string
*/
	if (!fgets (line, MAXDATALNLEN, param_file)) {
		if (param_file != NULL) {
		   fclose (param_file);
		   param_file = NULL;
		}

		(void)sprintf (buf, "ERROR: problems reading info line in Parameter File");
		return (buf);
	}

	if (Mparaminfo) {
		free (Mparaminfo);
	}
	Mparaminfo = strdup (line);

/*
**	See if version number is set
*/
	if (!fgets (line, MAXDATALNLEN, param_file)) {
		if (param_file != NULL) {
		   fclose (param_file);
		   param_file = NULL;
		}

		(void)sprintf (buf, "ERROR: problems reading version number in Parameter File");
		return (buf);
	}

	if (!fgets (line, MAXDATALNLEN, param_file)) {
		if (param_file != NULL) {
		   fclose (param_file);
		   param_file = NULL;
		}
		(void)sprintf (buf, "ERROR: problems reading dimension label in Parameter File");
		return (buf);
	}

/*
 *  Read in comments -- everything between version line and
 *  "** Dimensions **" line is a comment
 */

	Comments = (char **)malloc (1000 * sizeof (char *));
	nComments = 0;

	while (strncmp (line, "** Dimensions **", 16)) {
		if (!fgets (line, MAXDATALNLEN, param_file)) {
		   if (param_file != NULL) {
		      fclose (param_file);
		      param_file = NULL;
		   }
			(void)sprintf (buf, "ERROR: problems skipping comments in Parameter File");
			return (buf);
		}

		if (strncmp (line, "** Dimensions **", 16)) {
			printf ("Comment line = %s\n", line);
			Comments[nComments++] = strdup (line);
		}
	}
	//}

/*
**	Check dimension label
*/
	if (strncmp (line, "** Dimensions **", 16)) {
		if (param_file != NULL) {
		   fclose (param_file);
		   param_file = NULL;
		}
		(void)sprintf (buf, "ERROR: ** Dimensions ** label not found in Parameter File %s.",
		param_file_name);
		return (buf);
	}
  
	if (!fgets (line, MAXDATALNLEN, param_file)) {
		if (param_file != NULL) {
		   fclose (param_file);
		   param_file = NULL;
		}
		(void)sprintf (buf, "ERROR: unexpected end of Parameter File");
		return (buf);
	}

/*
* read in dimensions
*/
	while (strncmp (line, "** Parameters **", 16)) {

		if (strncmp (line, "####", 4)) {
		   if (param_file != NULL) {
		      fclose (param_file);
		      param_file = NULL;
		   }
			(void)sprintf (buf, "ERROR: expecting '####' found %s in Parameter File %s", line, param_file_name);
			return (buf);
		}

/*
**	Read dimension name from parameter file.
*/
		if (fgets (key, MAXDATALNLEN, param_file) == NULL) {
		   if (param_file != NULL) {
		      fclose (param_file);
		      param_file = NULL;
	       }
			(void)sprintf (buf, "ERROR: trying to read dimension name %s in Parameter File %s.", key, param_file_name);
			return (buf);
		}

		key[strlen(key)-1] = '\0';

		dim = dim_addr (key);
		if (dim) {
/*
**	Read dimension size from parameter file.
*/
			if (fgets (line, MAXDATALNLEN, param_file) == NULL) {
		       if (param_file != NULL) {
		          fclose (param_file);
		          param_file = NULL;
	           }
				(void)sprintf (buf, "ERROR: can't read dimension size for %s in Parameter File %s.", key, param_file_name);
				return (buf);
			}

			errno = 0;
			dim_size = strtol(line, &endptr, 10);
			if (errno != 0) {
		       if (param_file != NULL) {
		          fclose (param_file);
		          param_file = NULL;
		       }
				(void)sprintf (buf, "ERROR: size problem with %s in Parameter File %s", key, param_file_name);
				return (buf);
			}

/*
**	If necessary, reset dimension to value read from file.
*/
			if (dim->value != dim_size) {
				reset_dim (dim, dim_size);
			}

/*
* check if there are index names below
*/
			if (fgets (line, MAXDATALNLEN, param_file)) {
				if (strncmp (line, "** Parameters **", 16)) {
					if (dim->names) {
				//        free (dim->names);
						dim->names = NULL;
					}

					if (dim->notes) {
					//        free (dim->notes);
						dim->notes = NULL;
					}

					if (strncmp (line, "####", 4)) {
						dim->names = (char **)calloc (dim_size, sizeof (char *));
						dim->notes = (char **)calloc (dim_size, sizeof (char *));

						done = FALSE;
						i = 0;
						while (!done) {
							if (!strncmp (line, "####", 4)) {
								for (j = i; j < dim_size; j++) {
									dim->names[j] = NULL;
									dim->notes[j] = NULL;
								}
								done = TRUE;

							} else if (line[0] == '@') {
								i--;
								nch = (char *)strchr (line, '\n');
								if (nch) {
									*nch = '\0';
								}
								dim->notes[i] = strdup (&(line[1]));
								fgets (line, MAXDATALNLEN, param_file);
								i++;

							} else {
								nch = (char *)strchr (line, '\n');
								if (nch) {
									*nch = '\0';
								}
								dim->names[i] = strdup (line);
								fgets (line, MAXDATALNLEN, param_file);
								i++;
							}

							if ((i > dim_size) || ((i == dim_size) && (line[0] != '@'))) {
								done = TRUE;
							}
						}
					} else {
						dim->names = NULL;
						dim->files = NULL;
						dim->notes = NULL;
					}
				}
			} else {
		       if (param_file != NULL) {
		          fclose (param_file);  // EOL was returned -- done reading dimensions from this file;
		          param_file = NULL;
		       }
				return (NULL);
			}
		} else {
			(void)fprintf (stderr,"\nWARNING: dimension '%s' is not required; set in parameter file:\n         %s\n", key, param_file_name);
			fgets (line, MAXDATALNLEN, param_file);
			fgets (line, MAXDATALNLEN, param_file);
		}
	}

	if (param_file != NULL) {
	   fclose (param_file);
	   param_file = NULL;
	}
	return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : rp
 | COMMENT	:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *rp (char *param_file_name, int index) {

  FILE *param_file;
  PARAM *param;

  char line[MAXDATALNLEN];
  static char buf[256], *buf_ptr;


/*
* get param name, open file
*/
	if ((param_file = fopen (param_file_name, "r")) == NULL) {
		if (param_file_name)
			(void)sprintf (buf, "ERROR: cannot open Parameter File: %s", param_file_name);
		else
			(void)sprintf (buf, "ERROR: cannot open Parameter File");

		return (buf);
	}

	fgets (line, MAXDATALNLEN, param_file);
	if (index == 0) {  // if index equals zero, than this parameter file has dimension stuff and we need to skip over it.
		while (strncmp (line, "** Parameters **", 16)) {
			if (!fgets (line, MAXDATALNLEN, param_file)) {  // return if hits eol
		       if (param_file != NULL) {
		          fclose (param_file);
		          param_file = NULL;
		       }
			   return (NULL);
			}
		}
		fgets (line, MAXDATALNLEN, param_file);
	}

/*
**	Read in parameters.
*/
	while (!feof (param_file)) {
		buf_ptr = READ_param_head (&param, &param_file, param_file_name, line);
		if (buf_ptr) {
			if (buf_ptr == (char *)-1) {
		      if (param_file != NULL) {
		         fclose (param_file);
		         param_file = NULL;
		      }
			  return (NULL);
			} else {
		       if (param_file != NULL) {
		          fclose (param_file);
		          param_file = NULL;
		       }
			   return (buf_ptr);
			}
		}

		if (param != NULL) {
			buf_ptr = READ_param_values (param, param_file, line);
			if (buf_ptr) {
		        if (param_file != NULL) {
		           fclose (param_file);
		           param_file = NULL;
		        }
				return (buf_ptr);
			}
			updateparam (param->name);
		}
	}
    if (param_file != NULL) {
	   fclose (param_file);
	   param_file = NULL;
    }

	return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : READ_param_head
 | COMMENT		: Read the preliminary stuff for the parameter.  This is
 |                 the stuff between the ####s and where the data actually
 |                 starts.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *READ_param_head (PARAM **param_ptr, FILE **param_file, char *param_file_name, char line[]) {
  char key[MAXDATALNLEN];
  char dimen[MAXDATALNLEN];
  static char buf[256];
  char *temp, *npos, *tempfmt;
  int tempwidth, i, param_size, type;
  
/*
* space fwd to #### header
*/
  while (strncmp (line, "####", 4))
    if (!fgets (line, MAXDATALNLEN, *param_file)) {
		if (*param_file != NULL) {
		   fclose (*param_file);
		   *param_file = NULL;
		}
      return ((char *)-1);
    }
  
/*
* get key, column width and format
*/
  if (fgets (line, MAXDATALNLEN, *param_file) == NULL) {
	  (void)sprintf (buf, "\nERROR: Early end of Parameter File: %s", param_file_name);
    return (buf);
  }

/*
**	get the key
*/
  temp = (char *)strtok(line," ");


  npos = strchr(temp,'\n');
  if (npos) *npos = '\0';

  (void)strcpy(key,temp);
  key[strlen(temp)] = '\0';
	
/*
**	get the column width
*/
  temp = (char *)strtok(NULL," ");
  if (temp)
    tempwidth = atoi(temp);
  else
    tempwidth = 0;

/*
**	get the format
*/
  tempfmt = (char *)strtok(NULL," ");

/*
** markstro -- this check is added so that if there is just a space
**             after the width the parameter will not have a blank
**             format.
*/
  if (tempfmt && (strlen (tempfmt) < 2)) {
     tempfmt = NULL;
  }

/*
**  param is allocated by calls from the modules to declparam.
*/
	*param_ptr = param_addr (key);
	if (*param_ptr) {
	  /*
	  **  Set the read_in flag to true
	  */
		(*param_ptr)->read_in = 1;
/*
* save format and column width
*/
		(*param_ptr)->column_width = tempwidth;
		if (tempfmt) {
			tempfmt[strlen(tempfmt)-1] = '\0';
			if(!(*param_ptr)->format) {
				(*param_ptr)->format = (char *)(malloc(strlen(tempfmt)+1));
			} else {
				(*param_ptr)->format = (char *)(realloc((*param_ptr)->format, strlen(tempfmt) + 1));
			}   
			(void)strcpy((*param_ptr)->format, tempfmt);
		} else {
			(*param_ptr)->format = NULL;
		}
/*
* get number of dimensions
*/
		if(fgets(line, MAXDATALNLEN, *param_file) == NULL) {
			(void)sprintf (buf,"ERROR: reading param number of dimensions for %s in Parameter File %s", key, param_file_name);
			return buf;
		}

		if (isdigit(*line)) {
			if ((*param_ptr)->ndimen != atol(line)) {
				sprintf (buf, "\nERROR: number of dimensions for parameter %s doesn't match parameter declaration.\nParameter File: %s\n", key, param_file_name);
				return buf;
			}

			if((*param_ptr)->ndimen == 0) {
				(void)sprintf (buf, "\nERROR: number of dimensions is 0 for %s in Parameter File %s", key, param_file_name);
				return (buf);
			}
/*
* get number of dimensions if file format supports 2D arrays. Otherwise
* get dimension name.
*/
			for (i = 0; i < (*param_ptr)->ndimen; i++) {
				if(fgets(dimen, MAXDATALNLEN, *param_file) == NULL) {
					(void)sprintf (buf, "\nERROR: number of dimensions is wrong for %s in Parameter File %s", key, param_file_name);
					return (buf);
				}

				dimen[strlen(dimen) - 1] = '\0';
				if (strcmp(dimen, (*param_ptr)->dimen[i]->name)) {
					(void)sprintf (buf, "\nERROR: dimension specification is wrong for %s in Parameter File %s", key, param_file_name);
					return (buf);
				}
			} /* i */

/*
* get param size
*/
			if(fgets(line, MAXDATALNLEN, *param_file) == NULL) {
				(void)sprintf (buf, "ERROR: incorrect parameter size for %s in Parameter File %s", key, param_file_name);
				return (buf);
			}

			if((param_size = atol(line)) == 0) {
				(void)sprintf (buf, "\nERROR: incorrect parameter size for %s in Parameter File %s", key, param_file_name);
				return (buf);
			}

			if(param_size != (*param_ptr)->size) {
				(void)sprintf (buf, "\nERROR: incorrect parameter size for %s in Parameter File %s", key, param_file_name);
				return (buf);
			}

		} else {
			(*param_ptr)->ndimen = 1;
			strncpy(dimen, line, strlen(line));
			dimen[strlen(line)-1] = '\0';

			if (strcmp(dimen, (*param_ptr)->dimen[0]->name)) {
				(void)sprintf (buf, "\nERROR: incorrect dimension specified for parameter %s in Parameter File %s",
				  key, param_file_name);
				return (buf);
			}
			(*param_ptr)->size = getdim(dimen);
			param_size = (*param_ptr)->size;
		}
/*
* get type
*/

		if(fgets(line, MAXDATALNLEN, *param_file) == NULL) {
			(void)sprintf (buf, "\nERROR: incorrect data type specified for parameter %s in Parameter File %s", key, param_file_name);
			return (buf);
		}

		if((type = atol(line)) == 0) {
			sprintf (buf, "\nERROR: incorrect data type specified for parameter %s in Parameter File %s", key, param_file_name);
			return (buf);
		}

		if(type != (*param_ptr)->type) {
			sprintf (buf, "\nERROR: incorrect data type specified for parameter %s in Parameter File %s", key, param_file_name);
			return (buf);
		}
  
	} else {
		(void)printf ("\nWARNING: parameter %s is ignored as it is not required.\n", key);
		(void)printf ("         Read from Parameter File: %s\n", param_file_name);
	}

	return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : READ_param_values
 | COMMENT		: Read the values and comments for the parameter.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *READ_param_values (PARAM *param, FILE *param_file, char line[]) {
	int i, j;
	char *nch;
	int l1, l2, done;
	int	desc_count = 0;
	int repeat_count;
	char delims[] = ",";
	char *result = NULL;
	char *comp_ptr = NULL;
	static char crap[MAXDATALNLEN], crap2[MAXDATALNLEN];
	static char buf[256];
	float foo;
	double d;
	char *endp;
	long l;


/*
**  Space for the values and value_desc are allocated in declparam
*/
	done = FALSE;
	i = 0;
	while (!done) {
		if (!fgets (line, MAXDATALNLEN, param_file)) {
			done = TRUE;

		} else if (!strncmp (line, "####", 4)) {
			done = TRUE;

		} else if (!param) {
			;

		} else if (line[0] == '@') {
			i--;

			nch = (char *)strchr (line, '\n');
			if (nch) *nch = '\0';

			if (desc_count) {
				if (param->value_desc[i]) {
					l1 = strlen (param->value_desc[i]);
					l2 = strlen (line);
					param->value_desc[i] = (char *)realloc
					    (param->value_desc[i],
					    (l1 + l2 + 2) * sizeof (char));
					strcat (param->value_desc[i], "\n");
					strcat (param->value_desc[i], &(line[1]));
				} else {
					param->value_desc[i] = strdup (&(line[1]));
				}
			} else {

				param->value_desc[i] = strdup (&(line[1]));
			}
			i++;

		} else {
			desc_count = 0;
			result = NULL;
			//printf ("READ_param_values: line is %s\n", line);
			strncpy (crap, line, MAXDATALNLEN);
			//printf ("crap is %s\n", crap);

			result = strtok (crap, delims);
			while (result != NULL && !done) {
				//printf ("   READ_param_values: result is |%s|\n", result);

				strncpy (crap2, result, MAXDATALNLEN);
				//printf ("crap2 is %s\n", crap2);
				comp_ptr = strchr (crap2, '*');
				//printf ("comp_ptr is %s\n", comp_ptr);
				if (comp_ptr == NULL){
					repeat_count = 1;
					comp_ptr = crap2;
					//printf ("comp_ptr is %s\n", comp_ptr);
				} else {
					*comp_ptr = '\0';
					repeat_count = atol(crap2);
					comp_ptr++;
					//printf ("comp_ptr is %s\n", comp_ptr);
					foo = (float) atof(comp_ptr);
				}

				for (j = 0; j < repeat_count && !done; j++) {
					if (i < param->size) {
						switch (param->type) {

							case M_STRING:
                                comp_ptr[strlen(comp_ptr)-1] = '\0';
                                *((char **)param->value + i) = strdup (comp_ptr);
                                i++;

								//if (comp_ptr != endp && *endp == '\n') {

								//} else {
								//	sprintf (buf, "There is a parameter format error. Parameter name: %s Index = %d\n   The data type should be a character string or there could be white spaces after the values on the line.", param->name, (i+1));
									//printf ("%s", buf);
								//	return (buf);
								//}

								//((double *)(param->value))[i++] = atof(comp_ptr);
								break;

							case M_DOUBLE:
								d = strtod(comp_ptr, &endp);
								if (comp_ptr != endp && *endp == '\n') {
									((double *)(param->value))[i++] = d;
								} else {
									sprintf (buf, "\nERROR: parameter format error. Parameter name: %s Index = %d\n   The data type should be a double precision float or there could be white spaces after the values on the line.", param->name, (i+1));
									return (buf);
								}
								break;

							case M_FLOAT:
								d = strtod(comp_ptr, &endp);
								if (comp_ptr != endp && *endp == '\n') {
									((float *)(param->value))[i++] = (float)d;
								} else {
									sprintf (buf, "\nERROR: parameter format error. Parameter name: %s Index = %d\n   The data type should be a float or there could be white spaces after the values on the line.", param->name, (i+1));
									return (buf);
								}
								break;

							case M_LONG:
								l = strtol(comp_ptr, &endp, 0);
								if (comp_ptr != endp && *endp == '\n') {
									((int *)(param->value))[i++] = (int)l;
								} else {
									sprintf (buf, "\nERROR: parameter format error. Parameter name: %s Index = %d\n   The data type should be an integer or there could be white spaces after the values on the line.", param->name, (i+1));
									return (buf);
								}
								break;
						} // switch
				 
					} else { // if (i < param->size)
						done = TRUE;
						i++;
					} // if (i < param->size)
				}
				result = strtok(NULL, delims);
			} // while
		}
	}

	if (i < param->size) {
		sprintf (buf, "\nERROR: too few values read for paramter %s in Parameter File", param->name);
		return (buf);
	} else if (i > param->size) {
		sprintf (buf, "\nERROR: too many values read for paramter %s in Parameter File", param->name);
		return (buf);
	}
	return (NULL);
}
