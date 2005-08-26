#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* printf */

#define EXTERNAL
#include "phreeqc/global.h"
#undef EXTERNAL

/*
 *  Routines
 */
extern void add_all_components(void);
extern int build_tally_table(void);
extern int calc_dummy_kinetic_reaction(struct kinetics *kinetics_ptr);
extern int diff_tally_table(void);
extern int elt_list_to_tally_table(struct buffer *buffer_ptr);
extern int entity_exists (char *name, int n_user);
extern int extend_tally_table(void);
extern int free_tally_table(void);
extern int fill_tally_table(int *n_user, int n_buffer);
extern int get_tally_table_rows_columns(int *rows, int *columns);
extern int get_tally_table_column_heading(int column, int *type, char *string);
extern int get_tally_table_row_heading(int column, char *string);
extern int set_reaction_moles(int n_user, LDBLE moles);
extern int store_tally_table(double *array, int row_dim, int col_dim);
extern int warning_msg (const char *err_str);
extern int zero_tally_table(void);


#include "../include/IPhreeqc.h"

/*******************************
When using GNU gcc/g++/g77
compile fortran with:
	g77 -fno-second-underscore
********************************/
#if defined(__GNUC__)
#define LoadDatabaseF                  loaddatabasef_
#define AccumulateLineF                accumulatelinef_
#define RunF                           runf_
#define RunFileF                       runfilef_
#define GetSelectedOutputValueF        getselectedoutputvaluef_
#define OutputLastErrorF               outputlasterrorf_
#define OutputLinesF                   outputlinesf_
#define GetSelectedOutputRowCountF     getselectedoutputrowcountf_
#define GetSelectedOutputColumnCountF  getselectedoutputcolumncountf_
#define GetSelectedOutputValueF        getselectedoutputvaluef_
#endif

char *
f2cstring(char* fstring, int len)
{
    char *cstr, *str;
    int  i;

    str = fstring;
    for (i = len - 1; i >= 0 && !isgraph((int)str[i]); i--);
    cstr = (char *) malloc((size_t) (i + 2));
    if (!cstr) return 0;
	
    cstr[i + 1] = '\0';
    memcpy(cstr,str,i+1);
    return cstr;
}

void
padfstring(char *dest, char *src, unsigned int len)
{
    unsigned int sofar;

    for (sofar = 0; (sofar < len) && (*src != '\0'); ++sofar)
        *dest++ = *src++;

    while (sofar++ < len)
        *dest++ = ' ';
}

int 
LoadDatabaseF(char* filename, unsigned int filename_length)
{
	int n;
	char* cfilename;

	cfilename = f2cstring(filename, filename_length);
	if (!cfilename) {
		AddError("LoadDatabase: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	n = LoadDatabase(cfilename);

	free(cfilename);

	return n;
}

VRESULT
AccumulateLineF(char *line, unsigned int line_length)
{
	VRESULT n;
	char* cline;
	
	cline = f2cstring(line, line_length);
	if (!cline) {
		AddError("AccumulateLine: Out of memory.\n");
		return VR_OUTOFMEMORY;
	}

	n = AccumulateLine(cline);

	free(cline);
	
	return n;
}

int
RunF(int* output_on, int* error_on, int* log_on, int* selected_output_on)
{
	return Run(*output_on, *error_on, *log_on, *selected_output_on);
}

int
RunFileF(int* output_on, int* error_on, int* log_on, int* selected_output_on, char* filename, unsigned int filename_length)
{
	char* cline;
	
	cline = f2cstring(filename, filename_length);
	if (!cline) {
		AddError("RunFile: Out of memory.\n");
		return (int)VR_OUTOFMEMORY;
	}

	return RunFile(cline, *output_on, *error_on, *log_on, *selected_output_on);
}

int
GetSelectedOutputRowCountF(void)
{
	return GetSelectedOutputRowCount();
}

int
GetSelectedOutputColumnCountF(void)
{
	return GetSelectedOutputColumnCount();
}

VRESULT
GetSelectedOutputValueF(int *row, int *col, int *vtype, double* dvalue, char* svalue, unsigned int svalue_length)
{
	VRESULT result;
	VAR v;
	VarInit(&v);
	result = GetSelectedOutputValue(*row, *col, &v);

	switch (v.type) {
	case TT_EMPTY:
		*vtype = v.type;
		break;
	case TT_ERROR:
		*vtype = v.type;
		break;
	case TT_LONG:
		*vtype = TT_DOUBLE;
		*dvalue = (double)v.lVal;
		break;
	case TT_DOUBLE:
		*vtype = v.type;
		*dvalue = v.dVal;
		break;
	case TT_STRING:
		*vtype = v.type;
		padfstring(svalue, v.sVal, svalue_length);
		break;
	default:
		assert(0);
	}
	VarClear(&v);
	return result;
}

void
OutputLastErrorF(void)
{
	OutputLastError();	
}

void
OutputLinesF(void)
{
	OutputLines();	
}


