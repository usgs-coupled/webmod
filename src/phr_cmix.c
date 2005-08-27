#include <ctype.h>   /* isgraph */
#include <stdlib.h>  /* malloc */
#include <memory.h>  /* memcpy */
#include <assert.h>  /* assert */
#include <stdio.h>   /* printf */

#define EXTERNAL
#include "../IPhreeqc/src/phreeqc/global.h"
#undef EXTERNAL

#include "../IPhreeqc/include/IPhreeqc.h"

struct MixVars {
  int      count;
  int     *solutions;
  double  *fracs;
  int      index_conserv;
  double   fill_factor;
  int     *n_user;
  double   rxnmols;
  double   tempc;
  double   tsec;
  double  *array;
  int	   row_dim;
  int	   col_dim;
  int      files_on;
};

struct TallyVars {
	int  *rows;
	int  *columns;
	int  row;
	int  column;
	int  *type;
	char *string;
	size_t string_length;
};


extern void add_all_components(void);
extern int build_tally_table(void);
extern int calc_dummy_kinetic_reaction(struct kinetics *kinetics_ptr);
extern int diff_tally_table(void);
extern int entity_exists (char *name, int n_user);
extern int error_msg (const char *err_str, const int stop);
extern int extend_tally_table(void);
extern int free_tally_table(void);
extern int fill_tally_table(int *n_user, int index_conservative, int n_buffer);
extern int get_tally_table_rows_columns(int *rows, int *columns);
extern int get_tally_table_column_heading(int column, int *type, char *string, size_t string_length);
extern int get_tally_table_row_heading(int column, char *string, size_t string_length);
extern void malloc_error(void);
extern int set_reaction_moles(int n_user, LDBLE moles);
extern int set_reaction_temperature(int n_user, LDBLE tc);
extern int set_kinetics_time(int n_user, LDBLE step);
extern int store_tally_table(double *array, int row_dim, int col_dim, double fill_factor);
extern int warning_msg (const char *err_str);
extern int zero_tally_table(void);


extern char *f2cstring(char* fstring, int len);
extern void padfstring(char *dest, char *src, unsigned int len);


/*******************************
When using GNU gcc/g++/g77
compile fortran with:
	g77 -fno-second-underscore
********************************/
#if defined(__GNUC__)
#define build_tally_tableF              build_tally_tablef_
#define get_tally_table_column_headingF get_tally_table_column_headingf_
#define get_tally_table_row_headingF    get_tally_table_row_headingf_
#define get_tally_table_rows_columnsF   get_tally_table_rows_columnsf_
#define RunMixF                         runmixf_
#endif

int BuildTallyCallback(void* cookie)
{
	return build_tally_table();
}

int
build_tally_tableF(void)
{
	return CatchErrors(BuildTallyCallback, NULL);
}

int RowsColsCallback(void* cookie)
{
	struct TallyVars *pvars;

	if (!cookie) return ERROR;

	pvars = (struct TallyVars *)cookie;

	return get_tally_table_rows_columns(pvars->rows, pvars->columns);
}

int
get_tally_table_rows_columnsF(int *rows, int *columns)
{
	struct TallyVars vars;
	vars.rows    = rows;
	vars.columns = columns;
	return CatchErrors(RowsColsCallback, &vars);
}

int ColumnHeadingCallback(void* cookie)
{
	struct TallyVars *pvars;

	if (!cookie) return ERROR;

	pvars = (struct TallyVars *)cookie;

	return get_tally_table_column_heading(pvars->column, pvars->type, pvars->string, pvars->string_length);
}

int 
get_tally_table_column_headingF(int* column, int *type, char *string, unsigned int string_length)
{
	struct TallyVars vars;
	char* cstr;
	int rval;

    cstr = (char *) malloc((size_t) (string_length + 2));
    if (cstr == NULL) {
		malloc_error();
		return ERROR;
	}

	vars.column        = *column;
	vars.type          = type;
	vars.string        = cstr;
	vars.string_length = string_length;

	rval = CatchErrors(ColumnHeadingCallback, &vars);

	padfstring(string, cstr, string_length);
	free(cstr);
	return rval;
}

int RowHeadingCallback(void* cookie)
{
	struct TallyVars *pvars;

	if (!cookie) return ERROR;

	pvars = (struct TallyVars *)cookie;

	return get_tally_table_row_heading(pvars->row, pvars->string, pvars->string_length);
}

int
get_tally_table_row_headingF(int* row, char *string, unsigned int string_length)
{
	struct TallyVars vars;
	char* cstr;
	int rval;

    cstr = (char *) malloc((size_t) (string_length + 2));
    if (cstr == NULL) {
		malloc_error();
		return ERROR;
	}

	vars.row           = *row;
	vars.string        = cstr;
	vars.string_length = string_length;

	rval = CatchErrors(RowHeadingCallback, &vars);

	padfstring(string, cstr, string_length);
	free(cstr);
	return rval;
}

int PreMixCallback(void* cookie)
{
	struct MixVars *pvars;
	int i;
	char line[80];

	if (!cookie) return ERROR;

	pvars = (struct MixVars *)cookie;

	if (pvars->fill_factor <= 0.0) 
	{
		char buffer[80];
		sprintf(buffer, "Bad fill_factor in phr_mix for box %d %f.\n", pvars->n_user[Solution], pvars->fill_factor);
		error_msg(buffer, CONTINUE);
	}

	/* MIX */
	if (AccumulateLine("MIX") != VR_OK) {
		return ERROR;
	}
	for (i = 0; i < pvars->count; ++i) {
		sprintf(line, "\t%d %g", pvars->solutions[i], pvars->fracs[i]);
		if (AccumulateLine(line) != VR_OK) {
			return ERROR;
		}
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}

	/* COPY SOLUTION n_user[Solution] index_conserv */
	sprintf(line, "COPY SOLUTION %d %d", pvars->n_user[Solution], pvars->index_conserv);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}

	/* END */
	if (AccumulateLine("END") != VR_OK) {
		return ERROR;
	}


	if (pvars->n_user[Reaction]    < 0
		&&
		pvars->n_user[Exchange]    < 0
		&&
		pvars->n_user[Surface]     < 0
		&&
		pvars->n_user[Gas_phase]   < 0
		&&
		pvars->n_user[Pure_phase]  < 0
		&&
		pvars->n_user[Ss_phase]    < 0
		&&
		pvars->n_user[Kinetics]    < 0
		&&
		pvars->n_user[Temperature] < 0)
	{
	  if(pvars->files_on)	OutputLines();

	  zero_tally_table();
	  fill_tally_table(pvars->n_user, pvars->index_conserv, 0); /* initial */

	  return OK;
	}


	/* MIX */
	if (AccumulateLine("MIX") != VR_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], pvars->fill_factor);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}

	/* Reaction */
	if (pvars->n_user[Reaction] >= 0) {
		if (entity_exists("reaction", pvars->n_user[Reaction])) {
			set_reaction_moles(pvars->n_user[Reaction], pvars->rxnmols * pvars->fill_factor);
			sprintf(line, "USE reaction %d", pvars->n_user[Reaction]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION %d doesn't exist\n", pvars->n_user[Reaction]);
			warning_msg(line);
		}
	}


	/* Exchange */
	if (pvars->n_user[Exchange] >= 0) {
		if (entity_exists("exchange", pvars->n_user[Exchange])) {
			sprintf(line, "USE exchange %d", pvars->n_user[Exchange]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE exchange %d", pvars->n_user[Exchange]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EXCHANGE %d doesn't exist\n", pvars->n_user[Exchange]);
			warning_msg(line);
		}
	}
	/* Surface */
	if (pvars->n_user[Surface] >= 0) {
		if (entity_exists("surface", pvars->n_user[Surface])) {
			sprintf(line, "USE surface %d", pvars->n_user[Surface]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE surface %d", pvars->n_user[Surface]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SURFACE %d doesn't exist\n", pvars->n_user[Surface]);
			warning_msg(line);
		}
	}
	/* Gas_phase */
	if (pvars->n_user[Gas_phase] >= 0) {
		if (entity_exists("gas_phase", pvars->n_user[Gas_phase])) {
			sprintf(line, "USE gas_phase %d", pvars->n_user[Gas_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE gas_phase %d", pvars->n_user[Gas_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "GAS_PHASE %d doesn't exist\n", pvars->n_user[Gas_phase]);
			warning_msg(line);
		}
	}
	/* Pure_phase */
	if (pvars->n_user[Pure_phase] >= 0) {
		if (entity_exists("equilibrium_phases", pvars->n_user[Pure_phase])) {
			sprintf(line, "USE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EQUILIBRIUM_PHASES %d doesn't exist\n", pvars->n_user[Pure_phase]);
			warning_msg(line);
		}
	}
	/* Ss_phase */
	if (pvars->n_user[Ss_phase] >= 0) {
		if (entity_exists("solid_solution", pvars->n_user[Ss_phase])) {
			sprintf(line, "USE solid_solution %d", pvars->n_user[Ss_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE solid_solution %d", pvars->n_user[Ss_phase]);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SOLID_SOLUTION %d doesn't exist\n", pvars->n_user[Ss_phase]);
			warning_msg(line);
		}
	}
	/* Kinetics */
	if (pvars->n_user[Kinetics] >= 0) {
		if (entity_exists("kinetics", pvars->n_user[Kinetics])) {
			sprintf(line, "USE kinetics %d", pvars->n_user[Kinetics]);
			set_kinetics_time(pvars->n_user[Kinetics], pvars->tsec);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "KINETICS %d doesn't exist\n", pvars->n_user[Kinetics]);
			warning_msg(line);
		}
	}
	/* Temperature */
	if (pvars->n_user[Temperature] >= 0) {
		if (entity_exists("reaction_temperature", pvars->n_user[Temperature])) {
			sprintf(line, "USE reaction_temperature %d", pvars->n_user[Temperature]);
			set_reaction_temperature(pvars->n_user[Temperature], pvars->tempc);
			if (AccumulateLine(line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION_TEMPERATURE %d doesn't exist\n", pvars->n_user[Temperature]);
			warning_msg(line);
		}
	}

	/* END */
	if (AccumulateLine("END") != VR_OK) {
		return ERROR;
	}

	/* MIX */
	if (AccumulateLine("MIX") != VR_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], 1.0 / pvars->fill_factor);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (AccumulateLine(line) != VR_OK) {
		return ERROR;
	}


	/* OutputLines(); */

	zero_tally_table();
	fill_tally_table(pvars->n_user, pvars->index_conserv, 0); /* initial */

	return OK;
}

int PostMixCallback(void* cookie)
{
	struct MixVars *pvars;

	if (!cookie) return ERROR;

	pvars = (struct MixVars *)cookie;

	fill_tally_table(pvars->n_user, pvars->index_conserv, 1); /* final */
	store_tally_table(pvars->array, pvars->row_dim, pvars->col_dim, pvars->fill_factor);
	return OK;
}

/**
int
RunMixF(int *count, int *solutions, double *fracs, int *dest, double *conc, int *files_on,
	   int *n_user, double *rxnmols, double *tempc, double *tsec, double *array, int *row_dim, int *col_dim)


int
RunMixF(int *count, int *solutions, double *fracs, int *dest,
       /{ double *fill_factor }/, /{ int *index_rxn }/ double *conc, int *files_on,
	   int *n_user, double *rxnmols, double *tempc, double *tsec, double *array,
	   int *row_dim, int *col_dim)
**/

int
RunMixF(int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	struct MixVars vars;
	int result;

	vars.count         = *count;
	vars.solutions     = solutions;
	vars.fracs         = fracs;
	vars.index_conserv = *index_conserv;
	vars.fill_factor   = *fill_factor;
	vars.n_user        = n_user;
	vars.rxnmols       = *rxnmols;	
	vars.tempc         = *tempc;
	vars.tsec          = *tsec;
	vars.array         = array;
	vars.row_dim       = *row_dim;
	vars.col_dim       = *col_dim;	
	vars.files_on      = *files_on;	

	n_user[Solution] = *index_rxn;

	result = RunWithCallback(PreMixCallback, PostMixCallback, &vars, *files_on, *files_on, *files_on, *files_on);

	return result;
}
