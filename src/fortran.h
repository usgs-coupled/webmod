#pragma once

#include "IPhreeqc.h"                      // IPQ_RESULT

int CreateIPhreeqcMMSF(void);

IPQ_RESULT DestroyIPhreeqcMMSF(int *id);

int build_tally_tableF(int *id);

int get_tally_table_rows_columnsF(int *id, int *rows, int *columns);

int get_tally_table_column_headingF(int *id, int *column, int *type, char *string, unsigned int string_length);

int get_tally_table_row_headingF(int *id, int *row, char *string, unsigned int string_length);

int RunMixF(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
			double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
			int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
			int *row_dim, int *col_dim);



