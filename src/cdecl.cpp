#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

#include "fortran.h"

#if defined(__cplusplus)
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////
//
// CreateIPhreeqcMMS
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:default
IPQ_DLL_EXPORT int CREATEIPHREEQCMMS(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:default /names:default /assume:underscore
IPQ_DLL_EXPORT int CREATEIPHREEQCMMS_(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:default /names:lowercase
IPQ_DLL_EXPORT int createiphreeqcmms(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int createiphreeqcmms_(void)
{
	return CreateIPhreeqcMMSF();
}

///////////////////////////////////////////////////////////////////////////////
//
// DestroyIPhreeqcMMS
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:default
IPQ_DLL_EXPORT IPQ_RESULT DESTROYIPHREEQCMMS(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:default /names:default /assume:underscore
IPQ_DLL_EXPORT IPQ_RESULT DESTROYIPHREEQCMMS_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:default /names:lowercase
IPQ_DLL_EXPORT IPQ_RESULT destroyiphreeqcmms(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT IPQ_RESULT destroyiphreeqcmms_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}


///////////////////////////////////////////////////////////////////////////////
//
// build_tally_table
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:default
IPQ_DLL_EXPORT int BUILD_TALLY_TABLE(int *id)
{
	return build_tally_tableF(id);
}

// /iface:default /names:default /assume:underscore
IPQ_DLL_EXPORT int BUILD_TALLY_TABLE_(int *id)
{
	return build_tally_tableF(id);
}

// /iface:default /names:lowercase
IPQ_DLL_EXPORT int build_tally_table(int *id)
{
	return build_tally_tableF(id);
}

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int build_tally_table_(int *id)
{
	return build_tally_tableF(id);
}

///////////////////////////////////////////////////////////////////////////////
//
// build_tally_table
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:default
IPQ_DLL_EXPORT int GET_TALLY_TABLE_ROWS_COLUMNS(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:default /names:default /assume:underscore
IPQ_DLL_EXPORT int GET_TALLY_TABLE_ROWS_COLUMNS_(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:default /names:lowercase
IPQ_DLL_EXPORT int get_tally_table_rows_columns(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int get_tally_table_rows_columns_(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}


///////////////////////////////////////////////////////////////////////////////
//
// get_tally_table_column_heading
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int get_tally_table_column_heading_(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	return get_tally_table_column_headingF(id, column, type, string, string_length);
}

///////////////////////////////////////////////////////////////////////////////
//
// get_tally_table_row_heading
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int get_tally_table_row_heading_(int *id, int *row, char *string, unsigned int string_length)
{
	return get_tally_table_row_headingF(id, row, string, string_length);
}


///////////////////////////////////////////////////////////////////////////////
//
// RunMix
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int runmixf_(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	return RunMixF(id, count, solutions, fracs, index_conserv,
		fill_factor, index_rxn, conc_conserv, files_on,
		n_user, rxnmols, tempc, ph, tsec, array,
		row_dim, col_dim);
}

#if defined(__cplusplus)
}
#endif

