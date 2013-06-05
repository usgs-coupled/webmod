#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

#include "fortran.h"

#if defined(_WIN32) && !defined(_M_AMD64)

#if defined(__cplusplus)
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////
//
// CreateIPhreeqcMMS
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall CREATEIPHREEQCMMS(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall CREATEIPHREEQCMMS_(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall createiphreeqcmms(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall createiphreeqcmms_(void)
{
	return CreateIPhreeqcMMSF();
}

///////////////////////////////////////////////////////////////////////////////
//
// DestroyIPhreeqcMMS
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall DESTROYIPHREEQCMMS(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall DESTROYIPHREEQCMMS_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall destroyiphreeqcmms(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall destroyiphreeqcmms_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}
///////////////////////////////////////////////////////////////////////////////
//
// build_tally_table
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall BUILD_TALLY_TABLE(int *id)
{
	return build_tally_tableF(id);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall BUILD_TALLY_TABLE_(int *id)
{
	return build_tally_tableF(id);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall build_tally_table(int *id)
{
	return build_tally_tableF(id);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall build_tally_table_(int *id)
{
	return build_tally_tableF(id);
}

///////////////////////////////////////////////////////////////////////////////
//
// build_tally_table
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_ROWS_COLUMNS(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_ROWS_COLUMNS_(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall get_tally_table_rows_columns(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall get_tally_table_rows_columns_(int *id, int *rows, int *columns)
{
	return get_tally_table_rows_columnsF(id, rows, columns);
}

///////////////////////////////////////////////////////////////////////////////
//
// get_tally_table_column_heading
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_COLUMN_HEADING(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	return get_tally_table_column_headingF(id, column, type, string, string_length);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_COLUMN_HEADING_(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	return get_tally_table_column_headingF(id, column, type, string, string_length);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall get_tally_table_column_heading(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	return get_tally_table_column_headingF(id, column, type, string, string_length);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall get_tally_table_column_heading_(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	return get_tally_table_column_headingF(id, column, type, string, string_length);
}

///////////////////////////////////////////////////////////////////////////////
//
// get_tally_table_row_heading
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_ROW_HEADING(int *id, int *row, char *string, unsigned int string_length)
{
	return get_tally_table_row_headingF(id, row, string, string_length);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall GET_TALLY_TABLE_ROW_HEADING_(int *id, int *row, char *string, unsigned int string_length)
{
	return get_tally_table_row_headingF(id, row, string, string_length);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall get_tally_table_row_heading(int *id, int *row, char *string, unsigned int string_length)
{
	return get_tally_table_row_headingF(id, row, string, string_length);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall get_tally_table_row_heading_(int *id, int *row, char *string, unsigned int string_length)
{
	return get_tally_table_row_headingF(id, row, string, string_length);
}

///////////////////////////////////////////////////////////////////////////////
//
// RunMix
//
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall RUNMIXF(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	return RunMixF(id, count, solutions, fracs, index_conserv,
		fill_factor, index_rxn, conc_conserv, files_on,
		n_user, rxnmols, tempc, ph, tsec, array,
		row_dim, col_dim);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall RUNMIXF_(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	return RunMixF(id, count, solutions, fracs, index_conserv,
		fill_factor, index_rxn, conc_conserv, files_on,
		n_user, rxnmols, tempc, ph, tsec, array,
		row_dim, col_dim);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall runmixf(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	return RunMixF(id, count, solutions, fracs, index_conserv,
		fill_factor, index_rxn, conc_conserv, files_on,
		n_user, rxnmols, tempc, ph, tsec, array,
		row_dim, col_dim);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall runmixf_(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	return RunMixF(id, count, solutions, fracs, index_conserv,
		fill_factor, index_rxn, conc_conserv, files_on,
		n_user, rxnmols, tempc, ph, tsec, array,
		row_dim, col_dim);
}
///////////////////////////////////////////////////////////////////////////////
//
// MeltPack
//
///////////////////////////////////////////////////////////////////////////////

// /iface:default /names:default
IPQ_DLL_EXPORT int __stdcall MELTPACK(int *id, int *ipack, int *imelt, double *eps, double *ipf, double *fmelt, double *rstd)
{
	return MeltPackF(id, ipack, imelt, eps, ipf, fmelt, rstd);
}

// /iface:default /names:default /assume:underscore
IPQ_DLL_EXPORT int __stdcall MELTPACK_(int *id, int *ipack, int *imelt, double *eps, double *ipf, double *fmelt, double *rstd)
{
	return MeltPackF(id, ipack, imelt, eps, ipf, fmelt, rstd);
}

// /iface:default /names:lowercase
IPQ_DLL_EXPORT int __stdcall meltpack(int *id, int *ipack, int *imelt, double *eps, double *ipf, double *fmelt, double *rstd)
{
	return MeltPackF(id, ipack, imelt, eps, ipf, fmelt, rstd);
}

// /iface:default /names:lowercase /assume:underscore
// /iface:cref /assume:underscore
IPQ_DLL_EXPORT int __stdcall meltpack_(int *id, int *ipack, int *imelt, double *eps, double *ipf, double *fmelt, double *rstd)
{
	return MeltPackF(id, ipack, imelt, eps, ipf, fmelt, rstd);
}
#if defined(__cplusplus)
}
#endif

#endif // _WIN32