#include "fortran.h"

#include "IPhreeqcMMS.h"
#include "IPhreeqcMMS.hpp"
#include "IPhreeqcMMSLib.h"

int
CreateIPhreeqcMMSF(void)
{
	return CreateIPhreeqcMMS();
}

IPQ_RESULT
DestroyIPhreeqcMMSF(int *id)
{
	return DestroyIPhreeqcMMS(*id);
}

int
build_tally_tableF(int *id)
{
	IPhreeqcMMS* IPhreeqcMMSPtr = IPhreeqcMMSLib::GetInstance(*id);
	if (IPhreeqcMMSPtr)
	{
		return IPhreeqcMMSPtr->CatchBuildTally();		
	}
	return IPQ_BADINSTANCE;
}

int
get_tally_table_rows_columnsF(int *id, int *rows, int *columns)
{
	IPhreeqcMMS* IPhreeqcMMSPtr = IPhreeqcMMSLib::GetInstance(*id);
	if (IPhreeqcMMSPtr)
	{
		return IPhreeqcMMSPtr->CatchGetTallyRowsColumns(rows, columns);		
	}
	return IPQ_BADINSTANCE;
}

int
get_tally_table_column_headingF(int *id, int *column, int *type, char *string, unsigned int string_length)
{
	IPhreeqcMMS* IPhreeqcMMSPtr = IPhreeqcMMSLib::GetInstance(*id);
	if (IPhreeqcMMSPtr)
	{
		int n = IPhreeqcMMSPtr->CatchGetTallyColumnHeading(*column - 1, type, string, string_length);
		*type += 1;
		return n;
	}
	return IPQ_BADINSTANCE;
}

int
get_tally_table_row_headingF(int *id, int *row, char *string, unsigned int string_length)
{
	IPhreeqcMMS* IPhreeqcMMSPtr = IPhreeqcMMSLib::GetInstance(*id);
	if (IPhreeqcMMSPtr)
	{
		int n = IPhreeqcMMSPtr->CatchGetTallyRowHeading(*row - 1, string, string_length);
		return n;
	}
	return IPQ_BADINSTANCE;
}

int
RunMixF(int *id, int *count, int *solutions, double *fracs, int *index_conserv,
		double *fill_factor, int *index_rxn, double *conc_conserv, int *files_on,
		int *n_user, double *rxnmols, double *tempc, double *ph, double *tsec, double *array,
		int *row_dim, int *col_dim)
{
	IPhreeqcMMS* IPhreeqcMMSPtr = IPhreeqcMMSLib::GetInstance(*id);
	if (IPhreeqcMMSPtr)
	{
		struct MixVars vars;

		vars.id            = *id;
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

		IPhreeqcMMSPtr->PreMixCallback(&vars);
		int n = IPhreeqcMMSPtr->RunAccumulated();
		IPhreeqcMMSPtr->PostMixCallback(&vars);

		return n;
	}
	return IPQ_BADINSTANCE;
}
