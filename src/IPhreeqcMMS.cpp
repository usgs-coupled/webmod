#include "IPhreeqcMMS.hpp"

#include "Phreeqc.h"                // Phreeqc
#include "global_structures.h"      // OK, STOP
#include "IPhreeqc.h"

void padfstring(char *dest, const char *src, unsigned int len);

/*

Current IPhreeqc
	CreateIPhreeqc                                --    webb.f
		createiphreeqc_                           --    fwrap3.cpp
			CreateIPhreeqcF                       --    fwrap.cpp
				CreateIPhreeqc		              --    IPhreeqcLib.cpp
					IPhreeqcLib::CreateIPhreeqc   --    IPhreeqcLib.cpp


Old IPhreeqcMMS
	build_tally_table                             --    webb.f(test)
		build_tally_table                         --    TallyF.f(FLib)
			build_tally_tableF                    --    phr_cmix.c(IPhreeqcMMS)
				CatchErrors                       --    IPhreeqc.cpp(IPhreeqcMMS-IPhreeqc)
					BuildTallyCallback            --    phr_cmix.c(IPhreeqcMMS)
						build_tally_table         --    tally.c(IPhreeqcMMS-IPhreeqc)

*/

IPhreeqcMMS::IPhreeqcMMS(void)
// COMMENT: {6/6/2012 5:05:49 PM}: pfn_pre(0)
// COMMENT: {6/6/2012 5:05:49 PM}, pfn_post(0)
// COMMENT: {6/6/2012 5:05:49 PM}, cookie(0)
{
}

IPhreeqcMMS::~IPhreeqcMMS(void)
{
}

// COMMENT: {6/6/2012 5:05:45 PM}void IPhreeqcMMS::SetCallbackCookie(void* c)
// COMMENT: {6/6/2012 5:05:45 PM}{
// COMMENT: {6/6/2012 5:05:45 PM}	this->cookie = c;
// COMMENT: {6/6/2012 5:05:45 PM}}
// COMMENT: {6/6/2012 5:05:45 PM}
// COMMENT: {6/6/2012 5:05:45 PM}void IPhreeqcMMS::SetPostRunCallback(PFN_POSTRUN_CALLBACK pfn)
// COMMENT: {6/6/2012 5:05:45 PM}{
// COMMENT: {6/6/2012 5:05:45 PM}	this->pfn_post = pfn;
// COMMENT: {6/6/2012 5:05:45 PM}}
// COMMENT: {6/6/2012 5:05:45 PM}
// COMMENT: {6/6/2012 5:05:45 PM}void IPhreeqcMMS::SetPreRunCallback(PFN_PRERUN_CALLBACK pfn)
// COMMENT: {6/6/2012 5:05:45 PM}{
// COMMENT: {6/6/2012 5:05:45 PM}	this->pfn_pre = pfn;
// COMMENT: {6/6/2012 5:05:45 PM}}

// COMMENT: {6/6/2012 5:01:29 PM}void IPhreeqcMMS::do_run(const char* sz_routine, std::istream* pis, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie)
// COMMENT: {6/6/2012 5:01:29 PM}{
// COMMENT: {6/6/2012 5:01:29 PM}	ASSERT(FALSE);
// COMMENT: {6/6/2012 5:01:29 PM}	IPhreeqc::do_run(sz_routine, pis, this->pfn_pre, this->pfn_post, this->cookie);
// COMMENT: {6/6/2012 5:01:29 PM}	this->pfn_pre  = 0;
// COMMENT: {6/6/2012 5:01:29 PM}	this->pfn_post = 0;
// COMMENT: {6/6/2012 5:01:29 PM}	this->cookie   = 0;
// COMMENT: {6/6/2012 5:01:29 PM}}

// COMMENT: {6/4/2012 5:27:45 PM}int IPhreeqcMMS::CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie)
// COMMENT: {6/4/2012 5:27:45 PM}{
// COMMENT: {6/4/2012 5:27:45 PM}	int rvalue = OK;
// COMMENT: {6/4/2012 5:27:45 PM}	try
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		this->PhreeqcPtr->input_error = 0;
// COMMENT: {6/4/2012 5:27:45 PM}		this->io_error_count = 0;
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}		if (pfn)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			rvalue = pfn(cookie);
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (...)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		const char errmsg[] = "CatchErrors: Unhandled exception occured.\n";
// COMMENT: {6/4/2012 5:27:45 PM}		try
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			this->error_msg(errmsg, STOP); // throws PhreeqcStop
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}		catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	return rvalue;
// COMMENT: {6/4/2012 5:27:45 PM}}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}int IPhreeqcMMS::Catch(IPhreeqcMMSMemFunc pmfn)
// COMMENT: {6/4/2012 5:27:45 PM}{
// COMMENT: {6/4/2012 5:27:45 PM}	int rvalue = OK;
// COMMENT: {6/4/2012 5:27:45 PM}	try
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		this->PhreeqcPtr->input_error = 0;
// COMMENT: {6/4/2012 5:27:45 PM}		this->io_error_count = 0;
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}		if (pmfn)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			//rvalue = CALL_MEMBER_FN((*this), pmfn)();			
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (...)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		const char errmsg[] = "Catch: Unhandled exception occured.\n";
// COMMENT: {6/4/2012 5:27:45 PM}		try
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			this->error_msg(errmsg, STOP); // throws PhreeqcStop
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}		catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	return rvalue;
// COMMENT: {6/4/2012 5:27:45 PM}}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}int IPhreeqcMMS::Catch(PhreeqcMemFunc pmfn)
// COMMENT: {6/4/2012 5:27:45 PM}{
// COMMENT: {6/4/2012 5:27:45 PM}	int rvalue = OK;
// COMMENT: {6/4/2012 5:27:45 PM}	try
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		this->PhreeqcPtr->input_error = 0;
// COMMENT: {6/4/2012 5:27:45 PM}		this->io_error_count = 0;
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}		if (pmfn)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			rvalue = CALL_MEMBER_FN((*this->PhreeqcPtr), pmfn)();
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	catch (...)
// COMMENT: {6/4/2012 5:27:45 PM}	{
// COMMENT: {6/4/2012 5:27:45 PM}		const char errmsg[] = "Catch: Unhandled exception occured.\n";
// COMMENT: {6/4/2012 5:27:45 PM}		try
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			this->error_msg(errmsg, STOP); // throws PhreeqcStop
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}		catch (PhreeqcStop)
// COMMENT: {6/4/2012 5:27:45 PM}		{
// COMMENT: {6/4/2012 5:27:45 PM}			// do nothing
// COMMENT: {6/4/2012 5:27:45 PM}		}
// COMMENT: {6/4/2012 5:27:45 PM}	}
// COMMENT: {6/4/2012 5:27:45 PM}	return rvalue;
// COMMENT: {6/4/2012 5:27:45 PM}}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}
// COMMENT: {6/4/2012 5:27:45 PM}int IPhreeqcMMS::BuildTallyCallback(void *cookie)
// COMMENT: {6/4/2012 5:27:45 PM}{
// COMMENT: {6/4/2012 5:27:45 PM}	IPhreeqcMMS* pthis = (IPhreeqcMMS*)(cookie);
// COMMENT: {6/4/2012 5:27:45 PM}	return pthis->PhreeqcPtr->build_tally_table();
// COMMENT: {6/4/2012 5:27:45 PM}}

int IPhreeqcMMS::CatchBuildTally(void)
{
	int rvalue = OK;
	try
	{
		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		rvalue = this->PhreeqcPtr->build_tally_table();
	}
	catch (IPhreeqcStop)
	{
		rvalue = ERROR;
	}
	catch (...)
	{
		const char errmsg[] = "CatchBuildTally: Unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
		rvalue = ERROR;
	}
	return rvalue;
}

int IPhreeqcMMS::CatchGetTallyRowsColumns(int *rows, int *columns)
{
	int rvalue = OK;
	try
	{
		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		rvalue = this->PhreeqcPtr->get_tally_table_rows_columns(rows, columns);
	}
	catch (IPhreeqcStop)
	{
		rvalue = ERROR;
	}
	catch (...)
	{
		const char errmsg[] = "CatchGetTallyRowsColumns: Unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
		rvalue = ERROR;
	}
	return rvalue;
}

int IPhreeqcMMS::CatchGetTallyColumnHeading(int column, int *type, char *string, unsigned int string_length)
{
	int rvalue = OK;
	try
	{
		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		char* cstr = (char*) malloc((size_t) (string_length * 2));
		if (cstr == NULL)
		{
			this->PhreeqcPtr->malloc_error();
			return ERROR;
		}

		rvalue = this->PhreeqcPtr->get_tally_table_column_heading(column, type, cstr);

		padfstring(string, cstr, string_length);
		free(cstr);
	}
	catch (IPhreeqcStop)
	{
		rvalue = ERROR;
	}
	catch (...)
	{
		const char errmsg[] = "CatchGetTallyColumnHeading: Unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
		rvalue = ERROR;
	}
	return rvalue;
}

int IPhreeqcMMS::CatchGetTallyRowHeading(int row, char *string, unsigned int string_length)
{
	int rvalue = OK;
	try
	{
		this->PhreeqcPtr->input_error = 0;
		this->io_error_count = 0;

		char* cstr = (char*) malloc((size_t) (string_length * 2));
		if (cstr == NULL)
		{
			this->PhreeqcPtr->malloc_error();
			return ERROR;
		}

		rvalue = this->PhreeqcPtr->get_tally_table_row_heading(row, cstr);

		padfstring(string, cstr, string_length);
		free(cstr);
	}
	catch (IPhreeqcStop)
	{
		rvalue = ERROR;
	}
	catch (...)
	{
		const char errmsg[] = "CatchGetTallyColumnHeading: Unhandled exception occured.\n";
		try
		{
			this->error_msg(errmsg, STOP); // throws IPhreeqcStop
		}
		catch (IPhreeqcStop)
		{
			// do nothing
		}
		rvalue = ERROR;
	}
	return rvalue;
}

int IPhreeqcMMS::PreMixCallback(struct MixVars* pvars)
{
// COMMENT: {6/6/2012 5:18:37 PM}	struct MixVars *pvars;
	int i;
	char line[80];

	if (!pvars) return ERROR;

// COMMENT: {6/6/2012 5:14:41 PM}	pvars = (struct MixVars *)cookie;

	if (pvars->fill_factor <= 0.0) 
	{
		char buffer[80];
		sprintf(buffer, "Bad fill_factor in phr_mix for box %d %f.\n", pvars->n_user[Solution], pvars->fill_factor);
		this->PhreeqcPtr->error_msg(buffer, CONTINUE);
	}

	/* MIX */
	if (::AccumulateLine(pvars->id, "MIX") != VR_OK) {
		return ERROR;
	}
	for (i = 0; i < pvars->count; ++i) {
		sprintf(line, "\t%d %g", pvars->solutions[i], pvars->fracs[i]);
		if (::AccumulateLine(pvars->id, line) != VR_OK) {
			return ERROR;
		}
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* COPY SOLUTION n_user[Solution] index_conserv */
	sprintf(line, "COPY SOLUTION %d %d", pvars->n_user[Solution], pvars->index_conserv);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* END */
	if (::AccumulateLine(pvars->id, "END") != VR_OK) {
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
	  if (pvars->files_on) this->OutputAccumulatedLines();

	  this->PhreeqcPtr->zero_tally_table();
	  this->PhreeqcPtr->fill_tally_table(pvars->n_user, pvars->index_conserv, 0); /* initial */

	  return OK;
	}


	/* MIX */
	if (::AccumulateLine(pvars->id, "MIX") != VR_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* Reaction */
	if (pvars->n_user[Reaction] >= 0) {
		if (this->PhreeqcPtr->entity_exists("reaction", pvars->n_user[Reaction])) {
			this->PhreeqcPtr->set_reaction_moles(pvars->n_user[Reaction], pvars->rxnmols * pvars->fill_factor);
			sprintf(line, "USE reaction %d", pvars->n_user[Reaction]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION %d doesn't exist\n", pvars->n_user[Reaction]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}


	/* Exchange */
	if (pvars->n_user[Exchange] >= 0) {
		if (this->PhreeqcPtr->entity_exists("exchange", pvars->n_user[Exchange])) {
			sprintf(line, "USE exchange %d", pvars->n_user[Exchange]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE exchange %d", pvars->n_user[Exchange]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EXCHANGE %d doesn't exist\n", pvars->n_user[Exchange]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Surface */
	if (pvars->n_user[Surface] >= 0) {
		if (this->PhreeqcPtr->entity_exists("surface", pvars->n_user[Surface])) {
			sprintf(line, "USE surface %d", pvars->n_user[Surface]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE surface %d", pvars->n_user[Surface]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SURFACE %d doesn't exist\n", pvars->n_user[Surface]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Gas_phase */
	if (pvars->n_user[Gas_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("gas_phase", pvars->n_user[Gas_phase])) {
			sprintf(line, "USE gas_phase %d", pvars->n_user[Gas_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE gas_phase %d", pvars->n_user[Gas_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "GAS_PHASE %d doesn't exist\n", pvars->n_user[Gas_phase]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Pure_phase */
	if (pvars->n_user[Pure_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("equilibrium_phases", pvars->n_user[Pure_phase])) {
			sprintf(line, "USE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EQUILIBRIUM_PHASES %d doesn't exist\n", pvars->n_user[Pure_phase]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Ss_phase */
	if (pvars->n_user[Ss_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("solid_solution", pvars->n_user[Ss_phase])) {
			sprintf(line, "USE solid_solution %d", pvars->n_user[Ss_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE solid_solution %d", pvars->n_user[Ss_phase]);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SOLID_SOLUTION %d doesn't exist\n", pvars->n_user[Ss_phase]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Kinetics */
	if (pvars->n_user[Kinetics] >= 0) {
		if (this->PhreeqcPtr->entity_exists("kinetics", pvars->n_user[Kinetics])) {
			sprintf(line, "USE kinetics %d", pvars->n_user[Kinetics]);
			this->PhreeqcPtr->set_kinetics_time(pvars->n_user[Kinetics], pvars->tsec);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "KINETICS %d doesn't exist\n", pvars->n_user[Kinetics]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}
	/* Temperature */
	if (pvars->n_user[Temperature] >= 0) {
		if (this->PhreeqcPtr->entity_exists("reaction_temperature", pvars->n_user[Temperature])) {
			sprintf(line, "USE reaction_temperature %d", pvars->n_user[Temperature]);
			this->PhreeqcPtr->set_reaction_temperature(pvars->n_user[Temperature], pvars->tempc);
			if (::AccumulateLine(pvars->id, line) != VR_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION_TEMPERATURE %d doesn't exist\n", pvars->n_user[Temperature]);
			this->PhreeqcPtr->warning_msg(line);
		}
	}

	/* END */
	if (::AccumulateLine(pvars->id, "END") != VR_OK) {
		return ERROR;
	}
#ifdef SKIP
	/* MIX */
	if (::AccumulateLine(pvars->id, "MIX") != VR_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], 1.0 / pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}
#endif
	/* MIX */
	
	sprintf(line, "SOLUTION_MIX %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], 1.0 / pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != VR_OK) {
		return ERROR;
	}

	/* END */
	if (::AccumulateLine(pvars->id, "END") != VR_OK) {
		return ERROR;
	}

	/* OutputLines(); */

	this->PhreeqcPtr->zero_tally_table();
	this->PhreeqcPtr->fill_tally_table(pvars->n_user, pvars->index_conserv, 0); /* initial */

	if (pvars->files_on)
	{
		pvars->orig_out = ::GetOutputFileOn(pvars->id);
		::SetOutputFileOn(pvars->id, 1);
	}

	return OK;
}

int IPhreeqcMMS::PostMixCallback(struct MixVars* pvars)
{
	if (!pvars) return ERROR;

	if (pvars->files_on)
	{
		::SetOutputFileOn(pvars->id, pvars->orig_out);
	}


	this->PhreeqcPtr->fill_tally_table(pvars->n_user, pvars->index_conserv, 1); /* final */
	this->PhreeqcPtr->store_tally_table(pvars->array, pvars->row_dim, pvars->col_dim, pvars->fill_factor);
	return OK;
}

