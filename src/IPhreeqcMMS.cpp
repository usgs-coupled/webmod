#include "IPhreeqcMMS.hpp"

#include "Phreeqc.h"                // Phreeqc
#include "global_structures.h"      // OK, STOP
#include "IPhreeqc.h"
#include "Solution.h"

void padfstring(char *dest, const char *src, int *len);

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

		int strl = (int) string_length;
		padfstring(string, cstr, &strl);
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
		
		int strl = (int) string_length;
		padfstring(string, cstr, &strl);
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
	bool evap_or_sublimation;

	if (!pvars) return ERROR;

// COMMENT: {6/6/2012 5:14:41 PM}	pvars = (struct MixVars *)cookie;

	if (pvars->fill_factor <= 0.0) 
	{
		char buffer[80];
		sprintf(buffer, "Bad fill_factor in phr_mix for box %d %f.\n", pvars->n_user[Solution], pvars->fill_factor);
		this->error_msg(buffer, CONTINUE);
	}

	double pos_mix_h2o = 0.0;
	double pos_h2o_frac_sum = 0;
	double neg_mix_h2o = 0.0;
	double neg_factor = 1.0;

	evap_or_sublimation = false;
	for (i = 0; i < pvars->count; ++i) 
	{
		if (pvars->solutions[i] == 1 && pvars->fracs[i] < 0)
		{
			evap_or_sublimation = true;
			break;
		}
	}
	if (evap_or_sublimation )
	{
		for (i = 0; i < pvars->count; ++i) 
		{
			std::map<int, cxxSolution>::const_iterator it = this->PhreeqcPtr->Rxn_solution_map.find(pvars->solutions[i]);
			if (it != this->PhreeqcPtr->Rxn_solution_map.end())
			{
				// need to find mass of water for pvars->solutions[i]
				double h2o = it->second.Get_mass_water();
				if (pvars->fracs[i] >= 0.0)
				{
					pos_mix_h2o += pvars->fracs[i] * h2o;
					pos_h2o_frac_sum += pvars->fracs[i];
				}
				else
				{
					neg_mix_h2o = h2o;
				}
			}
			else
			{
				std::cerr << "Did not find solution " << pvars->solutions[i] << std::endl;
				return ERROR;
			}
		}
		if (neg_mix_h2o > 0.0)
		{
			double avg_h2o = pos_mix_h2o / pos_h2o_frac_sum;
			neg_factor = fabs(avg_h2o / neg_mix_h2o);
		}
	}

	/* MIX */
	if (::AccumulateLine(pvars->id, "MIX") != IPQ_OK) {
		return ERROR;
	}
	for (i = 0; i < pvars->count; ++i) {
		double factor;
		if (evap_or_sublimation)
		{
			factor = pvars->fracs[i] < 0.0 ? neg_factor : 1.0;
		}
		else
		{
			factor = 1.0;
		}
		sprintf(line, "\t%d %g", pvars->solutions[i], pvars->fracs[i]*factor);
		if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
			return ERROR;
		}
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

	/* COPY SOLUTION n_user[Solution] index_conserv */
	sprintf(line, "COPY SOLUTION %d %d", pvars->n_user[Solution], pvars->index_conserv);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

	/* END */
	if (::AccumulateLine(pvars->id, "END") != IPQ_OK) {
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
	if (::AccumulateLine(pvars->id, "MIX") != IPQ_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

	/* Reaction */
	if (pvars->n_user[Reaction] >= 0) {
		if (this->PhreeqcPtr->entity_exists("reaction", pvars->n_user[Reaction])) {
			this->PhreeqcPtr->set_reaction_moles(pvars->n_user[Reaction], pvars->rxnmols * pvars->fill_factor);
			sprintf(line, "USE reaction %d", pvars->n_user[Reaction]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION %d doesn't exist\n", pvars->n_user[Reaction]);
			this->warning_msg(line);
		}
	}


	/* Exchange */
	if (pvars->n_user[Exchange] >= 0) {
		if (this->PhreeqcPtr->entity_exists("exchange", pvars->n_user[Exchange])) {
			sprintf(line, "USE exchange %d", pvars->n_user[Exchange]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE exchange %d", pvars->n_user[Exchange]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EXCHANGE %d doesn't exist\n", pvars->n_user[Exchange]);
			this->warning_msg(line);
		}
	}
	/* Surface */
	if (pvars->n_user[Surface] >= 0) {
		if (this->PhreeqcPtr->entity_exists("surface", pvars->n_user[Surface])) {
			sprintf(line, "USE surface %d", pvars->n_user[Surface]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE surface %d", pvars->n_user[Surface]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SURFACE %d doesn't exist\n", pvars->n_user[Surface]);
			this->warning_msg(line);
		}
	}
	/* Gas_phase */
	if (pvars->n_user[Gas_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("gas_phase", pvars->n_user[Gas_phase])) {
			sprintf(line, "USE gas_phase %d", pvars->n_user[Gas_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE gas_phase %d", pvars->n_user[Gas_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "GAS_PHASE %d doesn't exist\n", pvars->n_user[Gas_phase]);
			this->warning_msg(line);
		}
	}
	/* Pure_phase */
	if (pvars->n_user[Pure_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("equilibrium_phases", pvars->n_user[Pure_phase])) {
			sprintf(line, "USE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE equilibrium_phases %d", pvars->n_user[Pure_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "EQUILIBRIUM_PHASES %d doesn't exist\n", pvars->n_user[Pure_phase]);
			this->warning_msg(line);
		}
	}
	/* Ss_phase */
	if (pvars->n_user[Ss_phase] >= 0) {
		if (this->PhreeqcPtr->entity_exists("solid_solution", pvars->n_user[Ss_phase])) {
			sprintf(line, "USE solid_solution %d", pvars->n_user[Ss_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
			sprintf(line, "SAVE solid_solution %d", pvars->n_user[Ss_phase]);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "SOLID_SOLUTION %d doesn't exist\n", pvars->n_user[Ss_phase]);
			this->warning_msg(line);
		}
	}
	/* Kinetics */
	if (pvars->n_user[Kinetics] >= 0) {
		if (this->PhreeqcPtr->entity_exists("kinetics", pvars->n_user[Kinetics])) {
			sprintf(line, "USE kinetics %d", pvars->n_user[Kinetics]);
			this->PhreeqcPtr->set_kinetics_time(pvars->n_user[Kinetics], pvars->tsec);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "KINETICS %d doesn't exist\n", pvars->n_user[Kinetics]);
			this->warning_msg(line);
		}
	}
	/* Temperature */
	if (pvars->n_user[Temperature] >= 0) {
		if (this->PhreeqcPtr->entity_exists("reaction_temperature", pvars->n_user[Temperature])) {
			sprintf(line, "USE reaction_temperature %d", pvars->n_user[Temperature]);
			this->PhreeqcPtr->set_reaction_temperature(pvars->n_user[Temperature], pvars->tempc);
			if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
				return ERROR;
			}
		}
		else {
			sprintf(line, "REACTION_TEMPERATURE %d doesn't exist\n", pvars->n_user[Temperature]);
			this->warning_msg(line);
		}
	}

	/* END */
	if (::AccumulateLine(pvars->id, "END") != IPQ_OK) {
		return ERROR;
	}

	/* MIX */
	if (::AccumulateLine(pvars->id, "MIX") != IPQ_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], 1.0 / pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

	/* SAVE solution n_user[Solution] */
	sprintf(line, "SAVE solution %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}

#ifdef SKIP
	/* MIX */
	
	sprintf(line, "SOLUTION_MIX %d", pvars->n_user[Solution]);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}
	sprintf(line, "\t%d %g", pvars->n_user[Solution], 1.0 / pvars->fill_factor);
	if (::AccumulateLine(pvars->id, line) != IPQ_OK) {
		return ERROR;
	}
#endif
	/* END */
	if (::AccumulateLine(pvars->id, "END") != IPQ_OK) {
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

int IPhreeqcMMS::Melt_pack(int ipack, int imelt, double eps, double ipf, double fmelt, double rstd)
{
	std::ostringstream strm;
	strm << "USE solution " << ipack << "\n";
	strm << "REACTION 1" << "\n";
	strm << "SELECTED_OUTPUT" << "\n";
	strm << "-reset false; -high\n";
	strm << "USER_PUNCH; -start\n";
	strm << "10 eps = " << -eps << "\n"; 
	strm << "20 rstd = " << rstd << "\n";
	strm << "30 fm = " << fmelt << "\n";
	strm << "40 ipf = " << ipf << "\n";
	strm << "50 fracp = 1/(1+ipf*fm/(1-fm))\n";
	strm << "60 fracm = 1-fracp\n";
	strm << "70 xfer_p = TOT(\"water\")*((1-fm)-fracp)/(0.001*GFW(\"H2O\"))\n";
	strm << "80 xfer_m = TOT(\"water\")*(fm-fracm)/(0.001*GFW(\"H2O\"))\n";
	strm << "90 PUNCH \"TITLE Makes melt water\", EOL$\n";
	strm << "100 PUNCH \"MIX\", EOL$\n";
	strm << "110 PUNCH " << "\"" << ipack <<"\"" <<" fracm, EOL$\n";
	strm << "120 PUNCH \"REACTION 1\", EOL$\n";
	strm << "130 PUNCH \"H2O \", xfer_m, EOL$\n";
	strm << "140 o_p = fracp*TOTMOLE(\"O\")+xfer_p\n";
	strm << "150 o_m = fracm*TOTMOLE(\"O\")+xfer_m\n";
	strm << "160 o18_m = (TOTMOLE(\"[18O]\")/o_p-eps/1000*rstd)/(1/o_m+1/o_p)\n";
	strm << "170 o18_p = (TOTMOLE(\"[18O]\")/o_m+eps/1000*rstd)/(1/o_p+1/o_m)\n";
	strm << "180 diff = o18_m-fracm*TOTMOLE(\"[18O]\")\n";
	strm << "190 PUNCH \"H2[18O] \", diff, EOL$\n";
	strm << "200 PUNCH 1, EOL$\n";
	strm << "210 PUNCH \"SAVE solution \" " << "\"" << imelt <<"\"" << ", EOL$\n";
	//strm << "220 PUNCH \"END\", EOL$\n";
	strm << "-end\n";
	//strm << "SAVE solution " << imelt  << "\n";
	strm << "END\n";

	this->SetSelectedOutputStringOn(true);
	//this->SetOutputStringOn(true);	
	//std::cerr << strm.str().c_str() << std::endl;
	this->RunString(strm.str().c_str());
	std::string selout = this->GetSelectedOutputString();
	//std::cerr << this->GetSelectedOutputString() << std::endl;
	//std::cerr << this->GetOutputString() << std::endl;
	//this->RunString("PRINT;-selected_output false");
	this->SetSelectedOutputStringOn(false);
	
	// Add selected output for mass of water
	VAR pvar1, pvar2;
	VarInit(&pvar1);
	VarInit(&pvar2);	
	selout.append("SELECTED_OUTPUT; -reset false\n");
	selout.append("USER_PUNCH; -heading mass_water\n");
	selout.append("10 PUNCH TOT(\"water\")\n");
	selout.append("END\n");
	
	// Make melt
	this->RunString(selout.c_str());
	this->GetSelectedOutputValue(1,0, &pvar1);	
	//std::cerr << this->GetOutputString() << std::endl;

	// Make snowpack
	std::ostringstream strm0;
	strm0 << "MIX 1 Make snowpack\n";
	strm0 << ipack << " 1\n";
	strm0 << imelt << " -1\n";
	strm0 << "SAVE solution " << ipack  << "\n";
	this->RunString(strm0.str().c_str());
	this->GetSelectedOutputValue(1,0, &pvar2);
	//std::cerr << this->GetOutputString() << std::endl;

	std::ostringstream strm2;
	strm2 << "SOLUTION_MIX " << imelt << "\n";
	strm2 << imelt << " " << 1./pvar1.dVal << "\n";
	strm2 << "END\n";
	strm2 << "SOLUTION_MIX " << ipack << "\n";
	strm2 << ipack << " " << 1./pvar2.dVal << "\n";
	strm2 << "END" << "\n";
	
	//std::cerr << strm2.str().c_str();
	this->RunString(strm2.str().c_str());
	//std::cerr << this->GetOutputString() << std::endl;
	return 1;
}