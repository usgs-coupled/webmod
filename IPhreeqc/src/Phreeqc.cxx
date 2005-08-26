// Phreeqc.cpp: implementation of the CPhreeqc class.
//
//////////////////////////////////////////////////////////////////////
#if defined(WIN32)
#include <windows.h>
const char NULL_DEVICE[] = "NUL";
#else
const char NULL_DEVICE[] = "/dev/null";
#endif

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

#include <stdio.h>

#include "Debug.h" // ASSERT

#include "Phreeqc.hxx"
#include "phreeqcns.hxx"
#include "PhreeqcParser.hxx"
#include "Output.hxx"
#include "ErrorReporter.hxx"

extern "C" int n_user_punch_index;
extern const std::string& GetAccumulatedLines(void);

CPhreeqc CPhreeqc::singleton;


CPhreeqc::CPhreeqc(void)
	: m_pPhreeqcParser(0)
	, m_pErrorReporter(0)
	, m_pOutput(0)
	, m_bDatabaseLoaded(false)
{
	ASSERT(phreeqc::input       == NULL);
	ASSERT(phreeqc::output      == NULL);
	ASSERT(phreeqc::log_file    == NULL);
	ASSERT(phreeqc::punch_file  == NULL);
	ASSERT(phreeqc::dump_file   == NULL);
	ASSERT(phreeqc::error_file  == NULL);
	ASSERT(phreeqc::input_error == 0);

	phreeqc::add_message_callback(phreeqc::default_handler, NULL);
	phreeqc::add_message_callback(CPhreeqc::MessageCallback, this);

	phreeqc::phast = FALSE;
}

CPhreeqc::~CPhreeqc(void)
{
	if (this->m_pPhreeqcParser) delete this->m_pPhreeqcParser;
	phreeqc::clean_up();
}

IErrorReporter* CPhreeqc::GetErrorReporter(void)
{
	return this->m_pErrorReporter;
}

int CPhreeqc::LoadDatabase(std::istream& is, IErrorReporter* pErrorReporter)
{
	try
	{
		ASSERT(this->m_pErrorReporter == 0);
		this->m_pErrorReporter = pErrorReporter;

		// phreeqc
		//
		phreeqc::clean_up();
		phreeqc::add_message_callback(phreeqc::default_handler, NULL);
		phreeqc::add_message_callback(CPhreeqc::MessageCallback, this);

		CPhreeqcParser parser(is);
		this->m_pPhreeqcParser = &parser;

		COutput output(false, false, false);
		if (output.Open(pErrorReporter)) {
			throw std::string("Cannot open output files.");
		}

		// database needs to be non-null
		ASSERT(phreeqc::database_file == NULL);
		phreeqc::database_file = fopen(NULL_DEVICE, "r");
		ASSERT(phreeqc::database_file != NULL);		

		phreeqc::input_error = 0;

		// phreeqc
		//
		phreeqc::state = INITIALIZE;
		phreeqc::initialize();
		phreeqc::phast = FALSE;

		phreeqc::simulation = 0;
		phreeqc::dup_print("Reading data base.", TRUE);

		phreeqc::read_input();
		phreeqc::tidy_model();
	}
	catch (std::string err)
	{
		++phreeqc::input_error;
	}
	catch (...)
	{
		std::string err_msg("LoadDatabase: Unhandled exception.\n");
		if (pErrorReporter) {
			pErrorReporter->AddError(err_msg.c_str());
		}
		else {
			std::cerr << err_msg.c_str() << std::endl;
		}

		++phreeqc::input_error;
	}

	if (phreeqc::database_file != NULL) ::fclose(phreeqc::database_file);
	phreeqc::database_file = NULL;

	this->m_pPhreeqcParser = 0;
	this->m_pErrorReporter = 0;

	this->m_bDatabaseLoaded = (phreeqc::input_error == 0);

	return phreeqc::input_error;
}

int CPhreeqc::Run(std::istream& is, bool output_on, bool error_on, bool log_on, bool selected_output_on, IErrorReporter* pErrorReporter)
{
	try
	{
		ASSERT(this->m_pErrorReporter == 0);
		this->m_pErrorReporter = pErrorReporter;

		if (!this->m_bDatabaseLoaded) {
			++phreeqc::input_error;
			throw std::string("phreeqc: No Database is loaded.\n");
		}
		phreeqc::input_error = 0;

		this->m_bSelectedOutputOn = selected_output_on;

		COutput output_files(output_on, error_on, log_on);
		if (output_files.Open(pErrorReporter)) {
			++phreeqc::input_error;
            throw std::string("phreeqc: Cannot open output files.\n");
		}
		ASSERT(this->m_pOutput == 0);
		this->m_pOutput = &output_files;

		ASSERT(this->m_pPhreeqcParser == 0);
		CPhreeqcParser parser(is);
		this->m_pPhreeqcParser = &parser;

		// run phreeqc
		//
		ASSERT(phreeqc::punch_file == NULL);

		char token[80];
		for (phreeqc::simulation = 1; ; ++phreeqc::simulation) {

			::sprintf(token, "Reading input data for simulation %d.", phreeqc::simulation);
			phreeqc::dup_print(token, TRUE);

			int save_punch_in = phreeqc::punch.in;

			if (phreeqc::read_input() == EOF) break;

			if (phreeqc::simulation > 1 && save_punch_in == TRUE && phreeqc::punch.new_def == TRUE)
			{
				this->ErrorMsg("phreeqc: Warning SELECTED_OUTPUT has been redefined.\n", CONTINUE);
			}
			if (phreeqc::simulation > 1 && phreeqc::keyword[39].keycount > 0)
			{
				this->ErrorMsg("phreeqc: Warning USER_PUNCH has been redefined.\n", CONTINUE);
			}

			// this->HandleSelectedOutput(selected_output_on);
			// selected_output truth table
			// if (punch.in == TRUE)
			//|
			//|pr    |punch   |sel_out_on|punch_file|tidy_punch|
			//|.punch|.new_def|          |          |          |
			//--------------------------------------------------
			//|   T  |   T    |    T     | file_name|.. NO     | 
			//|   T  |   T    |    F     | null_dev |.. NO     |
			//|   T  |   F    |    T     | file_name|   YES    |
			//|   T  |   F    |    F     | null_dev |  maybe   |
			//|   F  |   T    |    T     | null_dev |.. NO     |
			//|   F  |   T    |    F     | null_dev |.. NO     |
			//|   F  |   F    |    T     | null_dev |  maybe   |
			//|   F  |   F    |    F     | null_dev |  maybe   |
			if (phreeqc::punch.in == TRUE)
			{
				if (phreeqc::pr.punch == FALSE || !selected_output_on)
				{
					ASSERT(phreeqc::punch_file == NULL);
					phreeqc::punch_file = ::fopen(NULL_DEVICE, "w");
					if (phreeqc::punch_file == NULL)
					{
						++phreeqc::input_error;
						this->ErrorMsg("phreeqc: Unable to open null file.\n", CONTINUE);
					}
				}
				else 
				{
					if (phreeqc::punch.new_def == FALSE)
					{
						if (phreeqc::punch_file == NULL)
						{
							phreeqc::punch_file = ::fopen(this->m_sSelectedOutputFileName.c_str(), "w");
							if (phreeqc::punch_file == NULL)
							{
								std::string str("phreeqc: Unable to open ");
								str += this->m_sSelectedOutputFileName.c_str();
								str += ".\n";
								++phreeqc::input_error;
								this->ErrorMsg(str.c_str(), CONTINUE);
							}
							else
							{
								// output selected_output headings
								phreeqc::punch.new_def = TRUE;
								phreeqc::tidy_punch();
							}
						}
					}
					else 
					{
						// should have been opened in read_selected_output/OpenSelectedOutput
						ASSERT(phreeqc::punch_file != NULL);
					}
				}
				phreeqc::pr.punch = TRUE; // always TRUE when (phreeqc::punch.in == TRUE)
			}

			n_user_punch_index = -1;

			phreeqc::tidy_model();
			if (phreeqc::input_error > 0) {
				throw std::string("phreeqc: Input contains errors.\n");
			}

/*
 *   Calculate distribution of species for initial solutions
 */
			if (phreeqc::new_solution) phreeqc::initial_solutions(TRUE);

/*
 *   Calculate distribution for exchangers
 */
			if (phreeqc::new_exchange) phreeqc::initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
			if (phreeqc::new_surface) phreeqc::initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
			if (phreeqc::new_gas_phase) phreeqc::initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
			phreeqc::reactions();
/*
 *   Calculate inverse models
 */
			phreeqc::inverse_models();

/*
 *   Copy
 */
			if (phreeqc::new_copy) phreeqc::copy_entities();
/*
 *   End of simulation
 */
			phreeqc::dup_print("End of simulation.", TRUE);
		}

		phreeqc::dup_print ("End of run.", TRUE);

	}
	catch (std::string str)
	{
		if (pErrorReporter) {
			pErrorReporter->AddError(str.c_str());
		}
		else {
			std::cerr << str.c_str() << std::endl;
		}
	}
	catch (...)
	{
		const char err_msg[] = "phreeqc: Unhandled exception occured.\n";
		if (pErrorReporter) {
			pErrorReporter->AddError(err_msg);
		}
		else {
			std::cerr << err_msg << std::endl;
		}
	}

	this->m_pErrorReporter = 0;
	this->m_pPhreeqcParser = 0;
	this->m_pOutput        = 0;

	return phreeqc::input_error;
}

int CPhreeqc::GetLine(void)
{
	if (this->m_pPhreeqcParser) {
		return this->m_pPhreeqcParser->get_line();
	}
	return EOF;
}

int CPhreeqc::ErrorMsg(const char* err_str, const int stop)
{
	std::string str("ERROR: ");
	str += err_str;
	str += "\n";
	if (stop == STOP) {
		str += "Stopping.\n";
	}
	if (IErrorReporter* pErrorReporter = this->GetErrorReporter()) {
		pErrorReporter->AddError(str.c_str());
	}
	else {
		std::cerr << str;
		std::cerr << std::endl;
	}
	if (stop == STOP) {
		throw PhreeqcStop();
	}
	return OK;
}

int CPhreeqc::OpenSelectedOutput(const char* filename)
{
	this->m_sSelectedOutputFileName = filename;
	if (phreeqc::punch_file)
	{
		int close = ::fclose(phreeqc::punch_file);
		ASSERT(close == 0);
		phreeqc::punch_file = NULL;
	}

	if (this->m_bSelectedOutputOn)
	{
		ASSERT(phreeqc::punch_file == NULL);
		phreeqc::punch_file = ::fopen(filename, "w");
		if (phreeqc::punch_file == NULL)
		{
			std::string str("phreeqc: Unable to open ");
			str += this->m_sSelectedOutputFileName.c_str();
			str += ".\n";
			++phreeqc::input_error;
			this->ErrorMsg(str.c_str(), CONTINUE);
			return 1;
		}
	}
	else 
	{
		ASSERT(phreeqc::punch_file == NULL);
		phreeqc::punch_file = ::fopen(NULL_DEVICE, "w");
		if (phreeqc::punch_file == NULL)
		{
			++phreeqc::input_error;
			this->ErrorMsg("phreeqc: Unable to open null file.\n", CONTINUE);
			return 1;
		}
	}
	return 0;
}

int CPhreeqc::RunWithCallback(PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie, int output_on, int error_on, int log_on, int selected_output_on, IErrorReporter* pErrorReporter)
{
	try
	{
		ASSERT(this->m_pErrorReporter == 0);
		this->m_pErrorReporter = pErrorReporter;

		if (!this->m_bDatabaseLoaded) {
			++phreeqc::input_error;
			throw std::string("phreeqc: No Database is loaded.\n");
		}
		phreeqc::input_error = 0;

		this->m_bSelectedOutputOn = (selected_output_on != 0);

		COutput output_files((output_on != 0), (error_on != 0), (log_on != 0));
		if (output_files.Open(pErrorReporter)) {
			++phreeqc::input_error;
            throw std::string("phreeqc: Cannot open output files.\n");
		}
		ASSERT(this->m_pOutput == 0);
		this->m_pOutput = &output_files;

		if (pfn_pre) {
			pfn_pre(cookie);
		}

		std::istringstream iss(GetAccumulatedLines());

		ASSERT(this->m_pPhreeqcParser == 0);
		CPhreeqcParser parser(iss);
		this->m_pPhreeqcParser = &parser;

		// run phreeqc
		//
		ASSERT(phreeqc::punch_file == NULL);

		char token[80];
		for (phreeqc::simulation = 1; ; ++phreeqc::simulation) {

			::sprintf(token, "Reading input data for simulation %d.", phreeqc::simulation);
			phreeqc::dup_print(token, TRUE);

			int save_punch_in = phreeqc::punch.in;

			if (phreeqc::read_input() == EOF) break;

			if (phreeqc::simulation > 1 && save_punch_in == TRUE && phreeqc::punch.new_def == TRUE)
			{
				this->ErrorMsg("phreeqc: Warning SELECTED_OUTPUT has been redefined.\n", CONTINUE);
			}
			if (phreeqc::simulation > 1 && phreeqc::keyword[39].keycount > 0)
			{
				this->ErrorMsg("phreeqc: Warning USER_PUNCH has been redefined.\n", CONTINUE);
			}

			// this->HandleSelectedOutput(selected_output_on);
			// selected_output truth table
			// if (punch.in == TRUE)
			//|
			//|pr    |punch   |sel_out_on|punch_file|tidy_punch|
			//|.punch|.new_def|          |          |          |
			//--------------------------------------------------
			//|   T  |   T    |    T     | file_name|.. NO     | 
			//|   T  |   T    |    F     | null_dev |.. NO     |
			//|   T  |   F    |    T     | file_name|   YES    |
			//|   T  |   F    |    F     | null_dev |  maybe   |
			//|   F  |   T    |    T     | null_dev |.. NO     |
			//|   F  |   T    |    F     | null_dev |.. NO     |
			//|   F  |   F    |    T     | null_dev |  maybe   |
			//|   F  |   F    |    F     | null_dev |  maybe   |
			if (phreeqc::punch.in == TRUE)
			{
				if (phreeqc::pr.punch == FALSE || !selected_output_on)
				{
					ASSERT(phreeqc::punch_file == NULL);
					phreeqc::punch_file = ::fopen(NULL_DEVICE, "w");
					if (phreeqc::punch_file == NULL)
					{
						++phreeqc::input_error;
						this->ErrorMsg("phreeqc: Unable to open null file.\n", CONTINUE);
					}
				}
				else 
				{
					if (phreeqc::punch.new_def == FALSE)
					{
						if (!phreeqc::punch_file)
						{
							phreeqc::punch_file = ::fopen(this->m_sSelectedOutputFileName.c_str(), "w");
							if (phreeqc::punch_file == NULL)
							{
								std::string str("phreeqc: Unable to open ");
								str += this->m_sSelectedOutputFileName.c_str();
								str += ".\n";
								++phreeqc::input_error;
								this->ErrorMsg(str.c_str(), CONTINUE);
							}
							else
							{
								// output selected_output headings
								phreeqc::punch.new_def = TRUE;
								phreeqc::tidy_punch();
							}
						}
					}
					else 
					{
						// should have been opened in read_selected_output/OpenSelectedOutput
						ASSERT(phreeqc::punch_file != NULL);
					}
				}
				phreeqc::pr.punch = TRUE; // always TRUE when (phreeqc::punch.in == TRUE)
			}

			n_user_punch_index = -1;

			phreeqc::tidy_model();
			if (phreeqc::input_error > 0) {
				throw std::string("phreeqc: Input contains errors.\n");
			}

/*
 *   Calculate distribution of species for initial solutions
 */
			if (phreeqc::new_solution) phreeqc::initial_solutions(TRUE);

/*
 *   Calculate distribution for exchangers
 */
			if (phreeqc::new_exchange) phreeqc::initial_exchangers(TRUE);
/*
 *   Calculate distribution for surfaces
 */
			if (phreeqc::new_surface) phreeqc::initial_surfaces(TRUE);
/*
 *   Calculate initial gas composition
 */
			if (phreeqc::new_gas_phase) phreeqc::initial_gas_phases(TRUE);
/*
 *   Calculate reactions
 */
			phreeqc::reactions();
/*
 *   Calculate inverse models
 */
			phreeqc::inverse_models();

/*
 *   Copy
 */
			if (phreeqc::new_copy) phreeqc::copy_entities();
/*
 *   End of simulation
 */
			phreeqc::dup_print("End of simulation.", TRUE);
		}

		phreeqc::dup_print ("End of run.", TRUE);

		if (pfn_post) {
			pfn_post(cookie);
		}

	}
	catch (std::string str)
	{
		if (pErrorReporter) {
			pErrorReporter->AddError(str.c_str());
		}
		else {
			std::cerr << str.c_str() << std::endl;
		}
	}
	catch (...)
	{
		const char err_msg[] = "phreeqc: Unhandled exception occured.\n";
		if (pErrorReporter) {
			pErrorReporter->AddError(err_msg);
		}
		else {
			std::cerr << err_msg << std::endl;
		}
	}

	this->m_pErrorReporter = 0;
	this->m_pPhreeqcParser = 0;
	this->m_pOutput        = 0;

	return phreeqc::input_error;
}

int CPhreeqc::MessageCallback(const int type, const char *err_str, const int stop, void *cookie, const char *, va_list args)
{
	if (type == phreeqc::OUTPUT_ERROR) {
		CPhreeqc *self = (CPhreeqc*)cookie;
		self->ErrorMsg(err_str, stop);
	}
	return OK;
}

int CPhreeqc::CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie, IErrorReporter* pErrorReporter)
{
	int rvalue = OK;
	try
	{
		ASSERT(this->m_pErrorReporter == 0);
		this->m_pErrorReporter = pErrorReporter;

		phreeqc::input_error = 0;

		if (pfn) {
			rvalue = pfn(cookie);
		}

	}
	catch (...)
	{
		const char err_msg[] = "phreeqc: Unhandled exception occured.\n";
		if (pErrorReporter) {
			pErrorReporter->AddError(err_msg);
		}
		else {
			std::cerr << err_msg << std::endl;
		}
	}

	this->m_pErrorReporter = 0;

	return rvalue;
}
