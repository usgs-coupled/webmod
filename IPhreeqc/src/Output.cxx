#include <stdio.h> // fopen fclose

#include "Output.hxx"
#include "Debug.h" // ASSERT

#include "phreeqcns.hxx"
#include "ErrorReporter.hxx"

#if defined(_WIN32)
const char NULL_DEVICE[] = "NUL";
#else
const char NULL_DEVICE[] = "/dev/null";
#endif

const char OUTPUT_FILENAME[] = "phreeqc.out";
const char ERROR_FILENAME[]  = "phreeqc.err";
const char LOG_FILENAME[]    = "phreeqc.log";


COutput::COutput(bool bOutputOn, bool bErrorOn, bool bLogOn)
: m_bOutputOn(bOutputOn)
, m_bErrorOn(bErrorOn)
, m_bLogOn(bLogOn)
, m_sOutputFileName(OUTPUT_FILENAME)
, m_sErrorFileName(ERROR_FILENAME)
, m_sLogFileName(LOG_FILENAME)
{
}

COutput::~COutput(void)
{
	this->Close();
}

const char* COutput::GetOutputFileName(void)const
{
	return (this->m_bOutputOn) ? this->m_sOutputFileName.c_str() : NULL_DEVICE;
}

const char* COutput::GetErrorFileName(void)const
{
	return (this->m_bOutputOn) ? this->m_sErrorFileName.c_str() : NULL_DEVICE;
}

const char* COutput::GetLogFileName(void)const
{
	return (this->m_bOutputOn) ? this->m_sLogFileName.c_str() : NULL_DEVICE;
}

int COutput::Open(IErrorReporter* pErrorReporter)
{
	const char* name;
	const char err[] = "Unable to open ";
	int error_count = 0;

	ASSERT(phreeqc::output     == NULL);
	ASSERT(phreeqc::error_file == NULL);
	ASSERT(phreeqc::log_file   == NULL);
	ASSERT(phreeqc::punch_file == NULL);

	name = this->GetOutputFileName();
	phreeqc::output = ::fopen(name, "w");
	if (phreeqc::output == NULL)
	{
		if (pErrorReporter) {
			std::string str_err(err);
			str_err += name;
			str_err += "\n";
			pErrorReporter->AddError(str_err.c_str());
		}
		++error_count;
	}

	name = this->GetErrorFileName();
	phreeqc::error_file = ::fopen(name, "w");
	if (phreeqc::error_file == NULL)
	{
		if (pErrorReporter) {
			std::string str_err(err);
			str_err += name;
			str_err += "\n";
			pErrorReporter->AddError(str_err.c_str());
		}
		++error_count;
	}

	name = this->GetLogFileName();
	phreeqc::log_file = ::fopen(name, "w");
	if (phreeqc::log_file == NULL)
	{
		if (pErrorReporter) {
			std::string str_err(err);
			str_err += name;
			str_err += "\n";
			pErrorReporter->AddError(str_err.c_str());
		}
		++error_count;
	}

	return error_count;
}

int COutput::Close(IErrorReporter* pErrorReporter)
{
	const char err[] = "Unable to close ";
	int error_count = 0;

	ASSERT(phreeqc::output != NULL);
	ASSERT(phreeqc::error_file != NULL);
	ASSERT(phreeqc::log_file != NULL);

	if (phreeqc::output) {
		if (::fclose(phreeqc::output)) {
			if (pErrorReporter) {
				std::string str_err(err);
				str_err += this->GetOutputFileName();
				str_err += ".\n";
				pErrorReporter->AddError(str_err.c_str());
			}
			++error_count;
		}
		phreeqc::output = NULL;
	}

	if (phreeqc::error_file) {
		if (::fclose(phreeqc::error_file)) {
			if (pErrorReporter) {
				std::string str_err(err);
				str_err += this->GetErrorFileName();
				str_err += ".\n";
				pErrorReporter->AddError(str_err.c_str());
			}
			++error_count;
		}
		phreeqc::error_file = NULL;
	}

	if (phreeqc::log_file) {
		if (::fclose(phreeqc::log_file)) {
			if (pErrorReporter) {
				std::string str_err(err);
				str_err += this->GetLogFileName();
				str_err += ".\n";
				pErrorReporter->AddError(str_err.c_str());
			}
			++error_count;
		}
		phreeqc::log_file = NULL;
	}

	if (phreeqc::punch_file) {
		if (::fclose(phreeqc::punch_file)) {
			if (pErrorReporter) {
// COMMENT: {11/17/2003 4:58:42 PM}				std::string str_err(err);
// COMMENT: {11/17/2003 4:58:42 PM}				str_err += this->GetSelOutFileName();
// COMMENT: {11/17/2003 4:58:42 PM}				str_err += ".\n";
// COMMENT: {11/17/2003 4:58:42 PM}				pErrorReporter->AddError(str_err.c_str());
				pErrorReporter->AddError("Unable to close selected output file.\n");
			}
			++error_count;
		}
		phreeqc::punch_file = NULL;
	}

	return error_count;
}

// COMMENT: {11/17/2003 4:56:03 PM}int COutput::OpenSelectedOutput(const char* filename, IErrorReporter* pErrorReporter)
// COMMENT: {11/17/2003 4:56:03 PM}{
// COMMENT: {11/17/2003 4:56:03 PM}	const char* name;
// COMMENT: {11/17/2003 4:56:03 PM}	const char err[] = "Unable to open ";
// COMMENT: {11/17/2003 4:56:03 PM}	int error_count = 0;
// COMMENT: {11/17/2003 4:56:03 PM}
// COMMENT: {11/17/2003 4:56:03 PM}	assert(phreeqc::punch_file != NULL);
// COMMENT: {11/17/2003 4:56:03 PM}	if (phreeqc::punch_file) {
// COMMENT: {11/17/2003 4:56:03 PM}		if (::fclose(phreeqc::punch_file)) {
// COMMENT: {11/17/2003 4:56:03 PM}			if (pErrorReporter) {
// COMMENT: {11/17/2003 4:56:03 PM}				std::string str_err(err);
// COMMENT: {11/17/2003 4:56:03 PM}				str_err += this->GetSelOutFileName();
// COMMENT: {11/17/2003 4:56:03 PM}				str_err += ".\n";
// COMMENT: {11/17/2003 4:56:03 PM}				pErrorReporter->AddError(str_err.c_str());
// COMMENT: {11/17/2003 4:56:03 PM}			}
// COMMENT: {11/17/2003 4:56:03 PM}			++error_count;
// COMMENT: {11/17/2003 4:56:03 PM}		}
// COMMENT: {11/17/2003 4:56:03 PM}		phreeqc::punch_file = NULL;
// COMMENT: {11/17/2003 4:56:03 PM}	}
// COMMENT: {11/17/2003 4:56:03 PM}
// COMMENT: {11/17/2003 4:56:03 PM}	this->m_sSelOutName = filename;
// COMMENT: {11/17/2003 4:56:03 PM}
// COMMENT: {11/17/2003 4:56:03 PM}	name = this->GetSelOutFileName();
// COMMENT: {11/17/2003 4:56:03 PM}	phreeqc::punch_file = ::fopen(name, "w");
// COMMENT: {11/17/2003 4:56:03 PM}	if (phreeqc::punch_file == NULL)
// COMMENT: {11/17/2003 4:56:03 PM}	{
// COMMENT: {11/17/2003 4:56:03 PM}		if (pErrorReporter) {
// COMMENT: {11/17/2003 4:56:03 PM}			std::string str_err(err);
// COMMENT: {11/17/2003 4:56:03 PM}			str_err += name;
// COMMENT: {11/17/2003 4:56:03 PM}			str_err += ".\n";
// COMMENT: {11/17/2003 4:56:03 PM}			pErrorReporter->AddError(str_err.c_str());
// COMMENT: {11/17/2003 4:56:03 PM}		}
// COMMENT: {11/17/2003 4:56:03 PM}		++error_count;
// COMMENT: {11/17/2003 4:56:03 PM}	}
// COMMENT: {11/17/2003 4:56:03 PM}	return error_count;
// COMMENT: {11/17/2003 4:56:03 PM}}
