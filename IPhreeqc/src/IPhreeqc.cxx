#if defined(_DEBUG)
#pragma warning(disable : 4786)   // disable truncation warning
#endif
#include <iostream>  // std::cout
#include <fstream>   // std::ifstream
#include <memory>    // std::auto_ptr
#include <sstream>   // std::istringstream
#include <cstdarg>

#if defined(_WIN32)
#include <windows.h> // ::OutputDebugString
#endif

#include "phreeqcns.hxx"
#include "Phreeqc.hxx"
#include "ErrorReporter.hxx"
#include "SelectedOutput.hxx"
#include "../include/Var.h"
#include "../include/IPhreeqc.h"


const std::string& GetAccumulatedLines(void);
void ClearAccumulatedLines(void);


static std::string s_string_input;
static std::ostringstream s_oss;
static CErrorReporter<std::ostringstream> s_errorReporter;


int
LoadDatabase(const char* filename)
{
	try 
	{
		s_errorReporter.Clear();

		std::ifstream ifs;
		ifs.open(filename);

		if (!ifs.is_open()) {
			std::string str("LoadDatabase: Unable to open:");
			str += filename;
			str += ".\n";
			s_errorReporter.AddError(str.c_str());
			return 1;
		}

		CSelectedOutput::singleton.Clear();
		int errors = CPhreeqc::singleton.LoadDatabase(ifs, &s_errorReporter);
		return errors;
	}
	catch(...)
	{
		s_errorReporter.AddError("LoadDatabase: An unhandled exception occured.\n");
	}
	return 1;
}

void 
OutputLastError(void)
{
	std::cout << s_errorReporter.GetOS()->str().c_str() << std::endl;
}

#if defined(_WIN32)
void
DebugOutputLastError(void)
{
	std::istringstream iss(s_errorReporter.GetOS()->str());
	std::string line;
	while (std::getline(iss, line)) {
		::OutputDebugString(line.c_str());
		::OutputDebugString("\n");
	}
}
#endif


VRESULT
AccumulateLine(const char *line)
{
	try
	{
		s_errorReporter.Clear();
		s_string_input.append(line);
		s_string_input.append("\n");
		return VR_OK;
	}
	catch (...)
	{
		s_errorReporter.AddError("AccumulateLine: An unhandled exception occured.\n");
	}
	return VR_OUTOFMEMORY;
}

int
Run(int output_on, int error_on, int log_on, int selected_output_on)
{
	try 
	{
		s_errorReporter.Clear();
		CSelectedOutput::singleton.Clear();

		std::istringstream iss(GetAccumulatedLines());

		int errors = CPhreeqc::singleton.Run(
			iss, 
			(output_on != 0),
			(error_on != 0),
			(log_on != 0),
			(selected_output_on != 0),
			&s_errorReporter
			);

		if (errors > 0) {
			s_errorReporter.AddError("Start Input:\n");
			s_errorReporter.AddError(s_string_input.c_str());
			s_errorReporter.AddError("End Input:\n");
		}
		ClearAccumulatedLines();

		return errors;
	}
	catch(...)
	{
		s_errorReporter.AddError("Run: An unhandled exception occured.\n");
	}
	return 1;
}

int
RunFile(const char* filename, int output_on, int error_on, int log_on, int selected_output_on)
{
	try 
	{
		s_errorReporter.Clear();

		std::ifstream ifs;
		ifs.open(filename);
		if (!ifs.is_open()) {
			std::string str("RunFile: Unable to open:");
			str += filename;
			str += ".\n";
			s_errorReporter.AddError(str.c_str());
			return 1;
		}


		CSelectedOutput::singleton.Clear();

		int errors = CPhreeqc::singleton.Run(
			ifs, 
			(output_on != 0),
			(error_on != 0),
			(log_on != 0),
			(selected_output_on != 0),
			&s_errorReporter
			);

		return errors;
	}
	catch(...)
	{
		s_errorReporter.AddError("RunFile: An unhandled exception occured.\n");
	}
	return 1;
}


int
GetSelectedOutputRowCount(void)
{
	return CSelectedOutput::singleton.GetRowCount();
}

int
GetSelectedOutputColumnCount(void)
{
	return CSelectedOutput::singleton.GetColCount();
}

VRESULT
GetSelectedOutputValue(int row, int col, VAR* pVAR)
{
	s_errorReporter.Clear();
	if (!pVAR) {
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDARG pVar is NULL.\n");
		return VR_INVALIDARG;
	}

	VRESULT v = CSelectedOutput::singleton.Get(row, col, pVAR);
	switch (v) {
	case VR_OK:
		break;
	case VR_OUTOFMEMORY:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_OUTOFMEMORY Out of memory.\n");
		break;
	case VR_BADVARTYPE:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_BADVARTYPE pVar must be initialized(VarInit) and/or cleared(VarClear).\n");
		break;
	case VR_INVALIDARG:
		// not possible
		break;
	case VR_INVALIDROW:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDROW Row index out of range.\n");
		break;
	case VR_INVALIDCOL:
		s_errorReporter.AddError("GetSelectedOutputValue: VR_INVALIDCOL Column index out of range.\n");
		break;
	}
	return v;
}

size_t
AddError(const char* error_msg)
{
	return s_errorReporter.AddError(error_msg);
}

void
OutputLines(void)
{
	std::cout << s_string_input.c_str() << std::endl;
}

#if defined(WIN32)
void
DebugOutputLines(void)
{
	std::istringstream iss(s_string_input);
	std::string line;
	while (std::getline(iss, line)) {
		::OutputDebugString(line.c_str());
		::OutputDebugString("\n");
	}
}
#endif

const std::string&
GetAccumulatedLines(void)
{
	return s_string_input;
}

void
ClearAccumulatedLines(void)
{
	s_string_input.erase();
}

int
RunWithCallback(PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie, int output_on, int error_on, int log_on, int selected_output_on)
{
	try 
	{
		s_errorReporter.Clear();
		CSelectedOutput::singleton.Clear();

		std::istringstream iss(GetAccumulatedLines());

		int errors = CPhreeqc::singleton.RunWithCallback(
			pfn_pre,
			pfn_post,
			cookie,
			(output_on != 0),
			(error_on != 0),
			(log_on != 0),
			(selected_output_on != 0),
			&s_errorReporter
			);

		if (errors > 0) {
			s_errorReporter.AddError("Start Input:\n");
			s_errorReporter.AddError(s_string_input.c_str());
			s_errorReporter.AddError("End Input:\n");
		}
		ClearAccumulatedLines();

		return errors;
	}
	catch(...)
	{
		s_errorReporter.AddError("Run: An unhandled exception occured.\n");
	}
	return 1;
}

int
CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie)
{
	s_errorReporter.Clear();
	return CPhreeqc::singleton.CatchErrors(pfn, cookie, &s_errorReporter);
}
