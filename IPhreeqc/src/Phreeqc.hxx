// Phreeqc.h: interface for the CPhreeqc class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PHREEQC_H__9BBA7291_1046_44CB_AEE7_F4FF09F15FF0__INCLUDED_)
#define AFX_PHREEQC_H__9BBA7291_1046_44CB_AEE7_F4FF09F15FF0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdarg.h>
#include <stdio.h>
#include <string>

#include "../include/IPhreeqc.h"

class CPhreeqcParser;
class IErrorReporter;
class COutput;

class CPhreeqc  
{
protected:
	CPhreeqc(void);
public:
	static CPhreeqc singleton;
	virtual ~CPhreeqc(void);

	int CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie, IErrorReporter* pErrorReporter = 0);
	int LoadDatabase(std::istream& is, IErrorReporter* pErrorReporter = 0);
	int Run(std::istream& is, bool output_on, bool error_on, bool log_on, bool selected_output_on, IErrorReporter* pErrorReporter = 0);
	int RunWithCallback(PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie, int output_on, int error_on, int log_on, int selected_output_on, IErrorReporter* pErrorReporter = 0);


	IErrorReporter* GetErrorReporter(void);

	static int MessageCallback(const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args);


protected:
	CPhreeqcParser* m_pPhreeqcParser;
	IErrorReporter* m_pErrorReporter;
	COutput*        m_pOutput;
	bool            m_bDatabaseLoaded;
	std::string     m_sSelectedOutputFileName;
	bool            m_bSelectedOutputOn;

public:
	int GetLine(void);
	int ErrorMsg(const char* err_str, const int stop);
	int OpenSelectedOutput(const char* filename);
};

#endif // !defined(AFX_PHREEQC_H__9BBA7291_1046_44CB_AEE7_F4FF09F15FF0__INCLUDED_)
