#if !defined(__OUTPUT_HXX_INC)
#define __OUTPUT_HXX_INC


#include <string>
class IErrorReporter;


class COutput
{
public:
	COutput(bool output_on = false, bool error_on = false, bool log_on = false);
	virtual ~COutput(void);
public:
	int Open(IErrorReporter* pErrorReporter = 0);	
	int Close(IErrorReporter* pErrorReporter = 0);
	const char* GetOutputFileName(void)const;
	const char* GetErrorFileName(void)const;
	const char* GetLogFileName(void)const;
// COMMENT: {11/17/2003 4:56:52 PM}	const char* GetSelOutFileName(void)const;

// COMMENT: {11/17/2003 4:57:36 PM}	int OpenSelectedOutput(const char* filename, IErrorReporter* pErrorReporter = 0);
	
protected:
	bool m_bOutputOn;
	bool m_bErrorOn;
	bool m_bLogOn;
// COMMENT: {11/17/2003 4:57:04 PM}	bool m_bSelOutOn;
	std::string m_sOutputFileName;
	std::string m_sErrorFileName;
	std::string m_sLogFileName;
// COMMENT: {11/17/2003 4:57:09 PM}	std::string m_sSelOutName;
};


#endif // __OUTPUT_HXX_INC
