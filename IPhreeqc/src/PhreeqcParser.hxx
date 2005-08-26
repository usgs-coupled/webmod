#if !defined(__PHREEQC_PARSER_HXX_INC)
#define __PHREEQC_PARSER_HXX_INC


#include <string>
//#include <iosfwd>


struct PhreeqcStop{};


class CPhreeqcParser
{
public:
	CPhreeqcParser(std::istream& input);
	virtual ~CPhreeqcParser(void);

	// const char* GetErrorMsg(void);

	// overrides
	int get_logical_line(int *l);
	int get_line(void);
	// int error_msg (const char *err_str, const int stop);


protected:
	std::istream& m_input_stream;
	//std::ostream& m_output_stream;
	//std::ostream& m_error_stream;

	//std::string   m_errors;
};


#endif // __PHREEQC_PARSER_HXX_INC
