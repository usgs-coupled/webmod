#include <cctype>
#include <iostream>  // std::cerr std::cin std::cout

#include "PhreeqcParser.hxx"
#include "phreeqcns.hxx"
#include "Debug.h"   // ASSERT

CPhreeqcParser::CPhreeqcParser(std::istream& input)
: m_input_stream(input)
{
}

CPhreeqcParser::~CPhreeqcParser(void)
{
}

int CPhreeqcParser::get_logical_line(int *l)
{
/*
 *   Reads file fp until end of line, ";", or eof
 *   stores characters in line_save
 *   reallocs line_save and line if more space is needed
 *
 *   returns:
 *           EOF on empty line on end of file or
 *           OK otherwise
 *           *l returns length of line
 */
	int i, j;
	char c;
	i = 0;
	for (;;) {
		j = this->m_input_stream.get();
		if (j == EOF) break;
		c = (char) j;
		if (c == '\\') {
			j = this->m_input_stream.get();
			if (j == EOF) break;
			j = this->m_input_stream.get();
			if (j == EOF) break;
			c = (char) j;
		}
		if (c == ';' || c == '\n') break;
		if ( i + 20 >= phreeqc::max_line) {
			phreeqc::max_line *= 2;
			phreeqc::line_save = (char *) phreeqc::PHRQ_realloc (phreeqc::line_save, (size_t) phreeqc::max_line * sizeof(char));
			if (phreeqc::line_save == NULL) phreeqc::malloc_error();
			phreeqc::line = (char *) phreeqc::PHRQ_realloc (phreeqc::line, (size_t) phreeqc::max_line * sizeof(char));
			if (phreeqc::line == NULL) phreeqc::malloc_error();
		}
		phreeqc::line_save[i++] = c;
	}
	if (j == EOF && i == 0) {
		*l = 0;
		phreeqc::line_save[i] = '\0';
		return(EOF);
	} 
	phreeqc::line_save[i] = '\0';
	*l = i;
	return(OK);
}

/* ---------------------------------------------------------------------- */
int CPhreeqcParser::get_line(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line from input file put in "line".
 *   Copy of input line is stored in "line_save".
 *   Characters after # are discarded in line but retained in "line_save"
 *
 *   Arguments:
 *      fp is file name
 *   Returns:
 *      EMPTY,
 *      EOF,
 *      KEYWORD,
 *      OK,
 *      OPTION
 */
 	int i, j, return_value, empty, l;
	char *ptr;
	char token[MAX_LENGTH];

	return_value = EMPTY;
	while (return_value == EMPTY) {
/*
 *   Eliminate all characters after # sign as a comment
 */
		i=-1;
		j=0;
		empty=TRUE;
/*
 *   Get line, check for eof
 */
		if (this->get_logical_line(&l) == EOF) {
			phreeqc::next_keyword=0;
			return (EOF);
		}
/*
 *   Get long lines
 */
		j = l;
		ptr = ::strchr(phreeqc::line_save, '#');
		if (ptr != NULL) {
			j = (int)(ptr - phreeqc::line_save);
		}
		::strncpy(phreeqc::line, phreeqc::line_save, (unsigned) j);
		phreeqc::line[j] = '\0';
		for (i = 0; i < j; i++) {
			if (!isspace((int)phreeqc::line[i]) ) {
				empty = FALSE;
				break;
			}
		}
/*
 *   New line character encountered
 */

		if (empty == TRUE) {
			return_value=EMPTY;
		} else {
			return_value=OK;
		}
	}
/*
 *   Determine return_value
 */
	if (return_value == OK) {
		if ( phreeqc::check_key(phreeqc::line) == TRUE) {
			return_value=KEYWORD;
		} else {
			ptr = phreeqc::line;
			phreeqc::copy_token(token, &ptr, &i);
			if (token[0] == '-' && isalpha((int)token[1])) {
				return_value = OPTION;
			}
		}
	}
	return (return_value);
}

