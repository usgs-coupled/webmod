#include <cstdio>  // FILE
#include "Debug.h" // ASSERT


extern "C" {
int get_line(FILE *fp);
int error_checking_gui(void);
void swig_open_selected_output(const char* filename);
}

#include "Phreeqc.hxx"
#include "PhreeqcParser.hxx"


/* ---------------------------------------------------------------------- */
int get_line(FILE *fp)
/* ---------------------------------------------------------------------- */
{
	ASSERT(fp == NULL);
	return CPhreeqc::singleton.GetLine();
}

/* ---------------------------------------------------------------------- */
int error_checking_gui(void)
/* ---------------------------------------------------------------------- */
{
	return 1;
}

/* ---------------------------------------------------------------------- */
void swig_open_selected_output(const char* filename)
/* ---------------------------------------------------------------------- */
{
	CPhreeqc::singleton.OpenSelectedOutput(filename);
}
