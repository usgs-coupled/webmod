#pragma once

#include "IPhreeqc.hpp"
#include "Phreeqc.h"

#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

struct MixVars
{
  int          id;
  int          count;
  int         *solutions;
  double      *fracs;
  int          index_conserv;
  double       fill_factor;
  int         *n_user;
  double       rxnmols;
  double       tempc;
  double       tsec;
  double      *array;
  int	       row_dim;
  int	       col_dim;
  int          files_on;
  int          orig_out;
};


class IPQ_DLL_EXPORT IPhreeqcMMS : public IPhreeqc
{
public:
	IPhreeqcMMS(void);
	virtual ~IPhreeqcMMS(void);

public:
// COMMENT: {6/6/2012 5:06:16 PM}	void SetCallbackCookie(void* c);
// COMMENT: {6/6/2012 5:06:16 PM}	void SetPostRunCallback(PFN_POSTRUN_CALLBACK pfn);
// COMMENT: {6/6/2012 5:06:16 PM}	void SetPreRunCallback(PFN_PRERUN_CALLBACK pfn);
// COMMENT: {6/4/2012 5:27:58 PM}	int CatchErrors(PFN_CATCH_CALLBACK pfn, void *cookie);
// COMMENT: {6/4/2012 5:27:58 PM}
// COMMENT: {6/4/2012 5:27:58 PM}	typedef int (IPhreeqcMMS::*IPhreeqcMMSMemFunc)(void);  
// COMMENT: {6/4/2012 5:27:58 PM}	typedef int (Phreeqc::*PhreeqcMemFunc)(void);  
// COMMENT: {6/4/2012 5:27:58 PM}	
// COMMENT: {6/4/2012 5:27:58 PM}	int Catch(IPhreeqcMMSMemFunc pmfn);
// COMMENT: {6/4/2012 5:27:58 PM}
// COMMENT: {6/4/2012 5:27:58 PM}	int Catch(PhreeqcMemFunc pmfn);
// COMMENT: {6/4/2012 5:27:58 PM}
// COMMENT: {6/4/2012 5:27:58 PM}	static int BuildTallyCallback(void *cookie);
	int CatchBuildTally(void);
	int CatchGetTallyRowsColumns(int *rows, int *columns);
	int CatchGetTallyColumnHeading(int column, int *type, char *string, unsigned int string_length);
	int CatchGetTallyRowHeading(int row, char *string, unsigned int string_length);

	/*static*/
	int PreMixCallback(struct MixVars* pvars);
	/*static*/
	int PostMixCallback(struct MixVars* pvars);

// COMMENT: {6/6/2012 5:01:22 PM}protected:
// COMMENT: {6/6/2012 5:01:22 PM}	virtual void do_run(const char* sz_routine, std::istream* pis, PFN_PRERUN_CALLBACK pfn_pre, PFN_POSTRUN_CALLBACK pfn_post, void *cookie);

protected:
// COMMENT: {6/6/2012 5:05:58 PM}	PFN_PRERUN_CALLBACK pfn_pre;
// COMMENT: {6/6/2012 5:05:58 PM}	PFN_POSTRUN_CALLBACK pfn_post;
// COMMENT: {6/6/2012 5:05:58 PM}	void *cookie;

	friend class IPhreeqcMMSLib;
};

// COMMENT: {6/6/2012 5:13:01 PM}struct MixVars
// COMMENT: {6/6/2012 5:13:01 PM}{
// COMMENT: {6/6/2012 5:13:01 PM}  IPhreeqcMMS* ptr;
// COMMENT: {6/6/2012 5:13:01 PM}  int          id;
// COMMENT: {6/6/2012 5:13:01 PM}  int          count;
// COMMENT: {6/6/2012 5:13:01 PM}  int         *solutions;
// COMMENT: {6/6/2012 5:13:01 PM}  double      *fracs;
// COMMENT: {6/6/2012 5:13:01 PM}  int          index_conserv;
// COMMENT: {6/6/2012 5:13:01 PM}  double       fill_factor;
// COMMENT: {6/6/2012 5:13:01 PM}  int         *n_user;
// COMMENT: {6/6/2012 5:13:01 PM}  double       rxnmols;
// COMMENT: {6/6/2012 5:13:01 PM}  double       tempc;
// COMMENT: {6/6/2012 5:13:01 PM}  double       tsec;
// COMMENT: {6/6/2012 5:13:01 PM}  double      *array;
// COMMENT: {6/6/2012 5:13:01 PM}  int	       row_dim;
// COMMENT: {6/6/2012 5:13:01 PM}  int	       col_dim;
// COMMENT: {6/6/2012 5:13:01 PM}  int          files_on;
// COMMENT: {6/6/2012 5:13:01 PM}};
