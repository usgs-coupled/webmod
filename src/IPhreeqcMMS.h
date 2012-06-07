#pragma once

// COMMENT: {6/4/2012 9:51:11 PM}#include "IPhreeqc.hpp"
// COMMENT: {6/4/2012 9:51:11 PM}#include "Phreeqc.h"
#include "IPhreeqc.h"      // IPQ_RESULT

#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif


#if defined(__cplusplus)
extern "C" {
#endif

	int        CreateIPhreeqcMMS(void);
	IPQ_RESULT DestroyIPhreeqcMMS(int id);

#if defined(__cplusplus)
}
#endif

