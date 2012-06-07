#pragma once

#include "IPhreeqc.h"      // IPQ_RESULT


// COMMENT: {6/4/2012 9:46:18 PM}#if defined(__cplusplus)
// COMMENT: {6/4/2012 9:46:18 PM}extern "C" {
// COMMENT: {6/4/2012 9:46:18 PM}#endif
// COMMENT: {6/4/2012 9:46:18 PM}
// COMMENT: {6/4/2012 9:46:18 PM}	int        CreateIPhreeqcMMS(int id);
// COMMENT: {6/4/2012 9:46:18 PM}	IPQ_RESULT DestroyIPhreeqcMMS(int id);
// COMMENT: {6/4/2012 9:46:18 PM}
// COMMENT: {6/4/2012 9:46:18 PM}#if defined(__cplusplus)
// COMMENT: {6/4/2012 9:46:18 PM}}
// COMMENT: {6/4/2012 9:46:18 PM}#endif

class IPhreeqcMMS;

class IPhreeqcMMSLib
{
public:
	static int CreateIPhreeqcMMS(void);
	static IPQ_RESULT DestroyIPhreeqcMMS(int n);
	static IPhreeqcMMS* GetInstance(int n);
};
