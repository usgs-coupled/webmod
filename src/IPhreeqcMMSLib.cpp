#include "IPhreeqcMMSLib.h"

#include "IPhreeqcMMS.hpp"       // IPhreeqcMMS
#include "IPhreeqc.hpp"          // IPhreeqc
#include "IPhreeqc.h"            // IPQ_RESULT
#include "IPhreeqcMMS.h"


int
CreateIPhreeqcMMS(void)
{
	return IPhreeqcMMSLib::CreateIPhreeqcMMS();
}

// COMMENT: {6/4/2012 9:59:08 PM}int
// COMMENT: {6/4/2012 9:59:08 PM}CreateIPhreeqcMMSF(void)
// COMMENT: {6/4/2012 9:59:08 PM}{
// COMMENT: {6/4/2012 9:59:08 PM}	return CreateIPhreeqcMMS();
// COMMENT: {6/4/2012 9:59:08 PM}}


IPQ_RESULT
DestroyIPhreeqcMMS(int id)
{
	return IPhreeqcMMSLib::DestroyIPhreeqcMMS(id);
}

// COMMENT: {6/4/2012 9:59:13 PM}IPQ_RESULT
// COMMENT: {6/4/2012 9:59:13 PM}DestroyIPhreeqcMMSF(int *id)
// COMMENT: {6/4/2012 9:59:13 PM}{
// COMMENT: {6/4/2012 9:59:13 PM}	return DestroyIPhreeqcMMS(*id);
// COMMENT: {6/4/2012 9:59:13 PM}}

int
IPhreeqcMMSLib::CreateIPhreeqcMMS(void)
{
	int n = IPQ_OUTOFMEMORY;
	try
	{
		IPhreeqcMMS* IPhreeqcPtr = new IPhreeqcMMS;
		n = IPhreeqcPtr->Index;
	}
	catch(...)
	{
		return IPQ_OUTOFMEMORY;
	}
	return n;
}

IPQ_RESULT
IPhreeqcMMSLib::DestroyIPhreeqcMMS(int id)
{
	IPQ_RESULT retval = IPQ_BADINSTANCE;
	if (id >= 0)
	{
		if (IPhreeqcMMS *ptr = IPhreeqcMMSLib::GetInstance(id))
		{
			delete ptr;
			retval = IPQ_OK;
		}
	}
	return retval;
}

IPhreeqcMMS*
IPhreeqcMMSLib::GetInstance(int id)
{
	std::map<size_t, IPhreeqc*>::iterator it = IPhreeqc::Instances.find(size_t(id));
	if (it != IPhreeqc::Instances.end())
	{
		return dynamic_cast<IPhreeqcMMS*>((*it).second);
	}
	return 0;
}
