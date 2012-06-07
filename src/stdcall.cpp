#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

#include "fortran.h"

#if defined(__cplusplus)
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////
// CreateIPhreeqcMMS
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall CREATEIPHREEQCMMS(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall CREATEIPHREEQCMMS_(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall createiphreeqcmms(void)
{
	return CreateIPhreeqcMMSF();
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall createiphreeqcmms_(void)
{
	return CreateIPhreeqcMMSF();
}

///////////////////////////////////////////////////////////////////////////////
// DestroyIPhreeqcMMS
///////////////////////////////////////////////////////////////////////////////

// /iface:cvf
// /iface:stdcall /names:uppercase
IPQ_DLL_EXPORT int __stdcall DESTROYIPHREEQCMMS(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:cvf
// /iface:stdcall /names:uppercase /assume:underscore
IPQ_DLL_EXPORT int __stdcall DESTROYIPHREEQCMMS_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:stdcall
IPQ_DLL_EXPORT int __stdcall destroyiphreeqcmms(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}

// /iface:stdcall /assume:underscore
IPQ_DLL_EXPORT int __stdcall destroyiphreeqcmms_(int *id)
{
	return DestroyIPhreeqcMMSF(id);
}


#if defined(__cplusplus)
}
#endif
