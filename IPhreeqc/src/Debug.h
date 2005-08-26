#ifndef _INC_DEBUG_H
#define _INC_DEBUG_H

#if defined(_WIN32)
#include <crtdbg.h>
#define ASSERT _ASSERT
#else
#include <cassert>
#define ASSERT assert
#endif

#endif  /*_INC_DEBUG_H*/
