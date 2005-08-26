#include "TestCVar.h"


TestCVar::TestCVar()
{
}

TestCVar::~TestCVar()
{
}

void TestCVar::TestCVarCtor()
{
	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
}
