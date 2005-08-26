#include "TestVar.h"


TestVar::TestVar()
{
}

TestVar::~TestVar()
{
}

void TestVar::TestVarInit()
{
	VAR v;
	::VarInit(&v);
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
}
