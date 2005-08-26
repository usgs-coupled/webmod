#include "TestSelectedOutput.h"


// TODO remove these dependencies
extern "C" {
int warning_msg (const char *err_str);
extern int n_user_punch_index;
extern char **user_punch_headings;
extern int user_punch_count_headings;
}

extern "C" {
int EndRow(void);
}


int 
warning_msg(const char *err_str)
{
	return 0;
}

int n_user_punch_index;
char **user_punch_headings;
int user_punch_count_headings;


TestSelectedOutput::TestSelectedOutput()
{
}

TestSelectedOutput::~TestSelectedOutput()
{
}

void
TestSelectedOutput::TestEmpty()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);
}

void
TestSelectedOutput::TestSinglePushBack()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CVar v(7.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("pH", v) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	// row count doesn't change until EndRow is called
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestSinglePushBack");
#endif
}

void
TestSelectedOutput::TestMultiplePushBack()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CVar v1(7.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("pH", v1) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v2(8.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("pH", v2) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestMultiplePushBack");
#endif
}

void 
TestSelectedOutput::TestNewHeadingsPushBack()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CVar v1(7.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("pH", v1) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v2(8.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("pH", v2) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v3(9.0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBack("user_pH", v3) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 2);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);


	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 2);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestNewHeadingsPushBack");
#endif
}

void
TestSelectedOutput::TestPushBackDouble()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 7.0) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "pH") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_DOUBLE);
	CPPUNIT_ASSERT(vval.dVal == 7.0);
}

void
TestSelectedOutput::TestPushBackLong()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("Sim", 2) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading plus first row

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "Sim") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_LONG);
	CPPUNIT_ASSERT(vval.lVal == 2);
}

void
TestSelectedOutput::TestPushBackString()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackString("state", "i_soln") == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "state") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(vval.sVal, "i_soln") == 0);
}

void
TestSelectedOutput::TestPushBackEmpty()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("Empty") == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "Empty") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
}

void
TestSelectedOutput::TestDuplicateHeadings()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 7.0) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	// overwrite pH with 8.0
	//
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 8.0) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "pH") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_DOUBLE);
	CPPUNIT_ASSERT(vval.dVal == 8.0);
}

void
TestSelectedOutput::TestEndRow()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 7.0) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 1); // heading

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CVar v;
	CPPUNIT_ASSERT(v.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "pH") == 0);

	CVar vval;
	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
	CPPUNIT_ASSERT(vval.type == TT_DOUBLE);
	CPPUNIT_ASSERT(vval.dVal == 7.0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 8.0) == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);

	CVar vval3;
	CPPUNIT_ASSERT(vval3.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval3) == VR_OK);
	CPPUNIT_ASSERT(vval3.type == TT_DOUBLE);
	CPPUNIT_ASSERT(vval3.dVal == 7.0);

	CVar vval2;
	CPPUNIT_ASSERT(vval2.type == TT_EMPTY);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(2, 0, &vval2) == VR_OK);
	CPPUNIT_ASSERT(vval2.type == TT_DOUBLE);
	CPPUNIT_ASSERT(vval2.dVal == 8.0);
}

void
TestSelectedOutput::TestEndRow2()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 6.0) == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 7.0) == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 8.0) == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 9.0) == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

// COMMENT: {12/8/2003 6:06:03 PM}	CVar v;
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(v.type == TT_EMPTY);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(v.type == TT_STRING);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(::strcmp(v.sVal, "pH") == 0);
// COMMENT: {12/8/2003 6:06:03 PM}
// COMMENT: {12/8/2003 6:06:03 PM}	CVar vval;
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval.type == TT_EMPTY);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval) == VR_OK);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval.type == TT_DOUBLE);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval.dVal == 7.0);
// COMMENT: {12/8/2003 6:06:03 PM}
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackDouble("pH", 8.0) == 0);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);
// COMMENT: {12/8/2003 6:06:03 PM}
// COMMENT: {12/8/2003 6:06:03 PM}	CVar vval3;
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval3.type == TT_EMPTY);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &vval3) == VR_OK);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval3.type == TT_DOUBLE);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval3.dVal == 7.0);
// COMMENT: {12/8/2003 6:06:03 PM}
// COMMENT: {12/8/2003 6:06:03 PM}	CVar vval2;
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval2.type == TT_EMPTY);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(2, 0, &vval2) == VR_OK);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval2.type == TT_DOUBLE);
// COMMENT: {12/8/2003 6:06:03 PM}	CPPUNIT_ASSERT(vval2.dVal == 8.0);
}


void
TestSelectedOutput::TestTooManyHeadings()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	// USER_PUNCH
	// -headings 1.name 1.type 1.moles

	n_user_punch_index  = 0;
	user_punch_headings = NULL;
	user_punch_count_headings = 0;

	user_punch_headings = (char**)::realloc(user_punch_headings, (size_t) (user_punch_count_headings + 1) * sizeof(char *));
	user_punch_headings[user_punch_count_headings] = ::strdup("1.name");
	user_punch_count_headings++;

	user_punch_headings = (char**)::realloc(user_punch_headings, (size_t) (user_punch_count_headings + 1) * sizeof(char *));
	user_punch_headings[user_punch_count_headings] = ::strdup("1.type");
	user_punch_count_headings++;

	user_punch_headings = (char**)::realloc(user_punch_headings, (size_t) (user_punch_count_headings + 1) * sizeof(char *));
	user_punch_headings[user_punch_count_headings] = ::strdup("1.moles");
	user_punch_count_headings++;

	CPPUNIT_ASSERT(::EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 3);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestTooManyHeadings");
#endif

	// clean up headings
	//
	for (int i = 0; i < user_punch_count_headings; ++i) {
		::free(user_punch_headings[i]);
	}
	::free(user_punch_headings);
	user_punch_headings = NULL;
	user_punch_count_headings = 0;

	CVar head0, head1, head2;
	CVar val0, val1, val2;

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &head0) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 1, &head1) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 2, &head2) == VR_OK);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &val0) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 1, &val1) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 2, &val2) == VR_OK);

	CPPUNIT_ASSERT(head0.type == TT_STRING);
	CPPUNIT_ASSERT(head1.type == TT_STRING);
	CPPUNIT_ASSERT(head2.type == TT_STRING);

	CPPUNIT_ASSERT(val0.type == TT_EMPTY);
	CPPUNIT_ASSERT(val1.type == TT_EMPTY);
	CPPUNIT_ASSERT(val2.type == TT_EMPTY);

	CPPUNIT_ASSERT(::strcmp(head0.sVal, "1.name") == 0);
	CPPUNIT_ASSERT(::strcmp(head1.sVal, "1.type") == 0);
	CPPUNIT_ASSERT(::strcmp(head2.sVal, "1.moles") == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("sim", 1)            == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackString("state", "i_soln") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("soln", 22)          == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 6);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestTooManyHeadings");
#endif
}

void
TestSelectedOutput::TestNotEnoughHeadings()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	// USER_PUNCH
	// -headings 1.name 1.type 1.moles

	n_user_punch_index  = 0;
	user_punch_headings = NULL;
	user_punch_count_headings = 0;

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("sim", 1)            == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackString("state", "i_soln") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("soln", 22)          == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow()      == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 3);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);
#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestNotEnoughHeadings");
#endif

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("sim", 2)           == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackString("state", "react") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackLong("soln", 23)         == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("no_heading_1") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("no_heading_2") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("no_heading_3") == 0);

#if defined(_DEBUG)
	CSelectedOutput::singleton.Dump("TestNotEnoughHeadings");
#endif

	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 6);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 3);

	CVar head0, head1, head2, head3, head4, head5;
	CVar val0, val1, val2, val3, val4, val5;

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &head0) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 1, &head1) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 2, &head2) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 3, &head3) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 4, &head4) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 5, &head5) == VR_OK);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &val0) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 1, &val1) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 2, &val2) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 3, &val3) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 4, &val4) == VR_OK);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 5, &val5) == VR_OK);

	CPPUNIT_ASSERT(head0.type == TT_STRING);
	CPPUNIT_ASSERT(head1.type == TT_STRING);
	CPPUNIT_ASSERT(head2.type == TT_STRING);
	CPPUNIT_ASSERT(head3.type == TT_STRING);
	CPPUNIT_ASSERT(head4.type == TT_STRING);
	CPPUNIT_ASSERT(head5.type == TT_STRING);

	CPPUNIT_ASSERT(val0.type == TT_LONG);
	CPPUNIT_ASSERT(val1.type == TT_STRING);
	CPPUNIT_ASSERT(val2.type == TT_LONG);
	CPPUNIT_ASSERT(val3.type == TT_EMPTY);
	CPPUNIT_ASSERT(val4.type == TT_EMPTY);
	CPPUNIT_ASSERT(val5.type == TT_EMPTY);

	CPPUNIT_ASSERT(::strcmp(head0.sVal, "sim") == 0);
	CPPUNIT_ASSERT(::strcmp(head1.sVal, "state") == 0);
	CPPUNIT_ASSERT(::strcmp(head2.sVal, "soln") == 0);
	CPPUNIT_ASSERT(::strcmp(head3.sVal, "no_heading_1") == 0);
	CPPUNIT_ASSERT(::strcmp(head4.sVal, "no_heading_2") == 0);
	CPPUNIT_ASSERT(::strcmp(head5.sVal, "no_heading_3") == 0);

	CPPUNIT_ASSERT(val0.lVal == 1);
	CPPUNIT_ASSERT(::strcmp(val1.sVal, "i_soln") == 0);
	CPPUNIT_ASSERT(val2.lVal == 22);
}

void
TestSelectedOutput::TestInvalidRow()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CVar v;
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_INVALIDROW);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDROW);


	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(-1, -1, &v) == VR_INVALIDROW);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDROW);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("heading") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "heading") == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_EMPTY);


	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(2, 0, &v) == VR_INVALIDROW);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDROW);
}

void
TestSelectedOutput::TestInvalidCol()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CVar v;
	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_INVALIDROW);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDROW);


	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(-1, -1, &v) == VR_INVALIDROW);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDROW);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("heading") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v.sVal, "heading") == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(1, 0, &v) == VR_OK);
	CPPUNIT_ASSERT(v.type == TT_EMPTY);


	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, 1, &v) == VR_INVALIDCOL);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDCOL);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.Get(0, -1, &v) == VR_INVALIDCOL);
	CPPUNIT_ASSERT(v.type == TT_ERROR);
	CPPUNIT_ASSERT(v.vresult == VR_INVALIDCOL);
}

void
TestSelectedOutput::TestGet()
{
	CSelectedOutput::singleton.Clear();
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 0);

	CPPUNIT_ASSERT(CSelectedOutput::singleton.PushBackEmpty("heading") == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.EndRow() == 0);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetColCount() == 1);
	CPPUNIT_ASSERT(CSelectedOutput::singleton.GetRowCount() == 2);


	CVar v0 = CSelectedOutput::singleton.Get(0, 0);
	CPPUNIT_ASSERT(v0.type == TT_STRING);
	CPPUNIT_ASSERT(::strcmp(v0.sVal, "heading") == 0);

	CVar v1 = CSelectedOutput::singleton.Get(1, 0);
	CPPUNIT_ASSERT(v1.type == TT_EMPTY);
}
