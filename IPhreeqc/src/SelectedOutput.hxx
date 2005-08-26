// SelectedOutput.h: interface for the CSelectedOutput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined _INC_SELECTEDOUTPUT_H
#define _INC_SELECTEDOUTPUT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <map>
#include <list>
#include <vector>
#include "CVar.hxx"

// TODO: templatize
class CSelectedOutput  
{
protected:
	CSelectedOutput(void);
public:
	static CSelectedOutput singleton;
	virtual ~CSelectedOutput(void);

	int EndRow(void);
	void Clear(void);

	size_t GetRowCount(void)const;
	size_t GetColCount(void)const;

	CVar Get(int nRow, int nCol)const;
	VRESULT CSelectedOutput::Get(int nRow, int nCol, VAR* pVAR)const;

	int PushBack(const char* key, const CVar& var);

	int PushBackDouble(const char* key, double dVal);
	int PushBackLong(const char* key, long lVal);
	int PushBackString(const char* key, const char* sVal);
	int PushBackEmpty(const char* key);

#if defined(_DEBUG)
	void Dump(const char* heading);
	void AssertValid(void)const;
#endif

protected:
	friend std::ostream& operator<< (std::ostream &os, const CSelectedOutput &a);
// COMMENT: {12/8/2003 6:20:55 PM}	size_t m_nColIndex;
// COMMENT: {12/8/2003 6:22:17 PM}	size_t m_nRowIndex;

// COMMENT: {12/8/2003 6:18:00 PM}	size_t m_nColCount;
	size_t m_nRowCount;

	//std::map< std::string, std::vector<CVar> > m_mapHeadingToVarVector;
	//std::map< size_t, std::string > m_mapColToHeading;
	std::vector< std::vector<CVar> > m_arrayVar;
	std::vector<CVar> m_vecVarHeadings;
	std::map< std::string, size_t > m_mapHeadingToCol;

	//std::map< std::string, size_t > m_mapHeadingToColumn;
};

#endif // !defined(_INC_SELECTEDOUTPUT_H)
