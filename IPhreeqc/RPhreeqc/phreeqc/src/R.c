#include <R.h>
#include <Rdefines.h>
#include "Var.h"
#include "IPhreeqc.h"

SEXP
read_db(SEXP filename)
{
  const char* name;

  /* check args */
  if (!isString(filename) || length(filename) != 1) {
    error("read_db:filename is not a single string\n");
  }

  name = CHAR(STRING_ELT(filename, 0));
  if (LoadDatabase(name) != VR_OK) {
    error("read_db:Database contains errors.\n");
  }

  return(R_NilValue);
}

SEXP
read(SEXP string_input)
{
  const char* str_in;

  /* check args */
  if (!isString(string_input) || length(string_input) != 1) {
    error("read:string_input is not a single string\n");
  }

  str_in = CHAR(STRING_ELT(string_input, 0));

  if (AccumulateLine(str_in) != VR_OK) {
    error("read: Out of memory\n");
  }

  return(R_NilValue);
}

SEXP
run(void)
{
  if (Run(0, 0, 0, 0) != VR_OK) {
    error("run: input contains errors\n");
  }
  return(R_NilValue);
}

SEXP
runFile(SEXP filename)
{
  const char* name;

  /* check args */
  if (!isString(filename) || length(filename) != 1) {
    error("runFile: filename is not a single string\n");
  }

  name = CHAR(STRING_ELT(filename, 0));
  if (RunFile(name, 0, 0, 0, 0) != VR_OK) {
    error("runFile: input contains errors\n");
  }

  return(R_NilValue);
}


SEXP
getColumnCount()
{
  SEXP cols = R_NilValue;
  PROTECT(cols = allocVector(INTSXP, 1));
  INTEGER(cols)[0] = GetSelectedOutputColumnCount();
  UNPROTECT(1);
  return cols;
}

SEXP
getRowCount()
{
  SEXP rows = R_NilValue;
  PROTECT(rows = allocVector(INTSXP, 1));
  INTEGER(rows)[0] = GetSelectedOutputRowCount();
  UNPROTECT(1);
  return rows;
}

SEXP
getCol(int ncol)
{
  int r;
  int cols;
  int rows;
  VAR vn;
  VAR vv;

  SEXP ans = R_NilValue;

  cols = GetSelectedOutputColumnCount();
  rows = GetSelectedOutputRowCount();
  if (cols == 0 || rows == 0) {
    //error("getColumn: no data\n");
    return ans;
  }

  VarInit(&vn);
  GetSelectedOutputValue(0, ncol, &vn);

  VarInit(&vv);
  GetSelectedOutputValue(1, ncol, &vv);

  switch (vv.type) {
  case TT_LONG:
    PROTECT(ans = allocVector(INTSXP, rows-1));
    for (r = 1; r < rows; ++r) {
      VarInit(&vv);
      GetSelectedOutputValue(r, ncol, &vv);
      if (vv.lVal == -99) {
	INTEGER(ans)[r-1] = NA_INTEGER;
      }
      else {
	INTEGER(ans)[r-1] = vv.lVal;
      }
      VarClear(&vv);
    }
    UNPROTECT(1);
    break;
  case TT_DOUBLE:
    PROTECT(ans = allocVector(REALSXP, rows-1));
    for (r = 1; r < rows; ++r) {
      VarInit(&vv);
      GetSelectedOutputValue(r, ncol, &vv);
      if (vv.dVal == -999.999) {
	REAL(ans)[r-1] = NA_REAL;
      }
      else {
	REAL(ans)[r-1] = vv.dVal;
      }
      VarClear(&vv);
    }
    UNPROTECT(1);
    break;
  case TT_STRING:
    PROTECT(ans = allocVector(STRSXP, rows-1));
    for (r = 1; r < rows; ++r) {
      VarInit(&vv);
      GetSelectedOutputValue(r, ncol, &vv);
      SET_STRING_ELT(ans, r-1, mkChar(vv.sVal));
      VarClear(&vv);
    }
    UNPROTECT(1);
    break;
  case TT_EMPTY:
    break;
  case TT_ERROR:
    break;
  }
  return ans;
}


SEXP
getColumn(SEXP column)
{
  int ncol;

  PROTECT(column = AS_INTEGER(column));
  ncol = INTEGER_POINTER(column)[0];
  UNPROTECT(1);

  return getCol(ncol);
}

#define CONVERT_TO_DATA_FRAME

SEXP
getSelOut(void)
{
  int r;
  int c;
  int cols;
  int rows;
  VAR vn;

  SEXP list;
  SEXP attr;
  SEXP col;

#if defined(CONVERT_TO_DATA_FRAME)
  SEXP class;
  SEXP row_names;
#endif

  list = R_NilValue;

  cols = GetSelectedOutputColumnCount();
  rows = GetSelectedOutputRowCount();
  if (cols == 0 || rows == 0) {
    return list;
  }


  PROTECT(list = allocVector(VECSXP, cols));
  PROTECT(attr = allocVector(STRSXP, cols));
  for (c = 0; c < cols; ++c) {

    VarInit(&vn);
    GetSelectedOutputValue(0, c, &vn);

    PROTECT(col = getCol(c));

    SET_VECTOR_ELT(list, c, col);
    SET_STRING_ELT(attr, c, mkChar(vn.sVal));

    VarClear(&vn);
  }

  setAttrib(list, R_NamesSymbol, attr);


#if defined(CONVERT_TO_DATA_FRAME)
  /* Turn the data "list" into a "data.frame" */
  /* see model.c */

  PROTECT(class = mkString("data.frame"));
  setAttrib(list, R_ClassSymbol, class);
  UNPROTECT(1);

  PROTECT(row_names = allocVector(INTSXP, rows-1));
  for (r = 0; r < rows-1; ++r) INTEGER(row_names)[r] = r+1;
  setAttrib(list, R_RowNamesSymbol, row_names);
  UNPROTECT(1);
#endif

  UNPROTECT(2+cols);
  return list;
}

SEXP
getLastErrorString()
{
	SEXP ans = R_NilValue;
    PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(GetLastErrorString()));
	UNPROTECT(1);
	return ans;
}
