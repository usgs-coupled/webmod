!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM DRIVER
      IMPLICIT NONE
      INCLUDE '../IPhreeqc/include/IPhreeqc.f.inc'
      INCLUDE '../include/mms_phreeqc.inc'
      INTEGER iresult
      INTEGER      rows
      INTEGER      cols
      INTEGER      i
      INTEGER      j
      INTEGER      p
      INTEGER      len
      INTEGER      src(3)
      INTEGER      dest(3)
      CHARACTER(2) species(2)
      REAL*8       concs(2)
      INTEGER      solutions(4)
      REAL*8       fracs(4)
      REAL*8       conc(10)

      REAL*8        dvalue
      CHARACTER(50) heading
      CHARACTER(50) svalue
      INTEGER       vtype
      INTEGER       etype

! COMMENT: {3/3/2004 9:19:15 PM}      REAL*8        conc_conserv(2)
      REAL*8        conc_conserv(20)

      INTEGER       n_user(10)
      REAL*8        rxnmols
      REAL*8        tempc
      REAL*8        tsec
      REAL*8        ph
      REAL*8        array(50,50)

      !! load phreeqc.dat

      iresult = LoadDatabase('phreeqc.dat')
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors loading database:'
        CALL OutputLastError
        STOP
      ENDIF


      iresult = RunFile('test.pqi', .TRUE., .TRUE., .TRUE., .TRUE.)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors running test.pqi:'
        CALL OutputLastError
        STOP
      ENDIF


      cols = GetSelectedOutputColumnCount()

      ! headings
      DO 20 j=1,cols
        iresult = GetSelectedOutputValue(0, j, vtype, dvalue, svalue)
        len = INDEX(svalue, ' ')
        PRINT 50, svalue(1:len-1), ACHAR(9)
20    CONTINUE
      PRINT *

      rows = GetSelectedOutputRowCount()
      ! values
      DO 40 i=1,rows
        DO 30 j=1,cols
          iresult = GetSelectedOutputValue(i, j, vtype, dvalue, svalue)
          IF (iresult.EQ.VR_OK) THEN
            IF (vtype.eq.TT_EMPTY) THEN
              PRINT 50, ' ', ACHAR(9)
            ELSEIF(vtype.eq.TT_DOUBLE) THEN
              PRINT 60, dvalue, ACHAR(9)
            ELSEIF(vtype.eq.TT_STRING) THEN
              len = INDEX(svalue, ' ')
              PRINT 50, svalue(1:len-1), ACHAR(9)
            ENDIF
          ELSE
            IF (iresult.eq.VR_INVALIDROW) THEN
              PRINT 50, 'INVROW', ACHAR(9)
            ELSEIF (iresult.eq.VR_INVALIDCOL) THEN
              PRINT 50, 'INVCOL', ACHAR(9)
            ELSE
              PRINT 50, 'ERROR', ACHAR(9)
            ENDIF
          ENDIF
30      CONTINUE
        PRINT *
40    CONTINUE

50    FORMAT(A15,A,$)      
60    FORMAT(1PG15.7E2,A,$)

      !!}}
      
      !! run solns.pqi creating solutions 1-3
      iresult = RunFile('solns.pqi', .FALSE., .FALSE., .FALSE., .TRUE.)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors running solns.pqi:'
        CALL OutputLastError
        STOP
      ENDIF     


      iresult = build_tally_table()
      IF (iresult.NE.1) THEN
        PRINT *, 'Errors calling build_tally_table:'
        CALL OutputLastError
        STOP
      ENDIF
      PRINT *, "build_tally_table OK"

      iresult = get_tally_table_rows_columns(rows, cols)
      IF (iresult.NE.1) THEN
        PRINT *, 'Errors calling get_tally_table_rows_columns:'
        CALL OutputLastError
        STOP
      ENDIF
      PRINT *, "get_tally_table_rows_columns OK"
      PRINT *, "get_tally_table_rows_columns rows = ", rows
      PRINT *, "get_tally_table_rows_columns cols = ", cols

      DO 65 i = 1,cols
        iresult = get_tally_table_column_heading(i, etype, heading)
        IF (iresult.NE.1) THEN
          PRINT *, 'Errors during get_tally_table_column_heading:'
          CALL OutputLastError
          STOP
        ENDIF
        IF (etype.EQ.ET_SOLUTION) THEN
          PRINT *, 'Entity type: ET_SOLUTION'
        ELSEIF (etype.EQ.ET_REACTION) THEN
          PRINT *, 'Entity type: ET_REACTION'
        ELSEIF (etype.EQ.ET_SURFACE) THEN
          PRINT *, 'Entity type: ET_SURFACE'
        ELSEIF (etype.EQ.ET_GAS_PHASE) THEN
          PRINT *, 'Entity type: ET_GAS_PHASE'
        ELSEIF (etype.EQ.ET_PURE_PHASE) THEN
          PRINT *, 'Entity type: ET_PURE_PHASE'
        ELSEIF (etype.EQ.ET_SS_PHASE) THEN
          PRINT *, 'Entity type: ET_SS_PHASE'
        ELSEIF (etype.EQ.ET_KINETICS) THEN
          PRINT *, 'Entity type: ET_KINETICS'
        ELSEIF (etype.EQ.ET_MIX) THEN
          PRINT *, 'Entity type: ET_MIX'
        ELSEIF (etype.EQ.ET_TEMPERATURE) THEN
          PRINT *, 'Entity type: ET_TEMPERATURE'
        ELSE
          PRINT *, 'Entity type: ET_UNKNOWN'
        ENDIF
        len = INDEX(heading, ' ')
        WRITE(*, '(A,I1,A,A)'),'col_head(',i,')=',heading(1:len-1)
65    CONTINUE
      DO 68 i = 1,rows
        iresult = get_tally_table_row_heading(i, heading)
        IF (iresult.NE.1) THEN
          PRINT *, 'Errors during get_tally_table_row_heading:'
          CALL OutputLastError
          STOP
        ENDIF
        len = INDEX(heading, ' ')
        WRITE(*, '(A,I1,A,A)'),'row_head(',i,')=',heading(1:len-1)
68    CONTINUE



      !! copy solns 1-3 to 11-13
      DO 70 i=1,3
        src(i) = i
        dest(i) = i + 10
70    CONTINUE

      !! create solution 5
      species(1) = 'Cl'
      species(2) = 'Na'
      concs(1) = 1.0E-5
      concs(2) = 2.0E-5
      iresult = phr_precip(5, 2, species, concs)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors during phr_precip:'
        CALL OutputLastError
        STOP
      ENDIF

      do 1000 p=1,10

      iresult = phr_multicopy('solution', src, dest, 3)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors during phr_multicopy:'
        CALL OutputLastError
        STOP
      ENDIF

      
      !! create solution 100 and 200(conserv)
      solutions(1) = 1;
      solutions(2) = 2;
      solutions(3) = 3;
      solutions(4) = 5;
      fracs(1) = 0.20
      fracs(2) = 0.23
      fracs(3) = 0.27
      fracs(4) = 0.30
      
      n_user(ET_SOLUTION)    = -1
      n_user(ET_REACTION)    = 0
      n_user(ET_EXCHANGE)    = -1
      n_user(ET_SURFACE)     = -1
      n_user(ET_GAS_PHASE)   = -1
      n_user(ET_PURE_PHASE)  = -1
      n_user(ET_SS_PHASE)    = -1
      n_user(ET_KINETICS)    = -1
      n_user(ET_MIX)         = -1
      n_user(ET_TEMPERATURE) = -1

      iresult = AccumulateLine("SELECTED_OUTPUT")
      iresult = AccumulateLine("-reset false")
      iresult = AccumulateLine("-ph")
      iresult = AccumulateLine("-temperature")
      iresult = AccumulateLine("-totals Cl K Na")

      iresult = phr_mix(4,solutions,fracs,200,
     &                1.5d0,100,conc_conserv,.TRUE.,
     &                n_user,rxnmols,tempc,ph,tsec,array,
     &                50,50)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors during phr_mix:'
        CALL OutputLastError
        STOP
      ENDIF
1000  CONTINUE
      !{{
      PRINT *, 'pH:', ph
      PRINT *, 'Temp:', tempc
      !}}
      
      !! output solution 100
      cols = GetSelectedOutputColumnCount()
      PRINT *, 'Solution 100'
      DO 80 i = 3,cols
        iresult = GetSelectedOutputValue(0, i, vtype, dvalue, heading)
        IF (iresult.NE.0) THEN
          PRINT *, 'Errors during GetSelectedOutputValue:'
          CALL OutputLastError
          STOP
        ENDIF
        len = INDEX(heading, ' ')
        PRINT *, heading(1:len-1), ACHAR(9), conc_conserv(i-2)
80    CONTINUE

! COMMENT: {3/3/2004 9:32:11 PM}      iresult = build_tally_table()
! COMMENT: {3/3/2004 9:32:11 PM}      IF (iresult.NE.1) THEN
! COMMENT: {3/3/2004 9:32:11 PM}        PRINT *, 'Errors calling build_tally_table:'
! COMMENT: {3/3/2004 9:32:11 PM}        CALL OutputLastError
! COMMENT: {3/3/2004 9:32:11 PM}        STOP
! COMMENT: {3/3/2004 9:32:11 PM}      ENDIF
! COMMENT: {3/3/2004 9:32:11 PM}      PRINT *, "build_tally_table OK"
! COMMENT: {3/3/2004 9:32:11 PM}
! COMMENT: {3/3/2004 9:32:11 PM}
! COMMENT: {3/3/2004 9:32:11 PM}      iresult = get_tally_table_rows_columns(rows, cols)
! COMMENT: {3/3/2004 9:32:11 PM}      IF (iresult.NE.1) THEN
! COMMENT: {3/3/2004 9:32:11 PM}        PRINT *, 'Errors calling get_tally_table_rows_columns:'
! COMMENT: {3/3/2004 9:32:11 PM}        CALL OutputLastError
! COMMENT: {3/3/2004 9:32:11 PM}        STOP
! COMMENT: {3/3/2004 9:32:11 PM}      ENDIF
! COMMENT: {3/3/2004 9:32:11 PM}      PRINT *, "get_tally_table_rows_columns OK"
! COMMENT: {3/3/2004 9:32:11 PM}      PRINT *, "get_tally_table_rows_columns rows = ", rows
! COMMENT: {3/3/2004 9:32:11 PM}      PRINT *, "get_tally_table_rows_columns cols = ", cols
! COMMENT: {3/3/2004 9:32:11 PM}
! COMMENT: {3/3/2004 9:32:11 PM}
! COMMENT: {3/3/2004 9:32:11 PM}
! COMMENT: {3/3/2004 9:32:11 PM}      DO 85 i = 1,cols
! COMMENT: {3/3/2004 9:32:11 PM}        iresult = get_tally_table_column_heading(i, etype, heading)
! COMMENT: {3/3/2004 9:32:11 PM}        IF (iresult.NE.1) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Errors during get_tally_table_column_heading:'
! COMMENT: {3/3/2004 9:32:11 PM}          CALL OutputLastError
! COMMENT: {3/3/2004 9:32:11 PM}          STOP
! COMMENT: {3/3/2004 9:32:11 PM}        ENDIF
! COMMENT: {3/3/2004 9:32:11 PM}        IF (etype.EQ.ET_SOLUTION) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_SOLUTION'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_REACTION) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_REACTION'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_SURFACE) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_SURFACE'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_GAS_PHASE) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_GAS_PHASE'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_PURE_PHASE) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_PURE_PHASE'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_SS_PHASE) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_SS_PHASE'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_KINETICS) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_KINETICS'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_MIX) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_MIX'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSEIF (etype.EQ.ET_TEMPERATURE) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_TEMPERATURE'
! COMMENT: {3/3/2004 9:32:11 PM}        ELSE
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Entity type: ET_UNKNOWN'
! COMMENT: {3/3/2004 9:32:11 PM}        ENDIF
! COMMENT: {3/3/2004 9:32:11 PM}        len = INDEX(heading, ' ')
! COMMENT: {3/3/2004 9:32:11 PM}        WRITE(*, '(A,I1,A,A)'),'col_head(',i,')=',heading(1:len-1)
! COMMENT: {3/3/2004 9:32:11 PM}85    CONTINUE
! COMMENT: {3/3/2004 9:32:11 PM}      DO 88 i = 1,rows
! COMMENT: {3/3/2004 9:32:11 PM}        iresult = get_tally_table_row_heading(i, heading)
! COMMENT: {3/3/2004 9:32:11 PM}        IF (iresult.NE.1) THEN
! COMMENT: {3/3/2004 9:32:11 PM}          PRINT *, 'Errors during get_tally_table_row_heading:'
! COMMENT: {3/3/2004 9:32:11 PM}          CALL OutputLastError
! COMMENT: {3/3/2004 9:32:11 PM}          STOP
! COMMENT: {3/3/2004 9:32:11 PM}        ENDIF
! COMMENT: {3/3/2004 9:32:11 PM}        len = INDEX(heading, ' ')
! COMMENT: {3/3/2004 9:32:11 PM}        WRITE(*, '(A,I1,A,A)'),'row_head(',i,')=',heading(1:len-1)
! COMMENT: {3/3/2004 9:32:11 PM}88    CONTINUE

      
      PRINT *, 'All OK'      

      END PROGRAM DRIVER
