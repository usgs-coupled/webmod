!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_multicopy(id, keyword, srcarray, targetarray, count)
      IMPLICIT NONE
      INCLUDE       '../IPhreeqc/src/IPhreeqc.f.inc'
      INTEGER       id
      CHARACTER(*)  keyword
      INTEGER       srcarray(*)
      INTEGER       targetarray(*)
      INTEGER       count
      INTEGER       iresult
      INTEGER       i
      CHARACTER(80) line
      INTEGER       phr_multicopy

      phr_multicopy = 0
      DO 10 i=1,count
        IF (srcarray(i).GE.0.AND.targetarray(i).GE.0) THEN

          WRITE (line,100),'COPY ', keyword, srcarray(i), targetarray(i)
          phr_multicopy = AccumulateLine(id, line)
          IF (phr_multicopy.NE.0) THEN
            RETURN
          ENDIF

          phr_multicopy = RunAccumulated(id)
          IF (phr_multicopy.NE.0) THEN
            !!!CALL OutputErrorString(id)
            RETURN
          ENDIF
       ENDIF
10    CONTINUE
100   FORMAT(A,A,1X,2I10)
      END FUNCTION phr_multicopy
