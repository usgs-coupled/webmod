!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_multicopy(keyword, srcarray, targetarray, count)
      IMPLICIT NONE
      INCLUDE       '../IPhreeqc/include/IPhreeqc.f.inc'
      CHARACTER(*)  keyword
      INTEGER       srcarray(*)
      INTEGER       targetarray(*)
      INTEGER       count
      INTEGER       iresult
      INTEGER       i
      CHARACTER(80) line
      INTEGER       phr_multicopy

      DO 10 i=1,count
        IF (srcarray(i).GE.0.AND.targetarray(i).GE.0) THEN

          WRITE (line,100),'COPY ', keyword, srcarray(i), targetarray(i)
          iresult = AccumulateLine(line)

          !!!CALL OutputLines
          phr_multicopy = Run(.FALSE., .FALSE., .FALSE., .FALSE.)
          IF (phr_multicopy.NE.0) THEN
            !!!CALL OutputLastError
            RETURN
          ENDIF
        ENDIF
10    CONTINUE
100   FORMAT(A,A,1X,2I10)
      END FUNCTION phr_multicopy
