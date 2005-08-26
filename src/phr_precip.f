!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_precip(soln_id,count,aspecies,aconc)
      IMPLICIT NONE
      INCLUDE       '../IPhreeqc/include/IPhreeqc.f.inc'
      INTEGER       soln_id
      INTEGER       count
      CHARACTER*(*) aspecies(*)
      REAL*8        aconc(*)
      CHARACTER(80) line
      INTEGER       iresult
      INTEGER       i
      INTEGER       phr_precip

      WRITE (line,100),'SOLUTION ', soln_id
      iresult = AccumulateLine(line)

      WRITE (line,110),'-units  mol/kgw'
      iresult = AccumulateLine(line)

      DO 10 i=1,count
        WRITE (line,120), aspecies(i), aconc(i)
        iresult = AccumulateLine(line)
10    CONTINUE

      !!!CALL OutputLines
      phr_precip = Run(.FALSE., .FALSE., .FALSE., .FALSE.)

100   FORMAT(A,I10)
110   FORMAT(A)
120   FORMAT(A,1PG15.7E2)
      END FUNCTION phr_precip
