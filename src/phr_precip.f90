!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_precip(id,soln_id,count,aspecies,aconc,tempc,pH)
      USE IPhreeqc
      IMPLICIT NONE
      !INCLUDE       '../IPhreeqc/src/IPhreeqc.f90.inc'
      INTEGER       id
      INTEGER       soln_id
      INTEGER       count
      CHARACTER(*)  aspecies(*)
      DOUBLE PRECISION        aconc(*)
      DOUBLE PRECISION        tempc
      DOUBLE PRECISION        pH
      CHARACTER(80) line
      INTEGER       iresult
      INTEGER       i
      INTEGER       phr_precip

      WRITE (line,100) 'SOLUTION ', soln_id
      iresult = AccumulateLine(id, line)

      WRITE (line,110) '-units  mol/kgw'
      iresult = AccumulateLine(id, line)

      WRITE (line,120) '-temp ', tempc
      iresult = AccumulateLine(id, line)

      WRITE (line,120) '-pH ', pH
      iresult = AccumulateLine(id, line)

      DO 10 i=1,count
        WRITE (line,120)  aspecies(i), aconc(i)
        iresult = AccumulateLine(id, line)
10    CONTINUE

      phr_precip = RunAccumulated(id)

100   FORMAT(A,I10)
110   FORMAT(A)
120   FORMAT(A,1PG15.7E2)
      END FUNCTION phr_precip
