!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_mix(id,count,solutions,fracs,index_conserv, &
                     fill_factor,index_rxn,conc_conserv,files_on, &
                     n_user,rxnmols,tempc,ph,tsec,array, &
                     arr_rows,arr_cols)
#if defined(_WIN32)
      USE IFPORT   ! to enable 'SYSTEM' calls
#endif
      IMPLICIT NONE
      INCLUDE       'IPhreeqc.f.inc'
      INTEGER       id              ! 
      INTEGER       count           ! solution count
      INTEGER       solutions(*)    ! solution #'s
      REAL*8        fracs(*)        ! mixing fractions
      INTEGER       index_conserv   ! 
      REAL*8        fill_factor     ! 
      INTEGER       index_rxn       ! 
      REAL*8        conc_conserv(*) ! 
      LOGICAL       files_on        ! write to files
      INTEGER       n_user(*)       !
      REAL*8        rxnmols         !
      REAL*8        tempc           !
      REAL*8        ph              !
      REAL*8        tsec            !
      REAL*8        array(*)        !
      INTEGER       arr_rows        ! 
      INTEGER       arr_cols        ! 


      CHARACTER(80) line
      INTEGER       i
      INTEGER       cols
      INTEGER       vtype
      REAL*8        dvalue
      INTEGER       rows
      INTEGER       phr_mix

      INTEGER RunMixF

! Rick's debug/
!      integer          startmix(3), endmix(3), nstep, j, et_hyd, et_mix
!      integer, external  ::  getstep, elapsed_time
!      integer, save      ::  et_hold
!      LOGICAL          fil_temp, step1
!
!      data step1/.true./
!      save step1
!
!      fil_temp = files_on
!      startmix(1) = elapsed_time(2)
!      nstep=getstep()  !  Serial day number of run
!
!      if(step1) then
!        et_hold = 0
!        step1=.false.
!      endif
!      et_hyd = startmix(1) - et_hold

!      if(index_rxn.gt.208000000.and.index_rxn.lt.209000000) files_on = .true.
! /debug
      phr_mix = RunMixF(id,count,solutions,fracs,index_conserv, &
                     fill_factor,index_rxn,conc_conserv,files_on, &
                     n_user,rxnmols,tempc,ph,tsec,array, &
                     arr_rows,arr_cols)


      rows = GetSelectedOutputRowCount(id)
      IF (rows.NE.1.AND.rows.NE.3) THEN
        PRINT *, 'ERROR: Expected rows = 1 or 3'
        PRINT *, 'In phr_mix rows = ', rows
        PRINT *, 'No selected_ouput defined?'
        RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(id, 0, 1, vtype, dvalue, line)
      IF (phr_mix.EQ.IPQ_OK) THEN
         if ('pH' .NE. line(1:2)) THEN
            PRINT *, 'ERROR: Expected pH got:' , line
            PRINT *, 'In phr_mix'
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(id, 1, 1, vtype, dvalue, line)
      IF (phr_mix.EQ.IPQ_OK) THEN
         IF (vtype.eq.TT_DOUBLE) THEN
            pH = dvalue
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(id, 0, 2, vtype, dvalue, line)
      IF (phr_mix.EQ.IPQ_OK) THEN
         if ('temp(C)' .NE. line(1:7)) THEN
            PRINT *, 'ERROR: Expected temp(C) got:' , line
            PRINT *, 'In phr_mix'
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(id, 1, 2, vtype, dvalue, line)
      IF (phr_mix.EQ.IPQ_OK) THEN
         IF (vtype.eq.TT_DOUBLE) THEN
            tempc = dvalue
         ENDIF
      ELSE
         RETURN
      ENDIF


      cols = GetSelectedOutputColumnCount(id)
      DO 20 i=3,cols
        phr_mix = GetSelectedOutputValue(id, 1, i, vtype, dvalue, line)
        IF (phr_mix.EQ.IPQ_OK) THEN
          IF (vtype.eq.TT_DOUBLE) THEN
            conc_conserv(i - 2) = dvalue
          ENDIF
        ELSE
          RETURN
        ENDIF
20    CONTINUE
! Debug
!      files_on = fil_temp
!      endmix(1) = elapsed_time(2)
!      et_mix = endmix(1) - startmix(1)  
!      write(25,1000)nstep, et_hyd, et_mix, startmix(1), endmix(1), &
!         ' Rxn soln: ',index_rxn, &
!         ' # of inputs/solns: ',count, (solutions(j),j=1,count)
!      if(index_rxn.gt.208000000.and.index_rxn.lt.209000000) then
!        I = SYSTEM('type .\phreeqc.0.out >> .\output\select_mixes')
!      endif
!1000  format(I5,' ', 4(I8.8," "),A, I10, A, I3, 20I10) 
!      et_hold = endmix(1)
! /Debug            

100   FORMAT(A)
110   FORMAT(TR4,I10,1PG15.7E2)
120   FORMAT(A,I10)
      END FUNCTION phr_mix
