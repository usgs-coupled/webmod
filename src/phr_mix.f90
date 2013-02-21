!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_mix(id,count,solutions,fracs,index_conserv, &
                     fill_factor,index_rxn,conc_conserv,files_on, &
                     n_user,rxnmols,tempc,ph,ph_final,tsec,array, &
                     arr_rows,arr_cols)
#if defined(_WIN32)
      USE IFPORT   ! to enable 'SYSTEM' calls
#endif
      USE WEBMOD_IO, ONLY: nowtime
      USE WEBMOD_PHREEQ_MMS, ONLY: xdebug_start, xdebug_stop, nsolute 
      USE WEBMOD_OBSCHEM, ONLY : n_iso
      IMPLICIT NONE
      INCLUDE       'IPhreeqc.f90.inc'
      INTEGER       id              ! 
      INTEGER       count           ! solution count
      INTEGER       solutions(*)    ! solution #'s
      DOUBLE PRECISION        fracs(*)        ! mixing fractions
      INTEGER       index_conserv   ! 
      DOUBLE PRECISION        fill_factor     ! 
      INTEGER       index_rxn       ! 
      DOUBLE PRECISION        conc_conserv(*) ! 
      LOGICAL       files_on        ! write to files
      INTEGER       n_user(*)       !
      DOUBLE PRECISION        rxnmols         !
      DOUBLE PRECISION        tempc           !
      DOUBLE PRECISION        ph              !
      DOUBLE PRECISION        ph_final        !
      DOUBLE PRECISION        tsec            !
      DOUBLE PRECISION        array(*)        !
      INTEGER       arr_rows        ! 
      INTEGER       arr_cols        ! 


      CHARACTER(80) line
      INTEGER       i
      INTEGER       cols
      INTEGER       vtype
      DOUBLE PRECISION        dvalue
      INTEGER       rows
      INTEGER       phr_mix

      INTEGER RunMixF

! Rick's debug/
      integer          startmix(3), endmix(3), nstep, j, et_hyd, et_mix, iresult
      integer, external  ::  getstep, elapsed_time
      logical   fil_temp
      integer, save      ::  et_hold
!      LOGICAL          fil_temp, step1
      LOGICAL          step1, phr_print
!
      data step1/.true./
      save step1
!
      ! Initialize output variables
      do i = 1, nsolute+n_iso*3
        conc_conserv(i) = -987654321.
      enddo
      tempc = -987654321.
      ph = -987654321.
      ph_final = -987654321.
      
           
!      
      fil_temp = files_on
      phr_print = .false.
      startmix(1) = elapsed_time(2)
      nstep=getstep()  !  Serial day number of run
      et_hyd = startmix(1) - et_hold
!
      if(step1.and.files_on) then
!        write(25,"(A)")" nstep soln1 soln2 frac1 frac2 indx_cons index_rxn conc_cons"
        write(26,"(A)")"nstep ent_soln ent_rxn ent_exch ent_surf ent_gas ent_pure_ph "//&
           "ent_sld_soln ent_kin rxn indx_cons	indx_rxn	"//&
           "A	B	C	NoSolns	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10"
        I = SYSTEM('touch .\output\select_mixes')

        step1=.false.
      endif

      if(phr_print.or.files_on) then  ! set phr_print to true in watch window to print phreeq.out for this mix
            iresult = SetOutputFileOn(ID,.true.)
            files_on = .true.  ! stream
      endif
!      if(nstep.ge.xdebug_start.and.nstep.le.xdebug_stop) files_on = .true.
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

      phr_mix = GetSelectedOutputValue(id, 2, 1, vtype, dvalue, line)
      IF (phr_mix.EQ.IPQ_OK) THEN
         IF (vtype.eq.TT_DOUBLE) THEN
            pH_final = dvalue
         ENDIF
!      ELSE
!         RETURN
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
!      if(count.eq.-1) then
!        write(26,1000)nstep, solutions(1), solutions(2), fracs(1), fracs(2), &
!          index_conserv, index_rxn, (conc_conserv(i), i=1,15)
!       endif
!      if(index_rxn.gt.208000000.and.index_rxn.lt.209000000) then
      if(files_on) then
        endmix(1) = elapsed_time(2)
        et_mix = endmix(1) - startmix(1)  
!        write(26,1000)nstep, et_hyd, et_mix, startmix(1), endmix(1), &
!          ' Rxn: ',index_conserv, index_rxn, &
!          ' # of inputs/solns: ',count, (solutions(j),j=1,count)
        write(26,1000)nstep, (n_user(j),j=1,8), &
          ' Rxn: ',index_conserv, index_rxn, &
          ' # of inputs/solns: ',count, (solutions(j),j=1,count)
!        write(25,130)"Date:  ",(nowtime(i),i=1,3)        
!        iresult = SetOutputFileOn(ID,.FALSE.)
        I = SYSTEM('copy .\output\select_mixes + .\phreeqc.0.out .\output\select_mixes')
      endif
      files_on = fil_temp
      iresult = SetOutputFileOn(ID,files_on)
      et_hold = endmix(1)

! 1000 format(I5,' ', 4(I8.8," "),A, 2I10, A, I3, 20I10) 
 1000 format(I5,' ', 8I10 ,A, 2I10, A, I3, 22I10) 

!1000  format(I5,' ', 2(I9.9,' '),2F10.5, 2(I12," "), 15E30.5) 
!      et_hold = endmix(1)
! /Debug            

100   FORMAT(A)
110   FORMAT(TR4,I10,1PG15.7E2)
120   FORMAT(A,I10)
130   FORMAT(A,3I5)
      END FUNCTION phr_mix
