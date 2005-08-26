!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION phr_mix(count,solutions,fracs,index_conserv,
     &                fill_factor,index_rxn,conc_conserv,files_on,
     &                n_user,rxnmols,tempc,ph,tsec,array,
     &                arr_rows,arr_cols)
      IMPLICIT NONE
      INCLUDE       '../IPhreeqc/include/IPhreeqc.f.inc'
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

#if defined(_WIN32)
        INTERFACE 
          FUNCTION RunMixF(count,solutions,fracs,index_conserv,
     &                fill_factor,index_rxn,conc_conserv,files_on,
     &                n_user,rxnmols,tempc,ph,tsec,array,
     &                arr_rows,arr_cols)
            !DEC$ ATTRIBUTES C,REFERENCE::RunMixF
            !DEC$ ATTRIBUTES ALIAS:'_RunMixF'::RunMixF
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
          END FUNCTION RunMixF
        END INTERFACE
#else
        INTEGER RunMixF
#endif


      phr_mix = RunMixF(count,solutions,fracs,index_conserv,
     &                fill_factor,index_rxn,conc_conserv,files_on,
     &                n_user,rxnmols,tempc,ph,tsec,array,
     &                arr_rows,arr_cols)


      rows = GetSelectedOutputRowCount()
      IF (rows.NE.1.AND.rows.NE.3) THEN
        PRINT *, 'ERROR: Expected rows = 1 or 3'
        PRINT *, 'In phr_mix rows = ', rows
        PRINT *, 'No selected_ouput defined?'
        RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(0, 1, vtype, dvalue, line)
      IF (phr_mix.EQ.VR_OK) THEN
         if ('pH' .NE. line(1:2)) THEN
            PRINT *, 'ERROR: Expected pH got:' , line
            PRINT *, 'In phr_mix'
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(1, 1, vtype, dvalue, line)
      IF (phr_mix.EQ.VR_OK) THEN
         IF (vtype.eq.TT_DOUBLE) THEN
            pH = dvalue
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(0, 2, vtype, dvalue, line)
      IF (phr_mix.EQ.VR_OK) THEN
         if ('temp(C)' .NE. line(1:7)) THEN
            PRINT *, 'ERROR: Expected temp(C) got:' , line
            PRINT *, 'In phr_mix'
         ENDIF
      ELSE
         RETURN
      ENDIF

      phr_mix = GetSelectedOutputValue(1, 2, vtype, dvalue, line)
      IF (phr_mix.EQ.VR_OK) THEN
         IF (vtype.eq.TT_DOUBLE) THEN
            tempc = dvalue
         ENDIF
      ELSE
         RETURN
      ENDIF


      cols = GetSelectedOutputColumnCount()
      DO 20 i=3,cols
        phr_mix = GetSelectedOutputValue(1, i, vtype, dvalue, line)
        IF (phr_mix.EQ.VR_OK) THEN
          IF (vtype.eq.TT_DOUBLE) THEN
            conc_conserv(i - 2) = dvalue
          ENDIF
        ELSE
          RETURN
        ENDIF
20    CONTINUE      

100   FORMAT(A)
110   FORMAT(TR4,I10,1PG15.7E2)
120   FORMAT(A,I10)
      END FUNCTION phr_mix
