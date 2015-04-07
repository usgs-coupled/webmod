       INTERFACE
        FUNCTION CreateIPhreeqcMMS()
         IMPLICIT NONE
         INTEGER(KIND=4) :: CreateIPhreeqcMMS
        END FUNCTION CreateIPhreeqcMMS
       END INTERFACE


       INTERFACE
        FUNCTION DestroyIPhreeqcMMS(ID)
         IMPLICIT NONE
         INTEGER(KIND=4), INTENT(IN) :: ID
         INTEGER(KIND=4)             :: DestroyIPhreeqcMMS
        END FUNCTION DestroyIPhreeqcMMS
       END INTERFACE
       
       
       
       INTERFACE
        FUNCTION phr_mix(id,count,solutions,fracs,index_conserv, &
                     fill_factor,index_rxn,conc_conserv,files_on, &
                     n_user,rxnmols,tempc,ph,ph_final,tsec,array, &
                     arr_rows,arr_cols)
         IMPLICIT NONE
         INTEGER(KIND=4), INTENT(IN)    :: id              ! 
         INTEGER(KIND=4), INTENT(IN)    :: count           ! solution count
         INTEGER(KIND=4), INTENT(IN)    :: solutions(*)    ! solution #'s
         REAL(KIND=8),    INTENT(IN)    :: fracs(*)        ! mixing fractions
         INTEGER(KIND=4), INTENT(IN)    :: index_conserv   ! 
         REAL(KIND=8),    INTENT(IN)    :: fill_factor     ! 
         INTEGER(KIND=4), INTENT(IN)    :: index_rxn       ! 
         REAL(KIND=8),    INTENT(OUT)   :: conc_conserv(*) ! 
         LOGICAL(KIND=4), INTENT(IN)    :: files_on        ! write to files
         INTEGER(KIND=4), INTENT(IN)    :: n_user(*)       !
         REAL(KIND=8),    INTENT(IN)    :: rxnmols         !
         REAL(KIND=8),    INTENT(INOUT) :: tempc           !
         REAL(KIND=8),    INTENT(INOUT) :: ph              !
         REAL(KIND=8),    INTENT(INOUT) :: ph_final              !
         REAL(KIND=8),    INTENT(IN)    :: tsec            !
         REAL(KIND=8),    INTENT(INOUT) :: array(*)        !
         INTEGER(KIND=4), INTENT(IN)    :: arr_rows        ! 
         INTEGER(KIND=4), INTENT(IN)    :: arr_cols        ! 
         INTEGER(KIND=4)                :: phr_mix
        END FUNCTION phr_mix
       END INTERFACE

       INTERFACE
        FUNCTION phr_multicopy(id, keyword, srcarray, targetarray, count)
         IMPLICIT NONE
         INTEGER(KIND=4), intent(IN) ::  id
         CHARACTER(*),    intent(IN) ::  keyword
         INTEGER(KIND=4), intent(IN) ::  srcarray(*)
         INTEGER(KIND=4), intent(IN) ::  targetarray(*)
         INTEGER(KIND=4), intent(IN) ::  count
         INTEGER(KIND=4)             ::  phr_multicopy
        END FUNCTION phr_multicopy
       END INTERFACE

       INTERFACE
        FUNCTION phr_precip(id,soln_id,count,aspecies,aconc,tempc,ph)
         IMPLICIT NONE
         INTEGER(KIND=4), intent(IN) ::       id
         INTEGER(KIND=4), intent(IN) ::       soln_id
         INTEGER(KIND=4), intent(IN) ::       count
         CHARACTER(*),    intent(IN) ::       aspecies(*)
         DOUBLE PRECISION,intent(IN) ::       aconc(*)
         DOUBLE PRECISION,intent(IN) ::       tempc
         DOUBLE PRECISION,intent(IN) ::       ph
         INTEGER(KIND=4)             ::       phr_precip
        END FUNCTION phr_precip
       END INTERFACE

       INTERFACE
        FUNCTION build_tally_table(id)
         IMPLICIT NONE
         INTEGER(KIND=4), INTENT(IN) :: id
         INTEGER(KIND=4)             :: build_tally_table
        END FUNCTION build_tally_table
       END INTERFACE


       INTERFACE
        FUNCTION get_tally_table_rows_columns(id,rows,columns)
         IMPLICIT NONE
         INTEGER(KIND=4), INTENT(IN)  :: id
         INTEGER(KIND=4), INTENT(OUT) :: rows
         INTEGER(KIND=4), INTENT(OUT) :: columns
         INTEGER(KIND=4)              :: get_tally_table_rows_columns
        END FUNCTION get_tally_table_rows_columns
       END INTERFACE

       INTERFACE
        FUNCTION get_tally_table_column_heading(id,column,type,string)
         IMPLICIT NONE
         INTEGER(KIND=4),  INTENT(IN)  :: id
         INTEGER(KIND=4),  INTENT(IN)  :: column
         INTEGER(KIND=4),  INTENT(OUT) :: type
         CHARACTER(LEN=*), INTENT(OUT) :: string
         INTEGER(KIND=4)               :: get_tally_table_column_heading
        END FUNCTION get_tally_table_column_heading
       END INTERFACE


       INTERFACE
        FUNCTION get_tally_table_row_heading(id,row,string)
         IMPLICIT NONE
         INTEGER(KIND=4),  INTENT(IN)  :: id
         INTEGER(KIND=4),  INTENT(IN)  :: row
         CHARACTER(LEN=*), INTENT(OUT) :: string
         INTEGER(KIND=4)               :: get_tally_table_row_heading
        END FUNCTION get_tally_table_row_heading
       END INTERFACE
      ! get_tally_table_column_heading
      ! Entity types

      INTEGER    ET_SOLUTION
      PARAMETER (ET_SOLUTION=1)
      
      INTEGER    ET_REACTION
      PARAMETER (ET_REACTION=2)

      INTEGER    ET_EXCHANGE
      PARAMETER (ET_EXCHANGE=3)

      INTEGER    ET_SURFACE
      PARAMETER (ET_SURFACE=4)

      INTEGER    ET_GAS_PHASE
      PARAMETER (ET_GAS_PHASE=5)
      
      INTEGER    ET_PURE_PHASE
      PARAMETER (ET_PURE_PHASE=6)

      INTEGER    ET_SS_PHASE
      PARAMETER (ET_SS_PHASE=7)

      INTEGER    ET_KINETICS
      PARAMETER (ET_KINETICS=8)

      INTEGER    ET_MIX
      PARAMETER (ET_MIX=9)

      INTEGER    ET_TEMPERATURE
      PARAMETER (ET_TEMPERATURE=10)

      INTEGER    ET_UNKNOWN
      PARAMETER (ET_UNKNOWN=11)
