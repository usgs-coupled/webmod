      !PHREEQC_MMS prototypes
      
      INTEGER, EXTERNAL :: CreateIPhreeqcMMS
      INTEGER, EXTERNAL :: DestroyIPhreeqcMMS
      
      INTEGER, EXTERNAL :: phr_multicopy !!(keyword, srcarray, targetarray, count)
      INTEGER, EXTERNAL :: phr_precip    !!(soln_id, count, aspecies, aconc, tempc, ph)
      INTEGER, EXTERNAL :: phr_mix       !!(count, solutions, fracs, dest)

      INTEGER, EXTERNAL :: build_tally_table
      INTEGER, EXTERNAL :: get_tally_table_rows_columns   !!(rows, cols)
      INTEGER, EXTERNAL :: get_tally_table_column_heading !!(col,etype,string)	  
      INTEGER, EXTERNAL :: get_tally_table_row_heading    !!(row,string)

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
