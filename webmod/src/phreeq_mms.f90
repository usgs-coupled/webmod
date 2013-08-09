! *********************************************************************
!     phreeq_mms.f: compute solute concentrations and mass fluxes for 
!     hillslope reservoirs, ponds, and streams given fluxes calculated
!     in other modules
!
!     Rick Webb - August 2003
!
!     18 Mar 09 - Converted to use Fortran90 Module (WEBMOD_PHREEQ_MMS)
!
!     13 Mar 05 - Removed references to vmix_hold_ variables
!                 vmix_hold_well, and vmix_hold_diversion
!                 since these are transfered on the same time step
!                 as vmix_stream and vmix_diversion.
!                 Need to review mixing so that there is no double
!                 accounting (i.e. remove initial irrigation volumes
!                 at concentrations of t0 before mixing upstream
!                 t1 concentrations to make t1 final (us, sat, stream)
!
! **********************************************************************
! ***************** START MODULE *************
      MODULE WEBMOD_PHREEQ_MMS
      IMPLICIT NONE
      INCLUDE 'fmodules.f90.inc'
      INCLUDE 'IPhreeqc.f90.inc'
      INCLUDE 'mms_phreeqc.f90.inc'
      
!  Constants
!   M2mM is equal to 1000 to convert millimoles to moles. Duplicate into
!   a_thousand for other conversions
!
      real, PARAMETER :: almost_one = 0.9999 ! to avoid divide by zero in fractionation calculations
      real, PARAMETER :: zero = 0.0
      integer, PARAMETER :: one = 1  !  use to index evap water
      real, PARAMETER :: M2mM = 1000.0
      real, PARAMETER :: a_thousand = 1000.0
      real, PARAMETER :: a_million = 1000000.0
      real, PARAMETER :: inch2m = 0.0254
      real, PARAMETER :: ln_10 = 2.30258509299405 ! naperian log of 10 for converting log alpha to ln alpha
      integer, PARAMETER :: nentity = 11
      integer, PARAMETER :: nconvert = 3
      integer, PARAMETER :: nchemvar = 10
! maximum number deposition sources (precip + internal and external irrig)
      integer, PARAMETER :: maxdep = 3
! maximum number of species listed in phreeq_lut file or in tally_table
      integer, PARAMETER :: maxspecies = 1000
! maximum number of columns in tally table
      integer, PARAMETER :: maxtallycol = 100
      integer, save :: no_rxn(nentity)
      DATA no_rxn/nentity*-1/
      logical, save :: step1
      DATA step1/.true./
! Three coefficients for vapor liquid and vapor ice fractionation of isotopes.

!      INTERFACE
!         FUNCTION PHREEQ_ID(STRING, sol_lut)
!           CHARACTER(LEN=*) :: STRING 
!           CHARACTER(LEN=*) :: sol_lut(:,:)
!           INTEGER :: PHREEQ_ID
!         END FUNCTION PHREEQ_ID
!       END INTERFACE

!      INTERFACE
!         FUNCTION LENGTH(STRING)
!           CHARACTER(LEN=*) :: STRING
!           INTEGER :: LENGTH
!         END FUNCTION LENGTH
!      END INTERFACE

      INTERFACE
         FUNCTION solnnum(time,res_id,chemdat,mru,nac,hydro,stat)
           INTEGER :: time, res_id,chemdat,mru,nac,hydro,stat
           INTEGER :: solnnum
         END FUNCTION solnnum
       END INTERFACE

      INTERFACE
         FUNCTION isoln(soln_num,nchemdat,nmru,nac,nhydro,&
           ires,ichemdat,imru,inac,ihydro)
           INTEGER :: soln_num, nchemdat, nmru, nac, nhydro
           INTEGER :: ires, ichemdat, imru, inac, ihydro
           INTEGER :: isoln
         END FUNCTION isoln
       END INTERFACE

      INTERFACE 
         FUNCTION fill_ent(n_user,soln_number,nchemdat,nmru,&
           nac,clark_segs,src_init)
           INTEGER :: n_user(:),soln_number, nchemdat, nmru
           INTEGER :: nac, clark_segs
           INTEGER, INTENT(IN) :: src_init(:,:)
           INTEGER :: fill_ent
         END FUNCTION fill_ent
       END INTERFACE

      INTERFACE 
         FUNCTION wetbulb(pres,ta,rh)
           REAL, INTENT(IN) :: pres, ta, rh
           REAL :: wetbulb
         END FUNCTION wetbulb
       END INTERFACE

      INTERFACE
         FUNCTION checkfracs(count,solns,fracs,mixture)
           INTEGER :: count,solns(:),mixture
           DOUBLE PRECISION :: fracs(:)
           INTEGER :: checkfracs
         END FUNCTION checkfracs
       END INTERFACE

      INTERFACE
        INTEGER FUNCTION chemflag(index,chvar_lut, soln_number,nchemvar)
          INTEGER :: index, chvar_lut(:,:)
          INTEGER :: soln_number, nchemvar
        END FUNCTION chemflag
      END INTERFACE
      
      interface
         function update_chem(indx,totvol,imetric,&
              conc,ph,tempc,tally_table,n_ent,indxmru,&
              indxbas,restype)
           integer :: indx, imetric,n_ent(:)
           integer :: indxmru,indxbas,restype
           double precision :: totvol,conc(:),ph,tempc
           double precision :: tally_table(:,:)
           integer :: update_chem
         end function update_chem
      end interface
      
      INTERFACE
         FUNCTION fractionate(phase,indx,evap,totvol,ison,rh,tempevap)
           INTEGER :: phase, indx
           DOUBLE PRECISION :: evap,totvol,ison,rh,tempevap
           INTEGER :: fractionate
         END FUNCTION fractionate
      END INTERFACE

      INTERFACE
         FUNCTION reset_DI(ID,tempevap)
           INTEGER :: ID
           REAL :: tempevap
           INTEGER :: reset_DI
         END FUNCTION reset_DI
      END INTERFACE
      
!      INTERFACE
!         FUNCTION TEST(conc1)
!         double precision, intent (inout) :: conc1(:) 
!         INTEGER :: TEST
!         END FUNCTION TEST
!      END INTERFACE

       
!  Functions
!      INTEGER, EXTERNAL :: phreeq_id,chemflag, isoln, chemflag
!      INTEGER, EXTERNAL ::  update_chem, chem2var,fill_ent,checkfracs
!      INTEGER, EXTERNAL ::  chem2var, phreeq_id, length
      INTEGER, EXTERNAL ::  phreeq_id !!SRC, RunString

! strings to construct custom header in pqi file      
      integer, save :: isoh1_len, isoh2_len, isogl_len, isogs_len, sol_h2_len
      character*60, save :: isogs, isogl
      character*256 :: aline

      character*3000, save :: sol_header1, sol_header2, iso_header1,iso_header2
      character*3000, save :: inp_dir, phreeqmms_pqi
      character*3000, save :: cheminit_file, readline, phreeq_database



! Basin geometry - kPa is Air pressure in kiloPascals needed to derive wetbulb temperature for precip
      real, save :: basin_area
      real, save, allocatable :: ar_fill(:,:), mru_area(:),qdffrac(:)
      real, save, allocatable :: mru_elev(:), kPa(:)
      real, save, allocatable :: covden_sum(:), covden_win(:)
      real, save, allocatable :: snowcov_area(:),ac(:,:)
      real :: covden, acf

!  Dimensions, indices, and counters
      integer ::  i,j,ia,ires,imru,inac,ihydro,irip
      integer :: ichemdat,iunit,imet,unit_type
      integer, save :: nchem_ext, chem_ext,ppt_chem, nchem_sets
      integer, save ::  nmru, nac, nsolute, nchemobs,nchemdat
      integer, save ::  nchemdat_obs, nchemdep, nsoln, nchan
      integer, save ::  nresinp, nmru_res,nac_nmru_nresinp
      integer, save ::  chem_sim
      integer, save ::  clark_segs, nxkbin, nhydro
      integer, save, allocatable :: nacsc(:)

!  Indices describing chemistry of irrigation and regional ground water inputs
      integer, save, allocatable :: irrig_int_src(:), src_ext_irrig(:)
      integer, save, allocatable ::  src_gw1(:), src_gw2(:)
      
!  Initial solutions and reactants. Overwritten with solutions and reactants
!  to be tracked for each reservoir.

      integer, save, allocatable :: src_init(:,:)

! Booleans

      logical, save ::  chemdat_exists, phr_tf

! The c_ parameters determine what reservoirs and metrics are
! loaded into the ch_var variables

      integer, save ::  nphrsolns
! chvar_lut is an integer array describing the reservoirs and metrics 
! flagged in the nchemvar pars. The reservoir is described both by the
! unique solution number and the row number of the c_indx/c_chem tables.
! Note that reaction and net metrics only apply when the desired unit
! is mass and would not make sense for concentration 
!
!  Fields in chvar_lut: 
!  1) row number of solution, 
!  2) solution number, 
!  3) desired metric (c_metric: initial, in , out, reaction, final, net)
!  4) desired unit of output (c_units)
!  5) desired reservoir (c_ires)
!  6) desired observed water quality site (c_obs_indx)
!  7) desired model response unit (c_mru)
!  8) desired topographic index bin (c_stindx)
!  9) desired drainage segment (c_hyd_indx)
! 10) riparian index (for hillslope, riparian, or upland)

      integer, save, allocatable :: chvar_lut(:,:)
!
! chvar_conv contains values needed to convert mass and volumes in
! the c_chem matrix into the desired output variables described by
! the nchemvar parameters. The numerator is used to convert
! the moles in the c_chem matrix to the desired mass units. The
! mass conversions (to molecular weight and ion charge) for
! each solute are place in each chvar_conv(ivar, isol) where ivar
! is the index (chvar_01, _02, etc) and isol is the solute number.
! The denominator is equal to the area (canopy, snowpack, etc)  if
! a load is to be calculated or equal to one if mass or concentration
! is to be calculated. The area is stored in the last field
! (nsolute+1). At the end of each time step the volumes in each reservoir
! are used with these static conversion factors to produce the final
! values of mass, loads, or concentrations as indicated in table dimensioned
! by nchemvar.
! 
      double precision, save, allocatable :: chvar_conv(:,:)
      
!
! c_indx lists, in ascending order, all reservoir solution numbers
! representing solution compositions at the beginning (t0) and
! end (t+1) of each time step. The row number (or case) of the table
! is the lookup key common to c_indx and c_chem. The last rows of the
! c_indx and the c_chem tables describe pseudo-solutions compiled for
! indivual MRUs and the entire basin from volumes and masses computed
! for individual reservoirs.
! The second column will equal one ('1') if the reservoir has been
! selected for output to one of the 10 ch_vars.
!
! Since a single row in c_chem contains initial values, inputs,
! outputs, ET and reactions, and final values, there is no need to
! maintain this information in two solutions (t0 and t+1), so the summary
! fluxes will always be saved in the row associated with the initial t solution.
!
!
      integer, save, allocatable ::  c_indx(:,:) 
      integer, save :: indx, imetric

! The main data storage matrix, c_chem, is the master storage table of type geochem.
!
! TYPE (geochem) c_chem. Records moles and volumes init, in, out, rxn(ET), and final for
! all permanent and temporary reservoirs in the model. Also records final isotopic values
! (where applicable) along with temperature and pH.
!
      TYPE :: geochem
         double precision, allocatable :: M(:,:)   ! Moles of each solute for init, in, out, rxn, and final
         double precision :: vol(5)   ! volumes for init, in, out, ET, and final
         double precision, allocatable :: delta(:,:) ! delta values for any of nsolutes that are isotopes
         double precision :: Temp(5)   ! Temperature of reservoir inputs and outputs
         double precision :: pH(5)   ! pH of reservoir inputs and outputs
!         double precision :: Ent_Type(:)   ! pH of reservoir inputs and outputs
!         double precision :: Ent_M(:)   ! pH of reservoir inputs and outputs
      END TYPE geochem
     
      TYPE(geochem), save, allocatable :: c_chem(:) ! to be allocated by the total number of solutions, nphrsolns
      TYPE(geochem), save, allocatable :: c_chem_uzgen(:) ! additional rows of c_chem to hold uz composite mass and fluxes
      TYPE(geochem), save, allocatable :: c_chem_uzrip(:) ! additional rows of c_chem to hold uz composite mass and fluxes
      TYPE(geochem), save, allocatable :: c_chem_uzup(:) ! additional rows of c_chem to hold uz composite mass and fluxes
! These parameters enable better understanding of fluxes recorded in c_chem%M and c_chem_vol; index 4 is rxn for Moles or et for volumes.
      integer, PARAMETER :: init = 1
      integer, PARAMETER :: in = 2
      integer, PARAMETER :: out = 3
      integer, PARAMETER :: rxn = 4 ! for Moles
      integer, PARAMETER :: ET = 4  ! for volumes of water
      integer, PARAMETER :: fin = 5
      double precision :: vol  ! scalar
      logical :: not_isotope

! After each mix, the tally table is populated with the 
! net reactions (solution, kinetics, exchange, etc - tallycols)
! for all species (tally rows) that it knows about in a given run.
! The solutes of interest, nsolutes, are usually only a subset of all
! the species that phreeqc tracks. The row number index of the
! tally table that corresponds to each of the 'i' nsolutes will
! be stored in sol_id(i)%tally. mult contains a multiplier for
! each tally column to indicate whether the reaction (or other entity)
! added ions the the solution, mult()=1.0, or removed ions, mult()= -1.0

      double precision, allocatable ::   tally_table(:,:), mult(:)
      integer, save :: ntally_rows, ntally_cols
      integer, save, allocatable ::  n_user(:), n_ent(:)
      character*50, save, allocatable ::  tally_col_label(:)
      integer, save, allocatable ::  tally_col_type(:)
      character*12, save, allocatable ::  tally_row_label(:)

!  mixing variables

      integer, save :: validnacs, indx_rxn, mixture
      integer, save, allocatable :: src(:),srcdep(:),dest(:)
      double precision,save :: rxnmols,tempc,ph,ph_final,tsec,fill_factor,totvol,totvol_can
      double precision :: basin_in_vol, basin_out_vol, evap_frac
      double precision :: delta_res0, delta_res, delta_res_permil
      double precision :: log_a, eps_q, eps_diff
      double precision,allocatable:: mru_in_vol(:),mru_out_vol(:)
      double precision, save, allocatable :: fracs(:), conc(:)
      double precision, allocatable :: fracsdep(:)
      double precision :: evap, ison, rh  ! for raleigh fractionation

! convfactor times Moles per liter converts to user units

      double precision, save, allocatable :: convfactor(:)

! ch_vars provide a single-dimensioned variable for plotting or
! statistical reporting not possible with the multidimensioned
! model space (nmru, nac, nsolute, nmru_res, in, out, etc), The
! content of each ch_var is defined by the variables c_ires,
! c_metric,c_mru, c_stindx, c_hyd_indx, c_units, and c_rip
! all of which are dimensioned by nchemvar.

      double precision, save, allocatable :: ch_var_01_sol(:), ch_var_02_sol(:)
      double precision, save, allocatable :: ch_var_03_sol(:), ch_var_04_sol(:)
      double precision, save, allocatable :: ch_var_05_sol(:), ch_var_06_sol(:)
      double precision, save, allocatable :: ch_var_07_sol(:), ch_var_08_sol(:)
      double precision, save, allocatable :: ch_var_09_sol(:), ch_var_10_sol(:)

      double precision, save :: ch_var_01_m3, ch_var_02_m3
      double precision, save :: ch_var_03_m3, ch_var_04_m3
      double precision, save :: ch_var_05_m3, ch_var_06_m3
      double precision, save :: ch_var_07_m3, ch_var_08_m3
      double precision, save :: ch_var_09_m3, ch_var_10_m3
      
      double precision, save :: ch_var_01_tempc, ch_var_02_tempc
      double precision, save :: ch_var_03_tempc, ch_var_04_tempc
      double precision, save :: ch_var_05_tempc, ch_var_06_tempc
      double precision, save :: ch_var_07_tempc, ch_var_08_tempc
      double precision, save :: ch_var_09_tempc, ch_var_10_tempc

      double precision, save :: ch_var_01_pH, ch_var_02_pH
      double precision, save :: ch_var_03_pH, ch_var_04_pH
      double precision, save :: ch_var_05_pH, ch_var_06_pH
      double precision, save :: ch_var_07_pH, ch_var_08_pH
      double precision, save :: ch_var_09_pH, ch_var_10_pH


      integer, save, allocatable :: c_ires(:),c_metric(:)
      integer, save, allocatable :: c_mru(:), c_stindx(:), c_rip(:)
      integer, save, allocatable :: c_hyd_indx(:), c_units(:)
      integer, save, allocatable :: c_obs_indx(:)
      integer, save :: chemdat_flag ! indicates if file of time series of chem obs is available.

!
! chmru_soln indicates the index, or row number for each pseudo-solution
! that summarizes fluxes into and out of each MRU. Ditto for chbas_soln,
! the basin pseudo reservoir. These are reflected as indxm and indxm
! at any function call. indxuz will reference the riparian pseudosolution, ch_rip_solution,
! or the upland pseudosolution, ch_upland_soln as distinguised by riparian_thresh.
!

      integer, save, allocatable :: chmru_soln(:), chuz_soln(:), ch_rip_soln(:), ch_upland_soln(:)
      integer, save :: chbas_soln, indxm, indxb, indxuz

      double precision, save, allocatable ::&
           ch_basin_mass_g(:),ch_basin_conc_mgL(:), &
           ch_basin_rxn_g(:),&
           ch_basin_in_mgL(:), ch_basin_out_mgL(:),&
           ch_basin_in_permil(:), ch_basin_out_permil(:),&
           ch_basin_in_g(:), ch_basin_out_g(:),&
           ch_basin_net_g(:), ch_basin_in_gm2(:),&
           ch_basin_out_gm2(:), ch_basin_net_gm2(:)

      double precision, save :: ch_basin_out_tempC, ch_basin_out_pH
      double precision, save :: ch_basin_vol_m3

      double precision, save, allocatable ::&
           ch_mru_vol_m3(:),&
           ch_mru_mass_g(:,:),&
           ch_mru_conc_mgL(:,:),&
           ch_mru_rxn_g(:,:),&
           ch_mru_in_mgL(:,:),&
           ch_mru_out_mgL(:,:),&
           ch_mru_in_permil(:,:),&
           ch_mru_out_permil(:,:),&
           ch_mru_in_g(:,:),&
           ch_mru_out_g(:,:),&
           ch_mru_net_g(:,:),&
           ch_mru_in_gm2(:,:),&
           ch_mru_out_gm2(:,:),&
           ch_mru_net_gm2(:,:),&
           ch_mru_out_tempC(:),ch_mru_out_pH(:)

      integer, allocatable :: solns(:)
!
! Gordon-Craig evaporative fractionation model (relative humidity from obs and maximum diffusion in phreeq_lut)
!

      integer, PARAMETER :: snow = 0  ! flags sublimation from snow so use vapor/solid curve
      integer, PARAMETER :: water = 1 ! evaporation from canopy and UZ as water so use vapor/liquid curve
      integer, PARAMETER :: transp = 3 ! transpired water undergoes no fractionation
      double precision, save, allocatable :: iso_n(:)    ! iso_n
      double precision,allocatable:: evap_iso(:)    ! delta values for evap distinct from remaining water
      double precision :: eps_eq, delta_evap_permil
!
! Simple ionic pulse model with incongruous melting of snowpack
!
      integer, save, allocatable :: snow_ion_pulse(:)  ! flag to alert user of ionic pulse
      real :: percent_melt, max_factor, ion_factor   ! variables to compute pulse concentrations
! snow_ion_factor, and snow_iso_depl are parameters to simulate the concentration of solutes and
! depletion of heavy isotopes in snowmelt. snowpack_D and snowpack_18O are temporary variables describing
! the snowpack after sublimation. That value is used to assign delta_D and delta_18O to the DI used to
! create lighter isotopes in the melt.
      real, save, allocatable :: snow_ion_factor(:), snowmelt_D_depl(:), snowmelt_18O_depl(:)
      double precision :: evap_D_permil, evap_18O_permil ! permil deltas in evaporated vapor
      double precision :: res_D_permil, res_18O_permil ! permil deltas in reservoir after evaporation.
      double precision :: delta_D, delta_18O ! snowmelt_D_depl, snowmelt_18O_depl, and snow_ion_factor
                                             ! used to assign permil deltas for D and 18O after sublimation.
!
!  Initialization
!  These tables assign one or more typical 'sets' for the reservoirs in a MRU or 
!  hydrologic feature. Each MRU or hydro feature can be assigned to one of these sets
!  using the init_ parameters.
!
      integer, save :: init_soln_ppt
      integer, save, allocatable ::&
             eq_phset_table(:,:),&
             exchset_table(:,:),&
             init_eq_ph_hydro(:),&
             init_eq_ph_mru(:),&
             init_exch_hydro(:),&
             init_exch_mru(:),&
             init_kin_hydro(:),&
             init_kin_mru(:),&
             init_rxn_hydro(:),&
             init_rxn_mru(:),&
             init_soln_ext(:),&
             init_soln_hydro(:),&
             init_soln_mru(:),&
             init_surf_hydro(:),&
             init_surf_mru(:),&
             kinset_table(:,:),&
             rxnset_table(:,:),&
             solnset_table(:,:),&
             surfset_table(:,:),&
             uzindxinit(:,:,:),&
             uz_spread(:)
!
! debug and misc parameters
!
!
      double precision dvalue
      integer vtype, leng, id_len, ID
      character(15) heading
      character(60) line
      integer  nstep, datetime(6),xdebug_start,xdebug_stop
      integer iphrq_mru

      END MODULE WEBMOD_PHREEQ_MMS
!
!*********************************************************
!     Main phreeq_mms routine
!*********************************************************
!
      integer function phreeq_mms(arg)
      USE WEBMOD_IO, only: chemout_file_unit
      implicit none
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
      CHARACTER*256 SVN_ID
! Functions
      INTEGER, EXTERNAL :: phreeqmms_decl,phreeqmms_init,&
                           phreeqmms_run, phreeqmms_clean
      
      SVN_ID = &
           '$Id: phreeq_mms.f 36 2007-06-08 17:38:24Z rmwebb $ '
     
      phreeq_mms = 0

      IF ( Arg.EQ.'run' ) THEN
        phreeq_mms = phreeqmms_run()
      ELSEIF ( Arg.EQ.'declare' ) THEN
        phreeq_mms = phreeqmms_decl()
      ELSEIF ( Arg.EQ.'initialize' ) THEN
        phreeq_mms = phreeqmms_init()
      ELSEIF ( Arg.EQ.'clean' ) THEN
        phreeq_mms = phreeqmms_clean()
      ENDIF

      if(phreeq_mms.eq.1) then
         close (chemout_file_unit)
         close (14)
         close (16)
         close (17)
      end if

      END FUNCTION phreeq_mms

!**************************************************************
!
!     phreeqdecl - makes public variable declarations for the
!                     phreeq_mms module
!
!**************************************************************
 
      integer function phreeqmms_decl()

      USE WEBMOD_PHREEQ_MMS
      USE WEBMOD_OBSCHEM, ONLY : n_iso

#if defined(_WIN32)
      USE IFPORT
#endif
!      include 'fmodules.inc'
!      INCLUDE 'IPhreeqc.f90.inc'
!      INCLUDE 'mms_phreeqc.inc'

      phreeqmms_decl = 1
! Get dimensions

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nac = getdim('nac')
        if ( nac.eq.-1 ) return
      nxkbin = getdim('nxkbin')
        if ( nxkbin.eq.-1 ) return
      nhydro = getdim('nhydro')
        if ( nhydro.eq.-1 ) return
      nchan = getdim('nchan')
        if ( nchan.eq.-1 ) return
      nchemobs = getdim('nchemobs')
        if ( nchemobs.eq.-1 ) return
      nchem_ext = getdim('nchem_ext')
        if ( nchem_ext.eq.-1 ) return
      nsolute = getdim('nsolute')
        if ( nsolute.eq.-1 ) return
      nresinp = getdim('nresinp')
        if(nresinp.eq.-1) return
      nac_nmru_nresinp = getdim('nac_nmru_nresinp')
        if(nac_nmru_nresinp.eq.-1) return
      nmru_res = getdim('nmru_res')
        if ( nmru_res.eq.-1 ) return
      nchem_sets = getdim('nchem_sets')
        if ( nchem_sets.eq.-1 ) return
! nchemdat is equal to the number of deposition sources
! (nchemdep = one precip + nchem_ext), plus the number of unique
! water quality sampling sites, nchemobs.
!
      nchemdep = 1 + nchem_ext
!
! number of point concentrations
!
      nchemdat = nchemdep + nchemobs
!
! nsoln is the space to store the initial and final solute mass for each 
! hydrologic reservoir in the model
!
      nsoln = 1+nchemdat+2*nmru*(nac+14)+2*(nhydro+1)+1+nmru
!         old: 2*(((nmru+nac)*nmru)+nhydro)
!
! The ch_var_xx_sol and ch_var_xx_m3 variables below allow the user to track the
! solutes and volumes of any specific model reservoir. The reservoir
! is defined with the parameters associated with the nchemvar dimension
!
      allocate(ch_var_01_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_01_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(1)','c_units', ch_var_01_sol).ne.0) return

      allocate(ch_var_02_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_02_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(2)','c_units', ch_var_02_sol).ne.0) return

      allocate(ch_var_03_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_03_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(3)','c_units', ch_var_03_sol).ne.0) return

      allocate(ch_var_04_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_04_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(4)','c_units', ch_var_04_sol).ne.0) return

      allocate(ch_var_05_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_05_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(5)','c_units', ch_var_05_sol).ne.0) return

      allocate(ch_var_06_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_06_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(6)','c_units', ch_var_06_sol).ne.0) return

      allocate(ch_var_07_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_07_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(7)','c_units', ch_var_07_sol).ne.0) return

      allocate(ch_var_08_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_08_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(8)','c_units', ch_var_08_sol).ne.0) return

      allocate(ch_var_09_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_09_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(9)','c_units', ch_var_09_sol).ne.0) return

      allocate(ch_var_10_sol(nsolute))
      if(declvar('phreeqmms', 'ch_var_10_sol', 'nsolute', nsolute,&
           'double',&
           'Mass, flux or concentration of a solute for a modeled '//&
           'reservoir as defined by parameters with dimension '//&
           'nchemvar(10)','c_units', ch_var_10_sol).ne.0) return

      if(declvar('phreeqmms', 'ch_var_01_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(1)','cubic meters',&
           ch_var_01_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_02_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(2)','cubic meters',&
           ch_var_02_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_03_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(3)','cubic meters',&
           ch_var_03_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_04_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(4)','cubic meters',&
           ch_var_04_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_05_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(5)','cubic meters',&
           ch_var_05_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_06_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(6)','cubic meters',&
           ch_var_06_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_07_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(7)','cubic meters',&
           ch_var_07_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_08_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(8)','cubic meters',&
           ch_var_08_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_09_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(9)','cubic meters',&
           ch_var_09_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_10_m3', 'one', 1,&
           'double',&
           'The volume the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(10)','cubic meters',&
           ch_var_10_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_var_01_tempc', 'one', 1, 'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(1)','degrees Celsius',&
         ch_var_01_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_02_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(2)','degrees Celsius',&
         ch_var_02_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_03_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(3)','degrees Celsius',&
         ch_var_03_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_04_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(4)','degrees Celsius',&
         ch_var_04_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_05_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(5)','degrees Celsius',&
         ch_var_05_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_06_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(6)','degrees Celsius',&
         ch_var_06_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_07_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(7)','degrees Celsius',&
         ch_var_07_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_08_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(8)','degrees Celsius',&
         ch_var_08_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_09_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(9)','degrees Celsius',&
         ch_var_09_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_10_tempc', 'one', 1,'double',&
         'The temperature of the reservoir defined by the parameters '//&
         'associated with dimension nchemvar(10)','degrees Celsius',&
         ch_var_10_tempc).ne.0) return

      if(declvar('phreeqmms', 'ch_var_01_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(1)','pH units',&
           ch_var_01_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_02_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(2)','pH units',&
           ch_var_02_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_03_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(3)','pH units',&
           ch_var_03_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_04_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(4)','pH units',&
           ch_var_04_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_05_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(5)','pH units',&
           ch_var_05_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_06_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(6)','pH units',&
           ch_var_06_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_07_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(7)','pH units',&
           ch_var_07_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_08_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(8)','pH units',&
           ch_var_08_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_09_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(9)','pH units',&
           ch_var_09_pH).ne.0) return

      if(declvar('phreeqmms', 'ch_var_10_pH', 'one', 1,&
           'double',&
           'The pH of the reservoir defined by the parameters '//&
           'associated with dimension nchemvar(10)','pH units',&
           ch_var_10_pH).ne.0) return

!
! The variables that follow track the total storage and fluxes of solutes
! into and out of the basins and mrus. To track fluxes in specific 
! reservoirs described above, use the ch_var variables above.
!
      if(declvar('phreeqmms', 'ch_basin_vol_m3', 'one', 1,&
           'double','Total volume of water in the basin',&
           'm3', ch_basin_vol_m3).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_out_tempC', 'one', 1,&
           'double','Temperature of water leaving the basin',&
           'deg C', ch_basin_out_tempC).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_out_pH', 'one', 1,&
           'double','pH of water leaving the basin',&
           'pH units', ch_basin_out_pH).ne.0) return

      allocate(ch_basin_mass_g(nsolute))
      if(declvar('phreeqmms', 'ch_basin_mass_g', 'nsolute', nsolute,&
           'double','Total mass of solute in basin',&
           'g', ch_basin_mass_g).ne.0) return

      allocate(ch_basin_rxn_g(nsolute))
      if(declvar('phreeqmms', 'ch_basin_rxn_g', 'nsolute', nsolute,&
           'double','Total mass of solute created or consumed by '//&
           'reactions in basin',&
           'g', ch_basin_rxn_g).ne.0) return

      allocate(ch_basin_conc_mgL(nsolute))
      if(declvar('phreeqmms', 'ch_basin_conc_mgL', 'nsolute', nsolute,&
           'double','Average concentration of solute in basin',&
           'mg/L', ch_basin_conc_mgL).ne.0) return

      allocate(ch_basin_in_mgL(nsolute))
      if(declvar('phreeqmms', 'ch_basin_in_mgL', 'nsolute', nsolute,&
           'double','Average solute concentration of precip ,'//&
           'ground water influx and irrigation', 'mg/L',&
           ch_basin_in_mgL).ne.0) return

      allocate(ch_basin_out_mgL(nsolute))
      if(declvar('phreeqmms', 'ch_basin_out_mgL', 'nsolute', nsolute,&
           'double','Average solute concentration of stream at outlet',&
           'mg/L', ch_basin_out_mgL).ne.0) return

      allocate(ch_basin_in_permil(nsolute))
      if(declvar('phreeqmms', 'ch_basin_in_permil', 'nsolute', nsolute,&
           'double','Average permil value of isotope in precip ,'//&
           'ground water influx and irrigation', 'permil',&
           ch_basin_in_permil).ne.0) return

      allocate(ch_basin_out_permil(nsolute))
      if(declvar('phreeqmms', 'ch_basin_out_permil', 'nsolute', nsolute,&
           'double','Average permil value of isotope of stream '//&
           'at outlet','permil', ch_basin_out_permil).ne.0) return

      allocate(ch_basin_in_g(nsolute))
      if(declvar('phreeqmms', 'ch_basin_in_g', 'nsolute', nsolute,&
           'double','Total mass of solute into basin',&
           'g', ch_basin_in_g).ne.0) return

      allocate(ch_basin_out_g(nsolute))
      if(declvar('phreeqmms', 'ch_basin_out_g', 'nsolute', nsolute,&
           'double','Total mass of solute exported from basin',&
           'g', ch_basin_out_g).ne.0) return

      allocate(ch_basin_net_g(nsolute))
      if(declvar('phreeqmms', 'ch_basin_net_g', 'nsolute', nsolute,&
           'double','Net export of solute from basin.',&
           'g', ch_basin_net_g).ne.0) return

      allocate(ch_basin_in_gm2(nsolute))
      if(declvar('phreeqmms', 'ch_basin_in_gm2', 'nsolute', nsolute,&
           'double','Area-normalized flux of influents into basin.',&
           'g/m2', ch_basin_in_gm2).ne.0) return

      allocate(ch_basin_out_gm2(nsolute))
      if(declvar('phreeqmms', 'ch_basin_out_gm2', 'nsolute', nsolute,&
           'double','Area-normalized flux of effluent from basin.',&
           'g/m2', ch_basin_out_gm2).ne.0) return

      allocate(ch_basin_net_gm2(nsolute))
      if(declvar('phreeqmms', 'ch_basin_net_gm2', 'nsolute', nsolute,&
           'double','Net area-normalized flux of effluent from basin.',&
           'g/m2', ch_basin_net_gm2).ne.0) return
 
      allocate(snow_ion_pulse(nmru))
      if(declvar('phreeqmms', 'snow_ion_pulse', 'nmru', nmru,&
           'integer','Incongruous melting on this time step',&
           'none', snow_ion_pulse).ne.0) return

      allocate(ch_mru_vol_m3(nmru))
      if(declvar('phreeqmms', 'ch_mru_vol_m3', 'nmru', nmru,&
           'double','Total volume of water in the mru',&
           'm3', ch_mru_vol_m3).ne.0) return

      allocate(ch_mru_out_tempC(nmru))
      if(declvar('phreeqmms', 'ch_mru_out_tempC', 'nmru', nmru,&
           'double','Temperature of water leaving the mru',&
           'deg C', ch_mru_out_tempC).ne.0) return

      allocate(ch_mru_out_pH(nmru))
      if(declvar('phreeqmms', 'ch_mru_out_pH', 'nmru', nmru,&
           'double','pH of water leaving the mru',&
           'pH units', ch_mru_out_pH).ne.0) return

      allocate(ch_mru_mass_g(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_mass_g', 'nmru,nsolute',&
           nmru*nsolute,'double','Total mass of solute in mru',&
           'g', ch_mru_mass_g).ne.0) return

      allocate(ch_mru_rxn_g(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_rxn_g', 'nmru,nsolute', &
           nmru*nsolute,'double','Total mass of solute created '//&
           'or consumed by reactions in mru',&
           'g', ch_mru_rxn_g).ne.0) return

      allocate(ch_mru_conc_mgL(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_conc_mgL', 'nmru,nsolute',&
           nmru*nsolute,'double','Average concentration of '//&
           'solute in mru','mg/L', ch_mru_conc_mgL).ne.0) return

      allocate(ch_mru_in_mgL(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_in_mgL', 'nmru,nsolute',&
           nmru*nsolute,'double','Average solute concentration '//&
           'of precip, ground water influx, and irrigation',&
           'mg/L', ch_mru_in_mgL).ne.0) return

      allocate(ch_mru_out_mgL(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_out_mgL', 'nmru,nsolute',&
           nmru*nsolute,'double','Average solute concentration '//&
           'of mru discharge','mg/L', ch_mru_out_mgL).ne.0) return

      allocate(ch_mru_in_permil(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_in_permil', 'nmru,nsolute',&
           nmru*nsolute,'double','Average permil value of '//&
           'precip, ground water influx, and irrigation',&
           'permil', ch_mru_in_permil).ne.0) return

      allocate(ch_mru_out_permil(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_out_permil', 'nmru,nsolute',&
           nmru*nsolute,'double','Average permil value of '//&
           'mru discharge','permil', ch_mru_out_permil).ne.0) return

      allocate(ch_mru_in_g(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_in_g', 'nmru,nsolute',&
           nmru*nsolute,'double','Total mass of solute into mru',&
           'g', ch_mru_in_g).ne.0) return

      allocate(ch_mru_out_g(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_out_g', 'nmru,nsolute',&
           nmru*nsolute,'double','Total mass of solute '//&
           'exported from mru','g', ch_mru_out_g).ne.0) return

      allocate(ch_mru_net_g(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_net_g', 'nmru,nsolute',&
           nmru*nsolute,'double','Net export of solute from MRU.',&
           'g', ch_mru_net_g).ne.0) return

      allocate(ch_mru_in_gm2(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_in_gm2', 'nmru,nsolute',&
           nmru*nsolute,'double','Area-normalized flux of '//&
           'influents into mru.','g/m2', ch_mru_in_gm2).ne.0) return

      allocate(ch_mru_out_gm2(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_out_gm2', 'nmru,nsolute',&
           nmru*nsolute,'double','Area-normalized flux of '//&
           'effluent from mru.','g/m2', ch_mru_out_gm2).ne.0) return

      allocate(ch_mru_net_gm2(nmru,nsolute))
      if(declvar('phreeqmms', 'ch_mru_net_gm2', 'nmru,nsolute',&
           nmru*nsolute,'double','Net area-normalized flux '//&
           'of effluent from mru.','g/m2', ch_mru_net_gm2).ne.0) return

!
! parameters
!
      if(declparam('phreeqmms', 'xdebug_start', 'one', 'integer',&
         '0', '0', '1000000',&
         'Begin debugging on this step (0 if no debugging)',&
         'Begin debugging on this step (0 if no debugging)',&
         'integer').ne.0) return

      if(declparam('phreeqmms', 'xdebug_stop', 'one', 'integer',&
         '0', '0', '1000000',&
         'End debugging on this step',&
         'End debugging on this step',&
         'integer').ne.0) return

!      if(declparam('io', 'chemout_file_unit', 'one', 'integer',
!     +   '90', '50', '99',
!     +   'Unit number for file summarizing solute transport',
!     +   'Unit number for file summarizing solute transport',
!     +   'integer').ne.0) return

!
! A value of '0' for ppt_chem indicates that no time series of precipitation
! chemistry is available and therefore constant concentrations as defined by 
! init_soln_ppt will be used. If ppt_chem = 1 then a time series of precip
! chemistry will be read from the chemdat file.
!
      if(declparam('obs_chem', 'ppt_chem', 'one', 'integer',&
           '0', '0', '1',&
           '0 - Precipitation chemistry is constant as defined by '//&
           'init_soln_ppt and the pqi file; 1 - Precipitation '//&
           'chemistry is read every time step from the chemdep file',&
           '0 - Precipitation chemistry is constant as defined by '//&
           'init_soln_ppt and the pqi file; 1 - Precipitation '//&
           'chemistry is read every time step from the chemdep file',&
           'none') .ne.0) return

      if(declparam('obs_chem', 'chem_ext', 'one', 'integer',&
           '0', '0', '1',&
           '0 - External chemistry is constant as '//&
           'defined by init_soln_ext and the pqi file; '//&
           '1 - External chemistry is provided in '//&
           'the chemdep file for every time.',&
           '0 - External chemistry is constant as '//&
           'defined by init_soln_ext and the pqi file; '//&
           '1 - External chemistry is provided in '//&
           'the chemdep file for every time.',&
           'none') .ne.0) return

!
! parameters to describe water quality of irrigation and gw influx from outside
! of basin boundary. A maximum of one irrigation source and gw sources can be 
! defined for each MRU. GW sources could be a leaky irrigation canal crossing 
! the MRU and/or regional groundwater influx.
!

      allocate(src_ext_irrig(nmru))
      if(declparam('phreeqmms', 'src_ext_irrig', 'nmru',&
         'integer','0', 'bounded', 'nchem_ext',&
         'Chemical source for external irrigation',&
         'Chemical source for external irrigation',&
         'none').ne.0) return

      allocate(src_gw1(nmru))
      if(declparam('phreeqmms', 'src_gw1', 'nmru',&
         'integer','0', 'bounded', 'nchem_ext',&
         'Chemical source for 1st groundwater influx',&
         'Chemical source for 1st groundwater influx',&
         'none').ne.0) return

      allocate(src_gw2(nmru))
      if(declparam('phreeqmms', 'src_gw2', 'nmru',&
         'integer','0', 'bounded', 'nchem_ext',&
         'Chemical source for 2nd groundwater influx',&
         'Chemical source for 2nd groundwater influx',&
         'none').ne.0) return

!
! If there is an internal irrigation schedule for an MRU (irrig_sched_int>0),
! then irrig_int_src indicates if the MRU has recieves irrigations from a
! well in the MRU (value=0), or a stream segment (value > 0).
!
      ALLOCATE (irrig_int_src(nmru))
      if(declparam('irrig', 'irrig_int_src', 'nmru', 'integer',&
         '0', '0', '100',&
         ' 0 Irrigation from well in MRU; '//&
         '>0 Drainage segment ID that will provide irrigation water',&
         ' 0 Irrigation from well in MRU; '//&
         '>0 Drainage segment ID that will provide irrigation water',&
         'none')&
         .ne.0) return 

      if(declparam('phreeqmms', 'chem_sim', 'one', 'integer',&
         '1', '0', '1',&
         'Simulate solute fluxes (0=no; 1=yes)',&
         'Simulate solute fluxes (0=no; 1=yes)',&
         'integer').ne.0) return

      allocate(solnset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','solnset_table', 'nmru_res,nchem_sets',&
         'integer', '1', '1', '100',&
         'Solution sets available for initializing mru reservoir '//&
         'compositions',&
         'Solution sets for initializing mru reservoir compositions'//&
         'The solution numbers correspond to those read '//&
         'from the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(rxnset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','rxnset_table', 'nmru_res,nchem_sets',&
         'integer', '-1', '-1', '100',&
         'Sets of reaction IDs available for initializing MRU '//&
         'reservoir reactants (-1 if none).',&
         'Sets of reaction IDs available for initializing MRU '//&
         'reservoir reactants (-1 if none). '//&
         'The reactant IDs correspond to those read '//&
         'from the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(exchset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','exchset_table', 'nmru_res,nchem_sets',&
         'integer', '-1', '-1', '100',&
         'Sets of exchanger IDs available for initializing MRU '//&
         'reservoir exchange sites (-1 if none)',&
         'Sets of exchanger IDs available for initializing MRU '//&
         'reservoir exchange sites (-1 if none). '//&
         'The exchanger IDs correspond to those read '//&
         'from the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(surfset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','surfset_table', 'nmru_res,nchem_sets',&
         'integer', '-1', '-1', '100',&
         'Sets of surface species available for initializing MRU '//&
         'reservoirs (-1 if none)',&
         'Sets of surface species available for initializing MRU '//&
         'reservoir reactants (-1 if none). '//&
         'The suface species correspond to those read '//&
         'from the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(eq_phset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','eq_phset_table', 'nmru_res,nchem_sets',&
         'integer', '-1', '-1', '100',&
         'Sets of equilibrium phases available for initializing MRU '//&
         'reservoirs (-1 if none)',&
         'Sets of equilibrium phases available for initializing MRU '//&
         'reservoirs (-1 if none). '//&
         'The equilibrium phases correspond to those read '//&
         'from the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(kinset_table(nmru_res,nchem_sets))
      if(declparam('phreeqmms','kinset_table', 'nmru_res,nchem_sets',&
         'integer', '-1', '-1', '100',&
         'Sets of kinetic reactions available for initializing MRU '//&
         'reservoirs (-1 if none)',&
         'Sets of kinetic reactions available for initializing MRU '//&
         'reservoirs (-1 if none). '//&
         'The kinetics are defined in '//&
         'the phreeq input file, *.pqi.',&
         'none') .ne.0) return

      allocate(init_soln_ext(nchem_ext))
      if(declparam('phreeqmms','init_soln_ext', 'nchem_ext',&
         'integer', '1', '1', '100',&
         'Solution IDs to use for initializing the solute '//&
         'composition for each unique external source',&
         'Solution IDs to use for initializing the solute '//&
         'composition for each unique external source. '//&
         'The solution are described in the .pqi file.',&
         'none') .ne.0) return

      if(declparam('phreeqmms','init_soln_ppt', 'one',&
         'integer', '1', '1', '100',&
         'Solution ID describing typical precip chemistry.',&
         'Solution ID describing typical precip chemistry. '//&
         'The solution is described in the .pqi file.',&
         'none') .ne.0) return

      allocate(iso_n(nmru))
      if(declparam('phreeqmms','iso_n', 'nmru',&
         'double', '0.5', '0.5', '1.0',&
         'Isotopic fractionation factor "n"',&
         'Isotopic fractionation factor "n" of the '//&
         'Craig-Gordon evaporative fractionation model (1965)',&
         'none') .ne.0) return

      allocate(init_soln_mru(nmru))
      if(declparam('phreeqmms','init_soln_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Solution set to use for initializing the solute '//&
         'composition for each reservoir in an MRU',&
         'Solution set to use for initializing the solute '//&
         'composition for each reservoir in an MRU. The solution '//&
         'sets are described by the 2-dimensional solnset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_rxn_mru(nmru))
      if(declparam('phreeqmms','init_rxn_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Reactant set to use for initializing the reactants '//&
         'in each reservoir in an MRU',&
         'Reactant set to use for initializing the reactants '//&
         'in each reservoir in an MRU'//&
         'The reactant '//&
         'sets are described by the 2-dimensional rxnset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_exch_mru(nmru))
      if(declparam('phreeqmms','init_exch_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Exchanger set to use for initializing the solute '//&
         'composition for each reservoir in an MRU',&
         'Exchanger set to use for initializing the solute '//&
         'composition for each reservoir in an MRU. The exchanger '//&
         'sets are described by the 2-dimensional exchset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_surf_mru(nmru))
      if(declparam('phreeqmms','init_surf_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Surface-species set to use for initializing '//&
         'each reservoir in an MRU',&
         'Surface-species set to use for initializing '//&
         'each reservoir in an MRU'//&
         'The solution '//&
         'sets are described by the 2-dimensional surfset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_eq_ph_mru(nmru))
      if(declparam('phreeqmms','init_eq_ph_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Equilibrium-phase set to use for initializing '//&
         'each reservoir in an MRU',&
         'Equilibrium-phase species set to use for initializing '//&
         'each reservoir in an MRU'//&
         'The equilibrium-phase sets '//&
         'are described by the 2-dimensional eq_phset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_kin_mru(nmru))
      if(declparam('phreeqmms','init_kin_mru', 'nmru',&
         'integer', '1', '1', '100',&
         'Kinetic sets to use for initializing '//&
         'each reservoir in an MRU',&
         'Kinetic sets to use for initializing '//&
         'each reservoir in an MRU'//&
         'The kinetic '//&
         'sets are described by the 2-dimensional kinset_table '//&
         'parameter in conjuction with the .pqi file.',&
         'none') .ne.0) return

      allocate(init_soln_hydro(nhydro))
      if(declparam('phreeqmms','init_soln_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Solution IDs to use for initializing the solute '//&
         'composition for each stream segment or water body',&
         'Solution IDs to use for initializing the solute '//&
         'composition for each stream segment or water body. '//&
         'The solution are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(init_rxn_hydro(nhydro))
      if(declparam('phreeqmms','init_rxn_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Reactant IDs to use for initializing the reactants '//&
         'for each stream segment or water body (-1 if none)',&
         'Reactant IDs to use for initializing the reactants '//&
         'for each stream segment or water body (-1 if none). '//&
         'The solution are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(init_exch_hydro(nhydro))
      if(declparam('phreeqmms','init_exch_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Exchanger IDs to use for initializing exchange sites '//&
         'for each stream segment or water body (-1 if none)',&
         'Exchanger IDs to use for initializing exchange sites '//&
         'for each stream segment or water body (-1 if none). '//&
         'The exchangers are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(init_surf_hydro(nhydro))
      if(declparam('phreeqmms','init_surf_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Surface species to use for initializing '//&
         'each stream segment or water body (-1 if none)',&
         'Surface species to use for initializing '//&
         'each stream segment or water body (-1 if none). '//&
         'The surface species are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(init_eq_ph_hydro(nhydro))
      if(declparam('phreeqmms','init_eq_ph_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Equilibrium phases to use for initializing '//&
         'each stream segment or water body (-1 if none)',&
         'Equilibrium phases to use for initializing '//&
         'each stream segment or water body (-1 if none). '//&
         'The equilibrium phases are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(init_kin_hydro(nhydro))
      if(declparam('phreeqmms','init_kin_hydro', 'nhydro',&
         'integer', '1', '1', '100',&
         'Kinetics to use for initializing '//&
         'each stream segment or water body (-1 if none)',&
         'Kinetics to use for initializing '//&
         'each stream segment or water body (-1 if none). '//&
         'The kinetics are described in the .pqi file.',&
         'none') .ne.0) return

      allocate(C_ires(nchemvar))
      if(declparam('phreeqmms','C_ires', 'nchemvar', 'integer',&
         '6', '0', '99',&
         'Reservoir type (1-9, or 99) to save to ch_var',&
         'Reservoir type to save to ch_var: '//&
           ' 1= Canopy, '//&
           ' 2= Snowpack, '//&
           ' 3= Impermeable surface, '//&
           ' 4= O-horizon, '//&
           ' 5= Unsaturated zone, topo index '//&
           ' 6= Unsaturated zone, hill or riparian '//&
           ' 7= Unsaturated preferential flow, '//&
           ' 8= Saturated zone, '//&
           ' 9= Saturated preferential flow, '//&
           '99= Stream, pond, irrigation, or diversion.',&
         'none') .ne.0) return

      allocate(C_metric(nchemvar))
      if(declparam('phreeqmms','C_metric', 'nchemvar', 'integer',&
         '3', '1', '6',&
         'Chemical metric (1-6) for ch_var',&
         'Chemical metric for ch_var: '//&
           ' 1= Initial, '//&
           ' 2= Influent, '//&
           ' 3= Exfluent,  '//&
           ' 4= Reaction gain(+) or loss(-).'//&
           ' 5= Final, '//&
           ' 6= Net, ',&
         'none') .ne.0) return

      allocate(c_mru(nchemvar))
      if(declparam('phreeqmms','c_mru', 'nchemvar', 'integer',&
         '0', '0', '999',&
         'MRU for ch_var, 0 for basin sum for that reservoir',&
         'MRU for ch_var, 0 for basin sum for that reservoir',&
         'none') .ne.0) return

      allocate(c_rip(nchemvar))
      if(declparam('phreeqmms','c_rip', 'nchemvar', 'integer',&
         '0', '0', '2',&
         'Hillslope(0), Riparian(1), or Upland(2) for ch_var',&
         'Hillslope(0), Riparian(1), or Upland(2) for ch_var',&
         'none') .ne.0) return

      allocate(c_stindx(nchemvar))
      if(declparam('phreeqmms','c_stindx', 'nchemvar', 'integer',&
         '1', '0', '99',&
         'Specific topographic index for ch_var, 0 for average '//&
         'hillslope or riparian',&
         'Specific topographic index for ch_var, 0 for average '//&
         'hillslope or riparian',&
         'none') .ne.0) return

      allocate(c_hyd_indx(nchemvar))
      if(declparam('phreeqmms','c_hyd_indx', 'nchemvar', 'integer',&
         '1', '0', '99999',&
         'Specific drainage segment or node for ch_var',&
         'Specific drainage segment or node for ch_var. Unique '//&
         'IDs are assigned to each entity representing a stream, '//&
         'pond, irrigation source, or diversion.',&
         'none') .ne.0) return

      allocate(c_obs_indx(nchemvar))
      if(declparam('phreeqmms','c_obs_indx', 'nchemvar', 'integer',&
         '1', '0', '99',&
         'Index of solute chemistry observation station.',&
         'Index of solute chemistry observation station.',&
         'none') .ne.0) return

      allocate(C_units(nchemvar))
      if(declparam('phreeqmms','C_units', 'nchemvar', 'integer',&
         '7', '1', '13',&
         'Units for ch_var',&
         'Units for ch_var: '//&
         ' 1= mg, '//&
         ' 2= meq, '//&
         ' 3= mmol, '//&
         ' 4= mg/m², '//&
         ' 5= meq/m², '//&
         ' 6= mmol/m², '//&
         ' 7= mg/L, '//&
         ' 8= meq/L, '//&
         ' 9= mmol/L, '//&
         '10= convfactor(1),'//&
         '11= convfactor(2),'//&
         '12= convfactor(3),'//&
         '13= permil, '//&
         'Only mass or load units (1-6) allowed for '//&
         'Rx or Net metrics', 'none') .ne.0) return

!
! convfactor is available for up to 3 user-specified conversion factors
! to be used to convert reservoir concentrations to user-defined output
! units. The convfactor parameter may be used in either the .chemdat file
! to convert input units or here in the phreeqmms module to convert units
! for custom displays. For user-defined output units, specify 10, 11, or
! 12 in the c_units(nchemvar) param then specify the conversion factor
! in convfactor(1),(2), or (3).
!
      allocate(convfactor(nconvert))
      if(declparam('obs_chem', 'convfactor', 'nconvert', 'double',&
           '1', '0.000001', '1000000',&
           'User-defined conversion factor for solute '//&
           'inputs or outputs',&
           'User-defined conversion factors. Units are converted as '//&
           'follows: '//&
           'uconc(in)/convfactor = mol/L; or '//&
           'mol/L*convfactor = uconc(out).',&
           'user-defined') .ne.0) return

!
! Declare parameters established in other modules
!
! from basin_topg
!
      if(declparam('basin', 'basin_area', 'one', 'real',&
         '1.0', '0.01', '1e+09',&
         'Total basin area',&
         'Total basin area',&
         'km2').ne.0) return
      
      allocate(mru_area(nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',&
         '1.0', '0.01', '1e+09',&
         'MRU area',&
         'MRU area',&
         'km2').ne.0) return

      ALLOCATE (mru_elev(nmru))
      if(declparam('basin', 'mru_elev', 'nmru', 'real',&
         '0.', '-300.', '10000',&
         'Mean elevation for each MRU',&
         'Mean elevation for each MRU',&
         'meters').ne.0) return
     
!
! use cover density for canopy loads
!
      allocate(covden_sum(nmru))
      if(declparam('intcp', 'covden_sum', 'nmru', 'real',&
         '.5', '0.', '1.0',&
         'Summer vegetation cover density for major vegetation type',&
         'Summer vegetation cover density for the major '//&
         'vegetation type on each MRU',&
         'decimal percent')&
         .ne.0) return 

      allocate(covden_win(nmru))
      if(declparam('intcp', 'covden_win', 'nmru', 'real',&
         '.5', '0.', '1.0',&
         'Winter vegetation cover density for major vegetation type',&
         'Winter vegetation cover density for the major '//&
         'vegetation type on each MRU',&
         'decimal percent')&
         .ne.0) return


      allocate(nacsc(nmru))
      if(declparam('topc', 'nacsc', 'nmru', 'integer',&
         '1', '0', '100',&
         'Number of ln(a/tanB) increments in the subcatchment.',&
         'Number of ln(a/tanB) increments in the subcatchment.',&
         'none').ne.0) return

      allocate(ac(nac,nmru))
      if(declparam('topc', 'ac', 'nac,nmru', 'real',&
         '1', '0', '1',&
         'Fractional area for each ln(a/tanB) increment.',&
         'Fractional area for each ln(a/tanB) increment.',&
         'km2/km2').ne.0) return      

      allocate(qdffrac(nmru))
      if(declparam('topc', 'qdffrac', 'nmru', 'real',&
         '.3', '0', '1',&
         'Proportion of unsaturated zone drainage that runs off'//&
         ' as direct flow.','Fraction of unsaturated zone drainage'//&
         ' that runs off as direct flow.'//&
         'QDF=QDFFRAC*QUZ','Proportion')&
          .ne.0) return

      allocate(snow_ion_factor(nmru))
      if(declparam('topc', 'snow_ion_factor', 'nmru', 'real',&
         '1.0', '1.0', '10.0',&
         'Maximum ratio of solute concentrations in melt/pack', & 
         'Maximum ratio of solute concentrations in melt/pack', & 
         'ratio melt/pack').ne.0) return

      allocate(snowmelt_D_depl(nmru))
      if(declparam('topc', 'snowmelt_D_depl', 'nmru', 'real',&
         '-8.0', '-40', '0.0',&
         'Depletion of deuterium in melt versus pack', & 
         'Depletion of deuterium in melt versus pack', & 
         'permil per meltday').ne.0) return

      allocate(snowmelt_18O_depl(nmru))
      if(declparam('topc', 'snowmelt_18O_depl', 'nmru', 'real',&
         '-1.0', '-5.0', '0.0',&
         'Depletion of 18O in melt versus pack', & 
         'Depletion of 18O in melt versus pack', & 
         'permil per meltday').ne.0) return


! Other variables
      ALLOCATE(kPa(nmru))
      ALLOCATE (ar_fill(nhydro,nchan))
      ALLOCATE (snowcov_area(nmru))
      ALLOCATE (src_init(nsoln,nentity))
      ALLOCATE (chvar_lut(nchemvar,10))
      ALLOCATE (chvar_conv(nchemvar,nsolute+1))
      ALLOCATE (c_indx(nsoln,2))
      ALLOCATE (c_chem(nsoln))
      do i=1,nsoln
        ALLOCATE (c_chem(i)%M(nsolute,5))
        ALLOCATE (c_chem(i)%delta(nsolute,5))
      end do
      ALLOCATE (c_chem_uzgen(nmru))
      do i=1,nmru
        ALLOCATE (c_chem_uzgen(i)%M(nsolute,5))
        ALLOCATE (c_chem_uzgen(i)%delta(nsolute,5))
      end do
      ALLOCATE (c_chem_uzrip(nmru))
      do i=1,nmru
        ALLOCATE (c_chem_uzrip(i)%M(nsolute,5))
        ALLOCATE (c_chem_uzrip(i)%delta(nsolute,5))
      end do
      ALLOCATE (c_chem_uzup(nmru))
      do i=1,nmru
        ALLOCATE (c_chem_uzup(i)%M(nsolute,5))
        ALLOCATE (c_chem_uzup(i)%delta(nsolute,5))
      end do
      ALLOCATE (n_user(nentity))
      ALLOCATE (solns(nsoln))
      ALLOCATE (src(nsoln))
      ALLOCATE (srcdep(nresinp))
      ALLOCATE (dest(nsoln))
      ALLOCATE (fracs(nsoln))
      ALLOCATE (fracsdep(maxdep))
      ALLOCATE (chmru_soln(nmru))
      ALLOCATE (ch_rip_soln(nmru))
      ALLOCATE (ch_upland_soln(nmru))
      ALLOCATE (mru_in_vol(nmru))
      ALLOCATE (mru_out_vol(nmru))
      ALLOCATE (uz_spread(nentity))
      ALLOCATE (uzindxinit(nmru,nac,nentity))
      
      phreeqmms_decl = 0

      return
      end

!****************************************************************************
!
!     phreeqmms_init - initialize reservoir ids and solutions, and volumes
!
!****************************************************************************

      integer function phreeqmms_init()
#if defined(_WIN32)
      USE IFPORT
#endif
      USE WEBMOD_PHREEQ_MMS
      USE WEBMOD_OBSCHEM, ONLY :phq_lut,sol_id,sol_name,n_iso,iso_list
      USE WEBMOD_IO, only: chemout_file_unit


! Mixing variables from webmod_res
      USE WEBMOD_RESMOD, ONLY : vmix_can, vmix_snow, vmix_ohoriz, &
        vmix_uz, vmix_uz2can, vmix_uz2sat, vmix_sat2uz, vmix_uzgen,&
        vmix_uzrip, vmix_uzup, vmix_qdf, vmix_well, vmix_sat, vmix_satpref,&
        vmix_hill,  vmix_mru, vmix_hillexp, vmix_stream, vmix_diversion, &
        vmix_chan_loss, vmix_basin, uz2sat_vol, basin_qsim_cm
!      double precision vmix_imp(nmru,nresinp), vmix_rz(MAXMNR_3D) ! for later development
      USE WEBMOD_TOPMOD, ONLY : gw_loss,qpref_max, quz, st, riparian, uz_area
      IMPLICIT NONE
!      INCLUDE 'IPhreeqc.f90.inc'      
      integer, external ::  length
      interface
        function webmod_callback(x1, x2, str)
            double precision webmod_callback
            double precision, intent (in) :: x1
            double precision, intent (in) :: x2
            character(*), intent (in) :: str
        end function webmod_callback
      end interface

      logical filflg
      real dt

      integer k, l, ir, is, it, ih, path_len, file_len
      integer pqdat_len, res_id, ret, io, ivar, iresult
      integer uzwet, uzdry
!   moved to module
      !integer isoh1_len, isoh2_len, isogl_len, isogs_len, sol_h2_len
      !character*60 isogs, isogl
      !character*3000 sol_header1, sol_header2, iso_header1,iso_header2
      !character*3000 inp_dir, phreeqmms_pqi
      !character*3000 cheminit_file, readline, phreeq_database
      
      phreeqmms_init = 1

      if(getparam('phreeqmms', 'chem_sim', 1,&
           'integer',chem_sim) .ne.0) return

      if(getparam('phreeqmms', 'xdebug_start', 1,&
           'integer',xdebug_start) .ne.0) return

      if(getparam('phreeqmms', 'xdebug_stop', 1,&
           'integer',xdebug_stop) .ne.0) return

!
! Skip initialization (and later run) if no chemical simulations desired
!
      if(chem_sim.eq.1) then
!
      step1 = .true.
!
!     Phreeqc debug flag. If parameter x_debug_start.ge.1 then
!     print out mixing fractions, error files, etc during initialization
!     and then again on and after the time step x_debug_start.
!     If x_debug_start.eq.0 then never print out phreeq_debug information.
!
      if(xdebug_start.gt.0) then
         phr_tf=.true.
      else if(xdebug_start.eq.0) then
         phr_tf=.false.
      else
         print*,'parameter xdebug_start must be >= 0. Run terminated'
         return
      end if
      iresult = SetOutputFileOn(ID,phr_tf)
      iresult = SetErrorFileOn(ID,phr_tf)
      iresult = SetLogFileOn(ID,phr_tf)
      iresult = SetSelectedOutputFileOn(ID,phr_tf)

! get time step for reaction times
      dt = deltim()
      indx_rxn = 98 ! placekeeper for concentration of mixture after reactions

!
! Get user defined conversion factors
!
      if(getparam('obs_chem', 'convfactor', nconvert, 'double',&
           convfactor) .ne.0) return
!
! Get variables
!

      if(getvar('routec', 'clark_segs', 1, 'integer', clark_segs)&
         .ne.0) return

      if(getvar('routec', 'ar_fill', nhydro*nchan, 'real',&
           ar_fill) .ne.0) return

      if(getvar('obsc', 'chemdat_flag', 1, 'integer',&
           chemdat_flag) .ne.0) return

      if(chemdat_flag.eq.1) chemdat_exists = .true.

!
! Get parameters
!
!      if(getparam('io', 'chemout_file_unit', 1,&
!         'integer',chemout_file_unit ) .ne.0) return

      if(getparam('basin', 'basin_area', 1 , 'real', basin_area)&
         .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)&
         .ne.0) return

      if(getparam('basin', 'mru_elev', nmru, 'real', mru_elev)&
         .ne.0) return

      if(getparam('obs_chem', 'ppt_chem', 1, 'integer',&
           ppt_chem) .ne.0) return

      if(getparam('obs_chem', 'chem_ext', 1, 'integer',&
           chem_ext) .ne.0) return

      if(getparam('precip', 'src_ext_irrig', nmru, 'integer',&
           src_ext_irrig).ne.0) return

      if(getparam('precip', 'irrig_int_src', nmru, 'integer',&
           irrig_int_src).ne.0) return

      if(getparam('intcp', 'covden_sum', nmru, 'real', covden_sum)&
         .ne.0) return 

      if(getparam('intcp', 'covden_win', nmru, 'real', covden_win)&
         .ne.0) return 
       
      if(getparam('topc', 'ac', nac*nmru, 'real', AC)&
         .ne.0) return

      if(getparam('topc', 'src_gw1', nmru, 'integer', src_gw1)&
         .ne.0) return 

      if(getparam('topc', 'src_gw2', nmru, 'integer', src_gw2)&
         .ne.0) return 

      if(getparam('topc', 'qdffrac', nmru, 'real', qdffrac)&
         .ne.0) return

      if(getparam('phreeqmms', 'solnset_table', nmru_res*nchem_sets,&
           'integer', solnset_table) .ne.0) return

      if(getparam('phreeqmms', 'rxnset_table', nmru_res*nchem_sets,&
           'integer', rxnset_table) .ne.0) return

      if(getparam('phreeqmms', 'exchset_table', nmru_res*nchem_sets,&
           'integer', exchset_table) .ne.0) return

      if(getparam('phreeqmms', 'surfset_table', nmru_res*nchem_sets,&
           'integer', surfset_table) .ne.0) return

      if(getparam('phreeqmms', 'eq_phset_table', nmru_res*nchem_sets,&
           'integer', eq_phset_table) .ne.0) return

      if(getparam('phreeqmms', 'kinset_table', nmru_res*nchem_sets,&
           'integer', kinset_table) .ne.0) return

!      if(getparam('phreeqmms', 'riparian_thresh', nmru, 'real',&
!          riparian_thresh).ne.0) return

      if(getparam('phreeqmms', 'init_soln_ppt', 1, 'integer',&
           init_soln_ppt) .ne.0) return

      if(getparam('phreeqmms', 'init_soln_ext', nchem_ext,&
           'integer', init_soln_ext) .ne.0) return

      if(getparam('phreeqmms', 'init_soln_mru', nmru, 'integer',&
           init_soln_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_rxn_mru', nmru, 'integer',&
           init_rxn_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_exch_mru', nmru, 'integer',&
           init_exch_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_surf_mru', nmru, 'integer',&
           init_surf_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_eq_ph_mru', nmru, 'integer',&
           init_eq_ph_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_kin_mru', nmru, 'integer',&
           init_kin_mru) .ne.0) return

      if(getparam('phreeqmms', 'init_soln_hydro', nhydro, 'integer',&
           init_soln_hydro) .ne.0) return

      if(getparam('phreeqmms', 'init_rxn_hydro', nhydro, 'integer',&
           init_rxn_hydro) .ne.0) return

      if(getparam('phreeqmms', 'init_exch_hydro', nhydro, 'integer',&
           init_exch_hydro) .ne.0) return

      if(getparam('phreeqmms', 'init_surf_hydro', nhydro, 'integer',&
           init_surf_hydro) .ne.0) return

      if(getparam('phreeqmms', 'init_eq_ph_hydro', nhydro, 'integer',&
           init_eq_ph_hydro) .ne.0) return

      if(getparam('phreeqmms', 'init_kin_hydro', nhydro, 'integer',&
           init_kin_hydro) .ne.0) return

      if(getparam('phreeqmms', 'C_ires', nchemvar, 'integer',&
           c_ires) .ne.0) return

      if(getparam('phreeqmms', 'C_metric', nchemvar, 'integer',&
           c_metric) .ne.0) return

      if(getparam('phreeqmms', 'c_mru', nchemvar, 'integer',&
           c_mru) .ne.0) return

      if(getparam('phreeqmms', 'c_stindx', nchemvar, 'integer',&
           c_stindx) .ne.0) return

      if(getparam('phreeqmms', 'c_rip', nchemvar, 'integer',&
           c_rip) .ne.0) return

      if(getparam('phreeqmms', 'c_hyd_indx', nchemvar, 'integer',&
           c_hyd_indx) .ne.0) return

      if(getparam('phreeqmms', 'c_obs_indx', nchemvar, 'integer',&
           c_obs_indx) .ne.0) return

      if(getparam('phreeqmms', 'C_units', nchemvar, 'integer',&
           c_units) .ne.0) return

      if(getparam('basin', 'nacsc', nmru, 'integer', nacsc)&
         .ne.0) return

      if(getparam('phreeqmms', 'iso_n', nmru, 'double',&
           iso_n) .ne.0) return

      if(getparam('phreeqmms', 'snow_ion_factor', nmru, 'real',&
           snow_ion_factor) .ne.0) return

      if(getparam('phreeqmms', 'snowmelt_D_depl', nmru, 'real',&
           snowmelt_D_depl) .ne.0) return

      if(getparam('phreeqmms', 'snowmelt_18O_depl', nmru, 'real',&
           snowmelt_18O_depl) .ne.0) return


!
! Populate the chvar_lut table according to the parameters
! dimensioned by nchemvar
!
      do 22 i = 1, nchemvar
         chvar_lut(i,1) = -99
         chvar_lut(i,2) =&
              solnnum(0,c_ires(i),c_obs_indx(i),c_mru(i),c_stindx(i),&
              c_hyd_indx(i),c_rip(i))
         chvar_lut(i,3) = c_metric(i)
         chvar_lut(i,4) = c_units(i)
         chvar_lut(i,5) = c_ires(i)
         chvar_lut(i,6) = c_obs_indx(i)
         chvar_lut(i,7) = c_mru(i)
         chvar_lut(i,8) = c_stindx(i)
         chvar_lut(i,9) = c_hyd_indx(i)
         chvar_lut(i,10) = c_rip(i)  ! UZ: hillslope(0), riparian(1), or upland(2)
         
 22   continue

! Species and conversion units were read in obs_chem and are available in 
! phq_lut of type lut
!
! Open the initial solution file
      ret = getdataname (cheminit_file, '.pqi')
!
!  Kludge to run selected runs for GSA talk
!
!      cheminit_file='.\input\loch3.dat.pqi'
      open(unit=16, file=cheminit_file ,access='sequential',&
       form='formatted', status='old', iostat=io)
      file_len=index(cheminit_file,char(0))-1
      if(io.ne.0) then
         print*,'Cannot find the file of initial solutions ',&
              'and reservoir reactions: ', cheminit_file(1:file_len),&
              ' Run has been stopped.'
         return
      end if

!
! Prepend 'SELECTED OUTPUT' header to cheminit file and save
! in file, phreeqmms.pqi. This file will be created and used
! during each run.
!
      IF(control_string(inp_dir,'input_dir').NE.0) RETURN
      path_len = index(inp_dir,CHAR(0))-1   ! CHAR(0) is end of strings returned from control_string call
! Kludge for multi-CPU
      phreeqmms_pqi= inp_dir(1:path_len)//'phreeqmms.pqi'
      inquire(file=phreeqmms_pqi,exist=filflg)
      if (filflg) then
        open(unit=17,file=phreeqmms_pqi,status='old')
        close(unit=17,status='delete')
      endif
!----open the file.
      open (unit=17,file=phreeqmms_pqi,access='sequential',&
       form='formatted', status='new')

!      write(17,1000)(sol_name(i),i=1,nsolute)
! 1000 format('SELECTED_OUTPUT'/'-reset false'/'-pH'/'-temp'/&
!           '-tot  ',20(A10,1X))
!      iso_header1 = ""
!      iso_header2 = ""
! Create report of nsolute concentrations, as moles instead of molal to avoid conversion problems
! that arise from 18O being considered a solute
!
      write(sol_header1,*)(sol_name(i),i=1,nsolute)
      sol_header2 = ""
      i=1
      do while(i.le.nsolute)
        sol_h2_len = length(sol_header2)
        id_len = length(sol_name(i))
        if(sol_name(i)(1:id_len).eq."Alkalinity") then
          sol_header2 = sol_header2(1:sol_h2_len)//" ALK*TOT(""water"")"
        else
          sol_header2 = sol_header2(1:sol_h2_len)//" TOTMOL("""//sol_name(i)(1:id_len)//""")"
        endif
        i=i+1
      end do
!---- Append isotope ratios with user punch so that the permil 
!---- concentrations of the mixes, along with the calculated
!---- equilibrium fractionation between vapor and liquid or solid
!---- are sent back at the end of of the conc array.
      iso_header1 = sol_header1
      iso_header2 = sol_header2
      i = 1
!      if(n_iso.gt.0) then
      do while(i.le.n_iso)
       isoh1_len = length(iso_header1)
       isoh2_len = length(iso_header2)
       id_len = length(sol_name(iso_list(i)))
       isogl_len = length(phq_lut(sol_id(iso_list(i))%phq)%isofrac_gl)
       isogs_len = length(phq_lut(sol_id(iso_list(i))%phq)%isofrac_gs)
       isogl = phq_lut(sol_id(iso_list(i))%phq)%isofrac_gl(1:isogl_len)
       isogs = phq_lut(sol_id(iso_list(i))%phq)%isofrac_gs(1:isogs_len)
       iso_header1=iso_header1(1:isoh1_len)//' '//sol_name(iso_list(i)) &
        (1:id_len)//'_permil log_alpha_'//sol_name(iso_list(i &
         ))(1:id_len)// '_gl log_alpha_'//sol_name(iso_list(i &
         ))(1:id_len)//'_gs '
     
!    phq_lut(sol_id(iso_list(i))%phq)%isofrac_gl
!    phq_lut(sol_id(iso_list(i))%phq)%isofrac_gs
     
        iso_header2=iso_header2(1:isoh2_len)//' ISO("'//sol_name(iso_list(i)) &
        (1:id_len)//'") LK_NAMED("'//isogl(1:isogl_len)//'") LK_NAMED("'// &
        isogs(1:isogs_len)//'") '
        i=i+1
      end do
      isoh1_len = length(iso_header1)
      isoh2_len = length(iso_header2)
      write(17,1000)iso_header1(1:isoh1_len),iso_header2(1:isoh2_len)
 1000 format('SELECTED_OUTPUT'/'-reset false'/'-pH'/'-temp'/&
           'USER_PUNCH'/'-Headings ',A /'10 PUNCH ',A)
!      write(17,2000)iso_header1(1:isoh1_len),iso_header2(1:isoh2_len)
! 2000 format()
!     $ '10 PUNCH((TOT(""[18O]"")/TOT(""O"")) / 
!     $  2005.2e-6 - 1) * 1000 # VSMOW (Clark and Fritz, 1997)')")
!      end if
      do 20 while (io.eq.0)
         read(16,135,iostat=io,end=150,err=150)readline
         write(17,136)readline(1:length(readline))
 20   continue
 150  close(16)
      close(17)
      I = SYSTEM('type '//phreeqmms_pqi)
 135  format(A160)
 136  format(A)
!
! Open the geochemical database, phreeqc.dat  The phreeqc.dat file will be 
! placed in the data directory for now.
!
      IF(control_string(phreeq_database,'phreeq_database').NE.0) then
         PRINT *, 'phreeq_database needs to be defined in control file'
         STOP
      END IF
      pqdat_len = index(phreeq_database,CHAR(0))-1   ! CHAR(0) is end of strings returned from control_string call
!      phreeq_database = pqdat(1:pqdat_len)
      inquire(file=phreeq_database(1:pqdat_len),exist=filflg)
      if (filflg) then
         ID = CreateIPhreeqcMMS()
         iresult = LoadDatabase(ID, phreeq_database)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors loading database:'
            CALL OutputErrorString(id)
            STOP
         end if
        iresult = SetBasicFortranCallback(ID, webmod_callback)
      else
         print*,'Cannot open the phreeqc database ',&
              phreeq_database,' Run terminated.'
         return
      end if
!
! Allocate the conc as the number of solutes plus the number of isotopes
!
      ALLOCATE (conc(nsolute+n_iso*3))
!
      iphrq_mru = 1  ! point transpiration switch to MRU 1 for intial phreeqc callback function
!
! Provide initial solutions to phreeqc with the 
! initialization file, phreeqmms.pqi, that should
! now be available in the data directory
!
! Reserve solution 1 for pure water
!
      iresult = SetOutputFileOn(ID,phr_tf)
      iresult = SetErrorFileOn(ID,phr_tf)
      iresult = SetLogFileOn(ID,phr_tf)
      iresult = SetSelectedOutputFileOn(ID,phr_tf)
      iresult = RunFile(ID,phreeqmms_pqi)
      IF (iresult.NE.0) THEN
        PRINT *, 'Errors running solns.pqi:'
        CALL OutputErrorString(id)
        STOP
      ENDIF

      iresult = build_tally_table(ID)
      IF (iresult.NE.1) THEN
        PRINT *, 'Error building tally table:'
        CALL OutputErrorString(id)
        STOP
      ENDIF

      iresult = get_tally_table_rows_columns(ID,ntally_rows,ntally_cols)
      IF (iresult.NE.1) THEN
        PRINT *, 'Errors dimensioning rows and cols in tally table:'
        STOP
      ENDIF
      
! Allocate tally table arrays
      allocate (tally_table(ntally_rows + 1, ntally_cols))
      allocate (tally_row_label(ntally_rows + 1))
      allocate (tally_col_type(ntally_cols))
      allocate (tally_col_label(ntally_cols))
      ALLOCATE (n_ent(ntally_cols))
      allocate (mult(ntally_cols))
      sol_id%tally = -99 ! Initialize tally row of nsolute to -99
      do i = 1,ntally_rows
         iresult = get_tally_table_row_heading(ID,i,tally_row_label(i))
         IF (iresult.NE.1) THEN
            PRINT *, 'Errors configuring tally row headings:'
            STOP
         ENDIF
         do j = 1,nsolute
           if(sol_name(j).eq.tally_row_label(i))&
               sol_id(j)%tally = i
           if(sol_name(j).eq.'Alkalinity'.and.tally_row_label(i).eq.&
             'C')sol_id(j)%tally = i
         end do
      end do
      tally_row_label(ntally_rows + 1) = "Moles"
!
! Alert the user if solute not known to PHREEQC. i.e. not in any solution
! or pure phase in .pqi file and stop
!
      do j = 1,nsolute
        if(sol_id(j)%tally.eq.-99) then
           print *
           print *,'ERROR: Solute '// trim(sol_name(j)) //&
            'is not known to PHREEQC. Be sure that it described in '//&
              cheminit_file(1:file_len),&
            ' as one of the initial solutions or pure phases. ',&
            'Run Stopped.'
            STOP
        endif
      end do

!
! count how many of each entity (1, solution; 2, exchange, etc)
!
      do i = 1,ntally_cols
         n_ent(i) = 0
      end do

      do i = 1,ntally_cols
         iresult = get_tally_table_column_heading(ID,i,tally_col_type(i),&
             tally_col_label(i))
         IF (iresult.NE.1) THEN
            PRINT *, 'Errors configuring tally col headings:'
            STOP
         ENDIF
         n_ent(tally_col_type(i)) = n_ent(tally_col_type(i)) + 1
      end do
!
!  To mimimize sorting time in phreeqc, initialize all reservoirs in numerical order.
!
!  The run file above read solutions that should generally be less than 100.
!  The variable, nphrsolns will indicate the total number of solutions to be be modeled
!  Phreeq will also register the initial solutions read in the phreeqmms_pqi
!  file above, but other than solution 1, pure water, they will not be used in
!  the tables).
!
!  Store the DI water that is solution 1 as the first entry into the chemical tables
!
      src(1) = 1
      dest(1) = 1

      c_indx(1,1) = 1  ! order in c_chem tables
      c_indx(1,2) = 0  ! Is this solution tracked in one of the 10 chemvar variables (1=yes)
      nphrsolns = 1

!
!  In numerical order, solutions that may be active include the following:
!
!  chemical concentrations at a point (preciptation,irrigation,observations)
!  (1001, 1002, ...)
!     i = 1, nchem
!
      do 10001 i = 1, nchemdat
         nphrsolns = nphrsolns + 1
         src(nphrsolns) = 1
         dest(nphrsolns) =  solnnum(0,0,i,0,0,0,0)
         c_indx(nphrsolns,1) = dest(nphrsolns)
         c_indx(nphrsolns,2) =&
              chemflag(nphrsolns,chvar_lut,dest(nphrsolns),nchemvar)
10001 continue

!  Hillslope reservoirs at t=0 (100,000,000) and t=1 (200,000,000)
!    i = 0,1  start(0) or end(1) of time step
!     j = 1, nhillslope reservoirs (both real and temp; set at 14 for now)
!      k = 1, nmru
!       l = 1, nac (only for res_id = 7)
!        m = 1, nstat (last digit to be used for riparian or other discretization)
!                     (nstat not instituted yet - routine below, and the isoln
!                     function will have to be modified to include nstat)
!
! Also assign n_user reaction indices in preparation of initial mixes
!
!
      do 10050 i = 0,1
         do 10010 j = 1, 4 ! (canopy, snowpack, impermeable, and O-horizon)
            do 10010 k = 1, nmru
               nphrsolns = nphrsolns + 1
               src(nphrsolns) = 1
               dest(nphrsolns) =  solnnum(i,j,0,k,0,0,0)
               c_indx(nphrsolns,1) = dest(nphrsolns)
               c_indx(nphrsolns,2) =&
                    chemflag(nphrsolns,chvar_lut,&
                    dest(nphrsolns),nchemvar)
10010    continue
         do 10030 k = 1, nmru
            do 10030 l = 0, nacsc(k) 
               nphrsolns = nphrsolns + 1
               src(nphrsolns) = 1
               dest(nphrsolns) =  solnnum(i,5,0,k,l,0,0)
               c_indx(nphrsolns,1) = dest(nphrsolns)
               c_indx(nphrsolns,2) =&
                    chemflag(nphrsolns,chvar_lut,&
                    dest(nphrsolns),nchemvar)
               if(l.eq.0)ch_rip_soln(k)=nphrsolns ! zero index used psueudo solution of all riparian UZ (st >= riparian_thresh)
10030    continue
         do 10040 j = 6,14
            do 10040 k = 1, nmru
               nphrsolns = nphrsolns + 1
               src(nphrsolns) = 1
               dest(nphrsolns) =  solnnum(i,j,0,k,0,0,0)
               c_indx(nphrsolns,1) = dest(nphrsolns)
               c_indx(nphrsolns,2) =&
                    chemflag(nphrsolns,chvar_lut,&
                    dest(nphrsolns),nchemvar)
               if(j.eq.6.and.i.eq.0)ch_upland_soln(k)=nphrsolns ! Composite uz pseudo solution for UZ bins drier than riaparian_thresh
10040 continue
10050 continue
! 
!  Hydraulic reservoirs at t=0 (300,000,000) and t=1 (400,000,000)
!     i = 0,1
!      j = 1, nhydro (clark_segs in this initial existance)(res_id is
!                     equal to 99 for hydraulic segments)

      do 20050 i = 0,1
         do 20050 j = 0, clark_segs
            nphrsolns = nphrsolns + 1
            src(nphrsolns) = 1
            dest(nphrsolns) =  solnnum(i,99,0,0,0,j,0)
            c_indx(nphrsolns,1) = dest(nphrsolns)
            c_indx(nphrsolns,2) =&
                 chemflag(nphrsolns,chvar_lut,dest(nphrsolns),nchemvar)
20050 continue
!
! Basin and MRU reservoirs at t=0 (500,000,000). No need for t+1 since no mixing
!
! Note that the other than initilizing with DI, the basin and MRU 
! reservoirs are never processed by PHREEQC, they are just place holders
! in the c_chem matrix for summary calculations
! 

      nphrsolns = nphrsolns + 1
      src(nphrsolns) = 1
      dest(nphrsolns) =  solnnum(0,15,0,0,0,0,0)
      c_indx(nphrsolns,1) = dest(nphrsolns)
      c_indx(nphrsolns,2) =&
          chemflag(nphrsolns,chvar_lut,dest(nphrsolns),nchemvar)
      chbas_soln=nphrsolns
      indxb = nphrsolns
      do 30050 j = 1, nmru
        nphrsolns = nphrsolns + 1
        src(nphrsolns) = 1
        dest(nphrsolns) =  solnnum(0,15,0,j,0,0,0)
        c_indx(nphrsolns,1) = dest(nphrsolns)
        c_indx(nphrsolns,2) =&
          chemflag(nphrsolns,chvar_lut,dest(nphrsolns),nchemvar)
        chmru_soln(j)=nphrsolns
30050 continue
!
! The initial conditions for each reservoir are described with the
! init tables and init sets. These solutions and reaction entities are copied
! into src_init and individually renamed to the solnnum associated with that 
! nphrsoln row.
!
! src_init has eleven columns and nphrsoln rows. Not all entities
! are used in the current configuration of the model (see below).
! The solutions and other entities associated with a given reservoir 
! are unique but share the same ID number; that returned by the solnnum
! function. For example the solutions, reactants, surface reactors, etc,
! that are associated with the wettest unsaturated zone in the first MRU
! will have a number of 107001010. Before each mix that is
! static or a final export (imetric = -1 or 3), the fill_ent function will
! copy src_init values to the n_user vector. Each defined entity will evolve
! from its initial conditions during a run.
!
! Any non-negative value in columns 2 through 11 of n_user will indicate 
! inclusion in the geochemical mixing computed by PHREEQC; the first column
! is extraneous as the solutions being mixed are always included in each mix.
!
!           Entity                      n_user/keyword ! used
!           ==================          ==============   ====
! *         solution,                   1 ET_SOLUTION  ! explicit in mix call
! *         reaction,                   2 ET_REACTION    
! *         exchange,                   3 ET_EXCHANGE  
! *         surface,                    4 ET_SURFACE     
! *         gas_phase,                  5 ET_GAS_PHASE ! no
! *         equilibrium_phases,         6 ET_PURE_PHASE
! *         solid_solution,             7 ET_SS_PHASE  ! no
! *         kinetics,                   8 ET_KINETICS  
! *         mix,                        9 ET_MIX       ! no
! *         reaction_temperature        10 ET_TEMPERATURE ! explicit in mix call.
! *         unknown                     11 UnKnown     ! possible future use
!
! Populate UZ bin initialization for each entity by comparing initial 
! geochemistry (solutions, eq_ph, etc) for the wettest (ires=5) and driest (ires=6) bins.
!        - If both are -1, then all of that entity unknown (-1).
!          [uz_spread=0], *Not permitted for solutions*
!        - If they are the same and positive, all uz spaces initalized with that value.
!          [uz_spread=1]
!        - If they are different, use the wet value until the wetness index is less than 
!          riparian_thresh(nmru) after which point, the upland values are used (ires=6).
!          [uz_spread=2]
!        - If the value for the wet bin is not -1 and the dry bin is -1 then the initial
!          (soln, rxn, etc) will start with the value of the wet bin and increase by 1
!          for each nac bin (i.e. 10-19 or 30-39 for nacsc=10). For detailed assignments
!          like this, the user must verify that the entities are defined in the pqi file.
!          [uz_spread=3]
!
      do i=1,nmru
! solutions
        uzwet=solnset_table(5,init_soln_mru(i))
        uzdry=solnset_table(6,init_soln_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Solnset_table, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_SOLUTION)=0
            if(uzwet.gt.-1) uz_spread(ET_SOLUTION)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_SOLUTION)=1
          else 
            uz_spread(ET_SOLUTION)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_SOLUTION).eq.0) then  ! all solutions must be defined so this should never be true.
            uzindxinit(i,j,ET_SOLUTION)=-1
            PRINT *, 'Solnset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          elseif(uz_spread(ET_SOLUTION).eq.1) then
            uzindxinit(i,j,ET_SOLUTION)=uzwet
          elseif(uz_spread(ET_SOLUTION).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_SOLUTION)=uzwet
            else
              uzindxinit(i,j,ET_SOLUTION)=uzdry
            end if
          else if (uz_spread(ET_SOLUTION).eq.3) then
              uzindxinit(i,j,ET_SOLUTION)= uzwet+j-1  ! serial
          else
            PRINT *, 'Solnset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
           
        end do ! 1, nacsc
! reactions
        uzwet=rxnset_table(5,init_rxn_mru(i))
        uzdry=rxnset_table(6,init_rxn_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Exchset_table, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_REACTION)=0
            if(uzwet.gt.-1) uz_spread(ET_REACTION)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_REACTION)=1
          else 
            uz_spread(ET_REACTION)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_REACTION).eq.0) then 
            uzindxinit(i,j,ET_REACTION)=-1
          elseif(uz_spread(ET_REACTION).eq.1) then
            uzindxinit(i,j,ET_REACTION)=uzwet
          elseif(uz_spread(ET_REACTION).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_REACTION)=uzwet
            else
              uzindxinit(i,j,ET_REACTION)=uzdry
            end if
          else if (uz_spread(ET_REACTION).eq.3) then
              uzindxinit(i,j,ET_REACTION)= uzwet+j-1  ! serial
          else
            PRINT *, 'Exchset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
           
        end do ! 1, nacsc
! exchange
        uzwet=exchset_table(5,init_exch_mru(i))
        uzdry=exchset_table(6,init_exch_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Exchset_table, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_EXCHANGE)=0
            if(uzwet.gt.-1) uz_spread(ET_EXCHANGE)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_EXCHANGE)=1
          else 
            uz_spread(ET_EXCHANGE)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_EXCHANGE).eq.0) then
            uzindxinit(i,j,ET_EXCHANGE)=-1
          elseif(uz_spread(ET_EXCHANGE).eq.1) then
            uzindxinit(i,j,ET_EXCHANGE)=uzwet
          elseif(uz_spread(ET_EXCHANGE).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_EXCHANGE)=uzwet
            else
              uzindxinit(i,j,ET_EXCHANGE)=uzdry
            end if
          else if (uz_spread(ET_EXCHANGE).eq.3) then
              uzindxinit(i,j,ET_EXCHANGE)= uzwet+j-1  ! serial
          else
            PRINT *, 'Exchset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
           
        end do ! 1, nacsc
! surface complexation
        uzwet=surfset_table(5,init_surf_mru(i))
        uzdry=surfset_table(6,init_surf_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Surfset_table, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_SURFACE)=0
            if(uzwet.gt.-1) uz_spread(ET_SURFACE)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_SURFACE)=1
          else 
            uz_spread(ET_SURFACE)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_SURFACE).eq.0) then 
            uzindxinit(i,j,ET_SURFACE)=-1
          elseif(uz_spread(ET_SURFACE).eq.1) then
            uzindxinit(i,j,ET_SURFACE)=uzwet
          elseif(uz_spread(ET_SURFACE).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_SURFACE)=uzwet
            else
              uzindxinit(i,j,ET_SURFACE)=uzdry
            end if
          else if (uz_spread(ET_SURFACE).eq.3) then
              uzindxinit(i,j,ET_SURFACE)= uzwet+j-1  ! serial
          else
            PRINT *, 'Surfset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
        end do ! 1, nacsc
! equilibrium phases
        uzwet=eq_phset_table(5,init_eq_ph_mru(i))
        uzdry=eq_phset_table(6,init_eq_ph_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Eq_phtable, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_PURE_PHASE)=0
            if(uzwet.gt.-1) uz_spread(ET_PURE_PHASE)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_PURE_PHASE)=1
          else 
            uz_spread(ET_PURE_PHASE)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_PURE_PHASE).eq.0) then
            uzindxinit(i,j,ET_PURE_PHASE)=-1
          elseif(uz_spread(ET_PURE_PHASE).eq.1) then
            uzindxinit(i,j,ET_PURE_PHASE)=uzwet
          elseif(uz_spread(ET_PURE_PHASE).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_PURE_PHASE)=uzwet
            else
              uzindxinit(i,j,ET_PURE_PHASE)=uzdry
            end if
          else if (uz_spread(ET_PURE_PHASE).eq.3) then
              uzindxinit(i,j,ET_PURE_PHASE)= uzwet+j-1  ! serial
          else
            PRINT *, 'Eq_phset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
        end do ! 1, nacsc
! kinetics
        uzwet=kinset_table(5,init_kin_mru(i))
        uzdry=kinset_table(6,init_kin_mru(i))
          if(uzdry.lt.-1) then
            PRINT *, 'Kinset_table, row 6 (dry uz) must be '//&
                  'greater than or equal to -1. Run stopped'
            STOP
          endif
          if(uzdry.eq.-1) then
            if(uzwet.eq.-1) uz_spread(ET_KINETICS)=0
            if(uzwet.gt.-1) uz_spread(ET_KINETICS)=3
          elseif(uzdry.eq.uzwet) then
            uz_spread(ET_KINETICS)=1
          else 
            uz_spread(ET_KINETICS)=2
          end if
        do  j = 1, nacsc(i)
          if(uz_spread(ET_KINETICS).eq.0) then
            uzindxinit(i,j,ET_KINETICS)=-1
          elseif(uz_spread(ET_KINETICS).eq.1) then
            uzindxinit(i,j,ET_KINETICS)=uzwet
          elseif(uz_spread(ET_KINETICS).eq.2) then
            if(riparian(j,i))  then
              uzindxinit(i,j,ET_KINETICS)=uzwet
            else
              uzindxinit(i,j,ET_KINETICS)=uzdry
            end if
          else if (uz_spread(ET_KINETICS).eq.3) then
              uzindxinit(i,j,ET_KINETICS)= uzwet+j-1  ! serial
          else
            PRINT *, 'Kinset_table, row 5 (wet uz) has problems '//&
                 'Run stopped'
            STOP
          end if
        end do ! 1, nacsc
      enddo  ! 1, nmru
!
! populate src_init with initial solutions from pqi file
!
      do 40050 i = 1,nphrsolns
        k = isoln(c_indx(i,1),nchemdat,nmru,nac,clark_segs, &   ! returns the indices of the nphrsoln&
           ires,ichemdat,imru,inac,ihydro)
        do j = 1,11
           src_init(i,j) = -1  ! No reactions if not defined; solutions ignored
        end do
        if (i.eq.1) then  ! DI water
           src_init(i,1) = 1
        else if(ichemdat.eq.1) then ! Precip
           src_init(i,1) = init_soln_ppt
        else if(ichemdat.gt.1.and.ichemdat.le.(nchemdat-nchemobs))then
!                           Irrigation chemistry source. Ignore chemobs.
           src_init(i,ET_SOLUTION) = init_soln_ext(ichemdat-1)
        else if(ires.ge.1.and.ires.le.9) then ! hillslope reservoirs
           src_init(i,ET_SOLUTION) = &
               solnset_table(ires, init_soln_mru(imru))
           src_init(i,ET_REACTION) = &
               rxnset_table(ires,init_rxn_mru(imru))
           src_init(i,ET_EXCHANGE) = &
               exchset_table(ires,init_exch_mru(imru))
           src_init(i,ET_SURFACE) = &
                surfset_table(ires,init_surf_mru(imru))
           src_init(i,ET_PURE_PHASE) = &
                eq_phset_table(ires,init_eq_ph_mru(imru))
           src_init(i,ET_KINETICS) = &
                kinset_table(ires,init_kin_mru(imru))
!
!           This block sets initial geochemistry for the unsaturated zone bins.
!
           if(inac.gt.0.) then
            src_init(i,ET_SOLUTION)=uzindxinit(imru,inac,ET_SOLUTION)
            src_init(i,ET_REACTION)=uzindxinit(imru,inac,ET_REACTION)
            src_init(i,ET_EXCHANGE)=uzindxinit(imru,inac,ET_EXCHANGE)
            src_init(i,ET_SURFACE)=uzindxinit(imru,inac,ET_SURFACE)
            src_init(i,ET_PURE_PHASE)=&
                                    uzindxinit(imru,inac,ET_PURE_PHASE)
            src_init(i,ET_KINETICS)=uzindxinit(imru,inac,ET_KINETICS)
           end if
         else if(ihydro.gt.0) then  ! Draininge segment
           src_init(i,ET_SOLUTION) = init_soln_hydro(ihydro)
           src_init(i,ET_REACTION) = init_rxn_hydro(ihydro)
           src_init(i,ET_EXCHANGE) = init_exch_hydro(ihydro)
           src_init(i,ET_SURFACE) = init_surf_hydro(ihydro)
           src_init(i,ET_PURE_PHASE) = init_eq_ph_hydro(ihydro)
           src_init(i,ET_KINETICS) = init_kin_hydro(ihydro)
         else
            print*,'Reservoir ',c_indx(i,1),' row ',i,&
                   ' was not initialized.'
         end if
40050 continue
!$$$ *                 n_user/n_ent index
!$$$ *  Phreeq KEYWORD                    ET_Keyword  ! used in n_user 
!$$$ *         solution,                   1 Solution    ! no
!$$$ *         reaction,                   2 Reaction    
!$$$ *         exchange,                   3 Exchange    
!$$$ *         surface,                    4 Surface     
!$$$ *         gas_phase,                  5 Gas_phase   
!$$$ *         equilibrium_phases,         6 Pure_phase  
!$$$ *         solid_solution,             7 Ss_phase    
!$$$ *         kinetics,                   8 Kinetics    
!$$$ *         mix,                        9 Mix         ! nope
!$$$ *         reaction_temperature        10 Temperature 
!$$$ *         unknown                     11 UnKnown    ! not likely
!$$$ *
!$$$ */
! Use phr_multicopy here to initialize solution IDs.
! Use phr_mix subsequently so that concentrations are returned for
! initializing and updating the c_chem array
!
!$$$         iresult = phr_multicopy('solution',src, dest, nphrsolns)
         iresult = phr_multicopy(ID,'solution',src, dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy: solution'
            CALL OutputErrorString(id)
            STOP
         ENDIF

      print*,'Number of solutions tracked for this run: ',nphrsolns

         iresult = phr_multicopy(ID,'reaction',src_init(1,2),&
           dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy: reaction'
            CALL OutputErrorString(id)
            STOP
         ENDIF

         iresult = phr_multicopy(ID,'exchange',src_init(1,3),&
           dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy: exchange'
            CALL OutputErrorString(id)
            STOP
         ENDIF

         iresult = phr_multicopy(ID,'surface',src_init(1,4),&
           dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy: surface'
            CALL OutputErrorString(id)
            STOP
         ENDIF

         iresult = phr_multicopy(ID,'equilibrium_phases',src_init(1,6),&
           dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy:eq_phases'
            CALL OutputErrorString(id)
            STOP
         ENDIF

         iresult = phr_multicopy(ID,'kinetics',src_init(1,8),&
           dest, nphrsolns)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_multicopy: kinetics'
            CALL OutputErrorString(id)
            STOP
         ENDIF

      write(chemout_file_unit,'(A)')'row, chemrow, src, dest, '//&
                'chemvar?, src_init->'
      write(chemout_file_unit,123)(i,isoln(dest(i),nchemdat,nmru,nac,&
           clark_segs,ires,ichemdat,imru,inac,ihydro),&
           src(i),dest(i),c_indx(i,2),(src_init(i,j),j=1,11)&
              ,i=1,nphrsolns)
 123  format(16I10)
 
!
! update the src_init array with the unique reaction ids equal to the c_indx
!
      do 50050 i = 1,nphrsolns
         do j = 2,11  ! ignore solutions (solutions are j=1)
            if (src_init(i,j).ne.-1) src_init(i,j) = c_indx(i,1)
         end do
50050 continue

! Initialize c_chem volumes with 1 liter, Temp=25degC, ph=7, del value -1000 and solute mass = zero, 
      do i = 1,nphrsolns
         c_chem(i)%vol(init)=1.0e-3
         c_chem(i)%Temp(init)=25.D0
         c_chem(i)%pH(init)=7.D0
         do j = 1, nsolute
            do k = 1,5 ! Init, in, out, rxn|ET, final
              c_chem(i)%M(j,k) = 0.D0
              c_chem(i)%delta(j,k) = -1000.D0
            end do
          end do
      end do
! Check that all ch_vars linked to a known solution and assign factors
! to convert c_chem values to desired ch_var variables
!

      do 23 ivar = 1, nchemvar
         if (chvar_lut(ivar,1).eq. -99) then
            print*,'Could not find a solution for ch_var : ',ivar,&
                 '. Run terminated'
            return
         end if
         indx = chvar_lut(ivar,1)
         imet = chvar_lut(ivar,3)
         iunit = chvar_lut(ivar,4)
         unit_type = int((iunit-1)/3)+1 ! 1,mass; 2,load; 3,std conc; 4,user conc; 5,permil
         ires = chvar_lut(ivar,5)
         ichemdat = chvar_lut(ivar,6)
         imru = chvar_lut(ivar,7)
         inac = chvar_lut(ivar,8)
         ihydro = chvar_lut(ivar,9)
         irip = chvar_lut(ivar,10)

         if(imet.lt.1.or.imet.gt.6)then
            print*,'c_metric indicated in ch_var ',&
                 ivar,' is not valid. Select from: 1)  ',&
                 'Run terminated.'
            return
         else if (iunit.gt.6.and.imet.eq.4.or.&
                 iunit.gt.6.and.imet.eq.6) then
            print*,'Reaction or net metric for ',&
                 'ch_var ',ivar,' not permitted ',&
                 'since the c_units is a concentration. ',&
                 'Run terminated.'
            return
         else if (iunit.lt.1.or.iunit.gt.13) then
            print*,'Illegal units indicated for ch_var ',&
                 ivar,'. Run terminated'
            return
         else                   !Initial validation OK so proceed to
!                               !calculate partial conversion factors

!     Compute partial denominator to apply if a load. If the unit type
!     is any other than a load (mass or concentration) then the values
!     can be computed using the mass conversions in the numerator and
!     the data already in the c_chem matrix.
!     

            if(unit_type.eq.2) then ! load, so figure which area to apply
!
! if impermeable, use perc impermeable area (*mru_area) (implement later)
! if canopy, use mru_area and then multiply by canopy density later
! if snowpack, use mru_area and then multiply by snow covered area later
! if unsaturated zone, use fractional area (*mru) or composite areas if ires equals 6
! all other hillslope reservoirs use mru_area
! if stream segment or zero area, write error and set conversion
!   one square meter
!
! Note that for loads for canopy and snowpack will depend on the cover
! density or snow covered area on a given day. Therefore, the areal
! conversion factor, chvar_conv(*,nsolute+1) will be set to the mru area here
! in the init section and then adjusted in the run section before 
! assigning the values to the chvar variables.
!
               if (ires.eq.5) then ! unsaturated (5) or composite loni bin (6)
                  if(inac.eq.nacsc(imru)) then
                    ACF=0.5*AC(inac,imru)
                  else
                    ACF=0.5*(AC(inac,imru)+AC(inac+1,imru))
                  endif
                  print*,'chvar number ',ivar, 'uses an area of ',acf*mru_area(imru), &
                         ' sq km to compute loads for TWI ', inac
                  chvar_conv(ivar,nsolute+1) = acf*mru_area(imru)*a_million
               else if (ires.eq.6) then
                 if(irip.eq.0) then
                   chvar_conv(ivar,nsolute+1) = mru_area(imru)*a_million
                   print*,'chvar number ',ivar, 'uses the area of mru ',imru,' for the UZ area.'
                 elseif(irip.eq.1) then
                   chvar_conv(ivar,nsolute+1) = uz_area(1,imru)*a_million
                   print*,'chvar number ',ivar, 'uses a riparian area of ',uz_area(1,imru),' sq km to compute loads.'
                 elseif(irip.eq.2) then
                   chvar_conv(ivar,nsolute+1) = uz_area(2,imru)*a_million
                   print*,'chvar number ',ivar, 'uses an upland area of ',uz_area(2,imru),' sq km to compute loads.'
                 else
                   print*,'chvar number ',ivar, 'has c_ires of 6 so c_rip must equal 0,1,or 2, fix and restart'
                   return
                 endif
               else if (ires.eq.1) then ! canopy.
                  print*,'Loads calculated using canopy area ',&
                       'for chvar number ',ivar
                  chvar_conv(ivar,nsolute+1) =&
                       mru_area(imru)*a_million
               else if (ires.eq.2) then  ! Snowpack
                  print*,'Loads calculated using snow-covered area ',&
                       'for chvar number ',ivar
                  chvar_conv(ivar,nsolute+1) =&
                       mru_area(imru)*a_million
               else if (ires.lt.16) then ! all reservoirs other but streams
                 if(ires.eq.15.and.imru.eq.0) then  ! basin (should be done with basin variables but this allows a double check)
                  print*,'Loads calculated using basin area ',&
                       'for chvar number ',ivar
                  chvar_conv(ivar,nsolute+1) = basin_area*a_million
                 else
                  print*,'Loads calculated using mru area ',&
                       'for chvar number ',ivar
                  chvar_conv(ivar,nsolute+1) = mru_area(imru)*a_million
                 end if
               end if
               if(chvar_conv(ivar,nsolute+1).eq.0.or.ires.eq.99) then
                  chvar_conv(ivar,nsolute+1) = 1.0
                  print*,'Stream segment or zero area indicated ',&
                       'in chvar ',ivar,&
                       '. Loads computed using 1 sq.meter.'
               end if
            else                ! mass or concentration
               chvar_conv(ivar,nsolute+1) = 1.0
            end if
!     Compute numerators for unit conversions reflecting mass
!     and charge of each solute
            do 10 i = 1,nsolute
               if(unit_type.lt.4.or.unit_type.eq.5) then ! standard units or permil
                  if(mod((iunit+2),3).eq.0) then ! mg
                     chvar_conv(ivar,i) = phq_lut(sol_id(i)%phq)%M2mg
                  else if (mod((iunit+1),3).eq.0) then ! meq
                     chvar_conv(ivar,i) = phq_lut(sol_id(i)%phq)%M2meq
                  else          ! mmol
                     chvar_conv(ivar,i) = M2mM
                  end if
               else             ! custom conversions (10-12)
                  chvar_conv(ivar,i) = convfactor(iunit-9)
               end if
!             (iunit 13, permil has no conversion, it is reported directly)               
 10         continue
         end if
 23   continue



!
!      n_user(ET_SOLUTION) = -1

!      Initialize reaction parameters to 'mixing only'

!      do k=1,11
!         n_user(k) = -1
!      end do
!
      rxnmols = 0.0   ! No dry deposition
!      tempc = -1    ! These are returned values, not inputs. use TEMPERATURE block to assign temperature.
!      ph = 7.0      ! pH is initialized in the PQI solutions and updated with each mix and reaction.
      tsec = dt*3600   ! dt in hours
      fill_factor = 1.0
      iphrq_mru = 1    ! need at least one index for transpiration switch before phrees_mix is called

!
! Initialize nchemdat reservoirs (precip, irrigation).
! If no chemdat file was found these will be the concentrations
! used for the entire run, otherwise they will be overwritten
! with the time series observations in the chemdat file.
! All MRUs receive precipitation of the same concentration
!
      src(1) = init_soln_ppt
!     solnnum(time,reservoir ID,chemdat, mru,nac,hydro,stat)
      dest(1) = solnnum(0,0,1,0,0,0,0)
      fracs(1) = 1.0         
      totvol = 1D-3             ! in cubic meters

      iresult = fill_ent(n_user,dest(1),nchemdat,nmru,&
           nac,clark_segs,src_init)
      if(iresult.ne.0) then
         PRINT *, 'Errors assigning phreeq entities:'
         CALL OutputErrorString(id)
         STOP
      end if
! debug
!      conc = 0.0
!      iresult = test(conc)
      
!      write(chemout_file_unit,'(A)')'is, mixture, conc,' 
      iresult = phr_mix(ID,1, src(1), fracs,  dest(1),&
           fill_factor, dest(1), conc, phr_tf, n_user,rxnmols,tempc,&
           ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

      IF (iresult.NE.0) THEN
         PRINT *, 'Errors during phr_mix:'
         CALL OutputErrorString(id)
!!         STOP
      ENDIF
! debug
!      iresult = test(conc)

!
! Typical concentrations for each external source (irrigation and/or groundwater)
!
      do 182 ir = 1,nchem_ext
         src(1) = init_soln_ext(ir)
!     solnnum(time,reservoir ID,chemdat, mru,nac,hydro,stat)
         dest(1) = solnnum(0,0,ir+1,0,0,0,0) 
         fracs(1) = 1.0
         totvol = 1D-3 ! in cubic meters

!         iresult = fill_ent(n_user,dest(1),nchemdat,nmru,
!     $        nac,clark_segs,src_init)
!         if(iresult.ne.0) then
!            PRINT *, 'Errors assigning phreeq entities:'
!            CALL OutputErrorString(id)
!            STOP
!         end if
      
         iresult = phr_mix(ID,1, src(1), fracs,  dest(1),&
              fill_factor, dest(1), conc, phr_tf, no_rxn,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
 182  continue
 
! ET as DI (soln 1) included to be complete. Let all initial volumes = 1 kilogram (liter)
!
! Be careful not to assign anyting other than isotopic values for D and [18O] to solution 1
! (ET) or unexpected behavior such as a mass imbalance could result.
! If a stable isotope (D or 18O)is included as one of the nsolutes of interest, then
! evaporated water of specific delta value will be created in the fractionate subroutine
! and the isotopes will be tracked as basin exports in the update_chem routine.
!
         indx = 1
     
         totvol = 1D-3 ! in cubic meters
         c_chem(indx)%vol(init)=totvol
         fracs(1) = 1.0
         solns(1) = 1

!         iresult = fill_ent(n_user,dest(2),nchemdat,nmru,
!     $        nac,clark_segs,src_init)
!         if(iresult.ne.0) then
!            PRINT *, 'Errors assigning phreeq entities:'
!            CALL OutputErrorString(id)
!            STOP
!         end if

         iresult = phr_mix(ID,1, solns, fracs,  1,&
              fill_factor, 1, conc, phr_tf, no_rxn,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF

! Set basin index number
      indxb = chbas_soln
!
! Initialize individual hillslope reservoirs
!
      do 30 is = 1, nmru
! Establish MRU airpressure in kilopascals to derive wetbulb temperature for precipitation.
          kPa(is) = 101.325 * (1 - 2.25577e-5*mru_elev(is))**5.25588
!
! get the solution number for the pseudo-reservoir summarizing hillslope
! and unsaturated zone chemistry
!
         indxm = chmru_soln(is)
         iphrq_mru = is
!
!  Initialize hillslope reservoirs using the solnset_table parameter,
!  solnset_table(nmru_res,nmru) where reservoir ID:
!    1 = canopy interception. Note that the if there is any canopy
!        during winter or summer, then model maintains a residual
!        water content on the canopy at all times,
!    2 = snowpack (can melt completely),
!    3 = Impermeable surface (not implemented yet),
!    4 = O-horizon (fixed volume with seasonal canopy storage),
!    5 = unsaturated zone composition in wettest topographic index bin,
!    6 = unsaturated zone composition in driest topographic index bin,
!    7 = preferential flow through the unsaturated zone,
!    8 = saturated zone,
!    9 = preferential flow through the saturated zone.
!
!  Note special use of res_id 5 and 6, which represent the wet
!  and dry endmembers of the unsaturated zone reservoirs in the hillslope.
!  When the user is choosing a specific UZ reservoir for analysis using
!  the chem_var dimension, it will always be res_id 7 with further qualification
!  by nac. Res_id 8 will hold average values for hillslope, riparian
!  or upslope UZ reservoirs. Initial input concentrations are defined with
!  init_soln_ppt and init_soln_ext. Each MRU can then point to one of these
!  external concentrations using src_ext_irrig, src_gw1, or src_gw2.
!
! These reservoirs will also contains reactants and surfaces to
! simulate equilibrium and kinetic reactions as the nine indices only
! differ slightly from the solutions above in that the solutions between the
! the wet and dry wetness indices as described above.
!
!    1 = Canopy
!    2 = Snowpack
!    3 = Impermeable surface (not implemented)
!    4 = O-horizon
!    5 = Wettest unsaturated zone 
!    6 = Driest unsaturated zone
!    7 = Preferrential flow through the unsaturated zone (qdf)
!    8 = Saturated zone
!    9 = Preferential flow through the saturated zone (qpref)
!
! Geochemical activity in each of these reactive reservoirs is 
! described by a unique combination of reactants (i.e. dry deposition),
! exchange ions and exchangers (i.e. Ca and Mg on clay surfaces),
! surface species (i.e.. ferric oxide), equilibrium phases
! (i.e. mineral weathering), and kinetics (i.e. transformation
! rates of ammonia to nitrite)
!

! Initialize remaining reservoirs
!
         do 32 res_id = 1, 9
!
!              The solnset_table lists the initial solution sets for each of the 9
!              reservoir in a hillslope (see discussion of reservoir ID above).
            src(res_id) = solnset_table(res_id, init_soln_mru(is))
!              solnnum(time,reservoir ID,chemdat, mru,nac,hydro,stat)
            dest(res_id) = solnnum(0,res_id,0,is,0,0,0)
 32      continue



! Canopy
   
         indx = isoln(dest(1),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_can(is,1)
         c_chem(indx)%vol(init)=totvol
         
         solns(1) = src(1)

         iresult = fill_ent(n_user,dest(1),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(1),&
              fill_factor,  dest(1), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!
! Snowpack
!
         indx = isoln(dest(2),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_snow(is,1)
         c_chem(indx)%vol(init)=totvol
 
         solns(1) = src(2)

         iresult = fill_ent(n_user,dest(2),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(2),&
              fill_factor, dest(2), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!
! O-horizon
!
         indx = isoln(dest(4),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_ohoriz(is,1)
         c_chem(indx)%vol(init)=totvol
 
         solns(1) = src(4)
!         n_user(ET_REACTION) = 1
!         rxnmols = 1e-6

         iresult = fill_ent(n_user,dest(4),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(4),&
              fill_factor, dest(4), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
!         rxnmols = 0.0
!         n_user(ET_REACTION) = -1

! Update volume/mole matrix after a mix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!
! Unsaturated zone
!
!   The unsaturated zone geochemistry was assigned from pqi inputs
!   using the solution number of the wet (mru_res=5) and dry (mru_res=6)
!   endmember solutions as described above.
!
! Initialize 1 kg of solution for each UZ
!
         do 37 it = 1,nacsc(is)
!              solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)
! assign uzindx to riparian or upland
            if(riparian(it,is))  then
              indxuz=ch_rip_soln(is)
            else
              indxuz=ch_upland_soln(is)
            end if
            mixture = solnnum(0,5,0,is,it,0,0)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            solns(1) = uzindxinit(is,it,ET_SOLUTION)
            totvol = vmix_uz(it,is,1)
            c_chem(indx)%vol(init)=totvol
            fracs(1) = 1.0
            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if
! Check for negative fractions

               iresult=checkfracs(1,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,1, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, n_user,rxnmols,&
                 tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
            
! Update volume/mole matrix after a mix
            iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
! Accumulate initial volume of composite uz            
 37      end do
!     
!     Test the selected output for last loni
!                                !! output solution described in selected output
!         cols = GetSelectedOutputColumnCount()
!         PRINT *, 'Unsaturated zone solution for MRU: ',is,
!     $        ' and loni 30'
!            iresult = GetSelectedOutputValue
!     $           (0, 1, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), ph
!         
!            iresult = GetSelectedOutputValue
!     $           (0, 2, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), tempc
!         DO 80 i = 3,cols
!            iresult = GetSelectedOutputValue
!     $           (0, i, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), conc(i-2)
! 80      CONTINUE
!
! Preferential flow through the unsaturated zone (qdffrac)
!
         indx = isoln(dest(7),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_qdf(is,1)
         solns(1) = src(7)

         iresult = fill_ent(n_user,dest(7),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(7),&
              fill_factor, dest(7), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix
         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!
! Saturated zone
!
         indx = isoln(dest(8),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_sat(is,1)
         c_chem(indx)%vol(init)=totvol         
         solns(1) = src(8)

         iresult = fill_ent(n_user,dest(8),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(8),&
              fill_factor,  dest(8), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!************** print start *********
!     
!     output solution described in selected output
!         cols = GetSelectedOutputColumnCount()
!         PRINT *, 'Saturated zone solution for MRU: ',is
!            iresult = GetSelectedOutputValue
!     $           (0, 1, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), ph
!         
!            iresult = GetSelectedOutputValue
!     $           (0, 2, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), tempc
!         DO 80 i = 3,cols
!            iresult = GetSelectedOutputValue
!     $           (0, i, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), conc(i-2)
! 80      CONTINUE
!
!************** print stop **********
!
! Preferential flow through the saturated zone
!

         indx = isoln(dest(9),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_satpref(is,1)
         c_chem(indx)%vol(init)=totvol
         solns(1) = src(9)

         iresult = fill_ent(n_user,dest(9),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(9),&
              fill_factor, dest(9), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
 30   continue
!
! Initialize stream chemistry
!
 
      do 45 ih = 1,clark_segs
         dest(ih) = solnnum(0,99,0,0,0,ih,0)
         indx = isoln(dest(ih),nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         totvol = vmix_stream(ih)
         c_chem(indx)%vol(init)=totvol
         src(ih) =  init_soln_hydro(ih)
!                   solnnum(time,reservoir ID,chemdat, mru,nac,hydro,stat)
         solns(1) = src(ih)

         iresult = fill_ent(n_user,dest(ih),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

         iresult = phr_mix(ID,1, solns, fracs,  dest(ih),&
              fill_factor, dest(ih), conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
! Update volume/mole matrix after a mix: + basin

         iresult = update_chem(indx,totvol,1,conc,pH,tempc,&
                               tally_table,n_ent,0,indxb,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
 45   continue

!$$$c
!$$$c Initialize ch_vars to zero
!$$$      ch_var_01_m3 = 0
!$$$      ch_var_02_m3 = 0
!$$$      ch_var_03_m3 = 0
!$$$      ch_var_04_m3 = 0
!$$$      ch_var_05_m3 = 0
!$$$      ch_var_06_m3 = 0
!$$$      ch_var_07_m3 = 0
!$$$      ch_var_08_m3 = 0
!$$$      ch_var_09_m3 = 0
!$$$      ch_var_10_m3 = 0
!$$$      do 2 isol = 1,nsolute
!$$$      ch_var_01_sol(isol) = 0
!$$$      ch_var_02(isol) = 0
!$$$      ch_var_03(isol) = 0
!$$$      ch_var_04(isol) = 0
!$$$      ch_var_05(isol) = 0
!$$$      ch_var_06(isol) = 0
!$$$      ch_var_07(isol) = 0
!$$$      ch_var_08(isol) = 0
!$$$      ch_var_09(isol) = 0
!$$$      ch_var_10(isol) = 0
!$$$ 2    continue
!
!
! Transfer values from local variables
!
!$$$      iresult = chem2var(
!$$$     $     c_chem,chvar_lut,chvar_conv,nchemvar, nsolute,
!$$$     $     ch_var_01, ch_var_02, ch_var_03, ch_var_04,
!$$$     $     ch_var_05, ch_var_06, ch_var_07, ch_var_08,
!$$$     $     ch_var_09, ch_var_10, ch_var_01_m3,
!$$$     $     ch_var_02_m3, ch_var_03_m3, ch_var_04_m3,
!$$$     $     ch_var_05_m3, ch_var_06_m3, ch_var_07_m3,
!$$$     $     ch_var_08_m3, ch_var_09_m3, ch_var_10_m3)
!$$$
!$$$      IF (iresult.NE.0) THEN
!$$$         PRINT *, 'Errors assigning ch_var values:'
!$$$         STOP
!$$$      ENDIF
!
! If chem_sim indicates that no chemical simulations are to be completed
! notify the user.
!
      else
         print*,'Chem_sim parameter equals zero so all chemical ',&
           ' variables are uninitialized'
      end if

      phr_tf=.false. ! Make debug flag false until xdebug_start is reached
      iresult = SetOutputFileOn(ID,phr_tf)
      iresult = SetErrorFileOn(ID,phr_tf)
      iresult = SetLogFileOn(ID,phr_tf)
      iresult = SetSelectedOutputFileOn(ID,phr_tf)

      phreeqmms_init = 0
      return
      end

! ****************************************************************************
!
!     phreeqmms_run - run section
!
 
      integer function phreeqmms_run()

      USE WEBMOD_PHREEQ_MMS
!      USE WEBMOD_INTCP, ONLY : covden_win
      USE WEBMOD_OBSHYD, ONLY : relhum
      USE WEBMOD_OBSCHEM, ONLY : phq_lut, sol_id, sol_name,unit_lut,&
          n_iso, iso_list,c_precip_pH,c_precipT,cconc_precipM, &
          cconc_extM,cconc_obsM,c_ext_pH,c_extT,c_obs_pH,c_obsT

      USE WEBMOD_POTET, ONLY : transp_on
      ! Mixing variables from webmod_res
!      USE WEBMOD_PRECIP, ONLY : mru_ppt
      USE WEBMOD_TEMP1STA, ONLY: tmax_c
      USE WEBMOD_IRRIG, ONLY : irrig_ext_mru, irrig_hyd_mru, &
          irrig_frac_ext, irrig_frac_sat, irrig_frac_hyd, mru_ppt, mru_dep
      USE WEBMOD_TOPMOD, ONLY : gw_loss,qpref_max, st, quz, srzwet, riparian_thresh, uz_area, riparian
      USE WEBMOD_RESMOD, ONLY : vmix_can, vmix_snow, vmix_ohoriz, &
          vmix_uz, vmix_uz2can, vmix_uz2sat, vmix_sat2uz, vmix_uzgen,&
          vmix_uzrip, vmix_uzup, vmix_qdf, vmix_well, vmix_sat, vmix_satpref, &
          vmix_hill, vmix_mru, vmix_hillexp, vmix_stream, vmix_diversion, &
          vmix_chan_loss, vmix_basin, uz2sat_vol,basin_qsim_cm, &
          vmin_canopy
!      double precision vmix_imp(nmru,nresinp), vmix_rz(MAXMNR_3D) ! for later development
      implicit none

      integer, external ::  length

      integer endper
      logical end_run, end_yr, end_mo, end_dy, end_storm, minvol

      integer input_soln
      integer iresult
      integer res_id
      integer k, n, ib, ij, is,nmix ,ndep
      integer ibindx
      integer ib_start, ib_end
      integer chem2var

      real depvol, irrig_frac
      real dt
      
      double precision str_vol, vol_in, vol_out
      double precision chvar_conv_t(nchemvar,nsolute+1)
 
!      double precision cconc_precipM(nsolute)
!      double precision cconc_extM(nchem_ext,nsolute)
!      double precision cconc_obsM(nchemobs,nsolute)

      phreeqmms_run = 1

!/*
! For reactions
!
!      indx_rxn = 98 or destination for final reactions
! get time step for reaction times
      dt = deltim()
      rxnmols = 0.0   ! No dry deposition
      tempc = -1    ! Standard temperature (negative indicates use of temperatures in pqi file.
      ph = 7.0
      tsec = dt*3600   ! dt in hours
      fill_factor = 1.0
!
!
!*/

      nstep=getstep()

!
! Skip run time section if chem_sim=0
!
      if(chem_sim.eq.1) then
!
! Turn phreeq debugging on between xdbug_start and _stop
!
         phr_tf = .false.
         if(xdebug_start.ne.0) then
            if(nstep.le.xdebug_stop.and.nstep.ge.xdebug_start)&
                 phr_tf =.true.
         end if
         iresult = SetOutputFileOn(ID,phr_tf)
         iresult = SetErrorFileOn(ID,phr_tf)
         iresult = SetLogFileOn(ID,phr_tf)
         iresult = SetSelectedOutputFileOn(ID,phr_tf)

            
!
!  Get variables - No need provided by USE module above
!
! These time series are only available if the chemdat file was present
!
!      if(chemdat_exists) then
!         if(ppt_chem.ne.0) then
!            if(getvar('obs_c', 'cconc_precipM', nsolute, 'double', &
!                 cconc_precipM) .ne.0) return
!         end if
!
!         if(chem_ext.gt.0) then
!            if(getvar('obs_c', 'cconc_extM', nchem_ext*nsolute, &
!                 'double',cconc_extM) .ne.0) return
!         end if
!         if(nchemobs.gt.0) then
!            if(getvar('obs_c', 'cconc_obsM', nchemobs*nsolute, &
!                 'double',cconc_obsM) .ne.0) return
!         end if
!      end if


      if(getvar('io', 'endper', 1, 'integer', endper)&
         .ne.0) return

      if(getvar('potet', 'transp_on', nmru, 'integer', transp_on)&
         .ne.0) return

      if(getvar('snow', 'snowcov_area', nmru, 'real', snowcov_area)&
         .ne.0) return

      if(getvar('webr', 'vmix_can', nmru*nresinp, 'double', &
           vmix_can) .ne.0) return

      if(getvar('webr', 'vmix_snow', nmru*nresinp, 'double', &
           vmix_snow) .ne.0) return

      if(getvar('webr', 'vmix_ohoriz', nmru*nresinp, 'double',&
           vmix_ohoriz) .ne.0) return

      if(getvar('webr', 'vmix_uz', nac_nmru_nresinp, 'double', &
           vmix_uz) .ne.0) return

      if(getvar('webr', 'vmix_uz2can', nac*nmru, 'double', &
           vmix_uz2can) .ne.0) return

      if(getvar('webr', 'vmix_uz2sat', nac*nmru, 'double', &
           vmix_uz2sat) .ne.0) return

      if(getvar('webr', 'vmix_sat2uz', nac*nmru, 'double', &
           vmix_sat2uz) .ne.0) return

      if(getvar('webr', 'uz2sat_vol', nac*nmru, 'double', &
           uz2sat_vol) .ne.0) return

      if(getvar('webr', 'vmix_uzgen', nmru*nresinp, 'double', &
           vmix_uzgen) .ne.0) return

      if(getvar('webr', 'vmix_uzrip', nmru*nresinp, 'double', &
           vmix_uzrip) .ne.0) return

      if(getvar('webr', 'vmix_uzup', nmru*nresinp, 'double', &
           vmix_uzup) .ne.0) return

      if(getvar('webr', 'vmix_qdf', nmru*nresinp, 'double', &
           vmix_qdf) .ne.0) return

      if(getvar('webr', 'vmix_sat', nmru*nresinp, 'double', &
           vmix_sat) .ne.0) return

      if(getvar('webr', 'vmix_well', nmru, 'double', &
           vmix_well) .ne.0) return

      if(getvar('webr', 'vmix_satpref', nmru*nresinp, 'double', &
           vmix_satpref) .ne.0) return

      if(getvar('webr', 'vmix_hill', nmru*nresinp, 'double', &
           vmix_hill) .ne.0) return

      if(getvar('webr', 'vmix_mru', nmru*nresinp, 'double', &
           vmix_mru) .ne.0) return

      if(getvar('webr', 'vmix_hillexp', nhydro*nmru, 'double', &
           vmix_hillexp) .ne.0) return

      if(getvar('webr', 'vmix_stream', nhydro, 'double', &
           vmix_stream) .ne.0) return

      if(getvar('webr', 'vmix_diversion', nhydro, 'double', &
           vmix_diversion) .ne.0) return

      if(getvar('webr', 'vmix_chan_loss', nhydro, 'double', &
           vmix_chan_loss) .ne.0) return

      if(getvar('webr', 'vmix_basin', nresinp, 'double', &
           vmix_basin) .ne.0) return


      dt = deltim()

!
! Get date/time stamp for debugging
!
      call dattim('now', datetime)
! Zero basin isotope accumulators
      basin_in_vol=0D0
      basin_out_vol=0D0
      do n=1,nsolute
        ch_basin_in_permil(n)=0D0
        ch_basin_out_permil(n)=0D0
      end do
!
! If not the first time step, transfer final solutions from last step
! to initial solutions for this step. Concentrations for DI (soln 1)
! and observed input sources (2001,2002,...) are either constants
! or read in every time step and are not reset in this section. 
!
      if(.not.step1) then       ! copy ending solutions to initial solutions
         do 31 is = 1,nmru
            do 32 res_id = 1, 9
               src(res_id) = solnnum(1,res_id,0,is,0,0,0)
               dest(res_id) = solnnum(0,res_id,0,is,0,0,0)
 32         continue
            iresult = phr_multicopy(ID,'solution',src, dest, 9)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_multicopy:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     transfer UZ reservoirs
            do 3 ij = 1,nacsc(is)
               src(ij) = solnnum(1,5,0,is,ij,0,0)
               dest(ij) = solnnum(0,5,0,is,ij,0,0)
 3          continue
            iresult = phr_multicopy(id,'solution',src, dest, nacsc(is))
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_multicopy:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!    reset indicator of ionic pulse in snowmelt
            snow_ion_pulse(is)=0
 31      continue               ! end MRU loop

!
! Transfer stream chemistry
!
         do 4 ij = 1, clark_segs
               src(ij) = solnnum(1,99,0,0,0,ij,0)
               dest(ij) = solnnum(0,99,0,0,0,ij,0)
 4       continue
            iresult = phr_multicopy(id,'solution',src, dest, clark_segs)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_multicopy:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!
! Transfer final c_chem reservoirs values to inital values
! for this run. Temp and pH are now included and wull be 
! weighted by input and output volumes.
            do 203 i = 1,nphrsolns
               c_chem(i)%vol(init) = c_chem(i)%vol(fin)
               c_chem(i)%vol(in) = 0.0
               c_chem(i)%vol(out) = 0.0
               c_chem(i)%vol(ET) = 0.0
               c_chem(i)%vol(fin) = 0.0
               c_chem(i)%Temp(init) = c_chem(i)%Temp(fin)
               c_chem(i)%Temp(in) = 0.0
               c_chem(i)%Temp(out) = 0.0
               c_chem(i)%Temp(ET) = 0.0
               c_chem(i)%Temp(fin) = 0.0
               c_chem(i)%pH(init) = c_chem(i)%pH(fin)
               c_chem(i)%pH(in) = 0.0
               c_chem(i)%pH(out) = 0.0
               c_chem(i)%pH(ET) = 0.0
               c_chem(i)%pH(fin) = 0.0
               do 203 j = 1,nsolute
                 c_chem(i)%M(j,init) = c_chem(i)%M(j,fin)
                 c_chem(i)%M(j,in) = 0.0
                 c_chem(i)%M(j,out) = 0.0
                 c_chem(i)%M(j,rxn) = 0.0
                 c_chem(i)%M(j,fin) = 0.0
                 c_chem(i)%delta(j,init) = c_chem(i)%delta(j,fin)
                 c_chem(i)%delta(j,in) = -1000.0
                 c_chem(i)%delta(j,out) = -1000.0
                 c_chem(i)%delta(j,rxn) = -1000.0
                 c_chem(i)%delta(j,fin) = -1000.0
 203        continue
      else
         step1=.false.
      end if
!
! Decompose the endperiod variable
!
      end_run = .false.
      end_yr = .false.
      end_mo = .false.
      end_dy = .false.
      end_storm = .false.
!
! Save a few loops by testing the most common states first.
! Storms need to be tested independently. Others can just be tested
! for the end of the more frequent period. For example, if end_dy is
! not true then, by definition, neither is end_mo, end_yr, and end_run.
!
      if (endper.ne.0.and.mod(endper,2).ne.0) then
         endper = endper - 1
         end_storm = .true.
      end if
      
      if (endper.ne.0) then
         end_dy = .true.
         endper = endper - 2
         if (endper.ne.0) then
            end_mo = .true.
            endper = endper - 4
            if (endper.ne.0) then
               end_yr = .true.
               endper = endper - 8
               if (endper.ne.0) end_run = .true.
            end if
         end if
      end if

!      Initialize reaction parameters to 'mixing only'

!      do k=1,11
!         n_user(k) = -1
!      end do
!      rxnmols = 0.0   ! No default dry deposition
!      tempc = 25.0    ! Standard temperature
!      ph = 7          ! Dummy input designed to flag missing real input
!      tsec = dt*3600  ! convert dt, in hours, to seconds for phreeqc kinetic rx

!
! If no chemdat file was found use the initial solutions
! defined in the init section throughout the model run.
!
! If a chemdat file exists, create the solutions indicated
! by cconc_precipM(nsolute), cconc_extM(nchem_ext,nsolute),
! and cconc_obsM(nchemobs,nsolute) if ppt_chem=1, chem_ext=1,
! or nobs_chem>0, respectively.
!
! Also allow for varying chemical concentrations of irrigation and/or
! QW samples even if precip chemistry is to remain constant (ppt_chem=0).
! If precipitation chemistry is described in the chemdat file then
! ppt_chem = 1. Similarly if the composition of external irrigation
! is described in the chemdat file then chem_ext = 1
!
! nchemdat = nchemdep + nchemobs, 
! where
! nchemdep = 1(ppt chemistry either constant or varying)+nchem_ext
!
! nchemdat_obs equals the number of time varying chemical sources
! ppt_chem + chem_ext*nchem_ext + nchemobs
!
      nchemdat_obs = ppt_chem + chem_ext*nchem_ext + nchemobs

      if(chemdat_exists) then  ! create as many solutions as are in the chemdat file. Moved precip solution to MRU loop to assign dry bulb temperature.
         ib_start = 2-ppt_chem
         ib_end   = nchemdat_obs-ppt_chem+1
         do 482 ib = ib_start, ib_end
            ibindx = ib
            if(ib.gt.1.and.chem_ext.eq.0) & ! sample obs but no variable irrigation solutions
                 ibindx = ib+nchem_ext
            input_soln = solnnum(0,0,ibindx,0,0,0,0)
! pick up tempc and pH from data file
            if(ib.eq.1) then
              tempc = c_precipT
              if(tempc.eq.-999) tempc = 25  ! -999 is precip temperature flag to assign wetbulb temperature at elevation of MRU. Assign 25 to initial mix.
              ph = c_precip_pH
            else if (ib-1.le.nchem_ext*chem_ext) then
              tempc = c_extT(ib-1)
              ph = c_ext_pH(ib-1)
            else
              tempc = c_obsT(ibindx-nchemdep)
              ph = c_obs_pH(ibindx-nchemdep)
            end if
! pick up solute concentrations
         do 5 j = 1, nsolute
            if(ib.eq.1) then
               conc(j) = cconc_precipM(j)
            else if (ib-1.le.nchem_ext*chem_ext) then
               conc(j) = cconc_extM(ib-1,j)
            else
               conc(j) = cconc_obsM(ibindx-nchemdep,j)
            end if
 5       continue
! Create solution 
         iresult = phr_precip(ID,input_soln, nsolute, sol_name, conc, tempc, ph)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_precip:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
         if(ib.gt.nchemdep) then
! if there are observations (point samples) create a one liter solution with
! the same mass in inputs and outputs. precip and irrigation solutions
! will be created in the mru loop
            totvol=1e-3  ! one liter
            indx = isoln(input_soln,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,totvol,0,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         end if
 482    continue
      end if

! Begin MRU loop *********************************************

      do 10 is = 1, nmru
! Get the solution indices for tracking MRU and UZfluxes
!
         iphrq_mru = is
         mru_in_vol(is)=0D0
         mru_out_vol(is)=0D0
         delta_D = -1000D0
         delta_18O = -1000D0
         do n=1,nsolute
             ch_mru_in_permil(is,n)=0D0
             ch_mru_out_permil(is,n)=0D0
         end do
         indxm = chmru_soln(is)
         ndep = 0
         depvol = mru_dep(is)*mru_area(is)*inch2m*a_million
!
! Create input solution if precip or irrigation occured during this
! time step. Note that mru_dep is a combination of atmospheric precip,
! and internal and external irrigation sources.
!
         if(mru_dep(is).gt.0.0) then
            nmix = 1
! post input chemistry to c_chem matrix
! need to mix to same to retrieve concentration
            irrig_frac = irrig_frac_ext(is) + irrig_frac_sat(is)+&
                 irrig_frac_hyd(is)
!
! Atmospheric precipitation inputs **************************
!
            if(irrig_frac.lt.1.0) then !atmospheric dep > 0
               ndep = ndep + 1
               src(1) = solnnum(0,0,1,0,0,0,0)
               srcdep(ndep)=src(1)
               dest(1) = src(1)
               totvol=mru_dep(is)*(1-irrig_frac)*&
                    mru_area(is)*inch2m*a_million
               fracs(1) = 1.0
               fracsdep(ndep)=totvol/depvol


! Check for negative fractions

               iresult=checkfracs(nmix,src,fracs,dest(1))
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Assign wet bulb temperature to precip if indicated by temperature = -999 in chemdat file or if constant composition (ppt_chem=0)
               if(c_precipT.eq.-999.or.ppt_chem.eq.0) then
                   WRITE (line,100),'SOLUTION_MODIFY ', src(1)
                   iresult = AccumulateLine(id, line)
                   tempc=wetbulb(kPa(is),tmax_c(is), relhum(1)*100)
                   if(tempc.lt.0.) tempc=0.
                   WRITE (line,120),'-temp ', tempc
                   iresult = AccumulateLine(id, line)
                   WRITE (line,110),'END'
                   iresult = AccumulateLine(id, line)
                   iresult = RunAccumulated(id)
               end if
! Mix with no reactions (no_rxn instead of n_user)
               iresult = phr_mix(ID,nmix, src, fracs,  dest(1),&
                    fill_factor, dest(1), conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)


               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!
! Place model inputs with the input metric '2' such that the 5th index,
! the final volume or mass, works as an accumulator for the entire run
! Note that inputs from atmospheric precip is always the second index
! in c_chem so no need to use the isoln function here.
!
               iresult = update_chem(2,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
!
! external irrigation inputs ***********************************
!
            if(irrig_ext_mru(is).gt.0) then
               ndep = ndep + 1
               ib = src_ext_irrig(is)+1 ! first index (ib=1) is average basin precip concentration
               src(1)  = solnnum(0,0,ib,0,0,0,0)
               srcdep(ndep)=src(1)
               dest(1) = src(1)
               totvol=mru_dep(is)*irrig_frac_ext(is)*&
                    mru_area(is)*inch2m*a_million
               fracs(1) = 1.0
               fracsdep(ndep)=totvol/depvol

!               iresult = fill_ent(n_user,dest(1),nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if


! Check for negative fractions

               iresult=checkfracs(nmix,src,fracs,dest(1))
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,nmix, src, fracs,  dest(1),&
                    fill_factor, dest(1), conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!
! update c_chem irrigation inputs and add irrigation to inputs     
!
               indx = isoln(dest(1),nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
!
! irrigation from a well in the MRU ************************
!
            if(irrig_frac_sat(is).gt.0) then
               ndep = ndep + 1
               src(1)  = solnnum(0,8,0,is,0,0,0)
               srcdep(ndep)=src(1)
               dest(1) = src(1)
               totvol=mru_dep(is)*irrig_frac_sat(is)*&
                    mru_area(is)*inch2m*a_million
               fracs(1) = 1.0
               fracsdep(ndep)=totvol/depvol
!
!               iresult = fill_ent(n_user,dest(1),nchemdat,nmru,
!     $              nac,clark_segs,src_init)
               if(iresult.ne.0) then
                  PRINT *, 'Errors assigning phreeq entities:'
                  CALL OutputErrorString(id)
                  STOP
               end if


! Check for negative fractions

               iresult=checkfracs(nmix,src,fracs,dest(1))
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,nmix, src, fracs,  dest(1),&
                    fill_factor, dest(1), conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!
! update c_chem with irrigation outputs from the saturated zone
!
               indx = isoln(dest(1),nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires) 
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
               
!     add outputs from wells as outputs from MRU and basin pseudo solutions
               basin_out_vol=basin_out_vol+totvol
               mru_out_vol(is)=mru_out_vol(is)+totvol
               c_chem(indxm)%vol(out)= &
                          c_chem(indxm)%vol(out)+totvol
               c_chem(indxb)%vol(out)=&
                          c_chem(indxb)%vol(out)+totvol
               do k=1,nsolute
                 c_chem(indxm)%M(k,out)=&
                  c_chem(indxm)%M(k,out)+ c_chem(indx)%M(k,out)
                 c_chem(indxb)%M(k,out)=&
                  c_chem(indxb)%M(k,out)+ c_chem(indx)%M(k,out)
                 ch_mru_out_permil(is,k)=ch_mru_out_permil(is,k)+&
                       c_chem(indxm)%delta(k,out)*totvol
                 ch_basin_out_permil(k)=ch_basin_out_permil(k)+&
                       c_chem(indxb)%delta(k,out)*totvol
               end do
            end if
!
! irrigation from a stream diversion ***************************
!
            if(irrig_frac_hyd(is).gt.0) then
               ib = irrig_int_src(is)
               ndep = ndep + 1
               src(1)  = solnnum(0,99,0,0,0,ib,0)
               srcdep(ndep)=src(1)
               dest(1) = src(1)
               totvol=mru_dep(is)*irrig_frac_hyd(is)*&
                    mru_area(is)*inch2m*a_million
               fracs(1) = 1.0
               fracsdep(ndep)=totvol/depvol

!               iresult = fill_ent(n_user,dest(1),nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if


! Check for negative fractions

               iresult=checkfracs(nmix,src,fracs,dest(1))
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,nmix, src, fracs,  dest(1),&
                    fill_factor, dest(1), conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!
! update c_chem irrigation outputs from stream segments
!
               indx = isoln(dest(1),nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,indxb,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
!
! Mix precip and irrigation to obtain hillslope inputs
!
            dest(1) = solnnum(0,11,0,is,0,0,0)
!
! Check for negative fractions
!
            iresult=checkfracs(ndep,srcdep,fracsdep,dest(1))
            if(iresult.ne.0) then
              PRINT*,'Errors with mixing fractions'
            end if
! Mix with no_rxn
            iresult = phr_mix(ID,ndep, srcdep, fracsdep,  dest(1),&
                 fill_factor, dest(1), conc, phr_tf, no_rxn,rxnmols,&
                 tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF

!     Update volume/mole matrix with MRU combined inputs. As with
!     the individual atmospheric and irrigation inputs use
!     metric 2 (inputs) so that the final volumes and masses in
!     the c_chem matrix serve as accumulators.
!
!            totvol = mru_dep(is)
            totvol = depvol
            indx = isoln(dest(1),nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,indxm,indxb,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!
! Update isotope data
!
            basin_in_vol=basin_in_vol+totvol
            mru_in_vol(is)=mru_in_vol(is)+totvol
            do k=1,nsolute
              ch_mru_in_permil(is,k)=ch_mru_in_permil(is,k)+&
                c_chem(indxm)%delta(k,in)*totvol
              ch_basin_in_permil(k)=ch_basin_in_permil(k)+&
                c_chem(indxb)%delta(k,in)*totvol
            end do
        end if   ! mru_dep > 0
!
! Mix groundwater inputs (without reactions) into temporary reservoir
! and add to MRU and basin pseudo solutions. Actual inputs into groundwater
! will be tracked in the groundwater section below
!
        totvol = vmix_sat(is,20) + vmix_sat(is,21)
        if(totvol.gt.0.0) then
          dest(1) = solnnum(0,14,0,is,0,0,0)
          ib = src_gw1(is)+1
          src(1) = solnnum(0,0,ib,0,0,0,0) ! gw input from 1st source (channel leakage?)
          fracs(1) = vmix_sat(is,20)/totvol
          ib = src_gw2(is)+1
          src(2) = solnnum(0,0,ib,0,0,0,0) ! gw input from 2nd source (upgradient?)
          fracs(2) = vmix_sat(is,21)/totvol

!            iresult = fill_ent(n_user,dest(1),nchemdat,nmru,
!     $           nac,clark_segs,src_init)
!            if(iresult.ne.0) then
!               PRINT *, 'Errors assigning phreeq entities:'
!               CALL OutputErrorString(id)
!               STOP
!            end if


! Check for negative fractions

          iresult=checkfracs(2,src,fracs,dest(1))
          if(iresult.ne.0) then
             PRINT*,'Errors with mixing fractions'
          end if
! Mix with no_rxn
          iresult = phr_mix(ID,2, src, fracs,  dest(1),&
               fill_factor, dest(1), conc, phr_tf, no_rxn,rxnmols,&
               tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

          IF (iresult.NE.0) THEN
           PRINT *, 'Errors during phr_mix:'
           CALL OutputErrorString(id)
           STOP
          ENDIF
!
! Add inputs to MRU and basin pseudo solutions. Inputs to actual groundwater reservoir
! computed in groundwater section below
!
          basin_in_vol=basin_in_vol+totvol
          mru_in_vol(is)=mru_in_vol(is)+totvol
          c_chem(indxm)%vol(in)= &
            c_chem(indxm)%vol(in)+totvol
          c_chem(indxb)%vol(in)= &
            c_chem(indxb)%vol(in)+totvol
          do k=1,nsolute
            c_chem(indxm)%M(k,in)=&
              c_chem(indxm)%M(k,in)+ conc(k)*totvol*a_thousand   ! conc in Moles/L
            c_chem(indxb)%M(k,in)=&
              c_chem(indxb)%M(k,in)+ conc(k)*totvol*a_thousand
            ch_mru_in_permil(is,k)=ch_mru_in_permil(is,k)+&
                  c_chem(indxm)%delta(k,in)*totvol
            ch_basin_in_permil(k)=ch_basin_in_permil(k)+&
                  c_chem(indxb)%delta(k,in)*totvol
          end do
        end if ! totvol.gt.0
!ccccccccccc end inputs
! Begin Canopy
!
!        Record ET volume to canopy and make initial and final canopy volumes
!        available for reaction accounting
         mixture = solnnum(1,1,0,is,0,0,0) ! t+1 canopy use for throughfall
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_can(is,1)
               c_chem(indx)%vol(ET)=vmix_can(is,6)
               c_chem(indx)%vol(fin)=vmix_can(is,4)

! The following section executes on days with no transpiration
! (vmix_can(is,11) > 0), but when there has been some change in the
! interception chemistry because there was
!     throughfall (vmix_can(is,3) > 0),
!     precip (vmix_can(is,5) > 0), or
!     ET (vmix_can(is,6) > 0).
!
!

         if(vmix_can(is,11).eq.0.0.and.(vmix_can(is,3).ne.0.0&
           .or.vmix_can(is,5).ne.0.0.or.vmix_can(is,6).ne.0.0)) then
!
! Canopy scenario 1: ET from precipitation or canopy storage during
! this time step so no transpiration. Canopy scenario 2 describes days
! when canopy is dry so solutes are transpired from the unsaturated
! zone to the canopy. The code for scenario 2 including transpiration
! and canopy evaporation is included further on in this section
! with the rest of the UZ fluxes and transport.
!
!     Canopy inputs (use temporary input reservoir number 14)

            mixture = solnnum(0,14,0,is,0,0,0)

!     solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)
            totvol = vmix_can(is,2)

            if (totvol.gt.0.0) then ! Rain and/or leaves on
               
               solns(1) = solnnum(0,11,0,is,0,0,0)
               fracs(1) = (vmix_can(is,5)+vmix_can(is,17)+&
                  vmix_can(is,18)+vmix_can(is,19))/totvol
               solns(2) = solnnum(0,4,0,is,0,0,0) ! o-horiz on days of leaves on
               fracs(2) = vmix_can(is,10)/totvol

!               iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if


! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,2, solns, fracs,  mixture,&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)
               

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
            endif

! Recognize canopy reservoir, update c_chem, and prepare final solution
            mixture = solnnum(1,1,0,is,0,0,0) ! t+1 canopy after evaporation and throughfall

            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
! basin and mru inputs included in combined input section above
            iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     
! Mix intial and inputs to determine initial composition and volume to compute
! fractionation parameters if D and/or 18O are being tracked
!
            totvol = vmix_can(is,1)+vmix_can(is,2)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            if (totvol.gt.0.0) then
               solns(1) = solnnum(0,1,0,is,0,0,0) ! initial canopy chem
               fracs(1) = vmix_can(is,1)/totvol
               
               solns(2) = solnnum(0,14,0,is,0,0,0) ! canopy inputs
               fracs(2) = vmix_can(is,2)/totvol
               

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
!
! Mix initial and inputs storing the conservative mix back into the solution number for the initial chem (solns(1)
!
               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)


               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
! --------------  iso evap section
!
! If D or 18O are included in list of solutes of interest, compute fractionation of evaporated water.
! Evaporated D and 18O are tracked as basin exports in the function call. 
! if n_iso, the number of isotopes, =0, this section is skipped and DI is evaporated
               if(n_iso.ne.0.and.vmix_can(is,6).gt.0) then
                 evap=vmix_can(is,6)
                 ison=iso_n(is)
                 rh=relhum(1)
                 totvol_can = totvol - vmin_canopy(is)
! Rayleigh fractionation of evaporated water. Reduce volume 
                 iresult=fractionate(water,indx,evap,totvol_can,ison,rh,tempc)
                 IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during evaporation fractionation'
                  CALL OutputErrorString(id)
                  STOP
                 ENDIF
               end if
! --------------  //iso evap section

! remove evaporated water from canopy
               totvol = totvol-vmix_can(is,6)
               solns(1) = solnnum(0,1,0,is,0,0,0) ! mixture of initial and canopy inputs created in previous mix
               fracs(1) = (vmix_can(is,1)+vmix_can(is,2))/totvol
               solns(2) = solnnum(0,0,0,0,0,0,0)    ! soln 1: ET is DI water
               fracs(2) = -vmix_can(is,6)/totvol
! fill_ent to allow reactions in this final mix before export (uses final volume)
               iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
               if(iresult.ne.0) then
                 PRINT *, 'Errors assigning phreeq entities:'
                 CALL OutputErrorString(id)
                 STOP
               end if
! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                 fill_factor, mixture, conc, phr_tf, n_user,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
               IF (iresult.NE.0) THEN
                 PRINT *, 'Errors during phr_mix:'
                 CALL OutputErrorString(id)
                 STOP
               ENDIF

            endif
!     Update volume/mole matrix with output volumes
            totvol = vmix_can(is,3)
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                                  tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
! reset DI to pure water (no 18O or D)
            iresult = accumulateline(id, "solution 1")  ! solution 1 is always evaporative water
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors establishing DI solution 1'
               STOP
            ENDIF
            iresult = runaccumulated(ID)   ! assigns depleted deltas to evap using lines accumulated above         
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors establishing DI solution 1'
               STOP
            ENDIF

!     
! End of canopy mixes on days with rain or ET from intcp storage.
! When indicated, transpiration and canopy ET (canopy scenario 2)
! will take place in the UZ section below.
!
! The situation when there is transpiration on days of leaves off
! is not permitted in webmod_res (the transpiration rate is set 
! to zero. This avoids the necessity of waiting until the end of
! UZ zone mixing to deliver water to the canopy that will then be
! deposited on the surface as throughfall.
!
         else
            if(vmix_can(is,11).eq.0) then ! no canopy change. Allow reactions
               fracs(1)=1.0
               solns(1) = solnnum(0,1,0,is,0,0,0)
               mixture = solnnum(1,1,0,is,0,0,0)

               iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                    nac,clark_segs,src_init)
               if(iresult.ne.0) then
                  PRINT *, 'Errors assigning phreeq entities:'
                  CALL OutputErrorString(id)
                  STOP
               end if

               iresult = phr_mix(ID,1, solns, fracs,&
                    solns(1),fill_factor, mixture, conc, phr_tf,&
                    n_user,rxnmols,tempc,ph,ph_final,tsec,&
                    tally_table,ntally_rows,ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF

! update chem matrix - no volume needed if imetric = -1 since
! the final volume (same as init) is assessed in update_chem.
               indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,-1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
         end if
!
! Snowpack section
!
! The following section executes on days when there is
!     new snow,rain, or throughfall (vmix_snow(is,2) > 0), or
!     sublimation (vmix_snow(is,6) >0), or
!     melt (vmix_snow(is,3) > 0).
! Skip the section if none if these is true. This should work
! for days of no snow and also when the pack is unchanged
!
!        Record sublimation and make initial and final canopy volumes
!        available for reaction accounting

         mixture = solnnum(1,2,0,is,0,0,0) ! t+1 snowpack
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_snow(is,1)
               c_chem(indx)%vol(ET)=vmix_snow(is,6)
               c_chem(indx)%vol(fin)=vmix_snow(is,4)

         if(vmix_snow(is,2).gt.0.0.or.vmix_snow(is,6).gt.0.0 .or.&
          vmix_snow(is,3).gt.0.0) then

! Snowpack inputs (use temporary input reservoir number 14)

            mixture = solnnum(0,14,0,is,0,0,0)

!     solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)
            totvol = vmix_snow(is,2)

            if (totvol.gt.0.0) then ! Rain and/or throughfall
               
               solns(1) = solnnum(0,11,0,is,0,0,0) ! precip+irrig inputs
               fracs(1) = (vmix_snow(is,5)+vmix_snow(is,17)+&
                  vmix_snow(is,18)+vmix_snow(is,19))/totvol
               
               solns(2) = solnnum(1,1,0,is,0,0,0) ! throughfall computed above
               fracs(2) = vmix_snow(is,8)/totvol


! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,2, solns, fracs,  mixture,&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!     Update volume/mole matrix for snowpack inputs
               mixture = solnnum(1,2,0,is,0,0,0)

               indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
     
    
               iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            endif
!     
!     Mix and concentrate (if sublimation) solutes
!
            mixture = solnnum(1,2,0,is,0,0,0) ! t+1 snowpack solution
            totvol = vmix_snow(is,1)+vmix_snow(is,2)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            if (totvol.gt.0.0) then
               solns(1) = solnnum(0,2,0,is,0,0,0) ! initial snowpack chem
               fracs(1) = vmix_snow(is,1)/totvol
               
               solns(2) = solnnum(0,14,0,is,0,0,0) ! snowpack inputs
               fracs(2) = vmix_snow(is,2)/totvol

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
! --------------  iso evap section
!
! If D or 18O are included in list of solutes of interest, compute fractionation of sublimated water.
! Sublimated D and 18O are tracked as basin exports in the fractionate routine. 
! if n_iso, the number of isotopes, =0, this section is skipped and DI is evaporated
               if(n_iso.ne.0.and.vmix_snow(is,6).gt.0) then
                 evap=vmix_snow(is,6)
                 ison=iso_n(is)
                 rh=relhum(1)
!                 iresult = test(conc)
                 iresult=fractionate(snow,indx,evap,totvol,ison,rh,tempc)  ! Rayleigh fractionation of evapaporated (or sublimated) water
!                 iresult = test(conc)
                 IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during sublimation fractionation'
                  CALL OutputErrorString(id)
                  STOP
                 ENDIF
               elseif(n_iso.ne.0) then  ! No sublimation but isotopes tracked. Assign deltas to reservoirs for use in pulse computations
                 res_D_permil = -1000D0
                 res_18O_permil = -1000D0
                 i = 1
                 do while(i.le.n_iso)
                   id_len = length(sol_name(iso_list(i)))
                   if(sol_name(iso_list(i))(1:id_len).eq."D") then
                     res_D_permil = conc(nsolute+3*i-2) ! delta D of reservoir water 
                   elseif(sol_name(iso_list(i))(1:id_len).eq."[18O]") then
                     res_18O_permil = conc(nsolute+3*i-2)! delta 18O of reservoir water 
                   else
                     print *,"fractionated isotope not D or [18O]. Run stopped"
                     return
                   endif              
                   i=i+1
                 end do
               end if
! --------------  //iso evap section

! remove sublimated water from snow pack
               totvol = totvol-vmix_snow(is,6)
!               iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
!                               tally_table,n_ent,0,0,ires)
!               IF (iresult.NE.0) THEN
!                 PRINT *, 'Errors updating mole matrix:'
!                 STOP
!               ENDIF
               
               solns(1) = solnnum(1,2,0,is,0,0,0) ! mixture of initial and snowpack inputs created in previous mix
               fracs(1) = (vmix_snow(is,1)+vmix_snow(is,2))/totvol
               solns(2) = solnnum(0,0,0,0,0,0,0)    ! soln 1:  DI water
               fracs(2) = -vmix_snow(is,6)/totvol
               iresult = phr_mix(ID,2,solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
            endif
! Melt computations
!
!              if(snow_ion_factor(is).gt.1.0) then
!                mixture = solnnum(0,2,0,is,0,0,0) ! t0 snowpack solution for snowpack, use to store melt that has ionic pulse 
!              else
!                mixture = solnnum(1,2,0,is,0,0,0) ! t1 snowpack solution for snowpack and melt chemistry since no ion pulse simulated
!              endif
            if(vmix_snow(is,3).gt.0.0)then  ! Snow has melted
              totvol = vmix_snow(is,3)
              percent_melt = vmix_snow(is,3)/(vmix_snow(is,1)+ &
                             vmix_snow(is,2)-vmix_snow(is,6))
!              indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
!                  ires,ichemdat,imru,inac,ihydro)
! no ionic pulse because snow_ion_factor set to 1, or most (>90% of pack) if not all snow has melted.
              if(snow_ion_factor(is).le.1.0.or.vmix_snow(is,4).eq.zero.or.percent_melt.ge.0.9) then

! Snowmelt with no ion pulse being simulated. Track melt with same composition as
! pack. Remaining Pack composition set after possible ionic pulse block.
                   iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                                      tally_table,n_ent,0,0,ires)
                   IF (iresult.NE.0) THEN
                      PRINT *, 'Errors updating mole matrix:'
                      STOP
                   ENDIF
                   !iresult = update_chem(indx,totvol,-1,conc,pH,tempc,& ! totvol ignored as final volume used in update_chem
                   !                   tally_table,n_ent,0,0,ires)
                   !IF (iresult.NE.0) THEN
                   !   PRINT *, 'Errors updating mole matrix:'
                   !   STOP
                   !ENDIF

!
! Simulate ionic pulse and isotopic fractionation of melt if snow_ion_factor > 1. Throttle ion concentration factor (ICF) to 0.9*max_factor.
!
              else ! ionic pulse simulated
!                if(snow_ion_factor(is).gt.1.0.and.max_factor.gt.1.0) then ! no pulse if percent melt is so high that max factor > 1
                 mixture = solnnum(0,2,0,is,0,0,0) ! t0 snowpack solution for snowpack, use to store melt that has ionic pulse 
                 snow_ion_pulse(is)=1
                 ion_factor = snow_ion_factor(is)
                 !max_factor = 0.9/percent_melt
                 !if(ion_factor.gt.max_factor) ion_factor = max_factor
               
! Assign proper delta to concentrating DI to produce depletion in melt as reported in Taylor and others 2001.
! Parameters snowmelt_18O_depl and snowmelt_D_depl indicate the desired permil depletion of the heavier isotopes in the melt
! versus pack on each melt day. The DI delta should be (delta_pack + snowmelt_*_depl/(1-ICF)).
!
                 if(n_iso.ne.0) then   ! assign the proper delta to the DI, delta values for the pack 
                                       ! were stored in the fractionate subroutine above
                   ! int ipack, int imelt, double eps, double ipf, double fmelt, double rstd    
                   CALL meltpack(id,solns(1),mixture, dble(snowmelt_18O_depl(is)), dble(ion_factor), dble(percent_melt),phq_lut(sol_id(iso_list(1))%phq)%isoratio )
                   iresult = accumulateline(ID, "SELECTED_OUTPUT")
                   iresult = accumulateline(ID, "-reset false")
                   iresult = accumulateline(ID, "-pH")
                   iresult = accumulateline(ID, "-temp")
                   iresult = accumulateline(ID, "USER_PUNCH")
                   write(aline,1000)iso_header1(1:isoh1_len)
                   iresult = AccumulateLine(id, aline)
                   write(aline,1050)iso_header2(1:isoh2_len)
                   iresult = AccumulateLine(id, aline)
                   iresult = runaccumulated(ID)   ! assigns depleted deltas and temperature to evap using lines accumulated above
!                   if(snowmelt_D_depl(is).le.0.0.and.res_D_permil.ne.-1000) &
!                       delta_D = res_D_permil+snowmelt_D_depl(is)/(1-ion_factor)
!!                   if(snowmelt_18O_depl(is).lt.0.0.and.snowpack_18O.ne.-1000) &
!                   if(snowmelt_18O_depl(is).le.0.0.and.res_18O_permil.ne.-1000) &
!                       delta_18O = res_18O_permil+snowmelt_18O_depl(is)/(1-ion_factor)
!                   iresult = reset_DI(ID,zero)  ! This will create a solution '2' at zero degrees celsius to concentrate melt and produce the right delta 18O
!                   IF (iresult.NE.0) THEN
!                     PRINT *, 'Errors establishing deltas for DI solution 2'
!                     STOP
!                   ENDIF
                 end if
!                   solns(1) = solnnum(1,2,0,is,0,0,0) ! mixture of initial and snowpack inputs created in previous mix
!                   fracs(1) = ion_factor
!!                   solns(2) = solnnum(0,0,0,0,0,0,0)    ! soln 1:  DI water
!                   solns(2) = 2    ! soln 1:  concentrating solution
!                   fracs(2) = 1.0-ion_factor
!! Check for negative fractions
!                   iresult=checkfracs(2,solns,fracs,mixture)
!                   if(iresult.ne.0) then
!                      PRINT*,'Errors with mixing fractions'
!                   end if
! concentrated solutes in melt
                   !iresult = phr_mix(ID,2,solns, fracs,  mixture,&
                   !     fill_factor, mixture, conc, phr_tf, no_rxn,&
                   !     rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                   !     ntally_cols)
                   solns(1) = solnnum(0,2,0,is,0,0,0)
                   fracs(1) = 1.0
                   iresult = phr_mix(ID,1,solns, fracs,  mixture,&
                        fill_factor, mixture, conc, phr_tf, no_rxn,&
                        rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                        ntally_cols)
                   IF (iresult.NE.0) THEN
                      PRINT *, 'Errors during phr_mix:'
                      CALL OutputErrorString(id)
                      STOP
                   ENDIF
                   !     Update volume/mole matrix with snowmelt volumes and masses
                   iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                                      tally_table,n_ent,0,0,ires)
                   IF (iresult.NE.0) THEN
                      PRINT *, 'Errors updating mole matrix:'
                      STOP
                   ENDIF
!            endif


              endif  ! ionic pulse or not
            endif   ! did melt occur
! reconstitute diluted snowpack
            mixture = solnnum(1,2,0,is,0,0,0) 
!                   solns(1) = solnnum(1,2,0,is,0,0,0) ! mixture of initial and snowpack inputs created in previous mix
!                   fracs(1) = 1.0/(1-percent_melt)
!                   solns(2) = solnnum(0,2,0,is,0,0,0)    ! soln 2:  melt
!                   fracs(2) = 1-fracs(1)
!! Check for negative fractions
!                   iresult=checkfracs(2,solns,fracs,mixture)
!                   if(iresult.ne.0) then
!                      PRINT*,'Errors with mixing fractions'
!                   end if
! fill_ent to allow reactions in this final mix before export (uses final volume)
            iresult = fill_ent(n_user,mixture,nchemdat,nmru,nac,clark_segs,src_init)
            if(iresult.ne.0) then
              PRINT *, 'Errors assigning phreeq entities:'
              CALL OutputErrorString(id)
              STOP
            end if
! remove solutes from pack with mix
                   !iresult = phr_mix(ID,2,solns, fracs,  solns(1),&
                   !     fill_factor, mixture, conc, phr_tf, no_rxn,&
                   !     rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                   !     ntally_cols)
            solns(1) = solnnum(1,2,0,is,0,0,0)
            fracs(1) = 1.0
            iresult = phr_mix(ID,1,solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
            iresult = update_chem(indx,totvol,-1,conc,pH,tempc,&  ! totvol ignored as final volume is accessed in update_chem
                   tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
          else if(vmix_snow(is,4).gt.0) then   ! No change in snowpack so mix original to final soln so
                                              ! that reactions can take place. Skip this if snowpack is melted
            fracs(1)=1.0
            solns(1) = solnnum(0,2,0,is,0,0,0)
            mixture = solnnum(1,2,0,is,0,0,0)

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

            iresult = phr_mix(ID,1, solns, fracs,&
                 solns(1),fill_factor, mixture, conc, phr_tf,&
                 n_user,rxnmols,tempc,ph,ph_final,tsec,&
                 tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
! update chem matrix - no volume needed if imetric = -1 since
! the final volume (same as init) is assessed in update_chem.
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,totvol,-1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         end if
!     
! End of snowpack fluxes
!

!
! O-horizon section. Rinses with overland flow. Need to include
! DOC generation and other reactions.

! The following section executes on days when there is
!     overland flow (vmix_ohoriz(is,2)>0
!       from direct precip (vmix_ohoriz(is,5) > 0),
!       from snowmelt  (vmix_ohoriz(is,9) > 0), or
!       from residual canopy volume on day of leaves-off day
!         (vmix_ohoriz(is,8) > 0).
!     or to account for return of residual canopy volume
!         on days of leaves-on (vmix_can(is,10) > 0)
! Skip the section if none if these is true.
!
! Modify to include reaction chemistry
!
         mixture = solnnum(1,4,0,is,0,0,0) ! t+1 O-horizon
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_ohoriz(is,1)
               c_chem(indx)%vol(ET)=vmix_ohoriz(is,6)   ! This should be zero always
               c_chem(indx)%vol(fin)=vmix_ohoriz(is,4)

         if(vmix_ohoriz(is,2).gt.0.0.or.vmix_can(is,10).gt.0.0) then

!     O-horizon inputs (use temporary input reservoir number 14)

            mixture = solnnum(0,14,0,is,0,0,0)

!     solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)
            totvol = vmix_ohoriz(is,2)

            if (totvol.gt.0.0) then ! overland flow or leaves off
               
               solns(1) = solnnum(0,11,0,is,0,0,0) ! precip+irrig inputs
               fracs(1) = (vmix_ohoriz(is,5)+vmix_ohoriz(is,17)+&
                  vmix_ohoriz(is,18)+vmix_ohoriz(is,19))/totvol
               
               solns(2) = solnnum(1,1,0,is,0,0,0) ! throughfall at t+1
               fracs(2) = vmix_ohoriz(is,8)/totvol

               if(snow_ion_pulse(is).eq.1) then
                 solns(3) = solnnum(0,2,0,is,0,0,0) ! ionic snowmelt stored in t0 solution
               else
                 solns(3) = solnnum(1,2,0,is,0,0,0) ! no ionic pulse, use snowpack chemistry
               end if
               fracs(3) = vmix_ohoriz(is,9)/totvol

!               iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if

! Check for negative fractions

               iresult=checkfracs(3,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,3, solns, fracs,  mixture,&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!     Update volume/mole matrix for O-horizon inputs
               indx = isoln(solnnum(0,4,0,is,0,0,0),nchemdat,nmru,&
                    nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            endif
!     
!     The effluent to the canopy on days of leaves-on needs to be
!     tracked using the t0 concentrations. Use the t+1 solution
!     temporarily to return initial concentrations. This solution will
!     be overwritten in the next section
!     
            solns(1) = solnnum(0,4,0,is,0,0,0) ! t0 O-horizon solution
            mixture = solnnum(1,4,0,is,0,0,0) ! Temp t+1 O-horizon solution
            fracs(1) = 1.0
            totvol  = vmix_can(is,10) ! output to canopy
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)

!            iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $           nac,clark_segs,src_init)
!            if(iresult.ne.0) then
!               PRINT *, 'Errors assigning phreeq entities:'
!               CALL OutputErrorString(id)
!               STOP
!            end if

            iresult = phr_mix(ID,1, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
            
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     No ET volumes for this layer
!     Update volume/mole matrix with  volumes and masses
            totvol = vmix_can(is,10)
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     
!     Create t+1 O-horizon solution for export to transient hill
!     reservoir. Note that the initial volume is reduced by the 
!     volume that was passed to the canopy on days of leaves on.
!     This should avoid any double accounting when it rains on
!     days of leaves on.
!     
            mixture = solnnum(1,4,0,is,0,0,0) ! t+1 O-horizon solution
            totvol = vmix_ohoriz(is,1)+vmix_ohoriz(is,2)-&
                 vmix_can(is,10)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            if (totvol.gt.0.0) then
               solns(1) = solnnum(0,4,0,is,0,0,0) ! initial O-horiz chem
               fracs(1) = (vmix_ohoriz(is,1)-vmix_can(is,10))/totvol
               
               solns(2) = solnnum(0,14,0,is,0,0,0) ! O-horiz inputs
               fracs(2) = vmix_ohoriz(is,2)/totvol

! Check for negative fractions
               iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                    nac,clark_segs,src_init)
               if(iresult.ne.0) then
                  PRINT *, 'Errors assigning phreeq entities:'
                  CALL OutputErrorString(id)
                  STOP
               end if

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, n_user,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
            end if
!     Update volume/mole matrix with  volumes and masses
            totvol = vmix_ohoriz(is,3)-vmix_can(is,10)
!     Append outputs to hill reservoir to the outputs written 
!     the day of leaves-on
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     
!     End of O-horizon fluxes
!     
         else                   ! no change in o-horizon so allow reactions
            solns(1) = solnnum(0,4,0,is,0,0,0)
            fracs(1) = 1.0
            mixture = solnnum(1,4,0,is,0,0,0)

!$$$            n_user(ET_REACTION) = 1
!$$$            rxnmols = 1e-6

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

            iresult = phr_mix(ID,1, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, n_user,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     update chem matrix - no volume needed if imetric = -1 since
!     the final volume (same as init) is assessed in update_chem.
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,0D0,-1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!            rxnmols = 0.0
!            n_user(ET_REACTION) = -1
            
         end if
!
! Begin UZ mixing. Compute the solute flux from the saturated zone to the
! unsaturated zone resulting from root zone wetting and a lowering water table.
!
         solns(1) = solnnum(0,8,0,is,0,0,0)
         fracs(1) = 1.0
         mixture = solnnum(0,8,0,is,0,0,0)

!         iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $        nac,clark_segs,src_init)
!         if(iresult.ne.0) then
!            PRINT *, 'Errors assigning phreeq entities:'
!            CALL OutputErrorString(id)
!            STOP
!         end if

         iresult = phr_mix(ID,1, solns, fracs,  mixture,&
              fill_factor, mixture, conc, phr_tf, no_rxn,&
              rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
              ntally_cols)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF

!     update chem matrix with solute mass exported from the saturated zone to
!     the unsaturated zone.
         totvol = vmix_sat2uz(is)
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
              ires,ichemdat,imru,inac,ihydro)
         iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!     
!     Unsaturated zone
!     
!     Create new UZ solutions
!     
         do 93 j = 1, nacsc(is)
!
! Determine if UZ bin is in the riparian zone or upland using riparin_threh 
!
           if(st(j,is).ge.riparian_thresh(is))  then
              indxuz=ch_rip_soln(is)
           else
              indxuz=ch_upland_soln(is)
           end if

!         
!        Record evaporation and make initial and final UZ volumes
!        available for reaction accounting

         mixture = solnnum(0,5,0,is,j,0,0) ! t0 UZ bin
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
         c_chem(indx)%vol(init)=vmix_uz(j,is,1)
         c_chem(indx)%vol(ET)=vmix_uz(j,is,6)
         c_chem(indx)%vol(fin)=vmix_uz(j,is,4)
             
            if (vmix_uz(j,is,2).gt.0.0) then

               mixture = solnnum(0,14,0,is,j,0,0) ! Sum inputs in temp res 14
               totvol = vmix_uz(j,is,2)

               solns(1) = solnnum(0,11,0,is,0,0,0) ! precip+irrig
               fracs(1) = (vmix_uz(j,is,5)+vmix_uz(j,is,17)+&
                          vmix_uz(j,is,18)+vmix_uz(j,is,19))/totvol

               solns(2) = solnnum(1,1,0,is,0,0,0) ! throughfall
               fracs(2) = vmix_uz(j,is,8)/totvol

               if(snow_ion_pulse(is).eq.1) then
                 solns(3) = solnnum(0,2,0,is,0,0,0) ! ionic snowmelt stored in t0 solution
               else
                 solns(3) = solnnum(1,2,0,is,0,0,0) ! no ionic pulse, use snowpack chemistry
               end if
               fracs(3) = vmix_uz(j,is,9)/totvol

               solns(4) = solnnum(0,8,0,is,0,0,0) ! flux from wetting and water level change
               fracs(4) = vmix_uz(j,is,14)/totvol

! Check for negative fractions

               iresult=checkfracs(4,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,4, solns, fracs,  mixture,&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
               indx = isoln(solnnum(0,5,0,is,j,0,0),nchemdat,nmru,&
                    nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF

            end if

            mixture = solnnum(1,5,0,is,j,0,0) ! Create new UZ solution
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)

            if (vmix_uz(j,is,1).gt.0.0) then

               totvol = vmix_uz(j,is,1)+vmix_uz(j,is,2)

               solns(1) = solnnum(0,5,0,is,j,0,0) ! Initial UZ solution
               fracs(1) = vmix_uz(j,is,1)/totvol

               solns(2) = solnnum(0,14,0,is,0,0,0) ! UZ inputs
               fracs(2) = vmix_uz(j,is,2)/totvol

!               solns(3) = solnnum(0,0,0,0,0,0,0) ! ET (negative fraction)
!               fracs(3) = -vmix_uz(j,is,6)/totvol

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
!               iresult=checkfracs(3,solns,fracs,mixture)
!               if(iresult.ne.0) then
!                  PRINT*,'Errors with mixing fractions'
!               end if
! Mix
!               iresult = phr_mix(ID,3, solns, fracs,  solns(1),
!     $              fill_factor, mixture, conc, phr_tf, n_user,
!     $              rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,
!     $              ntally_cols)
!
                iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)
              IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
! --------------  iso evap section
!
! If D or 18O are included in list of solutes of interest, compute fractionation of evaporated water.
! Evaporated D and 18O are tracked as basin exports in the fractionate subroutine. 
! if n_iso, the number of isotopes, =0, this section is skipped and DI is evaporated
               if(n_iso.ne.0.and.vmix_uz(j,is,6).gt.0) then
                 evap=vmix_uz(j,is,6)
                 ison=iso_n(is)
                 rh=relhum(1)
!                   iresult=fractionate(water, indx,evap,totvol,ison,rh,tempc)  ! Rayleigh fractionation of evapaporated (or sublimated) water
!                 if(srzwet(j,is).gt.0) then ! water table is close to surface so permit fractionation
                 if(riparian(j,is)) then ! water table is close to surface so permit fractionation
                   iresult=fractionate(water, indx,evap,totvol,ison,rh,tempc)  ! Rayleigh fractionation of evapaporated (or sublimated) water
                 else
                   iresult=fractionate(transp, indx,evap,totvol,ison,rh,tempc)  ! transp keywork results in no fractionation for UZ bins with thick UZ
                 endif
                 IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during transpiration fractionation'
                  CALL OutputErrorString(id)
                  STOP
                 ENDIF
               end if
! --------------  //iso evap section

! remove evaporated water from soil
               totvol = totvol-vmix_uz(j,is,6)
               solns(1) = solnnum(1,5,0,is,j,0,0) ! mixture of initial UZ and canopy inputs created in previous mix
               fracs(1) = (vmix_uz(j,is,1)+vmix_uz(j,is,2))/totvol
               solns(2) = solnnum(0,0,0,0,0,0,0)    ! soln 1: ET is DI water
               fracs(2) = -vmix_uz(j,is,6)/totvol
! fill_ent to allow reactions in this final mix before export (uses final volume)
               iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                    nac,clark_segs,src_init)
               if(iresult.ne.0) then
                  PRINT *, 'Errors assigning phreeq entities:'
                  CALL OutputErrorString(id)
                  STOP
               end if
               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                    fill_factor, mixture, conc, phr_tf, n_user,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!! reset DI to pure water (no 18O or D)
!            iresult = accumulateline(id, "solution 1")  ! solution 1 is always evaporative water
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors establishing DI solution 1'
!               STOP
!            ENDIF
!            iresult = runaccumulated(ID)   ! assigns depleted deltas to evap using lines accumulated above         
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors establishing DI solution 1'
!               STOP
!            ENDIF

!     Update volume/mole matrix for UZ bin
!     Note that part of this output may be transpiration that is
!     tracked in the next section
               totvol = vmix_uz(j,is,3)
               iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
 93      continue               ! End of loni bin loop 
!*********************
! UZ/Canopy scenario 2: Transpire UZ water of final composition to canopy
! and evaporate.
!
! Create transient mixtures of unsaturated zone solutions for transpiration
! and recharge.
!
!        Record ET volume to canopy and make initial and final canopy volumes
!        available for reaction accounting
         mixture = solnnum(1,1,0,is,0,0,0) ! t+1 canopy use for throughfall
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_can(is,1)
               c_chem(indx)%vol(ET)=vmix_can(is,6)
               c_chem(indx)%vol(fin)=vmix_can(is,4)
         if (vmix_can(is,11).gt.0) then
            mixture = solnnum(0,12,0,is,0,0,0) ! transient reservoir
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            validnacs = 0
            totvol = vmix_can(is,11)
            do 26 ij = 1, nacsc(is)
               if(vmix_uz2can(ij,is).gt.0) then
                  validnacs = validnacs + 1
                  solns(validnacs)=solnnum(1,5,0,is,ij,0,0)
                  fracs(validnacs)=vmix_uz2can(ij,is)/totvol
               end if
 26         continue
!
! Check that vmix_can(is,11)>0 works for no canopy too
!

!            iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $           nac,clark_segs,src_init)
!            if(iresult.ne.0) then
!               PRINT *, 'Errors assigning phreeq entities:'
!               CALL OutputErrorString(id)
!               STOP
!            end if

! Check for negative fractions

               iresult=checkfracs(validnacs,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,validnacs, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     Update volume/mole matrix for uz2can transient reservoir
            iresult = update_chem(indx,totvol,0,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     Canopy inputs (use transient reservoir 14)
            mixture= solnnum(0,14,0,is,0,0,0)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            totvol=vmix_can(is,2)

            solns(1) = solnnum(0,11,0,is,0,0,0) ! Precip+irrig
            fracs(1) = (vmix_can(is,5)+vmix_can(is,17)+&
                  vmix_can(is,18)+vmix_can(is,19))/totvol

            solns(2) = solnnum(0,12,0,is,0,0,0) ! Transpiration inputs
            fracs(2) = vmix_can(is,11)/totvol

            solns(3) = solnnum(0,4,0,is,0,0,0) ! Transer from o-horiz, leaves on
            fracs(3) = vmix_can(is,10)/totvol
! Check for negative fractions

               iresult=checkfracs(3,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,3, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
! update inputs            
            mixture= solnnum(0,1,0,is,0,0,0)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
     
            iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                 tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     Mix and concentrate canopy solution
!$$$            write(*,396) (vmix_can(is,ik),ik=1,12)
            totvol = vmix_can(is,1)+vmix_can(is,2)
            mixture = solnnum(1,1,0,is,0,0,0) ! t+1 canopy 
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)

            solns(1) = solnnum(0,1,0,is,0,0,0) ! initial canopy chem
            fracs(1) = vmix_can(is,1)/totvol

            solns(2) = solnnum(0,14,0,is,0,0,0) ! inputs
            fracs(2) = vmix_can(is,2)/totvol


! Check for negative fractions

            iresult=checkfracs(2,solns,fracs,mixture)
            if(iresult.ne.0) then
               PRINT*,'Errors with mixing fractions'
            end if
! Mix
            iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
! --------------  iso evap section
!
! If D or 18O are included in list of solutes of interest, compute fractionation of evaporated water.
! Evaporated D and 18O are tracked as exports in the function call. 
! if n_iso, the number of isotopes, =0, this section is skipped and DI is evaporated
               if(n_iso.ne.0.and.vmix_can(is,6).gt.0) then
                 evap = vmix_can(is,6)
                 ison=iso_n(is)
                 rh=relhum(1)
                 totvol_can = totvol - vmin_canopy(is)
                 iresult=fractionate(water,indx,evap,totvol_can,ison,rh,tempc)  ! Rayleigh fractionation of evaporated water
                 IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during evaporation fractionation'
                  CALL OutputErrorString(id)
                  STOP
                 ENDIF
               end if
! --------------  //iso evap section

! remove evaporated water from canopy
               totvol = totvol-vmix_can(is,6)
               solns(1) = solnnum(0,1,0,is,0,0,0) ! conservative mixture of initial and canopy inputs created in previous mix
               fracs(1) = (vmix_can(is,1)+vmix_can(is,2))/totvol
               solns(2) = solnnum(0,0,0,0,0,0,0)    ! soln 1: ET is DI water
               fracs(2) = -vmix_can(is,6)/totvol
! fill_ent to allow reactions in this final mix before export (uses final volume)
               iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
               if(iresult.ne.0) then
                 PRINT *, 'Errors assigning phreeq entities:'
                 CALL OutputErrorString(id)
                 STOP
               end if
               ! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if

               iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                 fill_factor, mixture, conc, phr_tf, n_user,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
               IF (iresult.NE.0) THEN
                 PRINT *, 'Errors during phr_mix:'
                 CALL OutputErrorString(id)
                 STOP
               ENDIF
!! reset DI to pure water (no 18O or D)
!            iresult = accumulateline(id, "solution 1")  ! solution 1 is always evaporative water
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors establishing DI solution 1'
!               STOP
!            ENDIF
!            iresult = runaccumulated(ID)   ! assigns depleted deltas to evap using lines accumulated above         
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors establishing DI solution 1'
!               STOP
!            ENDIF

!     Update volume/mole matrix for canopy - Note totvol is set to zero because no
!     exported volume. The imetric 3 will insure that the reactions are accounted 
!     for and that the final masses after reaction are assigned.

            iresult = update_chem(indx,0D0,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         end if

!     No transpiration on days of leaves off 
!     Preferential flow through  unsaturated zone next.
!
! If there is recharge (quz>0), then generate a recharge solution from 
! contributions of individual UZ bins. Note that this will not execute
! for recharge to the saturated zone from water bypassing the unsaturated
! zone.
!
         if(quz(is).gt.0.0) then

            mixture = solnnum(1,13,0,is,0,0,0) ! transient recharge reservoir
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            validnacs = 0
            totvol = vmix_qdf(is,2)+vmix_sat(is,11)  ! 11, inputs from uz
            
            do 28 ij = 1, nacsc(is)
               if(vmix_uz2sat(ij,is).gt.0) then
                  validnacs = validnacs + 1
                  solns(validnacs)=solnnum(1,5,0,is,ij,0,0)
                  fracs(validnacs)=vmix_uz2sat(ij,is)/totvol
               end if
 28         continue

!
! New call will have 2 columns of solution concentration: 1st is conservative mix, 2nd
! is after reacting
!
!$$$            iresult = phr_mix(ID,validnacs, solns, fracs, index_conservative,
!$$$     $            fill_factor, index_reaction, *****conc_conservative****,
!$$$     $            phr_tf, n_user,rxnmols,tempc,ph,ph_final,tsec,
!$$$     $            tally_table,ntally_rows,ntally_cols)
!$$$
!            iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $           nac,clark_segs,src_init)
!            if(iresult.ne.0) then
!               PRINT *, 'Errors assigning phreeq entities:'
!               CALL OutputErrorString(id)
!               STOP
!            end if

! Check for negative fractions

               iresult=checkfracs(validnacs,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,validnacs, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,rxnmols,&
                 tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     Update volume/mole matrix for uz2sat transient reservoir
               iresult = update_chem(indx,totvol,0,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!     Update volume/mole matrix for UZ pref reservoir inputs
            indx = isoln(solnnum(0,7,0,is,0,0,0),nchemdat,nmru,&
                 nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
            totvol = vmix_qdf(is,2)
            iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         end if
!
! Preferential flow through unsaturated zone
! Check for fluxes if qdffrac > 0
!
! The following section executes on days when there is
!     preferential flow through the unsaturated zone
!     (vmix_qdf(is,2) > 0)
!
!
!        Record ET volume to canopy and make initial and final canopy volumes
!        available for reaction accounting
         mixture = solnnum(1,7,0,is,0,0,0) ! t+1 pref flow
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_qdf(is,1)
               c_chem(indx)%vol(ET)=vmix_qdf(is,6)
               c_chem(indx)%vol(fin)=vmix_qdf(is,4)

         if(vmix_qdf(is,2).gt.0.0) then ! flux through qdf flowpath
! Inputs chemistry from previous section so no need to
! create combination of inputs
            
            mixture = solnnum(1,7,0,is,0,0,0)
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            totvol = vmix_qdf(is,1)+vmix_qdf(is,2)

!     solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)

            solns(1) = solnnum(0,7,0,is,0,0,0) ! t0 solution
            fracs(1) = vmix_qdf(is,1)/totvol
            
            solns(2) = solnnum(1,13,0,is,0,0,0) ! recharge water
            fracs(2) = vmix_qdf(is,2)/totvol

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                 fill_factor, mixture, conc, phr_tf, n_user,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)


            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     Update volume/mole matrix for qdf outputs
               iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         else                   ! react the residual reservoir if static
            solns(1)=solnnum(0,7,0,is,0,0,0)
            mixture=solnnum(1,7,0,is,0,0,0)
            fracs(1) = 1.0

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

            iresult = phr_mix(ID,1, solns, fracs, solns(1),&
                 fill_factor, mixture, conc, phr_tf, n_user,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!     update chem matrix - no volume needed if imetric = -1 since
!     the initial volume is assessed in update_chem.
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,0D0,-1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         endif
!     
! End of flux calculations for preferential flow through the
! unsaturated zone.
!
! ******************************************
! Saturated Zone fluxes. Execute if saturated zone had
! inputs or outputs.
!     (vmix_sat(is,2) > 0), or
!     (vmix_sat(is,3) > 0).
!
!        Record initial and final volumes for reaction accounting

         mixture = solnnum(1,8,0,is,0,0,0) ! t+1 saturated zone water
         indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
               c_chem(indx)%vol(init)=vmix_sat(is,1)
               c_chem(indx)%vol(ET)=vmix_sat(is,6) ! Should be zero
               c_chem(indx)%vol(fin)=vmix_sat(is,4)

! TOPMODEL always has some baseflow but include check to be complete
!
         if(vmix_sat(is,3).gt.0.0.or.vmix_sat(is,2).gt.0.0) then
!
! Create combination of inputs including water that bypasses the UZ
!
            totvol = vmix_sat(is,2)
            if(totvol.gt.0) then
               mixture = solnnum(1,14,0,is,0,0,0)
               indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)

               solns(1) = solnnum(0,11,0,is,0,0,0) ! precip + irrig
               fracs(1) = vmix_sat(is,5)/totvol

               solns(2) = solnnum(1,1,0,is,0,0,0) ! throughfall
               fracs(2) = vmix_sat(is,8)/totvol

               if(snow_ion_pulse(is).eq.1) then
                 solns(3) = solnnum(0,2,0,is,0,0,0) ! ionic snowmelt stored in t0 solution
               else
                 solns(3) = solnnum(1,2,0,is,0,0,0) ! no ionic pulse, use snowpack chemistry
               end if
               fracs(3) = vmix_sat(is,9)/totvol

               solns(4) = solnnum(1,13,0,is,0,0,0) ! recharge from uz
               fracs(4) = vmix_sat(is,11)/totvol

               ib = src_gw1(is)+1
               solns(5) = solnnum(0,0,ib,0,0,0,0) ! gw input from 1st source (channel leakage?)
               fracs(5) = vmix_sat(is,20)/totvol

               ib = src_gw2(is)+1
               solns(6) = solnnum(0,0,ib,0,0,0,0) ! gw input from 2nd source (upgradient?)
               fracs(6) = vmix_sat(is,21)/totvol

!               iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if

! Check for negative fractions

               iresult=checkfracs(6,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
               iresult = phr_mix(ID,6, solns, fracs,  mixture,&
                    fill_factor, mixture, conc, phr_tf, no_rxn,&
                    rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                    ntally_cols)

               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors during phr_mix:'
                  CALL OutputErrorString(id)
                  STOP
               ENDIF
!     
!     Reset indx to the sat zone
!     
               mixture = solnnum(1,8,0,is,0,0,0) ! t+1 saturated zone water
               indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)

!
!     Update c_chem with input volumes and masses
!
                    iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if

! totvol for mixing now includes bypass volumes and is less the sat2uz flux and
! the amount of well irrigation 

            totvol = vmix_sat(is,1)-vmix_sat2uz(is)-&
                 vmix_well(is)+vmix_sat(is,2)

            solns(1) = solnnum(0,8,0,is,0,0,0) ! t0 solution
            fracs(1) = (vmix_sat(is,1)-vmix_sat2uz(is)-&
                 vmix_well(is))/totvol
            
            solns(2) = solnnum(1,14,0,is,0,0,0) ! inputs
            fracs(2) = vmix_sat(is,2)/totvol

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
            iresult = phr_mix(ID,2, solns, fracs, solns(1),fill_factor, &
                 mixture, conc, phr_tf, n_user,rxnmols,tempc,ph,ph_final,tsec,&
                 tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!
! Update c_chem with output volumes and masses (to stream) computed at new concentration.
! Reactions occur only in final volume: c_chem(indx)%vol(fin) = vmix_sat(is,4).
!
            totvol = vmix_sat(is,3)-vmix_sat2uz(is)-vmix_well(is)
            
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
!
! Update c_chem with satpref inputs. We are currently inside the loop
! where sat zone outputs, vmix_sat(is,3) > 0
!
            if(vmix_satpref(is,2).gt.0) then
               mixture=solnnum(0,9,0,is,0,0,0)
               totvol = vmix_satpref(is,2)
               indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                    ires,ichemdat,imru,inac,ihydro)
               iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
               IF (iresult.NE.0) THEN
                  PRINT *, 'Errors updating mole matrix:'
                  STOP
               ENDIF
            end if
         else                   ! react the reservoir if static
            solns(1) = solnnum(0,8,0,is,0,0,0)
            mixture= solnnum(1,8,0,is,0,0,0)
            fracs(1) = 1.0

            iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
                 nac,clark_segs,src_init)
            if(iresult.ne.0) then
               PRINT *, 'Errors assigning phreeq entities:'
               CALL OutputErrorString(id)
               STOP
            end if

            iresult = phr_mix(ID,1, solns, fracs,&
                 solns(1),fill_factor, mixture, conc, phr_tf, &
                 n_user,rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
! update chem matrix - no volume needed if imetric = -1 since
! the final volume (same as init) is assessed in update_chem.
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,0D0,-1,conc,pH,tempc,&
                 tally_table,n_ent,indxm,indxb,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
         endif
!
! Account for groundwater leakage, gw_loss, as exports from
! MRU and basin pseudo solutions.  Conc should reflect last
! non-reacted mix from the saturated zone above.
!
!!! Risky assumption
!
        totvol = gw_loss(is)*mru_area(is)*a_million ! cubic meters
        if(totvol.gt.0.0) then
          c_chem(indxm)%vol(out)= c_chem(indxm)%vol(out)+totvol
          c_chem(indxb)%vol(out)= c_chem(indxb)%vol(out)+totvol
          do k=1,nsolute
            c_chem(indxm)%M(k,out)=&
              c_chem(indxm)%M(k,out)+ conc(k)*totvol*a_thousand   ! conc in Moles/L
            c_chem(indxb)%M(k,out)=&
              c_chem(indxb)%M(k,out)+ conc(k)*totvol*a_thousand
            ch_mru_out_permil(is,k)=ch_mru_out_permil(is,k)+&
              c_chem(indxm)%delta(k,out)*totvol
            ch_basin_out_permil(k)=ch_basin_out_permil(k)+&
              c_chem(indxb)%delta(k,out)*totvol
          end do
        end if ! totvol.gt.0


!************** print to select.out *********
!     
!     Test the selected output for the GW of MRU 1
! output solution described in selected output
!      if(is.eq.1) then
!         cols = GetSelectedOutputColumnCount()
!         PRINT *, 'Saturated zone solution for MRU: ',is
!            iresult = GetSelectedOutputValue
!     $           (0, 1, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), ph
!         
!            iresult = GetSelectedOutputValue
!     $           (0, 2, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), tempc
!         DO 80 i = 3,cols
!            iresult = GetSelectedOutputValue
!     $           (0, i, vtype, dvalue, heading)
!            IF (iresult.NE.0) THEN
!               PRINT *, 'Errors during GetSelectedOutputValue:'
!               CALL OutputErrorString(id)
!               STOP
!            ENDIF
!            leng = INDEX(heading, ' ')
!            PRINT *, heading(1:leng-1), ACHAR(9), conc(i-2)
! 80      CONTINUE
!      endif
!************** print stop **********
! Preferential flow through the saturated zone. Execute if
! a preferential flow path exists (qpref_max > 0)
!
!       Record initial and final volumes for reaction accounting

        mixture = solnnum(1,9,0,is,0,0,0) ! t+1 preferential zone water
        indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                ires,ichemdat,imru,inac,ihydro)
        c_chem(indx)%vol(init)=vmix_satpref(is,1)
        c_chem(indx)%vol(ET)=vmix_satpref(is,6) ! Should be zero
        c_chem(indx)%vol(fin)=vmix_satpref(is,4)

        iresult = fill_ent(n_user,mixture,nchemdat,nmru,&
           nac,clark_segs,src_init)
        if(iresult.ne.0) then
           PRINT *, 'Errors assigning phreeq entities:'
           CALL OutputErrorString(id)
           STOP
        end if

        if(qpref_max(is).gt.0.0) then ! flowpath exists
          if(vmix_satpref(is,2).gt.0.0) then ! there is a flux
!     Inputs chemistry from previous section so no need to
!     create combination of inputs
            totvol = vmix_satpref(is,1)+vmix_satpref(is,2)

!     solnnum(time,res_ID,chemdat, mru, nac,hydro, stat)

            solns(1) = solnnum(0,9,0,is,0,0,0) ! t0 solution
            fracs(1) = vmix_satpref(is,1)/totvol

            solns(2) = solnnum(1,8,0,is,0,0,0) ! sat zone water
            fracs(2) = vmix_satpref(is,2)/totvol

! Check for negative fractions

            iresult=checkfracs(2,solns,fracs,mixture)
            if(iresult.ne.0) then
               PRINT*,'Errors with mixing fractions'
            end if
! Mix
            iresult = phr_mix(ID,2, solns, fracs,  solns(1),&
                 fill_factor, mixture, conc, phr_tf, n_user,rxnmols,&
                 tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!
! Update c_chem with output volumes and masses
            totvol = vmix_satpref(is,3)
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
          else                   ! react the reservoir if static
            solns(1) = solnnum(0,9,0,is,0,0,0)
!            mixture=solnnum(1,9,0,is,0,0,0)
            fracs(1) = 1.0

!            iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $           nac,clark_segs,src_init)
!            if(iresult.ne.0) then
!               PRINT *, 'Errors assigning phreeq entities:'
!               CALL OutputErrorString(id)
!               STOP
!            end if
!
            iresult = phr_mix(ID,1, solns, fracs,mixture,&
                 fill_factor, mixture, conc, phr_tf, n_user,rxnmols,&
                 tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
! update chem matrix - no volume needed if imetric = -1 since
! the final volume (same a init) is assessed in update_chem.
!            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,
!     $           ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,0D0,-1,conc,pH,tempc,&
                               tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
          endif ! vmix_satpref(is,2) > 0
        end if ! qprefmax > 0
!
! End of saturated zone flux calculations
!

!
! Combine outflows into transient hillslope reservoir.
!     
! The following section executes on days when there is
!     flux from the hillslope to the drainage system
!     (vmix_hill(is,3) >0). This should be every step.
!
        if(vmix_hill(is,3).gt.0.0) then

! Inputs: No need to use temporary input reservoir number 14
! since this is a transient reservoir where inputs = output
! and there is no residual volume. This may be modified in the
! future to become a riparian reservior. In that case the reservoir
! would have an area so precip and ET would exist along with a
! residual volume so res 15 and its associated two-part mixing
! would be needed.
!
          mixture = solnnum(1,10,0,is,0,0,0)
          totvol = vmix_hill(is,2)

          if (totvol.gt.0.0) then ! Outputs to stream
               
!$$$               solns(1) = solnnum(0,1,0,is,0,0,0) ! precip inputs
!$$$               fracs(1) = vmix_hill(is,5)/totvol
               
            solns(1) = solnnum(1,4,0,is,0,0,0) ! O-horizon
            fracs(1) = vmix_hill(is,10)/totvol

            solns(2) = solnnum(1,7,0,is,0,0,0) ! UZ pref flow
            fracs(2) = vmix_hill(is,12)/totvol

            solns(3) = solnnum(1,8,0,is,0,0,0) ! saturated zone
            fracs(3) = (vmix_hill(is,13)+vmix_hill(is,14))/totvol ! baseflow +exfiltration

            solns(4) = solnnum(1,9,0,is,0,0,0) ! Sat pref flow
            fracs(4) = vmix_hill(is,15)/totvol
!
!               iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $              nac,clark_segs,src_init)
!               if(iresult.ne.0) then
!                  PRINT *, 'Errors assigning phreeq entities:'
!                  CALL OutputErrorString(id)
!                  STOP
!               end if

! Check for negative fractions

            iresult=checkfracs(4,solns,fracs,mixture)
            if(iresult.ne.0) then
               PRINT*,'Errors with mixing fractions'
            end if
! Mix
            iresult = phr_mix(ID,4, solns, fracs,  mixture,&
                 fill_factor, mixture, conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)

            IF (iresult.NE.0) THEN
               PRINT *, 'Errors during phr_mix:'
               CALL OutputErrorString(id)
               STOP
            ENDIF
!$$$         write(*,198)(datetime(j),j=1,3), mixture,
!$$$     $        (solns(ik),fracs(ik),ik=1,5)
!$$$         write(*,199)(conc(ik),ik=1,nsolute)

!     Update volume/mole matrix in two calls to allow mru and basin accounting of transient hillslope outputs
            indx = isoln(mixture,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
            iresult = update_chem(indx,totvol,2,conc,pH,tempc,&
                            tally_table,n_ent,0,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
            iresult = update_chem(indx,totvol,3,conc,pH,tempc,&
                            tally_table,n_ent,indxm,0,ires)
            IF (iresult.NE.0) THEN
               PRINT *, 'Errors updating mole matrix:'
               STOP
            ENDIF
            mru_out_vol(is)=mru_out_vol(is)+totvol
            do k=1,nsolute
              ch_mru_out_permil(is,k)=ch_mru_out_permil(is,k)+&
                   c_chem(indxm)%delta(k,out)*totvol
            end do
          endif ! totvol.gt.0.0

! Add section here if the transient hill
! reservoir becomes a riparian reservoir.
        end if
!     
! End of transient hillslope fluxes
!
 10   continue ! end MRU loop

! Distribute output to stream segments -
! In webmod_res, the hillslope solutions are delivered to
! the stream segments and then the segments are moved one
! time step closer to the outlet. That is why there are
! operations between the hillslope inputs and the segments
! closer to the outlet in the routine below.
!
! The solute mass described by the stream solution at
! t0 receives contributions from all MRUs and then
! is moved one hydro segment closer to the outlet during
! each time step. There is no dispersion. So the t+1 volume
! and mass reflect upstream influent exclusively.
!
! The iterations are completed from the farthest segment first
! to the basin outlet last. Contributions to the segment ij=1
! were exported from the basin on this time step in webmod_res.
! The exported volume is the simulated exportvmix_basin(3) minus the stream and gw loss.
! This equals basin_qsim_cm

      do 20 ij = clark_segs,1,-1

         if(ij.gt.1) then
            str_vol = vmix_stream(ij-1)   ! previous volume + inputs
         else
            str_vol = basin_qsim_cm*1e4*basin_area
         end if

         mixture = solnnum(0,14,0,1,0,0,0) ! use temp res 14, MRU 1
         totvol = 0.0
         do 22 is = 1,nmru
            totvol = totvol + vmix_hillexp(ij,is)  ! hillslope inputs
 22      continue
         do 24 is = 1,nmru
            solns(is) = solnnum(1,10,0,is,0,0,0)
            fracs(is) = vmix_hillexp(ij,is)/totvol
            if(fracs(is).lt.0)fracs(is)=0.0
 24      continue
!         iresult = fill_ent(n_user,mixture,nchemdat,nmru,
!     $        nac,clark_segs,src_init)
!         if(iresult.ne.0) then
!            PRINT *, 'Errors assigning phreeq entities:'
!            CALL OutputErrorString(id)
!            STOP
!         end if

! Check for negative fractions

               iresult=checkfracs(nmru,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
         iresult = phr_mix(ID,nmru, solns, fracs,  mixture,&
              fill_factor, mixture, conc, phr_tf, no_rxn,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
!$$$         write(*,198)(datetime(j),j=1,3), mixture,
!$$$     $        (solns(ik),fracs(ik),ik=1,nmru)
!$$$         write(*,199)(conc(ik),ik=1,nsolute)

! Record stream inputs in c_chem. Append with upstream later
         indx = isoln(solnnum(0,99,0,0,0,ij,0),nchemdat,nmru,&
              nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
         iresult = update_chem(indx,totvol,2,conc,pH_final,tempc,&
                               tally_table,n_ent,0,0,ires)
         IF (iresult.NE.0) THEN
            PRINT *, 'Errors updating mole matrix:'
            STOP
         ENDIF
!
!     Mix hillslope with previous solution
!     
         mixture = solnnum(1,99,0,0,0,ij-1,0) ! defaults to basin at ij=1

         solns(1) = solnnum(0,14,0,1,0,0,0) ! hillslope inputs
         fracs(1) = totvol/str_vol
         if(fracs(1).lt.0)then
            print*,'negative fraction for hillslope input ',&
                 ' set to zero. It was equal to ',fracs(1)
            fracs(1)=0.0
         end if
         
         totvol = str_vol-totvol ! total - inputs = orig vol of segment
         
         solns(2) = solnnum(0,99,0,0,0,ij,0) ! orig chem in upstream section
         fracs(2) = totvol/str_vol
         if(fracs(2).lt.0)then
            print*,'negative fraction for stream seg ',&
                 ij,' set to zero. It was equal to ',fracs(2)
            fracs(2)=0.0
         end if

! fill entities with those of the upstream hydro segment (solns(2))
         iresult = fill_ent(n_user,solns(2),nchemdat,nmru,&
              nac,clark_segs,src_init)
         if(iresult.ne.0) then
            PRINT *, 'Errors assigning phreeq entities:'
            CALL OutputErrorString(id)
            STOP
         end if

! Check for negative fractions

               iresult=checkfracs(2,solns,fracs,mixture)
               if(iresult.ne.0) then
                  PRINT*,'Errors with mixing fractions'
               end if
! Mix
         iresult = phr_mix(ID,2, solns, fracs,  solns(2),&
              fill_factor, mixture, conc, phr_tf, n_user,rxnmols,&
              tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)

         IF (iresult.NE.0) THEN
            PRINT *, 'Errors during phr_mix:'
            CALL OutputErrorString(id)
            STOP
         ENDIF
!$$$         write(*,198)(datetime(j),j=1,3), mixture,
!$$$     $        (solns(ik),fracs(ik),ik=1,2)
!$$$         write(*,199)(conc(ik),ik=1,nsolute)

!
! Add stream diversions and channel losses to basin exports
!
        totvol = vmix_diversion(ij)+vmix_chan_loss(ij)
        if(totvol.gt.0.0) then
          c_chem(chbas_soln)%vol(out)= &
            c_chem(chbas_soln)%vol(out)+totvol
            basin_out_vol=basin_out_vol+totvol
          do k=1,nsolute
            c_chem(chbas_soln)%M(k,out)=&
              c_chem(chbas_soln)%M(k,out)+ conc(k)*totvol*a_thousand
            ch_basin_out_permil(k)=ch_basin_out_permil(k)+&
                  c_chem(indxb)%delta(k,out)*totvol
          end do
        end if ! totvol.gt.0

!
! Update these outputs and record as inputs to the next segment.
! If the outlet segment (ij=1) record as basin output.
!
        if (ij.gt.1) then
!
! Update c_chem to reflect outputs
!
          iresult = update_chem(indx,str_vol,3,conc,pH_final,tempc,&
                                tally_table,n_ent,0,0,ires)
          IF (iresult.NE.0) THEN
             PRINT *, 'Errors updating mole matrix:'
             STOP
          ENDIF
          indx = isoln(solnnum(0,99,0,0,0,ij-1,0),nchemdat,nmru,&
                 nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
          iresult = update_chem(indx,str_vol,2,conc,pH_final,tempc,&
                                   tally_table,n_ent,0,0,ires)
          IF (iresult.NE.0) THEN
             PRINT *, 'Errors updating mole matrix:'
             STOP
          ENDIF
        else   ! export last segment
! Mix once to erase reactions from last stream segment  
!            mixture = isoln(solnnum(0,99,0,0,0,1,0),nchemdat,nmru,
!     $           nac,clark_segs,ires,ichemdat,imru,inac,ihydro)
!            indx = chbas_soln
!            solns(1) = mixture
!            fracs(1) = 1.0
!            iresult = phr_mix(ID,1, solns, fracs,  mixture,
!     $        fill_factor, mixture, conc, phr_tf, no_rxn,rxnmols,
!     $        tempc,ph,ph_final,tsec,tally_table,ntally_rows,ntally_cols)
          iresult = update_chem(indx,str_vol,3,conc,pH_final,tempc,&
                             tally_table,n_ent,0,indxb,ires)
          IF (iresult.NE.0) THEN
             PRINT *, 'Errors updating mole matrix:'
             STOP
          ENDIF
          basin_out_vol=basin_out_vol+str_vol
          do k=1,nsolute
            ch_basin_out_permil(k)=ch_basin_out_permil(k)+&
                  c_chem(indxb)%delta(k,out)*str_vol
          end do
        end if

 20   continue                  ! end clark_segs loop


! Assign composite uz volumes (uzgen), then compute final values for each 
! reservoir, (ignoring temp and precip)
!
      do 200 k = 1,nmru ! composite uz variables
        c_chem_uzgen(k)%vol(init) = vmix_uzgen(k,1)
        c_chem_uzgen(k)%vol(in) = vmix_uzgen(k,2)
        c_chem_uzgen(k)%vol(out) = vmix_uzgen(k,3)
        c_chem_uzgen(k)%vol(ET) = vmix_uzgen(k,6)
!        c_chem(i)%vol(fin) = vmix_uzgen(k,4) ! Calculated below
        c_chem_uzrip(k)%vol(init) = vmix_uzrip(k,1)
        c_chem_uzrip(k)%vol(in) = vmix_uzrip(k,2)
        c_chem_uzrip(k)%vol(out) = vmix_uzrip(k,3)
        c_chem_uzrip(k)%vol(ET) = vmix_uzrip(k,6)
!        c_chem(i)%vol(fin) = vmix_uzrip(k,4) ! Calculated below
        c_chem_uzup(k)%vol(init) = vmix_uzup(k,1)
        c_chem_uzup(k)%vol(in) = vmix_uzup(k,2)
        c_chem_uzup(k)%vol(out) = vmix_uzup(k,3)
        c_chem_uzup(k)%vol(ET) = vmix_uzup(k,6)
!        c_chem(i)%vol(fin) = vmix_uzup(k,4) ! Calculated below
 200  continue
        k=nchemdat+2 !DI water and concentrations in the chemdat file do not have volumes associated
!                     with them so track volumes beginning with the canopy of the first MRU
      do 205 i = k,nphrsolns-nmru-1 ! last indices (nmru+1) are basin and mru assigned below
        minvol = .TRUE.
        c_chem(i)%vol(fin) = c_chem(i)%vol(init)+&
          c_chem(i)%vol(in)-c_chem(i)%vol(out)-&
          c_chem(i)%vol(ET)
        if(abs(c_chem(i)%vol(fin)).lt.0.001D0) minvol = .FALSE.  ! c_chem volumes in cubic meters
!        if(.not.water) print*,'final volume for ',c_indx(i,1),
!     $    ' less than 1 liter on time step ',nstep,'. final '//
!     $    'concentrations and volumes set to zero.'
        do 205 j = 1, nsolute
          if(.not.minvol) then
            c_chem(i)%M(j,fin) = 0D0
            c_chem(i)%vol(fin) = 0D0
          else if(c_chem(i)%M(j,fin).eq.0D0) then ! This updates final masses for conservative mixing. If reactions involved, update_chem has already set the final value.
            c_chem(i)%M(j,fin)=c_chem(i)%M(j,init) + c_chem(i)%M(j,in) -&
            c_chem(i)%M(j,out)+ c_chem(i)%M(j,rxn)
          end if
 205  continue
!
! Basin and MRU summaries
!

      i = chbas_soln
      c_chem(i)%vol(init) = vmix_basin(1)
      c_chem(i)%vol(in) = vmix_basin(2)
      c_chem(i)%vol(out) = vmix_basin(3)
      c_chem(i)%vol(ET) = vmix_basin(6)
      c_chem(i)%vol(fin) = vmix_basin(4)
      do 206 j = 1, nsolute
        c_chem(i)%M(j,fin) =&
         c_chem(i)%M(j,init) + c_chem(i)%M(j,in)-c_chem(i)%M(j,out)+&
         c_chem(i)%M(j,rxn)
206   continue
      vol_in = c_chem(i)%vol(in)
      vol_out = c_chem(i)%vol(out)
      ch_basin_vol_m3 = c_chem(i)%vol(fin)
      ch_basin_out_tempC = c_chem(i)%Temp(fin)
      ch_basin_out_pH = c_chem(i)%pH(fin)

      do n = 1,nsolute
         ch_basin_mass_g(n) =&
              c_chem(i)%M(n,fin)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
         ch_basin_rxn_g(n) =&
              c_chem(i)%M(n,rxn)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
         ch_basin_in_g(n) =&
              c_chem(i)%M(n,in)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
         ch_basin_in_gm2(n) =&
              ch_basin_in_g(n)/basin_area/a_million
         ch_basin_out_g(n) =&
              c_chem(i)%M(n,out)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
         ch_basin_out_gm2(n) =&
              ch_basin_out_g(n)/basin_area/a_million
         ch_basin_net_g(n) =&
              ch_basin_out_g(n)-ch_basin_in_g(n)
         ch_basin_net_gm2(n) =&
              ch_basin_out_gm2(n)-ch_basin_in_gm2(n)
         if(ch_basin_vol_m3.gt.0.0)then
            ch_basin_conc_mgL(n) = ch_basin_mass_g(n)/ch_basin_vol_m3
         else
            ch_basin_conc_mgL(n) = 0.0
         end if
         if(vol_in.gt.0.0)then
!           print*, 'Diff in basin vol_in and totvols = ',
!     $    vol_in - basin_in_vol, (vol_in - basin_in_vol)/basin_in_vol
           ch_basin_in_mgL(n) = ch_basin_in_g(n)/vol_in
! permil inputs and outputs were weighted by volumes so divide by total volume to get mean    
           ch_basin_in_permil(n)= ch_basin_in_permil(n)/basin_in_vol
         else
            ch_basin_in_mgL(n) = 0.0
            ch_basin_in_permil(n)=0D0
         end if

         if(vol_out.gt.0.0)then
!             print*, 'Diff in basin vol_out and totvols = ',
!     $  vol_out - basin_out_vol,(vol_out - basin_out_vol)/basin_out_vol
            ch_basin_out_mgL(n) =&
             ch_basin_out_g(n)/vol_out
! permil inputs and outputs were weighted by volumes so divide by total volume to get mean    
            ch_basin_out_permil(n)= ch_basin_out_permil(n)/basin_out_vol
         else
            ch_basin_out_mgL(n) = 0.0
            ch_basin_out_permil(n)=0D0
         end if
      end do


      do 207 j = 1,nmru ! mru variables
        i = chmru_soln(j)
        c_chem(i)%vol(init) = vmix_mru(j,1)
        c_chem(i)%vol(in) = vmix_mru(j,2)
        c_chem(i)%vol(out) = vmix_mru(j,3)
        c_chem(i)%vol(ET) = vmix_mru(j,6)
        c_chem(i)%vol(fin) = vmix_mru(j,4)
        ch_mru_vol_m3(j) = c_chem(i)%vol(fin)
        ch_mru_out_tempC(j) = c_chem(i)%Temp(fin)
        ch_mru_out_pH(j) = c_chem(i)%pH(fin)
        vol_in = c_chem(i)%vol(in)
        vol_out = c_chem(i)%vol(out)

        do 207 n = 1, nsolute
          c_chem(i)%M(n,fin) =&
           c_chem(i)%M(n,init) + c_chem(i)%M(n,in)-c_chem(i)%M(n,out)+&
             c_chem(i)%M(n,rxn)
           ch_mru_mass_g(j,n) =&
             c_chem(i)%M(n,fin)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
           ch_mru_rxn_g(j,n) =&
             c_chem(i)%M(n,rxn)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
           ch_mru_in_g(j,n) =&
             c_chem(i)%M(n,in)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
           ch_mru_in_gm2(j,n) =&
             ch_mru_in_g(j,n)/mru_area(j)/a_million
           ch_mru_out_g(j,n) =&
             c_chem(i)%M(n,out)*phq_lut(sol_id(n)%phq)%M2mg/a_thousand
           ch_mru_out_gm2(j,n) =&
             ch_mru_out_g(j,n)/mru_area(j)/a_million
           ch_mru_net_g(j,n) =&
             ch_mru_out_g(j,n)-ch_mru_in_g(j,n)
           ch_mru_net_gm2(j,n) =&
             ch_mru_out_gm2(j,n)-ch_mru_in_gm2(j,n)
           if(ch_mru_vol_m3(j).gt.0.0)then
             ch_mru_conc_mgL(j,n)=ch_mru_mass_g(j,n)/ch_mru_vol_m3(j)
           else
             ch_mru_conc_mgL(j,n)=0.0
           end if
           if(vol_in.gt.0.0)then
!             print*, 'MRU Diff in vol_in and totvols = ',
!     $    vol_in - mru_in_vol(j),(vol_in - mru_in_vol(j))/mru_in_vol(j)
             ch_mru_in_mgL(j,n) = ch_mru_in_g(j,n)/vol_in
! permil inputs and outputs were weighted by volumes so divide by total volume to get mean    
             ch_mru_in_permil(j,n)= ch_mru_in_permil(j,n)/mru_in_vol(j)
           else
             ch_mru_in_mgL(j,n) = 0.0
             ch_mru_in_permil(j,n)=0D0
            end if
           if(vol_out.gt.0.0)then
!            print*, 'MRU Diff in vol_out and totvols = ',
!     $      vol_out - mru_out_vol(j),(vol_out - mru_out_vol(j))/vol_out
            ch_mru_out_mgL(j,n) = ch_mru_out_g(j,n)/vol_out
            ch_mru_out_permil(j,n)=ch_mru_out_permil(j,n)/mru_out_vol(j)
           else
             ch_mru_out_mgL(j,n) = 0.0
             ch_mru_out_permil(j,n)=0D0
           end if
 207  continue

!
! Update chvar area conversions to loads for canopy and/or snowpack reservoirs
! Use temporary copy of conversion table to avoid rounding errors when converting
! back
!
      do 208 i = 1,nchemvar
         if(chvar_lut(i,4).gt.2.and.chvar_lut(i,4).lt.7) then ! load
            if (chvar_lut(i,5).eq.1) then ! canopy
               covden = covden_win(chvar_lut(i,7))
               if(transp_on(chvar_lut(i,7)).eq.1)&
                   covden = covden_sum(chvar_lut(i,7))
               if(covden.gt.0.0) then
                 chvar_conv_t(i,nsolute+1) = &
                  chvar_conv(i,nsolute+1) * covden
               else
                 chvar_conv_t(i,nsolute+1) = &
                  chvar_conv(i,nsolute+1)
               end if
             else if (chvar_lut(i,5).eq.2) then ! snowpack
               if(snowcov_area(chvar_lut(i,7)).gt.0.0) then
                 chvar_conv_t(i,nsolute+1)=chvar_conv(i,nsolute+1)*&
                    snowcov_area(chvar_lut(i,7))
               else
                 chvar_conv_t(i,nsolute+1) = &
                  chvar_conv(i,nsolute+1)
               end if
             else
               chvar_conv_t(i,nsolute+1) = chvar_conv(i,nsolute+1)
             end if
         else
            chvar_conv_t(i,nsolute+1) = chvar_conv(i,nsolute+1)
         end if
         do 208 j = 1,nsolute
            chvar_conv_t(i,j) = chvar_conv(i,j)
 208  continue
!
! Populate ch_vars
!
      iresult = chem2var()
      IF (iresult.NE.0) THEN
         PRINT *, 'Errors assigning ch_var values:'
         STOP
      ENDIF

      end if                    ! landing point if no chemical simulations to be done

      phreeqmms_run = 0
 
      return

 100  FORMAT(A,I10)
 110  FORMAT(A)
 120  FORMAT(A,1PG15.7E2)

 198  format(/,3i4,i10,6(/,i10,1x,f8.4))
 199  format(3(e10.4,1X))
 295  format(63(e10.4,1X))
 390  format(10f10.4)
 395  format(3i5,36(e12.5,1X))
 396  format(/,2(6e12.5))
 398  format(/,(I10,4(e12.5,1X)))
 1000 format('-Headings ',A)
 1050 format('10 PUNCH ',A)
      end
!***********************************************************************
!
!     phreeqmms_clean - Close the chemout file
!

      integer function phreeqmms_clean()

      USE WEBMOD_PHREEQ_MMS
      USE WEBMOD_IO, only: chemout_file_unit

      phreeqmms_clean = 1

      close (chemout_file_unit)
      close (unit=14)
      close (unit=16)
      close (unit=17)

      phreeqmms_clean = 0

      return
      end


!
! Lookup the PHREEQC internal ID
!
!      integer function phreeq_id(STRING, sol_lut)
!      CHARACTER*12 STRING 
!      character*12 sol_lut(50,2)
!      integer i, length, len1, len2
!      logical test1
!      phreeq_id = 0
!      i = 0
!      do while (phreeq_id.eq.0)
!         i = i+1
!         len1 = length(string)
!         len2 = length(sol_lut(i,1))
!         test1 = (string(1:len1) .eq. (sol_lut(i,1)(1:len2)))
!         if(STRING(1:len1) .eq. (sol_lut(i,1)(1:len2))) phreeq_id = i
!         if(i.gt.50) then
!            print*,"Species ", string(1:len1) ," cannot be "//
!     $           "found in the species lookup table"
!            print*," "
!            phreeq_id = -99
!         end if
!      end do
!      RETURN
!      END
!!
!! Function to return the length of the string using the standard
!! intrinsic FORTRAN LEN function
!!
!      INTEGER FUNCTION LENGTH(STRING)
!      INTEGER I
!      CHARACTER*(*) STRING 
!      DO 15, I = LEN(STRING), 1, -1 
!         IF(STRING(I:I) .NE. ' '.and.STRING(I:I)
!     $        .ne.'\f') GO TO 16 
! 15   CONTINUE 
! 16   LENGTH = I
!      RETURN
!      END
!
!
! Function to return an integer for the solution number.
!
! solnnum(time,reservoir ID,chemdat, mru,nac,hydro,stat)
!
! where
! time
!   = 0. solution at beginning of time step
!   = 1, solution at end of time step

! reservoir ID
!   =  1, canopy,
!   =  2, snowpack,
!   =  3, impermeable surface,
!   =  4, O-horizon,
!   =  5, unsaturated zone,
!   =  6, generalized unsaturated zone,
!   =  7, preferential flow (lateral) in the unsaturated zone,
!   =  8, saturated zone,
!   =  9, preferential flow in through the saturated zone,
!   = 10, combined MRU or basin when MRU=0,
!   = 11, combined hillslope inputs (transient),
!   = 12, transpired unsaturated zone water (transient),
!   = 13, recharge to saturated zone (transient),
!   = 14, mixture of inputs before mixing with reservoir
!   = 99, stream, pond, lake, diversion, or other drainage feature.

! chemdat
!   = unique ID for source of chemical deposition or an observation site

! mru
!   = model response unit

! nac
!   = topographic index bin (stindx in the ch_var list dimensioned by nchemvar)

! hydro
!   = unique ID for drainage feature

! stat
!   = summary flag for riparian, upslope, etc
!
! Numbering scheme
!
!
! Create hillslope reservoir ID's using 100,000,000 series for inital
! solutions and 200,000,000 series for solution concentrations at 
! end of time step. The reservoir type determines where the qualifying
! parameters will be obtained based on the dimension of the reservoir.
! the riparian dimension, nrip, is introduced for later discretization if needed.
!
! Reservoir solution ID (no commas in int*4 vector)  
!
! trrmmmnns where
! t = Time stamp, where
! 1 = beginning solution, and
! 2 = ending solution
!
! trrmmmnns where
!  rr = Hillslope reservoir type, where
!  00 = Input (precip, external chemistry, or observation)
!  01 = canopy interception (nmru),
!  02 = snowpack (nmru),
!  03 = impermeable surfaces (nmru, not active),
!  04 = O-horizon (nmru,nrip)
!  05 = unsaturated zone (nmru, nac),
!  06 = combined unsaturated zone
!  07 = unsaturated preferential flow (nmru,nrip),
!  08 = saturated zone (nmru,nrip),
!  09 = saturated zone preferential flow (nmru,nrip),
!  10 = pseudo reservoir for MRU or basin when mru = 0
!  11 = combined transient reservoir for all hillslopes inputs (nmru).
!  12 = transient reservoir for transpiration
!  13 = transient reservoir for recharge
!  14 = transient reservoir for mixing inputs

!  99 = stream, pond, lake, diversion, or other drainage feature.
!
! trrmmmnns where
!    mmm = MRU id
!
! trrmmmnns where
!       nn  = topographic index id
!
! trrmmmnns where
!         s   = optional dimension flag 
!               (0=basin,1=riparian, 3=upslope, etc)
!
!
! Create hydraulic segment ID's using 300,000,000 series for inital
! solutions and 400,000,000 series for solution concentrations at 
! end of time step.
! txxxhhhhh
! t   = 3, beginning solution
!       4, ending solution
!  xxx = placeholder for additional drainage dimension (000)
!  hh,hhh = hydraulic segment (stream, lake, diversion, etc)
!
! Basin and MRU summary ID's will be 500,000,000 series. No t+1 is
! needed since phreeqc never computes mru or basin compositions;
! they are simply the mathematical summation of distinct phreeqc reservoirs
! t00000mmm
! t   = 5, composite solution
! 00000, placeholder for additional mrus
! mmm = MRU ID, zero if summary for entire basin.
!

      INTEGER FUNCTION solnnum&
           (time,res_id,chemdat,mru,nac,hydro,stat)
      implicit none
      INTEGER time, res_id,chemdat,mru,nac,hydro,stat
      integer itime
      solnnum = -99
      itime = time*int(1e8)
      if(res_id.lt.0.or.res_id.gt.99) then
         write(*,100)
         return
      end if
      if (res_id.eq.99)then ! drainage segment
         solnnum = itime + int(3e8) + hydro
      else if (res_id.eq.15)then !Basin or MRU pseudo solution
         solnnum = itime + int(5e8) + mru
      else if (res_id.eq.0)then ! DI, precip, external chem, or observation
            solnnum = 1e3 + chemdat
            if(chemdat.eq.0) solnnum = 1  ! DI water for ET
      else if (res_id.le.14) then !one of the hillslope reservoirs
         solnnum = itime + int(1e8) + res_id*int(1e6) + mru*int(1e3)
      else ! reservoir id not identified in model 
         write(*,100)time,res_id,chemdat,mru,nac,hydro,stat
         return
      end if
      if (res_id.eq.4) solnnum = solnnum+stat !O-horizon
      if (res_id.eq.5) solnnum = solnnum+nac*10 !individual unsat zones
      if (res_id.ge.6.and.res_id.le.9) solnnum=solnnum+stat !hillslope division
      RETURN
 100  format('Problem assigning solution ID for : ',&
        'Time,res,chemdat,mru,nac,hydro,stat',7I4)
      END
!
! The isoln function below will need to be revised to include the stat
! value indicating hillslope versus riparian split if implemented. For now
! the composite UZ is stored in reservoir ID 6
!
! Modified to return only indx values for t0 solutions
!
      INTEGER FUNCTION isoln(soln_num,nchemdat,nmru,nac,nhydro,&
           ires,ichemdat,imru,inac,ihydro)

      implicit none
      INTEGER soln_num
      integer nmru, nac, nhydro, nchemdat
      integer icut, ires,imru, inac, ihydro,ichemdat

      isoln = 0

      ires = -99
      imru = -99
      inac = -99
      ihydro = -99
      ichemdat = -99

      if(soln_num.lt.1.or.soln_num.ge.int(6e8)) then   ! out of range
         write(*,100)
         return
      end if
      if(soln_num.eq.1) then        ! DI water
         isoln = 1
      else if (soln_num.lt.int(1e6)) then     ! chemical data 
         ires = 0
         isoln = soln_num - 999
         ichemdat = isoln - 1
      else if (soln_num.ge.int(5e8)) then    ! mru or basin (MRU=0)
         ires = 15
         imru =  soln_num - int(5e8)
         isoln = 1+nchemdat+2*nmru*(nac+14)+2*(nhydro+1)+1+imru
      else if (soln_num.ge.int(4e8)) then   ! t+1 stream or other hydro feature
         ihydro = soln_num - int(4e8)
         isoln = 1+nchemdat+2*nmru*(nac+14)+ihydro+1
      else if (soln_num.ge.int(3e8)) then ! t0 stream or other hydro feature
         ihydro = soln_num - int(3e8)
         isoln = 1+nchemdat+2*nmru*(nac+14)+ihydro+1
      else if (soln_num.ge.206000000) then ! t+1 reservoirs after the UZ
         icut = soln_num - int(2e8)
         ires = int(icut*1.0/1e6)
         imru = int((icut*1.0 - ires*1e6)/1000.)
         isoln = 1+nchemdat-nmru+(nac*nmru)+ires*nmru+imru
      else if (soln_num.ge.205000000) then ! t+1 UZ
         icut = soln_num - int(2e8)
         ires = int(icut*1.0/1e6) 
         imru = int((icut*1.0 - ires*1e6)/1000.)
         inac = int((icut*1.0 - ires*1e6 -imru*1e3)/10)
         isoln = 1+nchemdat+nmru*4+(imru-1)*(nac+1)+inac+1
      else if (soln_num.gt.int(2e8)) then ! t+1 MRU surface res
         icut = soln_num - int(2e8)
         ires = int(icut*1.0/1e6) 
         imru = int((icut*1.0 - ires*1e6)/1000)
         isoln = 1+nchemdat+(ires-1)*nmru+imru
      else if (soln_num.ge.106000000) then  ! t0 reservoirs after the UZ
         icut = soln_num - int(1e8)
         ires = int(icut*1.0/1e6)
         imru = int((icut*1.0 - ires*1e6)/1000)
         isoln = 1+nchemdat-nmru+(nac*nmru)+ires*nmru+imru
      else if (soln_num.ge.105000000) then ! t0 UZ zone
         icut = soln_num - int(1e8)
         ires = int(icut*1.0/1e6) 
         imru = int((icut*1.0 - ires*1e6)/1000.)
         inac = int((icut*1.0 - ires*1e6 -imru*1e3)/10)
         isoln = 1+nchemdat+nmru*4+(imru-1)*(nac+1)+inac+1
      else                         ! t0 MRU surface res
         icut = soln_num - int(1e8)
         ires = int(icut*1.0/1e6) 
         imru = int((icut*1.0 - ires*1e6)/1000)
         isoln = 1+nchemdat+(ires-1)*nmru+imru
      end if
      RETURN
 100  format('Problem inverting solution number')
      END


      INTEGER FUNCTION fill_ent(n_user,soln_number,nchemdat,nmru,&
           nac,clark_segs,src_init)

      integer isoln ! function
      INTEGER src_init(:,:), n_user(:), i, index
      integer soln_number, nchemdat, nmru, nac, clark_segs
      integer ires,ichemdat,imru,inac,ihydro

      fill_ent = 1
      index = isoln(soln_number,nchemdat,nmru,nac,clark_segs,&
                 ires,ichemdat,imru,inac,ihydro)
      n_user(1) = -1
      do 10 i = 2, 11
         n_user(i) = src_init(index,i)
 10   continue
      fill_ent = 0
      RETURN
      END

      INTEGER FUNCTION checkfracs(count,solns,fracs,mixture)

      integer i,count,solns(:),mixture,snowmelt
      double precision fracs(:)
      double precision sum
      logical ET, SNOW
      integer, parameter :: a_million = 1000000
      integer, parameter :: melt = 102
      integer, parameter :: concDI = 2
      checkfracs = 1
      sum = 0
      ET = .false.
      SNOW = .false.
      do 10 i = 1, count
         if(fracs(i).lt.0) then
            if(solns(i).eq.1) ET=.true. ! check if ET solution (soln 1)
            snowmelt=solns(i)/a_million
            if(solns(i).eq.concDI.or.snowmelt.eq.melt) SNOW=.true.
!  Ignore correction if ET since ET fraction is always negative. Snowmelt also can be negative if ionic pulse simuilated
            if(.not.ET.and..not.SNOW) then
               print*,'Negative fraction encountered while ',&
                    'mixing solution ',mixture
               print*,'Value of ',fracs(i),' set to zero.'
               fracs(i)=0
            end if
         end if
         sum = sum+fracs(i)
 10   continue
      if(abs(1.0-sum).gt.1e-3) then
         print*,'Mixture ',mixture,' fractions differ from one ',&
              ' by ',abs(1.0-sum)
         return
      end if
      checkfracs = 0
      RETURN
      END

      INTEGER FUNCTION chemflag(index,chvar_lut, soln_number,nchemvar)

      INTEGER chvar_lut(:,:), i, index
      integer soln_number, nchemvar

      chemflag = 0
      do 10 i = 1, nchemvar
         if(chvar_lut(i,2).eq.soln_number) then
            chvar_lut(i,1) = index
            chemflag = 1
         end if
 10   continue
      RETURN
      END
!
! Function update_chem() updates the c_chem matrix after a mixture.
! The volume and conc values will be used to update one of the five fields
! that are the last index in c_chem (init, in, out, reaction or ET, and
! final). The argument 'imetric' dictates which of the fields is to be updated.
!
!           c_chem field(s)
!  imetric -    updated     - Comment
!    1            1           Equilibrium masses after init rxn in 1.
!    0           2,3          Transient reservoir so in=out, no rxn
!   -1           4,5          Static volume, reactions in 4 w/
!                              final mass and vols in 5.
!    3          3,4,5         Export initial concentrations to update 3;
!                              React (4) and store final mass and vols in 5.
!
! Reactions are documented in the tally table. Tally_table has one row
! for each solute species germane to the phreeq.pqi file. Note that
! this may exceed the number of solutes defined by nsolute since congruous
! weathering of minerals may produce a variety of elements. Each column
! in tally_table represents one or more phreeqc entity type. The first two
! cols are solution concentrations, in Moles per liter; the first after 
! conservative mixing, the second after reactions (so n_ent(1) always 
! equals 2). n_ent() then lists how many columns of each entity follow:
! for example, if n_ent(2)=3, that means that the solution concentrations
! will be followed by three columns of reactions. Other than solutions, 
! reactions are the only entity where a positive value adds mass to a 
! reservoir solution; all others (exchange, etc) add mass when the value
! is negative and remove mass
!
! Equilibrium phases, solid solutions, and kinetics may have more than
! one column assigned, again depending on what is listed in the .pqi file.
!
!                                       n_ent()       used in n_user
! *         solution,                   1 Solution    ! no
! *         reaction,                   2 Reaction    
! *         exchange,                   3 Exchange    
! *         surface,                    4 Surface     
! *         gas_phase,                  5 Gas_phase   
! *         equilibrium_phases,         6 Pure_phase  
! *         solid_solution,             7 Ss_phase    
! *         kinetics,                   8 Kinetics    
! *         mix,                        9 Mix         ! nope
! *         reaction_temperature        10 Temperature
! *         unknown                     11 UnKnown    ! not likely
!
!
!
! Phreeqc will keep track of a t=0 (begin with 1,3, or 5) and a t+1 
! solution composition (begins with 2,4,or6).
! However, the c_chem matrix will only populate the row pertaining to 
! the t=0 solution (begin with 1,3,or 5) since the fluxes and t+1
! compositions are tracked in the 5 fields of that solution row
!
! c_chem is the master storage table with the following 5 fields (metrics)
! describing mass (in moles) for each the nsolutes
! initial, inputs, outputs, reaction gains/losses, ending mass
! tempc and pH are recorded for the final mix (either
! static (-1) or final (3)
!
! Initial and final volumes for the index being updated 
! c_chem(indx)%vol(init) and c_chem(indx)%vol(fin)
! must be set before entering this routine in order to
! compute the reaction and final masses.
!
! Update solutes for mru and basin pseudo indices when indxm or
! indxb are nonzero respectively. Do not track volumes for mru or basin.
!
      integer function update_chem(indx,totvol,imetric,&
            conc1,pH,tempc,tally_table,n_ent,indxmru,indxbas,restype)

      USE WEBMOD_PHREEQ_MMS, ONLY: c_chem, nsolute,a_thousand,&
       maxentity, ntally_cols, ntally_rows, datetime, nstep, phr_tf, &
       indxb, indxm, indxuz, c_indx, mult, init, in, out, rxn, ET, fin
      USE WEBMOD_OBSCHEM, ONLY : sol_id
      USE WEBMOD_RESMOD, ONLY: vmix_sat, vmix_uz
      USE WEBMOD_IO, ONLY: chemout_file_unit


      implicit none

!      integer sol_id(ntally_rows,2)
      integer i,j,k,n, indx,im,imetric,n_ent(:), is  , test
      integer indxmru, indxbas, restype
      logical react,step1,notstream,snow
      double precision vol,totvol, conc1(:),pH,tempc,tally_table(:,:)


      data step1/.true./
      save step1

      update_chem = 1

!
! Add up the number of tally table columns the first time through,
! ignoring solution concentrations that are always in column one.
! nally_cols is static for each run. Assign a multiplier of 1 for mass
! produced by reactions (n_ent(2)) and negatives for all other 
! entities (exchange, surface, gas, etx)
!
      if(step1) then
         do i = 1,ntally_cols
            mult(i) =1D0
            if(i.gt.2+n_ent(2)) mult(i) = -1D0
! Debug
! /Debug
         end do
         if(phr_tf) write(26,'(A)')'nstep yr mo dy index In_Out metric indxmru indxbas '// &
           'totvol vol_init vol_in vol_out vol_rxn vol_fin Cl_init Cl_in Cl_out Cl_rxn Cl_fin'
         step1=.false.
      end if  ! step1

! Debug
      test =  0 ! Snowpack
     
      if (phr_tf) then
! for nsolute = 12
!        write(26,123) nstep,(datetime(i),i=1,3),indx, &
!           ' in ',imetric, indxmru, indxbas, totvol,(c_chem(indx)%vol(i),i=1,5), &
!           ((c_chem(indx)%M(i,j),j=1,5),i=8,12)
! for nsolute = 1
        write(26,123) nstep,(datetime(i),i=1,3),indx, &
           ' in ',imetric, indxmru, indxbas, totvol,(c_chem(indx)%vol(i),i=1,5), &
           (c_chem(indx)%M(1,j),j=1,5)
      end if
!
! / Debug
!
!
!
! Check if stream
!
      notstream = .TRUE.
      if(restype.eq.99) notstream = .FALSE.
!
! Check if snowpack, which has special treatment of melt and final pack to simulate incongruous melting
!
      snow = .FALSE.
      if(restype.eq.2) snow = .TRUE.
!
!     Update only one slot based on concentrations
!
      j=1
      im = imetric
      if(imetric.eq.0) then ! unless imetric = 0, then update both inputs and outputs
        im = 2
        j=2
      end if
!
! This loop duplicates imports and exports for transient reservoirs (imetric=0).
! Otherwise only one pass through.
      do 5 i=1,j
      
! Debug
!        write(chemout_file_unit,1000)c_indx(indx,1), indx, indxmru,
!     $     indxbas, totvol, imetric, im, restype 
! 1000 format(I10, 3I5,f12.2,3I5)
! \Debug
!     
!     Update volumes,
!
         if(im.eq.1) then
           vol = c_chem(indx)%vol(init)
         else if(im.eq.-1) then
           vol = c_chem(indx)%vol(fin)   ! check if this should in 'in' instead of 'fin'
         else
           c_chem(indx)%vol(im) =&
               c_chem(indx)%vol(im)+totvol
         end if
!     
!     masses and isotopes,
!
         is=0  ! counter for number of isotopes tracked
         do 10 n = 1,nsolute
           if(im.eq.1) then ! initial masses use the reaction concentrations in the 2nd column of the tally table
             if(sol_id(n)%iso) then  ! record inital isotopic composition
               is=is+1
               c_chem(indx)%delta(n,1)=conc1(nsolute+is)
             end if
             k=sol_id(n)%tally
             c_chem(indx)%M(n,1) = tally_table(k,2)*vol*a_thousand
           else if(im.ne.-1) then
             c_chem(indx)%M(n,im)= c_chem(indx)%M(n,im)+ & !Imports and exports use conservative mixing&
              conc1(n)*totvol*a_thousand !This should cover snowpack chem on days of no snowpack.
             if(sol_id(n)%iso) then ! track input and output deltas which shoul only occur on a single input or output.
               is=is+1
               c_chem(indx)%delta(n,im)=conc1(nsolute+is)
             end if
           end if
           if(restype.eq.5) & ! sum riparian or upland input as delivered in indxuz
            c_chem(indxuz)%M(n,im)=c_chem(indxuz)%M(n,im)+ &
                                   c_chem(indx)%M(n,im)
           if(indxmru.ne.0) & ! sum mru inputs and outputs&
            c_chem(indxmru)%M(n,im)=c_chem(indxmru)%M(n,im)+ &
                                    c_chem(indx)%M(n,im)
           if(indxbas.ne.0)  & ! sum basin inputs and outputs&
            c_chem(indxbas)%M(n,im)=c_chem(indxbas)%M(n,im)+ &
                                    c_chem(indx)%M(n,im)
 10      continue
!
!     and temp and pH.
!
! Record temp and pH if imetric > 0
         if(im.gt.0) then
           c_chem(indx)%Temp = tempc
           c_chem(indx)%pH = ph
         endif
         if(indxmru.ne.0) then
           c_chem(indxmru)%Temp = tempc
           c_chem(indxmru)%pH = ph
         endif
         if(indxbas.ne.0) then
           c_chem(indxbas)%Temp = tempc
           c_chem(indxbas)%pH = ph
         endif
!
! Duplicate in and out values (2&3) if imetric = 0   
!
         if(imetric.eq.0) im = im + 1
 5    continue
!     
!     Compute reaction totals for imetric = -1 or 3 (static or last export). Snow always has a 3 and a -1 to accomodate
!     incongruent melting so only count reactions during last -1 step. All other reservoirs have either a 3 OR a -1, not both.
!
!     Accumulate reactions, temp, and pH for unchanged reservoirs and final exports
!
      if(imetric.eq.-1.or.(imetric.eq.3.and..not.snow)) then
          vol = c_chem(indx)%vol(fin)  ! reactions take place in final volumes
          if(snow.and.imetric.eq.3)  vol = c_chem(indx)%vol(out) ! not needed since no_rxn in melt step
          react = .FALSE.
          if(indxmru.ne.0) then ! record temp, pH, and del values of mru exports
             c_chem(indxmru)%Temp = tempc
             c_chem(indxmru)%pH = ph
             is=0            
             do n = 1,nsolute
               if(sol_id(n)%iso) then
                 is=is+1
                 c_chem(indxmru)%delta(n,out)=conc1(nsolute+is)
               end if
             end do
          endif
          if(indxbas.ne.0) then ! record temp, pH, and del values of basin exports
             c_chem(indxbas)%Temp = tempc
             c_chem(indxbas)%pH = ph
             is=0            
             do n = 1,nsolute
               if(sol_id(n)%iso) then
                 is=is+1
                 c_chem(indxbas)%delta(n,out)=conc1(nsolute+is)
               end if
             end do
          endif
          do 12 n = 1, nsolute
            k=sol_id(n)%tally
!            c_chem(indx)%M(n,rxn) = 0D0
            do 11 j = 3, ntally_cols  ! first two cols are conc. of conservative and reactive solutions
              c_chem(indx)%M(n,rxn) = c_chem(indx)%M(n,rxn)+&
                tally_table(k,j)*mult(j)*vol*a_thousand
              if(abs(c_chem(indx)%M(n,rxn)).gt.0.01) react = .TRUE.
   11       continue
!  sum reactions for basin, mru, and uz; mru only, if not stream. Note that the
!  indices are those shared through the PHREEQMMS_MOD. Reactions do not occur
!  in transient mixes (mrures>9).
           c_chem(indxb)%M(n,rxn) = c_chem(indxb)%M(n,rxn)+&
             c_chem(indx)%M(n,rxn)
           if(notstream) c_chem(indxm)%M(n,rxn)=c_chem(indxm)%M(n,rxn)+ &
              c_chem(indx)%M(n,rxn)
           if(restype.eq.5) & ! sum uz inputs&
            c_chem(indxuz)%M(n,rxn)=c_chem(indxuz)%M(n,rxn)+ &
              c_chem(indx)%M(n,rxn)
   12     continue
      end if
!
!     Use final concentration to establish final masses and record isotopes
!     This will avoid drift errors in final masses from additions.
!     im for final mix is either a 3 or a -1. The last one will 
!     stick (important for the incongruous melting)
!
      if(im.ne.1.and.im.eq.3.or.im.eq.-1) then 
        is=0
        do 15 n = 1, nsolute
          if(sol_id(n)%iso) then
            is=is+1
            c_chem(indx)%delta(n,fin)=conc1(nsolute+(3*is-2))
          end if
          k=sol_id(n)%tally
          c_chem(indx)%M(n,fin) = tally_table(k,2)*vol*a_thousand
   15   continue
      end if  
!  debug
!      write(25,121)nstep,indx,indxb, imetric,
!     $     (c_chem(indxb)%M(5,im),c_chem(indxb)%vol(im),im = 1,5)
! 121  format(I4, 2I10, i4, 10F12.1)
!      if (indx.gt.200) write(25,122)nstep,indx,indxb, imetric,
!     $     (c_chem(indx)%M(5,im),c_chem(indx)%vol(im),im = 1,5)
! 122  format(I4, 2I10, i4, 10F12.1)

      
      if (phr_tf) then
! for nsolute = 12
!        write(26,123) nstep,(datetime(i),i=1,3),indx, &
!           ' out ',imetric, indxmru, indxbas, totvol,(c_chem(indx)%vol(i),i=1,5), &
!           ((c_chem(indx)%M(i,j),j=1,5),i=8,12)
!           i=1
! for nsolute = 1
        write(26,123) nstep,(datetime(i),i=1,3),indx, &
           ' out ',imetric, indxmru, indxbas, totvol,(c_chem(indx)%vol(i),i=1,5), &
           (c_chem(indx)%M(1,j),j=1,5)
           i=1
      end if
! // debug

      update_chem = 0
 123  format(5I5,A,3I5,1X,31E14.6)

      RETURN
      END

!
! Function to convert and assign and manipulate values in the
! c_chem matrix to one of the nchemvar variables (10 variables
! for now).
!
!
      integer function chem2var()
   
      USE WEBMOD_PHREEQ_MMS

      integer ivar, isol, c_unit
      logical been_warned_iso, been_warned_vol, been_warned_area
      
! local variable
      double precision ch_var_tmp(nchemvar,nsolute)
      double precision ch_var_tmp_m3(nchemvar)
      double precision ch_var_tmp_tempc(nchemvar)
      double precision ch_var_tmp_pH(nchemvar)
      double precision moles

      data been_warned_iso /.false./
      data been_warned_vol /.false./
      data been_warned_area /.false./

      save been_warned_iso, been_warned_vol, been_warned_area

      chem2var = 1
!
!  Computations based on metric and unit
!
      do 1 ivar = 1,nchemvar
         indx = chvar_lut(ivar,1)
         imet = chvar_lut(ivar,3)
         c_unit = chvar_lut(ivar,4)
         unit_type = int((c_unit-1)/3)+1 ! 1,mass; 2,load; 3,std conc; 4,user conc'; 5,permil
!  Warn about limitations of permil units
         if(unit_type.eq.5) then
           if(.not.been_warned_iso) then
             print*,'chvar ',ivar,' indicates units of permil. ',&
             'A c_metric of 6 (net) is not valid. ',&
             'A c_metric of 4 (rxn) will report the delta of 18O or D ',&
             'in the evapotranspiration from that reservoir. ',&
             'If the selected solute is not an isotope ',&
             'the value will remain at zero. All subsequent ',&
             'warnings about isotopes for other chvars will be ignored.'
             been_warned_iso = .true.
           end if
         end if
!     Volumes first
         if(imet.eq.6) then     ! Net, so final-init
            vol = c_chem(indx)%vol(fin)-c_chem(indx)%vol(init)
         else                   ! All other metrics point in the imet direction
            vol = c_chem(indx)%vol(imet)  ! this is ET if imet=4, no concentrations for solutes permitted
         end if
         ch_var_tmp_m3(ivar)=vol
!     then mass
         do 10 isol = 1,nsolute
          if(imet.eq.6) then  ! Net, so final-init
             moles = c_chem(indx)%M(isol,fin)-c_chem(indx)%M(isol,init)
          else                ! All other metrics point in the imet direction
             moles = c_chem(indx)%M(isol,imet)
          end if
!     Finally conversions
          if(unit_type.eq.1) then ! mass only (mg, meq, or mmol)
             ch_var_tmp(ivar,isol)=moles*chvar_conv(ivar,isol)
          else if (unit_type.eq.2) then ! load so divide mass by area
            if(chvar_conv(ivar,nsolute+1).gt.0.0) then
              ch_var_tmp(ivar,isol)= moles*chvar_conv(ivar,isol)/&
                                       chvar_conv(ivar,nsolute+1)
            else
               if(.not.been_warned_area) then
                 print*,'chvar ',ivar,' indicates a load ',&
                  'unit and there is zero area. This may be ',&
                  'occuring because you are requesting a ',&
                  'load for canopy or snowpack on a day with there ',&
                  'is none. The requested concentration on such ',&
                  'days will be set to zero. All subsequent ',&
                  'warnings of this condition will be ignored.'
                 been_warned_area = .true.
               end if
              ch_var_tmp(ivar,isol)= 0.0
            end if
          else if (unit_type.lt.5) then ! standard or user-defined concentration so divide mass by volume
            if(vol.gt.0.0) then
               ch_var_tmp(ivar,isol)=&
                    moles*chvar_conv(ivar,isol)/vol/a_thousand
            else
              if(.not.been_warned_vol) then
                 print*,'chvar ',ivar,' indicates a concentration ',&
                  'unit and there is zero volume. This may be ',&
                  'occuring because you are requesting a ',&
                  'concentration for a day with zero flux (vol=0). ',&
                  'The requested concentration on such days ',&
                  'will be set to zero. All subsequent warnings ',&
                  'of this condition will be ignored.'
                 been_warned_vol = .true.
              end if
              ch_var_tmp(ivar,isol)= 0.0
            end if
          else if(c_unit.eq.13) then   ! permil isotopes. Change to unit_type 4 if more iso units added (pmc, etc)
            ch_var_tmp(ivar,isol)= c_chem(indx)%delta(isol,imet)
          else 
            print*,'chvar ',ivar,' has invalid units indicated.'
          end if
 10      continue

! Update_chem will overwrite temp and pH into all 5 metrics (in, out) on
! each mix. Final value should be a export and react or static reaction (update_chem 3 or -1)

      ch_var_tmp_tempc(ivar)=c_chem(indx)%Temp(imet)
      ch_var_tmp_pH(ivar)=c_chem(indx)%pH(imet)

 1    continue
!
! Transfer values from local variables
!
      ch_var_01_m3 = ch_var_tmp_m3(1)
      ch_var_02_m3 = ch_var_tmp_m3(2)
      ch_var_03_m3 = ch_var_tmp_m3(3)
      ch_var_04_m3 = ch_var_tmp_m3(4)
      ch_var_05_m3 = ch_var_tmp_m3(5)
      ch_var_06_m3 = ch_var_tmp_m3(6)
      ch_var_07_m3 = ch_var_tmp_m3(7)
      ch_var_08_m3 = ch_var_tmp_m3(8)
      ch_var_09_m3 = ch_var_tmp_m3(9)
      ch_var_10_m3 = ch_var_tmp_m3(10)
      ch_var_01_tempc = ch_var_tmp_tempc(1)
      ch_var_02_tempc = ch_var_tmp_tempc(2)
      ch_var_03_tempc = ch_var_tmp_tempc(3)
      ch_var_04_tempc = ch_var_tmp_tempc(4)
      ch_var_05_tempc = ch_var_tmp_tempc(5)
      ch_var_06_tempc = ch_var_tmp_tempc(6)
      ch_var_07_tempc = ch_var_tmp_tempc(7)
      ch_var_08_tempc = ch_var_tmp_tempc(8)
      ch_var_09_tempc = ch_var_tmp_tempc(9)
      ch_var_10_tempc = ch_var_tmp_tempc(10)
      ch_var_01_pH = ch_var_tmp_pH(1)
      ch_var_02_pH = ch_var_tmp_pH(2)
      ch_var_03_pH = ch_var_tmp_pH(3)
      ch_var_04_pH = ch_var_tmp_pH(4)
      ch_var_05_pH = ch_var_tmp_pH(5)
      ch_var_06_pH = ch_var_tmp_pH(6)
      ch_var_07_pH = ch_var_tmp_pH(7)
      ch_var_08_pH = ch_var_tmp_pH(8)
      ch_var_09_pH = ch_var_tmp_pH(9)
      ch_var_10_pH = ch_var_tmp_pH(10)
      do 2 isol = 1,nsolute
      ch_var_01_sol(isol) = ch_var_tmp(1,isol)
      ch_var_02_sol(isol) = ch_var_tmp(2,isol)
      ch_var_03_sol(isol) = ch_var_tmp(3,isol)
      ch_var_04_sol(isol) = ch_var_tmp(4,isol)
      ch_var_05_sol(isol) = ch_var_tmp(5,isol)
      ch_var_06_sol(isol) = ch_var_tmp(6,isol)
      ch_var_07_sol(isol) = ch_var_tmp(7,isol)
      ch_var_08_sol(isol) = ch_var_tmp(8,isol)
      
      ch_var_09_sol(isol) = ch_var_tmp(9,isol)
      ch_var_10_sol(isol) = ch_var_tmp(10,isol)
 2    continue

      chem2var = 0

      RETURN
      END FUNCTION chem2var

!
! Compute delta of residual liquid then back calculate delta of evaporated water (solution 1) and account for those
! isotopes leaving the reservoir, the mru, and the basin.
! phase: 0, snow; 1, water; 3, transpiration (no fractionation 
!        in deep UZ or transpiration from dry canopy)
! If phase = 3, reset to 1 to find right evaporation index in conc array
!
      integer function fractionate(phase,resindx,evapvol,resv,ison,rh,tempevap)
   
      USE WEBMOD_PHREEQ_MMS, only : aline, one, a_thousand, rxn, iso_n,&
             almost_one, fill_factor, phr_tf, no_rxn, n_ent,&
             rxnmols, tempc, pH, ph_final, tsec, tally_table, ires,&
             ntally_rows, ntally_cols, conc, indxm, indxb, indxuz,&
             ln_10, c_chem,ID, SETOUTPUTFILEON, SETERRORFILEON, &
             SETLOGFILEON, SETSELECTEDOUTPUTFILEON, res_D_permil, evap_D_permil, &
             res_18O_permil,evap_18O_permil, delta_res_permil

!      USE WEBMOD_PHREEQ_MMS
      USE WEBMOD_RESMOD, ONLY : vmin_canopy
      USE WEBMOD_OBSCHEM, ONLY : phq_lut, n_iso, sol_id, iso_list,&
            nsolute, sol_name

      implicit none
      
      integer :: i, k, resindx, solns(1), iresult
      integer :: id_len, phase   ! phase = 1 for water and 2 for snow
      integer, external ::  length
      integer, external ::  phr_mix, fill_ent
!      integer, external ::  update_chem, accumulateline, run
      integer, external  ::  accumulateline, runaccumulated
! local variables use resv in argument list to avoid changing totvol in the calling routine
      double precision :: evapvol, resvol, resv, evap_frac, rh, tempevap, tot_16O, tot_H
      double precision :: log_a, eps_eq, eps_diff,ison,fracs(1)
      double precision :: delta_res0, delta_res, delta_evap_permil,evap_D_M,evap_18O_M
!      integer :: test
!
!      allocate(conc1(size(conc)))
!      conc1 = conc
!      
!      iresult = test(conc)
      fractionate = 1
      res_D_permil = -1000D0
      evap_D_permil = -1000D0
      res_18O_permil = -1000D0
      evap_18O_permil = -1000D0
      tot_H = tally_table(2,2)  ! total moles of protium in reservoir (normalized to near 1 kg)
      tot_16O = tally_table(3,2) ! total moles of 16O in reservoir (normalized to near 1 kg)
      resvol = resv
! if the phase = 3 (deep UZ), no fractiontion.
      if(phase.eq.3) then
        evap_frac = 1.0
      else
        evap_frac = evapvol/resvol  ! evaporation fraction for kinetic (diffusive) fractionation
      endif
! if phase = 3 (no fractionation) set index to 1(water) to find right offset
! check delta of mixture and compute fractionation of evaporation
! Compute delta of residual liquid then back calculate delta of evaporated water (solution 1)
! if n_iso, the number of isotopes, =0, this section is skipped and pure water is evaporated
!      evap_frac = evapvol/resvol  ! evaporation fraction for kinetic (diffusive) fractionation
      if(abs(1.0-evap_frac).lt.0.0001) evap_frac=1.0D0
      i = 1
      iresult = accumulateline(id, "solution 1")  ! solution 1 is always evaporative water
      do while(i.le.n_iso)
       delta_res0 = conc(nsolute+3*i-2)/a_thousand    ! delta of reservoir water before evaporation, permil converted to numeric value 
       if(evap_frac.lt.almost_one.and.phase.ne.3) then
         log_a = conc(nsolute+3*i-phase)    ! common log of alpha of 1st isotope water=1(vapor-liquid), snow=0 (vapor-ice), returned from phreeq mix based on temperature
         eps_eq = ln_10 * log_a  ! equilibrium fractionation (small number i.e. permil/1000)
         eps_diff = (1.0-rh)*ison*&
                    phq_lut(sol_id(iso_list(i))%phq)%isofrac_max  ! diffusive fractionation (small number i.e. permil/1000)
         delta_res = (delta_res0+1)*(1-evap_frac)**&
                            (eps_eq+eps_diff)-1
         delta_evap_permil = (delta_res0-delta_res*(1-&
                              evap_frac))/evap_frac*a_thousand
         delta_res_permil = delta_res*a_thousand
       else
         delta_evap_permil = delta_res0*a_thousand
         delta_res_permil = delta_evap_permil
       end if
       id_len = length(sol_name(iso_list(i)))
       if(sol_name(iso_list(i))(1:id_len).eq."D") then
          res_D_permil = delta_res_permil
          evap_D_permil = delta_evap_permil
! use to double check mix below
          evap_D_M =  (evap_D_permil/a_thousand+one)* phq_lut(sol_id(iso_list(i))%phq)%isoratio * tot_H  ! diffusive fractionation (small number i.e. permil/1000)
       elseif(sol_name(iso_list(i))(1:id_len).eq."[18O]") then
          res_18O_permil = delta_res_permil
          evap_18O_permil = delta_evap_permil
! use to double check mix below
          evap_18O_M =  (evap_18O_permil/a_thousand+one)* phq_lut(sol_id(iso_list(i))%phq)%isoratio * tot_16O  ! diffusive fractionation (small number i.e. permil/1000)
       else
          print *,"fractionated isotope not D or [18O]. Run stopped"
          return
       endif
       write(aline,*)sol_name(iso_list(i)&
                       )(1:id_len)//" ",delta_evap_permil
       iresult = accumulateline(id,aline)
       i=i+1
      end do

! Assign depleted deltas to evaporative water and track evaporated isotopes
!           run_(int *output_on, int *error_on, int *log_on, int *selected_on)
!         iresult = run(1,1,0,1)   ! assigns depleted deltas to evap using lines accumulated above
      iresult = SetOutputFileOn(ID,phr_tf)
      iresult = SetErrorFileOn(ID,phr_tf)
      iresult = SetLogFileOn(ID,phr_tf)
      iresult = SetSelectedOutputFileOn(ID,phr_tf)
!
! Set temperature, pH and mix
!
      write(aline,*)'-temp ', tempevap
      iresult = AccumulateLine(id, aline)
      write(aline,*)'-ph 7 charge'  ! reduces negative alkalinity during acid rain deposition
      iresult = AccumulateLine(id, aline)
      iresult = runaccumulated(ID)   ! assigns depleted deltas and temperature to evap using lines accumulated above
      solns(1) =  one  ! solnnum(0,0,0,0,0,0,0) , evap water with correct deltas
      fracs(1) = 1.0
      iresult = phr_mix(ID,1, solns, fracs,  solns(1),&
                 fill_factor, solns(1), conc, phr_tf, no_rxn,&
                 rxnmols,tempc,ph,ph_final,tsec,tally_table,ntally_rows,&
                 ntally_cols)
      IF (iresult.NE.0) THEN
         PRINT *, 'Errors during phr_mix:'
         CALL OutputErrorString(id)
         STOP
      ENDIF
!
! Account for D and 18O isotopes removed during
! evaporation as negative reactive contributions
! so that the lost mass is subtracted from the initial
! mass. The loss from the basin, MRU, and gen UZ reservoir are
! traced in the update_chem routine

        do 15 i = 1, nsolute
          if(sol_id(i)%iso) then
            k=sol_id(i)%tally
            c_chem(resindx)%M(i,rxn) = c_chem(resindx)%M(i,rxn) - tally_table(k,2)*evapvol*a_thousand
!            c_chem(indxm)%M(i,rxn) = c_chem(indxm)%M(i,rxn) - tally_table(k,2)*evapvol*a_thousand
!            c_chem(indxb)%M(i,rxn) = c_chem(indxb)%M(i,rxn) - tally_table(k,2)*evapvol*a_thousand
!            if(ires.eq.5) & ! sum uz inputs&
!            c_chem(indxuz)%M(i,rxn)=c_chem(indxuz)%M(i,rxn) - tally_table(k,2)*evapvol*a_thousand
          end if
   15   continue
      
!
!      iresult = test(resindx, evapvol, i, conc)
!     $, pH,tempc,
!     $                         tally_table,n_ent,indxm,indxb,ires)
!      iresult = update_chem(resindx,evapvol,4,conc,pH,tempc,&
!                               tally_table,n_ent,indxm,indxb,ires)
!      IF (iresult.NE.0) THEN
!        PRINT *, 'Errors updating mole matrix:'
!        STOP
!      ENDIF
!
      fractionate = 0
 120  FORMAT(A,1PG15.7E2)
      RETURN
      END FUNCTION fractionate
!
!      integer function test(conc1)
!!     $, pH,tempc,
!!     $                         tally_table,n_ent,indxm,indxb,ires)
!      double precision, intent (inout) :: conc1(:) ! , e !, pH, tempc, tally_table(:,:)
!      integer :: i, j !, n_ent(:), indxm, indxb, ires
!      test = 0
!      end function test
!
!
! Create pure water with specific isotopic signature at specified temperature
! after sublimation to create snowmelt with ionic pulse and lighter isotopes.
!
      integer function reset_DI(ID,tempc)

      USE WEBMOD_PHREEQ_MMS, only: delta_D,delta_18O

      implicit none
      integer, external  ::  accumulateline, runaccumulated
! local variables use resv in argument list to avoid changing totvol in the calling routine
      integer :: ID, iresult
      real :: tempc
      character*256 :: aline
!
      reset_DI = 1
!
! solution 2 will be used for concentrating melt and producing propoer deltas. DI, or solution 1 is mass corrected to more accurately track basin exports
!
      iresult = accumulateline(ID, "solution 2")
      WRITE (aline,*)'-temp ', tempc
      if(delta_D.ne.-1000) then
          write(aline,*)"D ",delta_D
          iresult = AccumulateLine(ID, aline)
      elseif(delta_18O.ne.-1000) then
          write(aline,*)"[18O] ",delta_18O
          iresult = AccumulateLine(ID, aline)
      endif
      reset_DI = runaccumulated(ID)
      RETURN
    END FUNCTION reset_DI
    
!     program to converge on wet-bulb temp, given dry-bulb
!     and relative humidity - ported by RW on 23 May 2013 from tw.f 
!     written by Dave Stannard
 
      real function wetbulb (pres, ta, rh)
      implicit double precision (a-h,o-z)
      real pres, ta, rh
      wetbulb = -999.
!     site pressure in kilopascals
!      pres=101.3

!     Coefficients for the lowe (1976) equations for saturation 
!     vapor pressure, and the slope of that curve

!     sat. vap. press. curve
      a0=6.107799961
      a1=4.436518521e-1
      a2=1.428945805e-2
      a3=2.650648471e-4
      a4=3.031240396e-6
      a5=2.034080948e-8
      a6=6.136820929e-11

!     slope of sat. vap. press. curve
      b0=4.438099984e-1
      b1=2.857002636e-2
      b2=7.938054040e-4
      b3=1.215215065e-5
      b4=1.036561403e-7
      b5=3.532421810e-10
      b6=-7.090244804e-13

!     sat. vap. press. over ice
      a0i=6.109177956
      a1i=5.03469897e-1
      a2i=1.886013408e-2
      a3i=4.176223716e-4
      a4i=5.824720280e-6
      a5i=4.838803174e-8
      a6i=1.838826904e-10

!     slope sat. vap. press. curve over ice
      b0i=5.030305237e-1
      b1i=3.773255020e-2
      b2i=1.267995369e-3
      b3i=2.477563108e-5
      b4i=3.005693132e-7
      b5i=2.158542548e-9
      b6i=7.131097725e-12

!     Read in dry bulb temp and relative humidity (0-100)
!10    read (1,*,end=999)ta,rh

!     For first tw estimate, and to compute actual VP,
!     compute es (saturated vapor pressure) and s (slope of saturated
!     vapor pressure) at dry bulb
      if (ta .gt. 0.)then
       es=(a0+ta*(a1+ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*a6))))))/10.
       s=(b0+ta*(b1+ta*(b2+ta*(b3+ta*(b4+ta*(b5+ta*b6))))))/10.
      else
       es=(a0i+ta*(a1i+ta*(a2i+ta*(a3i+ta*(a4i+ta*(a5i+ta*a6i))))))/10.
       s=(b0i+ta*(b1i+ta*(b2i+ta*(b3i+ta*(b4i+ta*(b5i+ta*b6i))))))/10.
      endif

!     Compute actual VP (emeas) to converge upon
      emeas=es*rh/100.0d0

!     For first tw estimate, use gam, es, and s at dry bulb
      gam=pres*.00066*(1.+.00115*ta)   ! psychrometric constant
      tw=ta-(es-emeas)/(s+gam)
      etemp=emeas+(ta-tw)*gam
 !     write (6,*)'   Dry bulb       RH       Wet bulb     Mod. VP    Meas. VP'

!     For iterations, use gam, es, and s at wet bulb
20    if (tw .gt. 0.0)then
       es=(a0+tw*(a1+tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*a6))))))/10.
       s=(b0+tw*(b1+tw*(b2+tw*(b3+tw*(b4+tw*(b5+tw*b6))))))/10.
      else
       es=(a0i+tw*(a1i+tw*(a2i+tw*(a3i+tw*(a4i+tw*(a5i+tw*a6i))))))/10.
       s=(b0i+tw*(b1i+tw*(b2i+tw*(b3i+tw*(b4i+tw*(b5i+tw*b6i))))))/10.
      endif
 
      gam=pres*.00066*(1.+.00115*tw)
      tw=tw-(es-etemp)/(s+gam)
      etemp=emeas+(ta-tw)*gam
      emod=es-gam*(ta-tw)
 
!      write (6,2000)ta,rh,tw,emod,emeas,etemp

!     the following is the convergence criterion in degrees C
      if (abs((emod-emeas)/(s+gam)) .gt. .0001)go to 20
 
!      write (2,2000)ta,rh,tw,emod,emeas
      wetbulb = tw
      return
 !      go to 10
!999   stop
!2000  format(6f12.6)
    end function wetbulb 
    
double precision function webmod_callback(x1, x2, str)
    use WEBMOD_POTET
    use WEBMOD_PHREEQ_MMS
    double precision x1, x2
    character (*) str
    if (str .eq. "transp_on") then
        webmod_callback = dble(transp_on(iphrq_mru))
        return 
    else if (str .eq. "year") then
        webmod_callback = dble(datetime(1))
        return
    else if (str .eq. "month") then
        webmod_callback = dble(datetime(2))
        return
    else if (str .eq. "day") then
        webmod_callback = dble(datetime(3))
        return
    else if (str .eq. "doy") then
        webmod_callback = 1.1
        return
    endif
    webmod_callback = 0
    return 
end function webmod_callback



      
