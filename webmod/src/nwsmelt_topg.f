!***********************************************************************
!     nwsmelt.f:   snowmelt routine according to the National Weather
!                  Service River Forecast System - Snow Accumulation
!                  and Ablation Model (Eric A. Anderson, 1973)
!                  Rewritten and slightly modified for the use in
!                  combination with the Enns-Modell by M. Fuchs
!                  May 1997 , modified to more match newest routine
!                  of Anderson. This is done in 1998 3/9-12 by jvccaro
!**********************************************************************
c 8 Sept 2003 - Added RCS version control - RMTW
c
c Modifications by Rick Webb 2000 August-September
c
c Added a few lines to handle daily values with mixed precip and no
c snowpack.
c
c Removed snow undercatch factor, SCF. Undercatch to be calculated in 
c the precip module using parameter snow_adj(nmru, nmonths). Need to uncomment
c out selected lines and re-add SCF to the init and run argument lists to reactivate.
c
c The local PXTEMP variable will be assigned the value of tmax_allsnow_c declared in precip.
c
c 7/9/02 - Added basin average snow and precip values basin_net_snow and basin_net_rain RMTW
c          Removed partition of net_rain and net_snow from intercept to avoid more problems
c          on mixed rain/snow days. Just start with pmru as net_ppt.
c          Reinitialize PPTMIX_NOPACK(I1) = 0 at each MRU loop.
c
c June 2003 - Corrected TAK and TAIR variables that were not initialized nor saved (in MELTRATE).
c             Corrected aerial extent of snow equiv (SB = SBI() * 25.4)(in AESNOW). Now depletion
c             curves work as expected. - RMTW
c
c 06may10 - Port to Fortran 90 with module and dynamic memory
c
c 01mar12 - Change LIQW, SBI, SBWSI, SBAESCI, and LAGROI to internal variables to simplify initialization
c         - Add AESCI for initial aerial extent, AESC, derived from WEI, ACUMX, SI, and ADC
c
c***********************************************************************
c 
      MODULE WEBMOD_SNOW
      IMPLICIT NONE
      include 'fmodules.inc'

C   Dimensions
      integer, save :: nmru, nexlag, ndepl, ndeplval
      
C   Declared Variables
      integer, save, allocatable:: PPTMIX_NOPACK(:)
      real, save :: BASIN_SNOWMELT, PSOILBASIN, COVERBASIN,BASIN_PWEQV
      real, save :: BASIN_SNOWEVAP, basin_net_snow, basin_net_rain
      real, save, allocatable :: SNOW_EVAP(:),CHNGPCKMRU(:)
      real, save, allocatable :: SNOWMELT(:), PSOILMRU(:)
      real, save, allocatable :: TINDXI(:),PKWATER_EQUIV(:)
      real, save, allocatable :: SNOWCOV_AREA(:)
      
C   Declared Parameters
      real, save :: BASAREA, perAREA, ALAT, tmax_allsnow_c
      integer, save, allocatable:: MRUDEPL(:)
      real, save, allocatable :: ACUMX(:),  NEGHSI(:)
c      real, save, allocatable :: LIQWI(:),  SBWSI(:)
c      real, save, allocatable :: SBI(:),    SBAESCI(:)
      real, save, allocatable :: WEI(:), SUBRATE(:)
      real, save, allocatable :: UADJ(:),   DAYGM(:)
      real, save, allocatable :: MFMAX(:),  MFMIN(:)
      real, save, allocatable :: NMF(:),  MBASE(:)
      real, save, allocatable :: PLWHC(:), TIPM(:), SI(:)
      real, save, allocatable :: ADC(:)
      real, save, allocatable :: ELEV(:), MRUAREA(:)
C
C    Internal variables that were previously parameters but
C    assigned here based on WEI, ACUMX, and SI to simplify initial conditions
C
      real, save, allocatable :: LIQWI(:), SBWSI(:)
      real, save, allocatable :: SBI(:), SBAESCI(:), AESCI(:)
      real, save, allocatable :: LAGROI(:,:), LAGRO(:)
      real, save, allocatable :: STOREI(:), TMXPREI(:)
      
C   Undeclared Static Variables gotten from from other modules
C   Note that vars can be renamed  (i.e. pmru is net_ppt from ppt routine)
C
!  WEATHER DATA: previous public variables also used for model
      real, save, allocatable :: pmru(:), net_rain(:), Net_snow(:)
      real, save, allocatable :: tmxmru(:), tmnmru(:)
! pxtemp() = tmax__allsnow_c
      real, save, allocatable :: pxtemp(:)
    
      END MODULE WEBMOD_SNOW
c
!***********************************************************************
!
!     Main nwsmelt routine
!     ------------
      integer function nwsmelt_topg(arg)

      character*(*) arg
      CHARACTER*256 SVN_ID

      integer nwsmdecl, nwsminit, nwsmrun

      save SVN_ID

      SVN_ID = 
     $     '$Id: nwsmelt_topg.f 41 2008-07-24 16:23:18Z rmwebb $ '

      if(arg.eq.'declare') then
        nwsmelt_topg = nwsmdecl()
      else if(arg.eq.'initialize') then
        nwsmelt_topg = nwsminit()
      else if(arg.eq.'run') then
        nwsmelt_topg = nwsmrun()
      end if

      return
      end

!***********************************************************************
! 
!     nwsmdecl - DECLARE variables and parameters for snowmelt computations
!     --------
      integer function nwsmdecl()

      USE WEBMOD_SNOW

      nwsmdecl = 1

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nexlag = getdim('nexlag')
        if ( nexlag.eq.-1 ) return
      ndepl = getdim('ndepl')    ! number of depletion curves
        if ( ndepl.eq.-1 ) return
      ndeplval = getdim('ndeplval')    ! each depletion curve has 11 values
        if ( ndeplval.eq.-1 ) return

      ALLOCATE (SNOWMELT(Nmru))
      if(declvar('snow', 'snowmelt', 'nmru', nmru, 'real',
     &    'mean value of snowmelt in each mru',
     &    'inches',
     &   SNOWMELT).ne.0) return     

      ALLOCATE (PSOILMRU(Nmru))
      if(declvar('snow', 'psoilmru', 'nmru', nmru, 'real',
     &    'precipitation, irrigation, and thrufall on non-'//
     &    'snow-covered areas in each mru', 'inches',
     &   PSOILMRU).ne.0) return   
 
      ALLOCATE (SNOWCOV_AREA(Nmru))
      if(declvar('snow', 'snowcov_area', 'nmru', nmru, 'real',
     &    'percentage of snowcovered ground in each mru',
     &    'none',
     &   SNOWCOV_AREA).ne.0) return

      if(declvar('snow', 'basin_snowmelt', 'one', 1, 'real',
     &    'Basin area-weighted average of snowmelt',
     &    'inches',
     &   BASIN_SNOWMELT).ne.0) return
     
      if(declvar('snow', 'basin_net_snow', 'one', 1, 'real',
     &    'Basin area-weighted average snowfall minus interception',
     &    'inches',
     &   basin_net_snow).ne.0) return
     
      if(declvar('snow', 'basin_net_rain', 'one', 1, 'real',
     &    'Basin area-weighted average of rainfall minus interception',
     &    'inches',
     &   basin_net_rain).ne.0) return
     
      if(declvar('snow', 'psoilbasin', 'one', 1, 'real',
     &    'Basin area-weighted average of precipitation, irrigation '//
     &    'and thrufall on non-snow-covered areas',
     &    'inches',
     &   PSOILBASIN).ne.0) return

      if(declvar('snow', 'coverbasin', 'one', 1, 'real',
     &    'Basin area-weighted average of percentage of snowcover',
     &    'none',
     &   COVERBASIN).ne.0) return

      if(declvar('snow', 'basin_pweqv',  'one', 1, 'real',
     &    'Basin area-weighted average of total water in the snowpack',
     &    'inches',
     &   BASIN_PWEQV).ne.0) return
!
      if(declvar('snow', 'basin_snowevap',  'one', 1, 'real',
     &    'Basin area-weighted average of sublimation of the snowpack',
     &    'inches',
     &   BASIN_SNOWEVAP).ne.0) return
!
      ALLOCATE (PKWATER_EQUIV(Nmru))
      if(declvar('snow', 'pkwater_equiv','nmru', nmru, 'real',
     &    'snow-water-equivalent of snowpack for this day for each MRU',
     &    'inches',
     &   PKWATER_EQUIV).ne.0) return

      ALLOCATE (SNOW_EVAP(Nmru))
      if(declvar('snow', 'snow_evap',    'nmru', nmru, 'real',
     &    'sublimation from the snowpack for this day for each MRU',
     &    'inches',
     &   SNOW_EVAP).ne.0) return

      ALLOCATE (CHNGPCKMRU(Nmru))
      if(declvar('snow', 'chngpckmru',   'nmru', nmru, 'real',
     &    'change in water in snowpack for this day for each MRU',
     &    'inches',
     &   CHNGPCKMRU).ne.0) return

      ALLOCATE (PPTMIX_NOPACK(Nmru))
      if(declvar('snow', 'pptmix_nopack','nmru', nmru, 'integer',
     &     'variable set =1 if rain and snow same day with no pack ',
     &     'none',
     &   PPTMIX_NOPACK).ne.0) return

      ALLOCATE (TINDXI(Nmru))
      if(declvar('snow', 'tindxi',   'nmru', nmru, 'real',
     &    'antecedent snowpack temperature index for each MRU',
     &    'degrees F',
     &   TINDXI).ne.0) return

! OLD PARAMETERS USED IN THIS MODULE
! Copied parameter declarations from basin_prms.f

      ALLOCATE (ELEV(Nmru))  ! renamed in getparam
      if(declparam('basin', 'mru_elev', 'nmru', 'real',
     +   '0.', '-300.', '10000',
     +   'Mean elevation for each MRU',
     +   'Mean elevation for each MRU',
     +   'meters').ne.0) return
     
      ALLOCATE (mruarea(Nmru)) ! renamed in getparam
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('precip', 'tmax_allsnow_c', 'one', 'real',
     +   '0', '-10.', '10.',
     +   'All snow if tmax_c< this value; all rain if tmin_c>this '//
     $   'value.',' If MRU maximum temperature is below this value, '//
     +   'precipitation is assumed to be snow; alternately, if MRU '//
     +   'minimum temperature is above this value, precipitation is '//
     +   'assumed to be all rain.','degrees celsius').ne.0) return

      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return
      
! NEW PARAMETERS USED IN THIS MODULE----WILL BE INITIALIZED

      ALLOCATE (SUBRATE(Nmru))
      if(declparam('snow', 'SUBRATE',    'nmru', 'real',
     &   '0.0', '0.0', '.2',
     &   'Average daily snowpack sublimation rate for each MRU',
     &   'Average daily snowpack sublimation rate for each MRU, '//
     &   'INITIALIZE=0.0 if you '//
     &   'do not want to calculate snow evaporation, rate is * 0.5 '//
     &   'for each of the 2 daylight 6-hour periods',
     &   'inches/day').ne.0) return

      if(declparam('snow', 'ALAT',    'one', 'real',
     &   '47.0', '0.0', '90.',
     &   'Average latitude of modeled region',
     &   'Average latitude of modeled region, only need to input > 54',
     &   'degrees').ne.0) return

      ALLOCATE (TMXPREI(Nmru))
!      if(declparam('snow', 'TMXPREI', 'nmru', 'real',
!     &    '7.0', '0.', '50.',
!     &    'Init. daily maximum temperature from previous day',
!     &    'Init. daily maximum temperature from previous day',
!     &    'degrees C').ne.0) return
!
      ALLOCATE (ACUMX(Nmru))
      if(declparam('snow', 'ACUMX',   'nmru', 'real',
     &    '0.', '0.', '1000.',
     &    'Init. state: max. water equiv. since snow began to accum.',
     &    'Init. state: max. water equiv. since snow began to accum.',
     &    'inches').ne.0) return
  
      ALLOCATE (WEI(Nmru))
      if(declparam('snow', 'WEI',     'nmru', 'real',
     &    '0.', '0.', '1000.',
     &    'Initial water equivalent for each MRU',
     &    'Initial water equivalent for each MRU',
     &    'inches').ne.0) return
    
      ALLOCATE (NEGHSI(Nmru))
      if(declparam('snow', 'NEGHSI',  'nmru', 'real',
     &    '0.', '0.', '100.',
     &    'Initial negative heat storage for each MRU',
     &    'Initial negative heat storage for each MRU',
     &    'inches').ne.0) return    
 
      ALLOCATE (LIQWI(Nmru))
!       if(declparam('snow', 'LIQWI',  'nmru', 'real',
!     &   '0.', '0.', '100.',
!     &   'Initial liquid water in the snowpack for each MRU',
!     &   'Initial liquid water in the snowpack for each MRU',
!     &   'none').ne.0) return

      ALLOCATE (SBWSI(Nmru))
!       if(declparam('snow', 'SBWSI',  'nmru', 'real',
!     &    '0.', '0.', '1000.',
!     &    'Init. areal water equ. after new snowfall - see depl. curve',
!     &    'Init. areal water equ. after new snowfall - see depl. curve',
!     &    'inches').ne.0) return
    
      ALLOCATE (SBI(Nmru))
!      if(declparam('snow', 'SBI',    'nmru', 'real',
!     &    '0.', '0.', '1000.',
!     &    'Initial areal water equivalent just prior to the snowfall'//
!     &    ' - see depl. curve',
!     &    'Initial areal water equivalent just prior to the snowfall'//
!     &    ' - see depl. curve',
!     &    'inches').ne.0) return    
 
      ALLOCATE (AESCI(Nmru))
      ALLOCATE (SBAESCI(Nmru))
!      if(declparam('snow', 'SBAESCI', 'nmru', 'real',
!     &    '0.', '0.', '1.',
!     &    'Init. areal extent of snow cover prior to the new '//
!     &    'snowfall - see depletion curve',
!     &    'Init. areal extent of snow cover prior to the new '//
!     &    'snowfall - see depletion curve',
!     &    'none').ne.0) return 

      ALLOCATE (STOREI(Nmru))
!      if(declparam('snow', 'STOREI',  'nmru', 'real',
!     &    '0.', '0.', '100.',
!     &    'Initial routed liquid water in storage in the pack'//
!     &    ' - see lag function',
!     &    'Initial routed liquid water in storage in the pack'//
!     &    ' - see lag function',
!     &    'inches').ne.0) return 
 
 
      ALLOCATE (LAGRO(Nexlag)) ! local copy in run
      ALLOCATE (LAGROI(Nmru,Nexlag))
!      if(declparam('snow', 'LAGROI','nmru,nexlag', 'real',
!     &    '0.', '0.', '100.',
!     &    'Initial amount of lagged inflow for the first period'//
!     &    ' see lag- and attenuation function respectively',
!     &    'Initial amount of lagged inflow for the first period'//
!     &    ' see lag- and attenuation function respectively',
!     &    'inches/six hours').ne.0) return

      ALLOCATE (UADJ(Nmru))
      if(declparam('snow', 'UADJ',    'nmru', 'real',
     &   '0.05', '0.03', '0.2',
     &   'Average 6-hour wind function during rain on snow events',
     &   'Average 6-hour wind function during rain on snow events',
     &   'mm/(mb 6hr)').ne.0) return

c
c  Removed SCF. Snow undercatch corrections taken care of with
c  snow_adj(nmru,nmonths) in precip.
c
C      if(declparam('snow', 'SCF',    'nmru', 'real',
C     &   '1.15', '0.', '10.',
C     &   'snowfall correction factor',
C     &   'Multiplying factor to correct for precipitation gage catch'//
C     &   'deficiency during periods of snowfall',
C     &   'none').ne.0) return
c
      ALLOCATE (DAYGM(Nmru))
      if(declparam('snow', 'DAYGM',   'nmru', 'real',
     &   '0.1', '0.', '1.',
     &   'Ave. daily melt at the snow-soil interface during winter',
     &   'Ave. daily melt at the snow-soil interface during winter',
     &   'inches').ne.0) return
 
      ALLOCATE (MFMAX(Nmru))
      if(declparam('snow', 'MFMAX',   'nmru', 'real',
     &   '0.6', '0.5', '2.',
     &   'Maximum non-rain melt factor which occurs on June 21',
     &   'Maximum non-rain melt factor which occurs on June 21',
     &   'mm/(6hr.degrees C)').ne.0) return

      ALLOCATE (MFMIN(Nmru))
      if(declparam('snow', 'MFMIN',   'nmru', 'real',
     &   '0.2', '0.2', '1.',
     &   'Minimum non-rain melt factor which occurs on Dec 21',
     &   'Minimum non-rain melt factor which occurs on Dec 21',
     &   'mm/(6hr.degrees C)').ne.0) return

      ALLOCATE (NMF(Nmru))
      if(declparam('snow', 'NMF',     'nmru', 'real',
     &   '0.15', '0.05', '0.5',
     &   'Maximum value of negative melt factor which'//
     &   ' occurs on June 21',
     &   'Maximum value of negative melt factor which'//
     &   ' occurs on June 21',
     &   'mm/(6hr.degrees C)').ne.0) return

      ALLOCATE (MBASE(Nmru))
      if(declparam('snow', 'MBASE',   'nmru', 'real',
     &   '32.', '32.', '40.',
     &   'Base temperature for melt computations during non-rain per.',
     &   'Base temperature for melt computations during non-rain per.',
     &   'degrees F').ne.0) return
c
c     Made PXTEMP a local variable to be assigned the monthly value of
c     tmax_allsnow_c from the precip module
c
c      if(declparam('snow', 'PXTEMP',  'nmru', 'real',
c     &   '32.', '32.', '40.',
c     &   'temperature above which precipitation is assumed to be rain',
c     &   'temperature above which precipitation is assumed to be rain',
c     &   'degrees F').ne.0) return
c

      ALLOCATE (PLWHC(Nmru))
      if(declparam('snow', 'PLWHC',   'nmru', 'real',
     &   '0.05', '0.01', '0.4',
     &   'Percent (decimal) liquid water holding capacity',
     &   'Percent (decimal) liquid water holding capacity',
     &   'none').ne.0) return

      ALLOCATE (TIPM(Nmru))
      if(declparam('snow', 'TIPM',    'nmru', 'real',
     &   '0.3', '0.2', '0.6',
     &   'Antecedent temperature index parameter. 0.2 for'//
     &   ' Thick pack; 0.5 for Thin Pack',
     &   'Antecedent temperature index parameter. 0.2 for'//
     &   ' Thick pack; 0.5 for Thin Pack',
     &   'none').ne.0) return

      ALLOCATE (SI(Nmru))
      if(declparam('snow', 'SI',     'nmru', 'real',
     &   '5.', '0.0', '100.0',
     &   'Snow water equivalence above which snowcovered area '//
     &   'reaches a maximum',
     &   'Snow water equivalence above which snowcovered area '//
     &   'reaches a maximum',
     &   'inches').ne.0) return    
 
      ALLOCATE (MRUDEPL(Nmru))
      if(declparam('snow', 'MRUDEPL', 'nmru', 'integer',
     &   '1', 'bounded', 'ndepl',
     &   'Identifier for which areal depletion curve to use',
     &   'Identifier for which areal depletion curve to use',
     &   'none').ne.0) return     

      ALLOCATE (ADC(ndeplval))
      if(declparam('snow', 'ADC', 'ndeplval', 'real',
     &   '0.5', '0.0', '1.0',
     &   'Points on depletion curve',
     &   'Points on depletion curve',
     &   'none').ne.0) return     

      ALLOCATE (pmru(Nmru))
      ALLOCATE (net_rain(Nmru))
      ALLOCATE (net_snow(Nmru))
      ALLOCATE (tmxmru(Nmru))
      ALLOCATE (tmnmru(Nmru))
      ALLOCATE (pxtemp(Nmru))

      nwsmdecl = 0

      return
      end

!***********************************************************************
!
!     nwsminit - Initialize snow module - get parameter values
!     --------    and initial snowpack state data 
!
      integer function nwsminit()

      USE WEBMOD_SNOW

      integer n, n1, icrve, npts
      real r, fn, ai
      
      nwsminit = 1

! GET BASIN AREA 
      if(getparam('basin', 'basin_area',    1, 'real', BASAREA)
     &   .ne.0) return
! GET MRU AREA 
      if(getparam('basin', 'mru_area', nmru, 'real', MRUAREA)
     &   .ne.0) return
! GET ELEVATION OF MRUS (in meters)
      if(getparam('basin', 'mru_elev', nmru, 'real', ELEV)
     &   .ne.0) return
! Get the tmax_allsnow_c parameter. In precip modules, all rain if tmin_c > tmax_allsnow_c;
! all snow if tmax_c < tmax_allsnow_c. Value is usually aroung 0 degrees Celsius and is
! used here to determine if the precip during a given 6-hr period is rain or snow.
! tmax_allsnow_c is converted to Fahrenheit, the units of PXTEMP.

      if(getparam('precip', 'tmax_allsnow_c', 1, 'real', 
     +   tmax_allsnow_c).ne.0) return

      if(getparam('snow', 'MRUDEPL',   nmru, 'integer',MRUDEPL)
     &   .ne.0) return

      if(getparam('snow', 'SUBRATE',    nmru, 'real', SUBRATE)
     &   .ne.0) return

      if(getparam('snow', 'ALAT',           1, 'real', ALAT)
     &   .ne.0) return

!      if(getparam('snow', 'TMXPREI',   nmru, 'real', TMXPREI)
!     &   .ne.0) return

      if(getparam('snow', 'ACUMX',     nmru, 'real', ACUMX)
     &   .ne.0) return

      if(getparam('snow', 'WEI',       nmru, 'real', WEI)
     &   .ne.0) return

      if(getparam('snow', 'NEGHSI',    nmru, 'real', NEGHSI)
     &   .ne.0) return
     
!      if(getparam('snow', 'LIQWI',     nmru, 'real', LIQWI)
!     &   .ne.0) return
!
!      if(getparam('snow', 'SBWSI',     nmru, 'real', SBWSI)
!     &   .ne.0) return
 
!      if(getparam('snow', 'SBI',       nmru, 'real', SBI)
!     &   .ne.0) return
!
!      if(getparam('snow', 'SBAESCI',   nmru, 'real', SBAESCI)
!     &   .ne.0) return

!      if(getparam('snow', 'STOREI',    nmru, 'real', STOREI)
!     &   .ne.0) return

      if(getparam('snow', 'UADJ',      nmru, 'real', UADJ)
     &   .ne.0) return

c      if(getparam('snow', 'SCF',       nmru, 'real', SCF)
c     &   .ne.0) return

      if(getparam('snow', 'DAYGM',     nmru, 'real', DAYGM)
     &   .ne.0) return

      if(getparam('snow', 'MFMAX',     nmru, 'real', MFMAX)
     &   .ne.0) return

      if(getparam('snow', 'MFMIN',     nmru, 'real', MFMIN)
     &   .ne.0) return

      if(getparam('snow', 'NMF',       nmru, 'real', NMF)
     &   .ne.0) return

      if(getparam('snow', 'MBASE',     nmru, 'real', MBASE)
     &   .ne.0) return

      if(getparam('snow', 'PLWHC',     nmru, 'real', PLWHC)
     &   .ne.0) return

      if(getparam('snow', 'TIPM',      nmru, 'real', TIPM)
     &   .ne.0) return

      if(getparam('snow', 'SI',        nmru, 'real', SI)
     &   .ne.0) return

!      if(getparam('snow', 'LAGROI',  nmru*nexlag, 'real', LAGROI)
!     &   .ne.0) return

      if(getparam('snow', 'ADC',     ndeplval, 'real', ADC)
     &   .ne.0) return

! ZERO INITIAL VARIABLES

      DO 10 N=1,nmru
      TINDXI(N)       =0.0
      SNOWMELT(N)     =0.0
      PSOILMRU(N)     =0.0
      SNOWCOV_AREA(N) =0.0
      SNOW_EVAP(N)    =0.0
      PPTMIX_NOPACK(N)=0
!
! No liquid in initial snowpack. Lagged intervals set at two
!
      LAGROI(N,1) = 0.0
      LAGROI(N,2) = 0.0
      LIQWI(N) = 0.0
      STOREI(N) = 0.0
!
! Establish initial point on curve using WEI, ACCMX, SI, and ADC. Assume no liquid water or storage - RMTW
!
      AI = ACUMX(N)
      if(AI.GT.SI(N)) AI=SI(N)
      SBWSI(N) = WEI(N)
      ICRVE = MRUDEPL(N)-1
      if(WEI(N).GE.ACUMX(N)) THEN  ! set as accumulation period with AESC
        ACUMX(N)=WEI(N)
        SBI(N)=WEI(N)
        AESCI(N)=1.0
        SBAESCI(N)=1.0
      ELSE                         ! set as melt period, retaining ACUMX and an AESC max of 1.0 if WEI/AI > 1.0
        SBI(N)= WEI(N) + 0.05      ! 1.27 mm converted to inches
        R = (WEI(N)/AI)*10.0+1.0   ! Find point on ADC curve
        if(R.GE.11.0) then
          R = 1.0
          N1 = 10
        else
          N1  = R
          FN = N1
          R  = R - FN
        endif
        NPTS=ndeplval/ndepl        
        N1 = N1+NPTS*ICRVE ! assign lineear array index
        AESCI(N) = ADC(N1) + (ADC(N1+1)-ADC(N1))*R  ! INITIAL AREAL EXTENT on Curve, Must be less than or equal to 1.0
        SBAESCI(N) = AESCI(N)

      ENDIF 
 10   CONTINUE
!
! debug for sca
!
!      write(25,'(A)')" In: Yr Mo Dy 6hr WE LIQW ACCMAX "//
!     &  "SIT SB SBWS AESC SBAESC Curve"
!
      nwsminit = 0

      return
      end

!***********************************************************************
!
!      nwsmrun - Compute the snowmelt, the snowcover and the rain
!      -------   falling on bare ground for the each mru (Hydrologic
!                Response Unit) and the whole AREA. For EACH mru,
!                do TIME step, generally assumed to be 6 HOURS.
!               
      integer function nwsmrun()

      USE WEBMOD_SNOW
      USE WEBMOD_OBSHYD, ONLY : datetime

      integer I1, I2, I, jday, whichcurve, im
      integer itpx,ndt, nstep
      logical nosnow
c
c Added the following line to assign tmax_allrain_c to PXTEMP
c
c commented out since tmax_allsnow_c is now used
c
c      integer nowtime(6), mo
c
c Conversion of mm to cubic meters
c
      real cm(nmru),ci(nmru)
! PRECIP INFORMATION

!  SNOW MODEL PARAMETERS

! PUBLIC VARIABLES: for this module

! INITIAL/CARRYOVER PARAMETERS (MAJOR SNOW VARIABLES)--PRIVATE FOR THIS MODULE
   
! AIR temperatures: four 6-hour periods,temp of rain,snow,air,precp, past Temps
      real T6HR(4), TR, TS, TAIR, TPX, TMNPST, TMXPRE
! Precipitation
      real PXI
! PARAMETERS FOR EQUATIONS FOR THIS DAY FOR MRU (many have 'T' appended to names given
!                                                above, or else 'I' at end is deleted)
c      real SCFACT
      real RNorSNO
      real MFMAXT, MFMINT,  NMFT,    BASEM,   MF,     RATIO
      real NMINDEX,NMRATE,  TINDEX,  TSPM
      real RMIN,   SNEW,    RFMIN,   SFNEW
      real UADJT,  ELEV100m, PA,     SBC,     SBCI
      real WE,     NEGHS,   LIQW,    STORAGE
      real SFALL,  PXSOIL,  PACKRO,  EVAP,    DELT
      real DSFALL, DRAIN,   DQNET,   BEGINPACK
      real MELT,   RAINM,   QNET,    CNHS,    CNHSPX
      real GM,     GMRO,    GMWLOS,  RAIN,    TWE
      real WCPLW,  EXCESS,  LAGGED
      real SIT,    ACCMAX,  AESC, SBAESC,  SBWS,   SB
c$$$      real tavgc(nmru)

      nwsmrun = 1
  
      nstep=getstep()  ! add for debugging

c     begin initial values, set some values 
c----------------------------------------------------------------------------------
      
      SBC  = 0.0612  ! STEFAN/BOLTZMAN CONSTANT--MM/(((DEGK/100)**4)*HR)
      RMIN = 0.25    ! IF RAIN EXCEEDS RMIN/HR, THEN USE RAIN-ON-SNOW EQUATION
      SNEW = 1.5     ! IF SNOWFALL > SNEW/HR--THEN TINDEX=TPX
!
      BASIN_SNOWMELT = 0.0     ! ZERO THE BASIN AREA-WEIGHTED VALUES
      PSOILBASIN     = 0.0
      COVERBASIN     = 0.0
      BASIN_PWEQV    = 0.0
      BASIN_SNOWEVAP = 0.0
      basin_net_snow = 0.0
      basin_net_rain = 0.0
!
!--------------------------------------------------------------------------------------
! GET MRU TEMPERATURES AND PRECIPITATION, julian day

      if(getvar('temp',   'tmax_c',    nmru, 'real', tmxmru) 
     &   .ne.0) return
      if(getvar('temp',   'tmin_c',    nmru, 'real', tmnmru)
     &   .ne.0) return
      if(getvar('intcp', 'net_dep',   nmru, 'real', pmru  )
     &   .ne.0) return

! Get net snow,rain 
! Nets have been calculated but will be reset in this program since calculated here
      if(getvar('intcp', 'net_snow',  nmru, 'real', net_snow)
     &   .ne.0) return
      if(getvar('intcp', 'net_rain',  nmru, 'real', net_rain)
     &   .ne.0) return

c$$$      if(get*var('temp', 'tavgc', nmru, 'real', tavgc)
c$$$     +   .ne.0) return
c$$$
c     The following two lines need to be used if tmax_allrain_c(nmonths) is to be used
c     instead of tmax_allsnow_c.
c
c      call dattim('now', nowtime)
c      mo = nowtime(2)

c The julian day is needed to calculate the meltfactor
c As documented by Rob Payne, Day one must be the first day of spring
c for the Melt Factor to be computed with the correct phase - RMTW
c

      jday = julian('now', 'spring')

!
!--------------------------------------------------------------------------------------
!                                       -*-*--*-*--*-*--*-*--*-*--*-*--*-*--*-*--*-*-
! LOOP THROUGH EACH MRU (1000 LOOP), FOR EACH MRU GO THROUGH TIME LOOP (900 LOOP)
!                        ---------                                      --------
      do 1000 I1=1,nmru

         cm(i1)=mruarea(i1)*1000
         ci(i1)=mruarea(i1)*1000*25.4

      SUBRATE(I1)  = SUBRATE(I1) * 25.4     ! CONVERT DAILY SUBLIMATION RATE TO mm,reset @ end
      SNOWMELT(I1)  = 0.0
      PSOILMRU(I1)  = 0.0
      SNOW_EVAP(I1) = 0.0
      CHNGPCKMRU(I1)= 0.0
      PPTMIX_NOPACK(I1) = 0

      pxtemp(I1) = tmax_allsnow_c

c
c Modify slightly so that pkwater_equivalence is set to zero so that 
c webmod knows to skip snow chemistry section
c

c$$$      if ((WEI(I1) .EQ. 0.0) .AND. (PMRU(I1) .EQ. 0.0)) GO TO 1000  !--No need to do ANY updates 
c$$$

      if ((WEI(I1) .EQ. 0.0) .AND. (PMRU(I1) .EQ. 0.0)) then !--No need to do ANY updates 
         pkwater_equiv(i1) = 0.0
         go to 1000
      end if
!                                                                   ! NO PRECIP OR PACK
!-CALCULATE PERCENT OF TOTAL AREA FOR AREA-WEIGHTING SUMS
      perAREA = mruAREA(I1) / basAREA
!
!--------------------------------------------------------------------------------------
!--Initial SNOW-MODEL parameter values

      UADJT  = UADJ(I1)        ! wind function for mru
c      SCFACT = SCF(I1)         ! correction factor for snow at gages
      GM     = DAYGM(I1)*25.4  ! daily average ground-snow interface melt
      MFMAXT = MFMAX(I1)       ! maximum melt coefficient
      MFMINT = MFMIN(I1)       ! minimum melt coefficient
      NMFT   = NMF(I1)         ! negative melt factor
      BASEM  = MBASE(I1)       ! temperature above which snow melts
      BASEM  = (BASEM - 32.)*0.5555556     !CONVERT TO DEGREES C
      RNorSNO= PXTEMP(I1)      ! temperature that determines if rain or snow
c      RNorSNO= (RNorSNO - 32.)*0.5555556   !CONVERT TO DEGREES C (RNorSNO now - tmax_allsnow_c
c                                            so no conversion to celsius needed
      WCPLW  = PLWHC(I1)       ! percent liquid water holding capacity of snow
      TSPM   = TIPM(I1)        ! antecendant temperature index, > .5 past few 6hrs
!                                                         .2 past 1-7 days
! Get which areal depletion curve to use for this MRU

      whichcurve = MRUDEPL(I1)
! Eliminate initial assignment of yesterdays tmax by assuming the same as today's tmax
      if(nstep.eq.1) TMXPREI(I1) = TMXMRU(I1)

!---------------------------------------------------------------------------------------
!--Initial snowpack conditions, from yesterday (or for the 1st initializing)
!                                   ! Metric units, so convert inches to mm, Farh to Celsis
      WE     = WEI(I1)*25.4    ! WATER-EQUIVALENT IN FROZEN STATE
      TMXPRE = TMXPREI(I1)     ! Yesterdays's daily maximum temperature
      SIT    = si(I1)    *25.4 ! water-equiv above which 100% cover always occurs
      ACCMAX = ACUMX(I1) *25.4 ! max. water equivalent since snow began accumulating
      NEGHS  = NEGHSI(I1)*25.4 ! current negative heat storage for mru
      LIQW   = LIQWI(I1) *25.4 ! liquid water in pack for mru
      SBWS   = SBWSI(I1) *25.4 ! AREAl water equivalent after NEW snowfall
      TINDEX = TINDXI(I1)      ! antecedent temperature index for mru
      TINDEX = (TINDEX-32.)*0.5555556              !CONVERT TO DEGREES C
      SB     = SBI(I1) *25.4   ! AREAl water equivalent PRIOR to snowfall
      SBAESC = SBAESCI(I1)     ! AREAl extent snowcover prior to new snowfall for mru
      AESC   = AESCI(I1)          ! Areal extent of snowcover at the end of the day  ! Added by RMTW
      STORAGE  = STOREI(I1)  *25.4  ! the ROUTED liquid water in storage in pack for mru
      DO 75 I=1,nexlag
   75 LAGRO(I) = LAGROI(I1,I)*25.4  ! lagged/attenuated water in snowpack for mru
!                                   ! GENERALLY SET MAXLAG--nexlag =2
!-----------------------------------------------------------------------------------------
! SET the 4 6-hour temperature values, calculate the average, SAVE TODAYS MAX
      
      TMNPST =TMNMRU(I1) - TMXMRU(I1) + TMXPRE

      T6HR(1)=0.95  * TMNMRU(I1 )+ 0.05  * TMXPRE
      T6HR(2)=0.40  * TMNMRU(I1) + 0.60  * TMXMRU(I1) 
      T6HR(3)=0.925 * TMXMRU(I1) + 0.025 * TMNMRU(I1) + 0.05 * TMNPST
      T6HR(4)=0.33  * TMXMRU(I1) + 0.67  * TMNPST
 
      TMXPREI(I1 )= TMXMRU(I1)        ! SAVE TODAY'S MAX FOR TOMOOROW'S CALCULATIONS



!--------------------------------------------------------------------------------------
!   BEGINPACK=initial total water in snowpack (WATER EQ + LIQUID + LAGGED + STOR)

      LAGGED = 0
      DO 90 I=1,nexlag
   90 LAGGED = LAGGED + LAGRO(I)

      BEGINPACK = WE + LIQW + LAGGED + STORAGE ! STORE BEGINNING OF DAY TOTAL WATER IN MRU
!                                                      ----------------------
!----------------------------------------------------------------------------------------
!  initial daily values, including melt coefficient for this day, AND ATMOSPHERIC PRESSURE

      DSFALL = 0.0      ! DAILY SNOWFALL
      DRAIN  = 0.0      ! DAILY RAINFALL
      DQNET  = 0.0      ! DAILY NET RADIATION, CALCULATED

!--CALCULATE the meltfactor for this day
      call MELTFACT(JDAY,ALAT,MFMAXT,MFMINT,MF,RATIO)
!--CALCULATE atmospheric pressure, 1st elevation in 100's of meters
c$$$      ELEV100m = ELEV(I1)*0.01*0.3048  
c modified to match mru elevation in meters
      ELEV100m = ELEV(I1)*0.01  
      PA = 1012.4 - 11.34 * ELEV100M + 0.00745 * (ELEV100M ** 2.4)  ! ATMOSPHERIC PRESSURE

! End initial values
! --GET TIME STEP, IN HOURS, USING MMS LIBRARY FUNCTION
      DELT=deltim()
      NDT =24./DELT

      ITPX=24./NDT        ! integer of time-step,eg, 6.0 hr = DELT => 6 = ITPX

!---------------------------------------------------------*-*--*-*--*-*--*-*--*-*--
! --Beginning of 6 hour loop (900 loop), ASSUMING DELT=6 HOURS
!---------------------------------------------------------
! SET TO 6-HOURS FOR RIGHT NOW
      if (NDT .NE.4)  NDT=4
      if (ITPX.NE.6) ITPX=6

! --PUT HOURLY VALUES TO TIME-STEP VALUES (ASSUMING 6HRS, OK FOR HOURLY HERE)
      SFNEW =SNEW * ITPX
      RFMIN =RMIN * ITPX
      SBCI  =SBC  * ITPX

! --PUT AVERAGE DAILY GROUND MELT INTO MELT PER DELT 
      GM = GM / NDT

      DO 900 I2=1,NDT            !---DO 900-LOOP GOES BY NDT, EG: DELT=6, NDT=4 PERIODS

!---------------------------------------------------------------------------------
!  Initial 6-hour values, STILL OK FOR HOURLY VALUES, WILL NOTE WHEN NOT

      EVAP   = 0.0
      MELT   = 0.0
      RAIN   = 0.0
      RAINM  = 0.0
      CNHS   = 0.0
      CNHSPX = 0.0
      SFALL  = 0.0
      GMRO   = 0.0
      PXSOIL = 0.0
      NOSNOW = .FALSE.

!------------------------------------------------------------------------------------
!--Convert the WEATHER DATA (PRECIP TO mm, TEMP to cent), and put to time interval

      PXI   =  PMRU(I1)*25.4 / NDT
!---------------------------------------------------------- ---------------------------
! 6hr temp value                             

      TAIR  = T6HR(I2)      !  if NDT > 4, then this will bomb, could modify for > 4
      TPX   = TAIR          !  TPX = temperature of precip
!---------------------------------------------------------------------------------------
! IF NO PRECIP (but there is a pack), SKIP TO GROUND-MELT COMPUTATIONS

      if (PXI. EQ. 0.0)     GO TO 120
     
      if (TPX .GT. RNorSNO) go to 100  ! the fraction (PCTS()) of ANDERSON is not
!                                        used here only if TPX > or <= RnorSNO (~ 32F)
!-----------------------------------------------------------------------------------------

!-SNOWFALL:  temp <= RNorSNO

      TS = TPX
      if (TS .GT. 0.0) TS=0.0  

c      PXI    = PXI * SCFACT              ! Adjust for snowfall catch deficiency
      SFALL  = PXI                       ! Anderson further modified based on IFUT
      DSFALL = DSFALL + SFALL            ! this adjustment also used for others, see below
!                                          ALSO, no checks w/ Anderson's 'snof', Fuch uses 1.27 mm
      SBWS   = WE + LIQW + 0.75*PXI      ! Establish 75% of new snow as the threshold to drop to
                                         ! before returning to depletion curve
      WE     = WE + SFALL                    ! add snowfall to water equivalent 

! IF we+liqw => 3* sb, assume new accumulation period, (not in from Fuchs)
      if( (WE + LIQW).GE. (3.0*SB) ) ACCMAX = WE + LIQW

      CNHSPX = - TS*SFALL/160.       !  PRECP (-) HEAT store VALUE, WHEN < 0C
!                                          
      if (PXI .GT. SFNEW) TINDEX=TS              ! new snow > rate, set temp index snow to air temp
      PXI = 0.0                      ! SET PRECIP=0, SNOW ACCOUNTED FOR, WILL SKIP RAIN

      go to 120       ! Anderson did following calculations even if all snow,
!                       the 'go to 120' skips this for computational speed
!-----------------------------------------------------------------------------------
!-RAIN:   (temp > RnorSNO)
   
  100 if (WE .EQ. 0.0) THEN                        ! NO SNOW, SO ALL PRECIP GOES TO SNOW-FREE 
          PSOILMRU(I1) = PSOILMRU(I1) + PXI/25.4   ! SOIL COMPONENT, SO CAN'T MELT, SKIP REST
          DRAIN        = DRAIN        + PXI
          GO TO 300
      END IF

      TR=TPX
      IF(TR.LT.0.0) TR=0.0
      RAINM = 0.0125 * PXI * TR
      DRAIN = DRAIN  + PXI
      RAIN  = PXI                      ! RAIN IS THE RAINFALL FOR INTERVAL, BUT ON SNOW&BAREGROUND
!                                      ! RAIN WILL LATER BE REDUCED FOR RAIN ON SOILS=PXSOIL
!--------------------------------------------------------------------------------------
!  SNOWPACK computations: GROUNDMELT, ENERGY EXCHANGE, MELT, SURFACE-EXCHANGE
!---------------------------------------------------------------------------------------
!
!--Groundmelt----

  120 if (WE. LE. GM) THEN         ! WATER EQIV WILL BE REDUCED BY GM, OTHERWISE ALL
!                                    OF THE SNOW (old + NEW) WILL GO SINCE WE <= THE GM
       GMRO   = WE+LIQW
       RAIN   = 0.0             ! PXI=0 IF SNOW SO ONLY MATTERS IF RAIN,
       PXSOIL = PXI             ! IE, SNOW IS GONE SO RAIN GOES TO BARE GROUND 
        GO to 300               ! ALL SNOW WENT TO GROUNDMELT RUNOFF, GMRO, GO TO 300
!                              
      ELSE

       GMWLOS = (GM/WE) * LIQW      ! WATER LOSS FROM PACK LIQUID WATER DUE TO GROUNDMELT
       GMRO   =  GM + GMWLOS           ! GROUNDMELT RUNOFF = GROUNDMELT + LIQUID WATER LOSS

      END IF

!---------------------------------------------------------------------------------
!--Snowpack energy exchange computations, INCLUDES rest-of of Anderson's 'melt19' subroutine
!   compute melt and negative heat exchange index temperatures, then NEG HEAT EXCHANGE

      CALL MELTRATE(TAIR,BASEM,NMINDEX,TINDEX, NMRATE, RATIO,NMFT, CNHS,
     &TSPM,RAIN, RFMIN,  MELT, MF,RAINM, SBCI, UADJT,PA,EVAP,
     &SUBRATE,I2, I1)

!---------------------------------------------------------------------------------
!-AREAL extent of SNOW cover, Anderson's suboutine 'AESC19' SLIGHTLY MODIFIED
!                                                           ACCOUNTS FOR NEW SNOW (IF ANY)
!
!!debug/ for sca
!!
!      if(i1.eq.1) then   
!      write(25,121)"In:",(datetime(im), im=1,3), i2, we, LIQW,ACCMAX,
!     &  SIT,SB,SBWS,AESC,SBAESC,whichcurve
!      end if
! 121  format(A5,4I5, 8F12.3, i2)
!! /debug
!!      if (indx.gt.200) write(25,122)nstep,indx,indxb, imetric,
!!     $     (c_chem(indx)%M(5,im),c_chem(indx)%vol(im),im = 1,5)
!! 122  format(I4, 2I10, i4, 10F12.1)
!! Move to to maintain SCA as last step
!!      CALL AESNOW(WE,LIQW,ACCMAX,SIT,SB,SBWS,AESC,SBAESC,ADC,whichcurve)
!! debug/
!!
!      if(i1.eq.1) then      
!      write(25,121)"Out:",(datetime(im), im=1,3), i2, we, LIQW,ACCMAX,
!     &  SIT,SB,SBWS,AESC,SBAESC,whichcurve
!      end if
!!  /debug
!--------------------------------------------------------------------------------
!--Adjust period values for AESC < 1.0

      if (AESC .LT. 1.0)        THEN
         GM     = GM     * AESC
         GMRO   = GMRO   * AESC
         GMWLOS = GMWLOS * AESC
         MELT   = MELT   * AESC
         CNHS   = CNHS   * AESC
         EVAP   = EVAP   * AESC

         PXSOIL = RAIN   * (1.0-AESC)  ! RAIN IS RAIN ON: BARE GROUND AND ON SNOW
         RAIN   = RAIN   -  PXSOIL     ! RAIN NOW BECOMES RAIN ON SNOW
      ELSE                             ! --------
          PXSOIL = 0.0
      END IF
!-----------------------------------------------------------------------------------
!-Adjust for surface and groundmelt
      
       WE   = WE   - GM                   ! DECREASE WATER EQUIV  BY GROUNDMELT
       LIQW = LIQW - GMWLOS               ! DECREASE  LIQUID WATER BY WATER-LOSS
!                                                    ! check chns, IF IT BRINGS NEGHS TO NEGATIVE
      if ((CNHS+NEGHS).LT. 0.0) CNHS = -1.0*NEGHS    ! THEN DECREASE BY ABSOLUTE DIFF,IE,= NEG NEGHS
!---------------------------------------------------------------------------------
!-Surface melt computations, AT THIS POINT MELT IS ALWAYS => 0.0 SO CLEANED CODE

      if (MELT. GE. WE)   THEN
        MELT  = WE    + LIQW            ! MELT => WE, SO MELT IS EQUAL TO TWE
        DQNET = DQNET + MELT            ! SET DQNET FOR THIS INTERVAL (DELT), AND FINISH DELT
        EVAP = 0                        ! no snow left for sublimation - RMTW
        Go to 300
      ELSE
        WE = WE - MELT                  ! MELT  < WE, REDUCE WE BY MELT, THEN SET DQNET
      END IF
!------------------------------------------------------------------------------
!-ADJUST FOR SUBLIMATION, if parameter SUBLATE=0.0, then results same as Anderson
!
! Divert pack-eliminating sublimation to melt to avoid convergence problems 
! in PHREEQC (i.e. remaining solutes cannot evaporate, the must stay dissolved) - RMTW
!

      if (EVAP. GE. WE)    THEN
        EVAP  = 0.0
        MELT  = MELT+ WE+LIQW
        DQNET = DQNET + WE+LIQW
        GO to 300               ! IE, SNOW IS GONE SO RAIN GOES TO BARE GROUND
!                               ! REMAINING SNOW WENT TO SNOW EVAPORATION (melt for geochem), GO TO 300
      ELSE
        WE = WE - EVAP
      END IF
!-------------------------------------------------------------------------------
!--Qnet: net surface energy exchange FOR DELT WHEN ALL SNOW DID NOT MELT
      QNET  = MELT - CNHS - CNHSPX  + EVAP ! MELT (ENERGY,RAIN) + NEG HEAT SURFCE + ENERGY RAIN + EVP
      DQNET = DQNET + QNET
!-------------------------------------------------------------------------------
!--HEAT and Water balance for snowpack

      CALL HEATWATR(MELT,RAIN,CNHS,CNHSPX,WCPLW,WE,NEGHS,
     &              EXCESS,LIQW,TINDEX)

!-------------------------------------------------------------------------------
!-IF reached this far, then ROUTE water, add GMRO to PACKRO
!---Route excess and/or lagged water through the pack, this is Anderson's 'ROUT19' subroutine

      CALL ROUTE(ITPX,EXCESS,WE,AESC,STORAGE,LAGRO,PACKRO,NEXLAG)

      PACKRO = PACKRO + GMRO                             ! ADD GROUND MELT TO PACK RUNOFF

      GO TO 400
!--------------------------------------------------------------------------------
!-SNOW gone - layer or total pack: Initialize all variables for no snow
!-----------------------------------------------------
  300 CONTINUE            ! ALL went to runoff from pack
      LAGGED  = 0.0
      DO 350 I=1,nexlag
      LAGGED  = LAGRO(I)    + LAGGED
  350 LAGRO(I)= 0.0
      PACKRO  = GMRO + MELT + LAGGED + STORAGE + RAIN

      WE      = 0.0
      NEGHS   = 0.0
      LIQW    = 0.0
      SB      = 0.0
      SBAESC  = 0.0
      SBWS    = 0.0
      ACCMAX  = 0.0
      AESC    = 0.0
      TINDEX  = 0.0
      STORAGE = 0.0
      MELT    = 0.0
      NOSNOW  = .TRUE.
!
!-END of snowpack computations, SET VALUES FOR MRU (as INCHES) that are calc each Time-Step
!-----------------------------------------------------
  400 CONTINUE
      SNOWMELT(I1)   = SNOWMELT(I1)  + PACKRO / 25.4
      PSOILMRU(I1)   = PSOILMRU(I1)  + PXSOIL / 25.4
      SNOW_EVAP(I1)  = SNOW_EVAP(I1) + EVAP   / 25.4
      LAGGED  = 0.0
      DO 351 I=1,nexlag
      LAGGED  = LAGRO(I)    + LAGGED
  351 continue
!debug/ for sca
!
!      if(i1.eq.1) then   
!      write(25,121)"In:",(datetime(im), im=1,3), i2, we, LIQW,ACCMAX,
!     &  SIT,SB,SBWS,AESC,SBAESC,whichcurve
!      end if
! 121  format(A5,4I5, 8F12.3, i2)
! /debug
!      if (indx.gt.200) write(25,122)nstep,indx,indxb, imetric,
!     $     (c_chem(indx)%M(5,im),c_chem(indx)%vol(im),im = 1,5)
! 122  format(I4, 2I10, i4, 10F12.1)
  
! Added logic to update snow covered area only if snowpack remains.
! Moved AESNOW to last step of time loop for continuity of SB and SBAESC values
      if(.NOT.NOSNOW) CALL AESNOW(WE,LIQW,ACCMAX,SIT,SB,
     & SBWS,AESC,SBAESC,ADC,whichcurve)
!
! debug/
!
!      if(i1.eq.1) then      
!      write(25,121)"Out:",(datetime(im), im=1,3), i2, we, LIQW,ACCMAX,
!     &  SIT,SB,SBWS,AESC,SBAESC,whichcurve
!      end if
!  /debug

  900 continue                                                    ! END OF TIME LOOP
!----------------------------------------------------*-*--*-*--*-*--*-*--*-*--*-*-
!-Snow carryover values: Calculate total water in pack then
!                        recalculate snow-cover using end of period values
      if (.NOT.NOSNOW) then
       LAGGED     = 0.0
       do 950 I=1,nexlag
       LAGGED       = LAGRO(I) + LAGGED 
       LAGROI(I1,I) = LAGRO(I) / 25.4
  950  continue                  

       TWE            =  WE+LIQW+STORAGE+LAGGED ! TOTAL WATER IN PACK

!      CALL AESNOW(WE,LIQW,ACCMAX,SIT,SB,SBWS,AESC,SBAESC,ADC,whichcurve) 

      else

       TWE = 0.0
       do 951 I=1,nexlag
       LAGROI(I1,I) = 0.0
  951  continue                  

      end if

      WEI(I1)    = WE      / 25.4
      NEGHSI(I1) = NEGHS   / 25.4
      LIQWI(I1)  = LIQW    / 25.4
      SBWSI(I1)  = SBWS    / 25.4
      ACUMX(I1)  = ACCMAX  / 25.4
      SBI(I1)    = SB      / 25.4
      AESCI(I1)  = AESC
      SBAESCI(I1)= SBAESC
      STOREI(I1) = STORAGE / 25.4

! END OF carryover (NO PUBLIC VARIABLES, saved parameters)
!---------------------------------------
!-Public variables

      TINDXI(I1)         = ( TINDEX * (9./5.) ) + 32.  ! convert and save temp of pack
      PKWATER_EQUIV(I1)  = TWE / 25.4
      SNOWCOV_AREA(I1)   = AESC
      CHNGPCKMRU(I1)     = (TWE - BEGINPACK) / 25.4

!-SET basin values
      BASIN_SNOWMELT= BASIN_SNOWMELT + PERAREA * SNOWMELT(I1)
      PSOILBASIN    = PSOILBASIN     + PERAREA * PSOILMRU(I1)
      COVERBASIN    = COVERBASIN     + PERAREA * SNOWCOV_AREA(I1)
c$$$      BASIN_SNOWEVAP= BASIN_SNOWEVAP + PERAREA * SNOW_EVAP(I1)
      BASIN_PWEQV   = BASIN_PWEQV    + PERAREA * TWE / 25.4
!-SET prior values and basin average snow and rain values
      NET_SNOW(I1)=DSFALL / 25.4  ! DSFALL=SNOW FOR THIS DAY
      NET_RAIN(I1)=DRAIN  / 25.4  ! DRAIN =RAIN ON SNOW AND OR BARE-GROUND
      IF (WE.EQ.0.0 .AND. DSFALL.GT.0.0 .AND. DRAIN.GT.0.0) 
     &      PPTMIX_NOPACK(I1)=1

c     If last of pack melted on this time step. Set estimated ET to zero - unknown author
c
c     4 Apr 09 - Comment out as water evaporated in previous time steps
c                goes unaccounted for. RMTW
c

c      if(TWE.eq.0) snow_evap(i1)=0.0

      BASIN_SNOWEVAP= BASIN_SNOWEVAP + PERAREA * SNOW_EVAP(I1)


C
C    Set basin average snow and rain values
C
      basin_net_snow = basin_net_snow + PERAREA * NET_SNOW(I1)
      basin_net_rain = basin_net_rain + PERAREA * NET_RAIN(I1)
      if(putvar('intcp', 'net_snow',  nmru, 'real', NET_SNOW)
     &   .ne.0) return
      if(putvar('intcp', 'net_rain',  nmru, 'real', NET_RAIN)
     &   .ne.0) return

c
c Debug
c
c$$$      if(i1.eq.3) then
c$$$         write (*,450) TWE*cm(i1),net_rain(i1)*ci(i1),
c$$$     $        net_snow(i1)*ci(i1),snowmelt(i1)*ci(i1),
c$$$     $        psoilmru(i1)*ci(i1),snow_evap(i1)*ci(i1),
c$$$     $        CHNGPCKMRU(I1)*ci(i1)
c$$$! 
c$$$      end if

1000  continue                                                     ! END OF MRU LOOP 


      DO 1010 I = 1,nmru
 1010   SUBRATE(I) = SUBRATE(I) / 25.4 ! RESET SUBLIMATION TO INCHES

      
!----------------------------------------------------*-*--*-*--*-*--*-*--*-*--*-*--
      nwsmrun = 0
!
      RETURN                                              ! RETURN TO MODEL
 450  format('snow',/,7(1x,f10.3))
      END

      SUBROUTINE MELTFACT(IDN,ALAT,MFMAX,MFMIN,MF,RATIO)
!..................................................................................
!     SUBROUTINE COMPUTES SURFACE MELT BASED ON 100 PERCENT
!        SNOW COVER AND NON-RAIN CONDITIONS.
!.......................................
!     INITIALLY WRITTEN BY...                           ----------------------
!        ERIC ANDERSON - HRL    MAY 1980, MODIFIED TO DO ONLY MELT COEFFICIENT
!                                                       ----------------------
!                           in MARCH 1998 BY J.J.VACCARO
!.....................................................................................
      real    ALAT,  MFMAX, MFMIN, MF, RATIO
      real    DIFF,  DAYN,  XX,    X,  ADJMF
      integer IDN
!.......................................
!  INITIAL VALUES
      DIFF = MFMAX-MFMIN
      DAYN = float(IDN)
      IF(ALAT.LT.54.0) GO TO 120  ! melt factor for latitudes < 54 N
!.......................................
!  MELT FACTOR VARIATION FOR ALASKA.
      IF(IDN.GE.275) GO TO 102
      IF(IDN.GE.92) GO TO 101
      X=(91.0+DAYN)/183.0
      GO TO 105
  101 X=(275.0-DAYN)/(275.0-92.0)
      GO TO 105
  102 X=(DAYN-275.0)/(458.0-275.0)
  105 XX=(SIN(DAYN*2.0*3.1416/366.0)*0.5)+0.5
      IF(X.LE.0.48) GO TO 111
      IF(X.GE.0.70) GO TO 112
      ADJMF=(X-0.48)/(0.70-0.48)
      GO TO 110
  111 ADJMF=0.0
      GO TO 110

  112 ADJMF=1.0
  110 MF=(XX*ADJMF)*DIFF+MFMIN
      GO TO 125
!.......................................
!  MELT FACTOR VARIATION FOR THE LOWER 48.

  120 MF=(SIN(DAYN*2.0*3.1416/366.0)*DIFF*0.5)+(MFMAX+MFMIN)*0.5
  125 RATIO=MF/MFMAX
!
      RETURN
      END 

      SUBROUTINE MELTRATE(TAIR,MBASE,NMINDEX,TINDEX,NMRATE,RATIO,NMFT,
     &                    CNHS,TSPM,RAIN,RFMIN,MELT,MF,RAINM,
     &                    SBCI,UADJT,PA,EVP,SUBR,ISTEP,MRU)
!------------------------------------------------------------------------------
! REST OF ANDERSON'S MELT19 "PLUS" ADDITIONAL STUFF THAT ALL DEALS W/ MELT
!                                  MAKES MORE SENSE TO KEEP AS ONE
!
      USE WEBMOD_SNOW, ONLY : nmru

      integer ISTEP, MRU
      real TAIR,    MBASE, TMX,   TSUR,  EA, TAK, TAK4
      real QN,      QE,    QH
      real NMINDEX, TINDEX,NMRATE,RATIO, MF, NMFT, TSPM
      real CNHS,    RAIN,  RFMIN,  MELT, RAINM
      real SBCI,    UADJT, PA
      real EVP,     SUBR(nmru)

      TMX = TAIR - MBASE                      ! Air temp as index of melt using MBASE
      if (TMX  .LT. 0.0) TMX = 0.0
      TSUR= TAIR
      if (TSUR .GT. 0.0) TSUR = 0.0

      NMINDEX = TINDEX - TSUR                  ! NEGATIVE HEAT EXCHANGE INDEX
!                                              ! TINDEX ALWAYS <= 0 SINCE TEMP OF SNOWPACK
      NMRATE = RATIO  * NMFT                   ! NEGATIVE MELT FACTOR ASSUMED TO VARY AS 'MF'
      CNHS   = NMRATE * NMINDEX
      TINDEX = TINDEX + TSPM * (TAIR-TINDEX)   ! CALTE NEW ANTECDNT TEMP INDEX FOR NEXT TIME
!---------------------------------------------------------------------------------
      if (RAIN .LE. RFMIN)  THEN            ! rain >  rmin/hr--use rain-on-snow  equation 
!                           ****            ! rain =< rmin/hr--use light/no rain equation
!-Melt during non-rain or light rain periods, AND ADD 'RAINM' FROM EARLIER

        MELT = MF * TMX + RAINM
!-Calculate sublimation rate for 6am-6pm of non-rain days, will be =0 if subrate set =0

        if (ISTEP.GT.1. AND. ISTEP.LT.4) EVP=SUBR(MRU) * 0.5

      ELSE                 ! ELSE DO THE RAIN ON SNOW EQUATION
!     ****
!-Melt during rain > rmin (0.25 mm * delt), designed for six hours

       EA  =2.7489E8 * EXP(-4278.63 / (TAIR+242.792) )
!                 ASSUME 90 PERCENT RELATIVE HUMIDITY DURING RAIN-ON-SNOW
       EA  =0.90*EA
       TAK  =(TAIR+273.) * 0.01
       TAK4=TAK*TAK*TAK*TAK
       QN  =SBCI * (TAK4-55.55)
       QE  =8.5*(EA-6.11) * UADJT
!      IF(IFUT.EQ.0) QE=QE*WINDC             IFUT USED AGAIN, COMMENTED OUT, NO WINDC FACTOR
!       QH  =7.5 * 0.000646 * PA * UADJT * TAIR
       QH  = 8.5 * 0.00057 * PA * UADJT * TAIR  ! changed back to real value of psychrometric constant of 5.7e-4 mb/degC, same result - RW
!      IF(IFUT.EQ.0) QH=QH*WINDC
! ADD COMPONENTS
       MELT = QN + QE + QH + RAINM
       if (MELT .LT. 0.0)  MELT = 0.0

      END IF
!     ******
      RETURN
      END


      SUBROUTINE ROUTE(IT,EXCESS,WE,AESC,STORGE,EXLAG,PACKRO,NEXLAG)
!...............................................................................
!     THIS SUBROUTINE ROUTES EXCESS WATER THROUGH THE SNOW COVER FOR
!     SUBROUTINE INITIALLY WRITTEN BY...
!        ERIC ANDERSON - HRL   MAY 1980, added to MMS by J. Vaccaro 3/98
!       accounts for changes in other information
!       changed exlag to dimension=2 as far as can see L2=2, L1=1, for most any cases
!......................................................................................
      integer NEXLAG
      real EXLAG(NEXLAG)
      integer IT, N, I, L1, L2
      real EXCESS, WE, AESC, STORGE, PACKRO, PRIOR, WES
      real FN, FI, FIT, CL, R1, ENDL1
      real TERM, FLAG, POR1,POR2, EL, ELS, OS
!.......................................
!     INITIAL VALUES
      FIT   = float(IT)                 ! floating point of time-step
      PACKRO= 0.0                       ! set pack runoff =0.0
      CL    = 0.03*FIT/6.0              ! factor, =0.03 for DELT=6hr, =0.005 for DELT=1hr
      PRIOR = STORGE + EXLAG(1)
!.......................................
!   lag excess water, function of amount of excess and water equiv (we)
!                   ! check if no new excess and PRIOR water is =  0.0 (go to 190)
      if (EXCESS .EQ. 0.0 .and. PRIOR .EQ. 0.0) GO TO 190
!                   ! check if no new excess and PRIOR water is => 0.1 (go to 160)
      if (EXCESS .EQ. 0.0 .and. PRIOR .GE. 0.1) GO TO 160
!                   ! check  < 0.1 mm EXCESS or WE < 1 mm, don't lag,  (go to 120)
      if (EXCESS .LT .0.1  .or. WE .LT. 1.0)  GO TO 120  
!                                               
!  COMPUTE LAG IN HOURS AND PRORATE EXCESS.
      N = ( (EXCESS*4.0)**0.3 ) + 0.5  ! n varies from 2--low excs, to about 8--250 mm exces
      if(N.EQ.0) N=1
      FN = float(N)
      DO 110 I=1,N
      FI = float(I)
      TERM=CL*WE*FN/(EXCESS*(FI-0.5))
      if (TERM .GT. 150. 0) TERM=150.0
      FLAG = 5.33 * (1.0-EXP(-TERM))   ! Lag of the excess liquid water:float of lag in hours
      L2=(FLAG+FIT)/FIT+1.0            ! always 2, next always=1 so WHY DO
      L1=L2-1    
      ENDL1=L1*IT                     ! always=it, so also WHY DO
      POR2     = (FLAG+FIT-ENDL1)/FIT
      POR1     = 1.0-POR2
      EXLAG(L2) = EXLAG(L2) + POR2*EXCESS/FN
      EXLAG(L1) = EXLAG(L1) + POR1*EXCESS/FN
  110 CONTINUE

      GO TO 150
!
!  EXCESS OR WE SMALL, THUS NO LAG.
  120 EXLAG(1)=EXLAG(1)+EXCESS
!.....................................................................................
! CHECK TO SEE AMOUNT OF WATER TO ATTENUATE LAGGED EXCESS WATER - FUNCTION OF STORGE AND WE.

  150 if( (STORGE+EXLAG(1)) .EQ. 0.0) GO TO 190   ! no H2O to attenuate or RunOff
      if( (STORGE+EXLAG(1)) .GE. 0.1) GO TO 160   ! H2O, enough to attenuate
!
!   NO ATTENUATION
      PACKRO=STORGE+EXLAG(1)
      STORGE=0.0
      GO TO 190
!--------------------------------------------------------------------------------------
!--ATTENUATE LAGGED EXCESS WATER - FUNCTION OF STORGE AND WE.
!  EFFECT OF ATTENUATION COMPUTED USING A ONE-HOUR TIME STEP, mm convrted to inches
!                                                             to calculate R1 factor

  160 EL  = EXLAG(1)/FIT
      ELS = EL/(25.4*AESC)             
      WES = WE/(25.4*AESC)             
      TERM= 500.0*ELS/(WES**1.3)       
      if(TERM.GT.150.0) TERM=150.0
      R1=1.0/(5.0*EXP(-TERM)+1.0)              ! R1 is almost always 1.0, why bother

      DO 170 I=1,IT
      OS=(STORGE+EL)*R1
      PACKRO=PACKRO+OS
      STORGE=STORGE+EL-OS
  170 CONTINUE

      if (STORGE .GT .0.001) GO TO 190
      PACKRO=PACKRO+STORGE
      STORGE=0.0
!
!--Downshift water in EXLAG().
  190 DO 195 I=2,NEXLAG
  195 EXLAG(I-1)   = EXLAG(I)
      EXLAG(NEXLAG)= 0.0
!.................................................................................
      RETURN
      END 
      
      SUBROUTINE AESNOW(WE,LIQW,ACCMAX,SIT,SB,SBWS,AESC,
     &                  SBAESC,ADC,ICRVE)
!................................................................................
!     THIS SUBROUTINE COMPUTES THE AREAL EXTENT OF SNOW COVER USING THE
!        AREAL DEPLETION CURVE FOR THE 'SNOW-17 ' OPERATION.
!.......................................
!     SUBROUTINE INITIALLY WRITTEN BY...
!        ERIC ANDERSON - HRL   MAY 1980
!       slightly modified by J. Vaccaro, March, 1998
!.................................................................................

      USE WEBMOD_SNOW, ONLY : ndepl
      integer ICRVE, N
      real    WE, LIQW, ACCMAX,  SIT, SB, SBWS, AESC, SBAESC
      real    TWE, FN,  R,   AI
      real    ADC(11,ndepl)
! 
!-------------------------------------------------------------------------------

      TWE = WE + LIQW              ! LIQW AS INITIAL, WE may have been INCREASED DUE TO SFALL
! 
      if (TWE .GT. ACCMAX) ACCMAX=TWE  ! TOTAL WATER > LAST MAX ACUMLTED, RESET MAX as larger
      AI = ACCMAX                      ! Set checking variable to max SNOW accum FOR SEASON
      if (ACCMAX .GT. SIT)  AI=SIT     ! If max accum  > H2O equiv when always 100% COVER
!                                        reset AI to smaller VALUE of SIT

      if (TWE    .GE. AI)   go to 202  ! we and liquid => AI, 100% cover, set TWE prior to new snow
      if (TWE    .LE. SB)   go to 201  ! we and liquid <= water prior to new snow, recompute 
      if (TWE    .GE. SBWS) go to 203  ! we and liquid => water after new snow, set cover=100
!                                                                           and don't reset SB
      AESC = SBAESC + ((1.0 - SBAESC)*((TWE-SB)/(SBWS-SB))) !ON SAME CURVE, <100%, CONTINUE ON 
! RMTW copied these 3 lines from the 201 section to reset after melt below SBWS
      SB   = TWE + 1.27                     ! UPDATE SB W/ NEW POINT, ANDERSON USED 'SNOF'
      SBWS = TWE                            ! USING 1.27 MM HERE (NOTE THAT SNOF ALSO USED
      SBAESC = AESC                         ! IN FIRST PART FOR SNOWFALL, AND WE'RE NOT USING)
!                                                      CURVE CANNOT GO TO 100% BASED ON CHECKS ABOVE
      go to 204

  201 R  = (TWE/AI)*10.0+1.0           ! AREAL EXTENT SMALLER, RECOMPUTE, FIND WHICH POINT (R=N)
      N  = R
      FN = N
      R  = R - FN
      AESC = ADC(N,ICRVE) + (ADC(N+1,ICRVE)-ADC(N,ICRVE))*R  ! NEW AREAL EXTENT, CHECK IF NOW 100%
!      if (AESC .GT. 1.00) AESC = 1.00           ! This is checked with the jump to 204

      SB   = TWE + 1.27                     ! UPDATE SB W/ NEW POINT, ANDERSON USED 'SNOF'
      SBWS = TWE                            ! USING 1.27 MM HERE (NOTE THAT SNOF ALSO USED
      SBAESC = AESC                         ! IN FIRST PART FOR SNOWFALL, AND WE'RE NOT USING)

      go to 204        

  202 SB   = TWE                                  ! reset SB, SBWS, SINCE INCREASE 
      SBWS = TWE                                  ! SBWS here since now includes all SFALL
  203 AESC = 1.0
  204 if (AESC .LT. 0.05) AESC = 0.05
      if (AESC .GT. 1.00) AESC = 1.00

      RETURN
      END
      

      SUBROUTINE HEATWATR(MELT,RAIN,CNHS,CNHSPX,WCPLW,WE,NEGHS,
     &                    EXCESS,LIQW,TINDEX)
!-----------------------------------------------------------------------------------
!
! SUBROUTINE DOES THE HEAT AND WATER BALANCE FOR THE SNWPACK FOR THIS MRU FOR THIS
! TIME-STEP. EXCESS, AMOUNT OF WATER PACK CANN'T HOLD IS ALSO CALCULATED AND VALUES
! FOR NEGATIVE HEAT (NEGHS), LIQUID WATER (LIQW), WATER-EQUIVALENT (WE) ARE UPDATED
! put this as separate subroutine to make more readable (main was way to long)
!----------------------------------------------------------------------------------
      real LIQWMX, TINDEX, MELT, WCPLW, RAIN, CNHS, WE
      real EXCESS, CNHSPX, LIQW, WATER, HEAT, NEGHS

      WATER  = MELT  + RAIN        ! water is set to melt + rain falling on snow
      HEAT   = CNHS  + CNHSPX
      LIQWMX = WCPLW * WE       ! WE TIMES PERCENT LIQUID WATER HOLDING CAPACITY--MAX LIQD H2O
      NEGHS  = NEGHS + HEAT     ! UPDATE NEGHS TO INCLUDE HEAT FROM PRECIP, MELT (THESE ARE '-')

      if (NEGHS .LT. 0.0    ) NEGHS=0.0 !CNHS AND CNHSPX BROUGHT NEGHS TO -, SET=0, OTHERWISE '+'
      if (NEGHS. GT. 0.33*WE) NEGHS=0.33*WE   ! NEGATIVE HEAT CANN'T EXCEED 33% OF WE

! need to check if the water (old liquid + newmelt and rain) exceeds capacity

      if ((WATER+LIQW).LT.(LIQWMX+NEGHS+WCPLW*NEGHS)) GO TO 230

!--There is excess water
      EXCESS = WATER + LIQW - LIQWMX - NEGHS  - WCPLW*NEGHS  !FIND EXCES,RESET VALUES

      LIQW   = LIQWMX +  WCPLW*NEGHS  ! SET LIQUID TO MAX CAPCITY, THEN PUT THE NEG HEAT INTO THE 
      WE    = WE     +  NEGHS        ! WE, SO HAVE SET THESE AND NEGHS IS ADDED TO WE SO SET=0.0
      NEGHS  = 0.0

      go to 232

!--There is NO excess water
  230 if (WATER .LT. NEGHS) go to 231

      LIQW   = LIQW + WATER - NEGHS          ! water > neghs, increase liquid by amnt above neghs
      WE     = WE   + NEGHS
      NEGHS  = 0.0
      EXCESS = 0.0

      go to 232

!--There is NO excess water AND water total is smaller than the negative heat
  231 WE = WE + WATER                          ! water <= neghs, SO REFREEZE THE WATER
!                                                reduce NEGHS by WATER refrozen
      NEGHS = NEGHS-WATER
      EXCESS = 0.0

  232 continue

      if (TINDEX .GT. 0.0) TINDEX = 0.0
      if (NEGHS  .EQ. 0.0) TINDEX = 0.0

!--FINISHED HEAT AND WATER BALANCE FOR THIS DELT----------------------------------

      RETURN
      END

