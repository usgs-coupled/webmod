#include "defines.h"
!***********************a************************************************
!     MMS module that replaces call_modules.c and setdims.f
!***********************************************************************
      MODULE WEB_MODULE
      IMPLICIT NONE

      CHARACTER(LEN=32) :: Precip_module, Temp_module, Et_module
      CHARACTER(LEN=32) :: Srunoff_module, Solrad_module, Soltab_module
      CHARACTER(LEN=32) :: Summary_module
      INTEGER :: Model, Climate_dist
      INTEGER :: ID
!   Control Parameters
      CHARACTER(LEN=9) :: Model_mode
      END MODULE WEB_MODULE

!***********************************************************************
!***********************************************************************
      INTEGER FUNCTION call_modules(Arg)
      USE WEB_MODULE
      IMPLICIT NONE
      INCLUDE 'fmodules.inc'
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
! Functions (define in calling order)
      INTEGER, EXTERNAL :: io
      INTEGER, EXTERNAL :: soltab_prms
      INTEGER, EXTERNAL :: basin_topg
      INTEGER, EXTERNAL :: obs_chem
      INTEGER, EXTERNAL :: obs_webmod
      INTEGER, EXTERNAL :: obs_webx
      INTEGER, EXTERNAL :: temp_1sta_prms
      INTEGER, EXTERNAL :: precip_web
      INTEGER, EXTERNAL :: xyz_dist_web
      INTEGER, EXTERNAL :: irrig_web
      INTEGER, EXTERNAL :: ccsolrad_web
      INTEGER, EXTERNAL :: potet_hamon_prms
      INTEGER, EXTERNAL :: intcp_prms
      INTEGER, EXTERNAL :: nwsmelt_topg
      INTEGER, EXTERNAL :: topmod_chem
      INTEGER, EXTERNAL :: top2clark
      INTEGER, EXTERNAL :: route_clark
      INTEGER, EXTERNAL :: webmod_res
      INTEGER, EXTERNAL :: phreeq_mms
      INTEGER, EXTERNAL :: web_sum
      !DOUBLE PRECISION :: dt
      !INTEGER :: storm, runflg
      INTEGER :: runflg
!***********************************************************************
      call_modules = 1
      runflg = 0

      IF ( Arg.EQ.'declare' ) THEN


        PRINT *, 'The following modules are included in WEBMOD_1.0:'
        PRINT *, ' io.f, soltab_prms, basin_topg, obs_chem, obs_webmod,'
        PRINT *, ' obs_webx, temp_1sta_prms, precip_web,xyz_dist_dist,'
        PRINT *, ' irrig_web, ccsolrad_web, potet_hamon_prms'
        PRINT *, ' intcp_prms, nwsmelt_topg, topmod_chem '
        PRINT *, ' top2clark, route_clark, webmod_res, '
        PRINT *, ' phreeq_mms, web_sum.f'

        call_modules = declmodule(
     +'$Id: call_modules_webmod.f 494 2013-07-15 18:00:00Z rmwebb $')

        IF ( control_string(Model_mode, 'model_mode').NE.0 ) RETURN

!       Model (1=daily; 2=storm; 99=reserved for name files)',
!       Climate_dist (0=dist2, ide_dist; 1=1sta, laps, 2sta; 2=xyz)

        IF ( Model_mode(:5).EQ.'DAILY' .OR. Model_mode.EQ.'        ')
     +       THEN
          Model = 1
        ELSEIF ( Model_mode(:5).EQ.'STORM' ) THEN
          Model = 2
        ELSEIF ( Model_mode(:5).EQ.'TRENT' ) THEN
          Model = 3
        ELSEIF ( Model_mode(:4).EQ.'CFGI' ) THEN
          Model = 4
        ELSEIF ( Model_mode(:5).EQ.'_name' ) THEN
          Model = 99
        ELSE
          Model = -1
          PRINT 9002, Model_mode, Arg
          RETURN
        ENDIF

        Climate_dist = 1
        IF ( control_string(Precip_module, 'precip_module')
     +       .NE.0 ) RETURN
        IF ( Precip_module(:10).NE.'precip_web' .AND.
     +       Precip_module(:12).NE.'xyz_dist_web' ) THEN
          PRINT *, 'WARNING: Invalid precip_module value, reset to:',
     +             ' precip_web, value was: ', Precip_module //
     +             ' Temp_module set to temp_1sta_prms.'
        ELSE
          Precip_module = 'precip_web'
          Temp_module = 'temp_1sta_prms'
        ENDIF
        IF ( Precip_module(:12)=='xyz_dist_web' ) Climate_dist = 2

! if xyz is specified for one, it must be specified for the other
        IF ( Temp_module(:12).EQ.'xyz_dist_web' ) THEN
          Precip_module = Temp_module
        ELSEIF ( Precip_module(:12).EQ.'xyz_dist_web' ) THEN
          Temp_module = Precip_module
        ENDIF
      ENDIF

      IF ( Model==99 ) THEN
        IF ( Arg/='declare' ) THEN
          call_modules = 1
          RETURN
        ENDIF
        call_modules = io(Arg)
        call_modules = soltab_prms(Arg)
        call_modules = basin_topg(Arg)
        call_modules = obs_chem(Arg)
        call_modules = obs_webmod(Arg)
        call_modules = temp_1sta_prms(Arg)
        call_modules = precip_web(Arg)
        call_modules = obs_webx(Arg)
        call_modules = xyz_dist_web(Arg)
        call_modules = irrig_web(Arg)
        call_modules = ccsolrad_web(Arg)
        call_modules = potet_hamon_prms(Arg)
        call_modules = intcp_prms(Arg)
        call_modules = nwsmelt_topg(Arg)
        call_modules = topmod_chem(Arg)
        call_modules = top2clark(Arg)
        call_modules = route_clark(Arg)
        call_modules = webmod_res(Arg)
        call_modules = phreeq_mms(Arg)
        call_modules = web_sum(Arg)
        call_modules = 0
        WRITE (*, 9003)
        RETURN
      ENDIF

! All modules must be called during declare for model modes 1 & 2
      call_modules = io(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'io', Arg, call_modules
        RETURN
      ENDIF

      call_modules = soltab_prms(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'soltab_prms', Arg, call_modules
        RETURN
      ENDIF

      call_modules = basin_topg(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'basin_topg', Arg, call_modules
        RETURN
      ENDIF

      call_modules = obs_chem(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'obs_chem', Arg, call_modules
        RETURN
      ENDIF

      IF(Climate_dist.EQ.1) THEN
            call_modules = obs_webmod(Arg)
            IF ( call_modules.NE.0 ) THEN
             PRINT 9001, 'obs_webmod', Arg, call_modules
             RETURN
            ENDIF
            call_modules = temp_1sta_prms(Arg)
            IF ( call_modules.NE.0 ) THEN
             PRINT 9001, 'temp_1sta_prms', Arg, call_modules
             RETURN
            ENDIF
            call_modules = precip_web(Arg)
            IF ( call_modules.NE.0 ) THEN
             PRINT 9001, 'precip_web', Arg, call_modules
             RETURN
            ENDIF
      ELSEIF(Climate_dist.EQ.2) THEN
            call_modules = obs_webx(Arg)
            IF ( call_modules.NE.0 ) THEN
             PRINT 9001, 'obs_webx', Arg, call_modules
             RETURN
            ENDIF
            call_modules = xyz_dist_web(Arg)
            IF ( call_modules.NE.0 ) THEN
             PRINT 9001, 'xyz_dist_web', Arg, call_modules
             RETURN
            ENDIF
      ELSE
            PRINT 9004, Climate_dist
            RETURN
      ENDIF

      call_modules = irrig_web(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'irrig_web', Arg, call_modules
        RETURN
      ENDIF

      call_modules = ccsolrad_web(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'ccsolrad_web', Arg, call_modules
        RETURN
      ENDIF

      call_modules = potet_hamon_prms(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'potet_hamon_prms', Arg, call_modules
        RETURN
      ENDIF

      call_modules = intcp_prms(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'intcp_prms', Arg, call_modules
        RETURN
      ENDIF

      call_modules = nwsmelt_topg(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'nwsmelt_topg', Arg, call_modules
        RETURN
      ENDIF

      call_modules = topmod_chem(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'topmod_chem', Arg, call_modules
        RETURN
      ENDIF

      call_modules = top2clark(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'top2clark', Arg, call_modules
        RETURN
      ENDIF

      call_modules = route_clark(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'route_clark', Arg, call_modules
        RETURN
      ENDIF

      call_modules = webmod_res(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'webmod_res', Arg, call_modules
        RETURN
      ENDIF

      call_modules = phreeq_mms(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'phreeq_mms', Arg, call_modules
        RETURN
      ENDIF

      call_modules = web_sum(Arg)
      IF ( call_modules.NE.0 ) THEN
        PRINT 9001, 'web_sum', Arg, call_modules
        RETURN
      ENDIF

      call_modules = 0

 9001 FORMAT ('Error in ', A, ' module, arg = ', A, /, 'Return val =',
     +        I4)
 9002 FORMAT ('Error in call_modules_webmod: model_mode', A,
     +        ' not implemented', /, A)
 9003 FORMAT (//, ' _name files were generated and written to files', /,
     +        ' in the directory containing the control file.', /,
     +        ' Ignore any messages that follow this one.', /,
     +        ' Note, no simulation was computed.', /)
 9004 FORMAT ('Error in call_modules_webmod: Climate_mode', I4,
     +        ' not implemented')

      END FUNCTION call_modules

!***********************************************************************
!     declare the dimensions
!***********************************************************************

      INTEGER FUNCTION setdims()
      IMPLICIT NONE
      INCLUDE 'fmodules.inc'
! Local Variables
      INTEGER :: ncurves
!***********************************************************************
      setdims = 1

      if(declfix('five', 5, 5, 'Dimension of obj_q').ne.0) return

      if(declfix('nresinp', 21, MAXRESINP, 'Vector of reservoir '//
     $     'volumes used to calculate mixing ratios').ne.0) return

      if(declfix('nmru_res', 9, MAXMRURES, 'Vector of hillslope '//
     $   'reservoirs for initial solutions and reactants').ne.0)return

      if(declfix('nchemvar', 10, 10,
     $   'Number of available user-defined outputs to track '//
     $   'concentrations in any of the model reservoirs').ne.0) return

      if(declfix('nconvert', 3, 3,
     $     'Number of user-defined factors to convert arbitrary '//
     $     'concentrations to Moles per Liter').ne.0) return

      if(decldim('nchem_sets', 10, MAXSOLNSETS,
     $     'Number of unique solution and/or reaction entities '//
     $     'describing variety of typical hillslopes').ne.0) return

      if(decldim('nac_nmru_nresinp', MAXNMR_3D, MAXNMR_3D,
     $     'Vector of volumes used to calculate mixing ratios needed '//
     $     'for simulating new topographic bin solutions').ne.0) return

c
c nhydro will be used in the chemical mixing modules to limit the
c number of unique segments in the drainage system. For the Clarke
c routing module, nhydro equals the number of routing delay bins, nrtdelay,
c which varies according to distance from the outlet, stream velocities,
c and the time step. A check that nrtdelay and nhydro are equal is done in
c route_clark. In the mixing modules, the variable vmix_stream
c will be declared of dimension 'nhydro,nchan'. Lakes can be added
c later once a more realistic routing routine is implemented.
c
      if(decldim('nhydro', 1, 1,
     +   'Number of time-delay ordinates for stream '//
     $   'routing').ne.0) return

      if(decldim('nmru', 100, MAXMRU,
     +   'Number of model response units').ne.0) return

      if(decldim('nobs', 1, MAXOBS,
     +   'Number of stream measurement stations').ne.0) return

      if(decldim('nrain', 1, MAXRAIN,
     +   'Number of rain gages').ne.0) return

      if(decldim('ntemp', 1, MAXTEMP,
     +   'Number of temperature stations').ne.0) return

      if(decldim('nsol', 1, MAXSOL,
     +   'Number of solar radiation stations').ne.0) return

      if(decldim('nsolute', 1, MAXSOLUTE,
     +   'Number of solute species to be tracked. Each of nsolute '//
     +   'must be named as in the 3rd field in the phreeq_lut file.'
     +   ).ne.0) return

      if(decldim('nphq_lut', 1, MAXSOLUTE,
     +   'Number of species in the phreeq_lut file').ne.0) return

      if(decldim('nchemobs', 1, MAXCHEMOBS,
     +     'Number of unique chemical sampling sites').ne.0) return

      if(decldim('nchem_ext', 1, MAXCHEMOBS,
     +   'Number of chemically unique external '//
     $   'sources').ne.0) return

      if(decldim('ngw_ext', 1, MAXIRRIG,
     +   'Number of unique flux rates from external '//
     $   'ground water sources').ne.0) return

      if(decldim('nirrig_ext', 1, MAXIRRIG,
     +   'Number of unique application rates from external '//
     $   'irrigation sources').ne.0) return

      if(decldim('nirrig_int', 1, MAXIRRIG,
     +   'Number of unique application rates from internal '//
     $   'irrigation sources').ne.0) return
!
!      if(decldim('ngnd', 1, MAXGND,
!     +   'Number of soil heat fluxstations').ne.0) return
!
!      if(decldim('nwind', 1, MAXWIND,
!     +   'Number of wind velocity stations').ne.0) return

      if(decldim('nhum', 1, MAXHUM,
     +   'Number of relative humidity stations').ne.0) return
!
!      if(decldim('npres', 1, MAXPRES,
!     +   'Number of atmospheric pressure stations').ne.0) return

      if(decldim('nevap', 1, MAXEVAP,
     +   'Number of pan evaporation data sets').ne.0) return

      if(decldim('nsnow', 1, MAXSNOPIL,
     +   'Number of snow measurement stations').ne.0) return

      ncurves = 2 ! number of snow depletion curves
      
      if(decldim('ndepl', ncurves, ncurves,
     +   'Number of snow depletion curves').ne.0) return

      if(declfix('nmonths', 12, MAXMO,
     +   'Number of months in a year').ne.0) return

      if(decldim('ndeplval', ncurves*11, ncurves*11,
     +   'ndepl * 11').ne.0) return

      if(decldim('ntopchan', 5, MAXTOPCHAN, 
     +   'Number of channel segments associated with a '//
     +   'Clarke hydrograph channel (TOPMODEL)').ne.0) return

      if(decldim('nchan', 1, MAXCHAN, 
     +   'Generic number of channel segments for flow '//
     +   'routing').ne.0) return

      if(declfix('ndays', 366, MAXDAY,
     +   'Number of days in a year').ne.0) return

      if(decldim('nac', 30, MAXNAC,
     +   'Number of TTI bins').ne.0)
     +    return

      if(decldim('nform', 1, MAXFORM,
     +   'Number of form_data input values. This has a max of 1 but '//
     +   'is an array to satisy a readvar requirement. Set to 0 if'//
     +   'no form_data are included in the input file.').ne.0) return

      if(decldim('nexlag', 1, MAXLAG,
     +   'Number of excess liquid lag increments').ne.0)
     +    return

      if(declfix('nxkbin', 9, MAXXKBIN,
     +   'Number of bins to calculate infiltration excess.').ne.0)
     +    return

      if(declfix('nlapse', 3, MAXLAPSE,
     +   'Number of lapse rates in X, Y, and Z directions.')
     +    .ne.0) return


      setdims = 0
      END FUNCTION setdims
