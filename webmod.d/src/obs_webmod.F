#include "defines.h"
c***********************************************************************
c     obs_webmod.f: modified from obs_prms. -RMTW
c         4 Mar 05 - changed tmin to tsta_min_c and _f.
c                    Same with tmax and tstemp.
c **********************************************************************
c                        added temp_units parameter
c                        Snow water equivalence
c                        and irrigation schedules,
c                        and ground water influx, - earlier
c                        and relative humidity (for isotopes) - 24feb11
c
c                        added RCS version control,
c    27 apr  2010 - Port to Fortran 90 with module and dynamic memory
c
c        08jul16 - Removed nform, form_data, and strain as no storms and 
c                  no observed form for each MRU
c **********************************************************************
      MODULE WEBMOD_OBSHYD
      IMPLICIT NONE
#include "fmodules.inc"

      integer, save :: temp_units, nobs, nsol
c      integer, save :: ntemp, nrain, nsnow, nform, ngw_ext, nevap
      integer, save :: ntemp, nrain, nsnow, ngw_ext, nevap
      integer, save :: nirrig_ext, nirrig_int, nhum, nsolute
      logical, save :: rh_tf, sh, warned_rh
c      integer, save :: form_data(1), datetime(6)
      integer, save :: datetime(6)

      real, save, allocatable :: runoff(:)
      real, save, allocatable :: precip(:)
      real, save, allocatable :: irrig_ext(:)
      real, save, allocatable :: irrig_int_next(:)
      real, save, allocatable :: gw_ext(:)
      real, save, allocatable :: tsta_min_c(:), tsta_max_c(:)
      real, save, allocatable :: tsta_min_f(:), tsta_max_f(:)
      real, save, allocatable :: tsta_temp_c(:), tsta_temp_f(:)
      real, save, allocatable :: solrad(:)
      real, save, allocatable :: swe(:), relhum(:), spechum(:)
      real, save, allocatable :: pan_evap(:)
            
      END MODULE WEBMOD_OBSHYD

c **********************************************************************
c
c     main obs routine
c

      integer function obs_webmod(arg)
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      integer obsdecl, obsinit, obsrun
      save SVN_ID

      SVN_ID = 
     $     '$Id$ '

      obs_webmod = 0

      if (arg.eq.'declare') then
         obs_webmod = obsdecl()

      else if (arg.eq.'initialize') then
         obs_webmod = obsinit()

      else if (arg.eq.'run') then
         obs_webmod = obsrun()

      end if

      return
      end
c
c **********************************************************************
c
c     obsdecl - makes public variable declarations for the
c                     obs module
c
 
      integer function obsdecl()

      USE WEBMOD_OBSHYD
      IMPLICIT NONE

      obsdecl = 1
      
!
! Get dimensions
!
c      nform = getdim('nform')
c        if ( nform.eq.-1 ) return
      nobs = getdim('nobs')
        if ( nobs.eq.-1 ) return
      nsol = getdim('nsol')
        if ( nsol.eq.-1 ) return
      nevap = getdim('nevap')
        if ( nevap.eq.-1 ) return
      nrain = getdim('nrain')
        if ( nrain.eq.-1 ) return
      nsnow = getdim('nsnow')
        if ( nsnow.eq.-1 ) return
      ntemp = getdim('ntemp')
        if ( ntemp.eq.-1 ) return
      ngw_ext = getdim('ngw_ext')
        if ( ngw_ext.eq.-1 ) return
      nirrig_int = getdim('nirrig_int')
        if ( nirrig_int.eq.-1 ) return
      nirrig_ext = getdim('nirrig_ext')
        if ( nirrig_ext.eq.-1 ) return
      nsolute = getdim('nsolute')
        if ( nsolute.eq.-1 ) return
      nhum = getdim('nhum')
        if ( nhum.eq.-1 ) return
c
c Relative humidity is required for computations of wetbulb temperature and isotopic 
c fractionation in phreeq_mms (chem_sim=1, nsolute.ne.0). The value is derived from
c either observations at mru_rhsta (should point to same mru_tsta) or specific
c humidity available in grids from climate models including the Idaho data set.
c
      if (nsolute.gt.0) then
        ALLOCATE(relhum(nhum))
        if(declvar('obs', 'relhum', 'nhum', nhum, 'real',
     +   'Observed relative humidity',
     +   'fraction',
     +   relhum).ne.0) return
        ALLOCATE(spechum(nhum))
        if(declvar('obs', 'spechum', 'nhum', nhum, 'real',
     +   'Specific humidity',
     +   'ratio',
     +   spechum).ne.0) return
      endif

      if(declparam('obs', 'temp_units', 'one', 'integer',
     +   '1', '0', '1',
     +   'Units for observed temperature',
     +   'Units for observed temperature, 0=F, 1=C',
     +   'none').ne.0) return

      if(nobs.ne.0) then
        ALLOCATE(runoff(nobs))
        if(declvar('obs', 'runoff', 'nobs', nobs, 'real',
     +     'Observed runoff for each gage',
     +     'cfs',
     +   runoff).ne.0) return
      end if
      
      ALLOCATE(precip(nrain))
      if(declvar('obs', 'precip', 'nrain', nrain, 'real',
     +     'Observed precip at each rain gage',
     +     'inches',
     +   precip).ne.0) return

      if(ngw_ext.ne.0) then
        ALLOCATE(gw_ext(ngw_ext))
        if(declvar('obs', 'gw_ext', 'ngw_ext', ngw_ext, 'real',
     +  'Influx of groundwater from sources external to basin',
     $  'cfs/mile', gw_ext).ne.0) return
      endif

      if(nirrig_ext.gt.0) then
        ALLOCATE(irrig_ext(nirrig_ext))
        if(declvar('obs','irrig_ext','nirrig_ext',nirrig_ext,'real',
     +  'Depth of applied irrigation from an external source',
     $  'inches', irrig_ext).ne.0) return
      endif

      if(nirrig_int.gt.0) then
        ALLOCATE(irrig_int_next(nirrig_int))
        if(declvar('obs', 'irrig_int_next', 'nirrig_int', nirrig_int,
     $     'real', 'Depth of irrigation from an internal source '//
     $     'to apply the next day',
     $     'inches', irrig_int_next).ne.0) return
      end if

      ALLOCATE(tsta_min_c(ntemp))
      if(declvar('obs', 'tsta_min_c', 'ntemp', ntemp, 'real',
     +     'Observed daily minimum temperature',
     +     'degrees celsius',
     +   tsta_min_c).ne.0) return

      ALLOCATE(tsta_max_c(ntemp))
      if(declvar('obs', 'tsta_max_c', 'ntemp', ntemp, 'real',
     +     'Observed daily maximum temperature',
     +     'degrees celsius',
     +   tsta_max_c).ne.0) return

      ALLOCATE(tsta_temp_c(ntemp))
      if(declvar('obs', 'tsta_temp_c', 'ntemp', ntemp, 'real',
     +     'Average observed daily temperature',
     +     'degrees celsius',
     +   tsta_temp_c).ne.0) return

      ALLOCATE(tsta_min_f(ntemp))
      if(declvar('obs', 'tsta_min_f', 'ntemp', ntemp, 'real',
     +     'Observed daily minimum temperature',
     +     'degrees fahrenheit',
     +   tsta_min_f).ne.0) return

      ALLOCATE(tsta_max_f(ntemp))
      if(declvar('obs', 'tsta_max_f', 'ntemp', ntemp, 'real',
     +     'Observed daily maximum temperature',
     +     'degrees fahrenheit',
     +   tsta_max_f).ne.0) return

      ALLOCATE(tsta_temp_f(ntemp))
      if(declvar('obs', 'tsta_temp_f', 'ntemp', ntemp, 'real',
     +     'Average observed daily temperature',
     +     'degrees fahrenheit',
     +   tsta_temp_f).ne.0) return

      if(nsol.ne.0) then
      ALLOCATE(solrad(nsol))
      if(declvar('obs', 'solrad', 'nsol', nsol, 'real',
     +     'Observed solar radiation',
     +     'langleys',
     +   solrad).ne.0) return
      endif ! nsol.ne.0

      if(nevap.ne.0) then
       ALLOCATE(pan_evap(nevap))
       if(declvar('obs', 'pan_evap', 'nevap', nevap, 'real',
     +     'Observed pan evaporation',
     +     'inches',
     +    pan_evap).ne.0) return
      endif ! nevap.ne.0

      if(nsnow.ne.0) then
       ALLOCATE(swe(nsnow))
       if(declvar('obs', 'swe', 'nsnow', nsnow, 'real',
     +     'Observed snow water equivalent of snowpack',
     +     'inches',
     +    swe).ne.0) return
      endif ! nsnow.ne.0

c      if(declvar('obs', 'form_data', 'one', 1, 'integer',
c     +     'Form of precip - 0=not known, 1=snow, 2=rain',
c     +     'none',
c     +   form_data).ne.0) return

c      if(declvar('obs', 'route_on', 'one', 1, 'integer',
c     +     'Kinematic routing switch - 1=storm period',
c     +     'none',
c     +   route_on).ne.0) return
 
      obsdecl = 0

      return
      end

c **********************************************************************
c
c     obsinit - initializes obs module
c

      integer function obsinit()

      USE WEBMOD_OBSHYD
      IMPLICIT NONE

c      integer datetime(6)

      obsinit = 1

      if(getparam('obs', 'temp_units', 1,
     $     'integer',temp_units) .ne.0) return
c expect both relative humidity and specific humidity. Only one allowed if geochem. Establish which
c in the run section
      rh_tf=.true.
      sh=.true.
      warned_rh=.false.
c
c     print out start and end times
c
c      call dattim('start', datetime)
c      call dpint4('Start        :', datetime, 6, 2)
c      call dattim('end', datetime)
c      call dpint4('End          :', datetime, 6, 2)
c      call dpstr(' ', 2)

      obsinit = 0
 
      return
      end

c **********************************************************************
c
c     obsrun - runs obs module
c
 
      integer function obsrun()

      USE WEBMOD_OBSHYD
      IMPLICIT NONE

      integer i
c      integer nstep, datetime(6)
      !integer nstep
      !integer julcalen, julsolar, julwater
      double precision dt

      obsrun = 1

      dt = deltim()
      call dattim('now', datetime)
      if(nobs.ne.0) then
        if(readvar('obs','runoff').ne.0) return
      endif
      if(readvar('obs','precip').ne.0) return
      if(nirrig_ext.gt.0) then
        if(readvar('obs','irrig_ext').ne.0) return
      endif
      if(nirrig_int.gt.0) then
        if(readvar('obs','irrig_int_next').ne.0) return
      endif
      if(ngw_ext.gt.0) then
        if(readvar('obs','gw_ext').ne.0) return
      endif
      if(nsol.gt.0) then
        if(nsol.ne.1) then
            print*,'Dimension nsol set to ',nsol,' which is '//
     +         'greater than the maximum of one. Fix and restart.'
            return
        endif
        if(readvar('obs','solrad').ne.0) return
      endif
      if(nevap.gt.0) then
        if(readvar('obs','pan_evap').ne.0) return
      endif
      if(nsnow.gt.0) then
        if(readvar('obs','swe').ne.0) return
      endif
      if(nsolute.ne.0) then
        if(sh) then
          if(readvar('obs','spechum').ne.0) sh=.false.
        endif
        if(rh_tf.and.nhum.ne.0) then
          if(readvar('obs','relhum').ne.0) rh_tf=.false.
        endif   
        if(rh_tf.and.sh) then
             print*,'Both relhum and spechum found ',
     +       'in webmod.hydro.dat. Only one allowed. Fix and restart.'
             return
        elseif(rh_tf.or.sh) then  ! all good, use the rh_tf flag to distribute relative humidity,
                                !assign -99 to other variable to ensure fail if used.
           if(.not.warned_rh) then
              if(rh_tf) then
                print*,'Humidity values to be read from variable ',
     +       'relhum in webmod.hydro.dat instead.'
              else
                print*,'Humidity values to be distributed using ',
     +       'spechum in webmod.hydro.dat instead.'
              endif ! rh_tf
              warned_rh=.true.
           endif   ! not warned_rh
           do i=1,nhum
               if(rh_tf) then
                   spechum(i) = -99       
               else
                   relhum(i) = -99
                   if(nhum.gt.1) then
                     print*,'Only one specific humidity field allowed'
                     return                     
                   endif
               endif
           enddo
         else ! both false
             print*,'Either relhum and spechum needed ',
     +       'in webmod.hydro.dat. Neither found'
             return
         endif
      endif

c kludge to get ICRW paper out      
c      relhum(1) = 0.34
c end kludge
c      form_data(1) = 0
c      if(nform.gt.0) then
c        if(readvar('obs','form_data').ne.0) return
c      endif
c      route_on = 0

      if(temp_units.eq.1) then
         if(readvar('obs','tsta_min_c').ne.0) then
            print*,'Set temp_units to 1 (deg C) to '//
     $           'agree with temps in data file'
            return
         end if
         if(readvar('obs','tsta_max_c').ne.0) return
         do 10 i = 1, ntemp
            tsta_min_f(i) = tsta_min_c(i)*1.8 + 32.
            tsta_max_f(i) = tsta_max_c(i)*1.8 + 32.
            tsta_temp_c(i) = (tsta_max_c(i) + tsta_min_c(i))/2.
            tsta_temp_f(i) = (tsta_max_f(i) + tsta_min_f(i))/2.
            if(dt.lt.24.0) then
               if(readvar('obs','tsta_temp_c').ne.0) return
            end if
 10      continue
      else
         if(readvar('obs','tsta_min_f').ne.0) then
            print*,'Set temp_units to 0 (deg F) to '//
     $           'agree with temps in data file'
            return
         end if
         if(readvar('obs','tsta_max_f').ne.0) return
         do 20 i = 1, ntemp
            tsta_min_c(i) = (tsta_min_f(i)-32) / 1.8
            tsta_max_c(i) = (tsta_max_f(i)-32) / 1.8
            tsta_temp_c(i) = (tsta_max_c(i) + tsta_min_c(i))/2.
            tsta_temp_f(i) = (tsta_max_f(i) + tsta_min_f(i))/2.
            if(dt.lt.24.0) then
               if(readvar('obs','tsta_temp_f').ne.0) return
            end if
 20      continue
      endif


c      if(dt.lt.24.) then
c      
c         if(readvar('obs','route_on').ne.0) return
c      
c      endif
c
c     print out nstep, deltim, time and julian dates
c
c      dt = deltim()
c      nstep = getstep()
c      call dattim('now', datetime)
c      julcalen = julian('now', 'calendar')
c      julsolar = julian('now', 'solar')
c      julwater = julian('now', 'water')

c      call dpint4('Nstep        :', nstep, 1, 2)
c      call dpdble('Delta t      :', dt, 1, 2)
c      call dpint4('Date/Time    :', datetime, 6, 2)
c      call dpint4('Jul Calen.   :', julcalen, 1, 2)
c      call dpint4('Jul Solar    :', julsolar, 1, 2)
c      call dpint4('Jul Water    :', julwater, 1, 2)

      obsrun = 0
 
      return
      end

