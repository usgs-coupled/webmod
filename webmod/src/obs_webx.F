#include "defines.h"
c***********************************************************************
c     obs_webx.f: modified from obs_prms -> obs_webmod to
c                 include XYZ data input. -RMTW
c         4 Mar 05 - changed tmin to tsta_min_c and _f.
c                    Same with tmax and tstemp.
c **********************************************************************
c                        added temp_units parameter
c                        Snow water equivalence
c                        and irrigation schedules,
c                        and ground water influx
c
c                        added RCS version control,
c
c **********************************************************************
c
c     obsdecl - makes public variable declarations for the
c                     obs module
c
 
      integer function obsxdecl(runoff, precip,
     $     irrig_ext, irrig_int_next, gw_ext,
     $     tsta_min_c, tsta_max_c, tsta_temp_c,
     $     tsta_min_f, tsta_max_f, tsta_temp_f,
     $     solrad, swe, pan_evap, route_on, form_data,
     $     rain_day)
      use WEBMOD_IRRIG, only : nform
      implicit none
#include "fmodules.inc"
      real runoff(MAXOBS)
      real precip(MAXRAIN)
      real irrig_ext(MAXIRRIG)
      real irrig_int_next(MAXIRRIG)
      real gw_ext(maxirrig)
      real tsta_min_c(MAXTEMP), tsta_max_c(MAXTEMP)
      real tsta_min_f(MAXTEMP), tsta_max_f(MAXTEMP)
      real tsta_temp_c(MAXTEMP), tsta_temp_f(MAXTEMP)
      real solrad(MAXSOL)
      real swe(MAXSNOPIL)
      real pan_evap(MAXEVAP)
      integer form_data(MAXFORM), route_on, rain_day

 

      obsxdecl = 1

      if(declparam('obs', 'rain_code', 'nmonths', 'integer',
     +   '2', '1', '5',
     +   'Code indicating rule for precip station use',
     +   'Code indicating rule for precip station use: '//
     +   '1 = only precip if the regression stations have precip,'//
     +   '2 = only precip if any station in the basin has precip,'//
     +   '3 = precip if xyz says so,'//
     +   '4 = only precip if rain_day variable is set to 1,'//
     +   '5 = only precip if psta_freq_nuse stations see precip',
     +   'none')
     +   .ne.0) return

      if (declvar ('obs', 'rain_day', 'one', 1, 'integer',
     +   'Flag to force rain day',
     +   'none', rain_day).ne.0) return

      if(declvar('obs', 'runoff', 'nobs', MAXOBS, 'real',
     +     'Observed runoff for each gage',
     +     'cfs',
     +   runoff).ne.0) return

      if(declvar('obs', 'precip', 'nrain', MAXRAIN, 'real',
     +     'Observed precip at each rain gage',
     +     'inches',
     +   precip).ne.0) return

      if(declvar('obs', 'gw_ext', 'ngw_ext', MAXIRRIG, 'real',
     +  'Influx of groundwater from sources external to basin',
     $  'cfs/mile', gw_ext).ne.0) return

      if(declvar('obs', 'irrig_ext', 'nirrig_ext', MAXIRRIG, 'real',
     +  'Depth of applied irrigation from an external source',
     $  'inches', irrig_ext).ne.0) return

      if(declvar('obs', 'irrig_int_next', 'nirrig_int', MAXIRRIG,
     $     'real', 'Depth of irrigation from an internal source '//
     $     'to apply the next day',
     $     'inches', irrig_int_next).ne.0) return

      if(declparam('obs', 'temp_units', 'one', 'integer',
     +   '1', '0', '1',
     +   'Units for observed temperature',
     +   'Units for observed temperature, 0=F, 1=C',
     +   'none').ne.0) return

      if(declvar('obs', 'tsta_min_c', 'ntemp', MAXTEMP, 'real',
     +     'Observed daily minimum temperature',
     +     'degrees celsius',
     +   tsta_min_c).ne.0) return

      if(declvar('obs', 'tsta_max_c', 'ntemp', MAXTEMP, 'real',
     +     'Observed daily maximum temperature',
     +     'degrees celsius',
     +   tsta_max_c).ne.0) return

      if(declvar('obs', 'tsta_temp_c', 'ntemp', MAXTEMP, 'real',
     +     'Average observed daily temperature',
     +     'degrees celsius',
     +   tsta_temp_c).ne.0) return

      if(declvar('obs', 'tsta_min_f', 'ntemp', MAXTEMP, 'real',
     +     'Observed daily minimum temperature',
     +     'degrees fahrenheit',
     +   tsta_min_f).ne.0) return

      if(declvar('obs', 'tsta_max_f', 'ntemp', MAXTEMP, 'real',
     +     'Observed daily maximum temperature',
     +     'degrees fahrenheit',
     +   tsta_max_f).ne.0) return

      if(declvar('obs', 'tsta_temp_f', 'ntemp', MAXTEMP, 'real',
     +     'Average observed daily temperature',
     +     'degrees fahrenheit',
     +   tsta_temp_f).ne.0) return

      if(declvar('obs', 'solrad', 'nsol', MAXSOL, 'real',
     +     'Observed solar radiation',
     +     'langleys',
     +   solrad).ne.0) return

      if(declvar('obs', 'pan_evap', 'nevap', MAXEVAP, 'real',
     +     'Observed pan evaporation',
     +     'inches',
     +   pan_evap).ne.0) return

      if(declvar('obs', 'swe', 'nsnow', MAXSNOPIL, 'real',
     +     'Observed snow water equivalent of snowpack',
     +     'inches',
     +   swe).ne.0) return

      if(nform.eq.1) then
        if(declvar('obs', 'form_data', 'nform', MAXFORM, 'integer',
     +     'Form of precip - 0=not known, 1=snow, 2=rain',
     +     'none',
     +   form_data).ne.0) return
      endif
      if(declvar('obs', 'route_on', 'one', 1, 'integer',
     +     'Kinematic routing switch - 1=storm period',
     +     'none',
     +   route_on).ne.0) return
 

      obsxdecl = 0

      return
      end

c **********************************************************************
c
c     obsinit - initializes obs module
c

      integer function obsxinit(temp_units, ntemp, rain_code)
      IMPLICIT NONE

#include "fmodules.inc"

      !integer datetime(6), temp_units, ntemp, rain_code(maxmo)
      integer temp_units, ntemp, rain_code(maxmo)

      obsxinit = 1

      ntemp =  getdim('ntemp')

      if(getparam('obs', 'rain_code', MAXMO, 'integer',
     +   rain_code)
     +   .ne.0) return

      if(getparam('obs', 'temp_units', 1,
     $     'integer',temp_units) .ne.0) return
c
c     print out start and end times
c
c      call dattim('start', datetime)
c      call dpint4('Start        :', datetime, 6, 2)
c      call dattim('end', datetime)
c      call dpint4('End          :', datetime, 6, 2)
c      call dpstr(' ', 2)

      obsxinit = 0
 
      return
      end

c **********************************************************************
c
c     obsrun - runs obs module
c
 
      integer function obsxrun(route_on, form_data, tsta_min_c,
     $     tsta_max_c, tsta_temp_c, tsta_min_f, tsta_max_f,
     $     tsta_temp_f, temp_units, ntemp, rain_day, rain_code)
      IMPLICIT NONE

#include "fmodules.inc"

      integer i, rain_code(maxmo), rain_day
      !integer nstep, datetime(6), route_on, temp_units, ntemp
      integer datetime(6), route_on, temp_units, ntemp
c$$$      integer julcalen, julsolar, julwater
      integer form_data(MAXFORM)
      real tsta_min_c(maxtemp), tsta_max_c(maxtemp)
      real tsta_temp_c(maxtemp)
      real tsta_min_f(maxtemp), tsta_max_f(maxtemp)
      real tsta_temp_f(maxtemp)
      double precision dt

      obsxrun = 1

      dt = deltim()
      call dattim ('now', datetime)

      if(readvar('obs','runoff').ne.0) return
      if(readvar('obs','precip').ne.0) return
      if(readvar('obs','irrig_ext').ne.0) return
      if(readvar('obs','irrig_int_next').ne.0) return
      if(readvar('obs','gw_ext').ne.0) return
      if(readvar('obs','solrad').ne.0) return
      if(readvar('obs','pan_evap').ne.0) return
      if(readvar('obs','swe').ne.0) return

      form_data(1) = 0
      if(readvar('obs','form_data').ne.0) return

      route_on = 0

      if(temp_units.eq.1) then
         if(readvar('obs','tsta_min_c').ne.0) return
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
         if(readvar('obs','tsta_min_f').ne.0) return
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

      rain_day = 1
      if (rain_code(datetime(2)).eq.4) then
         if(readvar('obs','rain_day').ne.0) return
      end if

      if(dt.lt.24.) then
         if(readvar('obs','route_on').ne.0) return
      endif
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

      obsxrun = 0
 
      return
      end

c **********************************************************************
c
c     main obs routine
c

      integer function obs_webx(arg)
      IMPLICIT NONE

#include "fmodules.inc"

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      integer retval, obsxdecl, obsxinit, obsxrun
      integer route_on, form_data(MAXFORM)
      integer temp_units,ntemp
      integer rain_code(maxmo), rain_day

      real runoff(MAXOBS)
      real precip(MAXRAIN)
      real irrig_ext(MAXIRRIG)
      real irrig_int_next(MAXIRRIG)
      real gw_ext(maxirrig) 

      real tsta_min_c(MAXTEMP), tsta_max_c(MAXTEMP)
      real tsta_min_f(MAXTEMP), tsta_max_f(MAXTEMP)
      real tsta_temp_c(MAXTEMP), tsta_temp_f(MAXTEMP)
      real solrad(MAXSOL)
      real swe(MAXSNOPIL)
      real pan_evap(MAXEVAP)

      save runoff, precip
      save irrig_ext, irrig_int_next, gw_ext
      save tsta_min_c, tsta_max_c, tsta_temp_c
      save tsta_min_f, tsta_max_f, tsta_temp_f
      save solrad, swe, pan_evap
      save temp_units, ntemp
      save route_on, form_data, rain_day, rain_code
      save SVN_ID

      SVN_ID = 
     $     '$Id: obs_webx.f 29 2006-07-06 23:03:45Z rmwebb $ '

      retval = 0

      if (arg.eq.'declare') then
         retval = obsxdecl(runoff, precip,
     $     irrig_ext, irrig_int_next, gw_ext,
     $     tsta_min_c, tsta_max_c, tsta_temp_c,
     $     tsta_min_f, tsta_max_f, tsta_temp_f,
     $     solrad, swe, pan_evap, route_on, form_data,
     $     rain_day)

      else if (arg.eq.'initialize') then
         retval = obsxinit(temp_units, ntemp, rain_code)

      else if (arg.eq.'run') then
         retval = obsxrun(route_on, form_data, tsta_min_c,
     $     tsta_max_c, tsta_temp_c, tsta_min_f, tsta_max_f,
     $     tsta_temp_f, temp_units, ntemp, rain_day, rain_code)

      end if

      obs_webx = retval

      return
      end
