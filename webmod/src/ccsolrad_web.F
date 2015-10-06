#include "defines.h"
c***********************************************************************
c    ccsolrad_web.f: Subroutine to calculate daily solar radiation from
c                max-min temperature-cloud cover relationship
c
c              1.2 (rmwebb) Estimate cloud cover using basin
c                     tmin and tmax. This eliminates Basin_temp_mru.
c     version: 1.1 (rmwebb) eliminate radiation planes
c                     in favor of nmru (instead of hru)
c                     Also added RCS version control
c     version: 2.0 (rsregan)
c     version: 1.3 (markstro)
c    27 apr  2010 - (rmwebb) 
c                   Port to Fortran 90 with module and dynamic memory
c **********************************************************************
      MODULE WEBMOD_CSOL
      IMPLICIT NONE
#include "fmodules.inc"

C   Dimensions
      integer, save :: nmru, ndays, nmonths, nsol, nrain
C   Declared Variables
      real,save :: orad
      real, save, allocatable :: swrad(:), drad(:), mru_potsw(:)
      real, save, allocatable :: mru_sunhrs(:)
C   Declared Private Variables
      integer, save :: lp(24), nper, lsv, lday, kck
      real, save :: sper, cnt
      real, save, allocatable :: hssv(:), psv(:), ssv(:), hpsv(:)
C   Declared Parameters
      real, save :: radmax, radj_sppt, radj_wppt, crad_coef, crad_exp
      real, save, allocatable :: ccov_slope(:), ccov_intcp(:)
      real, save, allocatable :: ppt_rad_adj(:)
C   Undeclared Static Variables gotten from from other modules - soltab
      real, save, allocatable :: mru_soltab(:,:)    ! from soltab
      real, save, allocatable :: sunhrs_soltab(:,:) ! from soltab
      real, save, allocatable :: mru_cossl(:)       ! from soltab
      real, save, allocatable :: obsrad(:)          ! solrad from obs
      real, save, allocatable :: obs_rain(:)        ! precip from obs
      real,save :: basin_tmin_c, basin_tmax_c, basin_ppt ! from temp and precip
      
      END MODULE WEBMOD_CSOL
c***********************************************************************
c***********************************************************************
c     Main ccsolrad routine
c***********************************************************************

      integer function ccsolrad_web(Arg)

C Arguments
      character(len=*) arg
      CHARACTER(len=256) SVN_ID

C Functions
      integer csoldecl, csolinit, csolrun
      external csoldecl, csolinit, csolrun

      save SVN_ID

      SVN_ID = 
     $     '$Id: ccsolrad_web.f 29 2006-07-06 23:03:45Z rmwebb $ '

c***********************************************************************
      ccsolrad_web = 0

      if ( Arg.eq.'declare' ) then
         ccsolrad_web = csoldecl()

      elseif ( Arg.eq.'initialize' ) then
         ccsolrad_web = csolinit()

      elseif ( Arg.eq.'run' ) then
         ccsolrad_web = csolrun()

      endif

C******Debug level print
      call dpint4('End of ccsolrad, retval = ', ccsolrad_web, 1, 2)

      return
      end
c
c***********************************************************************
c     csoldecl - set up parameters for actual solar radiation
c                computations
c***********************************************************************

      integer function csoldecl()

      USE WEBMOD_CSOL
      IMPLICIT NONE

c***********************************************************************
      csoldecl = 1


!
! Get dimensions
!
      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      ndays = getdim('ndays')
        if ( ndays.eq.-1 ) return
      nmonths = getdim('nmonths')
        if ( nmonths.eq.-1 ) return
      nsol = getdim('nsol')
        if ( nsol.eq.-1 ) return
      nrain = getdim('nrain')
        if ( nrain.eq.-1 ) return

      if ( declpri('csoldecl_lp', 24, 'integer', lp).ne.0 ) return
      if ( declpri('csoldecl_nper', 1, 'integer', nper).ne.0 ) return
      if ( declpri('csoldecl_kck', 1, 'integer', kck).ne.0 ) return
      if ( declpri('csoldecl_lsv', 1, 'integer', lsv).ne.0 ) return
      if ( declpri('csoldecl_lday', 1, 'integer', lday).ne.0 ) return
      if ( declpri('csoldecl_sper', 1, 'real', sper).ne.0 ) return
      if ( declpri('csoldecl_cnt', 1, 'real', cnt).ne.0 ) return

      ALLOCATE (hpsv(nmru))
      if ( declpri('csoldecl_hpsv', nmru, 'real', hpsv)
     +     .ne.0 ) return
      ALLOCATE (hssv(nmru))
      if ( declpri('csoldecl_hssv', nmru, 'real', hssv)
     +     .ne.0 ) return
      ALLOCATE (psv(nmru))
      if ( declpri('csoldecl_psv', nmru, 'real', psv)
     +     .ne.0 ) return
      ALLOCATE (ssv(nmru))
      if ( declpri('csoldecl_ssv', nmru, 'real', ssv)
     +     .ne.0 ) return

      ALLOCATE (swrad(nmru))
      if ( declvar('solrad', 'swrad', 'nmru', nmru, 'real',
     +     'Computed shortwave radiation for each MRU',
     +     'langleys',
     + swrad).ne.0 ) return

      if ( declvar('solrad', 'orad', 'one', 1, 'real',
     +     'Observed or computed solar radiation on a horizontal'//
     +     ' surface',
     +     'langleys',
     + orad).ne.0 ) return

      ALLOCATE (drad(nmru))
      if ( declvar('solrad', 'drad', 'nmru', nmru, 'real',
     +     'Potential shortwave radiation for each MRU',
     +     'langleys',
     + drad).ne.0 ) return

c$$$      if ( decl*var('solrad', 'daily_swrad', 'nmru', nmru, 'real',
c$$$     +     'Computed daily shortwave radiation for each MRU,'//
c$$$     +     ' equals swrad if daily timestep',
c$$$     +     'langleys',
c$$$     + daily_swrad).ne.0 ) return

c$$$      if ( decl*var('solrad', 'daily_potsw', 'nmru', nmru, 'real',
c$$$     +     'Potential shortwave radiation for each radiation'//
c$$$     +     ' plane for each day',
c$$$     +     'langleys',
c$$$     + daily_potsw).ne.0 ) return
c$$$
      ALLOCATE (mru_potsw(nmru))
      if ( declvar('solrad', 'mru_potsw', 'nmru', nmru, 'real',
     +     'Potential shortwave radiation for each MRU '//
     +     'for each day',
     +     'langleys',
     + mru_potsw).ne.0 ) return

      ALLOCATE (mru_sunhrs(nmru))
      if ( declvar('solrad', 'mru_sunhrs', 'nmru', nmru, 'real',
     +     'Hours between sunrise and sunset for each MRU',
     +     'units of 12-hours',
     + mru_sunhrs).ne.0 ) return

      ALLOCATE (ccov_slope(nmonths))
      if ( declparam('solrad', 'ccov_slope', 'nmonths', 'real',
     +     '-.13', '-0.5', '-.01',
     +     'Slope in temperature cloud cover  relationship',
     +     'Coefficient in relationship: cloudcover ='//
     +     ' ccov_intcp + ccov_slope*(tmax_c-tmin_c)',
     +     'ccov per degC').ne.0 ) return

      ALLOCATE (ccov_intcp(nmonths))
      if ( declparam('solrad', 'ccov_intcp', 'nmonths', 'real',
     +     '1.83', '0.0', '5.0',
     +     'Intercept in temperature cloud cover relationship',
     +     'Intercept in relationship: cloudcover ='//
     +     ' ccov_intcp + ccov_slope*(tmax_c-tmin_c)',
     +     'ccov').ne.0 ) return

      if ( declparam('solrad', 'radj_sppt', 'one', 'real',
     +     '0.44', '0.0', '1.0',
     +     'Adjustment to solar radiation on precip day - summer',
     +     'Adjustment factor for computed solar radiation for'//
     +     ' summer day with greater than ppt_rad_adj inch precip',
     +     'none').ne.0 ) return

      if ( declparam('solrad', 'radj_wppt', 'one', 'real',
     +     '0.5', '0.0', '1.0',
     +     'Adjustment to solar radiation on precip day - winter',
     +     'Adjustment factor for computed solar radiation for'//
     +     ' winter day with greater than ppt_rad_adj inch precip',
     +     'none').ne.0 ) return

      if ( declparam('solrad', 'crad_coef', 'one', 'real',
     +     '0.4', '0.1', '0.7',
     +     'Coefficient in cloud cover-solar radiation relationship',
     +     'Coefficient(B) in Thompson(1976) equation:'//
     +     ' Solar radiation = B + (1.-B)*(1-cloudcover)**P'//
     +     ' Varies by region, contour map of values in reference.',
     +     'none').ne.0 ) return

      if ( declparam('solrad', 'crad_exp', 'one', 'real',
     +     '0.61', '0.2', '0.8',
     +     'Exponent in cloud cover-solar radiation relationship',
     +     'Exponent(P) in Thompson(1976) equation:'//
     +     ' Solar radiation = B +(1.-B)*(1-cloudcover)**P'//
     +     ' Author suggests value of 0.61.',
     +     'none')
     +     .ne.0 ) return

      if ( declparam('solrad', 'radmax', 'one', 'real',
     +     '0.8', '0.1', '1.0',
     +     'Maximum percent of potential solar radiation (decimal)',
     +     'The maximum portion of the potential solar radiation'//
     +     ' that may reach the ground due to haze, dust, smog, etc.',
     +     'decimal percent').ne.0 ) return

      ALLOCATE (ppt_rad_adj(nmonths))
      if ( declparam('solrad', 'ppt_rad_adj', 'nmonths', 'real',
     +     '0.02', '0.0', '0.5',
     +     'Radiation reduced if basin precip above this value',
     +     'If basin precip exceeds this value, radiation is'//
     +     ' mutiplied by summer or winter precip adjustment ',
     +     'inches').ne.0 ) return
!
! Allocate copies of variables from other modules
!
      ALLOCATE (mru_soltab(ndays,nmru))
      ALLOCATE (sunhrs_soltab(ndays,nmru))
      ALLOCATE (mru_cossl(nmru))
      if ( Nsol.gt.0 ) then
        ALLOCATE (obsrad(nsol))
      else
        ALLOCATE (obsrad(1))
      endif
      ALLOCATE (obs_rain(nrain))

      csoldecl = 0

      return 
      end

c***********************************************************************
c     csolinit - Initialize ccsolrad module - get parameter values,
c***********************************************************************

      integer function csolinit()


      USE WEBMOD_CSOL
      IMPLICIT NONE

C Local Variables
      integer i, nstep

C Save Variables
      integer lp_init(24)
      save lp_init

      data lp_init/19, 13, 15, 13, 15, 14, 14, 15, 14, 15, 14, 21, 20,
     +             15, 14, 15, 15, 14, 15, 14, 14, 14, 14, 18/

c***********************************************************************
      csolinit = 1

      nmru = getdim('nmru')
      if ( nmru.eq.-1 ) return

c$$$      Nradpl = getdim('nradpl')
c$$$      if ( Nradpl.eq.-1 ) return

      Nsol = getdim('nsol')
      if ( Nsol.eq.-1 ) return

      if ( getparam('solrad', 'ccov_slope', nmonths, 'real',
     +     ccov_slope).ne.0 ) return

      if ( getparam('solrad', 'ccov_intcp', nmonths, 'real',
     +     ccov_intcp).ne.0 ) return

      if ( getparam('solrad', 'radj_sppt', 1, 'real', radj_sppt)
     +     .ne.0 ) return

      if ( getparam('solrad', 'radj_wppt', 1, 'real', radj_wppt)
     +     .ne.0 ) return

      if ( getparam('solrad', 'crad_coef', 1, 'real', crad_coef)
     +     .ne.0 ) return

      if ( getparam('solrad', 'crad_exp', 1, 'real', crad_exp)
     +     .ne.0 ) return

      if ( getparam('solrad', 'radmax', 1, 'real', radmax)
     +     .ne.0 ) return

      if ( getparam('solrad', 'ppt_rad_adj', nmonths, 'real',
     +      ppt_rad_adj).ne.0 ) return

C**get soltab variables once at beginning of simulation, not done in run
      if ( getvar('soltab', 'mru_soltab', ndays*nmru, 'real', 
     + mru_soltab).ne.0 ) return

      if ( getvar('soltab', 'sunhrs_soltab', ndays*nmru, 'real', 
     + sunhrs_soltab).ne.0 ) return

      if ( getvar('soltab', 'mru_cossl', nmru, 'real', 
     + mru_cossl).ne.0 ) return

      nstep = getstep()
      if ( nstep.eq.0 ) then
 
         do 20 i = 1, 24
            Lp(i) = lp_init(i)
   20    continue

         do 22 i = 1, Nmru
            Drad(i) = 0.0
            Mru_potsw(i) = 0.0
            Mru_sunhrs(i) = 0.0
   22    continue

         do 30 i = 1, Nmru
            Swrad(i) = 0.0
   30    continue

         Orad = 0.0
         Lday = 0
         Nper = 1
         Kck = 0
         Lsv = 1
         Cnt = 0.0
         Sper = 0.0

         do 40 i = 1, Nmru
            Hpsv(i) = 0.0
            Hssv(i) = 0.0
            Psv(i) = 0.0
            Ssv(i) = 0.0
   40    continue

      endif

      csolinit = 0

      return
      end

c***********************************************************************
c     csolrun - Computes actual solar radiation on horizontal surface,
c               then determines values for each MRU.
c***********************************************************************

      integer function csolrun()

      USE WEBMOD_CSOL
      USE WEBMOD_IO, ONLY: nowtime
      IMPLICIT NONE

C Local Variables
      logical ppt
      integer mo, j, jkeep
      integer ip, jp, ipday, jpday, day
      integer myr, jsol, k
      integer nstep
c$$$      real horad, plrad(nmru), ccov, pptadj, radadj
      real horad, ccov, pptadj, radadj
      real hdy, c, r, ckrad

C Save Variables
      integer sday(13)
      save sday
      real nearzero
      parameter ( nearzero = 1.E-30 )
      data sday/356, 10, 23, 38, 51, 66, 80, 94, 109, 123, 138, 152,173/

c***********************************************************************
      csolrun = 1

c      call dattim('now', nowtime)
      myr = nowtime(1)
      mo = nowtime(2)
      day = nowtime(3)

      if ( Lday.eq.day ) goto 400

      Lday = day

      jsol = julian('now', 'solar')

      obsrad(1) = 0.0
      if ( Nsol.gt.0 ) then
         if ( getvar('obs', 'solrad', nsol, 'real', obsrad)
     +        .ne.0 ) return
      endif
      ckrad = obsrad(1)

c  Maintain degrees C so that slope remains constant

      if ( getvar('temp', 'basin_tmin_c', 1, 'real', basin_tmin_c)
     +     .ne.0 ) return

      if ( getvar('temp', 'basin_tmax_c', 1, 'real', basin_tmax_c)
     +     .ne.0 ) return

      if ( getvar('obs', 'precip', nrain, 'real', obs_rain)
     +     .ne.0 ) return

      if ( getvar('precip', 'basin_ppt', 1, 'real', basin_ppt)
     +     .ne.0 ) return

c **** begin new code from solrad.f

      nstep = getstep()
      if ( nstep.eq.1 ) then
c***
c*** initialize variables
c***
         k = myr
         if ( jsol.lt.11 ) k = myr + 1
         Lp(5) = 15
         if ( isleap(k).ne.0 ) Lp(5) = 16
         Nper = 1

         do 20 j = 1, 24
            jkeep = j
            Nper = Nper + Lp(j)
            if ( jsol.lt.Nper ) goto 30
   20    continue

   30    Lsv = jkeep
         Sper = Lp(Lsv)
         c = Nper - jsol
         Cnt = Sper - c
         Kck = 1
         if ( Cnt.lt.nearzero ) Kck = 0
      else
c***
c*** compute location of jsol in solar table rad
c***
         if ( Nper.ne.jsol ) goto 140
         if ( Lsv.eq.24 ) then
            Nper = 1
            Lsv = 0
            Lp(5) = 15
            k = myr + 1
            if ( isleap(k).ne.0 ) Lp(5) = 16
            goto 140
         else

            Lsv = Lsv + 1
            Sper = Lp(Lsv)
            Nper = Nper + Lp(Lsv)
         endif
      endif

      k = 13 - Lsv

      if ( k.lt.0 ) then
         ip = 13 + k
         jp = ip - 1
      elseif ( k.eq.0 ) then
         ip = 13
         jp = 12
      else
         ip = 13 - k
         jp = ip + 1
      endif

      ipday = sday(ip)
      jpday = sday(jp)

      do 110 j = 1, Nmru
         Hpsv(j) = Sunhrs_soltab(ipday, j)
         Hssv(j) = Sunhrs_soltab(jpday, j) - Hpsv(j)
         Psv(j) = Mru_soltab(ipday, j)
         Ssv(j) = Mru_soltab(jpday, j) - Psv(j)
  110 continue

      if ( Kck.eq.0 ) then

         do 120 j = 1, Nmru
            hdy = Hpsv(j)
            Drad(j) = Psv(j)
            Mru_potsw(j) = Drad(j)
            Mru_sunhrs(j) = hdy
  120    continue

         Cnt = 1.
         goto 160
      else
         Kck = 0
      endif
c***
c*** linearly interpolate between table values
c***
  140 r = Cnt/Sper

      do 150 j = 1, Nmru
         hdy = Hpsv(j) + (r*Hssv(j))
         Drad(j) = Psv(j) + (r*Ssv(j))
         Mru_potsw(j) = Drad(j)
         Mru_sunhrs(j) = hdy
  150 continue

      Cnt = Cnt + 1.

c **** end new code from solrad.f

  160 ppt = .FALSE.

      if ( basin_ppt.gt.Ppt_rad_adj(mo) ) ppt = .TRUE.

      horad = Drad(1)
      pptadj = 1.
      if ( ckrad.lt.nearzero .or. ckrad.gt.10000. ) then
         ccov = (Ccov_slope(mo)*(basin_tmax_c-basin_tmin_c))
     +          + Ccov_intcp(mo)
         if ( ccov.lt.nearzero ) ccov = 0.
         if ( ccov.gt.1. ) ccov = 1.
         if ( ppt ) then
            pptadj = Radj_wppt
            if ( mo.gt.4 .and. mo.lt.10 ) pptadj = Radj_sppt
         endif
         radadj = Crad_coef + (1.-Crad_coef)*((1.-ccov)**Crad_exp)
         if ( radadj.gt.Radmax ) radadj = Radmax
         radadj = radadj*pptadj
         obsrad(1) = radadj*horad
      endif
      Orad = obsrad(1)

      do 200 j = 1, nmru
         swrad(j) = (Drad(j)/horad)*obsrad(1)/Mru_cossl(j)
  200 continue

  400 csolrun = 0

      return
      end

