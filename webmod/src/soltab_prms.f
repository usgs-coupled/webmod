#include "defines.h"
c***********************************************************************
c     soltab_prms.f: compute potential solar radiation for timesteps
c               from 15 minutes to 1 day for each model response unit.
c               Computes day length, sunrise and sunset times.
c  
c***********************************************************************
c
c   Modify to compute for each model response unit instead of limiting
c   the computation to 'radiation plane' categories. Rick Webb 11Sept02.
c
c   Added version control 04Sept03 - RMTW
c    27 apr  2010 - Port to Fortran 90 with module and dynamic memory
c **********************************************************************
      MODULE WEBMOD_SOLTAB
      IMPLICIT NONE
      include 'fmodules.inc'

      integer, save :: nmru, ndays

      real, save :: l0, l1, l2, r1, dim(13)
      real, save, allocatable :: mru_aspect(:)
      real, save, allocatable :: mru_lat(:)
      real, save, allocatable :: mru_slope(:)

      real, save, allocatable :: mru_soltab(:,:)
      real, save, allocatable :: sunhrs_soltab(:,:)
      real, save, allocatable :: mru_cossl(:)
      real, save, allocatable :: mru_sunrise(:)
      real, save, allocatable :: mru_sunset(:)
      real, save, allocatable :: mru_sunrise2(:)
      real, save, allocatable :: mru_sunset2(:)

      END MODULE WEBMOD_SOLTAB

c***********************************************************************
c
c     Main soltab routine
c

      integer function soltab_prms(arg)

      USE WEBMOD_SOLTAB
      IMPLICIT NONE
      
      character(len=*) arg
      CHARACTER(len=256) SVN_ID

      integer stdecl, stinit
      save SVN_ID

      SVN_ID = 
     $     '$Id: soltab_prms.f 29 2006-07-06 23:03:45Z rmwebb $ '

      soltab_prms = 0

      if(arg.eq.'declare') then
        soltab_prms = stdecl()

      else if(arg.eq.'initialize') then
        soltab_prms = stinit()

      end if

C******Debug level print
      call dpint4('End of soltab, retval = ', soltab_prms, 1, 2)
      return
      end

c***********************************************************************
c 
c     stdecl - set up parameters for solar radiation computations
c

      integer function stdecl()

      USE WEBMOD_SOLTAB
      IMPLICIT NONE

      stdecl = 1
!
! Get dimensions
!
      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      ndays = getdim('ndays')
        if ( ndays.eq.-1 ) return

      ALLOCATE(mru_soltab(ndays,nmru))
      if(declvar('soltab', 'mru_soltab', 'ndays,nmru',
     +      ndays*nmru,'real',
     +     'Potential shortwave radiation for each model response '//
     +     'unit for each timestep',
     +     'langleys',
     +   mru_soltab).ne.0) return

      ALLOCATE(sunhrs_soltab(ndays,nmru))
      if(declvar('soltab', 'sunhrs_soltab', 'ndays,nmru',
     +      ndays*nmru,'real',
     +     'Hours between sunrise and sunset for model response unit',
     +     'hours',
     +   sunhrs_soltab).ne.0) return

      ALLOCATE(mru_cossl(nmru))
      if(declvar('soltab', 'mru_cossl', 'nmru', nmru,
     +     'real',
     +     'Cosine of the slope of the model response unit',
     +     'none',
     +   mru_cossl).ne.0) return

      ALLOCATE(mru_sunrise(nmru))
      if(declvar('soltab', 'mru_sunrise', 'nmru', nmru,
     +     'real',
     +     'Time of sunrise on the model response unit',
     +     'hours',
     +   mru_sunrise).ne.0) return

      ALLOCATE(mru_sunset(nmru))
      if(declvar('soltab', 'mru_sunset', 'nmru', nmru,
     +     'real',
     +     'Time of sunset on model response unit',
     +     'hours',
     +   mru_sunset).ne.0) return

      ALLOCATE(mru_sunrise2(nmru))
      if(declvar('soltab', 'mru_sunrise2', 'nmru', nmru,
     +     'real',
     +     'Time of second sunrise on the model response unit if '//
     +     'applicable',
     +     'hours',
     +   mru_sunrise2).ne.0) return

      ALLOCATE(mru_sunset2(nmru))
      if(declvar('soltab', 'mru_sunset2', 'nmru', nmru,
     +     'real', 'Time of second sunset on the model response '//
     +     'unit if applicable', 'hours',
     +   mru_sunset2).ne.0) return

      ALLOCATE(mru_slope(nmru))
      if(declparam('soltab', 'mru_slope', 'nmru', 'real',
     +   '0.0', '0.0', '10.0',
     +   'MRU slope in decimal vertical feet/horizontal feet',
     +   'MRU slope in decimal vertical feet/horizontal feet',
     +   'decimal percent')
     +   .ne.0) return

      ALLOCATE(mru_aspect(nmru))
      if(declparam('soltab', 'mru_aspect', 'nmru', 'real',
     +   '0.0', '0.0', '360.0',
     +   'Aspect of model response unit',
     +   'Aspect for each model response unit',
     +   'degrees')
     +   .ne.0) return

      ALLOCATE(mru_lat(nmru))
      if(declparam('soltab', 'mru_lat', 'nmru', 'real',
     +   '40.0', '20.0', '60.0',
     +   'Latitude of model response unit',
     +   'Latitude of model response unit',
     +   'degrees')
     +   .ne.0) return

      stdecl = 0

      return
      end

c***********************************************************************
c
c     stinit - Initialize soltab module - get parameter values,
c                compute mru_soltab (potential shortwave radiation)
c                and sunhrs_soltab (hours between sunrise and sunset)
c                for each model response unit for each day of the year.
c                
c

      integer function stinit()

      USE WEBMOD_SOLTAB
      IMPLICIT NONE

      integer n, jday(13), jd, is

      real v, w, x, y, func3
      real a
      real r0
      real d1, t, tx, dm(13)
      real t0, t1, t2, t3, t6, t7, t8, t9
      real e(13), d, i, day 

      data e/2.06699,2.06317,2.05582,2.04520,2.03243,2.01706,2.00080,
     +1.98553,1.96990,1.95714,1.94689,1.94005,1.93616/

      data dm/-.410152,-.383391,-.337430,-.27198,-.190532,-.09832,0.,
     +.09832,.190532,.27198,.33743,.383391,.410152/

      data jday/356,10,23,38,51,66,80,94,109,123,138,152,173/

c*****statement function
      func3(v,w,x,y)=r1*(sin(d)*sin(w)*(x-y)*3.8197+cos(d)*cos(w)*
     +(sin(x+v)-sin(y+v))*12./3.14159)

      stinit = 1

      if(nmru.eq.-1) return

      if(getparam('soltab', 'mru_slope', nmru, 'real',
     +   mru_slope).ne.0) return

      if(getparam('soltab', 'mru_aspect', nmru, 'real',
     +   mru_aspect).ne.0) return

      if(getparam('soltab', 'mru_lat', nmru, 'real',
     +   mru_lat).ne.0) return

      r0 = 2.0 


      do 210 n=1,nmru

       i = mru_slope(n)
       a = mru_aspect(n)
       l0 = mru_lat(n)
       i=atan(i)
       a=a/57.2958
       l0=l0/57.2958
       mru_cossl(n) = cos(i)
       r0=2.0

       l1=asin(cos(i)*sin(l0)+sin(i)*cos(l0)*cos(a))
       d1=cos(i)*cos(l0)-sin(i)*sin(l0)*cos(a)
       if(d1.eq.0.) d1=.0000000001
       l2=atan(sin(i)*sin(a)/d1)
       if(d1.lt.0.) l2=l2+3.14159


      do 200 is=1,13
       jd = jday(is)
       day=jday(is)
       d=dm(is)
       r1=60.*e(is)
       t=0.
       tx=-tan(l1)*tan(d)
       if(tx.lt.-1.0) t=3.14159
       if(tx.gt.1.0)  t=0.
       if(abs(tx).le.1.) t=acos(tx)
       t7=t-l2
       t6=-t-l2
       t=0.
       tx=-tan(l0)*tan(d)
       if(tx.lt.-1.0) t=3.14159
       if(tx.gt.1.0) t=0.
       if(abs(tx).le.1.) t=acos(tx)
       t1=t
       t0=-t
       t3=t7
       if(t7.gt.t1) t3=t1
       t2=t6
       if(t6.lt.t0) t2=t0
       if(i.ne.0) go to 150
       mru_soltab(jd,n)=func3 (0.,l0,t1,t0)
       sunhrs_soltab(jd,n)=(t1-t0)*3.8197/12.

       go to 200
  150 if(t3.ge.t2) go to 160
      t2=0.
      t3=0.
  160 t6=t6+6.28318
      if(t6.lt.t1) go to 190
      t7=t7-6.28318
      if(t7.gt.t0) go to 180
      mru_soltab(jd,n)=func3 (l2,l1,t3,t2)
      sunhrs_soltab(jd,n)=(t3-t2)*3.8197/12.

      go to 200
  180 t8=t0
      t9=t7
      go to 195
  190 t8=t6
      t9=t1
  195 mru_soltab(jd,n)=func3 (l2,l1,t3,t2) + func3 (l2,l1,t9,t8)
      sunhrs_soltab(jd,n)=(t3-t2+t9-t8)*3.8197/12.


 200  continue 
 210  continue

c      do 300 j=1,nradpl
c      write(*,8000) j, sunhrs_soltab(356,j),sunhrs_soltab(10,j),
c     + sunhrs_soltab(23,j), sunhrs_soltab(38,j),
c     + sunhrs_soltab(51,j), sunhrs_soltab(66,j),
c     + sunhrs_soltab(80,j), sunhrs_soltab(94,j),
c     + sunhrs_soltab(109,j), sunhrs_soltab(123,j),
c     + sunhrs_soltab(138,j), sunhrs_soltab(152,j),
c     + sunhrs_soltab(173,j)

c 8000 format(1x,i5,13f7.4)
c 300  continue

      stinit = 0

      return
      end


