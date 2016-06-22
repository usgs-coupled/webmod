#include "defines.h"
c***********************************************************************
c     potet_hamon.f:  Determine whether transpiration is occurring and
c                     compute the potential evapotranspiration for each
c                     MRU
c     version: 1.3 (markstro) modifed for xtop_prms by Rick Webb
c                  note areas in km2 and use of MRU instead of HRU
c
c                  5 Sept 03 - Added version control - RMTW
c        06may10 - Port to Fortran 90 with module and dynamic memory
c
c***********************************************************************
c 
      MODULE WEBMOD_POTET
      IMPLICIT NONE
#include "fmodules.inc"

C   Dimensions
      integer, save :: nmru, nmonths

C   Declared Variables
      integer, save, allocatable:: transp_on(:)
      real, save :: basin_potet
      real, save, allocatable :: potet(:)
C   Declared Parameters
      real, save :: basin_area
      integer, save, allocatable:: transp_beg(:), transp_end(:)
      real, save, allocatable :: mru_area(:), hamon_coef(:) 
      real, save, allocatable :: transp_tmax_c(:)

C   Private Variables
      integer, save :: lday
      integer, save,  allocatable:: transp_check(:)
      real, save, allocatable :: daily_potet(:), tmax_sum(:)

C   Undeclared Static Variables gotten from from other modules
      real, save, allocatable :: mru_sunhrs(:), tavgc(:), tmax_c(:)
    
      END MODULE WEBMOD_POTET
c
c***********************************************************************
c***********************************************************************
c
c     Main potet_hamon routine
c

      integer function potet_hamon_prms(arg)
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      integer petdecl, petinit, petrun
 
      save SVN_ID

      SVN_ID = 
     $     '$Id: potet_hamon_prms.f 29 2006-07-06 23:03:45Z rmwebb $ '
    
      potet_hamon_prms = 0

      if(arg.eq.'declare') then
        potet_hamon_prms = petdecl()

      else if(arg.eq.'initialize') then
        potet_hamon_prms = petinit()

      else if(arg.eq.'run') then
        potet_hamon_prms = petrun()

      end if

      return
      end


c***********************************************************************
c 
c     petdecl - set up parameters for potential et computations
c

      integer function petdecl()

     
      USE WEBMOD_POTET
      IMPLICIT NONE

      petdecl = 1

! Get dimensions
      nmru = getdim('nmru')
      if(nmru.eq.-1) return
      nmonths = getdim('nmonths')
      if(nmonths.eq.-1) return


      if(declpri('potet_lday', 1, 'integer', lday).ne.0) return

      ALLOCATE (daily_potet(Nmru))
      if(declpri('potet_daily_potet', nmru, 'real',
     +           daily_potet).ne.0) return
 
      ALLOCATE (transp_on(Nmru))
      if(declvar('potet', 'transp_on', 'nmru', nmru, 'integer',
     +    'Switch indicating whether transpiration is occurring,'//
     +    ' 0=no 1=yes',
     +    'none',
     +   transp_on).ne.0) return
 
      ALLOCATE (potet(Nmru))
      if(declvar('potet', 'potet', 'nmru', nmru, 'real',
     +    'Potential evapotranspiration on an MRU',
     +    'inches',
     +   potet).ne.0) return
 
      if(declvar('potet', 'basin_potet', 'one', 1, 'real',
     +    'Basin area-weighted average of potential et',
     +    'inches',
     +   basin_potet).ne.0) return
 
c      if(declparam('potet_hamon', 'mru_radpl', 'nmru', 'integer',
c     +   '1', 'bounded', 'nradpl',
c     +   'Index of radiation plane for MRU',
c     +   'Index of radiation plane used to compute solar '//
c     +   'radiation for an MRU',
c     +   'none').ne.0) return
 
      ALLOCATE (mru_area(Nmru))
      if(declparam('potet_hamon', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('potet_hamon', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return
      
      ALLOCATE (transp_beg(Nmru))
      if(declparam('potet', 'transp_beg', 'nmru', 'integer',
     +   '4', '1', '12',
     +   'Month to begin testing for transpiration',
     +   'Month to begin summing tmax_c for each MRU; when sum is '//
     +   '>= to transp_tmax_c, transpiration begins',
     +   'month').ne.0) return

      ALLOCATE (transp_end(Nmru))
      if(declparam('potet', 'transp_end', 'nmru', 'integer',
     +   '10', '1', '12',
     +   'End month of transpiration period',
     +   'Last month for transpiration computations; '//
     +   'Transpiration is computed through end of month',
     +   'month').ne.0) return
 
      ALLOCATE (transp_tmax_c(Nmru))
      if(declparam('potet', 'transp_tmax_c', 'nmru', 'real',
     +   '500.', '0.', '1000.',
     +   'Tmax index to determine start of transpiration',
     +   'Temperature index to determine the specific date of the '//
     +   'start of the transpiration period.  Subroutine sums tmax_c '//
     +   'for days above freezing for each MRU starting with the '//
     $   'first day of month transp_beg.  When the sum exceeds '//
     $   'this index, transpiration begins',
     +    'degrees celsius').ne.0) return

      ALLOCATE (hamon_coef(Nmonths))
      if(declparam('potet', 'hamon_coef', 'nmonths', 'real',
     +   '.0055', '.004', '.015',
     +   'Monthly air temp coefficient - Hamon',
     +   'Monthly air temperature coefficient used in Hamon '//
     +   'potential evapotranspiration computations, see '//
     +   'PRMS manual ',
     +   '????').ne.0) return

      ALLOCATE(transp_check(nmru))
      ALLOCATE(tmax_sum(nmru))
! getvars      
      ALLOCATE(mru_sunhrs(nmru))
      ALLOCATE(tavgc(nmru))
      ALLOCATE(tmax_c(nmru))

      petdecl = 0

      return
      end

c***********************************************************************
c
c     petinit - Initialize potet module - get parameter values,
c                set initial transp_on switch
c

      integer function petinit()


      USE WEBMOD_POTET
      IMPLICIT NONE

      integer starttime(6), mo, day, i

      petinit = 1

      if(getparam('potet', 'transp_beg', nmru, 'integer', transp_beg)
     +   .ne.0) return

      if(getparam('potet', 'transp_end', nmru, 'integer', transp_end)
     +   .ne.0) return

      if(getparam('potet', 'transp_tmax_c', nmru,
     $     'real', transp_tmax_c) .ne.0) return

      if(getparam('potet', 'hamon_coef', nmonths, 'real', hamon_coef)
     +   .ne.0) return

c      if(getparam('soltab', 'mru_radpl', nmru, 'integer', mru_radpl)
c     +   .ne.0) return

      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return


      
      call dattim('start',starttime)
      mo = starttime(2)
      day = starttime(3)
      
      do 10 i = 1,nmru
        transp_on(i) = 0
        tmax_sum(i) = 0.
        transp_check(i) = 0 

      if(mo.eq.transp_beg(i)) then
          if(day.gt.10) then
           transp_on(i) = 1
          else
            transp_check(i) = 1
          end if

        else if((transp_end(i)-transp_beg(i)).gt.0) then
          if(mo.gt.transp_beg(i).and.mo.lt.transp_end(i))
     +       transp_on(i) = 1

        else
          if((mo.gt.transp_beg(i).and.mo.le.12).or.
     +       (mo.ge.1.and.mo.lt.transp_end(i))) transp_on(i) = 1

        end if
   10 continue

      lday = 0

      petinit = 0

      return
      end

c***********************************************************************
c
c      petrun - Keeps track of transpiration on or off and computes
c               potential et for each MRU each day
c

      integer function petrun()

      USE WEBMOD_POTET
      USE WEBMOD_IO, ONLY: nowtime
      IMPLICIT NONE

      integer mo, i, jday, day
      real dyl, vpsat, vdsat, factor
      double precision dt

      petrun = 1

      dt = deltim()

c      call dattim('now', nowtime)
      mo = nowtime(2)
      day = nowtime(3)
      jday = julian('now', 'calendar')
      basin_potet = 0.

C******Set switch for active transpiration period
c      at the start of each new day (last day, lday.ne.day)
c

      if(lday.ne.day) then

        lday = day

        if(getvar('temp', 'tmax_c', nmru, 'real', tmax_c)
     +   .ne.0) return

        do 10 i= 1,nmru

C******If in checking period, then for each day
C******sum max temp until greater than temperature index parameter,
C******at which time, turn transpiration switch on, check switch off

           if(transp_check(i).eq.1) then
              if(tmax_c(i).gt.0.) tmax_sum(i) = tmax_sum(i) + tmax_c(i)

              if(tmax_sum(i).gt.transp_tmax_c(i)) then
                 transp_on(i) = 1
                 transp_check(i) = 0
                 tmax_sum(i) = 0.
              end if
              
C******Otherwise, check for month to turn check switch on or
C******transpiration switch off

         else
          if(day.eq.1) then
               if(mo.eq.transp_beg(i)) then
                  transp_check(i) = 1
                  if(tmax_c(i).gt.0.) tmax_sum(i) = tmax_sum(i) 
     +                 + tmax_c(i)
C******If transpiration switch on, check for end of period
 
            else
                if(transp_on(i).eq.1) then
                  if(mo.eq.transp_end(i)) transp_on(i) = 0
                end if
              end if
            end if
         end if
   10   continue

C******Compute potential et for each mru using Hamon formulation

        if(getvar('soltab', 'mru_sunhrs', nmru, 'real',
     +     mru_sunhrs).ne.0) return

        if(getvar('temp', 'temp_c', nmru, 'real', tavgc)
     +    .ne.0) return


        do 20 i=1,nmru
c          ir = mru_radpl(i)
          dyl = mru_sunhrs(i)
          vpsat = 6.108*EXP(17.26939*tavgc(i)/(tavgc(i)+237.3))
          vdsat = 216.7*vpsat/(tavgc(i)+273.2)
          daily_potet(i) = hamon_coef(mo)*dyl*dyl*vdsat

c          write(87,8020) mo,day, i, vpsat, vdsat, dyl,daily_potet(i)
c 8020     format(3i4, 5f8.4)

   20   continue

      end if

      do 30 i=1,nmru
        if(dt.lt.24.) then
          factor = dt/24.
        else
          factor = 1.
        end if
        potet(i) = factor * daily_potet(i)
        basin_potet = basin_potet + potet(i) * mru_area(i) 
   30 continue

      basin_potet = basin_potet / basin_area

      petrun = 0

      return
      end

