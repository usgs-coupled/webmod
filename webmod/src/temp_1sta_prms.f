c*********************************************************************
c     temp_1sta.f: Distribute temperatures to MRU's using 1 station
c                  and monthly lapse rate parameters
c
c     04sep03 - Added RCS version control - RMTW
c     01mar05 - Removed basin_tsta and temp_units.
c               Basin temps will be area weighted
c     07mar05 - Made lapse rate units deg celsius per kilometer
c     06may10 - Port to Fortran 90 with module and dynamic memory
c ********************************************************************
      MODULE WEBMOD_TEMP1STA
      IMPLICIT NONE
      include 'fmodules.inc'

C   Dimensions
      integer, save :: nmru, nmonths, ntemp
      
C   Declared Variables
      real, save :: basin_tmax_f, basin_tmax_c
      real, save :: basin_tmin_f, basin_tmin_c
      real, save :: basin_temp_f, basin_temp_c
      real, save, allocatable :: tmax_f(:),tmax_c(:)
      real, save, allocatable :: tmin_f(:),tmin_c(:)
      real, save, allocatable :: temp_f(:),temp_c(:)
C   Declared Parameters
      real, save :: basin_area
      integer, save, allocatable ::  mru_tsta(:)
      real, save, allocatable :: tsta_elev(:)
      real, save, allocatable :: tmax_lapse(:),tmin_lapse(:)
      real, save, allocatable :: tmax_adj(:), tmin_adj(:)
      real, save, allocatable :: mru_elev(:), mru_area(:)

C   Declared Private Variables
      real, save, allocatable :: elfac(:), tcrn(:), tcrx(:), tcr(:)

C   Undeclared Static Variables gotten from from other modules - soltab
      real, save, allocatable :: tsta_max_c(:),tsta_min_c(:)  ! from obs
      real, save, allocatable :: tsta_temp_c(:)               ! from obs
      
      END MODULE WEBMOD_TEMP1STA
c***********************************************************************
c***********************************************************************
c
c     Main temp_1sta routine
c

      integer function temp_1sta_prms(arg)

      character*(*) arg
      CHARACTER*256 SVN_ID
      integer t1decl, t1init, t1run, t1clean
      save SVN_ID

      SVN_ID = 
     $     '$Id: temp_1sta_prms.f 29 2006-07-06 23:03:45Z rmwebb $ '

      temp_1sta_prms = 0

      if(arg.eq.'declare') then
        temp_1sta_prms = t1decl()

      else if(arg.eq.'initialize') then
        temp_1sta_prms = t1init()

      else if(arg.eq.'run') then
        temp_1sta_prms = t1run() 

      end if

C******Debug level print
      call dpint4('End of temp_1sta, retval = ', temp_1sta_prms, 1, 2)

      end function temp_1sta_prms


c***********************************************************************
c 
c     t1decl - set up parameters for temperature computations
c

      integer function t1decl()


      USE WEBMOD_TEMP1STA

      t1decl = 1

!
! Get dimensions
!
      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      ntemp = getdim('ntemp')
        if ( ntemp.eq.-1 ) return
      nmonths = getdim('nmonths')
        if ( nmonths.eq.-1 ) return

      ALLOCATE (tmax_f(nmru))
      if(declvar('temp', 'tmax_f', 'nmru', nmru, 'real',
     +     'MRU adjusted daily maximum temperature',
     +     'degrees F',
     +   tmax_f).ne.0) return

      ALLOCATE (tmin_f(nmru))
      if(declvar('temp', 'tmin_f', 'nmru', nmru, 'real',
     +     'MRU adjusted daily minimum temperature',
     +     'degrees F',
     +   tmin_f).ne.0) return

      ALLOCATE (temp_f(nmru))
      if(declvar('temp', 'temp_f', 'nmru', nmru, 'real',
     +     'MRU adjusted temperature for timestep',
     +     'degrees F',
     +   temp_f).ne.0) return
 
      ALLOCATE (tmax_c(nmru))
      if(declvar('temp', 'tmax_c', 'nmru', nmru, 'real',
     +     'MRU adjusted daily maximum temperature',
     +     'degrees C',
     +   tmax_c).ne.0) return

      ALLOCATE (tmin_c(nmru))
      if(declvar('temp', 'tmin_c', 'nmru', nmru, 'real', 
     +     'MRU adjusted daily minimum temperature',
     +     'degrees C',
     +   tmin_c).ne.0) return

      ALLOCATE (temp_c(nmru))
      if(declvar('temp', 'temp_c', 'nmru', nmru, 'real',
     +     'MRU adjusted temperature for timestep',
     +     'degrees C',
     +   temp_c).ne.0) return
 
      if(declvar('temp', 'basin_tmax_c', 'one', 1, 'real',
     +     'Basin area-weighted daily maximum temperature',
     +     'degrees C',
     +   basin_tmax_c).ne.0) return

      if(declvar('temp', 'basin_tmin_c', 'one', 1, 'real', 
     +     'Basin area-weighted daily minimum temperature',
     +     'degrees C',
     +   basin_tmin_c).ne.0) return

      if(declvar('temp', 'basin_temp_c', 'one', 1, 'real',
     +     'Basin area-weighted temperature for timestep',
     +     'degrees C',
     +   basin_temp_c).ne.0) return

      if(declvar('temp', 'basin_tmax_f', 'one', 1, 'real',
     +     'Basin area-weighted daily maximum temperature',
     +     'degrees F',
     +   basin_tmax_f).ne.0) return

      if(declvar('temp', 'basin_tmin_f', 'one', 1, 'real', 
     +     'Basin area-weighted daily minimum temperature',
     +     'degrees F',
     +   basin_tmin_f).ne.0) return

      if(declvar('temp', 'basin_temp_f', 'one', 1, 'real',
     +     'Basin area-weighted temperature for timestep',
     +     'degrees F',
     +   basin_temp_f).ne.0) return

      ALLOCATE (tmax_lapse(nmonths))
      if(declparam('temp', 'tmax_lapse', 'nmonths', 'real',
     +   '9.8', '-15.0', '15.0',
     +   'Monthly maximum temperature lapse rate',
     +   'Array of twelve values representing the increase in '//
     +   'maximum temperature with decreasing elevation '//
     +   'for each month, January to December',
     +   'degrees celsius per kilometer')
     +    .ne.0) return

      ALLOCATE (tmin_lapse(nmonths))
      if(declparam('temp', 'tmin_lapse', 'nmonths', 'real',
     +   '9.8', '-15.0', '15.0',
     +   'Monthly minimum temperature lapse rate',
     +   'Array of twelve values representing the increase in '//
     +   'minimum temperature with decreasing elevation '//
     +   'for each month, January to December',
     +   'degrees celsius per kilometer')
     +    .ne.0) return

      ALLOCATE (tsta_elev(ntemp))
      if(declparam('temp', 'tsta_elev', 'ntemp', 'real',
     +   '0', '-300.', '10000.',
     +   'Temperature station elevation',
     +   'Elevation of each temperature '//
     +   'measurement station',
     +   'meters').ne.0) return

      ALLOCATE (tmax_adj(nmru))
      if(declparam('temp', 'tmax_adj', 'nmru', 'real',
     +   '0.0', '-10.', '10.0',
     +   'MRU maximum temperature adjustment',
     +   'Adjustment, in deg C, '//
     +   'to MRU maximum temperature '//
     +   'based on slope and aspect of MRU',
     +   'degrees C').ne.0) return

      ALLOCATE (tmin_adj(nmru))
      if(declparam('temp', 'tmin_adj', 'nmru', 'real',
     +   '0.0', '-10.0', '10.0',
     +   'MRU minimum temperature adjustment',
     +   'Adjustment, in deg C '//
     +   'to MRU minimum temperature '//
     +   'based on slope and aspect of MRU',
     +   'degrees C').ne.0) return

      ALLOCATE (mru_tsta(nmru))
      if(declparam('temp', 'mru_tsta', 'nmru', 'integer',
     +   '1', 'bounded', 'ntemp',
     +   'Index of temperature station for MRU',
     +   'Index of temperature station used to compute '//
     +   'MRU temperatures',
     +   'none').ne.0) return

      ALLOCATE (mru_area(nmru))
      if(declparam('temp_1sta', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('temp_1sta', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return
      
      ALLOCATE (mru_elev(nmru))
      if(declparam('temp_1sta', 'mru_elev', 'nmru', 'real',
     +   '0.', '-300.', '10000',
     +   'Mean elevation for each MRU',
     +   'Mean elevation for each MRU',
     +   'meters').ne.0) return
!
! Allocate local variables
!
      ALLOCATE (elfac(nmru))
      ALLOCATE (tcrn(nmru))
      ALLOCATE (tcrx(nmru))
      ALLOCATE (tcr(nmru))
!
! Allocate copies of variables from other modules
!
      ALLOCATE (tsta_max_c(ntemp))  ! from obs
      ALLOCATE (tsta_min_c(ntemp))  ! from obs
      ALLOCATE (tsta_temp_c(ntemp))  ! from obs
     
      t1decl = 0

      return
      end

c***********************************************************************
c
c     t1init - Initialize temp_1sta module - get parameter values,
c                compute elfac
c

      integer function t1init()

      USE WEBMOD_TEMP1STA
      
      integer j, k

      t1init = 1

      if(getparam('temp', 'tmin_lapse', nmonths, 'real', tmin_lapse)
     +   .ne.0) return

      if(getparam('temp', 'tmax_lapse', nmonths, 'real', tmax_lapse)
     +   .ne.0) return

      if(getparam('temp', 'tsta_elev', ntemp, 'real', tsta_elev)
     +   .ne.0) return

      if(getparam('basin', 'mru_elev', nmru, 'real', mru_elev)
     +   .ne.0) return

      if(getparam('temp', 'tmax_adj', nmru, 'real', tmax_adj)
     +   .ne.0) return

      if(getparam('temp', 'tmin_adj', nmru, 'real', tmin_adj)
     +   .ne.0) return

      if(getparam('temp', 'mru_tsta', nmru, 'integer', mru_tsta)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return


      do 10 j=1,nmru
        k = mru_tsta(j)
        elfac(j) = (mru_elev(j)-tsta_elev(k))/1000.
   10 continue

      t1init = 0

      return
      end

c***********************************************************************
c
c     t1run - Computes maximum, minumum and average temperature
c               for each MRU based on monthly lapse rate
c

      integer function t1run()
 
      USE WEBMOD_TEMP1STA
      
      integer j,mo,k,yr
      integer nowtime(6), nstep, day
      
      double precision dt

      t1run = 1

      call dattim('now',nowtime)
      nstep = getstep()
      day = nowtime(3)
      mo = nowtime(2)
      yr = nowtime(1)
      dt = deltim()

      if(getvar('obs', 'tsta_max_c', ntemp, 'real', tsta_max_c)
     +  .ne.0) return 
      if(getvar('obs', 'tsta_min_c', ntemp, 'real', tsta_min_c)
     +  .ne.0) return
      if(getvar('obs', 'tsta_temp_c', ntemp, 'real', tsta_temp_c)
     +     .ne.0) return

      basin_tmax_c = 0.
      basin_tmin_c = 0.
      basin_temp_c = 0.

      if(nstep.eq.1.or.day.eq.1) then
        do 10 j= 1,nmru
          tcrx(j) = (tmax_lapse(mo)*elfac(j))-tmax_adj(j)
          tcrn(j) = (tmin_lapse(mo)*elfac(j))-tmin_adj(j)
          tcr(j) = (tcrx(j) + tcrn(j))/2.
   10   continue
      end if

      do 20 j = 1,nmru
        k = mru_tsta(j)
        tmax_c(j) = tsta_max_c(k)-tcrx(j)
        tmin_c(j) = tsta_min_c(k)-tcrn(j)
        temp_c(j) = tsta_temp_c(k) -tcr(j)
        basin_tmax_c = basin_tmax_c + (tmax_c(j) * mru_area(j))
        basin_tmin_c = basin_tmin_c + (tmin_c(j) * mru_area(j))
        basin_temp_c = basin_temp_c + (temp_c(j) * mru_area(j))
        tmax_f(j) = tmax_c(j)*1.8+32.
        tmin_f(j) = tmin_c(j)*1.8+32.
        temp_f(j) = temp_c(j)*1.8+32.
   20 continue

      basin_tmax_c = basin_tmax_c / basin_area
      basin_tmin_c = basin_tmin_c / basin_area
      basin_temp_c = basin_temp_c / basin_area
      basin_tmax_f = basin_tmax_c*1.8+32.
      basin_tmin_f = basin_tmin_c*1.8+32.
      basin_temp_f = basin_temp_c*1.8+32.

      t1run = 0

      return
      end


