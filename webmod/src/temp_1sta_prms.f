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

C Constant
      real, PARAMETER :: nearzero = 1.E-10 ! Trap when trxn_**_days is actually equal to an integer
C   Dimensions
      integer, save :: nmru, nmonths, ntemp
      
C   Declared Variables
      real, save :: basin_tmax_f, basin_tmax_c
      real, save :: basin_tmin_f, basin_tmin_c
      real, save :: basin_temp_f, basin_temp_c
      real, save, allocatable :: tmax_f(:),tmax_c(:)
      real, save, allocatable :: tmin_f(:),tmin_c(:)
      real, save, allocatable :: trxn_ohoriz_c(:),trxn_uz_c(:)
      real, save, allocatable :: trxn_sat_c(:)
      real, save, allocatable :: temp_f(:),temp_c(:)
C   Declared Parameters
      real, save :: basin_area
      real, save :: trxn_ohoriz_days
      integer, save, allocatable :: trxn_ohoriz_stat(:)
      real, save, allocatable :: trxn_ohoriz_c_adj(:)
      real, save :: trxn_uz_days
      integer, save, allocatable :: trxn_uz_stat(:)
      real, save, allocatable :: trxn_uz_c_adj(:)
      real, save :: trxn_sat_days
      integer, save, allocatable :: trxn_sat_stat(:)
      real, save, allocatable :: trxn_sat_c_adj(:)
      integer, save, allocatable ::  mru_tsta(:)
      real, save, allocatable :: tsta_elev(:)
      real, save, allocatable :: tmax_lapse(:),tmin_lapse(:)
      real, save, allocatable :: tmax_adj(:), tmin_adj(:)
      real, save, allocatable :: mru_elev(:), mru_area(:)

C   Private Variables
      real, save, allocatable :: elfac(:), tcrn(:), tcrx(:), tcr(:)
      integer,save :: it_oh_days, it_uz_days, it_sat_days ! integer >= to real days. CEILING()
      real, save, allocatable :: trxn_c_array_oh(:,:)
      real, save, allocatable :: trxn_c_array_uz(:,:)
      real, save, allocatable :: trxn_c_array_sat(:,:)
      real :: trxn_frac ! it_**_days - trxn_**_days
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
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      !integer t1decl, t1init, t1run, t1clean
      integer t1decl, t1init, t1run
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
      IMPLICIT NONE

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

!
! parameters controlling moving averages of temperature to use for reactions in the hillslope
!
      ALLOCATE (trxn_ohoriz_stat(nmru))
      if(declparam('temp', 'trxn_ohoriz_stat', 'nmru', 'integer',
     +   '2', '1', '3', 
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average',
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average',
     +   'none').ne.0) return

      ALLOCATE (trxn_uz_stat(nmru))
      if(declparam('temp', 'trxn_uz_stat', 'nmru', 'integer',
     +   '2', '1', '3', 
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average', 
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average', 
     +   'none').ne.0) return
 
      ALLOCATE (trxn_sat_stat(nmru))
      if(declparam('temp', 'trxn_sat_stat', 'nmru', 'integer',
     +   '2', '1', '3', 
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average', 
     +   'Use tmin(1), tavg(2), or tmax(3) for moving average', 
     +   'none').ne.0) return
      
      ALLOCATE (trxn_ohoriz_c_adj(nmru))
      if(declparam('temp', 'trxn_ohoriz_c_adj', 'nmru', 'real',
     +  '0.0', '-5.0', '5.0',
     +  'Temperature adjustment to moving average of tmin, '//
     +  'tavg, or tmax for O-horizon','Temperature adjustment '//
     +  'to moving average of tmin, tavg, or tmax for O-horizon',
     +  'degrees C').ne.0) return

      ALLOCATE (trxn_uz_c_adj(nmru))
      if(declparam('temp', 'trxn_uz_c_adj', 'nmru', 'real',
     +  '0.0', '-5.0', '5.0',
     +  'Temperature adjustment to moving average of tmin, '//
     +  'tavg, or tmax for UZ','Temperature adjustment to moving '//
     +  'average of tmin, tavg, or tmax for UZ',
     +  'degrees C').ne.0) return

      ALLOCATE (trxn_sat_c_adj(nmru))
      if(declparam('temp', 'trxn_sat_c_adj', 'nmru', 'real',
     +  '0.0', '-5.0', '5.0',
     +  'Temperature adjustment to moving average of tmin, '//
     +  'tavg, or tmax for saturated zone','Temperature adjustment '//
     +  'to moving average of tmin, tavg, or tmax for saturated zone',
     +  'degrees C').ne.0) return

      if(declparam('temp', 'trxn_ohoriz_days', 'one', 'real',
     +   '7.0', '1.0', '10000.',
     +   'Temperature averaging window for O-horizon',
     +   'Temperature averaging window for O-horizon',
     +   'days')
     +    .ne.0) return
      
      if(declparam('temp', 'trxn_uz_days', 'one', 'real',
     +   '90.', '1.', '10000.',
     +   'Temperature averaging window for UZ',
     +   'Temperature averaging window for UZ',
     +   'days')
     +    .ne.0) return
      
      if(declparam('temp', 'trxn_sat_days', 'one', 'real',
     +   '365.', '1.', '10000.',
     +   'Temperature averaging window for saturated zone',
     +   'Temperature averaging window for saturated zone',
     +   'days')
     +    .ne.0) return
      
      ALLOCATE (trxn_ohoriz_c(nmru))
      if(declvar('temp', 'trxn_ohoriz_c', 'nmru', nmru, 'real', 
     +   'Reaction temperature for O-horizon','degrees C',
     +   trxn_ohoriz_c).ne.0) return

      ALLOCATE (trxn_uz_c(nmru))
      if(declvar('temp', 'trxn_uz_c', 'nmru', nmru, 'real', 
     +   'Reaction temperature for UZ','degrees C',
     +   trxn_uz_c).ne.0) return

      ALLOCATE (trxn_sat_c(nmru))
      if(declvar('temp', 'trxn_sat_c', 'nmru', nmru, 'real', 
     +   'Reaction temperature for saturated zone','degrees C',
     +   trxn_sat_c).ne.0) return
!
! Parameters from other modules
!
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
      IMPLICIT NONE
      
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

      if(getparam('temp', 'trxn_ohoriz_days', 1, 'real',
     +   trxn_ohoriz_days).ne.0) return

      if(getparam('temp', 'trxn_uz_days', 1, 'real',
     +   trxn_uz_days).ne.0) return

      if(getparam('temp', 'trxn_sat_days', 1, 'real',
     +   trxn_sat_days).ne.0) return

      if(getparam('temp', 'trxn_ohoriz_stat', nmru, 'integer',
     +   trxn_ohoriz_stat).ne.0) return

      if(getparam('temp', 'trxn_uz_stat', nmru, 'integer',
     +   trxn_uz_stat).ne.0) return

      if(getparam('temp', 'trxn_sat_stat', nmru, 'integer',
     +   trxn_sat_stat).ne.0) return

      if(getparam('temp', 'trxn_ohoriz_c_adj', nmru, 'real',
     +   trxn_ohoriz_c_adj).ne.0) return

      if(getparam('temp', 'trxn_uz_c_adj', nmru, 'real',
     +   trxn_uz_c_adj).ne.0) return

      if(getparam('temp', 'trxn_sat_c_adj', nmru, 'real',
     +   trxn_sat_c_adj).ne.0) return

      
      it_oh_days=ceiling(trxn_ohoriz_days)
      it_uz_days=ceiling(trxn_uz_days)
      it_sat_days=ceiling(trxn_sat_days)
      ALLOCATE(trxn_c_array_oh(nmru,it_oh_days))
      ALLOCATE(trxn_c_array_uz(nmru,it_uz_days))
      ALLOCATE(trxn_c_array_sat(nmru,it_sat_days))
      trxn_c_array_oh = 0.0
      trxn_c_array_uz = 0.0
      trxn_c_array_sat = 0.0

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
c             for each MRU based on monthly lapse rate. Also computes
c             running averages of tmin, tavg, or tmax to use in assigning
c             reaction temperatures to the o-horizon, the UZ, and the
c             saturated zone.

      integer function t1run()
 
      USE WEBMOD_TEMP1STA
      IMPLICIT NONE
      
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
! ohorizon reaction temperatures
        select case (trxn_ohoriz_stat(j))
        case (1)
          trxn_c_array_oh(j,1) = tmin_c(j)
        case (2)
          trxn_c_array_oh(j,1) = temp_c(j)
        case (3)
          trxn_c_array_oh(j,1) = tmax_c(j)
        case default
          print*, 'trxn_ohoriz_stat for MRU ',j,' is outside the'//
     $            ' range 1-3. Run terminated. Correct and restart.'
          return
        end select
        trxn_ohoriz_c(j) = 0.0
        do k = 1, it_oh_days
            trxn_ohoriz_c(j) = trxn_ohoriz_c(j) + trxn_c_array_oh(j,k)
        end do
! reduce sum by fraction of oldest temperature not to be included in average
        trxn_frac = it_oh_days-trxn_ohoriz_days
        if(trxn_frac.gt.nearzero) then ! will not execute if trxn_**_days is very close to an integer
            trxn_ohoriz_c(j) = trxn_ohoriz_c(j) - 
     +      trxn_frac*trxn_c_array_oh(j,it_oh_days)
        endif
        trxn_ohoriz_c(j) =trxn_ohoriz_c(j)/trxn_ohoriz_days+
     +                    trxn_ohoriz_c_adj(j)
! UZ reaction temperatures
        select case (trxn_uz_stat(j))
        case (1)
          trxn_c_array_uz(j,1) = tmin_c(j)
        case (2)
          trxn_c_array_uz(j,1) = temp_c(j)
        case (3)
          trxn_c_array_uz(j,1) = tmax_c(j)
        case default
          print*, 'trxn_uz_stat for MRU ',j,' is outside the'//
     $            ' range 1-3. Run terminated. Correct and restart.'
          return
        end select
        trxn_uz_c(j) = 0.0
        do k = 1, it_uz_days
            trxn_uz_c(j) = trxn_uz_c(j) + trxn_c_array_uz(j,k)
        end do
! reduce sum by fraction of oldest temperature not to be included in average
        trxn_frac = it_uz_days-trxn_uz_days
        if(trxn_frac.gt.nearzero) then ! will not execute if trxn_**_days is very close to an integer
            trxn_uz_c(j) = trxn_uz_c(j) - 
     +      trxn_frac*trxn_c_array_uz(j,it_uz_days)
        endif
        trxn_uz_c(j) =trxn_uz_c(j)/trxn_uz_days+
     +                    trxn_uz_c_adj(j)
! saturated zone reaction temperatures
        select case (trxn_sat_stat(j))
        case (1)
          trxn_c_array_sat(j,1) = tmin_c(j)
        case (2)
          trxn_c_array_sat(j,1) = temp_c(j)
        case (3)
          trxn_c_array_sat(j,1) = tmax_c(j)
        case default
          print*, 'trxn_sat_stat for MRU ',j,' is outside the'//
     $            ' range 1-3. Run terminated. Correct and restart.'
          return
        end select
        trxn_sat_c(j) = 0.0
        do k = 1, it_sat_days
            trxn_sat_c(j) = trxn_sat_c(j) + trxn_c_array_sat(j,k)
        end do
! reduce sum by fraction of oldest temperature not to be included in average
        trxn_frac = it_sat_days-trxn_sat_days
        if(trxn_frac.gt.nearzero) then ! will not execute if trxn_**_days is very close to an integer
            trxn_sat_c(j) = trxn_sat_c(j) - 
     +      trxn_frac*trxn_c_array_sat(j,it_sat_days)
        endif
        trxn_sat_c(j) =trxn_sat_c(j)/trxn_sat_days+
     +                    trxn_sat_c_adj(j)
!
        basin_tmax_c = basin_tmax_c + (tmax_c(j) * mru_area(j))
        basin_tmin_c = basin_tmin_c + (tmin_c(j) * mru_area(j))
        basin_temp_c = basin_temp_c + (temp_c(j) * mru_area(j))
        tmax_f(j) = tmax_c(j)*1.8+32.
        tmin_f(j) = tmin_c(j)*1.8+32.
        temp_f(j) = temp_c(j)*1.8+32.
   20       continue
            
      trxn_c_array_oh=cshift(trxn_c_array_oh,-1,2)
      trxn_c_array_uz=cshift(trxn_c_array_uz,-1,2)
      trxn_c_array_sat=cshift(trxn_c_array_sat,-1,2)
      basin_tmax_c = basin_tmax_c / basin_area
      basin_tmin_c = basin_tmin_c / basin_area
      basin_temp_c = basin_temp_c / basin_area
      basin_tmax_f = basin_tmax_c*1.8+32.
      basin_tmin_f = basin_tmin_c*1.8+32.
      basin_temp_f = basin_temp_c*1.8+32.

      t1run = 0

      return
      end


