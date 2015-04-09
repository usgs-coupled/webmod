c***********************************************************************
c     intcp.f:  Computes amount of intercepted rain and snow, and
c               determines evaporation from interception 
c
c               8 Sept 2003 - Added RCS version control - RMTW
c        06may10 - Port to Fortran 90 with module and dynamic memory
c
c***********************************************************************
c 
      MODULE WEBMOD_INTCP
      USE mmf, only : declparam
      IMPLICIT NONE
      include 'fmodules.inc'

C   Dimensions
      integer, save :: nmru, nevap, nmonths
      
C   Declared Variables
      integer, save, allocatable:: int_snow(:)
      integer, save, allocatable :: intcp_form(:), intcp_on(:)
      real, save :: basin_net_dep, basin_intcp_stor, basin_intcp_evap
      real, save, allocatable :: net_rain(:), net_snow(:), net_dep(:)
      real, save, allocatable :: intcp_evap(:), intcp_stor(:)
C   Declared Parameters
      real, save :: basin_area
      integer, save, allocatable:: cov_type(:)
      real, save, allocatable :: snow_intcp(:), srain_intcp(:)
      real, save, allocatable :: covden_sum(:), covden_win(:)
      real, save, allocatable :: epan_coef(:),mru_area(:),wrain_intcp(:)

C   Private Variables
      integer, save,  allocatable:: intcp_transp_on(:)
C   Undeclared Static Variables gotten from from other modules
      integer, save :: basin_dep
      integer, save, allocatable:: transp_on(:)
      real, save, allocatable :: mru_dep(:), mru_perv(:), mru_rain(:)
      real, save, allocatable :: mru_snow(:), pan_evap(:)
      real, save, allocatable :: pkwater_equiv(:), potet(:), swrad(:)
      real, save, allocatable :: temp_c(:)
    
      END MODULE WEBMOD_INTCP
c***********************************************************************
c***********************************************************************
c
c     Main intcp routine
c

      integer function intcp_prms(arg)
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID

      integer intdecl, intinit, intrun

      SVN_ID = 
     $     '$Id: intcp_prms.f 29 2006-07-06 23:03:45Z rmwebb $ '
    
      intcp_prms = 0

      if(arg.eq.'declare') then
        intcp_prms = intdecl()

      else if(arg.eq.'initialize') then
        intcp_prms = intinit()

      else if(arg.eq.'run') then
        intcp_prms = intrun()

      end if

      return
      end

c***********************************************************************
c 
c     intdecl - set up parameters for interception computations
c

      integer function intdecl()

      USE WEBMOD_INTCP
      IMPLICIT NONE

      intdecl = 1

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nevap = getdim('nevap')
        if ( nevap.eq.-1 ) return
      nmonths = getdim('nmonths')
        if ( nmonths.eq.-1 ) return

      ALLOCATE (net_rain(Nmru))
      if(declvar('intcp', 'net_rain', 'nmru', nmru, 'real',
     +     'mru_rain minus interception',
     +     'inches',
     +   net_rain).ne.0) return

      ALLOCATE (net_snow(Nmru))
      if(declvar('intcp', 'net_snow', 'nmru', nmru, 'real',
     +     'mru_snow minus interception',
     +     'inches',
     +   net_snow).ne.0) return

      ALLOCATE (net_dep(Nmru))
      if(declvar('intcp', 'net_dep', 'nmru', nmru, 'real',
     +     'MRU deposition (ppt and irrigation) with '//
     +     'interception removed',
     +     'inches',
     +   net_dep).ne.0) return

      if(declvar('intcp', 'basin_net_dep', 'one', 1, 'real',
     +     'Basin area-weighted average net_dep',
     +     'inches',
     +   basin_net_dep).ne.0) return

      ALLOCATE (intcp_stor(Nmru))
      if(declvar('intcp', 'intcp_stor', 'nmru', nmru, 'real',
     +     'Current interception storage on an MRU',
     +     'inches',
     +   intcp_stor).ne.0) return

      if(declvar('intcp', 'basin_intcp_stor', 'one', 1, 'real',
     +     'Basin area-weighted average interception storage',
     +     'inches',
     +   basin_intcp_stor).ne.0) return

      ALLOCATE (intcp_evap(Nmru))
      if(declvar('intcp', 'intcp_evap', 'nmru', nmru, 'real',
     +     'Evaporation from interception on each MRU',
     +     'inches',
     +   intcp_evap).ne.0) return

      if(declvar('intcp', 'basin_intcp_evap', 'one', 1, 'real',
     +     'Basin area-weighted evaporation from interception',
     +     'inches',
     +   basin_intcp_evap).ne.0) return

      ALLOCATE (intcp_form(Nmru))
      if(declvar('intcp', 'intcp_form', 'nmru', nmru, 'integer',
     +     'Form (rain or snow) of interception',
     +     'none',
     +   intcp_form).ne.0) return

      ALLOCATE (intcp_on(Nmru))
      if(declvar('intcp', 'intcp_on', 'nmru', nmru, 'integer',
     +     'Whether there is interception in the canopy, 0=no'//
     +     ' 1=yes ',
     +     'none',
     +   intcp_on).ne.0) return      

      ALLOCATE (int_snow(Nmru))
      if(declvar('intcp', 'int_snow', 'nmru', nmru, 'integer',
     +     'Whether snow has fallen from the canopy, 0=no'//
     +     ' 1=yes ',
     +     'none',
     +   int_snow).ne.0) return


      ALLOCATE (epan_coef(Nmonths))
      if(declparam('intcp', 'epan_coef', 'nmonths', 'real',
     +   '1.0', '0.2', '3.0',
     +   'Evaporation pan coefficient',
     +   'Evaporation pan coefficient',
     +   'none').ne.0) return

      ALLOCATE (mru_area(Nmru))
      if(declparam('intcp', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('intcp', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return
      
      ALLOCATE (snow_intcp(Nmru))
      if(declparam('intcp', 'snow_intcp', 'nmru', 'real',
     +   '.1', '0.', '5.',
     +   'Snow interception storage capacity',
     +   'Snow interception storage capacity for the major '//
     +   'vegetation type in the MRU',
     +   'inches')
     +   .ne.0) return
      
      ALLOCATE (srain_intcp(Nmru))
      if(declparam('intcp', 'srain_intcp', 'nmru', 'real',
     +   '.1', '0.', '5.',
     +   'Summer rain interception storage capacity',
     +   'Summer rain interception storage capacity for the major '//
     +   'vegetation type in the MRU',
     +   'inches')
     +   .ne.0) return
      
      ALLOCATE (wrain_intcp(Nmru))
      if(declparam('intcp', 'wrain_intcp', 'nmru', 'real',
     +   '.1', '0.', '5.',
     +   'Winter rain interception storage capacity',
     +   'Winter rain interception storage capacity for the major '//
     +   'vegetation type in the MRU',
     +   'inches')
     +   .ne.0) return

      ALLOCATE (cov_type(Nmru))
      if(declparam('intcp', 'cov_type', 'nmru', 'integer',
     +   '3', '0', '3',
     +   'Cover type designation for MRU',
     +   'Vegetation cover type designation for MRU:  '//
     +   '0 = bare soil, 1 = grasses, 2 = shrubs, 3 = trees',
     +   'none')
     +   .ne.0) return

      ALLOCATE (covden_sum(Nmru))
      if(declparam('intcp', 'covden_sum', 'nmru', 'real',
     +   '.5', '0.', '1.0',
     +   'Summer vegetation cover density for major vegetation type',
     +   'Summer vegetation cover density for the major '//
     +   'vegetation type on each MRU',
     +   'decimal percent')
     +   .ne.0) return 

      ALLOCATE (covden_win(Nmru))
      if(declparam('intcp', 'covden_win', 'nmru', 'real',
     +   '.5', '0.', '1.0',
     +   'Winter vegetation cover density for major vegetation type',
     +   'Winter vegetation cover density for the major '//
     +   'vegetation type on each MRU',
     +   'decimal percent')
     +   .ne.0) return 

!
! getvars from other modules
!
      ALLOCATE (mru_rain(Nmru))
      ALLOCATE (mru_snow(Nmru))
      ALLOCATE (mru_dep(Nmru))
      ALLOCATE (transp_on(Nmru))
      ALLOCATE (temp_c(Nmru))
      ALLOCATE (swrad(Nmru))
      ALLOCATE (mru_perv(Nmru))
      ALLOCATE (pan_evap(Nevap))
      ALLOCATE (potet(Nmru))
      ALLOCATE (pkwater_equiv(Nmru))
!
! Private variable
!
      ALLOCATE (intcp_transp_on(nmru))

      intdecl = 0

      return
      end

c***********************************************************************
c
c     intinit - Initialize intcp module - get parameter values,
c               set initial values.
c

      integer function intinit()

      USE WEBMOD_INTCP
      IMPLICIT NONE
      integer i

      intinit = 1

      if(getparam('intcp', 'snow_intcp', nmru, 'real', snow_intcp)
     +   .ne.0) return

      if(getparam('intcp', 'wrain_intcp', nmru, 'real', wrain_intcp)
     +   .ne.0) return

      if(getparam('intcp', 'srain_intcp', nmru, 'real', srain_intcp)
     +   .ne.0) return
 
      if(getparam('intcp', 'cov_type', nmru, 'integer', cov_type)
     +   .ne.0) return

      if(getparam('intcp', 'covden_sum', nmru, 'real', covden_sum)
     +   .ne.0) return 

      if(getparam('intcp', 'covden_win', nmru, 'real', covden_win)
     +   .ne.0) return 

      if(getparam('intcp', 'epan_coef', nmonths, 'real', epan_coef)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return

      do 10 i=1,nmru
        intcp_stor(i) = 0.
        intcp_on(i) = 0
        intcp_form(i) = 0
        intcp_evap(i) = 0.
   10 continue

      intinit = 0

      return
      end

c***********************************************************************
c
c     intrun - Computes and keeps track of intercepted precipitation
c              and evaporation for each MRU 
c

      integer function intrun()

      USE WEBMOD_INTCP
      IMPLICIT NONE

      integer i, nstep, mo, nowtime(6)
      real cov, stor, evcan, ta, sw, tc, air, can, bal, xmlt, z, xlos, d
      double precision dt

      if(getvar('precip', 'mru_rain', nmru, 'real', mru_rain)
     +  .ne.0) return

      if(getvar('precip', 'mru_snow', nmru, 'real', mru_snow)
     +   .ne.0) return

      if(getvar('precip', 'mru_dep', nmru, 'real', mru_dep)
     +   .ne.0) return

      if(getvar('precip', 'basin_dep', 1, 'real', basin_dep)
     +   .ne.0) return

      if(getvar('potet', 'transp_on', nmru, 'integer', transp_on)
     +   .ne.0) return

      if(getvar('temp', 'temp_c', nmru, 'real', temp_c)
     +   .ne.0) return

      if(getvar('solrad', 'swrad', nmru, 'real', swrad)
     +   .ne.0) return

      if(getvar('basin', 'mru_perv', nmru, 'real', mru_perv)
     +   .ne.0) return
     
      if(getvar('obs', 'pan_evap', nevap, 'real', pan_evap)
     +   .ne.0) return

      if(getvar('potet', 'potet', nmru, 'real', potet)
     +   .ne.0) return

      if(getvar('snow', 'pkwater_equiv', nmru, 'real', pkwater_equiv)
     +   .ne.0) return


      call dattim('now', nowtime)
      mo = nowtime(2)

      dt = deltim()
 
      nstep = getstep()
      if(nstep.eq.1) then
        do 10 i = 1,nmru
          intcp_transp_on(i) = transp_on(i)
   10   continue
      end if

      basin_net_dep = 0.
      basin_intcp_stor = 0.
      basin_intcp_evap = 0.
 
      do 1000 i= 1,nmru

C******Adjust interception amounts for changes in summer/winter cover
C******density

        int_snow(i) = 0
        cov = 0
        intcp_form(i) = 0      ! This was not set before and so lagged into next day - RMTW

c first step of non transpiration period

        if(transp_on(i).eq.0.and.intcp_transp_on(i).eq.1) then
           intcp_transp_on(i) = 0 
           if(intcp_stor(i).gt.0.) then
              intcp_evap(i) = intcp_stor(i) *
     +             (covden_sum(i)-covden_win(i))
           end if

c first step of transpiration period

        else if(transp_on(i).eq.1.and.intcp_transp_on(i).eq.0) then 
           intcp_transp_on(i) = 1
           if(intcp_stor(i).gt.0.) then
              if(covden_sum(i).eq.0.) then
                 intcp_stor(i) = 0.
              else
                 intcp_stor(i) = intcp_stor(i) *
     +                (covden_win(i)/covden_sum(i))
              end if
           end if
        else
           intcp_evap(i) = 0.
        end if

******Determine the amount of interception from rain

        if(mru_rain(i).gt.0.) then
           intcp_form(i) = 0
           if(cov_type(i).gt.1.or.
     +          (cov_type(i).eq.1.and.pkwater_equiv(i).eq.0.)) then
              if(transp_on(i).eq.1) then
                 cov = covden_sum(i)
                 stor = srain_intcp(i)
              else
                 cov = covden_win(i)
                 stor = wrain_intcp(i)
              end if
              call intercept(mru_rain(i), stor, cov, intcp_on(i),
     +             intcp_stor(i), net_rain(i))
           else
              cov = 0.
              intcp_stor(i) = 0.
              net_rain(i) = mru_rain(i)
              intcp_on(i) = 0
           end if
      else
         net_rain(i) = 0.
      end if
c******Determine amount of interception from snow

        if(mru_snow(i).gt.0.) then
           intcp_form(i) = 1
           if(cov_type(i).gt.1) then
              cov = covden_win(i)
              if(transp_on(i).eq.1) cov = covden_sum(i)
              stor = snow_intcp(i)
              call intercept(mru_snow(i), stor, cov, intcp_on(i),
     +             intcp_stor(i), net_snow(i))
           else
              cov = 0.
              intcp_stor(i) = 0.
              net_snow(i) = mru_snow(i)
              intcp_on(i) = 0
           end if
        else
           net_snow(i) = 0.
        end if

        net_dep(i) = net_rain(i) + net_snow(i)

C******compute evaporation or sublimation of interception

        if(intcp_on(i).eq.1) then

          intcp_evap(i) = 0.
          if(dt.lt.24.. and . mru_dep(i).gt.0.) go to 900

          evcan = potet(i) / epan_coef(mo)
          if(nevap.gt.0) then
            if(pan_evap(1).ne.-999) evcan = pan_evap(1)
         end if
            
c******Compute snow interception loss

          if(intcp_form(i).eq.1) then
            if(basin_dep.le.0.) then
              cov = covden_win(i)
              if(transp_on(i).eq.1) cov = covden_sum(i)
              ta = temp_c(i)
              if(ta.ge.0.) then
                sw = swrad(i) * .19
                tc = .585e-7 * ((ta +273.16)**4)
                air = tc*(1.-cov)*.85
                can = cov * tc
                bal = sw + air + can
                if(bal.gt.0.) then
                   xmlt = bal/203.2
                   z = intcp_stor(i) - xmlt
                   if(z.gt.0.) then
                      intcp_on(i) = 1
                      intcp_stor(i) = z
                      xlos = xmlt
                   else
                      xlos = intcp_stor(i)
                      intcp_stor(i) = 0.
                      intcp_on(i) = 0 
                   end if
                   d = evcan - xlos
                   if(d.lt.0) then
                      intcp_evap(i) = intcp_evap(i) + evcan
                      net_snow(i) = -d * cov 
                      net_dep(i) = net_snow(i)
                      int_snow(i) = 1
                   else
                      intcp_evap(i) = intcp_evap(i) + xlos
                   end if
                end if
c                intcp_on(i) = 1
             end if
          end if

       else if(intcp_form(i).eq.0) then
            cov = covden_sum(i)
            if(transp_on(i).eq.0) cov = covden_win(i)
            d = evcan - intcp_stor(i)
            if(d.lt.0.) then
              intcp_stor(i) = -d
              intcp_evap(i) = intcp_evap(i) + evcan
              intcp_on(i) = 1
            else
              intcp_evap(i) = intcp_evap(i) + intcp_stor(i)
              intcp_stor(i) = 0.
              intcp_on(i) = 0
           end if 
        end if
      end if 

 900    basin_net_dep = basin_net_dep + (net_dep(i) * mru_area(i))
        basin_intcp_stor = basin_intcp_stor + 
     +                     (cov * intcp_stor(i) * mru_area(i))
        basin_intcp_evap = basin_intcp_evap + 
     +                     (cov * intcp_evap(i) * mru_area(i))
        intcp_evap(i) = cov * intcp_evap(i)


 1000 continue

  
      basin_net_dep = basin_net_dep / basin_area
      basin_intcp_stor = basin_intcp_stor / basin_area
      basin_intcp_evap = basin_intcp_evap / basin_area

      intrun = 0

      return
      end


c***********************************************************************
c      Subroutine to compute interception of rain or snow
c***********************************************************************

      subroutine intercept(precip, stor_max, cov, intcp_on, intcp_stor,
     +                     net_precip)

      integer intcp_on
      real precip, stor_max, cov, intcp_stor, net_precip, avail_stor
      real thrufall

      intcp_on = 1
      avail_stor = stor_max - intcp_stor
      if(precip.gt.avail_stor) then
        intcp_stor = stor_max
        thrufall = precip - avail_stor
      else
        intcp_stor = intcp_stor + precip
        thrufall = 0.
      end if
      net_precip = (precip*(1.-cov)) + (thrufall * cov)

      return
      end
