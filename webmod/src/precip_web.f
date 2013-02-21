c***********************************************************************
c   precip_web.f: Determine form (rain, snow, mix) of precipitation 
c                 and distribute to mru's.
c***********************************************************************
c
c RMTW - 3-3-05  - Pruned from precip_irrig ver 1.2
c                  Now there are separate precip and irrigation
c                  routines so that XYZ integration is easier and
c                  more clear.
c
c         3-7-05 - Changing all temperature variables to both degC and
c                  degF. All temperature parameters will be in SI units
c                  of degrees C and meters.
c
c        06may10 - Port to Fortran 90 with module and dynamic memory
c
c***********************************************************************
c 
      MODULE WEBMOD_PRECIP
      IMPLICIT NONE
      include 'fmodules.inc'

C   Dimensions
      integer, save :: nmru, nmonths, nrain, nform
      
C   Declared Variables
      real, save :: basin_dep, basin_ppt
      real, save, allocatable :: mru_ppt(:)
      integer, save, allocatable :: newsnow(:), pptmix(:)
      real, save, allocatable :: mru_dep(:), mru_rain(:), mru_snow(:)
      real, save, allocatable :: prmx(:),strain_adj(:,:)
C   Declared Parameters
      real, save :: tmax_allsnow_c
      integer, save, allocatable ::  mru_psta(:)
      real, save, allocatable :: tmax_allrain_c(:)
      real, save, allocatable :: adjmix_rain(:), rain_adj(:,:) 
      real, save, allocatable :: snow_adj(:,:),  mru_area_frac(:)

C   Declared Private Variables
      integer, save :: step1

C   Undeclared Static Variables gotten from from other modules - soltab
      integer, save :: route_on   ! from obs
      integer, save ::  form_data(1)  ! from obs
      real, save, allocatable :: precip(:), tmax_c(:), tmin_c(:) ! from obs and temp
      
      END MODULE WEBMOD_PRECIP
c***********************************************************************
c
c     Main precip routine
c

      integer function precip_web(arg)

      include 'fmodules.inc'

      character*(*) arg
      CHARACTER*256 SVN_ID

      integer pptdecl, pptinit, pptrun
      save SVN_ID

      SVN_ID = 
     $     '$Id: precip_web.f 29 2006-07-06 23:03:45Z rmwebb $ '

      precip_web = 0

      if(arg.eq.'declare') then
         precip_web = pptdecl()

      else if(arg.eq.'initialize') then
        precip_web = pptinit()

      else if(arg.eq.'run') then
        precip_web = pptrun()

      end if
      
      return
      end
c
c     pptdecl - set up parameters for precipitation computations
c

      integer function pptdecl()

      USE WEBMOD_PRECIP

      pptdecl = 1
!
! Get dimensions
!
      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nrain = getdim('nrain')
        if ( nrain.eq.-1 ) return
      nform = getdim('nform')
        if ( nform.eq.-1 ) return
      nmonths = getdim('nmonths')
        if ( nmonths.eq.-1 ) return
        
       if(declpri('precip_step1', 1, 'integer', step1) .ne. 0) return

       if(declvar('precip', 'basin_ppt', 'one', 1, 'real',
     +     'Area-weighted adjusted average precip for basin',
     +     'inches',
     +   basin_ppt).ne.0) return

       if(declvar('precip', 'basin_dep', 'one', 1, 'real',
     +     'Area-weighted adjusted average precip+irrig for basin',
     +     'inches',
     +   basin_dep).ne.0) return

      ALLOCATE (mru_ppt(nmru))
      if(declvar('precip', 'mru_ppt', 'nmru', nmru, 'real',
     +     'Adjusted precip on each MRU',
     +     'inches',
     +   mru_ppt).ne.0) return

      ALLOCATE (mru_dep(nmru))
      if(declvar('precip', 'mru_dep', 'nmru', nmru, 'real',
     +     'Adjusted precip+irrig on each MRU',
     +     'inches',
     +   mru_dep).ne.0) return

      ALLOCATE (mru_rain(nmru))
      if(declvar('precip', 'mru_rain', 'nmru', nmru, 'real',
     +     'Computed rain on each MRU',
     +     'inches',
     +   mru_rain).ne.0) return

      ALLOCATE (mru_snow(nmru))
      if(declvar('precip', 'mru_snow', 'nmru', nmru, 'real',
     +     'Computed snow on each MRU',
     +     'inches',
     +   mru_snow).ne.0) return

      ALLOCATE (prmx(nmru))
      if(declvar('precip', 'prmx', 'nmru', nmru, 'real',
     +     'Proportion of rain in a mixed event',
     +     'decimal percent; equal to -0.1 on days with no rain',
     +   prmx).ne.0) return

      ALLOCATE (pptmix(nmru))
      if(declvar('precip', 'pptmix', 'nmru', nmru, 'integer',
     +     'Precip mixture - 0=no, 1=yes',
     +     'none',
     +   pptmix).ne.0) return

      ALLOCATE (newsnow(nmru))
      if(declvar('precip', 'newsnow', 'nmru', nmru, 'integer',
     +     'New snow on MRU, 0=no, 1=yes',
     +     'none',
     +   newsnow).ne.0) return

      ALLOCATE (tmax_allrain_c(nmonths))
      if(declparam('precip', 'tmax_allrain_c', 'nmonths', 'real',
     +   '5.', '0.', '15.',
     +   'Precip all rain if tmax_c above this value',
     +   'If MRU maximum temperature exceeds this value, '//
     +   'precipitation is assumed to be rain.','degrees celsius')
     +   .ne.0) return

      if(declparam('precip', 'tmax_allsnow_c', 'one', 'real',
     +   '0', '-10.', '10.',
     +   'All snow if tmax_c< this value; all rain if tmin_c>this '//
     $   'value.',' If MRU maximum temperature is below this value, '//
     +   'precipitation is assumed to be snow; alternately, if MRU '//
     +   'minimum temperature is above this value, precipitation is '//
     +   'assumed to be all rain.','degrees celsius').ne.0) return

      ALLOCATE (mru_psta(nmru))
      if(declparam('precip', 'mru_psta', 'nmru', 'integer',
     +   '1', 'bounded', 'nrain',
     +   'Index of precipitation station for MRU',
     +   'Index of precipitation station used to compute '//
     +   'rain and snow on MRU.',
     +   'none').ne.0) return

      ALLOCATE (adjmix_rain(nmonths))
      if(declparam('precip', 'adjmix_rain', 'nmonths', 'real',
     +   '1.', '0.', '3.',
     +   'Adjustment factor for rain in a rain/snow mix',
     +   'Monthly factor to adjust rain proportion in a mixed '//
     +   'rain/snow event',
     +   'none').ne.0) return

      ALLOCATE (rain_adj(nmru,nmonths))
      if(declparam('precip', 'rain_adj', 'nmru,nmonths', 'real',
     +   '1.0', '0.2', '5.0',
     +   'Rain adjustment factor, by month for each MRU',
     +   'Monthly factor to adjust measured precipitation to '//
     +   'each mru to account for differences in elevation, etc',
     +   'none')
     +   .ne.0) return

      ALLOCATE (snow_adj(nmru,nmonths))
      if(declparam('precip', 'snow_adj', 'nmru,nmonths', 'real',
     +   '1.0', '0.2', '5.0',
     +   'Snow adjustment factor, by month for each mru',
     +   'Monthly factor to adjust measured precipitation to '//
     +   'each mru to account for differences in elevation, etc',
     +   'none')
     +   .ne.0) return

      ALLOCATE (strain_adj(nmru,nmonths))
      if(declparam('precip', 'strain_adj', 'nmru,nmonths', 'real',
     +   '1.0', '0.2', '5.0',
     +   'Storm rain adjustment factor, by month for each mru',
     +   'Monthly factor to adjust measured precipitation to '//
     +   'each mru to account for differences in elevation, etc.'//
     +   'This factor is for the rain gage used for kinematic or'//
     +   'storm routing',
     +   'none').ne.0) return

      ALLOCATE (mru_area_frac(nmru))
      if(declparam('topc', 'mru_area_frac', 'nmru', 'real',
     +   '1', '0', '1',
     +   'Subcatchment area/total area',
     +   'Subcatchment area/total area',
     +   'none').ne.0) return
     
      ALLOCATE (precip(nrain))
      ALLOCATE (tmax_c(nmru))
      ALLOCATE (tmin_c(nmru))

      pptdecl = 0

      return
      end

c***********************************************************************
c
c     pptinit - Initialize precip module - get parameter values
c

      integer function pptinit()

      USE WEBMOD_PRECIP

      pptinit = 1

      step1 = 1

      if(getparam('precip', 'tmax_allrain_c', nmonths, 'real', 
     +   tmax_allrain_c).ne.0) return

      if(getparam('precip', 'tmax_allsnow_c', 1, 'real', tmax_allsnow_c)
     +   .ne.0) return

      if(getparam('precip', 'mru_psta', nmru, 'integer', mru_psta)
     +   .ne.0) return

      if(getparam('precip', 'adjmix_rain', nmonths, 'real', adjmix_rain)
     +   .ne.0) return

      if(getparam('precip', 'rain_adj', nmru*nmonths, 'real', rain_adj)
     +   .ne.0) return

      if(getparam('precip', 'snow_adj', nmru*nmonths, 'real', snow_adj)
     +   .ne.0) return

      if(getparam('precip', 'strain_adj', nmru*nmonths, 'real',
     +   strain_adj).ne.0) return

      if(getparam('basin', 'mru_area_frac', nmru, 'real',
     $     mru_area_frac) .ne.0) return

      form_data(1) = 0

      pptinit = 0

      return
      end

c***********************************************************************
c
c     pptrun - Computes precipitation form (rain, snow or mix) and
c               depth for each MRU, and basin weighted avg. precip
c

      integer function pptrun()

      USE WEBMOD_PRECIP
      
      integer i, ip, jday, mo
      integer nowtime(6)
      real ppt, pcor
      double precision dt
      integer  nstep

      pptrun = 1

      nstep = getstep()
      dt = deltim()
     
      if(getvar('obs', 'precip', nrain, 'real', precip)
     +   .ne.0) return

      if(nform.eq.1) then
        if(getvar('obs', 'form_data', 1, 'integer', form_data)
     +   .ne.0) return
      endif
      
      if(getvar('obs', 'route_on', 1, 'integer', route_on)
     +   .ne.0) return

      if(getvar('temp', 'tmax_c', nmru, 'real', tmax_c)
     +     .ne.0) return

      if(getvar('temp', 'tmin_c', nmru, 'real', tmin_c)
     +     .ne.0) return

      call dattim('now', nowtime)
      jday = julian('now', 'calendar')
      mo = nowtime(2)
      basin_dep = 0.
      basin_ppt = 0.

      call dpreal('precip - precip', precip, 1, 2)
c
c Initialize precip

      do 1000 i = 1,nmru
        newsnow(i) = 0
        ip = mru_psta(i)
        ppt = precip(ip)
        prmx(i) = -0.1

        mru_ppt(i) = precip(ip)

        if(ppt.eq.0) then    ! No deposition
          pptmix(i) = 0
          mru_dep(i) = 0.
          mru_rain(i) = 0.
          mru_snow(i) = 0.
          go to 1000          ! so return
       end if

c
c Since there is deposition, compute corrections to quantity and
c the ratios of rain to snow.
       
C******If within storm period for kinematic routing, adjust precip
C******by storm adjustment factor. Storm routing assumes all rain since
c      only summer storms are simulated to avoid complications with snowpack
c      processes.

        if(route_on.eq.1) then
          pcor = strain_adj(i,mo)
          pptmix(i) = 0
          prmx(i) = 1
          mru_ppt(i) = ppt * pcor
          mru_rain(i) = mru_ppt(i)
          mru_snow(i) = 0.

C******If observed temperature data are not available or if observed
C******form data are available and rain is explicitly specified then
C******precipitation is all rain.

        else if(form_data(1).eq.2) then
          pcor = rain_adj(i,mo)
          pptmix(i) = 0
          prmx(i)=1
          mru_ppt(i) = ppt * pcor
          mru_rain(i) = mru_ppt(i)
          mru_snow(i) = 0.

C******If form data are available and snow is explicitly specified or if
C******maximum temperature is below or equal to the base temperature for
C******snow then precipitation is all snow

        else if(form_data(1).eq.1.or.tmax_c(i).le.tmax_allsnow_c) then
          pcor = snow_adj(i,mo)
          pptmix(i) = 0
          mru_ppt(i) = ppt * pcor
          mru_rain(i) = 0
          mru_snow(i) = mru_ppt(i)
          if(mru_snow(i).gt.0) prmx(i) = 0
          newsnow(i) = 1

C******If minimum temperature is above base temperature for snow or
C******maximum temperature is above all_rain temperature then
C******precipitation is all rain

        else if(tmin_c(i).gt.tmax_allsnow_c.or.
     +          tmax_c(i).ge.tmax_allrain_c(mo)) then
          pcor = rain_adj(i,mo)
          pptmix(i) = 0
          mru_ppt(i) = ppt * pcor
          mru_rain(i) = mru_ppt(i)
          if(mru_rain(i).gt.0) prmx(i) = 1
          mru_snow(i) = 0.

C******Otherwise precipitation is a mixture of rain and snow

        else
          prmx(i) = ((tmax_c(i)-tmax_allsnow_c)/(tmax_c(i)-tmin_c(i)))*
     +             adjmix_rain(mo)

C******Unless mixture adjustment raises the proportion of rain to
C******greater than or equal to 1.0 in which case it all rain

          if(prmx(i).gt.1.) then
            pcor = rain_adj(i,mo)
            pptmix(i) = 0
            prmx(i) = 1.0
            mru_ppt(i) = ppt * pcor 
            mru_rain(i) = mru_ppt(i)
            mru_snow(i) = 0.

C******If not, it is a rain/snow mixture

c
          else
            pcor = snow_adj(i,mo)
            pptmix(i) = 1
c
c           RMTW - Changed the following logic so that snow_adj only applies
c                  to the snow portion. Otherwise rain amounts have snow_adj
c                  applied even if there is only a trace of snow. Update prmx.
c
c            mru_ppt(i) = ppt * pcor
c            mru_rain(i) = prmx(i) * mru_ppt(i)
c            mru_snow(i) = mru_ppt(i) - mru_rain(i)

            mru_snow(i) = (1-prmx(i))*ppt * pcor
            pcor = rain_adj(i,mo)
            mru_rain(i) = prmx(i)* ppt * pcor
            mru_ppt(i) = mru_rain(i) + mru_snow(i)
            prmx(i) = mru_rain(i)/mru_ppt(i)
c
            newsnow(i) = 1
          end if
       end if

       mru_dep(i) = mru_ppt(i)

       basin_ppt = basin_ppt + mru_ppt(i) * mru_area_frac(i)
       basin_dep = basin_dep + mru_dep(i) * mru_area_frac(i)

 1000 continue

      pptrun = 0

      return
      end

