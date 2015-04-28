#include "defines.h"
c***********************************************************************
c          irrig.f: Apply irrigation to mru's.
c***********************************************************************
c
c  
c RMTW - 3-9-05  - Pruned from precip_irrig module. Initial check-in
c                  Uses putvar to update mru_dep values so this module
c                  is not standalone and must be accompanied by a
c                  precipitation module (lapse or XYZ).
c
c        7-6-06  - Made mru_dep = rain + snow for all cases.
c
c        06may10 - Port to Fortran 90 with module and dynamic memory
c
c***********************************************************************
c 
      MODULE WEBMOD_IRRIG
      IMPLICIT NONE
      include 'fmodules.inc'

C   Dimensions
      integer, save :: nmru, nirrig_ext, nirrig_int
      integer, save :: nhydro, nmonths, nform
      
C   Declared Variables
      real, save :: basin_irr_ext,basin_irr_sat, basin_irr_hyd
      real, save, allocatable :: irrig_frac_ext(:), irrig_frac_sat(:)
      real, save, allocatable :: irrig_ext_mru(:), irrig_sat_mru(:)
      real, save, allocatable :: irrig_hyd_mru(:), irrig_frac_hyd(:)
      double precision, save, allocatable :: irrig_hyd_seg(:)
C   Declared Parameters
      real, save :: tmax_allsnow_c
      integer,save,allocatable::irrig_sched_ext(:),irrig_sched_int(:)
      integer, save, allocatable :: irrig_int_src(:)
      real, save, allocatable :: tmax_allrain_c(:), adjmix_rain(:)
      real, save, allocatable :: mru_area(:), mru_area_frac(:)
      real, save, allocatable :: irrig_int_init(:)

C   Declared Private Variables
      integer, save :: step1
      
      logical, save :: irrig_warn
      data irrig_warn/.false./

C   Undeclared Static Variables gotten from from other modules
      integer, save :: route_on, clark_segs
      real, save, allocatable :: irrig_ext(:), irrig_sat_next(:)
      real, save, allocatable :: mru_dep(:), mru_ppt(:), mru_rain(:)
      real, save, allocatable :: mru_snow(:),newsnow(:)
      integer, save, allocatable ::  prmx(:), pptmix(:)
      integer, save ::  form_data(1)
      double precision, save, allocatable :: irrig_hyd_next(:)
      double precision, save, allocatable :: irrig_hyd_seg_next(:)
      real, save, allocatable :: tmax_c(:), tmin_c(:)

      REAL, PARAMETER :: a_million = 1e6 , inch2m = 0.0254     
      
      END MODULE WEBMOD_IRRIG
c***********************************************************************
c 
c     irrdecl - set up parameters for precipitation computations
c

      integer function irrdecl()

      USE WEBMOD_IRRIG
      IMPLICIT NONE
      INTEGER, EXTERNAL :: declpri

      irrdecl = 1

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nhydro = getdim('nhydro')
        if ( nhydro.eq.-1 ) return
      nform = getdim('nform')
        if ( nform.eq.-1 ) return
      nmonths = getdim('nmonths')
        if ( nmonths.eq.-1 ) return
      nirrig_ext = getdim('nirrig_ext')
        if ( nirrig_ext.eq.-1 ) return
      nirrig_int = getdim('nirrig_int')
        if ( nirrig_int.eq.-1 ) return

       if(declpri('irrig_step1', 1, 'integer', step1) .ne. 0) return

       if(declvar('irrig', 'basin_irr_ext', 'one', 1, 'real',
     +     'Area weighted adjusted irrigation from external sources',
     +     'inches',
     +   basin_irr_ext).ne.0) return

       if(declvar('irrig', 'basin_irr_sat', 'one', 1, 'real',
     +     'Area weighted adjusted irrigation from shallow wells',
     +     'inches',
     +   basin_irr_sat).ne.0) return

       if(declvar('irrig', 'basin_irr_hyd', 'one', 1, 'real',
     +     'Area weighted adjusted irrigation from stream diversions',
     +     'inches',
     +   basin_irr_hyd).ne.0) return

      ALLOCATE (irrig_ext_mru(Nmru))
      if(declvar('irrig', 'irrig_ext_mru', 'nmru', nmru, 'real',
     +     'Irrigation applied to MRU from external source: ',
     $     'inches', irrig_ext_mru).ne.0) return

      ALLOCATE (irrig_sat_mru(Nmru))
      if(declvar('irrig', 'irrig_sat_mru', 'nmru', nmru, 'real',
     +     'Irrigation applied to MRU from a well.',
     $     'inches', irrig_sat_mru).ne.0) return

      ALLOCATE (irrig_hyd_mru(Nmru))
      if(declvar('irrig', 'irrig_hyd_mru', 'nmru', nmru, 'real',
     +     'Irrigation applied to MRU from a stream.',
     $     'inches', irrig_hyd_mru).ne.0) return

      ALLOCATE (irrig_hyd_seg(Nhydro))
      if(declvar('irrig', 'irrig_hyd_seg', 'nhydro', nhydro,
     $     'double','Volume of water diverted from a stream '//
     $     'segment to be applied to one or more MRUs.',
     $     'm3', irrig_hyd_seg).ne.0) return

      ALLOCATE (irrig_frac_ext(Nmru))
      if(declvar('irrig', 'irrig_frac_ext', 'nmru', nmru, 'real',
     +     'Fraction of deposition that was irrigation from an '//
     $     'external source; equal to -0.1 on days '//
     $     'with no rain or irrig','decimal percent',
     +   irrig_frac_ext).ne.0) return

      ALLOCATE (irrig_frac_sat(Nmru))
      if(declvar('irrig', 'irrig_frac_sat', 'nmru', nmru, 'real',
     +     'Fraction of deposition that was irrigation from a well; '//
     $     'equal to -0.1 on days with no rain or irrig',
     $     'decimal percent',
     +   irrig_frac_sat).ne.0) return

      ALLOCATE (irrig_frac_hyd(Nmru))
      if(declvar('irrig', 'irrig_frac_hyd', 'nmru', nmru, 'real',
     +     'Fraction of deposition that was irrigation from a '//
     $     'stream; equal to -0.1 on days '//
     $     'with no rain or irrig','decimal percent',
     +   irrig_frac_hyd).ne.0) return

      ALLOCATE (tmax_allrain_c(Nmonths))
      if(declparam('irrig', 'tmax_allrain_c', 'nmonths', 'real',
     +   '5.', '0.', '15.',
     +   'Precip all rain if tmax_c above this value',
     +   'If MRU maximum temperature exceeds this value, '//
     +   'precipitation is assumed to be rain.','degrees celsius')
     +   .ne.0) return

      if(declparam('irrig', 'tmax_allsnow_c', 'one', 'real',
     +   '0', '-10.', '10.',
     +   'All snow if tmax_c< this value; all rain if tmin_c>this '//
     $   'value.',' If MRU maximum temperature is below this value, '//
     +   'precipitation is assumed to be snow; alternately, if MRU '//
     +   'minimum temperature is above this value, precipitation is '//
     +   'assumed to be all rain.','degrees celsius').ne.0) return

      ALLOCATE (adjmix_rain(nmonths))
      if(declparam('irrig', 'adjmix_rain', 'nmonths', 'real',
     +   '1.', '0.', '3.',
     +   'Adjustment factor for rain in a rain/snow mix',
     +   'Monthly factor to adjust rain proportion in a mixed '//
     +   'rain/snow event',
     +   'none').ne.0) return

      ALLOCATE (irrig_sched_ext(Nmru))
      if(declparam('irrig', 'irrig_sched_ext', 'nmru', 'integer',
     +   '0',  '0', '100',
     +   'Index of external irrigation schedule for MRU; 0 if none',
     +   'Index of external irrigation schedule for MRU; 0 if none',
     +   'none').ne.0) return

      ALLOCATE (irrig_sched_int(Nmru))
      if(declparam('irrig', 'irrig_sched_int', 'nmru', 'integer',
     +   '0', '0', '100',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'none').ne.0) return
c
c If there is an internal irrigation schedule for an MRU (irrig_sched_int>0),
c then irrig_int_src indicates if the MRU has recieves irrigations from a
c well in the MRU (value=0), or a stream segment (value > 0).
c
      ALLOCATE (irrig_int_src(Nmru))
      if(declparam('irrig', 'irrig_int_src', 'nmru', 'integer',
     +   '0', '0', '100',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   'none')
     +   .ne.0) return 

      ALLOCATE (irrig_int_init(Nmru))
      if(declparam('irrig', 'irrig_int_init', 'nmru', 'real',
     + '0', '0', '100',
     + 'Irrigation from an internal source to be applied on first day',
     + 'Irrigation from an internal source to be applied on first day',
     + 'inches').ne.0) return

      ALLOCATE (mru_area(Nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      ALLOCATE (mru_area_frac(Nmru))
      if(declparam('topc', 'mru_area_frac', 'nmru', 'real',
     +   '1', '0', '1',
     +   'Subcatchment area/total area',
     +   'Subcatchment area/total area',
     +   'none').ne.0) return


      ALLOCATE(irrig_ext(nirrig_ext), irrig_sat_next(nmru)) 
      ALLOCATE(mru_dep(nmru), mru_ppt(nmru), mru_rain(nmru))
      ALLOCATE(mru_snow(nmru),newsnow(nmru))
      ALLOCATE(prmx(nmru), pptmix(nmru))
      ALLOCATE(irrig_hyd_next(nmru))
      ALLOCATE(irrig_hyd_seg_next(nhydro))
      ALLOCATE(tmax_c(nmru), tmin_c(nmru))

      irrdecl = 0

      return
      end

c***********************************************************************
c
c     irrinit - Initialize irrigation module - get parameter values
c

      integer function irrinit()

      USE WEBMOD_IRRIG
      IMPLICIT NONE

      irrinit = 1

      step1 = 1


      if(getparam('precip', 'tmax_allrain_c', nmonths, 'real', 
     +   tmax_allrain_c).ne.0) return

      if(getparam('precip', 'tmax_allsnow_c', 1, 'real',
     +   tmax_allsnow_c).ne.0) return

      if(getparam('precip', 'adjmix_rain', nmonths, 'real',
     +   adjmix_rain).ne.0) return

      if(getparam('precip', 'irrig_sched_ext', nmru, 'integer',
     $     irrig_sched_ext).ne.0) return

      if(getparam('precip', 'irrig_sched_int', nmru, 'integer',
     $     irrig_sched_int).ne.0) return

      if(getparam('precip', 'irrig_int_src', nmru, 'integer',
     $     irrig_int_src).ne.0) return

      if(getparam('precip', 'irrig_int_init', nmru, 'integer',
     $     irrig_int_init).ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real',
     $     mru_area) .ne.0) return

      if(getparam('basin', 'mru_area_frac', nmru, 'real',
     $     mru_area_frac) .ne.0) return

      form_data(1) = 0

      irrinit = 0

      return
      end

c***********************************************************************
c
c     irrrun - Computes form (rain, snow or mix) and amount of irrigation
c              from wells or diversions.
c

      integer function irrrun()

      USE WEBMOD_IRRIG
      IMPLICIT NONE
      INTEGER, EXTERNAL :: julian
      INTEGER, EXTERNAL :: putvar

      integer i, jday, mo
      integer nowtime(6)
      real basin_dep, dep
      double precision dt, diff
      integer  nstep,iext, iint, isrc

      irrrun = 1

      nstep = getstep()
      dt = deltim()
c
c Get the atmospheric precip values
c
      if(getvar('precip', 'mru_ppt', nmru, 'real',
     +    mru_ppt).ne.0) return

      if(getvar('precip', 'mru_dep', nmru, 'real',
     +    mru_dep).ne.0) return

      if(getvar('precip', 'mru_rain', nmru, 'real',
     +    mru_rain).ne.0) return

      if(getvar('precip', 'mru_snow', nmru, 'real',
     +    mru_snow).ne.0) return

      if(getvar('precip', 'prmx', nmru, 'real',
     +    prmx).ne.0) return

      if(getvar('precip', 'pptmix', nmru, 'integer',
     +    pptmix).ne.0) return

      if(getvar('precip', 'newsnow', nmru, 'integer',
     +    newsnow).ne.0) return
c
c Get irrigation amounts
c
      if (nirrig_ext.gt.0) then
         if(getvar('obs', 'irrig_ext', nirrig_ext, 'real', irrig_ext)
     +        .ne.0) return
      end if

      if(getvar('topc', 'irrig_sat_next', nmru, 'real',
     +    irrig_sat_next).ne.0) return

      if(getvar('routec', 'irrig_hyd_next', nmru, 'real',
     +    irrig_hyd_next).ne.0) return

      if(getvar('routec', 'clark_segs', 1, 'integer',
     +    clark_segs).ne.0) return

      if(getvar('routec', 'irrig_hyd_seg_next', nhydro, 'double',
     +    irrig_hyd_seg_next).ne.0) return

      if(nform.gt.0) then
        if(getvar('obs', 'form_data', nform, 'integer', form_data)
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

c precip variables have been initialized and/or computed in
c the precip module

      basin_dep = 0.    ! basin_dep is recomputed with irr added
c$$$      basin_ppt = 0. This goes unchanged
      basin_irr_ext = 0.
      basin_irr_sat = 0.
      basin_irr_hyd = 0.

c
c Set diversions from all segments to 0.0
c
      do i = 1, nhydro
         irrig_hyd_seg(i) = 0.0
      end do

c
c Initialize internal irrigation (from wells or streams internal to the watershed)
      if(step1.eq.1) then
         
         do i=1,nmru
            iint = irrig_sched_int(i)
            isrc = irrig_int_src(i)
            if(iint.ge.1) then
               if(isrc.eq.0) then
                  irrig_sat_next(i) = irrig_int_init(i)
                  irrig_hyd_next(i) = 0.0
               else
                  irrig_hyd_next(i) = irrig_int_init(i)
                  irrig_hyd_seg(isrc) = irrig_hyd_seg(isrc) +
     $                 (irrig_hyd_next(i)*inch2m*mru_area(i)*
     $                 a_million)
                  irrig_sat_next(i) = 0.0
               end if
            end if
         end do
         step1 = 0
      end if


      do 1000 i = 1,nmru

        irrig_frac_ext(i) = -0.1
        irrig_frac_sat(i) = -0.1
        irrig_frac_hyd(i) = -0.1
        
        iext = irrig_sched_ext(i)

        if(iext.eq.0) then
           irrig_ext_mru(i) = 0
        else if (iext.gt.0) then
           irrig_ext_mru(i) = irrig_ext(iext)
        else
           print*,'Irrig_sched_ext must be equal to  or '//
     $          'greater than zero. The value for MRU',i,' is ',iext
           return
        end if

        iint = irrig_sched_int(i)
        isrc = irrig_int_src(i)

        if(iint.eq.0) then
           irrig_sat_mru(i) = 0
           irrig_hyd_mru(i) = 0
        else if(iint.gt.0) then
           if(isrc.eq.0) then
              irrig_sat_mru(i) = irrig_sat_next(i)
              irrig_hyd_mru(i) = 0
           else if(isrc.gt.0) then
              irrig_hyd_mru(i) = irrig_hyd_next(i)
              irrig_hyd_seg(isrc) = irrig_hyd_seg(isrc) +
     $             (irrig_hyd_mru(i)*inch2m*mru_area(i)*
     $             a_million)
              irrig_sat_mru(i) = 0
           else
              print*,'Irrig_int_src for MRU ',i ,'must be '//
     $             'greater than or equal to zero. '//
     $             'The value for MRU ',i,' is ',irrig_int_src(i)
           print*,'Run terminated '
              return
           end if
        else
           print*,'Irrig_sched_int must be equal to or '//
     $          'greater than zero. The value for MRU',i,' is ',iint
           print*,'Run terminated '
           return
        end if

        dep = mru_ppt(i) + irrig_ext_mru(i)+irrig_sat_mru(i)+
     $           irrig_hyd_mru(i)

        if(dep.ne.0) then ! If no deposition, jump section
c                      and go to next MRU
c
c Since there is deposition, compute corrections to quantity and
c the ratios of rain to snow. Depths are only adjusted for
c atmospheric precipitation. The depth of irrigation is assumed
c to be correct. The final ratio of rain/snow computed for
c the precip will be that for irrigation on days of no atmospheric
c precip.
c
       
C******If within storm period for kinematic routing, adjust precip
C******by storm adjustment factor. Storm routing assumes all rain since
c      only summer storms are simulated to avoid complications with snowpack
c      processes.

        if(route_on.eq.1) then
          pptmix(i) = 0
          prmx(i) = 1
          mru_rain(i) = dep
          mru_snow(i) = 0.

C******If observed temperature data are not available or if observed
C******form data are available and rain is explicitly specified then
C******precipitation is all rain.

       else if(form_data(1).eq.2) then
          pptmix(i) = 0
          prmx(i)=1
          mru_rain(i) = dep
          mru_snow(i) = 0.

C******If form data are available and snow is explicitly specified or if
C******maximum temperature is below or equal to the base temperature for
C******snow then precipitation is all snow

        else if(form_data(1).eq.1.or.tmax_c(i).le.tmax_allsnow_c) then
          pptmix(i) = 0
          mru_rain(i) = 0
          mru_snow(i) = dep
          if(mru_snow(i).gt.0) prmx(i) = 0
          newsnow(i) = 1

C******If minimum temperature is above base temperature for snow or
C******maximum temperature is above all_rain temperature then
C******precipitation is all rain

        else if(tmin_c(i).gt.tmax_allsnow_c.or.
     +          tmax_c(i).ge.tmax_allrain_c(mo)) then
          pptmix(i) = 0
          mru_rain(i) = dep
          if(mru_rain(i).gt.0) prmx(i) = 1
          mru_snow(i) = 0.

C******Otherwise precipitation is a mixture of rain and snow

        else
          prmx(i) = ((tmax_c(i)-tmax_allsnow_c)/(tmax_c(i)-tmin_c(i)))*
     +             adjmix_rain(mo)
C******Unless mixture adjustment raises the proportion of rain to
C******greater than or equal to 1.0 in which case it all rain

          if(prmx(i).gt.1.) then
            pptmix(i) = 0
            prmx(i) = 1.0
            mru_rain(i) = dep
            mru_snow(i) = 0.

C******If not, it is a rain/snow mixture

c
          else
            pptmix(i) = 1
c
c           RMTW - Changed the following logic so that snow_adj only applies
c                  to the snow portion. Otherwise rain amounts have snow_adj
c                  applied even if there is only a trace of snow. Update prmx.
c
c            mru_ppt(i) = ppt * pcor
c            mru_rain(i) = prmx(i) * mru_ppt(i)
c            mru_snow(i) = mru_ppt(i) - mru_rain(i)

            mru_snow(i) = (1-prmx(i))*dep
            mru_rain(i) = prmx(i)* dep
            prmx(i) = mru_rain(i)/(mru_rain(i)+mru_snow(i))
            newsnow(i) = 1
         end if
       end if

c
c record fractions (may be eliminated later)
c
       
       mru_dep(i) = mru_rain(i) + mru_snow(i)
       irrig_frac_ext(i)=irrig_ext_mru(i)/dep
       irrig_frac_sat(i)=irrig_sat_mru(i)/dep
       irrig_frac_hyd(i)=irrig_hyd_mru(i)/dep

       basin_dep = basin_dep + mru_dep(i) * mru_area_frac(i)
       basin_irr_ext = basin_irr_ext + 
     $      irrig_ext_mru(i) * mru_area_frac(i)
       basin_irr_sat = basin_irr_sat + 
     $      irrig_sat_mru(i)* mru_area_frac(i)
       basin_irr_hyd = basin_irr_hyd + 
     $      irrig_hyd_mru(i)* mru_area_frac(i)
 
      end if        ! Landing pad for no dep in an MRU

 1000 continue
c
c Check that the diversions were equal to those set last time
c step in the route_clark module
c
      if(nstep.gt.1) then
         do i=1,nhydro
            diff = abs(irrig_hyd_seg_next(i)-irrig_hyd_seg(i))
            if(.not.irrig_warn.and.diff.gt.0.001) then
               print*,'Diversions set in route_clark do not equal '//
     $              'those deposited in the irrigation module.'
               print*,' This warning will not be repeated.'
               irrig_warn = .true.
! No need for return, this is just an alert to the user that irrigation
! applied was less than that requested by irrig_hyd_seg_next
!               return
            end if
         end do
      end if
c
c Use this for updating the deposition variables. 
c
      if(putvar('precip', 'basin_dep',  1, 'real', basin_dep)
     &   .ne.0) return

      if(putvar('precip', 'mru_dep',  nmru, 'real', mru_dep)
     &   .ne.0) return

      if(putvar('precip', 'mru_rain',  nmru, 'real', mru_rain)
     &   .ne.0) return

      if(putvar('precip', 'mru_snow',  nmru, 'real', mru_snow)
     &   .ne.0) return

      if(putvar('precip', 'prmx',  nmru, 'real', prmx)
     &   .ne.0) return

      if(putvar('precip', 'pptmix',  nmru, 'integer', pptmix)
     &   .ne.0) return

      if(putvar('precip', 'newsnow',  nmru, 'integer', newsnow)
     &   .ne.0) return

      irrrun = 0

      return
      end

c***********************************************************************
c
c     Main irrigation routine
c
      integer function irrig_web(arg)
      IMPLICIT NONE

      integer irrdecl, irrinit, irrrun
      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      save SVN_ID

      SVN_ID = 
     $     '$Id: irrig_web.f 33 2007-06-08 17:26:29Z rmwebb $ '

      irrig_web = 0

      if(arg.eq.'declare') then
        irrig_web = irrdecl()

      else if(arg.eq.'initialize') then
        irrig_web = irrinit()

      else if(arg.eq.'run') then
        irrig_web = irrrun()

      end if
      
      return
      end
