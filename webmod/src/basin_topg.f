c***********************************************************************
c     basin_topg.f: Declares basin and HRU physical parameters and 
c              reservoirs 
c
c     4 sept 2003 - added RCS version control. Removed the temperature
c              unit (temp_unit) parameter declaration since it is not
c              an intrinsic property of the basin. - RMTW
c     2 feb  2004 - Change all hru references to mru - RMTW
c    27 apr  2010 - Port to Fortran 90 with module and dynamic memory
c***********************************************************************

      MODULE WEBMOD_BASIN
      IMPLICIT NONE
      include 'fmodules.inc'
!   Dimensions and Local Variables
      INTEGER, SAVE :: Nmru
!   Declared Variables
      REAL, SAVE, ALLOCATABLE :: mru_perv(:), mru_imperv(:)
!   Declared Parameters
      REAL, SAVE :: basin_area
      REAL, SAVE, ALLOCATABLE :: mru_elev(:),mru_area(:),mru_slope(:)
      REAL, SAVE, ALLOCATABLE :: mru_percent_imperv(:)
      END MODULE WEBMOD_BASIN
c***********************************************************************
c
c     Main basin routine
c

      integer function basin_topg(arg)
      implicit none

! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
      CHARACTER*256 SVN_ID
! Functions
      integer, external:: basdecl, basinit
      save SVN_ID

      SVN_ID = 
     $     '$Id: basin_topg.f 31 2007-06-08 17:22:15Z rmwebb $ '

      basin_topg = 0

      if(arg.eq.'declare') then
        basin_topg = basdecl()
      else if(arg.eq.'initialize') then
        basin_topg = basinit()
      end if

C******Debug level print
      call dpint4('End of basin, retval = ', basin_topg, 1, 2)

      return
      end function basin_topg

c***********************************************************************
c 
c     basdecl - set up parameters
c

      integer function basdecl()
      USE WEBMOD_BASIN
      implicit none

      basdecl = 1

! Get Dimensions
      Nmru = getdim('nmru')
      IF ( Nmru.EQ.-1 ) RETURN

! Declared variables

      ALLOCATE (mru_perv(Nmru))
      if(declvar('basin', 'mru_perv', 'nmru', Nmru, 'real',
     +     'MRU pervious area',
     +     'km2',
     +  mru_perv) .ne.0) return

      ALLOCATE (mru_imperv(Nmru))
      if(declvar('basin', 'mru_imperv', 'nmru', Nmru, 'real',
     +     'MRU impervious area',
     +     'km2',
     +  mru_imperv) .ne.0) return

! Declared Parameters

      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return
      
      ALLOCATE (mru_elev(Nmru))
      if(declparam('basin', 'mru_elev', 'nmru', 'real',
     +   '0.', '-300.', '10000',
     +   'Mean elevation for each MRU',
     +   'Mean elevation for each MRU',
     +   'meters').ne.0) return
     
      ALLOCATE (mru_area(Nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      ALLOCATE (mru_slope(Nmru))
      if(declparam('basin', 'mru_slope', 'nmru', 'real',
     +   '0.0', '0.0', '10.0',
     +   'MRU slope in decimal vertical feet/horizontal feet',
     +   'MRU slope in decimal vertical feet/horizontal feet',
     +   'decimal percent')
     +   .ne.0) return

      ALLOCATE (mru_percent_imperv(Nmru))
      if(declparam('basin', 'mru_percent_imperv', 'nmru', 'real',
     +   '0.', '.00', '1+e09',
     +   'MRU impervious area in decimal percent',
     +   'MRU impervious area as a decimal percent of the total '//
     +   'MRU area',
     +   'decimal percent').ne.0) return

      basdecl = 0

      return
      end


c**********************************************************************
c     basinit - check for validity of basin parameters
c               and compute reservoir areas
c


      integer function basinit()

      USE WEBMOD_BASIN
      integer i
      real totarea, diff

      basinit = 1

      nmru = getdim('nmru')

      if(nmru.eq.-1) return

      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', Nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_percent_imperv', Nmru, 'real',
     +    mru_percent_imperv) .ne.0) return

      totarea = 0.
      do 100 i=1,nmru
        mru_imperv(i) = mru_percent_imperv(i) * mru_area(i)
        mru_perv(i) = mru_area(i) - mru_imperv(i)
        totarea = totarea + mru_area(i)
  100 continue
      diff = (totarea - basin_area)/basin_area
      if(ABS(diff).ge..01) then
        call dpstr('Sum of mru areas is not equal to basin area',0)
        return
      end if

      basinit = 0

      return
      end


