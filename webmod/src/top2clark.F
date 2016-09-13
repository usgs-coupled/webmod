#include "defines.h"
!**********************************************************************
!   top2clark.f
!
!   Module to couple topmodel outputs to stream channel segments
!   for the standard topmodel routing (based on Clark, 1945). This
!   coupling will permit distinct left and right banks to flow into
!   single channel segments.
!
!   Rick Webb 25 Oct 2002
!
!      8 Sept 2003 - Added RCS version control - RMTW
!     
!     References
!
!     Beven, K.J., and Kirkby, M.J., 1979, A physically based
!     variable contributing area model of basin hydrology:
!     Hydrology Science Bulletin 24(1), p. 43-69.
!
!     Clark, C.O., 1945, Storage and the unit hydrograph,
!     American Society of Civil Engineers Transactions, 
!     v.110, p.1419-1488.
!
!   18may09 - Add Fortran90 Module: WEBMOD_CLARK - RMTW
!
!
!***********************************************************************


!***********************************************************************
! ***************** START MODULE *************
      MODULE WEBMOD_CLARK
      IMPLICIT NONE
#include "fmodules.inc"
      

      
!   Parameters and initial flags
      real a_million
      DATA a_million / 1e6 /
!   Dimensions
      integer, save :: nmru, nchan
!   Declared Parameters
!      integer, save :: topout_file_unit
      integer, save, allocatable :: mru2chan(:)
      real, save :: basin_area
      real, save, allocatable :: mru_area(:), mru_area_frac(:)
!   Declared Variables
      real, save, allocatable :: qchanin(:)
      real, save, allocatable :: chan_area_frac(:), chan_area(:)

!   Undeclared Static Variables gotten from from other modules
      real, save, allocatable ::  qout(:)
      
      END MODULE WEBMOD_CLARK

!***********************************************************************


!***********************************************************************
! 
!     top2cdecl - declare variables and parameters
!

      integer function top2cdecl()

      USE WEBMOD_CLARK
      IMPLICIT NONE

      top2cdecl = 1

! Get dimensions

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nchan = getdim('nchan')
        if ( nchan.eq.-1 ) return

      ALLOCATE (qchanin(Nchan))
      if(declvar('top2c', 'qchanin', 'nchan', nchan, 'real', 
     + 'Sum of mru inputs into each channel reach for each'//
     + 'time step.','m^3',qchanin).ne.0) return

      ALLOCATE (chan_area_frac(Nchan))
      if(declvar('top2c', 'chan_area_frac', 'nchan', nchan, 'real', 
     + 'Decimal fraction of total basin area drained by a given '//
     + 'channel segment. Equal to the sum of contibuting mru_area_frac',
     + 'none',chan_area_frac).ne.0) return

      ALLOCATE (chan_area(Nchan))
      if(declvar('top2c', 'chan_area', 'nchan', nchan, 'real', 
     + 'Sum of MRU areas contributing to the channel reach',
     + 'km2',chan_area).ne.0) return

      ALLOCATE (mru2chan(Nmru))
      if(declparam('top2c', 'mru2chan', 'nmru', 'integer',
     +   '1', 'bounded', 'nchan',
     +   'Index of channel receiving discharge from MRU',
     +   'Index of channel receiving discharge from MRU',
     +   'none').ne.0) return

      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

      ALLOCATE (mru_area_frac(Nmru))
      if(declparam('topc', 'mru_area_frac', 'nmru', 'real',
     +   '1', '0', '1',
     +   'Subcatchment area/total area',
     +   'Subcatchment area/total area',
     +   'none').ne.0) return

      ALLOCATE (mru_area(Nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

!      if(decl*param('io', 'topout_file_unit', 'one', 'integer',
!     +   '80', '50', '99',
!     +   'Unit number for TOPMODEL output file',
!     +   'Unit number for TOPMODEL output file',
!     +   'integer').ne.0) return
!
! getvar qout from top_chem
!
      ALLOCATE (qout(nmru))

      top2cdecl = 0

      return
      end

!***********************************************************************
!
!     top2cinit - Initialize coupling module - get parameter values,
!

      integer function top2cinit()

      USE WEBMOD_CLARK
      USE WEBMOD_IO, ONLY : topout
      IMPLICIT NONE
      
      integer i, j
      real area_chk

!***  local variables

!      character*135 output_path
!      logical filflg

      top2cinit = 1

!----- set name for topmod unique output file 
!     Output file now set in mudule io.f
!
!      ret = getoutname (output_path, '.topout')
!      inquire(file=output_path,exist=filflg)
!      if (filflg) then
!        open(unit=80,file=output_path,status='old')
!        close(unit=80,status='delete')
!      endif

!-----open the file.
!      open (unit=80,file=output_path,access='sequential',
!     * form='formatted', status='old')


      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area_frac', nmru, 'real',
     + mru_area_frac) .ne.0) return

      if(getparam('top2c', 'mru2chan', nmru, 'integer', mru2chan)
     +   .ne.0) return

!      if(get*param('io', 'topout_file_unit', 1, 'integer',
!     +   topout_file_unit).ne.0) return

      do 10 i = 1, nchan
         chan_area(i)=0.0
         chan_area_frac(i)=0.0
 10   continue

      area_chk = 0.0
      do 20 i = 1, nmru
         j = mru2chan(i)
         chan_area(j) = chan_area(j) + mru_area(i)
         chan_area_frac(j) = chan_area_frac(j) + mru_area_frac(i)
         area_chk = area_chk + mru_area_frac(i)
 20   continue
      if(abs(area_chk - 1.0).gt.0.0001) then
         print*,'The sum of the mru fractional areas '//
     $        'do not add up to 1.0; Instead they add '//
     $        'to ',area_chk,'; Run terminated'
         return
      end if

      WRITE(topout%lun,604)(i, mru_area(i), mru2chan(i), i=1,nmru)
  604 FORMAT(//'HILLSLOPE TO CHANNEL ASSIGNMENTS'/
     $     'Hillslope(MRU)  Area(km2)  Channel'/
     $     '==============  =========  ======='/
     $     (3X,I3,11X,f7.2,6X, I3))

      WRITE(topout%lun,605)basin_area
 605  FORMAT(/'     Basin Area  ', f7.2,' sq.km')

      top2cinit = 0
      

      return
      end

!***********************************************************************
!
!     top2crun - pass hillslope discharge to channel segments
!

      integer function top2crun()

      USE WEBMOD_CLARK
      IMPLICIT NONE
      
      integer i, j

!***  local variables

!      character*135 output_path
!      logical filflg

!----- set name for topmod unique output file 

      top2crun = 1
!
!     Get the total discharge from each hillslop, in meters.
!      
      if(getvar('topc', 'qout', nmru, 'real', qout)
     +   .ne.0) return

      do 40 i = 1, nchan
         qchanin(i) = 0
 40   continue

      do 50 i = 1, nmru

      j = mru2chan(i)
!
!     Convert meters times area in sq.km to cubic meters
!
      qchanin(j) = qchanin(j) + (qout(i)* mru_area(i) * a_million)

 50   continue

      top2crun = 0

      return
      end

!***********************************************************************
!
!     top2cclean - No clean up routine needed since no I/O used
!

!***********************************************************************
!
!     Main top2clark routine
!
      integer function top2clark(arg)
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      integer top2cdecl,top2cinit, top2crun
      save SVN_ID

      SVN_ID = 
     $     '$Id$ '

      top2clark = 0

      if(arg.eq.'declare') then
        top2clark = top2cdecl()

      else if(arg.eq.'initialize') then
        top2clark = top2cinit()

      else if(arg.eq.'run') then
        top2clark = top2crun()

      end if

      return
      end
