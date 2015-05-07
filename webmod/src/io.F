#include "defines.h"
c***********************************************************************
c     Open and close the file for the objective functions.
c     Modified for xtopprms by rmwebb from obj_func_xml.f
c     version: 1.1.3 (markstro)
c
c     5/01/03 - added time step integer, endper, to flag
c       the end of storms, days, months, years, and run - RMTW
c
c     9/01/03 - Added RCS version control - RMTW
c
c     9/17/03 - Add file unit for outputting solute chemistry
c       from the webmod_res and phreeq_mms modules
c
c    27 apr  2010 - Port to Fortran 90 with module and dynamic memory
c
c
c    17 jan 2014 - Added files for volumes, solutes, and entities. LUNs for 
c     topout_file_unit, chemout_file_unit, and all out output files
c     now assigned using NEWUNIT.
c
c***********************************************************************
      MODULE WEBMOD_IO
      IMPLICIT NONE
#include "fmodules.inc"
!   Dimensions and Local Variables
!      integer, SAVE:: topout_file_unit, chemout_file_unit, phreeqout ! standard output, in addition to MMS files *.out, *.statvar, etc
      integer, save:: endper, yrdays, modays(12),nowtime(6)
      integer, save:: xdebug_start, xdebug_stop
      integer, save, allocatable :: vse_lun(:) ! LUNs for volume, solute, and entitity files
      integer, save:: print_vse, nf  ! In addition to the four standard files, additional detailed files
                                      ! may be output: nf files include volume files output in webmod_res and 
                                      ! chemistry files (solutes and entities) output in phreeq_mms
      data nf/0/ ! no extra files if print_vse = 0. This number increments for each detailed file opened in webmod_res and phreeq_mms
      TYPE :: outfiles   ! file names, shortnames, and logical unit numbers for input and output files.
         character(60) :: file   ! Output file
         integer       :: lun        ! integer returned by NEWUNIT
      END TYPE outfiles
!
      TYPE(outfiles), save :: topout, chemout, phreeqout,debug
c     ! 
c     ! TYPE(outfiles),save,allocatable :: c_mru(:), c_uzgen(:), These
c     !$ c_uzrip(:), c_uzup(:), c_can(:), c_snow(:), c_inperv(:),
c     !$ c_transp(:), c_ohoriz(:), c_uz(:,:), c_qdf(:), c_sat(:),
c     !$ c_satpref(:), c_hill(:), c_uz2sat(:), c_hyd(:)
c     ! 
c     ! integer, save, allocatable :: v_lun(:), c_lun(:)
c     ! integer, save :: nvf, ncf
      logical, save:: step1
      double precision, save:: endjday
      data modays/31,28,31,30,31,30,31,31,30,31,30,31/
      END MODULE WEBMOD_IO
c***********************************************************************
c     Main basin_sum routine c

      integer function io(arg)
      IMPLICIT NONE

! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
      CHARACTER(len=256) SVN_ID
      integer, EXTERNAL :: iodecl, ioinit, ioclean, iorun
      save SVN_ID

      SVN_ID = 
     $     '$Id: io.f 32 2007-06-08 17:24:28Z rmwebb $ '
      io = 0
      if(arg.eq.'declare') then
        io = iodecl()
      else if(arg.eq.'initialize') then
        io = ioinit()
      else if(arg.eq.'run') then
        io = iorun()
      else if(arg.eq.'cleanup') then
        io = ioclean ()
       end if
      return
      end
c***********************************************************************
c
c     iodecl
c

      integer function iodecl()

      USE WEBMOD_IO
      IMPLICIT NONE

      iodecl = 1

!      if(decl*param('io', 'topout_file_unit', 'one', 'integer',
!     +   '80', '50', '99',
!     +   'Unit number for TOPMODEL output file',
!     +   'Unit number for TOPMODEL output file',
!     +   'integer').ne.0) return

!      if(decl*param('io', 'chemout_file_unit', 'one', 'integer',
!     +   '90', '50', '99',
!     +   'Unit number for file summarizing solute transport',
!     +   'Unit number for file summarizing solute transport',
!     +   'integer').ne.0) return

c
c Single end period flag that can be decomposed to logical variables
c (for example, end_run=.true.) in modules where needed - RMTW
c
c Note: Step1 in the init section is not assigned a decl*pri
c since we want it to revert to .true.(from the data statement)
c upon a restart(init file). Data section is only read on first
c run so reset to true in the init section. Can be conditioned
c on savevar condition later.
c

      if(declparam('sumb', 'print_vse', 'one', 'integer', 
     +   '1', '0', '2',
     +   'Scope of geochemical output files (vse)',
     +   'Scope of geochemical output files (vse): 0 = none;'//
     +   '                  1 = basin and MRUs;            '//
     +   '                  2 = all reservoirs.                ',
     +   'none').ne.0)return

      if(declvar('io', 'endper', 'one', 1, 'integer', 
     $     'Composite flag indicating the end of a given period: '//
     $     'Run(16)+Year(8)+Month(4)+Day(2)+storm(1)',
     +     'none',endper).ne.0) return

      if(declparam('io', 'xdebug_start', 'one', 'integer',
     +    '0', '0', '1000000',
     +    'Begin debugging on this step (0 if no debugging)',
     +    'Begin debugging on this step (0 if no debugging)',
     +    'integer').ne.0) return

      if(declparam('io', 'xdebug_stop', 'one', 'integer',
     +    '0', '0', '1000000',
     +    'End debugging on this step',
     +    'End debugging on this step',
     +    'integer').ne.0) return

!
      if(declpri('io_modays', 12, 'integer', modays)
     + .ne. 0) return

      if(declpri('io_yrdays', 1, 'integer', yrdays ).ne. 0) return

      if(declpri('io_endjday', 1, 'double', endjday ).ne. 0) return

      iodecl = 0

      return
      end

c***********************************************************************
c
c     ioinit - 
c

      integer function ioinit()

!#if defined(_WIN32)
!      USE IFPORT
!#endif
      USE WEBMOD_IO
      IMPLICIT NONE
!!      INTEGER, EXTERNAL :: getoutname
      INTEGER, EXTERNAL ::open_res_files
      integer ret, chem_sim, nmru, nac, nhydro
      logical filflg

c End-period variables

      ioinit = 1

      step1 = .true.


! Get Dimensions
      nmru = getdim('nmru')
      IF (nmru.EQ.-1) RETURN
      nac = getdim('nac')
      IF (nac.EQ.-1) RETURN
      nhydro = getdim('nhydro')
      IF (nhydro.EQ.-1) RETURN


      if(getparam('sumb', 'print_vse', 1, 'integer', print_vse)
     +   .ne.0) return

      if(getparam('phreeqmms', 'chem_sim', 1,
     +     'integer',chem_sim) .ne.0) return
     
      if(getparam('io', 'xdebug_start', 1,
     +     'integer',xdebug_start) .ne.0) return

      if(getparam('io', 'xdebug_stop', 1,
     +     'integer',xdebug_stop) .ne.0) return

c
c Open topmodel output file
c
      ret = getoutname (topout%file, '.topout')
c
c Kludge to allow running on 2nd CPU
c
c      output_path='.\output\loch3.topout'
      inquire(file=topout%file,exist=filflg)
      if (filflg) then
        open(newunit=topout%lun,file=topout%file,status='old')
        close(topout%lun,status='delete')
      endif

c-----open the file.
      open (newunit=topout%lun,file=topout%file,
     +    access='sequential', form='formatted', status='new')

c
c Open solute model output file
c
      ret = getoutname (chemout%file, '.chemout')
c Kludge for running on 2nd CPU
c      output_path='.\output\loch3.chemout'
c End Kludge
      inquire(file=chemout%file,exist=filflg)
      if (filflg) then
        open(newunit=chemout%lun,file=chemout%file,status='old')
        close(unit=chemout%lun,status='delete')
      endif

c-----open the file.
      open (newunit=chemout%lun,file=chemout%file,
     +    access='sequential', form='formatted', status='new')
c
c Convert the end of run time to an absolute julian date
c
      endjday = djulian('end', 'absolute')
!
! select_mixes file
!
!      phreeqout%file='./output/select_mixes'
!      inquire(file=phreeqout%file,exist=filflg)
!      if (filflg) then
!        open(newunit=phreeqout%lun,file=phreeqout%file,status='old')
!        close(unit=phreeqout%lun,status='delete')
!      endif
!!!----open the file.
!      open (newunit=phreeqout%lun,file=phreeqout%file,access=
!     * 'sequential',form='formatted', status='new')
!      close (unit = phreeqout%lun)  ! Output for selected mixes. Empty at first then appended with each mix while debug is on
!
! Debug file
!
      !Debug%file='./output/Debug'
      Debug%file='./Debug'
      inquire(file=Debug%file,exist=filflg)
      if (filflg) then
        open(newunit=debug%lun,file=debug%file,status='old')
        close(unit=debug%lun,status='delete')
      endif
!----open the file.
      open (newunit=debug%lun,file=debug%file,access='sequential',
     * form='formatted', status='new')
!
! Detailed output in addition to the four above.
! When print_vse equals 1 or 2, calculate the number of volume files written in webmod_res and
! the number of solute and entity files written in phr3eeq_mms. 'nf' is the the total number of
! files to be open, written, and  
! when print_vse = 1 nf includes a summary file for the basin and each mru
! when print_vse = 2 nf includes addition files for each hillslope and stream reservoir
! the initial prefixes are for volume (v_), solutes (s_), and entities (e_).
! If chemistry in not being simulated (chem_sim=0), no solute or entity files will be written.
!
      if(print_vse.eq.1) then
          nf = nmru+1
          if(chem_sim.eq.1) then
              nf = 3*(nmru+1)
          endif
      elseif(print_vse.eq.2) then
          nf=1+nmru*(13+nac)+nhydro+1 ! first is basin, last is stream volumes, others are hillslope reservoirs
          if(chem_sim.eq.1) then
              nf =3*(1+nmru*(13+nac)+nhydro+1)
          endif
      endif
      if(nf.gt.0) allocate(vse_lun(nf))
      nf=0 ! reset to allow incremental assignments in webmod_res and phreeq_mms
      ioinit = 0

      return
      end

c***********************************************************************
c
c     iorun - File handling takes place in the init and cleanup 
c             sections of this module. The run section will establish
c             the flags indicating the end of periods.
c

      integer function iorun()

      USE WEBMOD_IO
      IMPLICIT NONE

      logical end_run, end_yr, end_mo, end_dy, end_storm
      integer iend_run, iend_yr, iend_mo, iend_dy, iend_storm
      integer year, mo, day, wyday, cday
      double precision jday, stepcheck, timestep
 
      iorun = 1

      end_run = .false.
      end_yr = .false.
      end_mo = .false.
      end_dy = .false.
      end_storm = .false.

      iend_run = 0
      iend_yr = 0
      iend_mo = 0
      iend_dy = 0
      iend_storm = 0

      call dattim('now', nowtime)
      wyday = julian('now', 'water')
      cday = julian('now','calendar')
      jday = djulian('now', 'absolute')
      year = nowtime(1)
      mo = nowtime(2)
      day = nowtime(3)
c
c Check for leap years on the the first iteration and on Jan 1
c
      if(step1.or.cday.eq.1) then
         if(isleap(year).eq.1) then
            yrdays = 366
            modays(2) = 29
         else
            yrdays = 365
            modays(2) = 28
         end if
         step1=.false.
      end if

      timestep = deltim()
      stepcheck = (1.1 * timestep/24.)+jday
c
c Check if its the end of the day. The condition tested is whether
c another 1.1 time steps would take you into another day
c
      if (timestep.ge.24.) then
         end_dy = .true.
         iend_dy = 1
      else if (dint(stepcheck).gt.dint(jday)) then
         end_dy = .true.
         iend_dy = 1
      end if

c
c Set all summary flags at end of run so that summaries will be
c printed even if run ends in the middle of the day and/or month.
c
c The end of year is flagged with the last data point
c on September 30, the end of the water year.
c      
      if(jday.ge.endjday) then
         iend_run = 1
         iend_yr = 1
         iend_mo = 1
         iend_dy = 1
         iend_storm = 1
         end_run = .true.
         end_yr = .true.
         end_mo = .true.
         end_dy = .true.
         end_storm = .true.
      else if(wyday.eq.yrdays.and.end_dy.eqv..true.) then
         iend_yr = 1
         iend_mo = 1
         end_yr = .true.
         end_mo = .true.
      else if(day.eq.modays(mo).and.end_dy.eqv..true.) then
         end_mo = .true.
         iend_mo = 1
      end if
c
c Compose the endper integer
c Run(16)+Year(8)+Month(4)+Day(2)+storm(1)
c
      endper = iend_run*16+iend_yr*8+iend_mo*4+iend_dy*2+iend_storm

      iorun = 0

      return
      end

c***********************************************************************
c
c     ioclean

      integer function ioclean ()
      USE WEBMOD_IO
      implicit none
      integer i

      ioclean = 1

      !close (unit = topout_file_unit)
      !
      !close (unit = chemout_file_unit)
      !
      !close (unit = 25)  ! Output for selected mixes
      !
      !close (unit = 26)  ! Debug file
      close (unit = topout%lun)
      
      close (unit = chemout%lun)
      
!      close (unit = phreeqout%lun)  ! Output for selected mixes
      
      close (unit = debug%lun)  ! Debug file
      
      if(nf.ne.0) then
          do i = 1, nf
            close (unit = vse_lun(i))  ! Close volume, solute, and entity files.
          end do
      endif
      ioclean = 0
      return
      end

c$$$      integer function WhatTime(IsIt)  ! used for debugging
c$$$!#if defined(_WIN32)
c$$$!        USE IFPORT
c$$$!#endif
c$$$      integer, intent(OUT) :: IsIt(3)
c$$$      call itime(IsIt)        ! now(1)=hour, (2)=minute, (3)=second
c$$$      return
c$$$      END FUNCTION WhatTime

