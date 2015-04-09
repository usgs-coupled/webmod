***********************************************************************
c   route_clark.f
c
c   Extracted from topmodg.f program and modified to include
c   two-dimensional parameter nchan by ntopchan. This allows
c   the standard topmodel routing (based on Clark, 1945) to be
c   applied to discharges from distinct left and right banks.
c
c   Receives input from the channel coupling module top2clark.f
c
c   Rick Webb 11 Sept 2002
c     
c   9 Sept 2003 - Added RCS version control - RMTW
c
c     References
c     Beven, K.J., and Kirkby, M.J., 1979, A physically based
c     variable contributing area model of basin hydrology:
c     Hydrology Science Bulletin 24(1), p. 43-69.
c
c     Clark, C.O., 1945, Storage and the unit hydrograph,
c     American Society of Civil Engineers Transactions, 
c     v.110, p.1419-1488.
c
c   25 August 2004 - New irrigation routines
c
c   21 Sept 2004 - Added stream loss in cfs/mi. Taken from channel
c     sections from water remaining after irrigation.
c
c     Also eliminated Reach velocity, RV(nmru), and channel
c     velocity CV(nmru), in favor of a single channel velocity
c     for the entire basin. It was incongruous to have the possibility
c     of different stream velocities assigned to the same conceptual
c     time delay channel segment.
c
c   17apr09 - Add Fortran90 Module: WEBMOD_ROUTE
c
c
c***********************************************************************


c***********************************************************************
c! ***************** START MODULE *************
      MODULE WEBMOD_ROUTE
      IMPLICIT NONE
      INCLUDE 'fmodules.inc'

C   Parameters and initial flags
      real inch2m,gpm2m3ph,a_million, qlin2dep
      data inch2m /0.0254/
      data gpm2m3ph /0.227125/
      data a_million / 1e6 /
      data qlin2dep /6.3343e-8/        ! conversion factors: inch to meter and gallons
c                                        per minute to cubic meters per hour.
c                                        Multiply stream loss  in cfs/mile *
c                                        qlin2dep*dt(hr)*channel segment(m)/basin_area(km2)
c                                        to compute area-normalized loss to deep aquifer,
c                                        in meters
C   Dimensions and index
      integer, save :: nchan, nmru, nobs, nhydro
      integer, save :: clark_segs, ntopchan, nirrig_int
C   Declared Parameters
      integer, save :: qobsta  !, topout_file_unit
      integer, save, allocatable :: irrig_int_src(:), irrig_sched_int(:)
      integer, save, allocatable :: nchan_d(:), mru2chan(:)
      real, save :: basin_area, CHV, Q0, dtinit,chan_loss_rate
      real, save, allocatable :: irrig_int_init(:), irrig_int_max(:)
      real, save, allocatable :: mru_area(:)
      real, save, allocatable :: mru_area_frac(:)
C   Declared Variables
      real, save :: chan_loss_basin, qbasincfs, qbasinm, QTOT
      real, save, allocatable :: Q(:,:),ar_fill(:,:), qchanin_m(:)
      real, save, allocatable :: maxchanloss(:), chan_loss(:)
      real, save, allocatable :: sum_chan_loss(:)
      double precision, save, allocatable ::  irrig_hyd_next(:)
      double precision, save, allocatable ::  irrig_hyd_seg_next(:)
C   local variables
      logical, save :: div_exists
      integer, save :: div_warn
      data div_warn/0/
      logical, save :: demand_warn
      data demand_warn/.false./
      integer, save, allocatable ::  nrsc(:), ndsc(:)
      integer, save ::  NTSTEP, irrsrc, irrsched
      real, save, allocatable :: ach(:,:), D(:,:)
      real, save, allocatable :: AR(:,:),TCH(:),tchmax(:)
      real, save, allocatable :: QSCM(:),qsccfs(:)
      real, save :: SUMQOBS,SSQ, F1, F2, Q0DT
      double precision, save, allocatable ::  qavail(:)
C   Undeclared variables gotten from from other modules
      integer, save :: endper
      real, save, allocatable  :: chan_area(:), chan_area_frac(:)
      real, save, allocatable :: qchanin(:), irrig_int_next(:)
      real, save, allocatable :: qobs(:), irrig_hyd_mru(:)

      END MODULE WEBMOD_ROUTE
c
c***********************************************************************
c
c     Main route_clark routine
c

      integer function route_clark(arg)
      IMPLICIT NONE

      character(len=*) arg
      CHARACTER(len=256) SVN_ID
      integer routecdecl,routecinit, routecrun

      SVN_ID = 
     $     '$Id: route_clark.f 37 2007-06-08 17:57:19Z rmwebb $ '
      
      route_clark = 0

      if(arg.eq.'declare') then
        route_clark = routecdecl()
      else if(arg.eq.'initialize') then
        route_clark = routecinit()
      else if(arg.eq.'run') then
        route_clark = routecrun()
      end if

      return
      end

c***********************************************************************
c 
c     routecdecl - declare variables and parameters for route_clark
c

      integer function routecdecl()


      USE WEBMOD_ROUTE
      USE mmf, only : declparam
      IMPLICIT NONE

      Routecdecl = 1

! Get dimensions

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nchan = getdim('nchan')
        if ( nchan.eq.-1 ) return
      nhydro = getdim('nhydro')
        if ( nhydro.eq.-1 ) return
      ntopchan = getdim('ntopchan')
        if ( ntopchan.eq.-1 ) return
      nirrig_int = getdim('nirrig_int')
        if ( nirrig_int.eq.-1 ) return
      nobs = getdim('nobs')
        if ( nobs.eq.-1 ) return

      ALLOCATE (q(nhydro,nchan))
      if(declvar('routec', 'q','nhydro,nchan', nhydro*nchan,
     $     'real', 'discharge discretized by nhydro and nchan',
     $     'm',q) .ne. 0) return
      Q = 0.0 ! Initialize to zero

      ALLOCATE (ar_fill(nhydro,nchan))
      if(declvar('routec', 'ar_fill','nhydro,nchan', nhydro*nchan,
     $     'real', 'Matrix of time delay ordinates for stream routing',
     $     'fraction',ar_fill) .ne. 0) return 

      if(declvar('routec', 'clark_segs','one', 1,
     $     'integer', 'number of time delay bins for this model '//
     $     'run.','none',clark_segs) .ne. 0) return 

      if(declparam('routec', 'chan_loss_rate', 'one', 'real',
     +   '0', '0', '100',
     +   'Channel loss rate',
     +   'Channel loss rate',
     +   'cfs/mi').ne.0) return

c
c The following channel loss variables are computed using the
c   chan_loss_rate, the channel velocity, chv, the time step, DT,
c   and area-distance metrics ach and d. Since the channels in the
c   last time-delay segments may be of lesser length than those
c   in segments closer to the outlet, the maxchanloss needs to be
c   dimensioned by nhydro.
c
      ALLOCATE (maxchanloss(nhydro))
      if(declvar('routec', 'maxchanloss','nhydro',nhydro, 'real',
     $     'Maximum rate of channel loss to deep aquifer '//
     $     'normalized to basin area.','m',maxchanloss).ne. 0) return 

      ALLOCATE (chan_loss(nhydro))
      if(declvar('routec', 'chan_loss','nhydro',nhydro,
     $   'real', 'Channel loss to deep aquifer normalized '//
     $   'to basin area.','m',chan_loss) .ne. 0) return 

      ALLOCATE (sum_chan_loss(nhydro))
      if(declvar('routec', 'sum_chan_loss','nhydro',nhydro,
     $   'real', 'Sum of channel loss to deep aquifer normalized '//
     $   'to basin area.','m',sum_chan_loss) .ne. 0) return 

      if(declvar('routec', 'chan_loss_basin','one', 1,
     $     'real', 'Average channel loss to deep aquifer '//
     $     'for entire bain','m',chan_loss_basin) .ne. 0) return 

      if(declpri('routec_sumqobs', 1, 'real', SUMQOBS)
     + .ne. 0) return 

      if(declpri('routec_SSQ', 1, 'real', SSQ)
     + .ne. 0) return 

      if(declpri('routec_f1', 1, 'real', F1)
     + .ne. 0) return 

      if(declpri('routec_f2', 1, 'real', F2)
     + .ne. 0) return 

      if(declpri('routec_ntstep', 1, 'integer', NTSTEP)
     + .ne. 0) return 

      ALLOCATE (AR(nhydro,nchan))
      if(declpri('routec_ar', nhydro*nchan, 'real', AR)
     + .ne. 0) return 

      ALLOCATE (nrsc(nchan))
      if(declpri('routec_nrsc', nchan, 'integer', nrsc) 
     + .ne. 0) return 

      ALLOCATE (ndsc(nchan))
      if(declpri('routec_ndsc', nchan, 'integer', ndsc)
     + .ne. 0) return 

      ALLOCATE (qchanin_m(nchan))
      if(declvar('routec', 'qchanin_m', 'nchan', nchan, 'real', 
     + 'Average depth of discharge received by channel from '//
     + 'contributing MRUs',
     +  'meters',qchanin_m).ne.0) return
  
      if(declvar('routec', 'qtot', 'one', 1, 'real', 
     + 'Sum of predicted runoff for all subcatchments',
     +  'meters',QTOT).ne.0) return
  
      if(declvar('routec', 'qbasinm', 'one', 1, 'real', 
     + 'Predicted basin runoff, in meters.',
     + 'meters',qbasinm).ne.0) return

      if(declvar('routec', 'qbasincfs', 'one', 1, 'real', 
     + 'Predicted basin runoff, in cfs.',
     + 'cfs',qbasincfs).ne.0) return

      if(declparam('topc', 'qobsta', 'one', 'integer',
     +   '1', 'bounded', 'nobs',
     +   'Index of streamflow station for calculating '//
     +   'objective function.','Index of streamflow station '//
     +   'for calculating objective function.','none').ne.0) return

      ALLOCATE (mru2chan(nmru))
      if(declparam('top2c', 'mru2chan', 'nmru', 'integer',
     +   '1', 'bounded', 'nchan',
     +   'Index of channel receiving discharge from MRU',
     +   'Index of channel receiving discharge from MRU',
     +   'none').ne.0) return

      if(declparam('routec', 'chv', 'one', 'real',
     +   '100', '1', '10000',
     +   'Main stream channel routing velocity.',
     +   'Main stream channel routing velocity.',
     +   'm/h').ne.0) return

      if(declparam('routec', 'q0', 'one', 'real',
     +   '1.', '0', '100000.',
     +   'Initial stream discharge at basin outlet.',
     
     +   'Initial stream discharge at basin outlet',
     +   'm^3/sec').ne.0) return

      ALLOCATE (ach(ntopchan,nchan))
      if(declparam('routec', 'ach', 'ntopchan,nchan', 'real',
     +   '1', '0', '1',
     +   'Fraction of MRU area draining to channel between '//
     +   'distance, d, and subcatchment outlet (0->1).',
     +   'Fraction of MRU area draining to channel between '//
     +   'distance, d, and subcatchment outlet (0->1).',
     +   'none').ne.0) return
     
      ALLOCATE (d(ntopchan,nchan))
      if(declparam('routec', 'd', 'ntopchan,nchan', 'real',
     +   '500', '0', '50000',
     +   'Distance from basin outlet.',
     +   'Distance from basin outlet. d(1) should be the ' //
     +   'distance upstream from the basin outlet to '//
     +   'the most downstream point of the MRU such that ach(1)=0.0',
     +   'm').ne.0) return

      ALLOCATE (nchan_d(nchan))
      if(declparam('routec', 'nchan_d', 'nchan', 'integer',
     +   '5', '0', '50',
     +   'Number of topmodel channel routing increments '//
     +   'for each channel.',
     +   'Number of topmodel channel routing increments '//
     +   'for each channel.',
     +   'none').ne.0) return

      if(declparam('topc', 'dtinit', 'one', 'real',
     +   '24', '0', '24',
     +   'Initial timestep for initialize function.',
     +   'Initial timestep for initialize function.',
     +   'hours').ne.0) return

c
c irrig_hyd_next will contain the depth of irrigation to be pumped from
c a river segment on the next time step.
c
      ALLOCATE (irrig_hyd_next(nmru))
      if(declvar('routec', 'irrig_hyd_next', 'nmru', nmru, 'double', 
     + 'Diverted stream water to be applied as irrigation on '//
     $ 'the next time step','inches',irrig_hyd_next).ne.0) return
c
c irrig_hyd_seg_next will contain the total volume of water to be diverted
c from a stream segment to be applied to one or more MRUs on the next
c time step
c
      ALLOCATE (irrig_hyd_seg_next(nhydro))
      if(declvar('routec', 'irrig_hyd_seg_next', 'nhydro', nhydro,
     $ 'double', 
     + 'Volume of water diverted from a stream segment to '//
     $ 'be applied to one or more MRUs on the next time step.',
     $ 'm3',irrig_hyd_seg_next).ne.0) return
c
c irrig_int_max is the maximum discharge rate for a pump placed in
c a well or stream providing irrigation to an MRU
c
      ALLOCATE (irrig_int_max(nmru))
      if(declparam('topc', 'irrig_int_max', 'nmru', 'real',
     +   '100','0.', '1000000',
     +   'Maximum discharge rate for a pump placed in a well or '//
     $   'stream to provide irrigation to an MRU',
     +   'Maximum discharge rate for a pump placed in a well or '//
     $   'stream to provide irrigation to an MRU',
     +   'gallons per minute').ne.0) return 
c
c If greater than zero, irrig_sched_int points to a schedule (irrig_int_next)
c describing depths of irrigation to be applied from an internal source.
c The composition of the irrigated water depends on the source, irrig_int_src.
c
      ALLOCATE (irrig_sched_int(nmru))
      if(declparam('precip', 'irrig_sched_int', 'nmru', 'integer',
     +   '0', '0', '100',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'none').ne.0) return

c
c If there is an internal irrigation schedule for an MRU (irrig_sched_int>0),
c then irrig_int_src indicates if the MRU has recieves irrigations from a
c well in the MRU (value=0), or a stream segment (value > 0).
c
      ALLOCATE (irrig_int_src(nmru))
      if(declparam('precip', 'irrig_int_src', 'nmru', 'integer',
     +   '0', '0', '100',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   'none')
     +   .ne.0) return 

      ALLOCATE (irrig_int_init(nmru))
      if(declparam('precip', 'irrig_int_init', 'nmru', 'real',
     + '0', '0', '100',
     + 'Irrigation from an internal source to be applied on first day',
     + 'Irrigation from an internal source to be applied on first day',
     + 'inches').ne.0) return

      ALLOCATE (mru_area(nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      ALLOCATE (mru_area_frac(nmru))
      if(declparam('topc', 'mru_area_frac', 'nmru', 'real',
     +   '1', '0', '1',
     +   'Subcatchment area/total area',
     +   'Subcatchment area/total area',
     +   'none').ne.0) return

      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

!      if(decl*param('io', 'topout_file_unit', 'one', 'integer',
!     +   '80', '50', '99',
!     +   'Unit number for TOPMODEL output file',
!     +   'Unit number for TOPMODEL output file',
!     +   'integer').ne.0) return
!     
! local variables
      ALLOCATE (qscm(nchan))
      ALLOCATE (qsccfs(nchan))
      ALLOCATE (chan_area(nchan))
      ALLOCATE (chan_area_frac(nchan))
      ALLOCATE (qchanin(nchan))
      ALLOCATE (irrig_int_next(nirrig_int))
      ALLOCATE (qobs(nobs))
      ALLOCATE (irrig_hyd_mru(nmru))
      ALLOCATE (tch(ntopchan))
      ALLOCATE (tchmax(nchan))
      ALLOCATE (qavail(nhydro))

      routecdecl = 0

      return
      end

c***********************************************************************
c
c     routecinit - Initialize Clark routing module - get parameter values,
c

      integer function routecinit()

      USE WEBMOD_ROUTE
      USE WEBMOD_IO, ONLY : topout
      IMPLICIT NONE
C***  local variables
      integer is, NCH, i, j, IN, NR, IR, ND, TIME
      integer nrtdelay
      real DT, maxchsegs
      real CHVDT, A1, A2
      real SUMAR, SUM, chan_area_tot

c      character*135 output_path
c      logical filflg

      routecinit = 1

c----- set name for topmod unique output file 
c      ret = getoutname (output_path, '.topout')
c      inquire(file=output_path,exist=filflg)
c      if (filflg) then
c        open(unit=80,file=output_path,status='old')
c        close(unit=80,status='delete')
c      endif

c-----open the file.
c      open (unit=80,file=output_path,access='sequential',
c     * form='formatted', status='old')


      if(getparam('precip', 'irrig_int_src', nmru, 'integer',
     $     irrig_int_src).ne.0) return

      if(getparam('precip', 'irrig_sched_int', nmru, 'integer',
     $     irrig_sched_int).ne.0) return

      if(getparam('precip', 'irrig_int_init', nmru, 'real',
     $     irrig_int_init).ne.0) return

      if(getparam('topc', 'irrig_int_max', nmru, 'real',
     $     irrig_int_max).ne.0) return

      if(getparam('topc', 'dtinit', 1 , 'real', dtinit)
     +   .ne.0) return
           
      if(getparam('top2c', 'mru2chan', nmru, 'integer',
     $     mru2chan).ne.0) return

      if(getparam('routec', 'chv', 1, 'real', CHV)
     +   .ne.0) return

      if(getparam('routec', 'qobsta', 1, 'integer', qobsta)
     +   .ne.0) return

      if(getparam('routec', 'q0', 1, 'real', Q0)
     +   .ne.0) return
c
c chan_loss_rate is used to determine the maximum channel loss rate
c for any time-delay bin (nclarksegs). Those values are stored in the
c variable maxchanloss so it is not necessary to pass chan_loss_rate
c  back out of the init routine.
c
      if(getparam('routec', 'chan_loss_rate', 1, 'real',
     $     chan_loss_rate).ne.0) return

      if(getparam('routec', 'nchan_d', nchan, 'integer', nchan_d)
     +   .ne.0) return

      if(getparam('routec', 'ach', ntopchan*nchan, 'real', ACH)
     +   .ne.0) return

      if(getparam('routec', 'd', ntopchan*nchan, 'real', D)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area_frac', nmru, 'real',
     $     mru_area_frac).ne.0) return

      if(getparam('basin', 'basin_area', 1 , 'real', basin_area)
     +   .ne.0) return

!      if(get*param('io', 'topout_file_unit', 1, 'integer',
!     +   topout_file_unit).ne.0) return


c
c    This use of get*var in the init loop is OK since the variables are
c    static compilations of mru parameters calculated in the linking
c    module top2clark.f -RW (note: get*var is separated to avoid
c    tripping up the xmbuild parser, that looks for the real word)
c    
c
      if(getvar('top2c', 'chan_area_frac', nchan, 'real', 
     +  chan_area_frac).ne.0) return

      if(getvar('top2c', 'chan_area', nchan, 'real', 
     +  chan_area).ne.0) return

      DT = dtinit

c
c Write a section title
c
      WRITE(topout%lun,603)DT, CHV
 603  FORMAT(//'CHANNEL ROUTING DATA'/
     $     ' (Time step, DT =',f8.4,' hours)'/
     $     ' (Channel velocity, CHV =',f8.4,' meters per hour)'/
     $     ' Chan  Area(km^2)  MaxDel(DTs) Histogram Ordinates'/
     $     ' ==== ===========  =========== ======================>')
c
c The following definitions of RVDT and CHVDT are consistant
c with their function. The terms are reversed in the original
c code.- RW
c
c Q0DT is the depth of water per time step that will
c yield Q0, the initial discharge in cubic meters per second
c
*  RVDT is the distance the water flows each time step
*    in the subcatchment.
*  CHVDT is the distance the water flows each time step
*    from the subcatchment outlet to the basin outlet.
*
* Sept 2004 - eliminated RV and CHV dimensioned by nmru in 
*    favor of a single basin channel velocity CHV - RMTW
*
* May 2005 - initialize maxchsegs = 0.0
*
      
      maxchsegs = 0.0      
      Q0DT = Q0*DT*0.0036/basin_area
      chan_area_tot = 0.0
c
c clark_segs will retain the maximum value of nrtdelay for use in 
c other modules (i.e. webmod_res.f).
c
      clark_segs = 0

      do 50 is = 1, nchan

      qchanin_m(is) = q0dt

      NCH = nchan_d(is)
      chan_area_tot = chan_area_tot + chan_area(is)

      CHVDT = CHV * DT

c      T0DT = T0(is) + ALOG(DT)
*
*  CONVERT DISTANCE/AREA TO TIME DELAY HISTOGRAM ORDINATES
*
      TCH(1) = D(1,is)/CHVDT
      DO 15 J = 2,NCH !NCH is the the number of (d,ach) pairs for each channel
      TCH(J) = TCH(1) + (D(J,is) - D(1,is))/CHVDT !sum time to farthest point
 15   CONTINUE
      tchmax(is)=tch(nch)
      if (tchmax(is).gt.maxchsegs) maxchsegs = tchmax(is)
      NR = INT(TCH(NCH))
      IF(FLOAT(NR).LT.TCH(NCH))NR=NR+1
      ND = INT(TCH(1)) !ND are the number of DTs below the MRU outlet
      NR = NR - ND  !NR are the number of DTs in the MRU
      DO 20 IR=1,NR !Loop to drain water feeding stream in MRU
      TIME = ND+IR ! Step time by DTs beginning with ND+1
      IF(TIME.GT.TCH(NCH))THEN 
       AR(IR,is)=1.0
      ELSE
       DO 21 J=2,NCH
       IF(TIME.LE.TCH(J))THEN
          AR(IR,is)=ACH(J-1,is)+(ACH(J,is)-ACH(J-1,is))*(TIME-TCH(J-1))/
     1      (TCH(J)-TCH(J-1))
          GOTO 20
       ENDIF
   21    CONTINUE
      ENDIF
   20 CONTINUE
   
      nrsc(is) = NR
      ndsc(is) = ND
      
      A1= AR(1,is)
      SUMAR=AR(1,is)
      AR(1,is)=AR(1,is)*chan_area_frac(is)
      IF(NR.GT.1)THEN
       DO 22 IR=2,NR
       A2=AR(IR,is)
       AR(IR,is)=A2-A1
       A1=A2
       SUMAR=SUMAR+AR(IR,is)
       AR(IR,is)=AR(IR,is)*chan_area_frac(is)
   22    CONTINUE
      ENDIF


 50   continue

*
*  INITIALISE Q0 VALUES HERE
*  Q0 IS THE INITIAL DISCHARGE FOR THIS SUBCATCHMENT
*
*
* Changed Q0 to scalar parameter of discharge at outlet. Internally
* scaled to Q0DT, a constant depth of discharge across all MRUs - RW
*
c
c  Initialise discharge array and store the maximum value of nrtdelay
      do 32 is=1,nchan
       nrtdelay=ndsc(is)+nrsc(is)
       if(nrtdelay.gt.clark_segs) clark_segs = nrtdelay
!      do 28 I=1,nrtdelay
! 28   Q(I,is)=0.
      SUM=0.
      ND = ndsc(is)
      NR = nrsc(is)
      DO 29 I=1,ND
      Q(I,is) = Q(I,is) + Q0DT*chan_area_frac(is)
   29 continue
      DO 30 I=1,NR
c$$$      SUM=SUM+AR(I,is)     ! I think this needs to be below
c                                the Q() assignment, otherwise the
c                                farthest segment has no water
      IN = ND + I 
      Q(IN,is)=Q(IN,is)+Q0DT*(chan_area_frac(is)-SUM)
      SUM=SUM+AR(I,is)
   30 continue
c$$$   30 Q(IN,is)=Q(IN,is)+Q0DT*(chan_area_frac(is)-SUM)
   32 continue
c
c Ensure that user sets the dimension decribing the number of time delay
c ordinates to the number of time-delay ordinates
c
      if(nhydro.ne.clark_segs) then
         print*,'The dimension nhydro needs to be set to ',clark_segs,
     $        ', the number of stream time-delay ordinates. Change '//
     $        'the dimension value and rerun.'
         return
      end if
c
c Check to see if Diversions for irrigation exist.
c If they do, then increase initial water in each segment
c by the amount of the diversion to be applied on the first time
c step. Also check that no segment indicated as a source for the
c diversion exceeds the total number of time-delay segments.
c
      div_exists = .false. ! flag to see if diversion exists
c                           if not, the diversion calcs in the run
c                           section can be skipped
      do i = 1,nmru
         irrsched = irrig_sched_int(i)
         irrsrc = irrig_int_src(i)
         irrig_hyd_next(i) = 0.0
         if(irrsched.gt.0 .and. irrsrc.gt.0) then
            div_exists = .true.
c add initial extractions to channel volumes
            in = irrsrc
            is = mru2chan(i)
            Q(IN,is)=Q(IN,is)+irrig_int_init(i)*mru_area_frac(i)*
     $           inch2m
         end if
      end do
c
c Set initial irrigation from stream to zero (to avoid NaN,
c the value is ignored on the first time step in topmod anyway;
c it is assigned the value irrig_hyd_init.
c Compute the average depth of water lost to the deep aquifer
c system from each time-delay bin and also add that to the initial
c volume.
c
      do i = 1, clark_segs
         irrig_hyd_seg_next(i) = 0.0
         if(i.eq.clark_segs) then
            maxchanloss(i) =
     $        (maxchsegs-int(maxchsegs))*chan_loss_rate*chvdt*qlin2dep*
     $        dt/basin_area
         else
            maxchanloss(i) =
     $           chan_loss_rate*chvdt*qlin2dep*dt/basin_area
         end if
c$$$         Q(IN,is)=Q(IN,is)+irrig_int_init(i)*mru_area_frac(i)*
c$$$     $        inch2m
         sum_chan_loss(i) = 0.0
      end do
c
c Create matrix for printing AR values filled with zero values
c for the ND time step bins
c
      do 27 is= 1,nchan
      DO 26 I=1,clark_segs
         if(I.GT.ndsc(is).and.I.le.nrsc(is)+ndsc(is)) then
            ar_fill(i,is)=ar(i-ndsc(is),is)
         else
            ar_fill(I,is)=0.0
         end if
 26   continue
      WRITE(topout%lun,604)is, chan_area(is), TCHMAX(is),
     $     (AR_FILL(I,is),I=1,clark_segs)
 604  FORMAT(1X,I4,E12.5,(1X,10E12.5))
 27   continue

      write(topout%lun,605) chan_area_tot,
     $     sumar, chan_loss_rate, maxchsegs*chvdt/1609.34
 605  format(//,'Sum of areas drained to channels ', 
     $     f10.3,' sq.km. (should equal basin area)'/
     $     'Sum of histogram ordinates for all channels  ',
     $     f10.4, ' (should equal 1.0000)'/
     $     'Channel loss to deep aquifer ',
     $     f10.4, ' cfs/mile over a total ',
     $     'length of ',f10.4,' miles')


      routecinit = 0

      return
      end

c***********************************************************************
c
c     routecrun -
c

      integer function routecrun()

      USE WEBMOD_ROUTE
      USE WEBMOD_IO, ONLY : topout
      IMPLICIT NONE

C***  local variables
      integer NR, IR, ND, i, is, j, jj, k, nstep
      real irrig_dep_max
      double precision reduce, frac_divert,vol_m3
      logical end_run, end_yr, end_mo, end_dy, end_storm
      real DT, qbar, varq, vare, e

      routecrun = 1


      if(getvar('top2c', 'qchanin', nchan, 'real',
     +   qchanin).ne.0) return

      if(getvar('obsc', 'irrig_int_next', nirrig_int, 'real',
     $     irrig_int_next).ne.0) return

      if(getvar('obs', 'runoff', nobs, 'real',
     +   qobs).ne.0) return

      if(getvar('io', 'endper', 1, 'integer', endper)
     +   .ne.0) return

      if(getvar('precip', 'irrig_hyd_mru', nmru, 'real',
     $     irrig_hyd_mru).ne.0) return

      DT = deltim()

      nstep = getstep()

      if(nstep.eq.1) then

        QTOT = 0.
        F1 = 0.
        F2 = 0.
        SUMQOBS = 0.
        SSQ = 0. 
        NTSTEP = 0
!          do i = 1, clark_segs
!            do j = 1, nchan
!              q(i,j) = 0.0
!            end do
!          end do

      end if

c
c Decompose the endperiod variable
c
      end_run = .false.
      end_yr = .false.
      end_mo = .false.
      end_dy = .false.
      end_storm = .false.
c
c Save a few loops by testing the most common states first.
c Storms need to be tested independently. Others can just be tested
c for the end of the more frequent period. For example, if end_dy is
c not true then, by definition, neither is end_mo, end_yr, and end_run.
c
      if (endper.ne.0.and.mod(endper,2).ne.0) then
         endper = endper - 1
         end_storm = .true.
      end if
      
      if (endper.ne.0) then
         end_dy = .true.
         endper = endper - 2
         if (endper.ne.0) then
            end_mo = .true.
            endper = endper - 4
            if (endper.ne.0) then
               end_yr = .true.
               endper = endper - 8
               if (endper.ne.0) end_run = .true.
            end if
         end if
      end if

c
c Extract diversions proportional to depths and then
c reset channel exports and available water for irrigation
c
c$$$      do i = 1, clark_segs
c$$$         do is = 1,nchan
c$$$            q_tot(i) = q_tot(i) + q(i,is)
c$$$         end do
c$$$         q(i,is) = q(i,is) - q(i,is)/q_tot*
c$$$     $        (irrig_hyd_seg_next(i)/a_million/basin_area)
c$$$         q_tot(i) = 0.0
c$$$      end do


      do i = 1, clark_segs
         irrig_hyd_seg_next(i) = 0.0
         qavail(i) = 0.0
      end do

c
c Accumulate total irrigation requested for next time step from
c each stream segment, in cubic meters
c
      if(div_exists) then
         do i = 1, nmru
            irrsrc = irrig_int_src(i)
            if(irrsrc.gt.0) then
               irrsched = irrig_sched_int(i)
c     
c     Limit hillslope application by the pump capacity
c     
               irrig_dep_max = irrig_int_max(i)
     $              *gpm2m3ph/mru_area(i)/1e6/inch2m*dt ! max inches per time step
               if(irrig_int_next(irrsched).le.irrig_dep_max) then
                  irrig_hyd_next(i) = irrig_int_next(irrsrc)
               else
                 if(div_warn.eq.0) then
                   print*,'River diversion rate reduced because of '//
     $               'limited pumping capacity on time step', ntstep,
     $               ' for MRU ',i
                   print*,'This warning will not be repeated'
                   div_warn = 1
                 endif
                 irrig_hyd_next(i) = irrig_dep_max
               end if
               vol_m3 = irrig_int_next(irrsched)*inch2m*a_million*
     $              mru_area(i)
               irrig_hyd_seg_next(irrsrc) =
     $              irrig_hyd_seg_next(irrsrc)+vol_m3
            end if
         end do
      end if

C  START LOOP with clark_segs on outside and nchan on inside.
c  This is the opposite of the original order in TOPMODEL.
c  Note the usage of ar_fill instead of ar. This is more efficient
c  since we can distribute hillslope inputs and track total volumes
c  at the same time.

      chan_loss_basin = 0.0

      do 100 ir = 1,clark_segs

         chan_loss(ir) = 0.0
C
C  CHANNEL ROUTING CALCULATIONS
C  allow for time delay to catchment outlet ND as well as 
C  internal routing array. Place hillslope discharge
c  into stream isochrons
c
         DO 40 is=1,nchan
c     
c     Convert channel inputs in cubic meters to
c     average depth in meters - RW

            qchanin_m(is) = qchanin(is)/chan_area(is)/a_million
c     
            Q(ir,is)=Q(ir,is)+qchanin_m(is)*ar_fill(ir,is)
            qavail(ir) = qavail(ir)+q(ir,is)
 40      CONTINUE
c
c Remove channel leakage from stream first
c
         if(maxchanloss(ir).gt.0.0) then
            if(qavail(ir).gt.maxchanloss(ir))then
               chan_loss(ir) = maxchanloss(ir)
               sum_chan_loss(ir) = sum_chan_loss(ir) + chan_loss(ir)
               qavail(ir) = qavail(ir) - chan_loss(ir)
            else
               chan_loss(ir) = qavail(ir)
               sum_chan_loss(ir) = sum_chan_loss(ir) + chan_loss(ir)
               qavail(ir) = 0.0
            end if
         end if

         chan_loss_basin = chan_loss_basin + chan_loss(ir)
c
c Convert depth in channels to volume for irrigation calcs
c

         qavail(ir) = qavail(ir)*basin_area*a_million ! Total volume in m3

c
c  Limit diversions to availability
c
         if(qavail(ir).lt.irrig_hyd_seg_next(ir)) then
            reduce = qavail(ir)/irrig_hyd_seg_next(ir)
            if(.NOT.demand_warn) then
              print*,'Demands for irrigation from stream segment ',
     $           ir,' exceed available volume. Irrigation to MRUs '//
     $           'pumping from this segment have been reduced.'
              print*,' This warning will not be repeated.'
              demand_warn=.true.
            end if   
            do i = 1,nmru
               if(irrig_int_src(i).eq.ir)
     $              irrig_hyd_next(i)=reduce*irrig_hyd_next(i)
            end do
            irrig_hyd_seg_next(ir) = qavail(ir)
         end if

c
c reduce depths in each segment (stratified by nchan) in proportion to
c diversions to be removed on the next time step from each segment
c
         if(qavail(ir).gt.0.0) then
            frac_divert = irrig_hyd_seg_next(ir)/qavail(ir)
         else
            frac_divert = 0.0
         end if

         do 95 is =1,nchan
            q(ir,is) = q(ir,is)*(1.0-frac_divert)
 95      continue
 100  continue  ! End of clark_segs loop

      do 202 is = 1,nchan

         qscm(is) = q(1,is)/ chan_area_frac(is)
         qsccfs(is) = (qscm(is) * chan_area(is) * 9809.63)/dt
         qtot = qtot + qscm(is)*chan_area_frac(is)

c*****shift route timing matrix
         NR = nrsc(is)
         ND = ndsc(is)

         jj= NR + ND - 1
         k = 1
         if(jj.gt.0) then
            do 60 j=1, jj
               k= j + 1
               Q(j,is) = Q(k,is)
 60         continue
         endif
         Q(k,is)= 0.

 202  continue
c
c  Calculate basin area-weighted discharge
c
c  The 110 loop uses qbasinm to temporarily store the sum
c  of the subcatchment discharge, in millions of cubic meters.
c  Note that chan_area is in sq.km and dt in hours.
c
c  
      qbasinm = 0.0
      do 110 is = 1,nchan
         qbasinm = qbasinm + (qscm(is)*chan_area(is))
 110  continue
c
c  Now discharge in millions of cubic meters is converted to cfs,
c  and then to depth in meters by dividing by the total basin area
c  in sq.km (units of 1,000,000 sq.m)
c
      qbasincfs = (qbasinm * 9809.63) / dt
      qbasinm = qbasinm / basin_area
c
c Accumulate standard TOPMODEL objective function expressions
c
      SUMQOBS =SUMQOBS + QOBS(qobsta)
      SSQ = SSQ + (QOBS(qobsta))**2
      F1=F1 + (qbasincfs - QOBS(qobsta))**2
      F2=F2 + ABS(qbasincfs - QOBS(qobsta))

      NTSTEP = NTSTEP + 1
c
c Print run summary
c
      if(end_run) then
c
c-- compute statistics
c
         QBAR = SUMQOBS / NTSTEP
         VARQ = (SSQ/NTSTEP - QBAR*QBAR)
         VARE = F1/NTSTEP
         E=1-VARE/VARQ
c
c  add objective function values to output file
         write(topout%lun,621)f1,e,f2,qbar,varq,vare
c      write(10,621)f1,e,f2,qbar,varq,vare
 621     format(//1x,'Objective function values'/
     1        1x,'F1 ',e12.5,'   E ',f12.5,'   F2 ',e12.5//
     2        1x,'Mean Obs Q ',e12.5,'    Variance Obs Q ',e12.5/
     3        '     Error Variance',e12.5)

      end if
      
      routecrun = 0

      return
      end

c***********************************************************************
c
c     routecclean - Just use to close the .topout file since this
c                   is the last module that uses a cleanup routine
c

c$$$      integer function routecclean(topout_file_unit, sumqobs,
c$$$     +                  ntstep, ssq, f1, f2)
c$$$
c$$$      include 'fmodules.inc'
c$$$
c$$$      integer topout_file_unit
c$$$
c$$$      real sumqobs, ssq, f1, f2, qbar, varq, vare, e
c$$$      integer ntstep
c$$$
c$$$      routecclean = 1
c$$$
c$$$c
c$$$c-- compute statistics and print here in the cleanup for now
c$$$c
c$$$
c$$$      QBAR = SUMQOBS / NTSTEP
c$$$      VARQ = (SSQ/NTSTEP - QBAR*QBAR)
c$$$      VARE = F1/NTSTEP
c$$$      E=1-VARE/VARQ
c$$$c
c$$$c  add objective function values to output file
c$$$      write(topout_file_unit,621)f1,e,f2,qbar,varq,vare
c$$$c      write(10,621)f1,e,f2,qbar,varq,vare
c$$$  621 format(//1x,'Objective function values'/
c$$$     1 1x,'F1 ',e12.5,'   E ',f12.5,'   F2 ',e12.5//
c$$$     2 1x,'Mean Obs Q ',e12.5,'    Variance Obs Q ',e12.5/
c$$$     3 '     Error Variance',e12.5)
c$$$c
c$$$c
c$$$C  Close topmodel output file
c$$$c
c$$$      close(unit=topout_file_unit)
c$$$
c$$$      routecclean = 0
c$$$
c$$$      return
c$$$      end



