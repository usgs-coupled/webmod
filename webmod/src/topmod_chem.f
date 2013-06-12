***********************************************************************
c   topmod_chem.f
c
c   last modified by Rick Webb $Date: 2007-06-08 11:58:32 -0600 (Fri, 08 Jun 2007) $
c   - corrected routing routine (eliminated nrtdelay as dimension)
c   - corrected max infiltration check for infiltration excess
c   - added cumulative variables for different flow paths
c   - added direct flow varable, cumulative variable, 
c     and parameter coefficient per efforts of David Kinner.
c   - changed Q0 parameter from 'm/time step' to cubic meters per second.
c     Also reduced dimension of Q0 from nmru to one since most data only
c     includes the initial discharge at the outlet.
c   - made base flow, qb, a public variable
c   - made qout, the sum of the instantneous flow paths, a public variable
c   - created public variables, runoffm(meters) and precipm(meters) to 
c     permit comparison of model predictions (meters) with observed
c     rainfall (inches) and dicharge (cfs). Also made public the variable
c     p that indicates the influx into the topmodule routines.
c   - added additional basin summary variables for exporting to
c     the summary modules (qout_basin, qdf_basin, etc).
c   - modified the equation that calculates UZ by moving DT from the
c     denominator to the numerator.
c   - added saves for qddffrac, t0, chv, rv, s_root_depth, and sr0
c   30apr02 - added recharge to root zone from shallow ground water,
c     srzwet(nac,nmru). SRZWET is weighted by area and summed per
c     subcatchment per time step in qwet(nmru). Accumulated values are
c     in the variable sumqwet(nmru). SRZWET, the subcatchment wetting
c     per time step are weighted by area and reported in qwet_basin,
c     the basin average recharge of ground water to root zone, in meters.
c     QWET is removed from the saturated zone so it increases the
c     saturation deficit, SBAR.
c   11sep02 - Extracted channel routing section to route_clark.f.
c     Each MRU will now point to a specific nchan segment. The
c     velocities and ntopchan segments per nchan will now be
c     associated with the routing routine and not topmodg.f
c     Removed radiation plane in favor of mru-specific parameters
c     and variables. Removed initial discharge, q0, in favor of
c     sbar0(nmru), the initial saturation deficits in each MRU.
c   29oct02 - Move q0 parameter to routing module route_clark to 
c     initialize channel routing.
c   01nov02 - Opening and closing of the TOPMODEL output file will
c     take place in the io.f module since writing to that file occurs 
c     here and in the route_clark module.
c   4Mar03 - Added two additional arguments to be passed to the subroutine
c     EXPINF. They are CONST, the constant developed to solve the finite
c     integral, and TP, the time to ponding. These values were calculated
c     only on the first of consecutive periods of recharge and then lost
c     or confused on subsequent days if calculating for multiple basins.
c     Also initialized initial guesses for F and F1 in the EXPINF. The 
c     cumulative infiltration, F is set to CUMF, the total culative
c     infiltration from the previous step. F1 is set to 0 so that the
c     initial estimates for time to ponding are correctly calculated.
c     Finally, I eliminated ocassional negative values of REX that
c     resulted from machine precision and conversions of P/DT. The
c     absolute value of infiltration excess, REX, must now be greater
c     than 1e-8 meters for a given time step or else it is be set to zero.
c   6Mar03 - Added lognormal distribution for vertical hydraulic
c     conductivity, XK. Distribution described by XK0, the median
c     of the log normal population, and XK_CV, the coefficient
c     of variation in log space. If infiltration excess is to be
c     calculated (infex=1), then on days with rain or snow melt,
c     infiltration excess will be calculated for each of 9 lognormally
c     distributed values of XK (at integer values of standard deviations
c     from -4 to +4. If XK_CV=0.0 then only one value of X0 is used to
c     compute infiltration excess = XK0. This introduces the new parameter,
c     XK_CV and a new dimension, nxkbin, set to the fixed value of nine
c     bins (nxkbin).
c   10Mar03 - Added basin totals for infiltration excess, rex_basin, and
c     saturated overland flow, qofs_basin, in meters. Also added acm_basin
c     and afx_basin to indicate total proportion of basin with saturated
c     overland flow or infiltration excess, respectively, for a given 
c     time step. INFEX has been reduced from dimension nmru to dimension
c     one. This will require editing a previous parameter file to show
c     infex listed as below
c
c     ####
c     infex 8
c     1 
c     one
c     1
c     1 
c     1 
c     ####
c
c   20Jun03 - Scale actual ET according to cover density and whether
c     transpiration is occurring.(later use cover type). Cover density
c     calculated in intcp_prms.f responds to the transpiration periods
c     calculated in potet_hamon_prms.f - RMTW
c
c   04Sep03 - Added version control
c 
c   12Sep03 - Added soil parameters necessary for solute mass balances.
c             soil parameters dimensioned by nmru have a prefix s_
c     s_theta_wp = wilting point in cm3/cm3
c     s_theta_0 = initial soil moisture content of root zone (wp<theta<fc)
c     s_theta_fc = field capacity in cm3/cm3
c     s_porosity = effective porosity = saturated water content
c     s_root_depth = rooting depth in meters. Water available for ET
c                  is equal to theta_fc - theta_wp)*root_depth.
c          s_drain will be a variable containing a static value for
c             the readily drainable soil porosity equal to the effective
c             porosity, s_porosity, minus the field capacity, s_theta_fc.
c     s_satpref_zmin = Minimum water table elevation, in meters, at which
c                  the deep preferential flow paths (i.e. tile drains)
c                  become active.
c     s_satpref_zmax = Water table elevation, in meters, at which the
c                  discharge from preferential flow in the saturated zone
c                  will be at its maximum discharge.
c     s_satpref_K = conductivity, in cm/s of preferential flow in the
c                  saturated zone. Discharge from zone, in m/hr, is
c                  (wt_elev - zmin) * satpref_K * 36
c                    for zmin < wt_elev < zmax
c
c   22Sep03 - All precip on impermeable areas assumed to runoff during a given
c     timestep. Use hru_perv and hru_imperv, in sq.km to scale runoff before
c     making volumes available to webmod_res. Maybe later we will add TR55 or
c     something similar for better representations of impervious areas.
c
c     Add depth to bedrock parameter, s_rock_depth, so that preferential
c     flows through the unsaturated zone can be parameterized.
c
c     Local water tables, in meters above land surface (negative values
c     that may be added to a local datum to yield elevations) will be
c     reported in the new variable z_wt_local.
c
c   10Oct03 - Added uz_depth, the water content in the unsaturated zone
c     assuming field capacity in the UZ below the root zone.
c
c     Made EA, the local variable for ET flux, a public variable, sae_local,
c     of dimension (maxnac,nmru) so that webres can compute fluxes specific
c     to each topographic index bin.
c
c   21apr04 - Added the pmacro parameter to permit direct recharge of unsaturated
c     zone storage without the need to 1st reach field capacity. This was a
c     feature of some of the original TOPMODEL versions including the Dave
c     Wolock's topmongo code and Leon Kaufmann's New Jersey code.
c
c   27apr04 - Since pmacro is added to suz, the unsaturated zone storage,
c     add an addition fraction, pmac_sat, that will deliver an amount of
c     potential infiltration (after Hortonian overland flow removed) directly
c     to the saturated zone without mixing with unsaturated zone waters. The
c     new variable tracking the bypass is qvpref, to imply Vertical PREFerential
c     flow.
c
c   17jun04 - Add direct groundwater influx across MRU (and basin) boundary.
c     dimension ngw_ext, the number of unique flux rates. Two gw sources
c     can be assigned to simulate a leaking irrigation canal and upgradient
c     groundwater, if present. The influxes to a MRU is recorded in the variables
c     variable gw_in1/2(nmru), meters of inflow (volume normalized by area)
c     into each mru based on
c        parameter, sched_gw1/2(nmru) that points to 
c        a variable, gw_ext(ngw_ext), a unique flux time series in cfs per mile,
c        and parameter gwbnd_len1/2(nmru), the distance normal to gw inflow for
c        each MRU, in meters.
c     
c   16sep04 - Enable simulations of leakage from saturated zone. This permits
c     a catchment to be modeled with variable saturation above a hard pan with
c     irreversable leakage to a deep groundwater system.
c        The variable gw_loss(nmru) (and associated gwl_basin) is controlled
c        by the parameter gw_loss_k(nmru), the conductivity across the aquitard
c       (bedrock) and head above bedrock depth (rock_depth) for the MRU.
c
c   17apr09 - Add Fortran90 Module: WEBMOD_TOPMODULE
c
c
c***********************************************************************


c***********************************************************************
c! ***************** START MODULE *************
      MODULE WEBMOD_TOPMOD
      IMPLICIT NONE
      INCLUDE 'fmodules.inc'

C   Parameters, constants, and initial flags
c     Normalized areas for nine bins representing the lognormal
c     distribution of vertical conductivity (XK0)
      real, PARAMETER :: ack(9) = (/2.33000E-04, 5.97800E-03,
     $  6.05970E-02, 2.41730E-01, 3.82924E-01, 2.41730E-01,
     $  6.05970E-02, 5.97800E-03, 2.33000E-04/)
! conversion factors: 
! inch to meter
      real, PARAMETER :: inch2m = 0.0254
! gallons per minute to cubic meters per hour
      real, PARAMETER :: gpm2m3ph = 0.227125
!
!     Multiply gw inputs in cfs/mile * qlin2dep*dt(hr)*gwbnd_len(m)/mru_area(km2)
!     to area normalized deposition to the saturated zone, in meters 
!
      real, PARAMETER :: qlin2dep = 6.3343e-8
      logical step1
      integer irr_warn
      data irr_warn/0/
!   Dimensions
      integer, save :: nsc, nmru, nac, nobs
      integer, save :: nxkbin, ngw_ext, nirrig_int
!   Declared Parameters
      integer, save :: INFEX, iout, qobsta, topout_file_unit
      integer, save :: irrsched,irrsrc,gwsched
      integer, save, allocatable :: nacsc(:), cov_type(:)
      real, save :: basin_area, irrig_dep_max, pump_coeff
      real, save, allocatable :: SZM(:), TD(:), T0(:)
      real, save, allocatable :: XK0(:),HF(:),DTH(:)
      real, save, allocatable :: AC(:,:), TL(:), ST(:,:), SR0(:)
      real, save, allocatable :: sc_area(:), AREA(:) ! local copies of mru_area and mru_area_frac
      real, save, allocatable :: resp_hr_full(:), resp_hr_min(:)
      real, save, allocatable :: resp_coef_min(:), resp_coef(:)
      real, save, allocatable :: xk_cv(:),xk(:,:)
      real, save, allocatable ::  covden_sum(:), covden_win(:)
      real, save, allocatable ::  s_theta_fc(:), s_theta_wp(:)
      real, save, allocatable ::  s_theta_0(:), s_porosity(:)
      real, save, allocatable ::  s_root_depth(:), s_rock_depth(:)
      real, save, allocatable ::  s_satpref_zmin(:)
      real, save, allocatable ::  s_satpref_zmax(:), s_satpref_k(:)
      real, save, allocatable ::  pmacro(:), pmac_sat(:)
      integer, save, allocatable :: irrig_int_src(:), irrig_sched_int(:)
      integer, save, allocatable ::  sched_gw1(:),sched_gw2(:)
      real, save, allocatable ::   irrig_sat_next(:), sbar_max(:)
      real, save, allocatable ::   irrig_int_max(:), gw_loss_k(:)
      real, save, allocatable ::   gwbnd_len1(:),gwbnd_len2(:)
      real, save, allocatable :: riparian_thresh(:)  ! used to aggregate riparian and upland TWIs
      real, save, allocatable :: uz_area(:,:)  ! areas of riparian, and upland TWI bins
C   Declared Variables
      real, save, allocatable :: SBAR(:),SRMAX(:),SD(:,:),QB(:),REX(:)
      real, save, allocatable :: QOF(:), QUZ(:), QSCM(:),QOFS(:)
      real, save, allocatable :: SZQ(:), SUZ(:,:), SRZ(:,:)
      real, save, allocatable :: srzwet(:,:), sae(:), cumf(:,:)
      real, save, allocatable :: BAL(:),SUMP(:), SUMAE(:)
      real, save, allocatable :: SUMQ(:), ACMAX(:), acm(:), afxmax(:)
      real, save, allocatable :: qdf(:), qdffrac(:)
      real, save, allocatable :: qsccfs(:), qout(:)
      real, save, allocatable :: last_z_wt_local(:,:)
      real, save :: gw_in1_basin, gw_in2_basin, gwl_basin
      real, save :: qbasinm, qbasincfs, basin_actet, srz_basin
      real, save :: runoffm, precipm, suz_basin, smav_basin
      real, save :: smcont_basin, sbar_basin,rex_basin, qofs_basin
      real, save :: qof_basin, qdf_basin, qb_basin, qout_basin
      real, save :: basin_infil,basin_surfdep, qexfil_basin, qwet_basin
      real, save :: afx_basin, acm_basin, z_wt_basin,qvpref_basin
      real, save :: qpref_basin, basin_uz2sat
      real, save, allocatable :: sumqrex(:), sumqofs(:), sumqb(:)
      real, save, allocatable :: sumquz(:), sumqdf(:)
      real, save, allocatable :: srz_sc(:), suz_sc(:), infil(:)
      real, save, allocatable :: smav_sc(:),smcont_sc(:), surfdep(:)
      real, save, allocatable :: qexfil(:), sumqexfil(:)
      real, save, allocatable :: qwet(:), sumqwet(:)
      real, save, allocatable :: afx(:), z_wt_local(:,:),z_wt(:)
      real, save, allocatable :: qvpref(:), sumqvpref(:)
      real, save, allocatable :: qpref(:), sumqpref(:)
      real, save, allocatable :: qpref_max(:),s_drain(:),uz_depth(:,:)
      real, save, allocatable :: sae_local(:,:), quz_local(:,:)
      real, save, allocatable :: uz_infil(:,:), uz2sat(:,:)
      real, save, allocatable :: gw_in1(:),gw_in2(:)
      real, save, allocatable :: sumgw_in1(:),sumgw_in2(:)
      real, save, allocatable :: gw_loss(:), sumgw_loss(:)
C getvars from other modules
      real, save, allocatable :: QOBS(:), potet(:), gw_ext(:)
      real, save, allocatable :: irrig_int_next(:), irrig_sat_mru(:)
      real, save, allocatable :: snowmelt(:), psoilmru(:)
      real, save, allocatable :: intcp_evap(:), snow_evap(:)
      integer, save, allocatable :: transp_on(:)
C   Local for init section
      integer, save, allocatable :: atag(:)
      real, save, allocatable ::  irrig_int_init(:), sbar0(:)
      logical, save, allocatable :: riparian(:,:)
C   Local for run section
      real, save :: DT, XKAREA, sbar_norm, sat_head, rate
      real, save :: POBS, P
      real, save :: z_wt_quick, qprefwt
      real, save, allocatable :: last_uz_dep(:,:), last_srz(:,:)
      real, save, allocatable :: last_suz(:,:)
      real, save, allocatable :: EX(:)
      integer, save, allocatable :: nxk(:)
      integer, save, allocatable :: irof(:,:),ihour(:,:)
      real, save, allocatable :: tp(:,:)
      double precision, save :: scconst
      double precision, save, allocatable :: const(:,:)

      END MODULE WEBMOD_TOPMOD

c***********************************************************************
c
c     Main topcomp routine
c

      integer function topmod_chem(arg)
      implicit none
      include 'fmodules.inc'

      character(LEN=*), INTENT(IN) :: arg
      CHARACTER*256 SVN_ID

! Functions
      INTEGER, EXTERNAL :: topmdecl,topminit,topmrun

      topmod_chem = 0

      SVN_ID = 
     $     '$Id: topmod_chem.f 38 2007-06-08 17:58:32Z rmwebb $ '
      
      if(arg.eq.'declare') then
        topmod_chem = topmdecl()

      else if(arg.eq.'initialize') then
        topmod_chem = topminit()

      else if(arg.eq.'run') then
        topmod_chem = topmrun()

      end if

      END FUNCTION topmod_chem
c 
c     topmdecl - declare variables and parameters for TOPMOD_CHEM
c

      integer function topmdecl()

      USE WEBMOD_TOPMOD


! Get dimensions

      nmru = getdim('nmru')
        if ( nmru.eq.-1 ) return
      nsc = nmru
      nac = getdim('nac')
        if ( nac.eq.-1 ) return
      nxkbin = getdim('nxkbin')
        if ( nxkbin.eq.-1 ) return
      nobs = getdim('nobs')
        if ( nobs.eq.-1 ) return
      ngw_ext = getdim('ngw_ext')
        if ( ngw_ext.eq.-1 ) return
      nirrig_int = getdim('nirrig_int')
        if ( nirrig_int.eq.-1 ) return


      Topmdecl = 1

      ALLOCATE (BAL(Nmru))
      if(declpri('topmdecl_bal', nmru, 'real', BAL) .ne. 0) return

      ALLOCATE (SZQ(Nmru))
      if(declpri('topmdecl_szq', nmru, 'real', SZQ) .ne. 0) return 

      ALLOCATE (IROF(Nmru,nxkbin))
      if(declpri('topmdecl_irof', nmru*nxkbin, 'integer', IROF)
     + .ne. 0) return 

      ALLOCATE (CUMF(Nmru,nxkbin))
      if(declpri('topmdecl_cumf', nmru*nxkbin, 'real', CUMF)
     + .ne. 0) return 

      ALLOCATE (const(Nmru,nxkbin))
      if(declpri('topmdecl_const', nmru*nxkbin, 'double', const)
     + .ne. 0) return

      ALLOCATE (tp(Nmru,nxkbin))
      if(declpri('topmdecl_tp', nmru*nxkbin, 'real', tp)
     + .ne. 0) return 

      ALLOCATE (nxk(Nmru))
      if(declpri('topmdecl_nxk', nmru, 'integer', nxk)
     + .ne. 0) return 

      ALLOCATE (afxmax(Nmru))
      if(declpri('topmdecl_afxmax', nmru, 'real', afxmax)
     + .ne. 0) return 

      ALLOCATE (xk(Nmru, nxkbin))
      if(declpri('topmdecl_xk', nmru*nxkbin, 'real', xk)
     + .ne. 0) return

c
c qw_inflow is the depth of ground water inputs into an MRU from
c sources outside of the surface water basin, in meters. Its
c derived from gw_ext, flux data from the mms data file, in cfs/mile
c and the boundary length for the MRU, gwbnd_len
c
      ALLOCATE (gw_in1(Nmru))
      if(declvar('topc','gw_in1','nmru',nmru,'real',
     $     'Ground water influx from 1st gw source external to basin',
     $     'm',gw_in1) .ne.0) return

      ALLOCATE (gw_in2(Nmru))
      if(declvar('topc','gw_in2','nmru',nmru,'real',
     $     'Ground water influx from 2nd gw source external to basin',
     $     'm',gw_in2) .ne.0) return

      ALLOCATE (sumgw_in1(Nmru))
      if(declvar('topc', 'sumgw_in1', 'nmru', nmru, 'real', 
     + 'Sum of groundwater flux into MRU from 1st gw source',
     +  'meters',sumgw_in1).ne.0) return

      ALLOCATE (sumgw_in2(Nmru))
      if(declvar('topc', 'sumgw_in2', 'nmru', nmru, 'real', 
     + 'Sum of groundwater flux into MRU from 2nd gw source',
     +  'meters',sumgw_in2).ne.0) return

      if(declvar('topc', 'gw_in1_basin', 'one', 1, 'real', 
     + 'Area-weighted influx of groundwater into basin from '//
     $ '1st gw source','meters',gw_in1_basin).ne.0) return

      if(declvar('topc', 'gw_in2_basin', 'one', 1, 'real', 
     + 'Area-weighted influx of groundwater into basin from '//
     $ '2nd gw source','meters',gw_in2_basin).ne.0) return
c
c qprefmax is the maximum discharge from the preferential flow
c through the saturated zone, in meters per hour
c
      ALLOCATE (qpref_max(Nmru))
      if(declvar('topc','qpref_max','nmru',nmru,'real',
     $     'Maximum discharge for preferential flow '//
     $     'through the saturated zone for each MRU',
     $     'm/hr',qpref_max) .ne.0) return
c
c standard topmodel variables
c
      ALLOCATE (SBAR(Nmru))
      if(declvar('topc', 'sbar', 'nmru', nmru, 'real', 
     + 'average saturation deficit for the subcatchment',
     + 'meters',SBAR).ne.0) return

      if(declvar('topc', 'sbar_basin', 'one', 1, 'real', 
     + 'Area-weighted mean saturation deficit for the basin',
     + 'meters',sbar_basin).ne.0) return

      ALLOCATE (srz_sc(Nmru))
      if(declvar('topc', 'srz_sc', 'nmru', nmru, 'real', 
     + 'Mean root zone deficit for subcatchment',
     + 'meters',srz_sc).ne.0) return

      ALLOCATE (suz_sc(Nmru))
      if(declvar('topc', 'suz_sc', 'nmru', nmru, 'real', 
     + 'Mean unsaturated zone storage for subcatchment',
     + 'meters',suz_sc).ne.0) return

      if(declvar('topc', 'srz_basin', 'one', 1, 'real', 
     + 'Basin area-weighted mean root zone deficit',
     + 'meters',srz_basin).ne.0) return

      if(declvar('topc', 'suz_basin', 'one', 1, 'real', 
     + 'Basin area-weighted mean unsaturated zone storage',
     + 'meters',suz_basin).ne.0) return

      if(declvar('topc', 'smav_basin', 'one', 1, 'real', 
     + 'Basin area-weighted mean root zone moisture storage.'//
     + ' Equal to SRMAX - srz_basin + residual water','meters'
     $  ,smav_basin).ne.0) return

      if(declvar('topc', 'smcont_basin', 'one', 1, 'real', 
     + 'Basin area-weighted soil moisture content.'//
     + ' Equal to smav/s_root_depth','cm3/cm3',
     +  smcont_basin).ne.0) return

      ALLOCATE (smav_sc(Nmru))
      if(declvar('topc', 'smav_sc', 'nmru', nmru, 'real', 
     + 'Mean root zone moisture storage for each subcatchment.'//
     + ' Equal to SRMAX - srz_sc + residual water','meters',
     $ smav_sc).ne.0) return

      ALLOCATE (smcont_sc(Nmru))
      if(declvar('topc', 'smcont_sc', 'nmru', nmru, 'real', 
     + 'Mean root soil moisture content for each subcatchment.'//
     + ' Equal to smav_sc/s_root_depth','cm3/cm3',
     +  smcont_sc).ne.0) return

      ALLOCATE (sae_local(nac,Nmru))
      if(declvar('topc', 'sae_local', 'nac,nmru', nac*nmru, 'real', 
     + 'evapotranspiration from the root zone for '//
     + 'each topographic index bin',
     + 'meters',sae_local).ne.0) return

      ALLOCATE (sae(Nmru))
      if(declvar('topc', 'sae', 'nmru', nmru, 'real', 
     + 'evapotranspiration from the root zone for '//
     + 'each subcatchment',
     + 'meters',sae).ne.0) return

      if(declvar('topc', 'basin_actet', 'one', 1, 'real', 
     + 'evapotranspiration from the root zone for '//
     + 'the entire basin',
     + 'meters',basin_actet).ne.0) return

      ALLOCATE (SD(nac,Nmru))
      if(declvar('topc', 'sd', 'nac,nmru', nac*nmru, 'real', 
     + 'local saturation deficit for each ln(a/tanB) increment',
     + 'meters',
     + SD).ne.0) return

c
c new variables for revised summary module
c
      ALLOCATE (infil(Nmru))
      if(declvar('topc', 'infil', 'nmru', nmru, 'real', 
     + 'infiltration into soil for each subcatchment',
     + 'meters',infil).ne.0) return

      if(declvar('topc', 'basin_infil', 'one', 1, 'real', 
     + 'Average infiltration for the basin',
     + 'meters',basin_infil).ne.0) return

      ALLOCATE (surfdep(Nmru))
      if(declvar('topc', 'surfdep', 'nmru', nmru, 'real', 
     + 'Surface deposition (rain+snowmelt) for each subcatchment',
     + 'meters',surfdep).ne.0) return

      if(declvar('topc', 'basin_surfdep', 'one', 1, 'real', 
     + 'Average surface deposition (rain+snowmelt) for the basin',
     + 'meters',basin_surfdep).ne.0) return

c
c new or modified variables for chem module
c
c
c srmax will be calculated as (s_theta_fc - s_theta_wp) * s_root_depth
c
      ALLOCATE (srmax(Nmru))
      if(declvar('topc', 'srmax', 'nmru', nmru, 'real',
     $ 'Storage between field capacity and wilting point. '//
     $ 'Derived as (s_th_fc - s_th_wp)* s_root_depth.',
     + 'm',
     $ srmax).ne.0) return

      ALLOCATE (tl(Nmru))
      if(declvar('topc', 'tl', 'nmru', nmru, 'real', 
     + 'Mean topographic index for MRU',
     + 'ln(a/tanB)',
     $ tl).ne.0) return

      ALLOCATE (sr0(Nmru))
      if(declvar('topc', 'sr0', 'nmru', nmru, 'real', 
     + 'initial root-zone soil moisture deficit for each MRU',
     + 'm',
     $ sr0).ne.0) return

      ALLOCATE (resp_coef(Nmru))
      if(declvar('topc', 'resp_coef', 'nmru', nmru, 'real', 
     + 'static value for response coefficient that throttles '//
     + 'the redistribution of hillslope moisture content '//
     + '(sd adjustment). Derived from response parameters.',
     + 'none',
     $ resp_coef).ne.0) return

      ALLOCATE (s_drain(Nmru))
      if(declvar('topc', 's_drain', 'nmru', nmru, 'real', 
     + 'Readily drainable soil porosity. Equal to the '//
     $ 'the effective porosity, s_porosity, minus the '//
     $ 'moisture content at field capacity, s_theta_fc.',
     + 'cm3/cm3',
     $ s_drain).ne.0) return

      ALLOCATE (uz_depth(nac,Nmru))
      if(declvar('topc', 'uz_depth', 'nac,nmru', nac*nmru,
     $ 'real','Total moisture content in the unsaturated zone '//
     $ 'at the end of the previous time step (resid water + '//
     $ 'UZ storage + root zone water',
     $ 'meters',
     $ uz_depth).ne.0) return

      ALLOCATE (uz_infil(nac,Nmru))
      if(declvar('topc', 'uz_infil', 'nac,nmru', nac*nmru,
     $ 'real','Local infiltration of unsaturated zone '//
     $ 'by precipitation and snowmelt.Does not include '//
     $ 'pmacro flow that bypasses the uz to become recharge.',
     $ 'meters',uz_infil).ne.0) return

      ALLOCATE (quz_local(nac,Nmru))
      if(declvar('topc', 'quz_local', 'nac,nmru', nac*nmru,
     $ 'real','Local recharge from the unsaturated zone '//
     $ 'to the shallow preferential flow and the saturated '//
     $ 'zone', 'meters', quz_local).ne.0) return

      ALLOCATE (z_wt_local(nac,Nmru))
      if(declvar('topc', 'z_wt_local', 'nac,nmru', nac*nmru,
     $ 'real','local water table height above land surface for '//
     $ 'each ln(a/tanB) increment',
     + 'meters',
     + z_wt_local).ne.0) return

      ALLOCATE (z_wt(Nmru))
      if(declvar('topc', 'z_wt', 'nmru', nmru,
     $ 'real','average water table height above land surface for '//
     $ 'each MRU',
     + 'meters',
     + z_wt).ne.0) return

      if(declvar('topc', 'z_wt_basin', 'one', 1,
     $ 'real','average water table height above land surface for '//
     $ 'the entire basin',
     + 'meters',
     + z_wt_basin).ne.0) return

      ALLOCATE (SUZ(nac,Nmru))
      if(declvar('topc', 'suz', 'nac,nmru', nac*nmru, 'real', 
     + 'unsaturated zone storage for each ln(a/tanB) increment',
     + 'meters',
     + SUZ).ne.0) return

      ALLOCATE (SRZ(nac,Nmru))
      if(declvar('topc', 'srz', 'nac,nmru', nac*nmru, 'real', 
     + 'root zone deficit for each ln(a/tanB) increment',
     + 'meters',SRZ)
     +   .ne.0) return

      ALLOCATE (uz2sat(nac,Nmru))
      if(declvar('topc', 'uz2sat', 'nac,nmru', nac*nmru,
     $ 'real','Transfer of water from unsaturated zone to the '//
     $ 'saturated zone resulting from a change in water '//
     $ 'table.', 'meters', uz2sat).ne.0) return

      if(declvar('topc', 'basin_uz2sat', 'one', 1,
     $ 'real','Basin average transfer of water from unsaturated '//
     $ 'zone to the '//
     $ 'saturated zone resulting from a change in water '//
     $ 'table.', 'meters', basin_uz2sat).ne.0) return

      ALLOCATE (qb(Nmru))
      if(declvar('topc', 'qb', 'nmru', nmru, 'real', 
     + 'Discharge from saturated zone for each subcatchment',
     +  'meters',qb).ne.0) return

      if(declvar('topc', 'qb_basin', 'one', 1, 'real', 
     + 'Area-weighted discharge from saturated zone for the basin',
     +  'meters',qb_basin).ne.0) return

      ALLOCATE (gw_loss(Nmru))
      if(declvar('topc', 'gw_loss', 'nmru', nmru, 'real', 
     + 'Loss of water from saturated to deep aquifer '//
     $ 'for each subcatchment', 'meters', gw_loss).ne.0) return

      ALLOCATE (sumgw_loss(Nmru))
      if(declvar('topc', 'sumgw_loss', 'nmru', nmru, 'real', 
     + 'Sum of water lost from saturated zone to deep aquifer '//
     $ 'for each subcatchment', 'meters',sumgw_loss).ne.0) return

      if(declvar('topc', 'gwl_basin', 'one', 1, 'real', 
     + 'Area-weighted loss of water from saturated zone to '//
     $ 'the deep aquifer for the basin',
     $ 'meters',gwl_basin).ne.0) return

      ALLOCATE (qexfil(Nmru))
      if(declvar('topc', 'qexfil', 'nmru', nmru, 'real', 
     + 'Exfiltration discharge from oversaturated zones '//
     + 'for each subcatchment',
     +  'meters',qexfil).ne.0) return

      ALLOCATE (sumqexfil(Nmru))
      if(declvar('topc', 'sumqexfil', 'nmru', nmru, 'real', 
     + 'Sum of exfiltration discharge from oversaturated zones '//
     + 'for each subcatchment',
     +  'meters',sumqexfil).ne.0) return

      if(declvar('topc', 'qexfil_basin', 'one', 1, 'real', 
     + 'Area-weighted exfiltration discharge from saturated'//
     + 'zone for the basin',
     +  'meters',qexfil_basin).ne.0) return

c
c irrig_sat_next will contain the depth of irrigation to be pumped from
c the saturated zone on this time step to be applied on the mru during
c the next time step.
c
      ALLOCATE (irrig_sat_next(Nmru))
      if(declvar('topc', 'irrig_sat_next', 'nmru', nmru, 'real', 
     + 'Ground water to be pumped and applied as irrigation on '//
     $ 'the next time step','inches',irrig_sat_next).ne.0) return
c
c parameters to include gw influx from outside of basin boundary.
c A maximum of 2 sources can be defined for each MRU. This could be
c a leaky irrigation canal crossing the MRU and/or regional groundwater
c influx.
c
      ALLOCATE (sched_gw1(Nmru))
      if(declparam('topc', 'sched_gw1', 'nmru', 'integer',
     +   '0','bounded', 'ngw_ext',
     +   'Index of time series describing ground water inputs into '//
     $   'the MRU from the 1st gw source outside the basin; 0 if none',
     +   'Index of time series describing ground water inputs into '//
     $   'the MRU from the 1st gw source outside the basin; 0 if none',
     +   'none').ne.0) return 

      ALLOCATE (sched_gw2(Nmru))
      if(declparam('topc', 'sched_gw2', 'nmru', 'integer',
     +   '0','bounded', 'ngw_ext',
     +   'Index of time series describing ground water inputs into '//
     $   'the MRU from the 2nd gw source outside the basin; 0 if none',
     +   'Index of time series describing ground water inputs into '//
     $   'the MRU from the 2nd gw source outside the basin; 0 if none',
     +   'none').ne.0) return 

      ALLOCATE (gwbnd_len1(Nmru))
      if(declparam('topc', 'gwbnd_len1', 'nmru', 'real',
     +   '0','0', '1000000',
     +   'Length of MRU cross-section receiving ground water influx'//
     $   ' from 1st gw source outside the basin',
     +   'Length of MRU cross-section receiving ground water influx'//
     $   ' from 1st gw source outside the basin',
     +   'meters').ne.0) return 

      ALLOCATE (gwbnd_len2(Nmru))
      if(declparam('topc', 'gwbnd_len2', 'nmru', 'real',
     +   '0','0', '1000000',
     +   'Length of MRU cross-section receiving ground water influx'//
     $   ' from 2nd gw source outside the basin',
     +   'Length of MRU cross-section receiving ground water influx'//
     $   ' from 2nd gw source outside the basin',
     +   'meters').ne.0) return 
c
c irrig_int_max is the maximum discharge rate for a pump placed in
c a well or stream providing irrigation to an MRU
c
      ALLOCATE (irrig_int_max(Nmru))
      if(declparam('topc', 'irrig_int_max', 'nmru', 'real',
     +   '100','0.', '1000000',
     +   'Maximum discharge rate for a pump placed in a well or '//
     $   'stream to provide irrigation to an MRU',
     +   'Maximum discharge rate for a pump placed in a well or '//
     $   'stream to provide irrigation to an MRU',
     +   'gallons per minute').ne.0) return 
c
c irrig_int_src indicates if the MRU has recieves irrigations from a
c well in the MRU (value=0), or a stream segment (value > 0).
c If irrigation is received from a well or stream internal
c to the basin, then the irrigation schedule describing the depths to be
c applied during the subsequent time period is recorded in irrig_sched_int.
c
      ALLOCATE (irrig_int_src(Nmru))
      if(declparam('precip', 'irrig_int_src', 'nmru', 'integer',
     +   '0', '0', '100',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   ' 0 Irrigation from well in MRU; '//
     $   '>0 Drainage segment ID that will provide irrigation water',
     +   'none')
     +   .ne.0) return 

      ALLOCATE (irrig_sched_int(Nmru))
      if(declparam('precip', 'irrig_sched_int', 'nmru', 'integer',
     +   '0', '0', '100',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'Index of internal irrigation schedule for MRU; 0 if none',
     +   'none').ne.0) return

      ALLOCATE (irrig_int_init(Nmru))
      if(declparam('precip', 'irrig_int_init', 'nmru', 'real',
     + '0', '0', '100',
     + 'Irrigation from an internal source to be applied on first day',
     + 'Irrigation from an internal source to be applied on first day',
     + 'inches').ne.0) return

      ALLOCATE (srzwet(nac,Nmru))
      if(declvar('topc', 'srzwet', 'nac,nmru', nac*nmru, 'real', 
     + 'depth of water transferred from the ground water to the '//
     + 'root zone deficit for each ln(a/tanB) increment',
     + 'meters',srzwet).ne.0) return

      ALLOCATE (qwet(Nmru))
      if(declvar('topc', 'qwet', 'nmru', nmru, 'real', 
     + 'Root zone recharge from shallow ground water '//
     + 'for each subcatchment',
     +  'meters',qwet).ne.0) return

      ALLOCATE (sumqwet(Nmru))
      if(declvar('topc', 'sumqwet', 'nmru', nmru, 'real', 
     + 'Sum of root zone recharge from shallow ground water '//
     + 'for each subcatchment',
     +  'meters',sumqwet).ne.0) return

      if(declvar('topc', 'qwet_basin', 'one', 1, 'real', 
     + 'Area-weighted root zone recharge from shallow ground '//
     + 'water for the basin',
     +  'meters',qwet_basin).ne.0) return

      ALLOCATE (qout(Nmru))
      if(declvar('topc', 'qout', 'nmru', nmru, 'real', 
     + 'Sum of overland flow, direct flow, and discharge'//
     + ' from the saturated zone for each subcatchment',
     +  'meters',qout).ne.0) return

      if(declvar('topc', 'qout_basin', 'one', 1, 'real', 
     + 'Area-weighted sum of overland flow, direct flow, '//
     + 'and discharge from the saturated zone for the basin',
     +  'meters',qout_basin).ne.0) return

      if(declvar('topc', 'p', 'one', 1, 'real', 
     + 'Sum of precipitation and snowmelt available '//
     + 'for input into the TOPMOD routines.',
     +  'meters',p).ne.0) return

      if(declvar('topc', 'runoffm', 'one', 1, 'real', 
     + 'Observed discharge in metric units'//
     + ' for comparing with predicted flows',
     +  'meters',runoffm).ne.0) return

      if(declvar('topc', 'precipm', 'one', 1, 'real', 
     + 'Observed precipitation in metric units '//
     + 'for comparing with predicted flows',
     +  'meters',precipm).ne.0) return

      ALLOCATE (sumqrex(Nmru))
      if(declvar('topc', 'sumqrex', 'nmru', nmru, 'real', 
     + 'Sum of predicted infiltration excess for each subcatchment',
     +  'meters',sumqrex).ne.0) return

      ALLOCATE (sumqofs(Nmru))
      if(declvar('topc', 'sumqofs', 'nmru', nmru, 'real', 
     +'Sum of predicted saturated overland flow for each subcatchment',
     +  'meters',sumqofs).ne.0) return

      ALLOCATE (sumqb(Nmru))
      if(declvar('topc', 'sumqb', 'nmru', nmru, 'real', 
     + 'Sum of saturated zone discharge for each subcatchment',
     +  'meters',sumqb).ne.0) return

      ALLOCATE (sumquz(Nmru))
      if(declvar('topc', 'sumquz', 'nmru', nmru, 'real', 
     + 'Sum of predicted flux from unsaturated zone'//
     + 'to the saturated zone for each subcatchment',
     +  'meters',sumquz).ne.0) return

      ALLOCATE (sumqdf(Nmru))
      if(declvar('topc', 'sumqdf', 'nmru', nmru, 'real', 
     + 'Sum of predicted direct flow for each subcatchment',
     +  'meters',sumqdf).ne.0) return

      ALLOCATE (qdf(Nmru))
      if(declvar('topc', 'qdf', 'nmru', nmru, 'real', 
     + 'Predicted direct flow for each subcatchment',
     +  'meters',qdf).ne.0) return

      if(declvar('topc', 'qdf_basin', 'one', 1, 'real', 
     + 'Area-weighted predicted direct flow for the basin',
     +  'meters',qdf_basin).ne.0) return


      ALLOCATE (qvpref(Nmru))
      if(declvar('topc', 'qvpref', 'nmru', nmru, 'real', 
     + 'Predicted vertical preferential flow bypassing '//
     $ 'the unsaturated zone for each subcatchment.',
     +  'meters',qvpref).ne.0) return

      ALLOCATE (sumqvpref(Nmru))
      if(declvar('topc', 'sumqvpref', 'nmru', nmru, 'real', 
     + 'Sum of vertical preferential flow that bypasses the '//
     $ 'unsaturated zone.','meters',sumqvpref).ne.0) return

      if(declvar('topc', 'qvpref_basin', 'one', 1, 'real', 
     + 'Area-weighted vertical preferential flow for the basin',
     +  'meters',qvpref_basin).ne.0) return

      ALLOCATE (qpref(Nmru))
      if(declvar('topc', 'qpref', 'nmru', nmru, 'real', 
     + 'Predicted deep preferential flow for each subcatchment. '//
     $ 'A function of s_satpref_zmin, s_satpref_zmax, and s_satpref_k',
     +  'meters',qpref).ne.0) return

c$$$      if(decl*var('topc', 'qpref', 'nmru', nmru, 'real', 
c$$$     + 'Predicted deep preferential flow for each subcatchment.',
c$$$     +  'meters',qpref).ne.0) return

      ALLOCATE (sumqpref(Nmru))
      if(declvar('topc', 'sumqpref', 'nmru', nmru, 'real', 
     + 'Sum of deep preferential flow in the saturated zone for '//
     $ 'each subcatchment','meters',sumqpref).ne.0) return

      if(declvar('topc', 'qpref_basin', 'one', 1, 'real', 
     + 'Area-weighted deep preferential flow for the basin',
     +  'meters',qpref_basin).ne.0) return

      ALLOCATE (acm(Nmru))
      if(declvar('topc', 'acm', 'nmru', nmru, 'real', 
     + 'Contributing area for each subcatchment',
     +  'decimal percent',acm).ne.0) return

      if(declvar('topc', 'acm_basin', 'one', 1, 'real', 
     + 'Proportion of basin with saturated overland flow',
     +  'decimal percent',acm_basin).ne.0) return

      ALLOCATE (afx(Nmru))
      if(declvar('topc', 'afx', 'nmru', nmru, 'real', 
     + 'Area experiencing infiltration excess '//
     + 'for each subcatchment',
     +  'decimal percent',afx).ne.0) return

      if(declvar('topc', 'afx_basin', 'one', 1, 'real', 
     + 'Proportion of basin with infiltration excess',
     +  'decimal percent',afx_basin).ne.0) return

      ALLOCATE (QOF(Nmru))
      if(declvar('topc', 'qof', 'nmru', nmru, 'real', 
     + 'Predicted overland flow','meters',QOF)
     +   .ne.0) return

      if(declvar('topc', 'qof_basin', 'one', 1, 'real', 
     + 'Area-weighted predicted overland flow for the basin',
     + 'meters',qof_basin)
     +   .ne.0) return

      ALLOCATE (QOFS(Nmru))
      if(declvar('topc', 'qofs', 'nmru', nmru, 'real', 
     + 'Predicted saturation excess overland flow',
     + 'meters',QOFS).ne.0) return

      if(declvar('topc', 'qofs_basin', 'one', 1, 'real', 
     + 'Area-weighted predicted saturated overland flow '//
     + 'for the basin', 'meters',qofs_basin)
     +   .ne.0) return

      ALLOCATE (QUZ(Nmru))
      if(declvar('topc', 'quz', 'nmru', nmru, 'real', 
     + 'Predicted flow from unsaturated zone to saturated zone',
     + 'meters',QUZ).ne.0) return

      ALLOCATE (REX(Nmru))
      if(declvar('topc', 'rex', 'nmru', nmru, 'real', 
     + 'Predicted infiltration excess surface runoff.',
     + 'meters',REX).ne.0) return
     
      if(declvar('topc', 'rex_basin', 'one', 1, 'real', 
     + 'Area-weighted predicted infiltration excess for the basin',
     + 'meters',rex_basin)
     +   .ne.0) return
c
c declare parameters
c
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
c
c Reintroduce the percent vertical macropore flow, pmacro.
c pmacro is the percentage of potential infiltration after
c surface control (Hortonian) that bypasses the root zone.
c A fraction of pmacro, pmac_sat, completely bypasses the entire
c unsaturated zone such that solutes deposited in precip and/or
c irrigation mix immediately with saturated zone waters.
c pmac_sat is also a fraction of potential infiltration and
c has an upper limit of pmacro.
c
      ALLOCATE (pmacro(Nmru))
      if(declparam('topc','pmacro', 'nmru',
     +   'real', '0.0', '0.0', '1.0',
     +   'Fraction of potential infiltration that bypasses '//
     $   'the root zone', 'Fraction of potential infiltration '//
     $   'that bypasses the root zone',
     +   'fraction') .ne.0) return

      ALLOCATE (pmac_sat(Nmru))
      if(declparam('topc','pmac_sat', 'nmru',
     +   'real', '0.0', '0.0', '1.0',
     +   'Fraction of potential infiltration that bypasses '//
     $   'the root zone and mixes immediately with the saturated zone',
     $   'Fraction of potential infiltration '//
     $   'that bypasses the root zone and mixes immediately '//
     $   'with the saturated zone',
     +   'fraction') .ne.0) return

      if(declparam('topc', 'qobsta', 'one', 'integer',
     +   '1', 'bounded', 'nobs',
     +   'Index of streamflow station for calculating '//
     +   'objective function.','Index of streamflow station '//
     +   'for calculating objective function.','none').ne.0) return

c
c The following 10 parameters, enable volume calculations and soil moisture
c calculations needed to track water table heights, reservoir volumes and 
c solute fluxes (to be computed by the phreeq_mms module).
c
      ALLOCATE (s_theta_fc(Nmru))
      if(declparam('topc','s_theta_fc', 'nmru',
     +   'real', '0.23', '0.01', '0.7',
     +   'Volumetric soil moisture content at field capacity',
     +   'Volumetric soil moisture content at field capacity. '//
     +   'Field capacity is determined as the moisture content '//
     +   'at which the hydraulic conductivity is equal to '//
     +   '1E-8 cm/s.', 'cm3/cm3') .ne.0) return

      ALLOCATE (s_theta_0(Nmru))
      if(declparam('topc','s_theta_0', 'nmru',
     +   'real', '0.23', '0.01', '0.7',
     +   'Initial volumetric soil moisture content in the root zone.',
     +   'Initial volumetric soil moisture content in the root zone.',
     +   'cm3/cm3') .ne.0) return

      ALLOCATE (s_theta_wp(Nmru))
      if(declparam('topc','s_theta_wp', 'nmru',
     +   'real', '0.13', '0.01', '0.56',
     +   'Volumetric soil moisture content at wilting point',
     +   'Volumetric soil moisture content at wilting point. '//
     +   'The wilting point is determined as the mositure content '//
     +   'at a tension of 15,300 cm (15 bars). Also known as '//
     +   'residual soil moisture content.',
     +   'cm3/cm3') .ne.0) return

      ALLOCATE (s_porosity(Nmru))
      if(declparam('topc','s_porosity', 'nmru',
     +   'real', '0.4', '0.1', '0.8',
     +   'Soil porosity', 'Effective soil porosity, equal '//
     +   'to saturated soil moisture content.',
     +   'cm3/cm3') .ne.0) return

c
c  Make srmax a variable to be calculated as
c  (s_theta_fc - s_theta_wp) * s_root_depth
c
c$$$      if(decl*param('topc', 'srmax', 'nmru', 'real',
c$$$     +   '1', '0', '5',
c$$$     +   'Available water capacity of root zone.',
c$$$     +   'Available water capacity of root zone. '//
c$$$     +   'Equal to the difference between field capacity '//
c$$$     +   'and the wilting point','m') .ne.0) return
c$$$
      ALLOCATE (gw_loss_k(Nmru))
      if(declparam('topc', 'gw_loss_k', 'nmru', 'real',
     +   '10.0', '-1.0', '11.0',
     +   'Negative common log of hydraulic conductivity for '//
     $   'bedrock or aquitard.',
     +   'Negative common log of hydraulic conductivity for '//
     $   'bedrock or aquitard. The rate that water is '//
     $   'irreversably lost from the saturated zone '//
     $   'to a deep aquifer is controlled by this parameter and '//
     $   'the variable z_wt, the average height above bedrock. '//
     $   'Any value of 11 or greater results in zero loss.',
     +   '-log10(cm/s)') .ne.0) return

      ALLOCATE (s_rock_depth(Nmru))
      if(declparam('topc', 's_rock_depth', 'nmru', 'real',
     +     '6.0', '0.1', '100',
     +     'Average depth to bedrock for the MRU.','Average depth to '//
     $     'bedrock. Must be greater than the rooting depth, '//
     $     's_root_depth.',
     +     'm') .ne.0) return

      ALLOCATE (s_root_depth(Nmru))
      if(declparam('topc', 's_root_depth', 'nmru', 'real',
     +   '1.8', '0.1', '50',
     +   'Rooting depth.','Rooting depth from ground surface, '//
     +   'Available water capacity (moisture content at field '//
     +   'capacity minus that at wilting point) * root_depth '//
     +   'equals the maximum soil moisture deficit, srmax: '//
     +   'smcont_sc = theta_wp+((srmax - srz)/root_depth).'//
     +   'Be sure to set root_depth to a value less then '//
     +   'the depth to bedrock, s_rock_depth',
     +   'm') .ne.0) return

      ALLOCATE (s_satpref_zmin(Nmru))
      if(declparam('topc', 's_satpref_zmin', 'nmru', 'real',
     +   '-5.0', '-100.0', '0.0',
     +   'Water table elevation at which preferential flow in the '//
     +   'saturated zone begins.','Water table elevation (z=0 at '//
     $   'surface) at which preferential flow in the saturated '//
     $   'zone begins. Set equal to s_satpref_zmax if path does '//
     $   'not exist.', 'm') .ne.0) return

      ALLOCATE (s_satpref_zmax(Nmru))
      if(declparam('topc', 's_satpref_zmax', 'nmru', 'real',
     +   '-5.0', '-100.0', '0.0',
     +   'Water table elevation at which preferential flow in the '//
     +   'saturated zone reaches a maximum.','Water table elevation '//
     +   '(z=0 at surface) at which preferential flow in the '//
     $   'saturated zone reaches a maximum. Set equal to '//
     $   's_satpref_zmin if path does not exist.', 'm') .ne.0) return

      ALLOCATE (s_satpref_k(Nmru))
      if(declparam('topc', 's_satpref_k', 'nmru', 'real',
     +   '0.0001', '0.0', '10.0',
     +   'Hydraulic conductivity for preferential flow through the '//
     +   'saturated zone','Hydraulic conductivity for preferential '//
     +   'flow through the saturated zone. The maximum discharge '//
     +   'qpref_max, will be proportional to (zmax - zmin)*satpref_k.',
     +   'cm/s') .ne.0) return

      ALLOCATE (t0(Nmru))
      if(declparam('topc', 't0', 'nmru', 'real',
     +   '-2', '-6', '4',
     +   'Mean subcatchment value of ln(T0)',
     +   'Mean subcatchment value of ln(T0) '//
     +   'where local transmissivity T=T0*exp(T0-TL)',
     +   'ln(m^2/h)').ne.0) return

      ALLOCATE (szm(Nmru))
      if(declparam('topc', 'szm', 'nmru', 'real',
     +   '.03', '0', '10',
     +   'Value of M in recession equation.',
     +   'Value of M in recession equation.'//
     +   'QB = SZQ*EXP(SBAR/M)','m'
     +    ).ne.0) return

      ALLOCATE (xk0(Nmru))
      if(declparam('topc', 'xk0', 'nmru', 'real',
     +   '.2', '.01', '5',
     +   'Median value of subcatchment surface vertical '//
     +   'hydraulic conductivity.',
     +   'Median value of subcatchment surface vertical '//
     +   'hydraulic conductivity. This is allowed to vary '//
     +   'from its theoretical link with T0 to allow for '//
     +   'anisotropic soil conductivities/transmissivities.',
     +   'm/hr').ne.0) return

      ALLOCATE (xk_cv(Nmru))
      if(declparam('topc', 'xk_cv', 'nmru', 'real',
     +   '14.13', '0.0', '100.0',
     +   'Coefficient of variation for population ln(XK0).',
     +   'Coefficient of variation for population ln(XK0). '//
     +   'A value of 0.0 indicates homogeneous soils for '//
     +   'the entire subcatchment, resulting in the original '//
     +   'behavior of TMOD9502. A value of 14.13 approximates '//
     +   'a lognormal distribution with an order of magnitude '//
     +   'change in XK0 for each standard deviation unit in '//
     +   'lognormal space.', 'unitless(standard deviation/mean)')
     +   .ne.0) return

      ALLOCATE (td(Nmru))
      if(declparam('topc', 'td', 'nmru', 'real',
     +   '60.', '0.', '240',
     +   'Unsaturated zone time delay per unit of storage deficit.',
     +   'Unsaturated zone time delay per unit of storage deficit'//
     +   'in QUZ=SUZ/(SD*TD).',' h/m').ne.0) return

      ALLOCATE (qdffrac(Nmru))
      if(declparam('topc', 'qdffrac', 'nmru', 'real',
     +   '.3', '0', '1',
     +   'Proportion of unsaturated zone drainage that runs off'//
     +   ' as direct flow.','Fraction of unsaturated zone drainage'//
     +   ' that runs off as direct flow.'//
     +   'QDF=QDFFRAC*QUZ','Proportion')
     +    .ne.0) return

      ALLOCATE (sbar0(Nmru))
      if(declparam('topc', 'sbar0', 'nmru', 'real',
     +   '.001', '0', '10',  
     +   'Initial soil moisture deficit, SBAR, in MRU.',
     +   'Initial soil moisture deficit, SBAR, in MRU.',
     +   'm').ne.0) return
            
      if(declparam('topc', 'infex', 'one', 'integer',
     +   '0', '0', '1',  
     +   'Switch for infiltration excess computation (1=y, 0=n).',
     +   'Switch for infiltration excess computation (1=y, 0=n).',
     +   'none').ne.0) return

      ALLOCATE (hf(Nmru))
      if(declparam('topc', 'hf', 'nmru', 'real',
     +   '.01', '0', '1',
     +   'Wetting front suction.',
     +   'Wetting front suction.','m')
     +    .ne.0) return

      ALLOCATE (dth(Nmru))
      if(declparam('topc', 'dth', 'nmru', 'real',
     +   '0.35', '0', '0.7',
     +   'Water content change across the wetting front.',
     +   'Water content change across the wetting front',
     +   'vol/vol').ne.0) return

      ALLOCATE (nacsc(Nmru))
      if(declparam('topc', 'nacsc', 'nmru', 'integer',
     +   '1', '0', '100',
     +   'Number of ln(a/tanB) increments in the subcatchment.',
     +   'Number of ln(a/tanB) increments in the subcatchment.',
     +   'none').ne.0) return

      ALLOCATE (ac(nac,Nmru))
      if(declparam('topc', 'ac', 'nac,nmru', 'real',
     +   '1', '0', '1',
     +   'Fractional area for each ln(a/tanB) increment.',
     +   'Fractional area for each ln(a/tanB) increment.',
     +   'km2/km2').ne.0) return      

      ALLOCATE (st(nac,Nmru))
      if(declparam('topc', 'st', 'nac,nmru', 'real',
     +   '5', '0', '25',
     +   'The ln(a/tanB) value.',
     +   'The ln(a/tanB) value.',
     +   'none').ne.0) return

      allocate(riparian_thresh(nmru))
      if(declparam('topc', 'riparian_thresh', 'nmru',
     +   'real','10.0', '1.0', '40.0',
     +   'UZ bins wetter than this are riparian',
     +   'UZ bins wetter than this are riparian; '//
     +   'UZ bins drier than this are upslope.',
     +   'lambda').ne.0) return

      ALLOCATE (AREA(Nmru)) ! local copy for mru_area_frac
      if(declparam('topc', 'mru_area_frac', 'nmru', 'real',
     +   '1', '0', '1',
     +   'Subcatchment area/total area',
     +   'Subcatchment area/total area',
     +   'none').ne.0) return

      if(declparam('topc', 'dtinit', 'one', 'real',
     +   '24', '0', '24',
     +   'Initial timestep for initialize function.',
     +   'Initial timestep for initialize function.',
     +   'hours').ne.0) return

      ALLOCATE (resp_hr_full(Nmru))
      if(declparam('topc', 'resp_hr_full', 'nmru', 'real',
     +   '6', '0', '100',
     +   'Time step after which the saturation readjustment '//
     +   'coefficient becomes 1','Time step after which the '//
     +   'saturation radjustment coefficient become 1 (i.e. '//
     +   'original TOPMODEL)','hours').ne.0) return

      ALLOCATE (resp_hr_min(Nmru))
      if(declparam('topc', 'resp_hr_min', 'nmru', 'real',
     +   '.083', '0', '24',
     +   'Time step below which the saturation readjustment '//
     +   'coefficient is set equal to resp_coef_min',
     +   'Time step below which the saturation readjustment '//
     +   'coefficient is set equal to resp_coef_min',
     +   'hours').ne.0) return

      ALLOCATE (resp_coef_min(Nmru))
      if(declparam('topc', 'resp_coef_min', 'nmru', 'real',
     +   '.01', '0', '1',
     +   'Minimum saturation readjustment coefficient',
     +   'Minimum saturation readjustment coefficient',
     +   'none').ne.0) return

      if(declparam('topc', 'iout', 'one', 'integer',
     +   '1', '0', '2',
     +   'Printed output switch. 0=min,1=summary,2=each timestep',
     +   'Printed output switch. 0=min,1=summary,2=each timestep',
     +   'none').ne.0) return

      ALLOCATE (sc_area(Nmru)) ! local copy of mru_area
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

      if(declparam('io', 'topout_file_unit', 'one', 'integer',
     +   '80', '50', '99',
     +   'Unit number for TOPMODEL output file',
     +   'Unit number for TOPMODEL output file',
     +   'integer').ne.0) return

!
! From init section
       ALLOCATE (atag(nac))
       ALLOCATE (riparian(nac,nmru))
! areas for riparian(1), and upland(2) TWI bins. all UZ areas (irip=0) uses the mru_area
       ALLOCATE (uz_area(2,nmru))
!
! From run section
       ALLOCATE (last_uz_dep(nac,Nmru))
       ALLOCATE (last_srz(nac,Nmru))
       ALLOCATE (last_suz(nac,Nmru))
       ALLOCATE (last_z_wt_local(nac,Nmru))
       ALLOCATE (SUMP(nmru))
       ALLOCATE (SUMAE(nmru))
       ALLOCATE (SUMQ(nmru))
       ALLOCATE (ACMAX(nmru))
       ALLOCATE (EX(nac))
       ALLOCATE (ihour(nac,nmru))

! Getvars from other modules
!

       ALLOCATE (transp_on(Nmru))
       ALLOCATE (snowmelt(Nmru))
       ALLOCATE (psoilmru(Nmru))
       ALLOCATE (potet(Nmru))
       ALLOCATE (QOBS(Nobs))
       ALLOCATE (irrig_int_next(nirrig_int))
       ALLOCATE (gw_ext(Ngw_ext))
       ALLOCATE (irrig_sat_mru(Nmru))
       ALLOCATE (intcp_evap(Nmru))
       ALLOCATE (snow_evap(Nmru))
       ALLOCATE (sbar_max(Nmru))

      topmdecl = 0

      return
      end

c***********************************************************************
c
c     topcinit - Initialize topcomp module - get parameter values,
c

      integer function topminit()

      USE WEBMOD_TOPMOD
      
      logical acflag
      integer is, j
      real dtinit

C***  local variables

      real tarea, SUMAC, T0DT, E_XK, sigma

      integer IA, ik

c      character*135 output_path
c      logical filflg

      topminit = 1
      step1 = .true.


c----- set name for topmod unique output file 
c      ret = getoutname (output_path, '.topout')
c      inquire(file=output_path,exist=filflg)
c      if (filflg) then
c        open(unit=80,file=output_path,status='old')
c        close(unit=80,status='delete')
c      endif

c-----open the file.
c      open (unit=80,file=output_path,access='sequential',
c     * form='formatted', status='new')

      if(getparam('topc', 'dtinit', 1 , 'real', dtinit)
     +   .ne.0) return

      if(getparam('topc', 'sched_gw1', nmru , 'integer', sched_gw1)
     +   .ne.0) return
           
      if(getparam('topc', 'sched_gw2', nmru , 'integer', sched_gw2)
     +   .ne.0) return
           
      if(getparam('topc', 'gwbnd_len1', nmru , 'real', gwbnd_len1)
     +   .ne.0) return
           
      if(getparam('topc', 'gwbnd_len2', nmru , 'real', gwbnd_len2)
     +   .ne.0) return
           
      DT = dtinit
      
      if(getparam('intcp', 'cov_type', nmru, 'integer', cov_type)
     +   .ne.0) return

      if(getparam('intcp', 'covden_sum', nmru, 'real', covden_sum)
     +   .ne.0) return 

      if(getparam('intcp', 'covden_win', nmru, 'real', covden_win)
     +   .ne.0) return 

      if(getparam('topc', 'iout', 1 , 'integer', iout)
     +   .ne.0) return

      if(getparam('topc', 't0', nmru, 'real', T0)
     +   .ne.0) return

      if(getparam('topc', 'szm', nmru, 'real', SZM)
     +   .ne.0) return

      if(getparam('topc', 'td', nmru, 'real', TD)
     +   .ne.0) return
c
c  srmax to be made into a variable derived from field capacity,
c  s_theta_fc, wilting point, s_theta_wp, and rooting depth, s_root_depth
c
c$$$      if(get*param('topc', 'srmax', nmru, 'real', SRMAX)
c$$$     +   .ne.0) return
c$$$

      
      if(getparam('topc', 's_theta_wp', nmru, 'real', s_theta_wp)
     +   .ne.0) return

      if(getparam('topc', 's_theta_0', nmru, 'real', s_theta_0)
     +   .ne.0) return

      if(getparam('topc', 's_theta_fc', nmru, 'real', s_theta_fc)
     +   .ne.0) return

      if(getparam('topc', 's_porosity', nmru, 'real', s_porosity)
     +   .ne.0) return

      if(getparam('topc', 's_root_depth', nmru, 'real', s_root_depth)
     +   .ne.0) return

      if(getparam('topc', 's_rock_depth', nmru, 'real', s_rock_depth)
     +   .ne.0) return

      if(getparam('topc', 's_satpref_zmin', nmru, 'real',
     $     s_satpref_zmin) .ne.0) return

      if(getparam('topc', 's_satpref_zmax', nmru, 'real',
     $     s_satpref_zmax) .ne.0) return

      if(getparam('topc', 's_satpref_k', nmru, 'real',
     $     s_satpref_k) .ne.0) return

      if(getparam('topc', 'gw_loss_k', nmru, 'real',
     $     gw_loss_k) .ne.0) return

      if(getparam('topc', 'qdffrac', nmru, 'real', qdffrac)
     +   .ne.0) return

      if(getparam('topc', 'resp_hr_full', nmru, 'real', resp_hr_full)
     +   .ne.0) return

      if(getparam('topc', 'resp_hr_min', nmru, 'real', resp_hr_min)
     +   .ne.0) return

      if(getparam('topc', 'resp_coef_min', nmru, 'real',
     +   resp_coef_min).ne.0) return
c
c sr0 is now calculated from initial soil moisture content, s_theta_0
c
c$$$      if(get*param('topc', 'sr0', nmru, 'real', SR0)
c$$$     +   .ne.0) return
c$$$
      if(getparam('topc', 'sbar0', nmru, 'real', SBAR0)
     +   .ne.0) return

      if(getparam('topc', 'qobsta', 1, 'integer', qobsta)
     +   .ne.0) return

      if(getparam('precip', 'irrig_int_src', nmru, 'integer',
     $     irrig_int_src).ne.0) return

      if(getparam('precip', 'irrig_int_init', nmru, 'real',
     $     irrig_int_init).ne.0) return

      if(getparam('precip', 'irrig_sched_int', nmru, 'integer',
     $     irrig_sched_int).ne.0) return

      if(getparam('topc', 'irrig_int_max', nmru, 'real',
     $     irrig_int_max).ne.0) return

      if(getparam('topc', 'infex', 1, 'real', INFEX)
     +   .ne.0) return

      if(getparam('topc', 'pmacro', nmru, 'real', pmacro)
     +   .ne.0) return

      if(getparam('topc', 'pmac_sat', nmru, 'real', pmac_sat)
     +   .ne.0) return

      if(getparam('topc', 'xk0', nmru, 'real', XK0)
     +   .ne.0) return

      if(getparam('topc', 'xk_cv', nmru, 'real', XK_CV)
     +   .ne.0) return

      if(getparam('topc', 'hf', nmru, 'real', HF)
     +   .ne.0) return
     
      if(getparam('topc', 'dth', nmru, 'real', DTH)
     +   .ne.0) return

      if(getparam('topc', 'mru_area_frac', nmru, 'real', AREA)
     +   .ne.0) return

      if(getparam('basin', 'nacsc', nmru, 'integer', nacsc)
     +   .ne.0) return

      if(getparam('basin', 'ac', nac*nmru, 'real', AC)
     +   .ne.0) return

      if(getparam('basin', 'st', nac*nmru, 'real', ST)
     +   .ne.0) return
     
      if(getparam('basin', 'mru_area', nmru, 'real', sc_area)
     +   .ne.0) return

      if(getparam('basin', 'basin_area', 1 , 'real', basin_area)
     +   .ne.0) return

      if(getparam('io', 'topout_file_unit', 1, 'integer',
     +   topout_file_unit).ne.0) return

c$$$      write(topout_file_unit, 8000) is
c$$$ 8000 format(1x,'Subcatchment  ', i3)

      do 50 is = 1, nsc
c
c Test that both pmacro is >= pmac_sat and that
c both are <= 1.
c     
         if(pmacro(is).gt.1.0.or.pmacro(is).lt.0.0.or..not.
     $        pmacro(is).ge.pmac_sat(is).and.pmac_sat(is).ge.0.0) then
          write(*,*)'ERROR: PMACRO must be greater than PMAC_SAT and '//
     $         'both must be between zero and one ',
     $         'Check values for MRU ',is
            return
         end if


c$$$      write(topout_file_unit, 8000) is
c$$$ 8000 format(1x,'Subcatchment  ', i3)

*  NAC IS NUMBER OF A/TANB ORDINATES
c  nsc is number of subcatchments
*  AREA IS SUBCATCHMENT AREA AS PROPORTION OF TOTAL CATCHMENT 
c      READ(8,*)(AC(J),ST(J),J=1,NAC)
*  AC IS DISTRIBUTION OF AREA WITH LN(A/TANB)
*  ST IS LN(A/TANB) VALUE

       tarea = ac(1,is)
      if(tarea.ne.0.0)then
         write(*,*)'ERROR: The first area for the first ln(a/tanB) ',
     $        'for mru ', is,'must equal zero. It is currently set to ',
     $         tarea, '. Please correct and rerun.'
         return
      endif
! Identify first (and wettest) TWI bin as riparian 
      riparian(1,is) = .true.
      uz_area(1,is) = 0.5*(AC(1,is)+AC(2,is))*sc_area(is)! uz_riparian
      uz_area(2,is) = 0.0! uz_uplands
c
c Use average ac (acf calculation) to check for multiple zero areas and warn
c user about using those indices.

       acflag = .FALSE.

       do 10 j=2,nacsc(is)
        riparian(j,is)=.true.
        tarea = tarea + ac(j,is)
        if(st(j,is).lt.riparian_thresh(is)) riparian(j,is)=.false. ! identify upland TWI bins based on fixed threshold
        if((ac(j-1,is)+ac(j,is))/2.eq.0) then  ! getting past check above insures that ac(1,is) = 0
          atag(j-1) = 1
          acflag = .true.
        else
          atag(j-1) = 0
        end if
!  Assign uz composite areas for load conversions
        if(j.eq.nacsc(is)) then
          ACF=0.5*AC(j,is)
        else
          ACF=0.5*(AC(j,is)+AC(j+1,is))
        endif
        ! uz_all uses the mru_area
        if(riparian(j,is)) then
         uz_area(1,is) = uz_area(1,is)+ACF*sc_area(is)! uz_riparian
        else
         uz_area(2,is) = uz_area(2,is)+ACF*sc_area(is)! uz_uplands
        endif
   10  continue

      if(abs(tarea-1.0).gt.0.0005)then
         write(*,*)'ERROR: The fractional areas for ln(a/tanB) ',
     $        'do not sum to one. Check ac values for MRU ',is
         return
      endif
c
c  Flag topographic indices if they exist
c
      if(acflag) then
         print*,'Because of zero average areas for consecutive ',
     $        'topographic indices, Reference to root zone or ',
     $        'unsaturated zone variable with the following ',
     $        'indices (nac,nmru)would be nonsensical.'
         do 6 j=1,nacsc(is)
            if (atag(j).eq.1) print*,j,is
 6       continue
         print*,' '
      end if
*
*  CALCULATE AREAL INTEGRAL OF LN(A/TANB)
*  NB.  a/tanB values should be ordered from high to low with ST(1)
*  as an upper limit such that AC(1) should be zero, with AC(2) representing
*  the area between ST(1) and ST(2)
      TL(is)=0.
      SUMAC=AC(1,is)
      DO 11 J=2,nacsc(is)
        SUMAC=SUMAC+AC(J,is)
        TL(is)=TL(is)+AC(J,is)*(ST(J,is)+ST(J-1,is))/2
   11 CONTINUE

c$$$      write(topout_file_unit,8010) TL(is), SUMAC
c$$$ 8010 format(1x, 'TL = ', f8.2,/, 1x, 'SUMAC = ', f8.2)
c$$$
c      AC(NAC+1,is)=0.
*
*
c      If(IOUT.ge.1)Write(10,600)TL, SUMAC
c  600 Format(1x,'TL = ',f8.2,/'SUMAC = ', f8.2)



*
*  READ PARAMETER DATA
c      READ(9,"(A)")SUBCAT
c      READ(9,*)SZM,T0,TD,CHV,RV,SRMAX,Q0,SR0,INFEX,XK0,HF,DTH
*
*  Convert parameters to m/time step DT
*  with exception of XK0 which must stay in m/h
*                    Q0 is already in m/time step
*                    T0 is input as Ln(To)
c
c
c  Changed Q0 to from 'm/timestep' to m^3/s so conversion will
c  be needed. Also changed Q0 from dimension nmru to one - RW
c
c  Added a linear ramp function to scale the rate that the
c  saturation deficit will respond to the new SBAR at the
c  beginning of each time step. The saturation deficit
c  for each ln(a/tanB) bin was calculated as
c
c  SD(IA,is)=SBAR(is)+SZM(is)*(TL(is)-ST(IA,is))
c
c  will now be 
c
c  SD(IA,is)=SBAR(is)+(SZM(is)*(TL(is)-ST(IA,is)))*resp_coef
c
c  where resp_coef will ramp up from resp_coef_min to 1 in the
c  interval resp_hr_min < DT < resp_hr_max. As long as the time
c  step remains constant, so will resp_coef.
c
c  This allows for more realistic collection and delivery of
c  apparent exfiltrated water (i.e. SD<0) over various time-scales.
c
c  15 sep 2003 remove loop here since we are already inside the
c              DO 50 1, nsc loop
c
c$$$      do 12 j=1,nsc
c$$$      if(dt.le.resp_hr_min(j)) then
c$$$         resp_coef(j) = resp_coef_min(j)
c$$$      else if (dt.ge.resp_hr_full(j)) then
c$$$         resp_coef(j) = 1.0
c$$$      else
c$$$         resp_coef(j) = resp_coef_min(j) + ((1.0 - resp_coef_min(j))*
c$$$     +      (dt - resp_hr_min(j))/(resp_hr_full(j) - resp_hr_min(j)))
c$$$      end if
c$$$ 12   continue
c      do 12 j=1,nsc
      if(dt.le.resp_hr_min(is)) then
         resp_coef(is) = resp_coef_min(is)
      else if (dt.ge.resp_hr_full(is)) then
         resp_coef(is) = 1.0
      else
         resp_coef(is) = resp_coef_min(is) + ((1.0 - resp_coef_min(is))*
     $      (dt - resp_hr_min(is))/(resp_hr_full(is) - resp_hr_min(is)))
      end if
c 12   continue
c
      T0DT = T0(is) + ALOG(DT)

*  Calculate SZQ parameter - maximum baseflow discharge for each
*  subcatchment for each time step.
      SZQ(is) = EXP(T0DT-TL(is))
**
c      WRITE(80,604)TCH(NCH),(AR(IR,is),IR=1,NR)
c  604 FORMAT(1X,'SUBCATCHMENT ROUTING DATA'/
c     1  1X,'Maximum Routing Delay  ',E12.5/
c     2  1X,'Histogram ordinates  ',/(1X,5E12.5))
c
c Convert moisture contents into maximum root zone deficit (srmax)
c
      srmax(is) = (s_theta_fc(is) - s_theta_wp(is)) * s_root_depth(is)
c
c     Initialize lognormal distribution of vertical hydraulic
c     conductivity using the median value, XK0, and the
c     coefficient of variation, XK_CV. The XK array holds the
c     back-transformed z values from -4 to +4 standard
c     deviations. If the coefficient of variation, xk_cv,
c     is set equal to zero, then a homogeneous soil is assumed
c     with a vertical conductivity equal to XK0.   - RMTW
c
         if (xk_cv(is).ne.0.0) then
            e_xk = log(xk0(is))
            sigma = sqrt(log(xk_cv(is)**2+1))
            nxk(is) = nxkbin
            do 55 ik = 1,nxkbin
               xk(is,ik) = exp(e_xk+(ik-5)*sigma)
 55         continue
         else
            xk(is,1)=xk0(is)
            nxk(is) = 1
         endif
         

c
 50   continue
c      write(80,605) sumar
c  605 format(1X,'Sum of histogram ordinates for'//
c     1  ' all subcatchments  ', f10.4)
c
c Place infex loop below inside of 50 loop above
c
cc$$$      do 56 is = 1, nsc
c$$$         if (xk_cv(is).ne.0.0) then
c$$$            e_xk = log(xk0(is))
c$$$            sigma = sqrt(log(xk_cv(is)**2+1))
c$$$            nxk(is) = nxkbin
c$$$            do 55 ik = 1,nxkbin
c$$$               xk(is,ik) = exp(e_xk+(ik-5)*sigma)
c$$$ 55         continue
c$$$         else
c$$$            xk(is,1)=xk0(is)
c$$$            nxk(is) = 1
c$$$         endif
c$$$ 56   continue

      If(INFEX.ne.0) then
         Write(topout_file_unit,8020)(ack(ik),ik=1,9)
 8020    Format('TOPMODEL output file'//
     $       'Infiltration excess parameters'/
     +       'Proportional areas: '/3X,9(1x,e12.5)/
     +       'MRU Vertical Hydraulic Conductivities, XK'/
     +       '=== ',9('============ '))
         do 57 is = 1,nsc
            Write(topout_file_unit,8021)is,(xk(is,ik),ik=1,nxk(is))
 8021       Format(I3,9(1x,e12.5))
 57      continue
      else
         Write(topout_file_unit,8025)
 8025    Format('TOPMODEL output file'//
     $        'Infiltration excess not simulated')
      endif
*
*  INITIALISE SRZ AND Q0 VALUES HERE
*  SR0 IS INITIAL ROOT ZONE STORAGE DEFICIT BELOW FIELD CAPACITY
*


c
c  Write a subcatchment characteristics and initial water
c  balance for each MRU
c
      If(IOUT.ge.1) write(topout_file_unit, 8000)
 8000 format(//'Characteristics of subcatchments'/
     $     'MRU: Subcatchment ID'/
     $     'TL: Mean Topographic Index'/
     $     'SZM: Shape factor for exponential transmissivity'/
     $     'T0: Naperian log of maximum transmissivity'/
     $     'SZQ: Maximum baseflow discharge in m/hr'/
     $     'PMACRO: Fraction of infiltration that bypasses root zone'/
     $     'PMAC_SAT: Fraction of infiltration that bypasses'/
     $     '  the root zone and the unsaturated zone (Max=PMACRO)'/
     $     'QDFFRAC: Fraction of recharge diverted to stream'/
     $     'S_TH_WP: Soil moisture at wilting point(cm3/cm3)'/
     $     'S_TH_0: Initial soil moisture(cm3/cm3)'/
     $     'S_TH_FC: Soil moisture at field capacity(cm3/cm3)'/
     $     'S_POROSITY: Soil porosity(cm3/cm3)'/
     $     'S_ROOT_DEP: Rooting Depth (m)'/
     $     'S_ROCK_DEP: Depth to bedrock (m)'/
     $     'GW_LOSS_K: Conductivity of bedrock or '/
     $     '  aquitard (-log(cm/s))'/
     $     'S_SATP_ZMN: Minimum water table elevation, '/
     $     '  in meters, at which the deep preferential flow '/
     $     '  paths (i.e. tile drains) become active.'/
     $     'S_SATP_ZMX: Water table elevation, in meters, '/
     $     '  at which the discharge from preferential flow in '/
     $     '  the saturated zone will be at its maximum discharge.'/
     $     'S_SATP_K: Conductivity, in cm/s, of preferential '/
     $     '  flow in the saturated zone.'/
     $     'QPREF_MAX: Maximum discharge from preferential flow '/
     $     '  in the saturated zone (m/hr).'/
     $     'SRMAX: Storage between field capacity and wilting point'/
     $     '  Derived as (s_th_fc - s_th_wp)* s_root_dep.'/
     $     'SR0: Initial root zone deficit'/
     $     '  Derived as (s_th_fc - s_th_0) * s_root_dep.'/
     $     'SBAR0: Initial soil moisture deficit, in meters.'/
     $     '  [adjusted to include 1st day well pumping]'/
     $     'Z_0: Initial average water level above bedrock.'/
     $     '  Derived as (s_rock_dep - (sbar0/s_porosity)'/
     $     'BAL: Inital water balance, in meters ='/
     $     '  (S_ROCK_DEP * S_POROSITY) -(SR0+SBAR0)'//
     $     'MRU  TL    SZM        T0         SZQ        PMACRO     ',
     $     'PMAC_SAT   QDFFRAC    S_TH_WP    S_TH_0     S_TH_FC    ',
     $     'S_POROSITY S_ROOT_DEP S_ROCK_DEP GW_LOSS_K  S_SATP_ZMN ',
     $     'S_SATP_ZMX S_SATP_K   QPREF_MAX  SRMAX      SR0        ',
     $     'SBAR0      Z_0        INIT STORE'/
     $     '==== ===== ========== ========== ========== ========== ',
     $     '========== ========== ========== ========== ========== ',
     $     '========== ========== ========== ========== ========== ',
     $     '========== ========== ========== ========== ========== ',
     $     '========== ========== ==========')


*  INITIALISE STORES
*  and print initial conditions - RMTW
*

c Initialize unsaturated zone water content
      do 26 is=1,nsc
         sr0(is) = (s_theta_fc(is) - s_theta_0(is)) * s_root_depth(is)
         
c Readily drainable porosity is the difference getween saturated
c soil moisture content and field capacity

         s_drain(is) = s_porosity(is) - s_theta_fc(is)
c
c The maximum moisture deficit
c
         sbar_max(is) = s_rock_depth(is)*s_drain(is)

c
c If there is an amount of initial irrigation from a shallow well in this
c MRU then increase the amount of available water (decrease SBAR) in the
c saturated zone by that amount.
c
         if(irrig_int_init(is).ne.0.0 .and.
     $        irrig_int_src(is).eq.0)
     $        sbar0(is) = sbar0(is) -
     $        (irrig_int_init(is)*inch2m)

         DO 25 IA=1,nacsc(is)
            SUZ(IA,is)=0.       ! No recharge water in UZ
            SRZ(IA,is)=SR0(is)
C  CALCULATE LOCAL STORAGE DEFICIT
c
c
            SD(IA,is)=SBAR0(is)+(SZM(is)*(TL(is)-ST(IA,is)))
     $           *resp_coef(is)


c
c Add uz_depth to indicate the initial water content in the
c unsaturated zone. uz_depth is determined by the saturation deficit,
c the root zone deficit, and the unsaturated zone storage
c
            z_wt_local(IA,is) = - (sd(ia,is)/s_drain(is))

            if(z_wt_local(ia,is).ge.-s_root_depth(is)) then 
                    ! just account for the root zone deficit plus suz
               uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
     $              s_theta_wp(is)*s_root_depth(is) + suz(ia,is)
            else
                                ! add additional UZ water at field capacity
               uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
     $              s_theta_wp(is)*s_root_depth(is) + suz(ia,is) -
     $              (s_root_depth(is)+z_wt_local(ia,is))*s_theta_fc(is)
            end if
 25      continue
c
c calculate maximum discharge from preferential flows throught the
c unsaturated zone, in meters per hour
c
      qpref_max(is) = (s_satpref_zmax(is) - s_satpref_zmin(is)) *
     $            s_satpref_k(is) * 36.0
c Use soil characteristics to compute maximum root zone deficit,
c in meters
c      srmax(is) = (s_theta_fc(is) - s_theta_wp(is))* s_root_depth(is)

c The initial average water table elevation
      z_wt(is) =  -sbar0(is)/s_drain(is)
      SBAR(is)=sbar0(is)
c      SBAR(is)=-SZM(is)*ALOG(Q0DT/SZQ(is))
      sumqrex(is)=0
      sumqofs(is)=0
      sumqb(is)=0
      sumquz(is)=0
      sumqdf(is)=0
      sumqvpref(is)=0
      sumqpref(is)=0
      sumqexfil(is)=0
      sumqwet(is)=0
      sumgw_loss(is) = 0

      BAL(is)=(s_rock_depth(is) * s_porosity(is)) - SBAR(is) - SR0(is)


      If(IOUT.ge.1) write(topout_file_unit,603)is, TL(is),szm(is),
     $     T0(is), szq(is), pmacro(is), pmac_sat(is),qdffrac(is), 
     $     s_theta_wp(is),s_theta_0(is),
     $     s_theta_fc(is), s_porosity(is), s_root_depth(is),
     $     s_rock_depth(is), gw_loss_k(is),s_satpref_zmin(is),
     $     s_satpref_zmax(is),
     $     s_satpref_k(is), qpref_max(is), srmax(is), SR0(is), 
     $     SBAR(is), z_wt(is), BAL(is)
  603 format(I4,1X,F5.2,22(1X,E10.3))

c$$$      write(topout_file_unit,8040)
c$$$     +                   is, BAL(is),SBAR(is),SR0(is)
c$$$ 8040 Format(1x, 'Subcatchment = ', i4,/
c$$$     +       3x,'Initial Balance BAL ',e12.5/
c$$$     1       3x,'Initial SBAR        ',e12.5/
c$$$     2       3x,'Initial SR0         ',e12.5)



 26   continue
   

*
*  Initialise water balance.  BAL is positive for storage
c

c    move this print header to the run so it appears above the
c    detailed table output.

c$$$      If(IOUT.ge.2)Write(topout_file_unit,8050)
c$$$ 8050 format(1x,'scat  it       p        ep       q(it)       quz',
c$$$     1'      q       sbar       qof       acm       afx')
c$$$
      topminit = 0

      return
      end

c***********************************************************************
c
c     topcrun -
c

      integer function topmrun()

      USE WEBMOD_TOPMOD
      
      implicit none

C***  local variables

      real ACF, UZ, EA, OF, DF
      real RINT, EP
      real sccumf,schf, scdth, scszm, scxk0, sctp
      real SUMRZ, SUMUZ, sd_temp
      integer is, j, ik
      integer IA, IB, nstep
      integer scirof
      integer endper
      logical end_run, end_yr, end_mo, end_dy, end_storm

      topmrun = 1

      if(getvar('io', 'endper', 1, 'integer', endper)
     +   .ne.0) return

      if(getvar('potet', 'transp_on', nmru, 'integer', transp_on)
     +   .ne.0) return

      if(getvar('nwsmelt', 'snowmelt', nmru, 'real', snowmelt)
     +   .ne.0) return

      if(getvar('nwsmelt', 'psoilmru', nmru, 'real', psoilmru)
     +   .ne.0) return

      if(getvar('obs', 'potet', nmru, 'real', potet)
     +   .ne.0) return

      if(getvar('obs', 'runoff', nobs, 'real', QOBS)
     +   .ne.0) return

      if(getvar('obs', 'irrig_int_next', nirrig_int, 'real',
     $     irrig_int_next).ne.0) return

      if(getvar('obs', 'gw_ext', Ngw_ext, 'real',
     $     gw_ext).ne.0) return

      if(getvar('precip', 'basin_ppt', 1, 'real', POBS)
     +   .ne.0) return

      if(getvar('precip', 'irrig_sat_mru', nmru, 'real',
     $     irrig_sat_mru).ne.0) return

      if(getvar('intcp', 'intcp_evap', nmru, 'real',
     +   intcp_evap).ne.0) return

      if(getvar('snow', 'snow_evap', nmru, 'real',
     +   snow_evap).ne.0) return

      DT = deltim()
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
c  Alter the hillslope response coeffient in case of new dt
c
      do 120 j=1,nsc
      if(dt.le.resp_hr_min(j)) then
         resp_coef(j) = resp_coef_min(j)
      else if (dt.ge.resp_hr_full(j)) then
         resp_coef(j) = 1.0
      else
         resp_coef(j) = resp_coef_min(j) + ((1.0 - resp_coef_min(j))*
     +      (dt - resp_hr_min(j))/(resp_hr_full(j) - resp_hr_min(j)))
      end if
 120  continue
c
c     Calculate metric equivalents of runoff and precip to facilitate
c     comparison with the predicted flow paths.
c

      runoffm = qobs(qobsta)*DT/basin_area/9809.63
      precipm = pobs * .0254


      nstep = getstep()

      if(nstep.eq.1) then

        do 12 is=1,nsc

        SUMP(is) = 0.
        sumgw_in1(is) = 0.
        sumgw_in2(is) = 0.
        sumgw_loss(is) = 0.
        SUMAE(is) = 0.
        SUMQ(is) = 0.
        sumqexfil(is) = 0.
        sumqwet(is) = 0.
        acmax(is) = 0.
        afxmax(is)=0.
        do 13 ik=1,nxkbin
        IROF(is,ik) = 0
        CUMF(is,ik) = 0.
        const(is,ik) = 0.
        tp(is,ik) = 0.
 13     continue
 12   continue
c$$$        QTOT = 0.
c$$$        F1 = 0.
c$$$        F2 = 0.
c$$$        SUMQOBS = 0.
c$$$        SSQ = 0. 
c$$$        NTSTEP = 0
        
      end if

C
C  Reset basin time-step averages to zero
C

      basin_actet = 0
      srz_basin = 0
      suz_basin = 0
      smav_basin = 0
      smcont_basin = 0
      sbar_basin = 0
      rex_basin = 0
      afx_basin = 0
      qofs_basin = 0
      acm_basin = 0
      qof_basin = 0
      basin_infil = 0
      basin_surfdep = 0
      qdf_basin = 0
      qvpref_basin = 0
      qpref_basin = 0
      gw_in1_basin = 0
      gw_in2_basin = 0
      qb_basin = 0
      qexfil_basin = 0
      qwet_basin = 0
      qout_basin = 0
      z_wt_basin = 0
      gwl_basin = 0
C
C  Initialise contributing area counts, root zone deficit,
c  unsaturated zone storage, and withdrawels for irrigation

      do 5 is = 1,nsc
         srz_sc(is) = 0
         suz_sc(is) = 0
         smav_sc(is) = 0
         smcont_sc(is) = 0
         irrig_sat_next(is) = 0
         do 6 ia = 1, nacsc(is)
            ihour(ia,is)=0
    6    continue
    5 continue

c
c  START THE LOOP ON SUBCATCHMENTS
C  START LOOP ON SUBCATCHMENTS

      do 100 is = 1,nsc

C
C  START LOOP ON TIME STEPS
c
C      If(IOUT.ge.2)Write(80,101)
C  101 format(1x,'  it       p        ep       q(it)       quz',
C     1'      q       sbar       qof')
C
c 
c  Reset time-step accumulators
c
      sae(is)=0.
      QOF(is)=0.
      QUZ(is)=0.
      REX(is)=0.
      afx(is) = 0.
      acm(is) = 0.
      qexfil(is)=0.
      qwet(is) = 0.
c$$$      z_wt(is) = 0.
C

      EP= (potet(is)- intcp_evap(is) - snow_evap(is)) * .0254

c Since snow_evap is a parameter, the total ET could exceed PET

      if(ep.lt.0.0)ep=0.0


      P= (snowmelt(is) + psoilmru(is)) * .0254
      infil(is) = p
      SUMP(is) = SUMP(is) + P      
C
C  SKIP INFILTRATION EXCESS CALCULATIONS IF INFEX = 0
      IF(INFEX.EQ.1) THEN
C
C****************************************************************
C  INFILTRATION EXCESS CALCULATIONS USING EXPINF ROUTINE BASED ON
C  GREEN-AMPT INFILTRATION IN A SOIL WITH CONDUCTIVITY DECLINING
C  EXPONENTIALLY WITH DEPTH (REF. BEVEN, HSJ, 1984)
C
C  NOTE THAT IF INFILTRATION EXCESS DOES OCCUR IT WILL DO SO OVER
C  THE WHOLE SUBCATCHMENT BECAUSE OF HOMOGENEOUS SOIL ASSUMPTION
C
C  ALL PARAMETERS AND VARIABLES ON INPUT MUST BE IN M/H
C
C  THIS SECTION CAN BE OMITTED WITHOUT PROBLEM
c
c     Modified Mar 2003 to include lognormal distribution of XK0
c     values evaluated for nine bins. - RMTW
c
c
C***************************************************************
         DO 25 ik = 1,nxk(is)

            IF(P.GT.0.)THEN
C
C  Adjust Rainfall rate from m/time step to m/h
             RINT = P/DT
             scirof = IROF(is,ik)
             sccumf = CUMF(is,ik)
             scxk0 = XK(is,ik)
             schf = HF(is)
             scdth = DTH(is)
             scszm = SZM(is)
             scconst = const(is,ik)
             sctp = tp(is,ik)
             CALL EXPINF(scirof,nstep,RINT,DF,sccumf,scxk0,
     +                         schf,scdth,scszm, DT, scconst, sctp)
C   DF is volumetric increment of infiltration and is returned in m/DT
c
c   Eliminate precision artifacts with next statement
               if(abs(p-df).lt.1e-8) then
                  df=p
                  scirof = 0
                  sccumf = 0
               endif
               xkarea = ack(ik)
               if(nxk(is).eq.1) xkarea = 1.0
             REX(is) = REX(is)+(P - DF)*xkarea
c$$$             P= P - REX(is)
             IROF(is,ik) = scirof
             CUMF(is,ik) = sccumf
               CONST(is,ik) = scconst
               tp(is,ik) = sctp
             
               If(IROF(is,ik).EQ.1)afx(is) = afx(is) + xkarea
            ELSE
             REX(is)=0.
             IROF(is,ik)=0
             CUMF(is,ik)=0.
               CONST(is,ik) = 0.0
               tp(is,ik) = 0.0
            ENDIF
 25      CONTINUE
         P = P-REX(is)
      ENDIF
C****************************************************************
C
C P IS RAINFALL AVAILABLE FOR INFILTRATION AFTER SURFACE CONTROL
C   CALCULATION
c
c   Also remove from P the amount of water that bypasses both the
c   root zone and the unsaturated zone storage, qvpref - RMTW
c
      qvpref(is) = p * pmac_sat(is)
      p = p - qvpref(is)
C
c
C  START LOOP ON A/TANB INCREMENTS
c
      DO 30 IA=1,nacsc(is)
c
      if(ia.eq.nacsc(is)) then
        ACF=0.5*AC(IA,is)
      else
        ACF=0.5*(AC(IA,is)+AC(IA+1,is))
      endif
      UZ=0.
      EX(IA)=0.
      srzwet(ia,is)=0
      if(acf.ne.0.0) then
         last_z_wt_local(ia,is)=z_wt_local(ia,is)
         last_uz_dep(ia,is) = uz_depth(ia,is)
         last_srz(ia,is) = srz(ia,is)
         last_suz(ia,is) = suz(ia,is)
C
C  CALCULATE LOCAL STORAGE DEFICIT and associated fluxes
c
      SD(IA,is)=SBAR(is)+(SZM(is)*(TL(is)-ST(IA,is)))*resp_coef(is)
c      SD(IA,is)=SBAR(is)+SZM(is)*(TL(is)-ST(IA,is))
c
cccccccMove to end so that depth reflect ending water table after ET
c$$$c Add uz_depth to indicate the initial water content in the
c$$$c unsaturated zone. uz_depth is determined by the saturation deficit,
c$$$c the root zone deficit, and the unsaturated zone storage
c$$$c
c$$$      z_wt_local(IA,is) = - (sd(ia,is)/s_drain(is))
c$$$      z_wt(is) = z_wt(is) + z_wt_local(ia,is) * acf
c$$$
c$$$      if(z_wt_local(ia,is).ge.-s_root_depth(is)) then 
c$$$            ! just account for the root zone deficit plus suz
c$$$         uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
c$$$     $        s_theta_wp(is)*s_root_depth(is) + suz(ia,is)
c$$$      else
c$$$            ! add additional UZ water at field capacity
c$$$         uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
c$$$     $        s_theta_wp(is)*s_root_depth(is) + suz(ia,is) -
c$$$     $        (s_root_depth(is)+z_wt_local(ia,is))*s_theta_fc(is)
c$$$      end if
c
c
      IF(SD(IA,is).LT.0.) THEN
        if(abs(sd(ia,is)).ge.srz(ia,is)) then
           srzwet(ia,is) = srz(ia,is)
           srz(ia,is) = 0
           qexfil(is) = qexfil(is) - ((sd(ia,is)+srzwet(ia,is)) * acf)
        else
           srzwet(ia,is) = abs(sd(ia,is))
           srz(ia,is) = srz(ia,is) - srzwet(ia,is)
        endif
        qwet(is) = qwet(is) + (srzwet(ia,is) * acf)
c        SRZ(IA,is) = SRZ(IA,is) + (sd(ia,is)*(1-resp_coef(is)))
c         if (srz(ia,is).lt.0) then
c            srz(ia,is)=0
c         endif
        SD(IA,is)=0.
      ENDIF
C
C  ROOT ZONE CALCULATIONS
      SRZ(IA,is) = SRZ(IA,is) - (1.0-pmacro(is))*(P+qvpref(is))
      IF(SRZ(IA,is).LT.0.)THEN
        SUZ(IA,is) = SUZ(IA,is) - SRZ(IA,is)
        SRZ(IA,is) = 0.
      ENDIF
      suz(ia,is) = suz(ia,is)+ (pmacro(is)-pmac_sat(is))*(P+qvpref(is))
C
C  UZ CALCULATIONS
      IF(SUZ(IA,is).GT.SD(IA,is))THEN
         EX(IA) = SUZ(IA,is) - SD(IA,is)
         SUZ(IA,is)=SD(IA,is)
         QOF(is)=QOF(is)+EX(IA)*acf
         ACM(IS)=ACM(IS)+ACF
      ENDIF
c
c Compute total depth of water infiltrated into local wetness bin.
c This p is without infiltration excess and without the water
c that bypasses the unsaturated zone
c
      uz_infil(ia,is) = p-ex(ia)
     
C
C  CALCULATE DRAINAGE FROM SUZ
      IF(SD(IA,is).GT.0.)THEN
         UZ=SUZ(IA,is)*DT/(SD(IA,is)*TD(is))
         IF(UZ.GT.SUZ(IA,is))UZ=SUZ(IA,is)
         SUZ(IA,is)=SUZ(IA,is)-UZ
         IF(SUZ(IA,is).LT.0.0000001)SUZ(IA,is)=0.
         QUZ(is)=QUZ(is)+UZ*ACF
      ENDIF
c
c
c
c Track local recharge flux for weights in webmod_res.f
c
      quz_local(ia,is) = uz

C
C***************************************************************
C  CALCULATE EVAPOTRANSPIRATION FROM ROOT ZONE DEFICIT
C
      EA=0.
c
c  Modify to check if transpiration is on. If it is
c  then proceed to calculate actual ET:
c  If srz/srmax is 0.01 or less, then assume evaporation can occur 
c    over 100% of the basin (EA as originally calculated)
c  Else if transp_on then scale the predicted EA by the summer cover
c    density. This means that only transpiration can remove water from
c    the root zone.
c  Else if transpiration is off then scale EA by the ten percent of the
c    winter cover density since the transpiration rate is greatly reduced
c    during winter even though the vegetation can still intercept and
c    evaporate water.
c
      EA=EP*(1 - SRZ(IA,is)/SRMAX(is))
      IF(EA.GT.SRMAX(is)-SRZ(IA,is))EA=SRMAX(is)-SRZ(IA,is)
      IF(SRZ(IA,is)/SRMAX(is).GT.0.01) THEN
         IF(TRANSP_ON(IS).EQ.1) THEN
            EA = EA * COVDEN_SUM(IS)
         ELSE
            EA = EA * 0.1 * COVDEN_WIN(IS)
         ENDIF
      ENDIF

      SRZ(IA,is)=SRZ(IA,is)+EA
      sae_local(ia,is) = EA

c$$$      IF(EP.GT.0.)THEN
c$$$        EA=EP*(1 - SRZ(IA,is)/SRMAX(is))
c$$$        IF(EA.GT.SRMAX(is)-SRZ(IA,is))EA=SRMAX(is)-SRZ(IA,is)
c$$$        SRZ(IA,is)=SRZ(IA,is)+EA
c$$$      ENDIF
      SUMAE(is) = SUMAE(is) + EA * ACF
      sae(is) = sae(is) + EA *ACF
c
c     Calculate balance terms per step
c
      srz_sc(is) = srz_sc(is) + SRZ(IA,is)*ACF
      suz_sc(is) = suz_sc(is) + SUZ(IA,is)*ACF

C
C***************************************************************
C
C
C  CALCULATION OF FLOW FROM FULLY SATURATED AREA
C  This section assumes that a/tanB values are ordered from high to low
C
!      OF=0.
!      IF(IA.GT.1)THEN
!        IB=IA-1
!        IF(EX(IA).GT.0.)THEN
!c  Both limits are saturated
!          OF=AC(IA,is)*(EX(IB)+EX(IA))/2
!          ACM(IS)=ACM(IS)+ACF
!          ihour(ib,is) = ihour(ib,is) + 1
!         ELSE
!c  Check if lower limit saturated (higher a/tanB value)
!             IF(EX(IB).GT.0.)THEN
!               ACF=ACF*EX(IB)/(EX(IB)-EX(IA))
!               OF=ACF*EX(IB)/2
!               ACM(IS)=ACM(IS)+ACF
!               ihour(ib,is) = ihour(ib,is) + 1
!             ENDIF
!         ENDIF
!      ENDIF
!      QOF(is)=QOF(is)+OF
C
C  Set contributing area plotting array
c      CA(IT) = ACM
c

C      write(80,1009) nstep,ia, is, acf, sd(ia,is),sbar(is),szm(is),
C     +  tl(is), st(ia,is),srz(ia,is), suz(ia,is), uz
C 1009 format(3I4, 9e10.3)

cccccccMove to end so that depth reflect ending water table after ET
c Add uz_depth to indicate the initial water content in the
c unsaturated zone. uz_depth is determined by the saturation deficit,
c the root zone deficit, and the unsaturated zone storage
c
c$$$      z_wt_local(IA,is) = - (sd(ia,is)/s_drain(is))
c$$$      z_wt(is) = z_wt(is) + z_wt_local(ia,is) * acf
c$$$
c$$$      if(z_wt_local(ia,is).ge.-s_root_depth(is)) then 
c$$$            ! just account for the root zone deficit plus suz
c$$$         uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
c$$$     $        s_theta_wp(is)*s_root_depth(is) + suz(ia,is)
c$$$      else
c$$$            ! add additional UZ water at field capacity
c$$$         uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
c$$$     $        s_theta_wp(is)*s_root_depth(is) + suz(ia,is) -
c$$$     $        (s_root_depth(is)+z_wt_local(ia,is))*s_theta_fc(is)
c$$$      end if


C  END OF A/TANB LOOP
      end if
 30   CONTINUE

C
C  ADD INFILTRATION EXCESS
      QOFS(is)=QOF(is)
      QOF(is)=QOF(is)+REX(is)
c
c Compute surface deposition and total infiltration for mru and basin variables
c
      surfdep(is) = infil(is)
      infil(is) = infil(is) - qof(is)
c
c Variable infilitration excess areas now. Independent variables
c afx(is) and afxmax(is) will track area of ponded water
c for each subcatchment. So comment out next line - RMTW
c
c
c$$$      IF(IROF(is).EQ.1)ACMAX(is)=1.
      IF(ACM(is).GT.ACMAX(is))ACMAX(is)=ACM(is)

      if(afx(is).gt.afxmax(is))afxmax(is) = afx(is)
c
c  Compute loss from saturated zone to deep aquifer system using
c  the water table height from the previous time step.
c
      sat_head = s_rock_depth(is)+z_wt(is)
      if(sat_head.gt.0.0.and.gw_loss_k(is).lt.11) then
         rate = 10**(-gw_loss_k(is))
         gw_loss(is) = sat_head * rate *36.0* dt
      else
         gw_loss(is) = 0.0
      end if
c
c Warn user if simulations indicate water table near or below the bedrock depth      
      sbar_norm = sbar(is)/sbar_max(is)

      if(sbar_norm.gt.1.0) then
         sbar_norm = 1.0
         print*,'Water table in MRU ',is,
     $        ' is close to or below bedrock depth'
      end if

c
c Limit irrigation amounts requested for next day, irrig_int_next, to
c the pumping rate throttled using a quadratic function to account for
c decreased pumping capacity in depleted surface-water aquifers. Use
c the readily drainable soil porosity for estimating the maximum
c drawdown
c
      irrsched = irrig_sched_int(is)
      if(irrsched.gt.0) then    ! internal irrigation source present
         irrsrc = irrig_int_src(is)
         if(irrsrc.eq.0) then   ! saturated zone is source of internal irrigation
            pump_coeff = sqrt(1.-sbar_norm**2)
            irrig_dep_max = pump_coeff*irrig_int_max(is)*
     $           gpm2m3ph/sc_area(is)/1e6/inch2m*dt ! max inches per time step
            if(irrig_int_next(irrsched).le.irrig_dep_max) then
               irrig_sat_next(is) = irrig_int_next(irrsched)
            else
               irrig_sat_next(is) = irrig_dep_max
               if(irr_warn.eq.0) then
                  print*,'Irrigation rate reduced because of limited'//
     $                 ' pumping capacity on time step', nstep,
     $                 ' for MRU ',is
                  print*,' This warning will not be repeated'
                  irr_warn = 1
               end if
            end if
         else
            irrig_sat_next(is) = 0.0
         end if
      end if
c
c Add groundwater influx from external sources
c
      gw_in1(is)=0.0
      gw_in2(is)=0.0
      gwsched=sched_gw1(is)
      if(gwsched.ne.0) then
         gw_in1(is) = gw_ext(gwsched)*qlin2dep*gwbnd_len1(is)*
     $        dt/sc_area(is)
      end if
      gwsched=sched_gw2(is)
      if(gwsched.ne.0) then
         gw_in2(is) = gw_ext(gwsched)*qlin2dep*gwbnd_len2(is)*
     $        dt/sc_area(is)
      end if

C  CALCULATE SATURATED ZONE DRAINAGE
c
c  Added direct flow to model quick flow through macropores per
c  work by Dave Kinner - RW
c
      qdf(is)=qdffrac(is)*quz(is)
      quz(is)=quz(is)-qdf(is)

      QB(is)=SZQ(is)*EXP(-SBAR(is)/SZM(is))
c
c  Compute drainage from preferential flow in the unsaturated zone using
c  the water table height from the previous time step.
c  v2.0 add infiltration and subtract recharge and direct flow to avoid 
c  one day time lag in preferential flow.
      z_wt_quick = z_wt(is) +
     $     (+QUZ(is)-QB(is)-qwet(is)+qvpref(is)-qexfil(is)
     $     -(irrig_sat_mru(is)*inch2m)+
     $     gw_in1(is)+gw_in2(is)-gw_loss(is))/s_drain(is)
c      z_wt_quick = z_wt(is)
      if(z_wt_quick.le.s_satpref_zmin(is).or.
     $     qpref_max(is).lt.1e-12) then
         qpref(is) = 0.0
      else if(z_wt_quick.ge.s_satpref_zmax(is)) then
         qpref(is) = qpref_max(is) * dt
      else
         qpref(is) = (z_wt_quick - s_satpref_zmin(is)) /
     $        (s_satpref_zmax(is)-s_satpref_zmin(is))*qpref_max(is)
      end if
c
c  Calculate qpref such that it cannot drain to an elevation below the tile drain
c
      if(qpref(is).gt.z_wt_quick - s_satpref_zmin(is))
     $   qpref(is)=z_wt_quick-s_satpref_zmin(is)+(-QUZ(is)-qvpref(is)
     $     +QB(is)+qwet(is)+qexfil(is)+(irrig_sat_mru(is)*inch2m)-
     $     gw_in1(is)-gw_in2(is) + gw_loss(is))/s_drain(is)
      if(qpref(is).lt.0) qpref(is) = 0.0

      
c      QB(is)=SZQ(is)*EXP(-SBAR(is)/SZM(is))
c
c Adjust tile drain (qpref), QB and qexfil to allow the tile drain
c to be the dominant flow path when the water table is above the top of the
c drain. This implies that the conductivity of the tile drain should be higher than
c the Ksat at the depth of the tile drain.
c
      if(qpref_max(is).gt.0) then
         qprefwt = qpref(is)/qpref_max(is)
         qb(is) = qb(is) - qpref(is)*qprefwt
         qexfil(is)=qexfil(is) - qpref(is)*qprefwt
         if(qb(is).lt..00001) qb(is)=.00001
         if(qexfil(is).lt..00001) qexfil(is)=.00001
      end if
      
      SBAR(is)=SBAR(is)-QUZ(is)-qvpref(is)+QB(is)+qwet(is)+
     $     qexfil(is)+qpref(is)+(irrig_sat_mru(is)*inch2m)-
     $     gw_in1(is)-gw_in2(is) + gw_loss(is)
      QOUT(is)=QB(is)+QOF(is)+qdf(is)+qexfil(is)+qpref(is)

cccccccMove to end so that depth reflect ending water table after ET
c$$$c Add uz_depth to indicate the final water content in the
c$$$c unsaturated zone. uz_depth is determined by the saturation deficit,
c$$$c the root zone deficit, and the unsaturated zone storage
c$$$c
c
c      z_wt_local(ia,is) = 0.0
      z_wt(is) = 0.0
      do 250 ia = 1, nacsc(is)
        if(ia.eq.nacsc(is)) then
          ACF=0.5*AC(IA,is)
        else
          ACF=0.5*(AC(IA,is)+AC(IA+1,is))
        endif

         sd_temp = 
     $        SBAR(is)+(SZM(is)*(TL(is)-ST(IA,is)))*resp_coef(is)
c$$$         if(is.eq.1) print*,ia,sd_temp
         if(sd_temp.lt.0) sd_temp = 0.0
         z_wt_local(IA,is) = - sd_temp/s_drain(is)
         z_wt(is) = z_wt(is) + z_wt_local(ia,is) * acf

         if(acf.ne.0.0) then

            if(z_wt_local(ia,is).ge.-s_root_depth(is)) then 
                                !  account for the root zone deficit plus suz
               uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
     $              s_theta_wp(is)*s_root_depth(is) + suz(ia,is)
            else
                                !  additional UZ water at field capacity
               uz_depth(ia,is) = srmax(is) - srz(ia,is) + 
     $              s_theta_wp(is)*s_root_depth(is) + suz(ia,is) -
     $              (s_root_depth(is)+z_wt_local(ia,is))*s_theta_fc(is)
            end if
c
c Use fluxes and change in uz_depth to calculate the transfer of
c water from the unsaturated zone to the saturated zone that occurs
c when the water table changes. This is just a convoluted way of tracking
c the water table changes while eliminating variations in the root zone
c and the unsaturated zone storage.
c
c$$$         uz2sat(ia,is) = (uz_depth(ia,is) - last_uz_dep(ia,is))
c$$$     $        - (srz(ia,is) - last_srz(ia,is))


c$$$            uz2sat(ia,is) = (last_uz_dep(ia,is) - uz_depth(ia,is))
c$$$     $           + (last_srz(ia,is)-srz(ia,is))
c$$$     $           - (last_suz(ia,is)-suz(ia,is))
c
c Try just using the depth to water table for uz2sat. Since the control volume
c for the UZ is all porosity below field capacity in the rooting zone, no
c net flux resulting from water table changes (uz2sat) can occur unless
c either the current or previous water table was below the root zone depth.
c
            if(z_wt_local(ia,is).ge.-s_root_depth(is)) then
               if(last_z_wt_local(ia,is).ge.-s_root_depth(is)) then
                  uz2sat(ia,is)=0.0
               else
                  uz2sat(ia,is)=(-s_root_depth(is)-
     $                 last_z_wt_local(ia,is))*s_theta_fc(is)
               end if
            else
               if(last_z_wt_local(ia,is).ge.-s_root_depth(is)) then
                  uz2sat(ia,is)=(s_root_depth(is)+
     $                 z_wt_local(ia,is))*s_theta_fc(is)
               else
                  uz2sat(ia,is) = (z_wt_local(ia,is) -
     $                 last_z_wt_local(ia,is))*s_theta_fc(is)
               end if
                     
            end if

         else
            uz_depth(ia,is) = 0.0
            uz2sat(ia,is) = 0.0
         end if
         
 250  continue
C
C  Sum different flow paths
      SUMQ(is)=SUMQ(is)+QOUT(is)
      sumqexfil(is) = sumqexfil(is)+qexfil(is)
      sumqwet(is) = sumqwet(is) + qwet(is)
      sumqrex(is)=sumqrex(is)+rex(is)
      sumqofs(is)=sumqofs(is)+qofs(is)
      sumqb(is)=sumqb(is)+qb(is)
      sumquz(is)=sumquz(is)+quz(is)
      sumqdf(is)=sumqdf(is)+qdf(is)
      sumqvpref(is)=sumqvpref(is)+qvpref(is)
      sumqpref(is)=sumqpref(is)+qpref(is)
      sumgw_in1(is)=sumgw_in1(is)+gw_in1(is)
      sumgw_in2(is)=sumgw_in2(is)+gw_in2(is)
      sumgw_loss(is)=sumgw_loss(is)+gw_loss(is)

C
C  CHANNEL ROUTING CALCULATIONS
C  allow for time delay to catchment outlet ND as well as 
C  internal routing array
c$$$      DO 40 IR=1,NR
c$$$      IN= 1+ND+IR-1
c$$$c      IF(IN.GT.NSTEP)GO TO 10
c$$$      Q(IN,is)=Q(IN,is)+QOUT(is)*AR(IR,is)
c$$$   40 CONTINUE
   
c      qscm(is) = q(1,is)/ area(is)
c      qsccfs(is) = (qscm(is) * sc_area(is) * 9808.333) / dt
c      qtot = qtot + qscm(is)
C
c
c     Print topmodel fluxes every time step if iout = 2.
c
c     Note that step1 does not use declpri. because we want it
c     to reset to .true. on a restart with a init file. so that
c     the header gets printed. - RMTW
c
!      If(IOUT.ge.2) then
!         if (step1)then
!            Write(topout_file_unit,8050)
c$$$ 8050 format(//' MRU step     p       gw_in      ep       qout    ',
c$$$     $         '   quz        qb       sbar      qof    ',
c$$$     $         '   acm       afx'/
c$$$     $         ' === ==== ========= ========= ========= ========= ',
c$$$     $         '========= ========= ========= ========= ',
c$$$     $         '========= =========')
! 8050 format(//' MRU step    surfdep infil p       gw_in   soil_et   q'
!     $         'out    quz        qb      gwloss    sbar       qof    ',
!     $         '   suz       afx'/
!     $         ' === ==== ========= ========= ========= ========= ',
!     $         '========= ========= ========= ========= ========= ',
!     $         '========= =========')
!            step1=.false.
!         endif

!        write(topout_file_unit,1000) is, nstep,surfdep(is),infil(is),
!     $   p, gw_in1(is)+ gw_in2(is), sae(is), qout(is),
!     +   quz(is), qb(is), gw_loss(is),sbar(is), qof(is),
!     $   acm(is), afx(is)
! 1000   format(i4,i5, 13e10.3)
!        write(topout_file_unit,1000) is, nstep,surfdep(is),infil(is),
!     $   p, gw_in1(is)+ gw_in2(is), sae(is), qout(is),
!     +   quz(is), qb(is), gw_loss(is),sbar(is), qof(is),
!     $   acm(is), afx(is)
! 1000   format(i4,i5, 13e10.3)
!      endif

      basin_actet = basin_actet + (sae(is)*area(is))
c$$$      if (is.eq.1) print*,sae(is),basin_actet
  100 continue
c
c  Calculate basin area-weighted discharge and
c  root zone soil moisture 
c  
c      qbasinm = 0.0
      do 110 is = 1,nsc
c         qbasinm = qbasinm + (qscm(is)*sc_area(is))
         smav_sc(is) = srmax(is) - srz_sc(is) + 
     $        s_theta_wp(is)*s_root_depth(is)
         smcont_sc(is) = smav_sc(is)/s_root_depth(is)
         srz_basin = srz_basin + (srz_sc(is)*area(is))
         suz_basin = suz_basin + (suz_sc(is)*area(is))
         smav_basin = smav_basin + (smav_sc(is)*area(is))
         smcont_basin = smcont_basin + (smcont_sc(is)*area(is))
         sbar_basin = sbar_basin + (sbar(is)*area(is))
         qexfil_basin = qexfil_basin + (qexfil(is)*area(is))
         qwet_basin = qwet_basin + (qwet(is)*area(is))
         rex_basin = rex_basin + (rex(is)*area(is))
         afx_basin = afx_basin + (afx(is)*area(is))
         qofs_basin = qofs_basin + (qofs(is)*area(is))
         acm_basin = acm_basin + (acm(is)*area(is))
         qof_basin = qof_basin + (qof(is)*area(is))
         basin_infil = basin_infil + (infil(is)*area(is))
         basin_surfdep = basin_surfdep + (surfdep(is)*area(is))
         qdf_basin = qdf_basin + (qdf(is)*area(is))
         qvpref_basin = qvpref_basin + (qvpref(is)*area(is))
         qpref_basin = qpref_basin + (qpref(is)*area(is))
         gw_in1_basin = gw_in1_basin + (gw_in1(is)*area(is))
         gw_in2_basin = gw_in2_basin + (gw_in2(is)*area(is))
         qb_basin = qb_basin + (qb(is)*area(is))
         qout_basin = qout_basin + (qout(is)*area(is))
         z_wt_basin = z_wt_basin + (z_wt(is)*area(is))
         gwl_basin = gwl_basin +(gw_loss(is)*area(is))
 110  continue

c
c Print summary at end of run
c
      if (end_run.and.iout.ge.1) then

!         If(IOUT.ge.1) Write(topout_file_unit,650)
         Write(topout_file_unit,650)
 650     FORMAT(/'Water Balance for Subcatchments : '//
     1        '  MRU       SUMP    SUMGWIN     SUMAE    ',
     $        '  SUMQ     SUMGW_LOSS   SUMRZ      SUMUZ  ',
     $        '   SBAR      BAL       ACMAX AFXMAX'/
     $        '  ===   ========== ========== ========== ',
     $        '========== ========== ========== ========== ',
     $        '========== ==========  ===== ======')

C  CALCULATE BALANCE TERMS

         do 5000 is= 1, nsc
            SUMRZ = 0.
            SUMUZ = 0.
            DO 50 IA =1,nacsc(is)
               if(ia.eq.nacsc(is)) then
                 ACF=0.5*AC(IA,is)
               else
                 ACF=0.5*(AC(IA,is)+AC(IA+1,is))
               endif
               SUMRZ = SUMRZ + SRZ(IA,is)*ACF
               SUMUZ = SUMUZ + SUZ(IA,is)*ACF 
 50         CONTINUE

            BAL(is) = BAL(is) + SBAR(is) +SUMP(is)+sumgw_in1(is)+
     +          sumgw_in2(is)- SUMAE(is)- SUMQ(is) + SUMRZ - SUMUZ -
     $          sumgw_loss(is) - (s_rock_depth(is) * s_porosity(is))

!            If(IOUT.ge.1) Write(topout_file_unit,655) is,SUMP(is),
            Write(topout_file_unit,655) is,SUMP(is),
     $        sumgw_in1(is)+sumgw_in2(is),SUMAE(is),SUMQ(is),
     $        sumgw_loss(is),SUMRZ,SUMUZ,SBAR(is),BAL(is),
     $        acmax(is), afxmax(is)
 655        FORMAT(i4,3x,9e11.3,2(2x,f5.2))

 5000    continue

      end if
c
c End of Summary Info
c
c
c Leftovers from earlier versions
c
c      qbasincfs = (qbasinm * 9808.333) / dt
c      qbasinm = qbasinm / basin_area


c---  compute statistics.
c
c These calculations and the final printing were moved to
c the route_clarke routine since the observed values should
c be compared with the predictions of discharge at the outlet
c where stream discharge measurements are assumed to be made.
c

c      SUMQOBS =SUMQOBS + QOBS(qobsta)
c      SSQ = SSQ + (QOBS(qobsta))**2
c      F1=F1 + (qbasincfs - QOBS(qobsta))**2
c      F2=F2 + ABS(qbasincfs - QOBS(qobsta))
c      NTSTEP = NTSTEP + 1

      topmrun = 0

      return
      end

c$$$c***********************************************************************
c$$$c
c$$$c     sumbclean - Close the topout file
c$$$c
c$$$
c$$$      integer function topclean(SBAR, SUZ, SRZ, AC, BAL, nsc, nacsc, 
c$$$     +                            iout, SUMP, SUMAE, SUMQ, ACMAX,
c$$$     +                    sumqrex, sumqofs, sumqb, sumquz, sumqdf,
c$$$     +                    qdffrac, afxmax, topout_file_unit)
c$$$
c$$$      include 'fmodules.inc'
c$$$
c$$$      integer nsc, nac, iout, is, IA, nacsc(nmru)
c$$$c$$$      integer nsc, nac, iout, is, IA, NTSTEP, nacsc(nmru)
c$$$
c$$$      real SUZ(MAXNAC,nmru), SRZ(MAXNAC,nmru)
c$$$      real BAL(nmru), AC(MAXNAC,nmru), SBAR(nmru)
c$$$      real SUMP(nmru), SUMAE(nmru), SUMQ(nmru), ACMAX(nmru)
c$$$      real sumqrex(nmru), sumqofs(nmru), sumqb(nmru)
c$$$      real sumquz(nmru), sumqdf(nmru), qdffrac(nmru)
c$$$      real SUMRZ, SUMUZ, ACF, afxmax(nmru)
c$$$c$$$      real SUMQOBS, SSQ, F1, F2
c$$$c$$$      real QBAR, VARQ, VARE, E
c$$$      integer topout_file_unit
c$$$
c$$$      topclean = 1
c$$$
c$$$      If(IOUT.ge.1) Write(topout_file_unit,650)
c$$$  650 FORMAT(/'Water Balance for Subcatchments : '//
c$$$     1 '  MRU       SUMP      SUMAE       SUMQ      SUMRZ   ',
c$$$     2 '  SUMUZ      SBAR        BAL      ACMAX AFXMAX'/
c$$$     $ '  ===   ========== ========== ========== ========== ',
c$$$     $ '========== ========== ==========  ===== ======')
c$$$
c$$$C  CALCULATE BALANCE TERMS
c$$$
c$$$      do 100 is= 1, nsc
c$$$
c$$$        nac = nacsc(is)
c$$$        SUMRZ = 0.
c$$$        SUMUZ = 0.
c$$$        DO 50 IA =1,NAC
c$$$          ACF=0.5*(AC(IA,is)+AC(IA+1,is))
c$$$          SUMRZ = SUMRZ + SRZ(IA,is)*ACF
c$$$          SUMUZ = SUMUZ + SUZ(IA,is)*ACF 
c$$$   50   CONTINUE
c$$$
c$$$      BAL(is) = BAL(is) + SBAR(is) +SUMP(is) - SUMAE(is) - SUMQ(is) 
c$$$     +          + SUMRZ - SUMUZ
c$$$
c$$$      If(IOUT.ge.1) Write(topout_file_unit,655) is,SUMP(is),
c$$$     $     SUMAE(is),SUMQ(is),SUMRZ,SUMUZ,SBAR(is),BAL(is),
c$$$     $     acmax(is), afxmax(is)
c$$$ 655  FORMAT(i4,3x,7e11.4,2(2x,f5.2))
c$$$
c$$$c$$$      If(IOUT.ge.1)WRITE(topout_file_unit,651)ACMAX(is),
c$$$c$$$     $     afxmax(is), qdffrac(is)
c$$$c$$$ 651  FORMAT(9X,'Maximum contributing area ', e12.5/
c$$$c$$$     $     9X,'Maximum ponded area (infiltration excess) ', e12.5/
c$$$c$$$     1     9X,'Proportion QUZ delivered as direct flow', f12.2)
c$$$
c$$$ 100  continue
c$$$
c$$$c$$$c-- compute statistics
c$$$c$$$
c$$$c$$$
c$$$c$$$      QBAR = SUMQOBS / NTSTEP
c$$$c$$$      VARQ = (SSQ/NTSTEP - QBAR*QBAR)
c$$$c$$$      VARE = F1/NTSTEP
c$$$c$$$      E=1-VARE/VARQ
c$$$c$$$c
c$$$c$$$c  add objective function values to output file
c$$$c$$$      write(topout_file_unit,621)f1,e,f2,qbar,varq,vare
c$$$c$$$c      write(10,621)f1,e,f2,qbar,varq,vare
c$$$c$$$  621 format(//1x,'Objective function values'/
c$$$c$$$     1 1x,'F1 ',e12.5,'   E ',f12.5,'   F2 ',e12.5//
c$$$c$$$     2 1x,'Mean Obs Q ',e12.5,'    Variance Obs Q ',e12.5/
c$$$c$$$     3 '     Error Variance',e12.5)
c$$$c
c$$$c
c$$$c
c$$$
c$$$
c$$$c      close (unit=80)
c$$$
c$$$      topclean = 0
c$$$
c$$$      return
c$$$      end


c********************************************************************

C
      SUBROUTINE EXPINF(IROF,IT,RINT,DF,CUMF,XK0,HF,DTH,SZM,DT,
     +                  const,tp)
C
c      INCLUDE TMCOMMON.FOR
      DOUBLE PRECISION CONST,SUM,FC,FUNC,CD,SZF,XKF

C*************************************************************
C
C  SUBROUTINE TO CALCULATE INFILTRATION EXCESS RUNOFF USING THE
C  EXPONENTIAL GREEN-AMPT MODEL.
C
C**************************************************************
C
C
      real XK0, HF, DTH, SZM, E
      real RINT,DF,CUMF, F, F1, F2, R2
      real TP, DT, FAC, ADD
      real DFUNC
      
      integer IROF,IT, i, j



C  Note that HF and DTH only appear in product CD
      E = 0.001
      CD=HF*DTH
      SZF = 1./SZM
      XKF = XK0
      F = CUMF
      F1 = 0
      IF(IROF.EQ.1)GO TO 10
C  PONDING HAS ALREADY OCCURRED - GO TO EXCESS CALCULATION
C
      IF(CUMF.EQ.0.)GOTO 7
C  FIRST TIME STEP, OVERFLOW IF CUMF=0, GO DIRECT TO F2 CALCULATION
C  INITIAL ESTIMATE OF TIME TO PONDING
      F1=CUMF
      R2=-XKF*SZF*(CD+F1)/(1-EXP(SZF*F1))
      IF(R2.LT.RINT)THEN
C  PONDING STARTS AT BEGINNING OF TIME STEP
      TP=(IT-1.)*DT
      IROF=1
      F=CUMF
      GO TO 8
      ENDIF
    7 F2=CUMF+DT*RINT
      IF(F2.EQ.0.)GO TO 20
      R2=-XKF*SZF*(CD+F2)/(1-EXP(SZF*F2))
      IF(R2.GT.RINT)GO TO 20
      F=CUMF+R2*DT
      DO 9 I=1,20
      R2=-XKF*SZF*(CD+F)/(1-EXP(SZF*F))
      IF(R2.GT.RINT)THEN
        F1=F
        F=(F2+F)*0.5
      IF(ABS(F-F1).LT.E)GO TO 11
      ELSE
        F2=F
        F=(F1+F)*0.5
      IF(ABS(F-F2).LT.E)GO TO 11
      ENDIF
    9 CONTINUE
      WRITE(6,600)
  600 FORMAT(1X,'MAXIMUM NO OF ITERATIONS EXCEEDED')
   11 CONTINUE
      TP=(IT-1)*DT+(F-CUMF)/RINT
      IF(TP.GT.IT*DT)GO TO 20
C
C  SET UP DEFINITE INTEGRAL CONSTANT USING FP
C
    8 CONST =0
      FAC=1
      FC=(F+CD)
      DO 12 J=1,10
      FAC=FAC*J
      ADD=(FC*SZF)**J/(J*FAC)
      CONST=CONST+ADD
   12 CONTINUE
      CONST=DLOG(FC)-(DLOG(FC)+CONST)/DEXP(SZF*CD)
C
      IROF=1
      F=F+0.5*RINT*(IT*DT-TP)
   10 CONTINUE
C
C  NEWTON-RAPHSON SOLUTION FOR F(T)
      DO 14 I=1,40
C
C  CALCULATE SUM OF SERIES TERMS
      FC=(F+CD)
      SUM=0.
      FAC=1.
      DO 13 J=1,10
      FAC=FAC*J
      ADD=(FC*SZF)**J/(J*FAC)
      SUM=SUM+ADD
   13 CONTINUE
      FUNC=-(DLOG(FC)-(DLOG(FC)+SUM)/DEXP(SZF*CD)-CONST)/(XKF*SZF)
     1     -(IT*DT-TP)
      DFUNC=(EXP(SZF*F)-1)/(XKF*SZF*FC)
      DF=-FUNC/DFUNC
      F=F+DF
      IF(ABS(DF).LE.E)GO TO 15
   14 CONTINUE
      WRITE(6,600)
 15   CONTINUE
      IF(F.LT.CUMF+(RINT*DT))THEN
C      IF(F.LT.CUMF+RINT)THEN
      DF=F-CUMF
      CUMF=F
C  SET UP INITIAL ESTIMATE FOR NEXT TIME STEP
c  Commented out since F is not passed through argument list.
c  Reinitialize F to CUMF upon entry to subroutine. - RMTW
c     F=F+DF
      RETURN
      ENDIF
   20 CONTINUE
C  THERE IS NO PONDING IN THIS TIME STEP
      IROF=0
      DF = RINT*DT
      CUMF=CUMF+DF
      RETURN
      END
C
