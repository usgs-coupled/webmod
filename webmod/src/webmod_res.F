#include "defines.h"
***********************************************************************
c   webmod_res.f
c
c   last modified by Rick Webb $Date: 2009-04-24$
c
c   This module collects water flux information from other modules
c   linked in the xtop_prms model. Residual water contents and
c   reservoir volumes and mixing ratios are then calculated for
c   input into the phreeq_mms module.
c
c   webmod_res is therefore an interface module that can be used
c   as a template for preparing fluxes computed in other modules,
c   xprms for example, for input into phreeq_mms. All depths are
c   to be in meters. If they are not they must be converted for
c   volumes to be calculated correctly.
c
c   Fluxes begin with the canopy interception,
c   proceed through snowpack and hillslope reservoirs, and then
c   finish with inputs into the stream routing.
c
c   24Apr09 - Added FORTRAN 90 Module WEBMODRES
c
c   17Sep03 - Created with RCS version control
c 
c     
c  This module can be adapted to work with other models by performing the
c  necessary water tracking calculations to provide a list of reservoirs
c  with water fluxes into and out of each similar reservoir.
c
c  The inputs to phreeq_mms are series of volumes, in cubic meters,
c  of inputs for a forward-feeding constantly mixed reactor series.
c
c  For example,
c  inputs for the root zone would be
c  root_zone_1* = root_zone_0 + precip + snowmelt
c                 + recharge from the saturated zone
c
c  ET (pure water) removed from root_zone_1* results in root_zone_1 that
c  provides output chemistry for recharge, QUZ.
c
c  QUZ_1 = QUZ_0 + root_zone_1
c
c  and so on for the cascading reservoirs
c  
c The main structural element is a volume mixing matrix, for example,
c vmix_can to describe volumes entering and leaving the canopy.
c All mixing matrices are dimensioned (nmru,nresinp) except for the
c unsaturated zone, vmix_uz(nac,nmru,nresinp), the reservoirs
c that distribute hillslope exports to strema segments, 
c vmix_hillexp(nhydro,nmru), and the stream segments, vmix_stream(nhydro).
c
c A hillslope reservoir, vmix_hill, is provided for mixing all unique
c hillslope inputs before routing to the drainage system, vmix_stream.
c
c The beginning and ending volumes for the hillslope reservoir, and all
c other transient reservoirs, will always be zero for any time step.
c
c The final index of any hillslope mixing matrix contains nresinp values.
c The index nresinp begins with four volumes for the current reservoir
c type:

c     1) the initial volumes
c     2) the sum of input volumes,
c     3) the sum of descharge volumes,
c     4) and the reservoir volume at the end of the time step.
c These are followed by volumes from all potential inputs (zero if no 
c input from that source in that timestep). For WEBMOD_RES inputs and
c matrix placements include the following:
c     5) precipitation (rain or snow)
c     6) ET (considered a negative input for mixing)
c     7) impervious surface (not implemented yet)
c     8) canopy or vegetation
c     9) snow melt
c    10) O-horizon (flushed with overland flow, reservoir leaf-on water)
c    11) unsaturated zone (root zone, recharge, residual water)
c    12) shallow unsaturated zone (fixed volume flushed with direct flow)
c    13) saturated zone (base flow)
c    14) exfiltrated water (springs or wetlands)(inactive)
c    15) deep preferential flow (tile drains or other)
c    16) hillslope reservoir (mixing box that feeds stream segments)
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c    20) groundwater inflow from outside of basin from 1st source
c    21) groundwater inflow from outside of basin from 2nd source
c           Two sources of groundwater from external sources are
c           made available to simulate sources from a leaky drainage
c           canal and influx from ground water sources upgradient of
c           the canal.
c
c vmix_imp, impervious area inputs (ET is considered nill
c  as this water is expected to move quickly into storm sewers)
c     5) precip
c     9) melt
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c           impervious area outputs: uz, hill reservoir
c
c vmix_can, canopy inputs
c     5) precip
c     6) et
c    10) O-horizon, only on first day of tranpiration when canopy
c        increases from winter to summer cover density.
c    11) unsaturated zone (transpiration)
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c           canopy output (throughfall): snowpack, root zone, 
c                          O-horizon for overland flow and on the first
c                          day when transpiration stops (leaves-off)
c
c vmix_snow, snowpack inputs
c     5) precip
c     6) et
c     8) canopy
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c            snowpack output: root zone, o_horizon
c
c vmix_ohoriz, O-horizon inputs (fixed volume for now, provides solutes
c              for overland flow to pick up, and provides cyclical reservoir
c              for canopy moisture on days of leaves_on and leaves_off)
c     5) precip
c     8) canopy (on first day of transpiration)
c     9) melt
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c              O-horizon outputs: stream, canopy first day of no transpiration 
c
c vmix_uz, Unsaturated zone (combination of root zone and unsaturated zone -
c                This reservoir has three dimensions nac, nmru, and nresinp)
c     5) precip
c     6) et
c     8) canopy
c     9) melt
c    14) groundwater (uz2sat-srz_wet)
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c          UZ outputs: saturated zone, canopy(transpiration),shallow
c                      preferential flow (qdf)
c
c vmix_uzgen, Composite of unsaturated zone (combination of root zone and unsaturated zone -
c vmix_uzrip, riparian uz composite (as determined using riparian_thresh
c vmix_uzup, uplands uz composite
c                These reservoir has two dimensions nmru, and nresinp)
c     5) precip
c     6) et
c     8) canopy
c     9) melt
c    14) groundwater (uz2sat-srz_wet)
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c
c vmix_qdf, Direct flow (preferential flow through the unsaturated zone, qdf)
c     
c     5) precip
c     8) canopy
c     9) melt
c    11) quz, discharge from storage in the unsaturated zone
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c  other indices zero
c
c           Direct flow output: hill reservoir (mixing box before stream)
c
c vmix_sat, Saturated zone
c     5) precip
c     8) canopy
c     9) melt
c    11) quz, recharge from storage in the unsaturated zone
c    17) irrigation from well
c    18) irrigation from diversion
c    19) irrigation from external sources
c    20-21) ground water influx from channel(gw_in1) or upgradient(gw_in2)
c  other indices zero
c
c           Saturated zone output: deep preferential flow, 
c                                  hill reservoir,irrigation
c
c vmix_satpref, preferential flow in the saturated zone (i.e. tile drains)
c    13) saturated zone
c  other indices zero
c
c               satpref output: hill reservoir
c
c vmix_hill, transient storage for mixing all hillslope inputs
c      before stream routing
c     7) impermeable areas
c    10) O-horizon
c    12) qdf, preferential flow in the unsaturated zone
c    13) qb, baseflow discharge from the saturated zone
c    14) exfiltration (inactive)
c    15) satpref, preferential flow in the saturated zone
c  other indices zero
c
c            hill output: vmix_hillexp (mixing matrix to stream)
c
c vmix_hillexp, distributes each hillslop export to individual
c      drainage segments. It has indices of nhydro, set to 
c      clark_segs, and nmru.
c
c            hillexp output: vmix_stream
c
c vmix_stream, the drainage segments. The combined hillslope inputs are 
c      mixed with the existing stream segment contents and then translated
c      downstream one segment until exiting the basin. The exported volume
c      for each time step is recorded in vmix_basin(3).
c
c In addition to the transient hillslope reservoir, five more
c transient reservoirs will be used:
c
c vmix_uz2can(nac,nmru) will be used to mix unsaturated-zone water
c      from individual wetness index bins to provide water chemistry
c      for the flux transpired up to the canopy.
c
c vmix_sat2uz(nmru) will be used to track how much total water to
c      transfer from the saturated zone to the unsaturated zone as
c      a result of wetting the root zone and from stranded pore water
c      resulting from the water table lowering.
c
c vmix_uz2sat(nac,nmru) will be used to mix uz waters for recharge.
c      This water is apportioned to preferential flow and recharge
c      of the saturated zone.
c
c vmix_well(nmru) will be used to track water pumped from the
c       saturated zone. It will be equal to to irrig_sat_mru.
c
c vmix_diversion(nhydro) will be used to track water pumped or
c      diverted from stream segments.
c
c 27apr04 - Modified to include flux of potential infiltration
c      of precip, melt, and throughfall directly to saturated zone.
c
c 14jun04 - Added irrigation from internal and external sources
c
c 29sep04 - Added ground water loss from MRU and channel segments
c
c 22nov04 - Added basin summary variables (_cm) to be passed to
c      the summary module web_sum.f
c
c   may09 - Add Fortran90 Module: WEBMOD_ROUTE
c
c
c***********************************************************************c
c***********************************************************************
      MODULE WEBMOD_RESMOD

      IMPLICIT NONE

#include "fmodules.inc"

      real a_million, inch2m
      DATA a_million / 1e6 /
      DATA inch2m / 0.0254 /
! m3cm converts cubic meters to cm depth using basin area. Set in Init section
      double precision m3cm

      logical, save :: resstep1
      integer, save :: nmru, nac, nobs, nresinp, nsolute, nchemvar
      integer, save :: nchemobs, nchemdat, nchem_ext
c      integer, save :: clark_segs, nhydro, nchan,qobsta
      integer, save :: clark_segs, nhydro, nchan
c      integer, save :: nirrig_int, nirrig_ext
      integer, save, allocatable :: web_transp_on(:), transp_on(:)
      integer, save, allocatable :: nacsc(:), mru2chan(:)
      real, save :: c_can_depth, basin_area
      real, save, allocatable :: covden_sum(:),covden_win(:)
!      real, save, allocatable :: s_ohoriz_depth(:), s_root_depth(:)
      real, save, allocatable :: s_ohoriz_depth(:)
      real, save, allocatable :: s_porosity(:), s_theta_wp(:)
      real, save, allocatable :: s_theta_fc(:), uz_depth(:,:), ac(:,:)
      real, save, allocatable :: qdffrac(:), mru_area(:)
      real, save, allocatable :: mru_area_frac(:),chan_area(:)
      real, save, allocatable :: ar_fill(:,:)
      double precision, save, allocatable :: can_ohoriz_vol(:)
      double precision, save, allocatable :: uz2sat_vol(:,:)

c Mixing variables
c      double precision, save :: vmix_imp(:,:)
      double precision, save :: vmix_basin0
      double precision, save, allocatable :: vmix_can(:,:)
      double precision, save, allocatable :: vmin_canopy(:)
      double precision, save, allocatable :: vmix_snow(:,:)
      double precision, save, allocatable :: vmix_ohoriz(:,:)
c$$$      double precision, save, allocatable :: vmix_rz(MAXMNR_3D)  ! possible root zone for later
      double precision, save, allocatable :: vmix_uz(:,:,:)  ! combination of root zone and suz
      double precision, save, allocatable :: vmix_uz2can(:,:)  ! UZ to canopy
      double precision, save, allocatable :: vmix_uz2qdf(:,:)  ! Direct flow
      double precision, save, allocatable :: vmix_uz2sat(:,:)  ! Recharge water
      double precision, save, allocatable :: vmix_sat2uz(:)  ! Water transferred from sat to uz
      double precision, save, allocatable :: vmix_uzgen(:,:)
      double precision, save, allocatable :: vmix_uzrip(:,:)
      double precision, save, allocatable :: vmix_uzup(:,:)
      double precision, save, allocatable :: vmix_qdf(:,:)
      double precision, save, allocatable :: vmix_sat(:,:)
      double precision, save, allocatable :: vmix_well(:)
      double precision, save, allocatable :: vmix_satpref(:,:)
      double precision, save, allocatable :: vmix_hill(:,:)
      double precision, save, allocatable :: vmix_hillexp(:,:)
      double precision, save, allocatable :: vmix_stream(:)
      double precision, save, allocatable :: vmix_diversion(:)
      double precision, save, allocatable :: vmix_chan_loss(:)
      double precision, save, allocatable :: vmix_mru(:,:), vmix_mru0(:)
      double precision, save, allocatable :: vmix_basin(:)
c
c basin vars
c
      real, save :: basin_sto_cm,basin_in_cm,basin_out_cm,basin_et_cm
      real, save :: basin_ppt_cm,basin_irr_ext_cm
      real, save :: basin_irr_sat_cm,basin_irr_hyd_cm
      real, save :: basin_intcp_sto_cm, basin_intcp_cm
      real, save :: basin_intcp_evap_cm
      real, save :: basin_transp_cm, basin_thruf_cm,basin_net_dep_cm
      real, save :: basin_net_rain_cm,basin_net_snow_cm,basin_pweqv_cm
      real, save :: basin_snowmelt_cm,basin_snowevap_cm,basin_surfdep_cm
      real, save :: basin_soil_moist_cm,basin_sssto_cm,basin_ofhort_cm
      real, save :: basin_ofdunn_cm,basin_sroff_cm,basin_infil_cm
      real, save :: basin_qdf_cm
      real, save :: basin_uz_et_cm,basin_gw_sto_cm,basin_gw1_in_cm
      real, save :: basin_gw2_in_cm,basin_qwet_cm
      real, save :: basin_vpref_cm, basin_recharge_cm
      real, save :: basin_uz2sat_cm,basin_sat2uz_cm
      real, save :: basin_qwell_cm,basin_qpref_cm,basin_gwloss_cm
      real, save :: basin_gwflow_cm,basin_chan_sto_cm
      real, save :: basin_stflow_cm,basin_chan_div_cm,basin_chan_loss_cm
      real, save :: basin_qsim_cm, basin_qobs_cm
      real, save :: basin_qsim_m3s,basin_qobs_m3s
      real, save, allocatable :: runoff_m3s(:)
!
! Local variables in init
!
      real, save, allocatable :: maxchanloss(:), intcp_on(:)
!      real, save, allocatable :: s_theta_0(:),s_rock_depth(:)
      real, save, allocatable :: wei(:), q(:,:)
      real, save, allocatable :: s_satpref_zmax(:), s_satpref_zmin(:)
      real, save, allocatable :: sr0(:),sbar0(:)
!
! Local variables and getvars in run
!
      real, save, allocatable :: mru_ppt(:), mru_dep(:)
      real, save, allocatable :: net_dep(:)
      real, save, allocatable :: net_rain(:), net_snow(:)
c      real, save, allocatable :: irrig_sat_mru(:)
c      real, save, allocatable :: irrig_hyd_mru(:)
c      real, save, allocatable :: irrig_ext_mru(:)
      real, save, allocatable :: intcp_stor(:), intcp_evap(:)
      real, save, allocatable :: pkwater_equiv(:), psoilmru(:)
      real, save, allocatable :: snowmelt(:), snow_evap(:)
      real, save, allocatable :: srz_sc(:), suz_sc(:)
      real, save, allocatable :: sae(:)
      real, save, allocatable :: sae_local(:,:)
      real, save, allocatable :: sd(:,:)
!      real, save, allocatable :: z_wt_local(:,:)
      real, save, allocatable :: suz(:,:),srz(:,:)
!      real, save, allocatable :: srzwet(:,:)
      real, save, allocatable :: qb(:)
      real, save, allocatable :: qdf(:), gw_loss(:)
      real, save, allocatable :: qpref(:), acm(:), afx(:)
      real, save, allocatable :: qof(:),qofs(:), quz(:)
      real, save, allocatable :: rex(:)
      real, save, allocatable :: quz_local(:,:),qdf_local(:,:)
      real, save, allocatable :: uz_infil(:,:)
      real, save, allocatable :: uz2sat(:,:), qvpref(:)
c      double precision, save, allocatable :: irrig_hyd_seg(:)
      real, save, allocatable :: gw_in1(:), gw_in2(:)
      real, save, allocatable :: chan_loss(:)
      
! Construction of volume file names (v_*)    
      !integer, save :: isoh1_len, isoh2_len, isogl_len, isogs_len, sol_h2_len
      !character*60, save :: isogs, isogl
      !character*256 :: aline
      !
      !character*3000, save :: sol_header1, sol_header2, iso_header1,iso_header2
      logical filflg

      TYPE :: outfiles   ! file names, shortnames, and logical unit numbers for input and output files.
         character(60) :: file   ! Output file
         integer       :: lun        ! integer returned by NEWUNIT
      END TYPE outfiles
!
      TYPE(outfiles),save :: vf_bas, vf_hyd
      TYPE(outfiles),save,allocatable :: vf_mru(:), vf_uzgen(:),
     $ vf_uzrip(:), vf_uzup(:), vf_can(:), vf_snow(:), 
     $ vf_transp(:), vf_ohoriz(:), vf_uz(:,:), vf_qdf(:), vf_sat(:),
     $ vf_satpref(:), vf_hill(:), vf_uz2sat(:), vf_hillexp(:)
      !TYPE(outfiles),save, allocatable :: vf_imperv(:)

!
      character(3000), save :: out_dir, hdr, filename, mruid, nacid
      integer :: tmplun, path_len, j

      END MODULE WEBMOD_RESMOD
      
c***********************************************************************
c
c     Main webmod_res routine
c

      integer function webmod_res(arg)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: Arg
      CHARACTER(LEN=256) SVN_ID

      integer, external :: webrdecl,webrinit, webrrun

      SVN_ID =
     $     '$Id: webmod_res.f 40 2009-04-24 rmwebb $ '
      
      webmod_res = 0

      if(arg.eq.'declare') then
        webmod_res = webrdecl()
      else if(arg.eq.'initialize') then
        webmod_res = webrinit()
      else if(arg.eq.'run') then
        webmod_res = webrrun()
c
c     No cleanup routine used
c
       end if

      END FUNCTION webmod_res


c***********************************************************************
c 
c     webrdecl - declare variables and parameters for WEBMOD_RES
c

      integer function webrdecl()
      
      USE WEBMOD_RESMOD
      IMPLICIT NONE

      integer check_uz_dim

      webrdecl = 1

!
! Get dimensions
!



      nmru = getdim('nmru')
      if(nmru.eq.-1) return
      nac = getdim('nac')
      if(nac.eq.-1) return
      nresinp = getdim('nresinp')
      if(nresinp.eq.-1) return
      nobs = getdim('nobs')
      if(nobs.eq.-1) return
      nchan = getdim('nchan')
      if(nchan.eq.-1) return
      nhydro = getdim('nhydro')
      if(nhydro.eq.-1) return
      nsolute = getdim('nsolute')
      if(nsolute.eq.-1) return
!      nhcs = get*dim('nhcs')
      nchemvar = getdim('nchemvar')
      if(nchemvar.eq.-1) return
!      nmru_res = get*dim('nmru_res')
!      nirrig_ext = getdim('nirrig_ext')
!      if(nirrig_ext.eq.-1) return
!      nirrig_int = getdim('nirrig_int')
!      if(nirrig_int.eq.-1) return
      nchem_ext = getdim('nchem_ext')
      if(nchem_ext.eq.-1) return
      nchemobs = getdim('nchemobs')
      if(nchemobs.eq.-1) return

c number of point concentrations

      nchemdat = 1+nchem_ext+nchemobs
c
c The following check must be made to ensure that the value representing the 
c maximum matrix size for the three dimensional public variable, vmix_uz, is
c correct. If it isn't, then geochemical fluxes will not be simulated correctly.
c
      check_uz_dim = getdim('nac_nmru_nresinp')
      if(check_uz_dim.eq.-1) return

      if(check_uz_dim.ne.nac*nmru*nresinp) then
         print*,'The dimension size for nac_nmru_nresinp must ',
     $        'be set to ',nac*nmru*nresinp,
     $        ' for solute fluxes to be computed correctly'
         return
      end if
c
c Volume mixing matrices. A maximum of six possible mixing
c volumes for inputs will be use. This can be change by altering
c nresinp in setdim_web.f and nresinp in fmodules_web.inc (to
c be copied to fmodules.inc). This represents the maximum unique
c reservoir inputs for a new reservoir, specifically, the new 
c root_zone composition may be determined by the compositions of
c the old root_zone, precipitation, snowmelt, throughfall, and ET.
c
c Impermeable areas (precip); Use later once impermeable areas
c are included explicitly (i.e. TR55 or something)
c
c$$$      if(decl*var('webr', 'vmix_imp', 'nmru,nresinp', nmru*nresinp,
c$$$     $     'double', 'volumes to mix for making new solution for '//
c$$$     $     'impermeable areas', 'm3',vmix_imp).ne.0) return
c
c Canopy - (precip, ET, root_zone)
c
      allocate(vmix_can(nmru,nresinp))
      if(declvar('webr', 'vmix_can', 'nmru,nresinp', nmru*nresinp,
     $     'double', 'volumes to mix for making new canopy solution',
     + 'm3',vmix_can).ne.0) return
c
c Snowpack - (precip, canopy, ET)
c
      allocate(vmix_snow(nmru,nresinp))
      if(declvar('webr', 'vmix_snow', 'nmru,nresinp', nmru*nresinp,
     $ 'double', 'volumes to mix for making new snowpack solution',
     + 'm3',vmix_snow).ne.0) return
c
c O-horizon - (precip)
c
      allocate(vmix_ohoriz(nmru,nresinp))
      if(declvar('webr', 'vmix_ohoriz', 'nmru,nresinp', 
     $     nmru*nresinp, 'double', 'volumes to mix for making '//
     $     'new o-horizon solution', 'm3',vmix_ohoriz).ne.0) return
c
c Root zone and unsaturated zone to be combined for single UZ chemistry
c
c Root Zone (precip, canopy, ET, snowmelt, sat_zone)
c
c Might need this later as a separate reservoir
c
c$$$      if(decl*var('webr', 'vmix_rz', 'nmru_nac_nresinp',
c$$$     $  MAXMNR_3D, 'double', 'volumes to mix for making '//
c$$$     $  'new root zone solution', 'm3',vmix_rz).ne.0) return
c
c Unsaturated zone (precip, canopy, ET, snowmelt, sat_zon))
c
      allocate(vmix_uz(nac,nmru,nresinp))
      if(declvar('webr', 'vmix_uz', 'nac_nmru_nresinp',
     $     nac*nmru*nresinp, 'double', 'volumes to mix for making '//
     $     'new unsaturated zone storage, suz, solution',
     $     'm3',vmix_uz).ne.0) return
c
c Unsaturated zone mixing box for feeding canopy transpiration
c
      allocate(vmix_uz2can(nac,nmru))
      if(declvar('webr', 'vmix_uz2can', 'nac,nmru',
     $     nac*nmru, 'double', 'volumes to mix for making '//
     $     'water transpired to canopy',
     $     'm3',vmix_uz2can).ne.0) return
c
c Unsaturated zone mixing box for direct flow
      allocate(vmix_uz2qdf(nac,nmru))
      if(declvar('webr', 'vmix_uz2qdf', 'nac,nmru',
     $     nac*nmru, 'double', 'volumes to mix for direct flow',
     $     'm3',vmix_uz2qdf).ne.0) return
c
c Unsaturated zone mixing box for recharge
c
      allocate(vmix_uz2sat(nac,nmru))
      if(declvar('webr', 'vmix_uz2sat', 'nac,nmru',
     $     nac*nmru, 'double', 'volumes to mix for recharge',
     $     'm3',vmix_uz2sat).ne.0) return
c
c Flux from saturated zone for irrigating the land surface on the
c next time step
c
      allocate(vmix_well(nmru))
      if(declvar('webr', 'vmix_well', 'nmru',
     $     nmru, 'double', 'volume of saturated zone water '//
     $     'pumped to the surface for irrigation',
     $     'm3',vmix_well).ne.0) return
c
c Flux from saturated zone to uz for root zone wetting and stranded
c pore water.
c
      allocate(vmix_sat2uz(nmru))
      if(declvar('webr', 'vmix_sat2uz', 'nmru',
     $     nmru, 'double', 'volume transferred from saturated '//
     $     'zone to the unsaturated zone',
     $     'm3',vmix_sat2uz).ne.0) return
c
c Net flux (uz2sat - srzwet) where uz2sat is the transfer of
c water from the unsaturated zone to the saturated zone resulting
c from a change in water table, and srzwet is the flux of water
c from the saturated zone to the unsaturated zone that occurs where
c artesian conditions replenish moisture deficits in the root zone.
c
      allocate(uz2sat_vol(nac,nmru))
      if(declvar('webr', 'uz2sat_vol', 'nac,nmru',
     $     nac*nmru, 'double', 'Net volume of water transferred'//
     $     'from the unsaturated zone to the saturated zone.'//
     $     'resulting from a changing water table below the root zone',
     $     'm3',uz2sat_vol).ne.0) return
c
c Unsaturated zone composite - (suz, combined rz/suz)
c
      allocate(vmix_uzgen(nmru,nresinp))
      if(declvar('webr', 'vmix_uzgen', 'nmru,nresinp', nmru*nresinp,
     $ 'double', 'volumes to track to simulate composite solution '//
     $ 'of individual uz bins', 'm3',vmix_uzgen).ne.0) return
c
c Unsaturated zone composite for riparian zone - (suz, combined rz/suz)
c
      allocate(vmix_uzrip(nmru,nresinp))
      if(declvar('webr', 'vmix_uzrip', 'nmru,nresinp', nmru*nresinp,
     $ 'double', 'volumes to track to simulate composite solution '//
     $ 'of individual uz bins in riparian zone', 'm3',
     $ vmix_uzrip).ne.0) return
c
c Unsaturated zone composite for uplands - (suz, combined rz/suz)
c
      allocate(vmix_uzup(nmru,nresinp))
      if(declvar('webr', 'vmix_uzup', 'nmru,nresinp', nmru*nresinp,
     $ 'double', 'volumes to track to simulate composite solution '//
     $ 'of individual uz bins in uplands', 'm3',vmix_uzup).ne.0) return
c
c Unsaturated zone preferential flow - (suz,  combined rz/suz)
c
      allocate(vmix_qdf(nmru,nresinp))
      if(declvar('webr', 'vmix_qdf', 'nmru,nresinp', nmru*nresinp,
     $ 'double', 'volumes to mix for making new direct flow solution',
     $ 'm3',vmix_qdf).ne.0) return
c
c Saturated zone - (sat)
c
      allocate(vmix_sat(nmru,nresinp))
      if(declvar('webr', 'vmix_sat', 'nmru,nresinp', nmru*nresinp,
     $ 'double','volumes to mix for making new saturated zone '//
     $ 'solution', 'm3',vmix_sat).ne.0) return
c
c Saturated zone preferential flow - (sat zone)
c
      allocate(vmix_satpref(nmru,nresinp))
      if(declvar('webr', 'vmix_satpref', 'nmru,nresinp', 
     $     nmru*nresinp, 'double', 
     $     'volumes to mix for making new deep preferential '//
     $     'flow solution', 'm3',vmix_satpref).ne.0) return
c
c Hillslope -(combination of overland flow, direct flow, satpref,
c    baseflow)
c
      allocate(vmix_hill(nmru,nresinp))
      if(declvar('webr', 'vmix_hill', 'nmru,nresinp', 
     $     nmru*nresinp, 'double', 
     $     'volumes to mix for making new stream solution',
     +     'm3',vmix_hill).ne.0) return
c
c MRU inputs and outputs
c
      allocate(vmix_mru(nmru,nresinp))
      if(declvar('webr', 'vmix_mru', 'nmru,nresinp', 
     $     nmru*nresinp, 'double', 
     $     'volumes to track inputs and outputs for MRU',
     +     'm3',vmix_mru).ne.0) return
c
c Hillslope to stream mixing matrix -(Distributes hillslope exports
c to stream segments.
c
      allocate(vmix_hillexp(nhydro,nmru))
      if(declvar('webr', 'vmix_hillexp', 'nhydro,nmru', 
     $     nhydro*nmru, 'double', 
     $     'volumes to mix for making new stream solutions',
     +     'm3',vmix_hillexp).ne.0) return
c
c Streamflow -(Stream segments - segment 1, which was located closest
c to the outlet, is exported from the basin at the end of each time step)
c
      allocate(vmix_stream(nhydro))
      if(declvar('webr', 'vmix_stream', 'nhydro', 
     $     nhydro, 'double', 
     $     'volume in each stream segment',
     +     'm3',vmix_stream).ne.0) return
c
c Streamflow diversions to be applied as irrigation
c
      allocate(vmix_diversion(nhydro))
      if(declvar('webr', 'vmix_diversion', 'nhydro', 
     $     nhydro, 'double', 
     $     'volume of stream water pumped or diverted for'//
     $     'irrigation',
     +     'm3',vmix_diversion).ne.0) return

      allocate(vmix_chan_loss(nhydro))
      if(declvar('webr', 'vmix_chan_loss', 'nhydro', 
     $     nhydro, 'double', 
     $     'volume of stream water lost by seepage through '//
     $     'channel bed',
     +     'm3',vmix_chan_loss).ne.0) return
c
c Basin -(precipitation in. ET and streamflow out)
c
      allocate(vmix_basin(nresinp))
      if(declvar('webr', 'vmix_basin', 'nresinp', 
     $     nresinp, 'double', 
     $     'volumes to check overall water balance for basin.',
     +     'm3',vmix_basin).ne.0) return

c
c Declare basin variables needed for summary module
c
      if(declvar('webr', 'basin_sto_cm', 'one', 1, 'real',
     + 'Storage in basin including streams, groundwater, subsurface '//
     + 'storage, soil moisture, snowpack, and interception',
     + 'cm',
     + basin_sto_cm).ne.0) return

      if(declvar('webr', 'basin_in_cm', 'one', 1, 'real',
     + 'Basin inputs from precip, irrigation, and regional groundwater',
     + 'cm',
     +  basin_in_cm).ne.0) return

      if(declvar('webr', 'basin_out_cm', 'one', 1, 'real',
     + 'Basin discharge and losses to deep aquifers',
     + 'cm',
     +  basin_out_cm).ne.0) return

      if(declvar('webr', 'basin_et_cm', 'one', 1, 'real',
     + 'Basin evaporatranspiration',
     + 'cm',
     +  basin_et_cm).ne.0) return

      if(declvar('webr', 'basin_ppt_cm', 'one', 1, 'real',
     + 'Basin precipitation',
     + 'cm',
     +  basin_ppt_cm).ne.0) return

      if(declvar('webr', 'basin_irr_ext_cm', 'one', 1, 'real',
     + 'Basin irrigation from external sources (canals and deep wells)',
     + 'cm',
     +  basin_irr_ext_cm).ne.0) return

      if(declvar('webr', 'basin_irr_sat_cm', 'one', 1, 'real',
     + 'Basin irrigation from shallow wells',
     + 'cm',
     +  basin_irr_sat_cm).ne.0) return

      if(declvar('webr', 'basin_irr_hyd_cm', 'one', 1, 'real',
     + 'Basin irrigation from stream diversions',
     + 'cm',
     +  basin_irr_hyd_cm).ne.0) return

      if(declvar('webr', 'basin_intcp_sto_cm', 'one', 1, 'real',
     + 'Basin storage on leaves/canopy',
     + 'cm',
     +  basin_intcp_sto_cm).ne.0) return

      if(declvar('webr', 'basin_intcp_cm', 'one', 1, 'real',
     + 'Basin average interception',
     + 'cm',
     +  basin_intcp_cm).ne.0) return

      if(declvar('webr', 'basin_intcp_evap_cm', 'one', 1, 'real',
     + 'Basin average evaporation from leaves/canopy',
     + 'cm',
     +  basin_intcp_evap_cm).ne.0) return

      if(declvar('webr', 'basin_transp_cm', 'one', 1, 'real', 
     + 'Basin average transpiration',
     $ 'cm',basin_transp_cm).ne.0) return

      if(declvar('webr', 'basin_thruf_cm', 'one', 1, 'real',
     + 'Basin average canopy throughfall',
     + 'cm',
     +  basin_thruf_cm).ne.0) return

      if(declvar('webr', 'basin_net_dep_cm', 'one', 1, 'real',
     + 'Basin net deposition (direct plus throughfall)',
     + 'cm',
     +  basin_net_dep_cm).ne.0) return

      if(declvar('webr', 'basin_net_rain_cm', 'one', 1, 'real',
     + 'Basin net rain',
     + 'cm',
     +  basin_net_rain_cm).ne.0) return

      if(declvar('webr', 'basin_net_snow_cm', 'one', 1, 'real',
     + 'Basin net snow',
     + 'cm',
     +  basin_net_snow_cm).ne.0) return

      if(declvar('webr', 'basin_pweqv_cm', 'one', 1, 'real',
     + 'Basin packwater equivalence',
     + 'cm',
     +  basin_pweqv_cm).ne.0) return

      if(declvar('webr', 'basin_snowmelt_cm', 'one', 1, 'real',
     + 'Basin snowmelt',
     + 'cm',
     +  basin_snowmelt_cm).ne.0) return

      if(declvar('webr', 'basin_snowevap_cm', 'one', 1, 'real',
     + 'Basin snow evaporation',
     + 'cm',
     +  basin_snowevap_cm).ne.0) return

      if(declvar('webr', 'basin_surfdep_cm', 'one', 1, 'real',
     + 'Basin hillslope deposition (direct plus melt)',
     + 'cm',
     +  basin_surfdep_cm).ne.0) return

      if(declvar('webr', 'basin_soil_moist_cm', 'one', 1, 'real',
     + 'Basin soil moisture in root zone',
     + 'cm',
     +  basin_soil_moist_cm).ne.0) return

      if(declvar('webr', 'basin_sssto_cm', 'one', 1, 'real',
     + 'Basin subsurface storage',
     + 'cm',
     +  basin_sssto_cm).ne.0) return

      if(declvar('webr', 'basin_ofhort_cm', 'one', 1, 'real',
     + 'Basin infiltration excess (Hortonian overland flow)',
     + 'cm',
     +  basin_ofhort_cm).ne.0) return

      if(declvar('webr', 'basin_ofdunn_cm', 'one', 1, 'real',
     + 'Basin saturation excess (Dunnian overland flow)',
     + 'cm',
     +  basin_ofdunn_cm).ne.0) return

      if(declvar('webr', 'basin_sroff_cm', 'one', 1, 'real',
     + 'Basin overland flow',
     + 'cm',
     +  basin_sroff_cm).ne.0) return

      if(declvar('webr', 'basin_infil_cm', 'one', 1, 'real',
     + 'Basin infiltration',
     + 'cm',
     +  basin_infil_cm).ne.0) return

      if(declvar('webr', 'basin_qdf_cm', 'one', 1, 'real',
     + 'Basin lateral macropore flow',
     + 'cm',
     +  basin_qdf_cm).ne.0) return

      if(declvar('webr', 'basin_uz_et_cm', 'one', 1, 'real', 
     + 'Basin evaporation from the soil',
     $ 'cm',basin_uz_et_cm).ne.0) return

      if(declvar('webr', 'basin_gw_sto_cm', 'one', 1, 'real',
     + 'Basin groundwater storage',
     + 'cm',
     +  basin_gw_sto_cm).ne.0) return

      if(declvar('webr', 'basin_gw1_in_cm', 'one', 1, 'real',
     + 'Basin inputs from leaky canals',
     + 'cm',
     +  basin_gw1_in_cm).ne.0) return

      if(declvar('webr', 'basin_gw2_in_cm', 'one', 1, 'real',
     + 'Basin inputs from upgradient groundwater',
     + 'cm',
     +  basin_gw2_in_cm).ne.0) return

      if(declvar('webr', 'basin_qwet_cm', 'one', 1, 'real',
     + 'Basin root-zone replenishment from groundwater',
     + 'cm',
     +  basin_qwet_cm).ne.0) return

      if(declvar('webr', 'basin_vpref_cm', 'one', 1, 'real',
     + 'Basin infiltration delivered directly to groundwater',
     + 'cm',
     +  basin_vpref_cm).ne.0) return

      if(declvar('webr', 'basin_recharge_cm', 'one', 1, 'real',
     + 'Basin recharge to groundwater not including qvpref',
     + 'cm',
     +  basin_recharge_cm).ne.0) return

      if(declvar('webr', 'basin_uz2sat_cm', 'one', 1, 'real',
     + 'Basin UZ water engulfed into saturated zone as '//
     + 'water table rises','cm',
     +  basin_uz2sat_cm).ne.0) return

      if(declvar('webr', 'basin_sat2uz_cm', 'one', 1, 'real',
     + 'Basin water left in UZ as water table lowers',
     + 'cm',
     +  basin_sat2uz_cm).ne.0) return

      if(declvar('webr', 'basin_qwell_cm', 'one', 1, 'real',
     + 'Basin pumping from shallow wells',
     + 'cm',
     +  basin_qwell_cm).ne.0) return

      if(declvar('webr', 'basin_qpref_cm', 'one', 1, 'real',
     + 'Basin preferential flow in the saturated zone',
     + 'cm',
     +  basin_qpref_cm).ne.0) return

      if(declvar('webr', 'basin_gwloss_cm', 'one', 1, 'real',
     + 'Basin groundwater loss to deep aquifer',
     + 'cm',
     +  basin_gwloss_cm).ne.0) return

      if(declvar('webr', 'basin_gwflow_cm', 'one', 1, 'real',
     + 'Basin baseflow',
     + 'cm',
     +  basin_gwflow_cm).ne.0) return

      if(declvar('webr', 'basin_chan_sto_cm', 'one', 1, 'real',
     + 'Basin channel storage',
     + 'cm',
     +  basin_chan_sto_cm).ne.0) return

      if(declvar('webr', 'basin_stflow_cm', 'one', 1, 'real',
     + 'Basin hillslope runoff to channels',
     + 'cm',
     +  basin_stflow_cm).ne.0) return

      if(declvar('webr', 'basin_chan_div_cm', 'one', 1, 'real',
     + 'Basin channel diversions for irrigation',
     + 'cm',
     +  basin_chan_div_cm).ne.0) return

      if(declvar('webr', 'basin_chan_loss_cm', 'one', 1, 'real',
     + 'Basin channel losses to deep aquifer',
     + 'cm',
     +  basin_chan_loss_cm).ne.0) return

      if(declvar('webr', 'basin_qsim_cm', 'one', 1, 'real',
     + 'Basin simulated discharge',
     + 'cm',
     +  basin_qsim_cm).ne.0) return

      if(declvar('webr', 'basin_qsim_m3s', 'one', 1, 'real',
     + 'Basin simulated discharge',
     + 'm^3/s',
     +  basin_qsim_m3s).ne.0) return
c
c for clarity, use the nobs dimensioned runoff variable to
c create a nobs dimensioned runoff_m3s and the derive a single
c observed volume and equivalent depth of discharge to be passed
c to the summary module. No variable declared of nobs=0
c
      if(nobs.ne.0) then
        allocate(runoff_m3s(nobs))
        if(declvar('webr', 'runoff_m3s', 'nobs', nobs, 'real',
     +   'Observed discharge at one of the nobs stations',
     +   'm3s',
     +    runoff_m3s).ne.0) return
      end if

      if(declvar('webr', 'basin_qobs_cm', 'one', 1, 'real',
     + 'Basin observed discharge (at station qobsta)',
     + 'cm',
     +  basin_qobs_cm).ne.0) return

      if(declvar('webr', 'basin_qobs_m3s', 'one', 1, 'real',
     + 'Basin observed discharge (at station qobsta)',
     + 'm^3/s',
     +  basin_qobs_m3s).ne.0) return

c
c declare parameters need from other modules
c$$$c
c$$$c io_chem - output file number - Use this in phreeq_mms, not here
c$$$c
c$$$      if(decl*param('io', 'chemout_file_unit', 'one', 'integer',
c$$$     +   '90', '50', '99',
c$$$     +   'Unit number for solute chemistry output file',
c$$$     +   'Unit number for solute chemistry output file',
c$$$     +   'integer').ne.0) return
c
c from basin_topg
c
      if(declparam('basin', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

c      if(declparam('topc', 'qobsta', 'one', 'integer',
c     +   '1', 'bounded', 'nobs',
c     +   'Index of streamflow station for calculating '//
c     +   'objective function.','Index of streamflow station '//
c     +   'for calculating objective function.','none').ne.0) return

      allocate(mru_area(nmru))
      if(declparam('basin', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return
c
c potet_hamon_prms - no parameters needed
c
c intcp_prms
c
      allocate(covden_sum(nmru))
      if(declparam('intcp', 'covden_sum', 'nmru', 'real',
     +   '.5', '0.', '1.0',
     +   'Summer vegetation cover density for major vegetation type',
     +   'Summer vegetation cover density for the major '//
     +   'vegetation type on each MRU',
     +   'decimal percent')
     +   .ne.0) return 

      allocate(covden_win(nmru))
      if(declparam('intcp', 'covden_win', 'nmru', 'real',
     +   '.5', '0.', '1.0',
     +   'Winter vegetation cover density for major vegetation type',
     +   'Winter vegetation cover density for the major '//
     +   'vegetation type on each MRU',
     +   'decimal percent')
     +   .ne.0) return

c
c nwsmelt_topg
c
      allocate(WEI(nmru))
      if(declparam('snow', 'WEI', 'nmru', 'real',
     &    '0.', '0.', '1000.',
     &    'Initial water equivalent for each MRU',
     &    'Initial water equivalent for each MRU',
     &    'inches').ne.0) return
c
c topmodg_chem
c
      allocate(s_porosity(nmru))
      if(declparam('topc','s_porosity', 'nmru',
     +   'real', '0.4', '0.001', '0.8',
     +   'Soil porosity', 'Effective soil porosity, equal '//
     +   'to saturated soil moisture content.',
     +   'cm3/cm3') .ne.0) return

      allocate(s_theta_fc(nmru))
      if(declparam('topc','s_theta_fc', 'nmru',
     +   'real', '0.23', '0.001', '0.7',
     +   'Volumetric soil moisture content at field capacity',
     +   'Volumetric soil moisture content at field capacity. '//
     +   'Field capacity is determined as the moisture content '//
     +   'at which the hydraulic conductivity is equal to '//
     +   '1E-8 cm/s.', 'cm3/cm3') .ne.0) return

      allocate(s_theta_wp(nmru))
      if(declparam('topc','s_theta_wp', 'nmru',
     +   'real', '0.13', '0.001', '0.56',
     +   'Volumetric soil moisture content at wilting point',
     +   'Volumetric soil moisture content at wilting point. '//
     +   'The wilting point is determined as the mositure content '//
     +   'at a tension of 15,300 cm (15 bars). Also known as '//
     +   'residual soil moisture content.',
     +   'cm3/cm3') .ne.0) return


!      allocate(s_root_depth(nmru))
!      if(declparam('topc', 's_root_depth', 'nmru', 'real',
!     +   '1.8', '0.1', '100',
!     +   'Rooting depth.','Rooting depth from ground surface, '//
!     +   'Available water capacity (moisture content at field '//
!     +   'capacity minus that at wilting point) * root_depth '//
!     +   'equals the maximum soil moisture deficit, srmax: '//
!     +   'smcont_sc. smcont_sc = (srmax - srz)/root_depth.'//
!     +   'Be sure to set root_depth to a value less then '//
!     +   'the depth to bedrock, s_rock_depth',
!     +   'm') .ne.0) return

!      allocate(s_rock_depth(nmru))
!      if(declparam('topc', 's_rock_depth', 'nmru', 'real',
!     +     '6.0', '0.1', '300.0',
!     +     'Average depth to bedrock for the MRU.','Average depth to '//
!     $     'bedrock. Must be greater than the rooting depth, '//
!     $     's_rock_depth.',
!     +     'm') .ne.0) return

!      allocate(s_theta_0(nmru))
!      if(declparam('topc','s_theta_0', 'nmru',
!     +   'real', '0.23', '0.01', '0.7',
!     +   'Initial volumetric soil moisture content in the root zone.',
!     +   'Initial volumetric soil moisture content in the root zone.',
!     +   'cm3/cm3') .ne.0) return
            
      allocate(sbar0(nmru))
      if(declparam('topc', 'sbar0', 'nmru', 'real',
     +   '.001', '0', '10',  
     +   'Initial soil moisture deficit, SBAR, in MRU.',
     +   'Initial soil moisture deficit, SBAR, in MRU.',
     +   'm').ne.0) return

      allocate(qdffrac(nmru))
      if(declparam('topc', 'qdffrac', 'nmru', 'real',
     +    '.3', '0', '1',
     +   'Fraction of recharge that runs off as direct flow.',
     +   'Fraction of recharge that runs off as direct flow. '//
     +   'QDF=QDFFRAC*QUZ','Fraction')
     +    .ne.0) return

c     +   'Fraction of infiltration that becomes'//
c     +   ' lateral direct flow.','Fraction of infiltration'//
c     +   ' that becomes lateral direct flow.'//
c     +   'QDF=QDFFRAC*INFILTRATION','Proportion')
c     +    .ne.0) return
      
      allocate(nacsc(nmru))
      if(declparam('topc', 'nacsc', 'nmru', 'integer',
     +   '1', '0', '100',
     +   'Number of TTI bins in the subcatchment.',
     +   'Number of TTI bins in the subcatchment.',
     +   'none').ne.0) return

      allocate(ac(nac,nmru))
      if(declparam('topc', 'ac', 'nac,nmru', 'real',
     +   '1', '0', '1',
     +   'Fractional area for each TTI bin.',
     +   'Fractional area for each TTI bins',
     +   'km2/km2').ne.0) return      

      allocate(s_satpref_zmin(nmru))
      if(declparam('topc', 's_satpref_zmin', 'nmru', 'real',
     +   '-5.0', '-100.0', '0.0',
     +   'Water table elevation at which preferential flow in the '//
     +   'saturated zone begins.','Water table elevation (z=0 at '//
     $   'surface) at which preferential flow in the saturated '//
     $   'zone begins. Set equal to s_satpref_zmax if path does '//
     $   'not exist.', 'm') .ne.0) return

      allocate(s_satpref_zmax(nmru))
      if(declparam('topc', 's_satpref_zmax', 'nmru', 'real',
     +   '-5.0', '-100.0', '0.0',
     +   'Water table elevation at which preferential flow in the '//
     +   'saturated zone reaches a maximum.','Water table elevation '//
     +   '(z=0 at surface) at which preferential flow in the '//
     $   'saturated zone reaches a maximum. Set equal to '//
     $   's_satpref_zmin if path does not exist.', 'm') .ne.0) return

      allocate(mru_area_frac(nmru))
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

      allocate(mru2chan(nmru))
      if(declparam('top2c', 'mru2chan', 'nmru', 'integer',
     +   '1', 'bounded', 'nchan',
     +   'Index of channel receiving discharge from MRU',
     +   'Index of channel receiving discharge from MRU',
     +   'none').ne.0) return

c
c The following additional parameters, enable volume calculations and
c soil moisture calculations needed to track water table heights, 
c reservoir volumes and solute fluxes (to be computed by the 
c phreeq_mms module).
c
c
c srmax will be calculated as (s_theta_fc - s_theta_wp) * s_root_depth
c
c$$$      if(decl*pri('webr_srmax', 'srmax', nmru, 'real', srmax)
c$$$     + .ne.0) return

      if(declparam('webres', 'c_can_depth', 'one', 'real',
     +   '.001', '0.00001', '0.01',
     +     'Fixed depth for residual water on the canopy',
     $     'Fixed depth for residual water on the canopy. '//
     $     'Residual water will be tranferred to and from the '//
     $     'O-horizon on days of leaves off or leaves on, '//
     $     'respectively.', 'm') .ne.0) return

      allocate(s_ohoriz_depth(nmru))
      if(declparam('webres', 's_ohoriz_depth', 'nmru', 'real',
     +   '.005', '0.0001', '0.5',
     +   'Fixed depth for O-horizon','Fixed depth for O-horizon. '//
     +   'Carbon produced in the the O-horizon reservoir will be '//
     +   'available for transport by overland flow.',
     +   'm') .ne.0) return

! locals from init
      allocate (vmix_mru0(nmru))
      allocate (transp_on(nmru))
      allocate (web_transp_on(nmru))
      allocate (sr0(nmru))
      allocate (intcp_on(nmru))
      allocate (can_ohoriz_vol(nmru))
      allocate (vmin_canopy(nmru))
      allocate (uz_depth(nac,nmru))
      allocate (q(nhydro,nchan))
      allocate (ar_fill(nhydro,nchan))
      allocate (maxchanloss(nhydro))
      allocate (chan_area(nchan))

! locals from run
      allocate (mru_ppt(nmru))
      allocate (mru_dep(nmru))
      allocate (net_dep(nmru))
      allocate (net_rain(nmru))
      allocate (net_snow(nmru))
c      allocate (irrig_sat_mru(nmru))
c      allocate (irrig_hyd_mru(nmru))
c      allocate (irrig_ext_mru(nmru))
      allocate (intcp_stor(nmru))
      allocate (intcp_evap(nmru))
      allocate (pkwater_equiv(nmru))
      allocate (psoilmru(nmru))
      allocate (snowmelt(nmru))
      allocate (snow_evap(nmru))
      allocate (srz_sc(nmru))
      allocate (suz_sc(nmru))
      allocate (sae(nmru))
      allocate (sae_local(nac,nmru))
      allocate (sd(nac,nmru))
!      allocate (z_wt_local(nac,nmru))
      allocate (suz(nac,nmru))
      allocate (srz(nac,nmru))
!      allocate (srzwet(nac,nmru))
      allocate (qb(nmru))
      allocate (gw_loss(nmru))
      allocate (qdf(nmru))
      allocate (qpref(nmru))
      allocate (acm(nmru))
      allocate (afx(nmru))
      allocate (qof(nmru))
      allocate (qofs(nmru))
      allocate (quz(nmru))
      allocate (rex(nmru))
      allocate (quz_local(nac,nmru))
      allocate (qdf_local(nac,nmru))
      allocate (uz_infil(nac,nmru))
      allocate (uz2sat(nac,nmru))
      allocate (qvpref(nmru))
c      allocate (irrig_hyd_seg(nhydro))
      allocate (gw_in1(nmru))
      allocate (gw_in2(nmru))
      allocate (chan_loss(nhydro))

      webrdecl = 0

      return
      end

c***********************************************************************
c
c     webrinit - Initialize webmod_res module - get parameter values
c                and establish initial volumes
c

      integer function webrinit()

      USE WEBMOD_RESMOD
      USE WEBMOD_TOPMOD, ONLY : riparian, s_theta_0, 
     $                          s_rock_depth, s_root_depth
      USE WEBMOD_IO, ONLY : print_vse, nf, vse_lun
      IMPLICIT NONE
      integer, external :: length      
      integer is, ia, ih, filelen
      real acf
      logical v_alloc
      integer i

      webrinit = 1

      resstep1=.true.
      v_alloc=.TRUE. ! limit allocations to the first time through the mru loop


c  Get the areas of pervious versus impervious from the basin_topg module.
c  Since the values are static, its OK that we retrieve the values here.
c
c For now we will use the hortonian overland flow to simulate overland
c flow from impervious area. We can get more explicit later. Note need to
c use the * in the MMS get var function once the lines are commented out
c so as not to confuse the MMS parser.
c
c$$$      if(get*var('basin', 'mru_perv', nmru, 'real', mru_perv)
c$$$     +   .ne.0) return
c$$$
c$$$      if(get*var('basin', 'mru_imperv', nmru, 'real', mru_imperv)
c$$$     +   .ne.0) return
c$$$c
c Gather other static variables and parameters
c
c$$$      if(get*param('topc', 'dtinit', 1 , 'real', dtinit)
c$$$     +   .ne.0) return
           
      if(getparam('intcp', 'covden_sum', nmru, 'real', covden_sum)
     +   .ne.0) return 

      if(getparam('intcp', 'covden_win', nmru, 'real', covden_win)
     +   .ne.0) return 

      if(getparam('snow', 'WEI', nmru, 'real', wei)
     +   .ne.0) return 

c$$$      if(get*param('topc', 't0', nmru, 'real', T0)
c$$$     +   .ne.0) return
c$$$
c$$$      if(get*var('topc', 'tl', nmru, 'integer', tl)
c$$$     +   .ne.0) return

c$$$      if(get*var('topc', 'srmax', nmru, 'real', srmax)
c$$$     +   .ne.0) return
c$$$
      if(getvar('topc', 'uz_depth', nac*nmru, 'real', uz_depth)
     +   .ne.0) return

      if(getvar('potet', 'transp_on', nmru, 'integer', transp_on)
     +   .ne.0) return

c$$$      if(get*param('topc', 'szm', nmru, 'real', SZM)
c$$$     +   .ne.0) return
c$$$
c$$$      if(get*param('topc', 'td', nmru, 'real', TD)
c$$$     +   .ne.0) return
c$$$      
!      if(getparam('topc', 's_theta_0', nmru, 'real', s_theta_0)
!     +   .ne.0) return

      if(getparam('topc', 's_porosity', nmru, 'real', s_porosity)
     +   .ne.0) return

      if(getparam('topc', 's_theta_wp', nmru, 'real', s_theta_wp)
     +   .ne.0) return

      if(getparam('topc', 's_theta_fc', nmru, 'real', s_theta_fc)
     +   .ne.0) return

!      if(getparam('topc', 's_rock_depth', nmru, 'real', s_rock_depth)
!     +   .ne.0) return

!      if(getparam('topc', 's_root_depth', nmru, 'real', s_root_depth)
!     +   .ne.0) return

      if(getparam('topc', 's_satpref_zmax', nmru, 'real',
     $     s_satpref_zmax) .ne.0) return

      if(getparam('topc', 's_satpref_zmin', nmru, 'real',
     $     s_satpref_zmin) .ne.0) return

      if(getparam('topc', 'qdffrac', nmru, 'real', qdffrac)
     +   .ne.0) return

      if(getparam('topc', 'sbar0', nmru, 'real', SBAR0)
     +   .ne.0) return

c$$$      if(get*param('topc', 'infex', 1, 'real', INFEX)
c$$$     +   .ne.0) return
c$$$

      if(getparam('basin', 'nacsc', nmru, 'integer', nacsc)
     +   .ne.0) return

      if(getparam('basin', 'ac', nac*nmru, 'real', AC)
     +   .ne.0) return

      if(getparam('basin', 'mru_area', nmru, 'real', mru_area)
     +   .ne.0) return

      if(getparam('basin', 'mru_area_frac', nmru, 'real',
     $     mru_area_frac) .ne.0) return

      if(getparam('basin', 'basin_area', 1 , 'real', basin_area)
     +   .ne.0) return

      if(getvar('top2c', 'chan_area', nchan, 'real', 
     +  chan_area).ne.0) return

      if(getparam('top2c', 'mru2chan', nmru , 'integer', mru2chan)
     +   .ne.0) return

c      if(getparam('basin', 'qobsta', 1, 'integer', qobsta)
c     +   .ne.0) return

      if(getvar('routec', 'clark_segs', 1, 'integer', clark_segs)
     +   .ne.0) return

      if(getvar('routec', 'q', nhydro*nchan, 'real', q)
     +   .ne.0) return

      if(getvar('routec', 'ar_fill', nhydro*nchan, 'real',
     $     ar_fill) .ne.0) return

      if(getvar('routec', 'maxchanloss', nhydro, 'real',
     $     maxchanloss) .ne.0) return

      if(getparam('webr', 'c_can_depth', 1, 'real',
     +   c_can_depth).ne.0) return

      if(getparam('webr', 's_ohoriz_depth', nmru, 'real',
     +   s_ohoriz_depth).ne.0) return

c$$$      DT = dtinit
c
c Conversion factor for cubic meters to centimeters
c
      m3cm = 1e-4/basin_area
     
c
c Initialize storage variables and begin loop on MRU's
c
      basin_sto_cm = 0.0
      basin_intcp_sto_cm = 0.0
      basin_pweqv_cm = 0.0
      basin_soil_moist_cm = 0.0
      basin_sssto_cm = 0.0
      basin_gw_sto_cm = 0.0
      basin_chan_sto_cm = 0.0

      vmix_basin(1) = 0.0

      do is = 1, nmru
      
         vmix_mru(is,1) = 0.0

c
c Initialize impermeable area volume to 0.1 mm (1e-4 m) for the impermeable area.
c This minimum volume can be used to store dry deposition later.
c
c$$$         vmix_imp(is,1) = 0.0001* mru_imperv(is) * a_million
c
c Initialize canopy volume. Assume that depth is c_can_depth (defaults to 
c 0.1 mm = 0.0001 m) when the model sees the canopy as dry (intcp_stor=0.0).
c Solutes will be stored in this residual moisture in the canopy.
c
c Use the initial value of transp_on from the potet module to assign the
c initial canopy density.
c

         web_transp_on(is) = transp_on(is)
c
c calculate a fixed reservoir volume, in cubic meters, that will be exchanged
c back and forth between the canopy and the O-horizon on days of leaves-on and
c leaves-off
c
         can_ohoriz_vol(is) = 
     $        c_can_depth * (covden_sum(is)-covden_win(is))
     $        * mru_area(is) * a_million

         if(transp_on(is).eq.0) then
            vmix_can(is,1) = c_can_depth * covden_win(is) * mru_area(is)
     $           * a_million
         else
            vmix_can(is,1) = c_can_depth * covden_sum(is) * mru_area(is)
     $           * a_million
         end if
         
         basin_intcp_sto_cm = basin_intcp_sto_cm + vmix_can(is,1)

         vmix_mru(is,1) = vmix_mru(is,1) + vmix_can(is,1)
         vmix_basin(1) = vmix_basin(1) + vmix_can(is,1)
c     
c Initialize snowpack volume (WEI in inches)
c     
         vmix_snow(is,1) = wei(is) * mru_area(is) * a_million * inch2m

         basin_pweqv_cm = basin_pweqv_cm + vmix_snow(is,1)

         vmix_mru(is,1) = vmix_mru(is,1) + vmix_snow(is,1) 
         vmix_basin(1) = vmix_basin(1) + vmix_snow(is,1) 

c
c Initialize O-horizon volume, s_ohoriz_depth(nmru). The default depth
c is 0.5 mm. We might make s_ohoriz_depth a basin parameter or even static
c if it proves insensitive. This may be modified once we add explicit
c impervious areas. The o-horizon will also hold the residual water and
c solutes dropped from the canopy on the day of leaves off (transp_on=0)
c and provide the residual water back to the canopy on the first day
c of leaves on (transp_on=1). The residual canopy depth is set with the
c parameter c_can_depth (in meters).
c
         vmix_ohoriz(is,1) = s_ohoriz_depth(is) * mru_area(is)
     $        * a_million
         
         if(transp_on(is).eq.0)
     $        vmix_ohoriz(is,1) = 
     $        vmix_ohoriz(is,1) + can_ohoriz_vol(is)

         vmix_mru(is,1) = vmix_mru(is,1) + vmix_ohoriz(is,1) 
         vmix_basin(1) = vmix_basin(1) + vmix_ohoriz(is,1) 

c The variable uz_diff from the topmod_chem module contains the
c change in water content of the unsaturated zone during a given
c time step. The difference depth
c of water in the unsaturated zone at the end of this time step
c so that value will be loaded into the 1st nresinp index for the
c vmix_uz(nac,nmru,nresinp) matrix in the run routine. Here we will
c assign the initial volumes assigned in topmod_chem using the initial
c basin saturation deficit, the initial root zone deficit and the
c distribution of water table predicted using T0, szm, and the resp_coef.
c
            vmix_uzgen(is,1) = 0.0
            vmix_uzrip(is,1) = 0.0
            vmix_uzup(is,1) = 0.0
            do 75 ia = 1, nacsc(is)
              if(ia.eq.nacsc(is)) then
                ACF=0.5*AC(IA,is)
              else
                ACF=0.5*(AC(IA,is)+AC(IA+1,is))
              endif
              vmix_uz(ia,is,1) = uz_depth(ia,is) * acf *
     $             mru_area(is) * a_million
              vmix_uzgen(is,1) = vmix_uzgen(is,1) + vmix_uz(ia,is,1)
              if(riparian(ia,is))  then
                vmix_uzrip(is,1)=vmix_uzrip(is,1) + vmix_uz(ia,is,1)
              else
                vmix_uzup(is,1)=vmix_uzup(is,1) + vmix_uz(ia,is,1)
              end if
              vmix_mru(is,1) = vmix_mru(is,1) + vmix_uz(ia,is,1)
              vmix_basin(1) = vmix_basin(1) + vmix_uz(ia,is,1) 
              basin_sssto_cm =
     $             basin_sssto_cm + vmix_uz(ia,is,1) 
              vmix_uz2can(ia,is) = 0.0
              vmix_uz2sat(ia,is) = 0.0
              vmix_uz2qdf(ia,is) = 0.0

 75         continue
c
c Water transferred from the saturated zone to the unsaturated zone.
c
            vmix_sat2uz(is) = 0.0
c
c Water pumped at the end of the day for application the following morning
c

            vmix_well(is) = 0.0
c
c
c Initialize shallow preferential flow volume with residual volume
c equal wilting point moisture times root depth * percentage of
c flow attributed to direct flow
c
c Note that there is no residual volume if either the rooting depth
c or the percent preferential flow is equal to 0.
c
         vmix_qdf(is,1) = s_theta_wp(is)*s_root_depth(is)*qdffrac(is)
     $        * mru_area(is) * a_million

         vmix_mru(is,1) = vmix_mru(is,1) + vmix_qdf(is,1)
         vmix_basin(1) = vmix_basin(1) + vmix_qdf(is,1) 

c
c Initialize saturated zone volume with the initial storage in the
c subcatchment
c
         sr0(is) = (s_theta_fc(is) - s_theta_0(is)) * s_root_depth(is)

         basin_soil_moist_cm = basin_soil_moist_cm +
     $        (s_theta_0(is) * s_root_depth(is))*
     $        mru_area(is) * a_million
c
c No unsaturated zone storage (suz) to start, so the unsaturated zone
c storage will equal the root zone storage (just duplicate at the end
c of the init section)
c

         vmix_sat(is,1) = ((s_rock_depth(is) * s_porosity(is))
     $        - SBAR0(is) - SR0(is)) *mru_area(is) * a_million

         basin_gw_sto_cm = basin_gw_sto_cm + vmix_sat(is,1)

         vmix_mru(is,1) = vmix_mru(is,1) + vmix_sat(is,1)
         vmix_basin(1) = vmix_basin(1) + vmix_sat(is,1)
c
c Initialize storage in the preferential flow through the saturated zone
c with saturated soil moisture times the thickness of the preferential
c flow zone.
c
         vmix_satpref(is,1) = s_porosity(is) * (s_satpref_zmax(is) -
     $        s_satpref_zmin(is)) * mru_area(is) * a_million

c do not include this volume in the basin volume since it is accounted
c for in vmix_sat
c         vmix_basin(1) = vmix_basin(1) + vmix_satpref(is,1) 

c
c Initialize storage in the transient store that will mix all hillslope
c inputs. This transient reservoir will feed the stream reservoirs and
c might be eliminated in the future if we can feed the drainage directly.
c For now, we'll keep it separate to provide flexibility for coupling to
c a variety of routing routines.
c
         vmix_hill(is,1) = 0.0

c Initialize stream inputs to zero
         do 85 ih = 1, clark_segs
            vmix_hillexp(ih,is) = 0.0
 85      continue

c Assign initial MRU volume

         vmix_mru0(is) = vmix_mru(is,1)

c end of mru loop

      end do
c
c Initialize stream segment volumes. The initial and subsequent
c total volumes for the stream segments will be stored in the nmru+1 index
c Add volumes for channel leakage to initial volumes since it will be removed
c on the first time step
c

      do is = 1, clark_segs
         vmix_stream(is) = 0.0
         do ih = 1, nchan
            vmix_stream(is) = vmix_stream(is) + 
     $           q(is,ih)*chan_area(ih)*a_million
         end do
         vmix_stream(is) = vmix_stream(is) +
     $        maxchanloss(is)*basin_area*a_million
         vmix_diversion(is) = 0.0
         vmix_chan_loss(is) = 0.0
         basin_chan_sto_cm = basin_chan_sto_cm + vmix_stream(is) 
         vmix_basin(1) = vmix_basin(1) + vmix_stream(is) 
      end do

      vmix_basin0 = vmix_basin(1)

      print*,'vmix_basin0 = ',vmix_basin0
c
c Convert volumes to initial depths
c
      basin_sto_cm = vmix_basin0 * m3cm
      basin_intcp_sto_cm = basin_intcp_sto_cm * m3cm
      basin_pweqv_cm = basin_pweqv_cm * m3cm
      basin_soil_moist_cm = basin_soil_moist_cm * m3cm
      basin_sssto_cm = basin_sssto_cm * m3cm
      basin_gw_sto_cm = basin_gw_sto_cm * m3cm
      basin_chan_sto_cm = basin_chan_sto_cm * m3cm
c
c Open volume files if print_vse= 1 (basin and mru) or 2 (all reservoirs)
c
      if(print_vse.ge.1) then
!        nf=0
c composite basin volumes
        IF(control_string(out_dir,'output_dir').NE.0) RETURN
        path_len = index(out_dir,CHAR(0))-1   ! CHAR(0) is end of strings returned from control_string call
        vf_bas%file = out_dir(1:path_len)//'v_basin'
        inquire(file=vf_bas%file,exist=filflg)
        if (filflg) then
          tmplun = my_newunit()
          open(unit=tmplun,file=vf_bas%file,status='old')
          close(unit=tmplun,status='delete')
        endif
!----open the file.
        vf_bas%lun = my_newunit()
        open (unit=vf_bas%lun,file=vf_bas%file,access='sequential',
     $    form='formatted', status='new')
        nf=nf+1
        vse_lun(nf)=vf_bas%lun
        write(vf_bas%lun,10)
!
!open mru volume files (only compsite mru if print_vse = 1)
!
        allocate(vf_mru(nmru))
        do i = 1, nmru
! composite mru
          write(filename,20)i
          filelen=length(filename)
          vf_mru(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_mru(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_mru(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_mru(i)%lun = my_newunit()
          open (unit=vf_mru(i)%lun,file=vf_mru(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_mru(i)%lun
          write(vf_mru(i)%lun,10)
! additional reservoir files
          if(print_vse.eq.2) then
              if(v_alloc) then ! only allocate once
                allocate(vf_uzgen(nmru))
                allocate(vf_uzrip(nmru))
                allocate(vf_uzup(nmru))
                allocate(vf_can(nmru))
                allocate(vf_snow(nmru))
!        allocate(vf_imperv(nmru))
                allocate(vf_transp(nmru))
                allocate(vf_ohoriz(nmru))
                allocate(vf_uz(nmru,nac))
                allocate(vf_qdf(nmru))
                allocate(vf_sat(nmru))
                allocate(vf_satpref(nmru))
                allocate(vf_hill(nmru))
                allocate(vf_uz2sat(nmru))
                v_alloc=.FALSE.
              end if
! composite uz
          write(filename,30)i
          filelen=length(filename)
          vf_uzgen(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_uzgen(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_uzgen(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_uzgen(i)%lun = my_newunit()
          open (unit=vf_uzgen(i)%lun,file=vf_uzgen(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_uzgen(i)%lun
          write(vf_uzgen(i)%lun,10)
! composite riparian uz
          write(filename,40)i
          filelen=length(filename)
          vf_uzrip(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_uzrip(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_uzrip(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_uzrip(i)%lun = my_newunit()
          open (unit=vf_uzrip(i)%lun,file=vf_uzrip(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_uzrip(i)%lun
          write(vf_uzrip(i)%lun,10)
! composite upland uz
          write(filename,50)i
          filelen=length(filename)
          vf_uzup(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_uzup(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_uzup(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_uzup(i)%lun = my_newunit()
          open (unit=vf_uzup(i)%lun,file=vf_uzup(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_uzup(i)%lun
          write(vf_uzup(i)%lun,10)
! canopy
          write(filename,60)i
          filelen=length(filename)
          vf_can(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_can(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_can(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_can(i)%lun = my_newunit()
          open (unit=vf_can(i)%lun,file=vf_can(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_can(i)%lun
          write(vf_can(i)%lun,10)
! snowpack
          write(filename,70)i
          filelen=length(filename)
          vf_snow(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_snow(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_snow(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_snow(i)%lun = my_newunit()
          open (unit=vf_snow(i)%lun,file=vf_snow(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_snow(i)%lun
          write(vf_snow(i)%lun,10)
! impervious surface
    !      write(filename,80)i
    !      filelen=length(filename)
    !      vf_imperv(i)%file = out_dir(1:path_len)//filename(1:filelen)
    !        inquire(file=vf_imperv(i)%file,exist=filflg)
    !      if (filflg) then
    !        tmplun = my_newunit()
    !        open(unit=tmplun,file=vf_imperv(i)%file,status='old')
    !        close(unit=tmplun,status='delete')
    !      endif
    !!----open the file.
    !      vf_imperv(i)%lun = my_newunit()
    !      open (unit=vf_imperv(i)%lun,file=vf_imperv(i)%file,
    ! $      access='sequential',form='formatted', status='new')
    !      nf=nf+1
    !      vse_lun(nf)=vf_imperv(i)%lun
    !      write(vf_imperv(i)%lun,10)
! transpiration
          write(filename,90)i
          filelen=length(filename)
          vf_transp(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_transp(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_transp(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_transp(i)%lun = my_newunit()
          open (unit=vf_transp(i)%lun,file=vf_transp(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_transp(i)%lun
          write(vf_transp(i)%lun,12)('UZ',j,j=1,nac)
! o-horizon
          write(filename,100)i
          filelen=length(filename)
          vf_ohoriz(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_ohoriz(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_ohoriz(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_ohoriz(i)%lun = my_newunit()
          open (unit=vf_ohoriz(i)%lun,file=vf_ohoriz(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_ohoriz(i)%lun
          write(vf_ohoriz(i)%lun,10)
! individual unsaturated zone reservoirs
          do j = 1, nac
           write(filename,110)i,j
           filelen=length(filename)
           vf_uz(i,j)%file = out_dir(1:path_len)//filename(1:filelen)
             inquire(file=vf_uz(i,j)%file,exist=filflg)
           if (filflg) then
             tmplun = my_newunit()
             open(unit=tmplun,file=vf_uz(i,j)%file,status='old')
             close(unit=tmplun,status='delete')
           endif
    !----open the file.
           vf_uz(i,j)%lun = my_newunit()
           open (unit=vf_uz(i,j)%lun,file=vf_uz(i,j)%file,
     $      access='sequential',form='formatted', status='new')
           nf=nf+1
           vse_lun(nf)=vf_uz(i,j)%lun
           write(vf_uz(i,j)%lun,10)
          end do
! Direct flow
          write(filename,120)i
          filelen=length(filename)
          vf_qdf(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_qdf(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_qdf(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_qdf(i)%lun = my_newunit()
          open (unit=vf_qdf(i)%lun,file=vf_qdf(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_qdf(i)%lun
          write(vf_qdf(i)%lun,10)
! saturated zone
          write(filename,130)i
          filelen=length(filename)
          vf_sat(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_sat(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_sat(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_sat(i)%lun = my_newunit()
          open (unit=vf_sat(i)%lun,file=vf_sat(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_sat(i)%lun
          write(vf_sat(i)%lun,10)
! preferential flow through the saturated zone (tile drains)
          write(filename,140)i
          filelen=length(filename)
          vf_satpref(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_satpref(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_satpref(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_satpref(i)%lun = my_newunit()
          open (unit=vf_satpref(i)%lun,file=vf_satpref(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_satpref(i)%lun
          write(vf_satpref(i)%lun,10)
! Combined hillslope discharge (overland flow, qdf, and baseflow)
          write(filename,150)i
          filelen=length(filename)
          vf_hill(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_hill(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_hill(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_hill(i)%lun = my_newunit()
          open (unit=vf_hill(i)%lun,file=vf_hill(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_hill(i)%lun
          write(vf_hill(i)%lun,10)
! uz2sat - net flux of water from uz to sat from changes in water table
!          and sat water moving up to meet evap demand.
          write(filename,160)i
          filelen=length(filename)
          vf_uz2sat(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_uz2sat(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_uz2sat(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_uz2sat(i)%lun = my_newunit()
          open (unit=vf_uz2sat(i)%lun,file=vf_uz2sat(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_uz2sat(i)%lun
          write(vf_uz2sat(i)%lun,13)('UZ',j,j=1,nac)
          endif !print_vse=2, mru section
      enddo ! mru loop
! volumes of water exported to stream segments from each MRU on that day
      if(print_vse.ge.1) then
          allocate(vf_hillexp(nhydro))
          do i = 1, nhydro
          write(filename,170)i
          filelen=length(filename)
          vf_hillexp(i)%file = out_dir(1:path_len)//filename(1:filelen)
            inquire(file=vf_hillexp(i)%file,exist=filflg)
          if (filflg) then
            tmplun = my_newunit()
            open(unit=tmplun,file=vf_hillexp(i)%file,status='old')
            close(unit=tmplun,status='delete')
          endif
    !----open the file.
          vf_hillexp(i)%lun = my_newunit()
          open (unit=vf_hillexp(i)%lun,file=vf_hillexp(i)%file,
     $      access='sequential',form='formatted', status='new')
          nf=nf+1
          vse_lun(nf)=vf_hillexp(i)%lun
          write(vf_hillexp(i)%lun,15)i
          end do
! volumes of water in each stream segment at end of day
        filename = 'v_hyd'
        filelen=length(filename)
        vf_hyd%file = out_dir(1:path_len)//filename(1:filelen)
          inquire(file=vf_hyd%file,exist=filflg)
        if (filflg) then
          tmplun = my_newunit()
          open(unit=tmplun,file=vf_hyd%file,status='old')
          close(unit=tmplun,status='delete')
        endif
    !----open the file.
        vf_hyd%lun = my_newunit()
        open (unit=vf_hyd%lun,file=vf_hyd%file,
     $      access='sequential',form='formatted', status='new')
        nf=nf+1
        vse_lun(nf)=vf_hyd%lun
        write(vf_hyd%lun,17)('hyd',i,i=1,nhydro)
        endif ! print_vse=2, hydro section
      endif ! print_vse=1 
      
      webrinit = 0
!     '12345678911234567892123456783123456789412
!
! headers in files
!
 10   format('Volumes, in cubic meters, Specific inputs, ',
     $ '(except for ET) listed after "Final" column'/'nstep Year Mo Dy',
     $ '          Init         Inputs        Outputs',
     $ '          Final         Precip             ET',
     $ '     Impervious    Throughfall           Melt',
     $ '      O-horizon             UZ        UZ_Pref',
     $ '       Sat_zone       Sat_Pref',
     $ '      Hillslope           Well      Diversion',
     $ ' Outside_Source           GW_1           GW_2')

 12   format('Volumes in cubic meters pulled from each UZ by ',
     $ 'transpiration.'/'nstep Year Mo Dy',50(A10,I4.3))
 13   format('Volumes, in cubic meters, of recharge from each UZ',
     $ /'nstep Year Mo Dy',50(A10,I4.3))
 15   format('Volumes, in cubic meters'/'nstep Year Mo Dy Inputs',
     $       ' to stream segment ',I5,' beginning at MRU 1')
 17   format('Volumes, in cubic meters at end of day'/'nstep Year Mo Dy',
     $       '     Discharge',20(A10,I4.3))
!
! file names
!
 20   format('v_mru',I3.3)
 30   format('v_mru',I3.3,'_uzgen')
 40   format('v_mru',I3.3,'_uzrip')
 50   format('v_mru',I3.3,'_uzup')
 60   format('v_mru',I3.3,'_can')
 70   format('v_mru',I3.3,'_snow')
! 80   format('v_mru',I3.3,'_imperv')
 90   format('v_mru',I3.3,'_transp')
 100  format('v_mru',I3.3,'_ohoriz')
 110  format('v_mru',I3.3,'_uz',I2.2)
 120  format('v_mru',I3.3,'_qdf')
 130  format('v_mru',I3.3,'_sat')
 140  format('v_mru',I3.3,'_satpref')
 150  format('v_mru',I3.3,'_hill')
 160  format('v_mru',I3.3,'_uz2sat')
 170  format('v_hyd',I3.3)

      return
      end

c***********************************************************************
c
c     webrrun - convert model fluxes to storage changes and mixing
c               coefficients for input into phreeq_mms
c

      integer function webrrun()

      USE WEBMOD_RESMOD
      USE WEBMOD_IO, only: phreeqout, print_vse, debug
      USE WEBMOD_OBSHYD, ONLY: runoff,nirrig_int,nirrig_ext
      USE WEBMOD_IRRIG, ONLY: irrig_ext_mru,irrig_sat_mru,irrig_hyd_mru,
     $                        irrig_hyd_seg
      USE WEBMOD_TOPMOD, only: z_wt_local, srzwet, riparian,
     $                         s_theta_0, s_root_depth, qobsta
      IMPLICIT NONE

c variables and parameters from other modules
      double precision timestep
      real fac
      integer endper, datetime(6), nstep
      logical end_run, end_yr, end_mo, end_dy, end_storm
      logical basinq_found
      integer i, is, ia, ih, k
      real smav_basin, p
      real covden, transp_rate
      real acf
      real snowfrac, offrac
c
c Temporary volume variables for using within the mru loops.
c No storage variables used since these will be assigned
c using the vmix_??(is,4) indices
c
      double precision v_dep,v_ppt,v_ppt2can,v_ppt2gnd
      double precision v_irrsat,v_irrsat2can,v_irrsat2gnd
      double precision v_irrhyd,v_irrhyd2can,v_irrhyd2gnd
      double precision v_irrext,v_irrext2can,v_irrext2gnd
      double precision v_intcp, v_intcp_evap
      double precision v_transp, v_thruf
      double precision v_netdep, v_netrain, v_netsnow
      double precision v_pweqv, v_snowmelt, v_snowevap
      double precision v_surfdep, v_surfdep2, dep_ratio, v_sssto
      double precision v_ofhort, v_ofdunn, v_sroff  ! sroff = hort+dunn
      double precision v_infil_loc, v_infil
      double precision v_qdf, v_uz_et
      double precision v_gw1_in,v_gw2_in
      double precision v_qwet_loc, v_qwet
      double precision v_rech_loc, v_rech, v_qdf_loc
      double precision v_uz2sat,v_sat2uz
      double precision v_qwell, v_qpref,v_gw_loss
      double precision v_qb
      double precision v_chan_loss
      double precision ppt2hill,irrext2hill,irrsat2hill,irrhyd2hill
      double precision can2hill,melt2hill
      double precision uz_bypass,bypass_frac
      double precision vsat_tmp

      webrrun = 1

c
c Get date/time stamp for debugging
c
      call dattim('now', datetime)

      nstep=getstep()

      timestep = deltim()
c
c from io_chem
c
      if(getvar('io', 'endper', 1, 'integer', endper)
     +   .ne.0) return
c
c from basin_topg - no variables needed
c
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
c Get variables from other modules:
c
c Precipitation fluxes from precip_prms
c
      if(getvar('precip', 'mru_dep', nmru, 'real', mru_dep)
     +   .ne.0) return

      if(getvar('precip', 'mru_ppt', nmru, 'real', mru_ppt)
     +   .ne.0) return

c      if(getvar('irrig', 'irrig_sat_mru', nmru, 'real',
c     $     irrig_sat_mru).ne.0) return

c      if(getvar('irrig', 'irrig_hyd_mru', nmru, 'real',
c     $     irrig_hyd_mru).ne.0) return

c      if(getvar('irrig', 'irrig_ext_mru', nmru, 'real',
c     $     irrig_ext_mru).ne.0) return

c
c Fluxes from impermeable areas
c ;later

c Transpiration flag from potet_hamon_prms
c
      if(getvar('potet', 'transp_on', nmru, 'integer', transp_on)
     +   .ne.0) return
c
c Canopy fluxes and from intcp_prms
c
      if(getvar('intcp', 'net_dep', nmru, 'real', net_dep)
     +   .ne.0) return

      if(getvar('intcp', 'net_rain', nmru, 'real', net_rain)
     +   .ne.0) return

      if(getvar('intcp', 'net_snow', nmru, 'real', net_snow)
     +   .ne.0) return

      if(getvar('intcp', 'intcp_on', nmru, 'integer', intcp_on)
     +   .ne.0) return

      if(getvar('intcp', 'intcp_stor', nmru, 'real', intcp_stor)
     +   .ne.0) return

      if(getvar('intcp', 'intcp_evap', nmru, 'real', intcp_evap)
     +   .ne.0) return

c
c Snowpack fluxes from nwsmelt_topg
c
      if(getvar('nwsmelt', 'pkwater_equiv', nmru,
     $     'real', pkwater_equiv) .ne.0) return

      if(getvar('nwsmelt', 'snowmelt', nmru, 'real', snowmelt)
     +   .ne.0) return

      if(getvar('nwsmelt', 'psoilmru', nmru, 'real', psoilmru)
     +   .ne.0) return

      if(getvar('nwsmelt', 'snow_evap', nmru, 'real', snow_evap)
     +   .ne.0) return

c
c Hillslope fluxes from topmod_chem
c
      if(getvar('topc', 'srz_sc', nmru, 'real',srz_sc).ne.0) return

      if(getvar('topc', 'suz_sc', nmru, 'real',suz_sc).ne.0) return

      if(getvar('topc', 'smav_basin', 1, 'real', 
     $     smav_basin).ne.0) return

      if(getvar('topc', 'sae', nmru, 'real',
     $     sae).ne.0) return

      if(getvar('topc', 'sae_local', nac*nmru, 'real',
     $     sae_local).ne.0) return

      if(getvar('topc', 'sd', nac*nmru, 'real', SD).ne.0) return

!      if(getvar('topc', 'z_wt_local', nac*nmru, 'real', 
!     +     z_wt_local).ne.0) return

      if(getvar('topc', 'uz_depth', nac*nmru, 'real', 
     +     uz_depth).ne.0) return

      if(getvar('topc', 'suz', nac*nmru, 'real', 
     +     SUZ).ne.0) return

      if(getvar('topc', 'srz', nac*nmru, 'real', 
     + SRZ) .ne.0) return

      if(getvar('topc', 'qb', nmru, 'real', 
     +     qb).ne.0) return

!      if(get*var('topc', 'srzwet', nac*nmru, 'real', 
!     +     srzwet).ne.0) return

      if(getvar('topc', 'p', 1, 'real', p).ne.0) return

      if(getvar('topc', 'qdf', nmru, 'real', qdf).ne.0) return

      if(getvar('topc', 'qpref', nmru, 'real', qpref).ne.0) return

      if(getvar('topc', 'acm', nmru, 'real', acm).ne.0) return

      if(getvar('topc', 'afx', nmru, 'real', afx).ne.0) return

      if(getvar('topc', 'qof', nmru, 'real', QOF) .ne.0) return

      if(getvar('topc', 'qofs', nmru, 'real', QOFS).ne.0) return

      if(getvar('topc', 'quz_local', nac*nmru, 'real',
     $     quz_local).ne.0) return

      if(getvar('topc', 'qdf_local', nac*nmru, 'real',
     $     qdf_local).ne.0) return

      if(getvar('topc', 'uz_infil', nac*nmru, 'real',
     $     uz_infil).ne.0) return

      if(getvar('topc', 'uz2sat', nac*nmru, 'real',
     $     uz2sat).ne.0) return

      if(getvar('topc', 'quz', nmru, 'real', QUZ).ne.0) return

      if(getvar('topc', 'qvpref', nmru, 'real', qvpref).ne.0) return

      if(getvar('topc', 'rex', nmru, 'real', REX).ne.0) return

      if(getvar('topc', 'gw_in1', nmru, 'real', gw_in1).ne.0) return

      if(getvar('topc', 'gw_in2', nmru, 'real', gw_in2).ne.0) return

      if(getvar('topc', 'gw_loss', nmru, 'real', gw_loss).ne.0) return
c
c Stream routing info from top2clark and route_clark. The static
c variables, mru2chan, chan_area, clark_segs and ar_fill were read
c in the init section.

      if(getvar('routec', 'q', nhydro*nchan, 'real', 
     +  q).ne.0) return

      if(getvar('routec', 'chan_loss', nhydro, 'real', 
     +  chan_loss).ne.0) return

c      if(getvar('topc', 'irrig_hyd_seg', nhydro, 'double', 
c    +  irrig_hyd_seg).ne.0) return

c
c observed discharge. Convert to equivalent depth, basin_qobs_cm, below
c
c        if(getvar('obs', 'runoff', nobs, 'real', runoff_m3s)
c     +  .ne.0) return
c
c     Convert observed basin discharge from cfs to cubic meters per second
c
      basinq_found = .false.
      if(nobs.gt.0) then
        do 105 j=1,nobs
           runoff_m3s(j) = runoff(j) * 0.028316847
           if(j.eq.qobsta) then
              basin_qobs_m3s = runoff_m3s(j)
              basinq_found = .true.
           end if
 105    continue
      end if
      if (.not.basinq_found) then
         basin_qobs_m3s = 0.0
c         print*,'nobs equals zero or qobsta does not '//
c     $     'point to observed discharge data. basin_qobs_m3s '//
c     $     'and basin_qobs_cm set to zero.'
      end if

      fac = timestep/24

      basin_qobs_cm= basin_qobs_m3s *8.64*fac/basin_area

c
c If not the first time step, transfer final volumes from last step
c to initial volumes for this step
c
      if(.not.resstep1) then
         vmix_basin(1) = vmix_basin(4)
         do 2 is = 1,nmru
c     vmix_imp(is,1) = vmix_imp(is,1)
            vmix_mru(is,1) = vmix_mru(is,4)
            vmix_uzgen(is,1) = vmix_uzgen(is,4)
            vmix_uzrip(is,1) = vmix_uzrip(is,4)
            vmix_uzup(is,1) = vmix_uzup(is,4)
            vmix_can(is,1) = vmix_can(is,4) 
            vmix_snow(is,1) = vmix_snow(is,4)
            vmix_ohoriz(is,1) = vmix_ohoriz(is,4)
            vmix_qdf(is,1) = vmix_qdf(is,4) 
            do 3 j = 1,nacsc(is)
               vmix_uz(j,is,1) = vmix_uz(j,is,4)
 3          continue
            vmix_sat(is,1) = vmix_sat(is,4)
            vmix_satpref(is,1) = vmix_satpref(is,4)
            vmix_hill(is,1) = vmix_hill(is,4)
            do 4 k = 1, clark_segs
               vmix_hillexp(k,is) = 0.0
 4          continue
 2       continue
      else
         resstep1=.false.
c     
c     write debug header
!      write(chemout_file_unit,233)
      end if

c
c Zero basin volume accumulators
c
      do 5 is = 2, nresinp
         vmix_basin(is) = 0.0
 5    continue

      basin_sto_cm = 0.0
      basin_in_cm = 0.0
      basin_out_cm = 0.0
      basin_et_cm = 0.0
      basin_ppt_cm = 0.0
      basin_irr_ext_cm = 0.0
      basin_irr_sat_cm = 0.0
      basin_irr_hyd_cm = 0.0
      basin_intcp_sto_cm = 0.0
      basin_intcp_cm = 0.0
      basin_intcp_evap_cm = 0.0
      basin_transp_cm = 0.0
      basin_thruf_cm = 0.0
      basin_net_dep_cm = 0.0
      basin_net_rain_cm = 0.0
      basin_net_snow_cm = 0.0
      basin_pweqv_cm = 0.0
      basin_snowmelt_cm = 0.0
      basin_snowevap_cm = 0.0
      basin_surfdep_cm = 0.0
      basin_soil_moist_cm = 0.0
      basin_sssto_cm = 0.0
      basin_ofhort_cm = 0.0
      basin_ofdunn_cm = 0.0
      basin_sroff_cm = 0.0
      basin_infil_cm = 0.0
      basin_qdf_cm = 0.0
      basin_uz_et_cm = 0.0
      basin_gw_sto_cm = 0.0
      basin_gw1_in_cm = 0.0
      basin_gw2_in_cm = 0.0
      basin_qwet_cm = 0.0
      basin_vpref_cm = 0.0
      basin_uz2sat_cm = 0.0
      basin_recharge_cm = 0.0
      basin_sat2uz_cm = 0.0
      basin_qwell_cm = 0.0
      basin_qpref_cm = 0.0
      basin_gwloss_cm = 0.0
      basin_gwflow_cm = 0.0
      basin_chan_sto_cm = 0.0
      basin_stflow_cm = 0.0
      basin_chan_div_cm = 0.0
      basin_chan_loss_cm = 0.0
      basin_qsim_cm = 0.0
      basin_qsim_m3s = 0.0
c
c MRU Loop
c
      do 10 is = 1,nmru

c debug
c$$$            if (is.eq.3) then
!               write(26,295)(quz_local(ij,is),uz2sat(ij,is),
!     $              srzwet(ij,is),ij=1,11)
c$$$            end if
c$$$            if(datetime(3).eq.9.and.is.eq.3) return
c end debug
c
c Zero the volume accumulators. The first index (j=1) contains the storage
c from the previous time step. vmix_sat2uz trancks the transfer of water from
c the saturated zone to the unsaturated zone that results from root zone wetting
c and stranded pore water as the water table falls. Since there is only
c one source there is no need for the nresinp dimension.
c
         vmix_sat2uz(is) = 0.0
         do 100 j = 2,nresinp
c     vmix_imp(is,2) = 0.0
            vmix_mru(is,j) = 0.0
            vmix_can(is,j) = 0.0
            vmix_snow(is,j) = 0.0
            vmix_ohoriz(is,j) = 0.0
            vmix_qdf(is,j) = 0.0
            do 1000 k = 1,nacsc(is)
               vmix_uz(k,is,j) = 0.0
 1000       continue
            vmix_uzgen(is,j) = 0.0
            vmix_uzrip(is,j) = 0.0
            vmix_uzup(is,j) = 0.0
            vmix_sat(is,j) = 0.0
            vmix_satpref(is,j) = 0.0
            vmix_hill(is,j) = 0.0
 100     continue

c
c Zero the local volumes and average basin depths
c
      v_dep = 0.0
      v_ppt = 0.0
      v_ppt2can = 0.0
      v_ppt2gnd = 0.0
      v_irrsat = 0.0
      v_irrsat2can = 0.0
      v_irrsat2gnd = 0.0
      v_irrhyd = 0.0
      v_irrhyd2can = 0.0
      v_irrhyd2gnd = 0.0
      v_irrext = 0.0
      v_irrext2can = 0.0
      v_irrext2gnd = 0.0
      v_intcp = 0.0
      v_intcp_evap = 0.0
      v_transp = 0.0
      v_thruf = 0.0
      v_netdep = 0.0
      v_netrain = 0.0
      v_netsnow = 0.0
      v_pweqv = 0.0
      v_snowmelt = 0.0
      v_snowevap = 0.0
      v_surfdep = 0.0
      v_surfdep2 = 0.0
      v_sssto = 0.0
      v_ofhort = 0.0
      v_ofdunn = 0.0
      v_sroff = 0.0
      v_infil = 0.0
      v_qdf = 0.0
      v_uz_et = 0.0
      v_gw1_in = 0.0
      v_gw2_in = 0.0
      v_qwet = 0.0
      v_uz2sat = 0.0
      v_rech = 0.0
      v_sat2uz = 0.0
      v_qwell = 0.0
      v_qpref = 0.0
      v_gw_loss = 0.0
      v_qb = 0.0
      v_chan_loss = 0.0
      ppt2hill = 0.0
      irrext2hill = 0.0
      irrsat2hill = 0.0
      irrhyd2hill = 0.0
      can2hill = 0.0
      melt2hill = 0.0
      uz_bypass = 0.0
      bypass_frac = 0.0
c
c Leave out impervious areas for now.
c     vmix(is,2) = rain + snow
c     vmix_imp(is,5) = rain
c     vmix_imp(is,6) = ET
c     vmix_imp(is,9) = snow melt 
c
c Canopy Volumes:
c If there is interception on the canopy, or it has rained, Canopy ET
c (intcp_evap) will be removed from the storage and no transpiration
c nor mobilization of solutes from the root zone will occur on that day.
c
c If the canopy is dry as indicated by intcp_on=0, then the canopy will
c accumulate solutes drawn from the root zone at a rate of 75% of the
c total root zone evaporation (sae) during transpiration periods and
c 5% of the sae during non transpiration periods. The source of the
c water supplied to the canopy will be a well-mixed reactor of root
c zone contributions weighted by the local ET predicted for each wetness
c index. Transpiration will result in no change of canopy volume.

c Note: Canopy variables intcp_evap and incp_stor in inches
c
        if(transp_on(is).eq.1) then
           covden = covden_sum(is)
           transp_rate = 0.7
        else
           covden = covden_win(is)
           transp_rate = 0.1
        end if
c record minimum canopy volume for isotopic fractionation
        vmin_canopy(is) = c_can_depth*covden*mru_area(is)*a_million
c
c Compute volumes of precipitation and irrigation
c v_irrsat and v_irrhyd set to zero when nirrig_int=0
c v_irrext set to zero when nirrig_ext=0
c totals
        v_ppt = mru_ppt(is)*mru_area(is)*a_million*inch2m
        if(nirrig_int.ne.0) then
          v_irrsat = irrig_sat_mru(is)*mru_area(is)*a_million*inch2m
          v_irrhyd = irrig_hyd_mru(is)*mru_area(is)*a_million*inch2m
        else
          v_irrsat = 0
          v_irrhyd = 0
        endif
        if(nirrig_ext.ne.0) then
          v_irrext = irrig_ext_mru(is)*mru_area(is)*a_million*inch2m
        else
          v_irrext = 0
        endif
        v_dep = v_ppt + v_irrsat + v_irrhyd + v_irrext
c
c Track basin and mru inputs in their respective categories
c precip (5)
        vmix_mru(is,5) = v_ppt 
        vmix_basin(5) = vmix_basin(5)+ v_ppt
c irrsat (17)
        vmix_mru(is,17) = v_irrsat 
        vmix_basin(17) = vmix_basin(17)+ v_irrsat
c irrhyd (18)
        vmix_mru(is,18) = v_irrhyd 
        vmix_basin(18) = vmix_basin(18)+ v_irrhyd
c irrext (19)
        vmix_mru(is,19) = v_irrext 
        vmix_basin(19) = vmix_basin(19)+ v_irrext

        vmix_mru(is,2) = v_dep
        vmix_basin(2) = vmix_basin(2)+ v_dep
c
c to canopy
        v_ppt2can = v_ppt*covden
        v_irrsat2can = v_irrsat*covden
        v_irrhyd2can = v_irrhyd*covden
        v_irrext2can = v_irrext*covden

c unimpeded by canopy
        v_ppt2gnd = v_ppt*(1-covden)
        v_irrsat2gnd = v_irrsat*(1-covden)
        v_irrhyd2gnd = v_irrhyd*(1-covden)
        v_irrext2gnd = v_irrext*(1-covden)

c On the day that transpiration stops (leaves off), the concentration of
c the remaining leaves will remain constant. The mass of water and solutes
c on the leaves that fell will be delivered to the o-horizon
c
c Also, on days of leaves off, make the transpiration rate zero so that
c throughfall chemistry reflects the solutions at the start
c of the day and does not need to account for mixing in the unsaturated
c zone first.


        if(transp_on(is).eq.0.and.web_transp_on(is).eq.1) then !leaves off
           web_transp_on(is) = 0
           transp_rate = 0.0
           vmix_can(is,3) = vmix_can(is,1) -                ! move canopy water
     $          (covden_win(is)/covden_sum(is))* vmix_can(is,1) ! to throughfall
           vmix_ohoriz(is,8) = can_ohoriz_vol(is)  ! recharge o-horizon

c On the first day of transpiration, the area of new summer leaves will
c gain water and solutes from the O-horizon with a volume equal to the
c residual canopy storage times the increase from covden_win to covden_sum.


        else if(transp_on(is).eq.1.and.web_transp_on(is).eq.0) then !leaves on
           web_transp_on(is) = 1
           vmix_can(is,10) = can_ohoriz_vol(is)
           vmix_ohoriz(is,3) = can_ohoriz_vol(is)  ! collect water from o-horizon
        end if

c
c compute thrufall and intercepts ahead of inputs since transpiration will only
c occur on days with no thrufall 
c
c     Discharge volume (throughfall) (add additional throughfall on 
c     day of leaves off)

        vmix_can(is,3) = (net_dep(is)-(mru_dep(is)*(1-covden)))
     $       * inch2m * mru_area(is) * a_million + vmix_can(is,3)
c separate out true thufall from water flux on day of leaves off
        v_thruf = vmix_can(is,3) - vmix_ohoriz(is,8)
c zero out rounding errors. Throughfall less than 1 liter set to zero.
        if(abs(v_thruf).lt.1e-3) v_thruf=0.0
        if(abs(vmix_can(is,3)).lt.1e-3) vmix_can(is,3)=0.0
c
c Further development may consider distinguishing bare ground, grass, shrubs,
c or trees further than just the volume and density fluxes reported by
c the interception model.

c$$$        if(intcp_evap(is).gt.0.or.mru_dep(is).gt.0.0) then ! Canopy is wet
        vmix_can(is,5)  = v_ppt2can !precip
        vmix_can(is,17) = v_irrsat2can !well irrigation
        vmix_can(is,18) = v_irrhyd2can !river diversion
        vmix_can(is,19) = v_irrext2can !external irrigation

c          vmix_can(is,11) = 0        ! no root zone water if canopy is wet (already zero from init)

        vmix_can(is,6) = intcp_evap(is) * inch2m * mru_area(is) ! ET
     $       * a_million

c Allow transpiration on days when there is no throughfall
        if(v_thruf.eq.0.0) then
           vmix_can(is,11) = sae(is) * transp_rate !sae in meters
     $          *covden * mru_area(is) * a_million
c Add transpired water to canopy evap
           vmix_can(is,6) = vmix_can(is,6) + vmix_can(is,11) 
        end if
c
c Evaporate sae volume not transpired to canopy from the unsaturated zone
c
        v_uz_et = (sae(is)*mru_area(is)* a_million) - vmix_can(is,11)
c
c Sum of inputs, precip + uz (including recharge from 
c O-horizon for leaves-on day)
        vmix_can(is,2) = vmix_can(is,5)+ vmix_can(is,17)+
     $       vmix_can(is,18)+vmix_can(is,19)+
     $       vmix_can(is,10)+vmix_can(is,11)
c
c Set export to zero if negative result from rounding errors
c
        if(vmix_can(is,3).lt.0.0) vmix_can(is,3) = 0.0

c volume of interception, transpiration and throughfall - 
c corrected for transpiration and exchange with o-horizon
c on days of leaves off and leaves on

c$$$        v_thruf = vmix_can(is,3) - vmix_ohoriz(is,8)
        v_transp =  vmix_can(is,11)
        v_intcp = vmix_can(is,2)- v_transp - vmix_can(is,10) - v_thruf

c Final volume

        vmix_can(is,4) = vmix_can(is,1)+vmix_can(is,2)-vmix_can(is,3)
     $       - vmix_can(is,6)
c
c Collect evaporation in mru and basin variables
c
        vmix_mru(is,6) = vmix_mru(is,6) + vmix_can(is,6)
        vmix_basin(6) = vmix_basin(6) + vmix_can(is,6)

c volumes of net deposition, rain, and snow

        v_netrain = net_rain(is)*mru_area(is)*a_million*inch2m
        v_netsnow = net_snow(is)*mru_area(is)*a_million*inch2m
        v_netdep  = net_dep(is)*mru_area(is)*a_million*inch2m ! should equal sum

c
c - check final volume against intcp_stor
c
c$$$        can_depth_check = (vmix_can(is,4) - c_can_depth*a_million
c$$$     $     *mru_area(is)*covden)/a_million/mru_area(is)/covden/inch2m
c$$$
c$$$        if (vmix_can(is,4).gt.1000.or.(is.eq.1.and.
c$$$     $       vmix_can(is,3).lt.-0.1)) then
c$$$            print*,'depth check: ',can_depth_check, 
c$$$     $       '  intcp_stor: ',intcp_stor(is)
c$$$            print*,vmix_can(is,1),
c$$$     $       vmix_can(is,2),vmix_can(is,3),  vmix_can(is,4),
c$$$     $           vmix_can(is,5),vmix_can(is,11),  vmix_can(is,6)
c$$$        end if

c$$$        vmix_can(is,1) = vmix_can(is,4)

c$$$c
c$$$c test vmix_ohoriz
c$$$c
c$$$        vmix_ohoriz(is,4) = vmix_ohoriz(is,1) + vmix_ohoriz(is,8)
c$$$     $       - vmix_ohoriz(is,3)
c$$$        vmix_ohoriz(is,1) = vmix_ohoriz(is,4)
c
c The chemistry of the net_dep that contributes to snowpack or falls
c on bare ground (psoilmru) is a mixture of daily rainfall in mru_ppt
c and the throughfall, vmix_can(is,3).
c
c If there was no snowpack and no inputs, then 
c assign zeros to final volume and snowfrac and skip
c the rest of the snow section.

        if (vmix_snow(is,1).eq.0.0.and.net_dep(is).eq.0.0) then
           vmix_snow(is,4) = 0.0
           snowfrac = 0.0
        else

c The fraction of the net deposition that contributes to the
c snowpack, 1-(psoilmru/net_dep), will be placed in the variable snowfrac.
c 
           if(net_dep(is).gt.0.0)  then
              snowfrac = 1 - psoilmru(is)/net_dep(is)
           else
              snowfrac = 0
           end if
c 
c Calculate snowpack inputs
c ([5-precip, 17-19 irrigation deposition unimpeded by canopy],
c  6-ET, and 8-throughfall)
c Precip, irrigation and snowmelt variables are in inches
c 
c Snowpack fluxes are described by
c Change in pack storage = (net_dep-psoilmru) - snow_evap - snowmelt
c
c Snowpack
           vmix_snow(is,5)=v_ppt2gnd*snowfrac
           vmix_snow(is,17)=v_irrsat2gnd*snowfrac
           vmix_snow(is,18)=v_irrhyd2gnd*snowfrac
           vmix_snow(is,19)=v_irrext2gnd*snowfrac

           vmix_snow(is,6)=snow_evap(is)*mru_area(is)*a_million*inch2m
c 
c Canopy throughfall corrected for artifact of flux to 0_horizon on day
c of leaves off
c 
           vmix_snow(is,8) = v_thruf*snowfrac
c 
c Sum snow inputs
c 
           vmix_snow(is,2) = vmix_snow(is,5) + vmix_snow(is,17)+
     $          vmix_snow(is,18)+ vmix_snow(is,19) + vmix_snow(is,8)
c 
c Snowmelt output
c 
           vmix_snow(is,3) = snowmelt(is)*mru_area(is)*a_million*inch2m
c 
c Total snowpack fluxes (set to packwater equivalence to avoid rounding
c errors.
c 
            vmix_snow(is,4) = pkwater_equiv(is)*inch2m*mru_area(is)*
     $           a_million

c$$$           vmix_snow(is,4) = vmix_snow(is,1) + vmix_snow(is,2)
c$$$     $          - vmix_snow(is,3)-vmix_snow(is,6)

c 
c Debug section
c 
!           if(is.eq.1) then
!              write(*,396)nstep, vmix_snow(is,1),vmix_snow(is,2),
!     $             vmix_snow(is,3),vmix_snow(is,4),vmix_snow(is,5),
!     $             vmix_snow(is,6), vmix_snow(is,8)
!           end if
c
c /debug
c
 396       format(I8,7(1x,f12.3))
 397       format(12(1x,f10.3))

c 
c Collect evaporation in basin variable
c     
           vmix_mru(is,6) = vmix_mru(is,6) + vmix_snow(is,6)
           vmix_basin(6) = vmix_basin(6) + vmix_snow(is,6)

        end if

c
c Compute hillslope inputs before overland flow and before vertical bypass
c volumes are considered.
c
        ppt2hill = v_ppt2gnd * (1-snowfrac) ! atmospheric ppt
        irrext2hill = v_irrext2gnd * (1-snowfrac) ! external irrigation
        irrsat2hill = v_irrsat2gnd * (1-snowfrac) !atmospheric irrsat
        irrhyd2hill = v_irrhyd2gnd * (1-snowfrac) ! atmospheric irrhyd
        can2hill = v_thruf*(1-snowfrac) ! thrufall
        melt2hill = (snowmelt(is)*mru_area(is)*a_million*inch2m)  ! snowmelt
c
c Total surface depostion to be distributed as overland flow, infiltration,
c and vertical bypass.
c
        v_surfdep = ppt2hill + irrext2hill + irrsat2hill +
     $       irrhyd2hill + can2hill + melt2hill
c
c Calculate overland flow inputs as a fraction, offrac, of total inputs
c to the land surface(5-ppt, 17-19 irrigation, 8-throughfall, 9-melt).
c Throughfall inputs are increased by can_ohoriz_vol  on the day of 
c leaves off. The overland flow fraction does not include mass transferred
c between the canopy and o-horiz on days of leaves on or leaves off
c
        if(qof(is).gt.0) then
           offrac = qof(is)/(snowmelt(is)+psoilmru(is))/inch2m
           vmix_ohoriz(is,5)= ppt2hill*offrac
           vmix_ohoriz(is,19)= irrext2hill*offrac
           vmix_ohoriz(is,17)= irrsat2hill*offrac
           vmix_ohoriz(is,18)= irrhyd2hill*offrac
c (is,8)includes moisture gained on day of leaves off, Add to this,
c the fraction of throughfall when overland flow occurs.
           vmix_ohoriz(is,8) = vmix_ohoriz(is,8)+
     $          can2hill*offrac
           vmix_ohoriz(is,9) = melt2hill*offrac
c
c Sum O-horizon inputs. 
c
           vmix_ohoriz(is,2) = vmix_ohoriz(is,5)+
     $          vmix_ohoriz(is,19)+vmix_ohoriz(is,17)+
     $          vmix_ohoriz(is,18)+vmix_ohoriz(is,8)+
     $          vmix_ohoriz(is,9)
c
c Compute O-horizon outputs
c (same as inputs plus volume to canopy on day of leaves on and less
c  volume from canopy on day of leaves off)
c
           vmix_ohoriz(is,3) = vmix_ohoriz(is,3) + vmix_ohoriz(is,5)+
     $          vmix_ohoriz(is,19)+vmix_ohoriz(is,17)+
     $          vmix_ohoriz(is,18)+ vmix_ohoriz(is,9)+
     $          v_thruf*(1-snowfrac)*offrac
c     $       +vmix_can(is,10)
c
c Compute new O-horizon volume
c
        else
           offrac = 0.0
           vmix_ohoriz(is,2) = vmix_ohoriz(is,8)
        end if
        vmix_ohoriz(is,4) = vmix_ohoriz(is,1) + vmix_ohoriz(is,2)
     $       - vmix_ohoriz(is,3)
c
c accumulate overland flow volumes for basin summaries
c
        v_ofhort = rex(is)* mru_area(is) * a_million
        v_ofdunn = qofs(is)* mru_area(is) * a_million
        v_sroff = qof(is)* mru_area(is) * a_million

c Compute weighted inputs from atmosphere, thrufall, and snowfall
c
c 4/28/04 multiplied by 1-offrac?
c$$$        v_surfdep = 
c$$$     $          (ppt_irrig*(1-snowfrac)            ! atmos precip
c$$$     $        + (vmix_can(is,3)-vmix_ohoriz(is,8))*(1-snowfrac) ! through fall
c$$$     $        + snowmelt(is)*mru_area(is)*a_million*inch2m)*(1-offrac) ! snowmelt
c$$$c
c$$$c Compute volumes of surface inputs
c$$$c
c$$$        ppt2hill = v_ppt2gnd * (1-snowfrac)*(1-offrac) ! atmospheric ppt
c$$$        irrext2hill = v_irrext2gnd * (1-snowfrac)*(1-offrac) ! external irrigation
c$$$        irrsat2hill = v_irrsat2gnd * (1-snowfrac)*(1-offrac) ! atmospheric irrsat
c$$$        irrhyd2hill = v_irrhyd2gnd * (1-snowfrac)*(1-offrac) ! atmospheric irrhyd
c$$$        can2hill = (vmix_can(is,3)-vmix_ohoriz(is,8))*(1-snowfrac)*
c$$$     $       (1-offrac)                                       ! throughfall
c$$$        melt2hill = (snowmelt(is)*mru_area(is)*a_million*inch2m)*  ! snowmelt
c$$$     $       (1-offrac)
c
c note that v_surfdep is reduced by the uz_bypass volume in the next section - Not(2/16/05)
c
c On days of melt or other surface deposition:
c Compute the volume of water bypassing the unsaturated zone. Note that in topmod_chem,
c the amount bypassing the unsaturated zone is computed after removing hortonian
c overland flow, but before computing dunnian. The reasoning is that the SD and
c SBAR values indicating surface saturation and rejected infiltration are considered
c to describe the more homogeneous soils in the pedons and not necessarily the porous
c vertical preferenctial flow paths.
c
        if(qvpref(is).gt.0) then               ! if qvpref>0 then v_surfdep2 is also
           v_surfdep2 = v_surfdep-v_ofhort
           dep_ratio = v_surfdep2/v_surfdep ! equals 1.0 if no hortonian flow
           uz_bypass = qvpref(is) * mru_area(is)*a_million
           bypass_frac = uz_bypass/v_surfdep2
        else
           dep_ratio = 0.0
           uz_bypass = 0.0
           bypass_frac = 0.0
        end if

c
c Begin topographic index loop. Root zone and unsaturated zone combined.
c

c initialize accumulators

        do 150 ia = 1, nacsc(is)

c Area weights. Uses same topmod.f assumption that the first topographic
c index value is a maximum value with zero percent area having a greater
c value (i.e. st > max and ac = 0)
         if(ia.eq.nacsc(is)) then
           ACF=0.5*AC(IA,is)
         else
           ACF=0.5*(AC(IA,is)+AC(IA+1,is))
         endif
c
c Mixing volumes for water transpired to canopy weighted by ET from each
c wetness index bin.
c
           if(vmix_can(is,11).gt.0) then
              vmix_uz2can(ia,is) = transp_rate*sae_local(ia,is)*acf
     $             * covden *mru_area(is)*a_million
           else
              vmix_uz2can(ia,is) = 0.0
           end if
c$$$            vmix_uz2can(ia,is) = vmix_can(is,11)*sae_local(ia,is)
c$$$     $           /sae(is)
c
c The fluxes between the unsaturated zone and the saturated zone will be
c tracked as follows:
c 1. Root zone wetting from the saturated zone (srzwet) will be added to
c   any negative uz2sat values to compute a volume to be transfered from 
c   saturated zone. These fluxes will be added to the uz inputs and tracked
c   as outputs from the saturated zone but using the solute chemistry as
c   it was in the saturated zone at the beginning of the time step.
c 2. After mixing infiltration and saturated zone inputs, the recharge, 
c   quz, and any positive uz2sat values will be combined to compute the flux
c   of water from the uz to the saturated zone.
c 3. Finally the recharge will be mixed with saturated zone volume
c   considering the corrected volume for waters transferred to the uz
c   at the beginning of the time step.
c
        v_qwet_loc = srzwet(ia,is)* acf*mru_area(is) *a_million
        v_qwet = v_qwet + v_qwet_loc

c volumes of recharge and lateral macropore flow
        v_rech_loc = quz_local(ia,is)*acf*mru_area(is) *a_million
        v_qdf_loc = qdf_local(ia,is)*acf*mru_area(is) *a_million
c
c uz2sat and uz2sat_vol are positive if the water table rose (below
c the root zone) and negative if the water table was fell.
c
        uz2sat_vol(ia,is) = uz2sat(ia,is)*acf*mru_area(is) *a_million
c
c where uz2sat is a mass closure result for porewater exchange from unsaturated
c zone to the saturated zone resulting from falling or rising water table.
c As the water table rises, a positive value results as residual water
c below field capacity that was in the unsaturated zone now forms part
c of the saturated zone water. Conversely, a falling water table will
c result in a negative value as water that was in the saturated zone
c below field capacity now becomes part of the unsaturated zone. Therefore
c positive uz2sat values are included as outputs from the uz and negative
c values are recorded as outputs from the saturated zone (at concentrations
c of the beginning of the time step) and inputs into the unsaturated zone.
c

c
c Unsaturated zone inputs (5-precip, 17-19, irrigation,
c                          6-ET, 8-canopy throughfall,
c                          9-melt, 14- (srz_wet-uz2sat))
c
        v_infil_loc = uz_infil(ia,is)*acf *mru_area(is)*a_million

        if(v_infil_loc.gt.0.and.v_surfdep.gt.0.) then
           vmix_uz(ia,is,5) =
     $          v_infil_loc*ppt2hill/v_surfdep
           vmix_uz(ia,is,17) =
     $          v_infil_loc*irrsat2hill/v_surfdep
           vmix_uz(ia,is,18) =
     $          v_infil_loc*irrhyd2hill/v_surfdep
           vmix_uz(ia,is,19) =
     $          v_infil_loc*irrext2hill/v_surfdep
           vmix_uz(ia,is,8) =
     $          v_infil_loc*can2hill/v_surfdep
           vmix_uz(ia,is,9) =
     $          v_infil_loc*melt2hill/v_surfdep
        end if

        v_infil = v_infil + v_infil_loc
        vmix_uzgen(is,5)=vmix_uzgen(is,5)+vmix_uz(ia,is,5)
        vmix_uzgen(is,17)=vmix_uzgen(is,17)+vmix_uz(ia,is,17)
        vmix_uzgen(is,18)=vmix_uzgen(is,18)+vmix_uz(ia,is,18)
        vmix_uzgen(is,19)=vmix_uzgen(is,19)+vmix_uz(ia,is,19)
        vmix_uzgen(is,8)=vmix_uzgen(is,8)+vmix_uz(ia,is,8)
        vmix_uzgen(is,9)=vmix_uzgen(is,9)+vmix_uz(ia,is,9)
        if(riparian(ia,is))  then
            vmix_uzrip(is,5)=vmix_uzrip(is,5)+vmix_uz(ia,is,5)
            vmix_uzrip(is,17)=vmix_uzrip(is,17)+vmix_uz(ia,is,17)
            vmix_uzrip(is,18)=vmix_uzrip(is,18)+vmix_uz(ia,is,18)
            vmix_uzrip(is,19)=vmix_uzrip(is,19)+vmix_uz(ia,is,19)
            vmix_uzrip(is,8)=vmix_uzrip(is,8)+vmix_uz(ia,is,8)
            vmix_uzrip(is,9)=vmix_uzrip(is,9)+vmix_uz(ia,is,9)
        else
            vmix_uzup(is,5)=vmix_uzup(is,5)+vmix_uz(ia,is,5)
            vmix_uzup(is,17)=vmix_uzup(is,17)+vmix_uz(ia,is,17)
            vmix_uzup(is,18)=vmix_uzup(is,18)+vmix_uz(ia,is,18)
            vmix_uzup(is,19)=vmix_uzup(is,19)+vmix_uz(ia,is,19)
            vmix_uzup(is,8)=vmix_uzup(is,8)+vmix_uz(ia,is,8)
            vmix_uzup(is,9)=vmix_uzup(is,9)+vmix_uz(ia,is,9)
        endif
c
c ET - note that the water transpired to the canopy is drawn from
c the unsaturated zone as water with solutes and therefore appears
c as an output. The remaining ET tracked here will evaporate water
c only thereby concentrating the solutes left in the unsaturated
c zone.
c
        vmix_uz(ia,is,6) = 
     $    (sae_local(ia,is)*acf*mru_area(is)*a_million) ! total ET volume (sae in meters)
     $    - vmix_uz2can(ia,is)                          ! minus that translocated
                                                        ! to the canopy
c
c Collect evaporation in uzgen, mru, and basin variables
c
        vmix_uzgen(is,6)=vmix_uzgen(is,6)+vmix_uz(ia,is,6)
        if(riparian(ia,is)) then
          vmix_uzrip(is,6)=vmix_uzrip(is,6)+vmix_uz(ia,is,6)
        else  
          vmix_uzup(is,6)=vmix_uzup(is,6)+vmix_uz(ia,is,6)
        endif
        vmix_mru(is,6) = vmix_mru(is,6) + vmix_uz(ia,is,6)
        vmix_basin(6) = vmix_basin(6) + vmix_uz(ia,is,6)
c
c
c The flux uz2sat indicates the exchange of water from the unsaturated zone
c to the saturated zone resulting from a change in the water table. The
c unsaturated zone is defined as the wilting point moisture from the
c surface to the bottom of the root zone, the root zone moisture between
c wilting point and the field capacity for the same depth, and any readily
c draining water between field capacity and saturation that has not yet
c been delivered to the saturated zone (suz). Therefore, water table
c fluxtuations in the root zone will result in a uz2sat value of zero.
c negative as water below field capacity that was located in the saturated
c zone is now considered part of the unsaturated zone. Conversely, a rising
c water table below the root-zone depth will result in a positive uz2sat value
c as pore water that was included in the unsaturated zone is now in the saturated
c zone domain.
c
        if(uz2sat(ia,is).lt.0.0) then
           vmix_uz(ia,is,13) = v_qwet_loc-uz2sat_vol(ia,is)
           v_sat2uz = v_sat2uz - uz2sat_vol(ia,is)
        else
           vmix_uz(ia,is,13) = v_qwet_loc
        end if


c
c add uz inputs from saturated zone
c
        vmix_sat2uz(is) = vmix_sat2uz(is)+ vmix_uz(ia,is,13)
        vmix_uzgen(is,13)=vmix_uzgen(is,13)+vmix_uz(ia,is,13)
        if(riparian(ia,is)) then
          vmix_uzrip(is,13)=vmix_uzrip(is,13)+vmix_uz(ia,is,13)
        else  
          vmix_uzup(is,13)=vmix_uzup(is,13)+vmix_uz(ia,is,13)
        endif

c
c sum inputs
c
        vmix_uz(ia,is,2) = vmix_uz(ia,is,5) + vmix_uz(ia,is,8)
     $       + vmix_uz(ia,is,17) + vmix_uz(ia,is,18)
     $       + vmix_uz(ia,is,19) 
     $       + vmix_uz(ia,is,9) + vmix_uz(ia,is,13)
c
c collect inputs in uzgen
c
        vmix_uzgen(is,2)=vmix_uzgen(is,2)+vmix_uz(ia,is,2)
        if(riparian(ia,is)) then
          vmix_uzrip(is,2)=vmix_uzrip(is,2)+vmix_uz(ia,is,2)
        else  
          vmix_uzup(is,2)=vmix_uzup(is,2)+vmix_uz(ia,is,2)
        endif
c
c Unsaturated zone outputs (other than evaporation)
c
c Recharge to preferential flow and recharge weighted by recharge from
c each wetness index bin. Also add the volume of unsaturated zone pore
c water incorporated into the saturated zone from a rising water table.
c
!        vmix_uz2sat(ia,is) = v_rech_loc + v_qdf_loc  
        vmix_uz2sat(ia,is) = v_rech_loc
        vmix_uz2qdf(ia,is) = v_qdf_loc
        
        vmix_uz(ia,is,3) = vmix_uz2can(ia,is) +
     $       v_rech_loc + v_qdf_loc   !transpiration + local recharge + preferential flow

        if(uz2sat(ia,is).gt.0)  then
c track volume of recharge from UZ (not including bypass)
c   for reporting or debugging later
           vmix_uz2sat(ia,is) =
     $          vmix_uz2sat(ia,is) + uz2sat_vol(ia,is)
c               add volume engulfed by rising water to UZ output
           vmix_uz(ia,is,3) = vmix_uz(ia,is,3) +
     $        uz2sat_vol(ia,is)
           v_uz2sat = v_uz2sat + uz2sat_vol(ia,is)
        end if
c$$$        vmix_uz(ia,is,3) = vmix_uz2can(ia,is) + !transpiration
c$$$     $       vmix_uz2sat(ia,is) !plus recharge and volume from rising water table
c
c collect composite outputs
c
        vmix_uzgen(is,3)=vmix_uzgen(is,3)+vmix_uz(ia,is,3)
        if(riparian(ia,is)) then
          vmix_uzrip(is,3)=vmix_uzrip(is,3)+vmix_uz(ia,is,3)
        else  
          vmix_uzup(is,3)=vmix_uzup(is,3)+vmix_uz(ia,is,3)
        endif

c
c collect inputs for the preferential flow and the recharge from the
c unsaturated zone (source index 11)
c
        vmix_qdf(is,11) = vmix_qdf(is,11) + 
     $       v_qdf_loc
        vmix_sat(is,11) = vmix_sat(is,11) + 
     $       v_rech_loc

c track matrix recharge separate from uz-bypass or uz2sat volumes
        v_rech = v_rech + v_rech_loc

c add engulfed water as uz output to sat zone
        if(uz2sat(ia,is).gt.0.0) vmix_sat(is,11) =
     $       vmix_sat(is,11) + uz2sat_vol(ia,is)

c
c New storage for unsaturated zone
c
        vmix_uz(ia,is,4) = vmix_uz(ia,is,1) + vmix_uz(ia,is,2)
     $       - vmix_uz(ia,is,3) - vmix_uz(ia,is,6)

        v_sssto = v_sssto + vmix_uz(ia,is,4)

c
c Debugging output to chemout file
c
c Good thing I didn't try this with ia eq 3 because the acf would equal 0
c
!            if(ia.eq.23.and.is.eq.1)
!     $           write(chemout_file_unit,234)mru_ppt(is),
!     $           z_wt_local(ia,is), sae_local(ia,is),
!     $           sae(is),srz(ia,is),srzwet(ia,is),uz_depth(ia,is),
!     $           uz_infil(ia,is),quz_local(ia,is),uz2sat(ia,is),
!     $           vmix_uz(ia,is,1)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,2)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,3)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,4)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,5)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,6)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,8)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,9)/mru_area(is)/acf/a_million,
!     $           vmix_uz(ia,is,14)/mru_area(is)/acf/a_million,
!     $           vmix_uz2can(ia,is)/mru_area(is)/acf/a_million
! 233     format(' mru_ppt   z_wt_loc      sae_local sae       ',
!     $        'srz       srzwet    ',
!     $        'uz_depth  uz_infil  quz_local  uz2sat    ',
!     $        'vmix_uz01 vmix_uz02 vmix_uz03 vmix_uz04 vmix_uz05 ',
!     $        'vmix_uz06 vmix_uz08 vmix_uz09 vmix_uz14 vmix_uz2can')
! 234  format(20f10.5)
!

c$$$            if(ia.eq.10.and.is.eq.1)
c$$$     $          write(*,234)quz_local(ia,is),qdffrac(is),
c$$$     $          vmix_uz2sat(ia,is)/acf/mru_area(is)/a_million
c$$$ 233     format(' quz_loc   qdffrac    vmix_uz2sat')
c$$$ 234  format(3(e12.5,1X))

 150       continue                  ! end of loni loop

! Compute new composite volumes
           vmix_uzgen(is,4) = vmix_uzgen(is,1) + vmix_uzgen(is,2)
     $           - vmix_uzgen(is,3) - vmix_uzgen(is,6)
           vmix_uzrip(is,4) = vmix_uzrip(is,1) + vmix_uzrip(is,2)
     $           - vmix_uzrip(is,3) - vmix_uzrip(is,6)
           vmix_uzup(is,4) = vmix_uzup(is,1) + vmix_uzup(is,2)
     $           - vmix_uzup(is,3) - vmix_uzup(is,6)


c
c     Calculate shallow preferential flow fluxes
c
c Unsaturated zone preferential flow inputs (11 - unsaturated zone) calculated above
           vmix_qdf(is,2) = vmix_qdf(is,11)
c Unsaturated zone preferential flow outputs same as inputs
c 
           vmix_qdf(is,3) = vmix_qdf(is,11)
c
c New unsat zone preferential volume is static
c
           vmix_qdf(is,4) = vmix_qdf(is,1) + vmix_qdf(is,2)
     $          - vmix_qdf(is,3)
c
c Collect saturated zone inputs (5 - atmospheric precip; 8 - throughfall;
c 9 - snowmelt; 11 - unsaturated zone, 17 - irrsat, 18 - irrhyd,
c 19 - irrext, 20 - gw1, 21 - gw2)
c
           vmix_sat(is,5) = ppt2hill*dep_ratio*bypass_frac
           vmix_sat(is,17) = irrsat2hill*dep_ratio*bypass_frac
           vmix_sat(is,18) = irrhyd2hill*dep_ratio*bypass_frac
           vmix_sat(is,19) = irrext2hill*dep_ratio*bypass_frac
           vmix_sat(is,8) = can2hill*dep_ratio*bypass_frac
           vmix_sat(is,9) = melt2hill*dep_ratio*bypass_frac
c           vmix_sat(is,11) = recharge, calculated above
           vmix_sat(is,20) = gw_in1(is)*mru_area(is)*a_million
           vmix_sat(is,21) = gw_in2(is)*mru_area(is)*a_million

           vmix_sat(is,2) = vmix_sat(is,5) +
     $          vmix_sat(is,17)+vmix_sat(is,18)+
     $          vmix_sat(is,19)+ vmix_sat(is,8) +
     $          vmix_sat(is,9) + vmix_sat(is,11) +
     $          vmix_sat(is,20) + vmix_sat(is,21)

c
c Collect GW influx as mru and basin inputs
c

c gw1    (20)
           vmix_mru(is,20) = vmix_sat(is,20)
           vmix_basin(20) = vmix_basin(20)+ vmix_sat(is,20)
c gw2    (21)
           vmix_mru(is,21) = vmix_sat(is,21) 
           vmix_basin(21) = vmix_basin(21)+ vmix_sat(is,21)

           vmix_mru(is,2) = vmix_mru(is,2) + 
     $          vmix_sat(is,20) + vmix_sat(is,21)
           vmix_basin(2) = vmix_basin(2) +
     $          vmix_sat(is,20) + vmix_sat(is,21)
c
c Saturated zone outputs (deep preferential flow, 
c                       loss to deep aquifer, water removed for
c                       irrigation and stream/hill reservoir)
c
           vsat_tmp = vmix_sat(is,3)
           v_qpref = qpref(is)*mru_area(is)*a_million
           v_gw_loss = gw_loss(is)*mru_area(is)*a_million
           v_qb = qb(is)*mru_area(is)*a_million
           v_qwell = v_irrsat

           vmix_sat(is,3) = vmix_sat2uz(is)+
     $          v_qpref +              ! preferential flow
     $          v_qb +                 ! baseflow
     $          v_gw_loss +            ! loss to deep aquifer
     $          v_qwell                ! irrigation
c
c irrig_sat_mru was applied at the beginning of the time step.
c
c collect gw_loss and irrigation as mru and basin output
c
         vmix_mru(is,3) = vmix_mru(is,3) + v_gw_loss + v_qwell
         vmix_basin(3) = vmix_basin(3) + v_gw_loss + v_qwell

c debug
c           if(is.eq.1)then
c              write(*,295)vsat_tmp,vmix_sat(is,2),vmix_sat(is,3),
c     $          z_wt_local(10,is),
c     $          (uz2sat_net(ij,is),ij=1,11),
c     $          (srzwet(ia,is),ia=1,11),
c     $          (sd(ik,is),ik=1,11)
c           end if
c           if(datetime(1).eq.1993.and.
c     $          datetime(3).eq.13.and.is.eq.1) return
c /debug
c
c New saturated zone volume
c
           vmix_sat(is,4) =  vmix_sat(is,1) +  vmix_sat(is,2) 
     $          -  vmix_sat(is,3)
c
c Calculate deep preferential fluxes
c Inputs (13 - saturated zone); place in total inputs
c
           vmix_satpref(is,13) = v_qpref
           vmix_satpref(is,2) = vmix_satpref(is,13)
c
c Ouput equal volume to stream since volume is static
c
           vmix_satpref(is,3) = vmix_satpref(is,2)
c
c Sum to zero storage change
c
           vmix_satpref(is,4) = vmix_satpref(is,1) + vmix_satpref(is,2)
     $          - vmix_satpref(is,3)
c
c Calculate transient hillslope fluxes (zero residual storage)
c Inputs:
c    10) O-horizon
c    12) qdf, preferential flow in the unsaturated zone
c    13) qb, baseflow discharge from the saturated zone
c    14) exfiltration (inactive)
c    15) satpref, preferential flow in the saturated zone
c
c overland flow minus recharge to canopy on first day of transpiration
c
           vmix_hill(is,10) = vmix_ohoriz(is,3) - vmix_can(is,10)
c shallow preferential flow, qdf
           vmix_hill(is,12) = vmix_qdf(is,3)
c baseflow, qb
           vmix_hill(is,13) = v_qb
c exfiltration, qexfil : exfiltration in no longer simulated - RW
           vmix_hill(is,14) = 0
c deep preferential flow
           vmix_hill(is,15) = vmix_satpref(is,3)
c
c Sum hillslope inputs
c
           vmix_hill(is,2) = vmix_hill(is,10) + vmix_hill(is,12)
     $          + vmix_hill(is,13) + vmix_hill(is,14)
     $          + vmix_hill(is,15)
c
c Ouput the same
c
           vmix_hill(is,3) = vmix_hill(is,2)
c
c Record export from mru and set final volume
c
           vmix_mru(is,3) = vmix_mru(is,3) + vmix_hill(is,3)

           vmix_mru(is,4) = vmix_mru(is,1) + vmix_mru(is,2) -
     $        vmix_mru(is,3)- vmix_mru(is,6)

c
c Zero residual volume
c
           vmix_hill(is,4) = vmix_hill(is,1) + vmix_hill(is,2)
     $          - vmix_hill(is,3)
c
c Test basin discharge and ET
c
c           vmix_basin(3) = vmix_basin(3) + vmix_hill(is,3)
c
c Calculate basin water balance
c
c$$$            if(is.eq.1)
c$$$     $          write(*,234)
c$$$     $          qof(is),
c$$$     $          vmix_hill(is,10)/mru_area(is)/a_million,
c$$$     $          qdf(is),
c$$$     $          vmix_hill(is,12)/mru_area(is)/a_million,
c$$$     $          qb(is),
c$$$     $          vmix_hill(is,13)/mru_area(is)/a_million,
c$$$     $          vmix_hill(is,14)/mru_area(is)/a_million,
c$$$     $          qpref(is),
c$$$     $          vmix_hill(is,15)/mru_area(is)/a_million
c$$$            if(is.eq.4)
c$$$     $          write(*,234)
c$$$     $          (vmix_hill(is,j),j=1,10)
c$$$ 233        format('        qof     hill10        qdf     hill12',
c$$$     $             '         qb     hill13     qexfil     hill14',
c$$$     $             '      qpref     hill15')
c$$$ 234  format(10(1X,e10.2))

c This section accumulates the volume in the depth variable. The
c volumes will be converted to depths at the end of this run section

        basin_ppt_cm = basin_ppt_cm + v_ppt
        basin_irr_sat_cm = basin_irr_sat_cm + v_irrsat
        basin_irr_hyd_cm = basin_irr_hyd_cm + v_irrhyd
        basin_irr_ext_cm = basin_irr_ext_cm + v_irrext
        basin_intcp_sto_cm = basin_intcp_sto_cm + vmix_can(is,4)
        basin_intcp_cm = basin_intcp_cm + v_intcp
        basin_intcp_evap_cm = basin_intcp_evap_cm + vmix_can(is,6)
        basin_transp_cm = basin_transp_cm + v_transp
        basin_thruf_cm = basin_thruf_cm + v_thruf
        basin_net_rain_cm = basin_net_rain_cm + v_netrain
        basin_net_snow_cm = basin_net_snow_cm + v_netsnow
        basin_net_dep_cm = basin_net_dep_cm + v_netdep
        basin_pweqv_cm = basin_pweqv_cm + vmix_snow(is,4)
        basin_snowmelt_cm = basin_snowmelt_cm + vmix_snow(is,3)
        basin_snowevap_cm = basin_snowevap_cm + vmix_snow(is,6)
        basin_surfdep_cm = basin_surfdep_cm + v_surfdep
c soil moisture in the root zone not explicitly tracked in webmodres since
c it is assumed part of the overall unsaturated zone reservoir.
c So just assign the topmodel smav_basin to basin_soil_moist_cm
c at the end of this run section
        basin_sssto_cm = basin_sssto_cm + v_sssto
        basin_ofhort_cm = basin_ofhort_cm + v_ofhort
        basin_ofdunn_cm = basin_ofdunn_cm + v_ofdunn
        basin_sroff_cm = basin_sroff_cm + v_sroff
        basin_infil_cm = basin_infil_cm + v_infil
        basin_qdf_cm = basin_qdf_cm + vmix_qdf(is,3)
        basin_uz_et_cm = basin_uz_et_cm + v_uz_et
c$$$        write(*,*)'is,gw_sto,vmix_sat(is,4)',is,basin_gw_sto_cm,
c$$$     $       vmix_sat(is,4)
        basin_gw_sto_cm = basin_gw_sto_cm + vmix_sat(is,4)
        basin_gw1_in_cm = basin_gw1_in_cm + vmix_sat(is,20)
        basin_gw2_in_cm = basin_gw2_in_cm + vmix_sat(is,21)
        basin_qwet_cm = basin_qwet_cm + v_qwet
        basin_vpref_cm = basin_vpref_cm + uz_bypass
        basin_uz2sat_cm = basin_uz2sat_cm + v_uz2sat
        basin_recharge_cm = basin_recharge_cm + v_rech
        basin_sat2uz_cm = basin_sat2uz_cm + v_sat2uz
        basin_qwell_cm = basin_qwell_cm + v_qwell 
        basin_qpref_cm = basin_qpref_cm + v_qpref
        basin_gwloss_cm = basin_gwloss_cm + v_gw_loss
        basin_gwflow_cm = basin_gwflow_cm + v_qb
        basin_stflow_cm = basin_stflow_cm + vmix_hill(is,3)

 10   continue

c
c debug
c
c$$$      write(*,234) (vmix_stream(ih),ih=1,10)
c$$$ 233        format('        qof     hill10        qdf     hill12',
c$$$     $             '         qb     hill13     qexfil     hill14',
c$$$     $             '      qpref     hill15')
c$$$ 234  format(10(1X,e10.2))

c
c Calculate stream inputs and fluxes
c
      do 200 ih = 1, clark_segs
         do 20 is = 1,nmru
            vmix_hillexp(ih,is) = vmix_hill(is,3)
     $           *ar_fill(ih,mru2chan(is))
     $           /(chan_area(mru2chan(is))/basin_area)
c     add hillslope contribution to existing stream segment volume
c     and remove any_diversions or seepage through the stream bed
            vmix_stream(ih) =
     $           vmix_stream(ih)+vmix_hillexp(ih,is)

 20      continue
         vmix_diversion(ih) = irrig_hyd_seg(ih) ! for clarity
         vmix_chan_loss(ih) = chan_loss(ih)*basin_area*a_million 
c
c appears that chan_loss has already bean removed and that vmix_stream
c represents the residual
c
c$$$         vmix_stream(ih) = vmix_stream(ih) -
c$$$     $        vmix_diversion(ih)-vmix_chan_loss(ih)
c$$$         vmix_stream(ih) = vmix_stream(ih) - vmix_diversion(ih)
         basin_chan_div_cm = basin_chan_div_cm + vmix_diversion(ih)
         basin_chan_loss_cm = basin_chan_loss_cm + vmix_chan_loss(ih)
         vmix_basin(3) = vmix_basin(3) + vmix_chan_loss(ih)+ 
     $                   vmix_diversion(ih)
      
 200  continue
c Export hydro-segment nearest the outlet and shift such that
c the segment farthest from the outlet is set to zero.
      vmix_basin(3) = vmix_basin(3) + vmix_stream(1)
      basin_qsim_cm = vmix_stream(1)

      do 17  j = 1,clark_segs-1
         vmix_stream(j)=vmix_stream(j+1)
         basin_chan_sto_cm = basin_chan_sto_cm + vmix_stream(j)
 17   continue
      vmix_stream(clark_segs)=0.0

      vmix_basin(4) = vmix_basin(1) + vmix_basin(2) - vmix_basin(3)
     $     - vmix_basin(6)
c
c Convert volumes to average basin depths for summary module
c
      basin_sto_cm = vmix_basin(4)*m3cm
      basin_in_cm = vmix_basin(2) * m3cm
      basin_out_cm = vmix_basin(3) * m3cm
      basin_et_cm = vmix_basin(6) * m3cm
      basin_ppt_cm =  basin_ppt_cm * m3cm
      basin_irr_ext_cm = basin_irr_ext_cm * m3cm
      basin_irr_sat_cm = basin_irr_sat_cm * m3cm
      basin_irr_hyd_cm = basin_irr_hyd_cm * m3cm
      basin_intcp_sto_cm = basin_intcp_sto_cm * m3cm
      basin_intcp_cm = basin_intcp_cm * m3cm
      basin_intcp_evap_cm = basin_intcp_evap_cm * m3cm
      basin_transp_cm = basin_transp_cm * m3cm
      basin_thruf_cm = basin_thruf_cm * m3cm
      basin_net_dep_cm = basin_net_dep_cm * m3cm
      basin_net_rain_cm = basin_net_rain_cm * m3cm
      basin_net_snow_cm = basin_net_snow_cm * m3cm
      basin_pweqv_cm = basin_pweqv_cm * m3cm
      basin_snowmelt_cm = basin_snowmelt_cm * m3cm
      basin_snowevap_cm = basin_snowevap_cm * m3cm
      basin_surfdep_cm = basin_surfdep_cm * m3cm
      basin_soil_moist_cm = smav_basin * 100   ! convert to cm and pass on
      basin_sssto_cm = basin_sssto_cm * m3cm
      basin_ofhort_cm = basin_ofhort_cm * m3cm
      basin_ofdunn_cm = basin_ofdunn_cm * m3cm
      basin_sroff_cm = basin_sroff_cm *m3cm
      basin_infil_cm = basin_infil_cm * m3cm
      basin_qdf_cm = basin_qdf_cm * m3cm
      basin_uz_et_cm = basin_uz_et_cm * m3cm
      basin_gw_sto_cm = basin_gw_sto_cm * m3cm
      basin_gw1_in_cm = basin_gw1_in_cm * m3cm
      basin_gw2_in_cm = basin_gw2_in_cm * m3cm
      basin_qwet_cm = basin_qwet_cm * m3cm
      basin_vpref_cm = basin_vpref_cm * m3cm
      basin_uz2sat_cm = basin_uz2sat_cm * m3cm
      basin_recharge_cm = basin_recharge_cm * m3cm
      basin_sat2uz_cm = basin_sat2uz_cm * m3cm
      basin_qwell_cm = basin_qwell_cm * m3cm
      basin_qpref_cm = basin_qpref_cm * m3cm
      basin_gwloss_cm = basin_gwloss_cm * m3cm
      basin_gwflow_cm = basin_gwflow_cm * m3cm
      basin_chan_sto_cm = basin_chan_sto_cm * m3cm
      basin_stflow_cm = basin_stflow_cm * m3cm
      basin_chan_div_cm = basin_chan_div_cm * m3cm
      basin_chan_loss_cm = basin_chan_loss_cm * m3cm

      basin_qsim_cm = basin_qsim_cm * m3cm
      basin_qsim_m3s = basin_qsim_cm*basin_area/8.64/fac

c
c Write output volume files when print_vse >= 1
c
      if(print_vse.ge.1) then
       write(vf_bas%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_basin(i),i=1,21)
        do is=1,nmru
          write(vf_mru(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_mru(is,i),i=1,21)
          if(print_vse.eq.2) then
          write(vf_uzgen(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_uzgen(is,i),i=1,21)
          write(vf_uzrip(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_uzrip(is,i),i=1,21)
          write(vf_uzup(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_uzup(is,i),i=1,21)
          write(vf_can(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_can(is,i),i=1,21)
          write(vf_snow(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_snow(is,i),i=1,21)
!          write(vf_imperv(is)%lun,123) nstep,(datetime(i),i=1,3),
!     $      (vmix_imp(is,i),i=1,21)
          write(vf_transp(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_uz2can(i,is),i=1,nac)
          write(vf_ohoriz(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_ohoriz(is,i),i=1,21)
          write(vf_qdf(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_qdf(is,i),i=1,21)
          write(vf_sat(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_sat(is,i),i=1,21)
          write(vf_satpref(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_satpref(is,i),i=1,21)
          write(vf_hill(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_hill(is,i),i=1,21)
          write(vf_uz2sat(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_uz2sat(i,is),i=1,nac)
            do ia=1,nac
              write(vf_uz(is,ia)%lun,123) nstep,(datetime(i),i=1,3),
     $          (vmix_uz(ia,is,i),i=1,21)
            end do
          end if
        end do
        do is=1,nhydro
          write(vf_hillexp(is)%lun,123) nstep,(datetime(i),i=1,3),
     $      (vmix_hillexp(is,i),i=1,nmru)
        end do
        write(vf_hyd%lun,123) nstep,(datetime(i),i=1,3),
     $     basin_qsim_cm/m3cm, (vmix_stream(i),i=1,nhydro)
      end if
! temp print out  uz vols
       !if(end_run) then
       !    write(debug%lun,*)'O Horiz:',vmix_ohoriz(1,4)
       !    write(debug%lun,*)'UZ_pref:',vmix_qdf(1,4)
       !    write(debug%lun,*)'Sat:',vmix_sat(1,4)
       !    write(debug%lun,*)'UZ volumes_m3:',(vmix_uz(ia,1,4),ia=1,11)
       !endif
!        
 123  format(2I5,2I3,21(E14.6,1X))
!      if (end_dy) write(*,"(85E14.6)")
!     $ (vmix_can(1,i), i=1,6),(vmix_ohoriz(1,i), i=1,6),(vmix_snow(1,i),
!     $ i=1,6),
!     $ (vmix_uzgen(1,i), i=1,6),(vmix_qdf(1,i), i=1,6),(vmix_sat(1,i),
!     $  i=1,6),(vmix_satpref(1,i), i=1,6),(vmix_mru(1,i), i=1,6),
!     $ (vmix_stream(i), i=1,4),(vmix_basin(i), i=1,6),(vmix_hill(1,i),
!     $ i=1,15),v_ofhort,v_ofdunn,v_sroff,v_dep
!
c
c End of run
c
c$$$      if (end_run) print*,'Total depth of rainfall was ',
c$$$     $     (vmix_basin(1) - vmix_basin0)/basin_area/a_million,' ',
c$$$     $     vmix_basin0
      webrrun = 0

      return
 295  format(79(e10.4,1X))
      end

c
c No cleanup routine needed since output files are closed in the
c io_chem module
c

