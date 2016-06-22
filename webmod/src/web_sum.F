#include "defines.h"
**********************************************************************
* Modified from PRMS version by Rick Webb
*
* 5/15/02 - Corrected if/endif around year loop so that the user
*           can output only run totals
*
* 10/6/03 - Added column for deep preferential flow.
*           Changed gwdef to gwsto since model now has an absoute
*            volume dictated by depth to bedrock.
*
* 10/8/04 - Summations or min/max over time periods are now
*           consolidated into a single <var>_sum variable dimensioned by
*           five. The first index is a storm index (not currently
*           used except to mirror the unit values. This is followed by
*           daily sums <var>_sum(2) (previously _dy variables), then
*           monthly(3), yearly(4), and run totals(5). Summation vars
*           are now zeroed at the beginning of the run section if the
*           end of a time period occured on the last step.
*
*           To keep with this new scheme, the sum_obj_fun(5) will
*           now be obj_fun_sum(5,5) with the first index the objective
*           function type and the second index the time period.
*
*           Added irrigation and losses to deep aquifers
*
* 11/4/04 - Modified the way that average values are accumulated.
*           The time steps in the storm, day, month, year are tracked.
*           The cumulative total for each period is divided by the time
*           steps to compute a period average. Also added average solar
*           radiation values using this approach.
*
*11/26/04 - Moved computation of unit values to the webmod_res module
* 
*03/04/05 - Changed Basin min and max temperature to summarize the
*           area-weighted average. 
*
c 17apr09 - Add Fortran90 Module: WEBMOD_SUM
c
c Jan 2016 - Added chemical summaries
c
c***********************************************************************

**********************************************************************
c     Sums values for daily, monthly, yearly and total flow
c
c***********************************************************************

c***********************************************************************
c! ***************** START MODULE *************
      MODULE WEBMOD_SUM
      IMPLICIT NONE
#include "fmodules.inc"
      
C   Parameters and initial flags
c Conversion factors
c Ld2Wm2, Convert Langleys/day to Watts/sq.meter
c m3cm, Convert basin volumes in cubic meters to depth in centimeters
c       given the basin area (set in the  init)

      real, PARAMETER :: Ld2Wm2 = 0.4833333333333
      real m3cm

C   Dimensions and flags
      integer, save :: nmru
      logical, save :: uprt, dprt, mprt, yprt, tprt

C   Declared Parameters
      integer, save :: print_type, print_freq, print_objfunc
      integer :: header_prt, print_explanation
      real, save ::  dtinit, basin_area
C   Declared Variables
      integer, save :: lastper(5), tot_steps(5)
      real, save ::  swrad_W,swrad_W_avg(5),tot_swrad_W(5)
      real, save ::  basin_tmaxc_max(5), basin_tminc_min(5)
      real, save ::  basin_potet_cm, basin_potet_cm_sum(5)
      real, save ::  basin_in_cm_sum(5)
      real, save ::  basin_out_cm_sum(5)
      real, save ::  basin_et_cm_sum(5)
      real, save ::  basin_ppt_cm_sum(5)
      real, save ::  basin_irr_ext_cm_sum(5)
      real, save ::  basin_irr_sat_cm_sum(5)
      real, save ::  basin_irr_hyd_cm_sum(5)
      real, save ::  basin_intcp_cm_sum(5)
      real, save ::  basin_intcp_evap_cm_sum(5)
      real, save ::  basin_transp_cm_sum(5)
      real, save ::  basin_thruf_cm_sum(5)
      real, save ::  basin_net_dep_cm_sum(5)
      real, save ::  basin_net_rain_cm_sum(5)
      real, save ::  basin_net_snow_cm_sum(5)
      real, save ::  basin_snowmelt_cm_sum(5)
      real, save ::  basin_snowevap_cm_sum(5)
      real, save ::  basin_surfdep_cm_sum(5)
      real, save ::  basin_ofhort_cm_sum(5)
      real, save ::  basin_ofdunn_cm_sum(5)
      real, save ::  basin_sroff_cm_sum(5)
      real, save ::  basin_infil_cm_sum(5)
      real, save ::  basin_qdf_cm_sum(5)
      real, save ::  basin_uz_et_cm_sum(5)
      real, save ::  basin_gw1_in_cm_sum(5)
      real, save ::  basin_gw2_in_cm_sum(5)
      real, save ::  basin_qwet_cm_sum(5)
      real, save ::  basin_vpref_cm_sum(5)
      real, save ::  basin_uz2sat_cm_sum(5)
      real, save ::  basin_recharge_cm_sum(5)
      real, save ::  basin_sat2uz_cm_sum(5)
      real, save ::  basin_qwell_cm_sum(5)
      real, save ::  basin_qpref_cm_sum(5)
      real, save ::  basin_gwloss_cm_sum(5)
      real, save ::  basin_gwflow_cm_sum(5)
      real, save ::  basin_stflow_cm_sum(5)
      real, save ::  basin_chan_div_cm_sum(5)
      real, save ::  basin_chan_loss_cm_sum(5)
      real, save ::  basin_qsim_cm_sum(5)
      real, save ::  basin_qobs_cm_sum(5)
      real, save ::  basin_qsim_m3s_avg(5)
      real, save ::  basin_qobs_m3s_avg(5)
      real, save ::  basin_init_m3_tot(5)
      real, save ::  wat_bal_cm,wat_bal_cm_sum(5)
      real, save ::  obj_func(5),obj_func_sum(5,5)
      real, save ::  tot_basin_qsim_m3s(5),tot_basin_qobs_m3s(5)
      real, save :: last_basin_stor,last_chan_stor
      character (len=:), save, allocatable :: bufferclabels ! chemout buffer for labels
      character (len=:), save, allocatable :: buffercspan  ! chemout buffer for header and span lines
      character (len=:), save, allocatable :: buffercdat  ! chemout buffer for data lines

      character (len=100) :: span1, spansol

c total discharge of water and solutes at basin outlet
      double precision, save :: ch_outlet_m3_sum(5)
      double precision, save :: ch_outlet_tempC_avg(5)
      double precision, save :: ch_outlet_pH_avg(5)
      double precision, save, allocatable :: ch_outlet_permil_avg(:,:)
      double precision, save, allocatable :: ch_outlet_M_sum(:,:)

c in, out, ET or rxn, final for water and solutes
      double precision, save :: ch_basin_m3_in_sum(5)
      double precision, save :: ch_basin_m3_out_sum(5)
      double precision, save :: ch_basin_m3_ET_sum(5)
      double precision, save :: ch_basin_tempC_in_avg(5)
      double precision, save :: ch_basin_tempC_out_avg(5)
      double precision, save :: ch_basin_pH_in_avg(5)
      double precision, save :: ch_basin_pH_out_avg(5)
      double precision, save, allocatable:: ch_basin_permil_in_avg(:,:)
      double precision, save, allocatable:: ch_basin_permil_out_avg(:,:)
      double precision, save, allocatable:: ch_basin_permil_ET_avg(:,:)
      double precision, save, allocatable :: ch_basin_M_in_sum(:,:)
      double precision, save, allocatable :: ch_basin_M_out_sum(:,:)
      double precision, save, allocatable :: ch_basin_M_rxn_sum(:,:)

      END MODULE WEBMOD_SUM

c***********************************************************************
c
c     Main basin_sum routine
c

      integer function web_sum(arg)
      IMPLICIT NONE

      character(*) arg
      CHARACTER(256) SVN_ID

      integer sumbdecl, sumbinit, sumbrun

      save SVN_ID

      SVN_ID = 
     $     '$Id: web_sum.f 39 2007-06-08 18:01:49Z rmwebb $ '

      web_sum = 0

      if(arg.eq.'declare') then
        web_sum = sumbdecl()

      else if(arg.eq.'initialize') then
        web_sum = sumbinit()

      else if(arg.eq.'run') then
        web_sum = sumbrun()

      end if

C******Debug level print
c      call dpint4('End of basin_sum, retval = ', web_sum, 1, 2)
      return
      end
      
c***********************************************************************
c
c     sumbdecl - set up basin summary parameters
c
      integer function sumbdecl()

      USE WEBMOD_SUM
      USE WEBMOD_IO, only: debug
      USE WEBMOD_PHREEQ_MMS, only: nsolute

      IMPLICIT NONE
      sumbdecl = 1
      
      if(declparam('sumb', 'print_explanation', 'one', 'integer', 
     +   '1', '0', '1',
     +   'Print detailed explanation in output file (0,No;1,Yes)',
     +   'Print detailed explanation in output file (0,No;1,Yes)',
     +   'none').ne.0)return

      if(declparam('sumb', 'print_type', 'one', 'integer', 
     +   '1', '0', '2',
     +   'Type of output data file',
     +   'Output data file: 0 = observed and predicted flow only;'//
     +   '                  1 = water balance table;            '//
     +   '                  2 = detailed output.                ',
     +   'none').ne.0)return

      if(declparam('sumb', 'print_freq', 'one', 'integer',
     +   '1', '0', '31',
     +   'Frequency for output data file',
     +   'Output data file: 0 = no output file;                 '//
     +   '                  1 = output run totals;              '//
     +   '                  2 = output yearly totals;           '//
     +   '                  4 = output monthly totals;          '//
     +   '                  8 = output daily totals;            '//
     +   '                 16 = output unit totals;             '//
     +   'For combinations, add index numbers, for example, daily'//
     +   ' plus yearly output = 10; yearly plus total = 3',
     +   'none').ne.0)return      

      if(declparam('sumb', 'print_objfunc', 'one', 'integer',
     +   '0', '0', '1',
     +   'Switch to turn objective function printing off and on',
     +   'Output data file: 0 = no print of objective function;'//
     +   '                  1 = print objective function  ',
     +   'none').ne.0)return

      if(declparam('topc', 'dtinit', 'one', 'real',
     +   '24', '0', '24',
     +   'Initial timestep for initialize function.',
     +   'Initial timestep for initialize function.',
     +   'hours').ne.0) return

c      if(declparam('sumb', 'qobsta', 'one', 'integer',
c     +   '1', 'bounded', 'nobs',
c     +   'Index of streamflow station for calculating '//
c     +   'objective function.','Index of streamflow station '//
c     +   'for calculating objective function.','none').ne.0) return

      if(declparam('precip', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

c total discharge of water and solutes at basin outlet

      if(declvar('phreeqmms', 'ch_outlet_m3_sum', 'five', 5,
     +  'double','Total volume of water discharged at outlet '//
     $   'for time step of (1) unit, (2) day, (3) month, '//
     $   '(4) year, or (5) entire run.',
     +  'm3', ch_outlet_m3_sum).ne.0) return
      ch_outlet_m3_sum = 0.0d0

      if(declvar('phreeqmms', 'ch_outlet_tempC_avg', 'five', 5,
     +  'double','Temperature of water discharged at outlet '//
     $  'for time step of (1) unit, (2) day, (3) month, '//
     $  '(4) year, or (5) entire run.',
     +  'deg C', ch_outlet_tempC_avg).ne.0) return
      ch_outlet_tempC_avg = 0.0d0

      if(declvar('phreeqmms', 'ch_outlet_pH_avg', 'five', 5,
     +  'double','pH of water discharged at outlet '//
     $  'for time step of (1) unit, (2) day, (3) month, '//
     $  '(4) year, or (5) entire run.',
     +  'pH units', ch_outlet_pH_avg).ne.0) return
      ch_outlet_pH_avg = 0.0d0

      allocate(ch_outlet_permil_avg(nsolute,5))
      if(declvar('phreeqmms', 'ch_outlet_permil_avg', 'nsolute,five', 
     +  nsolute*5,'double','Average delta of isotope discharged at '//
     $  'outlet for time step of (1) unit, (2) day, (3) month, '//
     $  '(4) year, or (5) entire run.',
     +  'permil', ch_outlet_permil_avg).ne.0) return
      ch_outlet_permil_avg = 0.0d0

      allocate(ch_outlet_M_sum(nsolute,5))
      if(declvar('phreeqmms', 'ch_outlet_M_sum', 'nsolute,five', 
     +  nsolute*5,'double','Moles of solute discharged at outlet '//
     $  'for time step of (1) unit, (2) day, (3) month, '//
     $  '(4) year, or (5) entire run.',
     +  'Moles', ch_outlet_M_sum).ne.0) return
      ch_outlet_M_sum = 0.0d0

c in, out, ET or rxn, final for water and solutes
      
      if(declvar('phreeqmms', 'ch_basin_m3_in_sum', 'five', 5,
     +     'double','Total volume of water entering the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'm3', ch_basin_m3_in_sum).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_m3_out_sum', 'five', 5,
     +      'double','Total volume of water leaving the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'm3', ch_basin_m3_out_sum).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_m3_ET_sum', 'five', 5,'double',
     +      'Total volume of evapotranspiration from the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'm3', ch_basin_m3_ET_sum).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_tempC_in_avg', 'five', 5,
     +     'double','Average temperature of water entering the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'degC', ch_basin_tempC_in_avg).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_tempC_out_avg', 'five', 5,
     +     'double','Average temperature of water leaving the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'degC', ch_basin_tempC_out_avg).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_pH_in_avg', 'five', 5,
     +     'double','Average pH of water entering the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'pH units', ch_basin_pH_in_avg).ne.0) return

      if(declvar('phreeqmms', 'ch_basin_pH_out_avg', 'five', 5,
     +     'double','Average pH of water leaving the basin'//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +      'pH units', ch_basin_pH_out_avg).ne.0) return

      allocate(ch_basin_permil_in_avg(nsolute,5))
      if(declvar('phreeqmms','ch_basin_permil_in_avg','nsolute,five',
     +     nsolute*5,'double','Average delta of water entering the '//
     $     'basin for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'permil', ch_basin_permil_in_avg).ne.0) return

      allocate(ch_basin_permil_out_avg(nsolute,5))
      if(declvar('phreeqmms','ch_basin_permil_out_avg','nsolute,five',
     +     nsolute*5,'double','Average delta of water leaving the '//
     $     'basin for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'permil', ch_basin_permil_out_avg).ne.0) return

      allocate(ch_basin_permil_ET_avg(nsolute,5))
      if(declvar('phreeqmms', 'ch_basin_permil_ET_avg', 'nsolute,five',
     +      nsolute*5, 'double','Average delta of water evaporated '//
     $     'from the basin for time step of (1) unit, (2) day, '//
     $     '(3) month, (4) year, or (5) entire run.',
     +     'permil', ch_basin_permil_ET_avg).ne.0) return
      
      allocate(ch_basin_M_in_sum(nsolute,5))
      if(declvar('phreeqmms', 'ch_basin_M_in_sum', 'nsolute,five',
     +      nsolute*5,'double','Moles of solute entering the basin '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'Moles', ch_basin_M_in_sum).ne.0) return

      allocate(ch_basin_M_out_sum(nsolute,5))
      if(declvar('phreeqmms', 'ch_basin_M_out_sum', 'nsolute,five',
     +     nsolute*5,'double','Moles of solute leaving the basin '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'Moles', ch_basin_M_out_sum).ne.0) return

      allocate(ch_basin_M_rxn_sum(nsolute,5))
      if(declvar('phreeqmms', 'ch_basin_M_rxn_sum', 'nsolute,five',
     +  nsolute*5, 'double','Moles of solute produced (+) or '//
     +     'consumed (-) by reactions in the basin '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'Moles', ch_basin_M_rxn_sum).ne.0) return
c
c computed incoming shortwave radiation in Watts per square meters
c
      if(declvar('sumb', 'swrad_W', 'one', 1, 'real', 
     + 'Basin average incident shortwave radiation on a '//
     $ 'horizontal surface.',
     + 'W/m2', swrad_W).ne.0) return

      if(declvar('sumb', 'swrad_W_avg', 'five', 5, 'real',
     +     'Average shortwave radiation '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'W/m2',
     +   swrad_W_avg).ne.0) return

c temperature extremes

      if(declvar('sumb', 'basin_tmaxc_max', 'five', 5, 'real',
     +     'Maximum mean basin temp observed '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'deg C',
     +   basin_tmaxc_max).ne.0) return

      if(declvar('sumb', 'basin_tminc_min', 'five', 5, 'real',
     +     'Minimum mean basin temp observed '//
     $     'for time step of (1) unit, (2) day, (3) month, '//
     $     '(4) year, or (5) entire run.',
     +     'deg C',
     +   basin_tminc_min).ne.0) return

c Track potential ET directly as reported from the potet module

      if(declvar('sumb', 'basin_potet_cm', 'one', 1, 'real',
     +    'Basin potential ET ',
     +    'cm',
     +   basin_potet_cm).ne.0) return

      if(declvar('sumb', 'basin_potet_cm_sum', 'five', 5, 'real',
     +    'Sum of potential ET '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_potet_cm_sum).ne.0) return

c  -----------------------------
c  The following state variables are retrieved
c  from webmod_res with no summary computations:
c
c  basin_sto_cm, basin_intcp_sto_cm, basin_pweqv_cm,
c  basin_soil_moist_cm, basin_sssto_cm, basin_gw_sto_cm,
c  basin_chan_sto_cm
c  ------------------------
c  Unit values for fluxes and stores, expressed as average
c  basin depth in centimeters, are computed in webmod_res.
c  Only the summary vars_sum (sum and avg) will be declared
c  in this module. 
c
      if(declvar('sumb', 'basin_in_cm_sum', 'five', 5, 'real',
     +    'Sum of basin inputs from precip, irrigation, '//
     $    'and regional groundwater '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_in_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_out_cm_sum', 'five', 5, 'real',
     +    'Sum of basin discharge and losses to deep aquifers '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_out_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_et_cm_sum', 'five', 5, 'real',
     +    'Sum of basin evapotranspiration '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_et_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_ppt_cm_sum', 'five', 5, 'real',
     +    'Sum of basin precipitation '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_ppt_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_irr_ext_cm_sum', 'five', 5, 'real',
     +    'Sum of basin irrigation from external sources '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_irr_ext_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_irr_sat_cm_sum', 'five', 5, 'real',
     +    'Sum of basin irrigation from shallow wells '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_irr_sat_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_irr_hyd_cm_sum', 'five', 5, 'real',
     +    'Sum of basin irrigation from stream diversions '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_irr_hyd_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_intcp_cm_sum', 'five', 5, 'real',
     +    'Sum of basin interception '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_intcp_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_intcp_evap_cm_sum', 'five', 5, 'real',
     +    'Sum of evaporation from interception '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_intcp_evap_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_transp_cm_sum', 'five', 5, 'real',
     +    'Sum of transpiration from root zone '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_transp_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_thruf_cm_sum', 'five', 5, 'real',
     +    'Sum of canopy throughfall '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_thruf_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_net_dep_cm_sum', 'five', 5, 'real',
     +    'Sum of basin direct deposition plus throughfall '//
     +    '(precip+irrig) '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_net_dep_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_net_rain_cm_sum', 'five', 5, 'real',
     +    'Sum of basin net rain '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_net_rain_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_net_snow_cm_sum', 'five', 5, 'real',
     +    'Sum of basin net snow '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_net_snow_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_snowmelt_cm_sum', 'five', 5, 'real',
     +    'Sum of basin snowmelt '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_snowmelt_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_snowevap_cm_sum', 'five', 5, 'real',
     +    'Sum of basin snow evaporation '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_snowevap_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_surfdep_cm_sum', 'five', 5, 'real',
     +    'Sum of basin hillslope deposition '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_surfdep_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_ofhort_cm_sum', 'five', 5, 'real',
     +    'Sum of basin infiltration excess '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_ofhort_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_ofdunn_cm_sum', 'five', 5, 'real',
     +    'Sum of basin saturation excess '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_ofdunn_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_sroff_cm_sum', 'five', 5, 'real',
     +    'Sum of basin overland flow '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_sroff_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_infil_cm_sum', 'five', 5, 'real',
     +    'Sum of basin infiltration '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_infil_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qdf_cm_sum', 'five', 5, 'real',
     +    'Sum of basin lateral macropore flow '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qdf_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_uz_et_cm_sum', 'five', 5, 'real',
     +    'Sum of basin evaporation from soil '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_uz_et_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_gw1_in_cm_sum', 'five', 5, 'real',
     +    'Sum of inputs from leaky canals '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_gw1_in_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_gw2_in_cm_sum', 'five', 5, 'real',
     +    'Sum of inputs from upgradient groundwater '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_gw2_in_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qwet_cm_sum', 'five', 5, 'real',
     +    'Sum of basin root-zone replenishment by groundwater '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qwet_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_vpref_cm_sum', 'five', 5, 'real',
     +    'Sum of basin infiltration delivered directly to '//
     $    'groundwater for time step of (1) unit, (2) day, '//
     $    '(3) month, (4) year, or (5) entire run.',
     +    'cm',
     +   basin_vpref_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_uz2sat_cm_sum', 'five', 5, 'real',
     +    'Sum of basin UZ water engulfed into saturated zone '//
     $    'as water table rises '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_uz2sat_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_recharge_cm_sum', 'five', 5, 'real',
     +    'Sum of basin recharge to groundwater from UZ '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_recharge_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_sat2uz_cm_sum', 'five', 5, 'real',
     +    'Sum of basin water left in UZ as water table lowers '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_sat2uz_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qwell_cm_sum', 'five', 5, 'real',
     +    'Sum of basin water pumped from shallow wells '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qwell_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qpref_cm_sum', 'five', 5, 'real',
     +    'Sum of basin preferential flow through the saturated '//
     $    'zone for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qpref_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_gwloss_cm_sum', 'five', 5, 'real',
     +    'Sum of basin groundwater loss to deep aquifer '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_gwloss_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_gwflow_cm_sum', 'five', 5, 'real',
     +    'Sum of basin baseflow '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_gwflow_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_stflow_cm_sum', 'five', 5, 'real',
     +    'Sum of basin hillslope runoff to channels '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_stflow_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_chan_div_cm_sum', 'five', 5, 'real',
     +    'Sum of basin diversions for irrigation '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_chan_div_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_chan_loss_cm_sum', 'five', 5, 'real',
     +    'Sum of basin chan_lossersions for irrigation '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_chan_loss_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qsim_cm_sum', 'five', 5, 'real',
     +    'Sum of basin simulated discharge '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qsim_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qobs_cm_sum', 'five', 5, 'real',
     +    'Sum of basin observed discharge (at station qobsta) '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   basin_qobs_cm_sum).ne.0) return

      if(declvar('sumb', 'basin_qsim_m3s_avg', 'five', 5, 'real',
     +    'Average basin predicted runoff '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'm^3/s',
     +   basin_qsim_m3s_avg).ne.0) return

      if(declvar('sumb', 'basin_qobs_m3s_avg', 'five', 5, 'real',
     +    'Average basin observed discharge (at station qobsta) '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'm^3/s',
     +   basin_qobs_m3s_avg).ne.0) return

      if(declvar('sumb', 'obj_func', 'five', 5, 'real',
     +     'Objective function for each time step.\n'//
     +     'Index 1: |Obs-Pred.|\n'//
     +     '      2: (Obs-Pred)^2.\n'//
     +     '      3: |(ln(obs+1)-ln(Pred+1)|.\n'//
     +     '      4: (ln(Obs+1)-ln(Pred+1))^2.\n'//
     +     '      5: Obs-Pred.','none', obj_func).ne.0) return

      if(declvar('sumb', 'obj_func_sum', 'five,five', 5*5, 'real',
     +     'Objective function for each time period.\n'//
     +     'Index 1: |Obs-Pred.|'//
     +     '      2: (Obs-Pred)^2.'//
     +     '      3: |(ln(obs+1)-ln(Pred+1)|'//
     +     '      4: (ln(Obs+1)-ln(Pred+1))^2'//
     +     '      5: Obs-Pred.'//
     +     'Index 1: storm \n'//
     +     '      2: day \n'//
     +     '      3: month \n'//
     +     '      4: year \n'//
     +     '      5: run',
     $     'none', obj_func_sum).ne.0) return

      if(declvar('sumb', 'wat_bal_cm', 'one', 1, 'real', 
     + 'Error in basin water balance',
     +  'cm', wat_bal_cm).ne.0) return

      if(declvar('sumb', 'wat_bal_cm_sum', 'five', 5, 'real',
     +    'Cumulative error in basin water balance '//
     $    'for time step of (1) unit, (2) day, (3) month, '//
     $    '(4) year, or (5) entire run.',
     +    'cm',
     +   wat_bal_cm_sum).ne.0) return

      if(declpri('sumb_last_basin_stor', 1, 'real', last_basin_stor)
     + .ne. 0) return 

      if(declpri('sumb_last_chan_stor', 1, 'real', last_chan_stor)
     + .ne. 0) return 
c
c lastper flags (1) whether the last time step was an end of
c   (1) storm, (2) day, (3) month, and (4) year. Needed so that
c   last value of time period is not reset before leaving routine.
c
      if(declpri('sumb_lastper', 5, 'integer', lastper).ne. 0) return
c
c Running averages for instantaneous discharge observations will be 
c   computed for storms, day, month, year and run periods using the
c   number of time steps tot_steps and the vectors of cumulative totals.
c   The discharge vectors contain total values for storm, day, month,
c   year, and run. This approach will need modification if different
c   time steps are premitted in the future (i.e. day to hours during storms)
c
c Add solar radiation as a rate variable also even though it is in Langleys
c   per day.
c
      if(declpri('sumb_tot_steps', 5, 'integer', tot_steps)
     + .ne. 0) return

      if(declpri('sumb_tot_swrad_W', 5, 'real',
     $     tot_swrad_W).ne. 0) return

      if(declpri('sumb_tot_basin_qsim_m3s', 5, 'real',
     $     tot_basin_qsim_m3s).ne. 0) return

      if(declpri('sumb_tot_basin_qobs_m3s', 5, 'real',
     $     tot_basin_qobs_m3s).ne. 0) return

      sumbdecl = 0

      return
      end

c***********************************************************************
c
c     sumbinit - Initialize basinsum module - get parameter values
c                set to zero
c

      integer function sumbinit()

#if defined(__INTEL_COMPILER)
      USE IFPORT
#endif
      USE WEBMOD_SUM
      USE WEBMOD_IO, only: chemout
      USE WEBMOD_OBSCHEM, ONLY :phq_lut,sol_id,sol_name,n_iso,iso_list
      USE WEBMOD_RESMOD, ONLY : basin_sto_cm, basin_chan_sto_cm,
     $  basin_intcp_sto_cm,basin_pweqv_cm,basin_sssto_cm,
     $  basin_soil_moist_cm, basin_gw_sto_cm
      USE WEBMOD_PHREEQ_MMS, ONLY : cheminit_file, c_chem_basin,
     $  c_chem_hyd, chem_sim, nsolute

      IMPLICIT NONE
      !integer, external :: getdataname  ! (defined in obs_chem)

      integer, external :: length  ! (defined in obs_chem)
     
      integer pftemp, nstep
      integer i, lengdat, lengpar, lengpqi
      double precision timestep

      integer DATE_TIME(8)
      character(8) yyyymmdd
      character(10)hhmmssdsss
      integer Exe_yr, Exe_mo, Exe_dy, Exe_hr, Exe_min
      real Exe_sec
      
      integer ret
c      integer today(3), now(3)

      character(90) datafile, paramfile

      character(80) buffer
      character(120) bufferw
      character(4000) bufferd
      character(4000) buffer_txt

 
      sumbinit = 1
      
      if (.not. allocated(bufferclabels)) then
          allocate(character(len=40000) :: bufferclabels)
      endif      
      if (.not. allocated(buffercspan)) then
          allocate(character(len=40000) :: buffercspan)
      endif      
      if (.not. allocated(buffercdat)) then
          allocate(character(len=40000) :: buffercdat)
      endif    
      
      nmru = getdim('nmru')
      if(nmru.eq.-1) return

      if(getparam('topc', 'dtinit', 1 , 'real', dtinit)
     +   .ne.0) return

      timestep = dtinit

      if(getparam('basin', 'basin_area', 1, 'real', basin_area)
     +   .ne.0) return

      if(getparam('sumb', 'print_explanation', 1,
     $     'integer', print_explanation).ne.0) return

      if(getparam('sumb', 'print_type', 1, 'integer', print_type)
     +   .ne.0) return

      if(getparam('sumb', 'print_freq', 1, 'integer', print_freq)
     +   .ne.0) return

      if(getparam('sumb', 'print_objfunc', 1, 'integer', print_objfunc)
     +   .ne.0) return

c Set conversion factor for converting basin volumes in cubic meters
c to depth in centimeters

      m3cm = 1e-4/basin_area

      if(print_freq.gt.0) then

        nstep = getstep()

        if(nstep.eq.0) then


           last_basin_stor = 0.
           last_chan_stor = 0.
c
c set all lastper flags to 1 so that variables are
c zeroed at the beginning of the time step
c

          do 3 i = 1,5
             lastper(i) = 1
 3        end do


C******Set unit print switch      
        pftemp = print_freq 
        if(pftemp.ge.16)  then
          uprt = .true.
          pftemp = pftemp-16
        else
          uprt = .false.
        end if

C******Set daily print switch      
        if(pftemp.ge.8)  then
          dprt = .true.
          pftemp = pftemp-8
        else
          dprt = .false.
        end if

C******Set monthly print switch        
        if(pftemp.ge.4) then
          mprt = .true.
          pftemp = pftemp-4
        else
          mprt = .false.
        end if

C******Set yearly print switch
        if(pftemp.ge.2) then
          yprt = .true.
          pftemp =pftemp-2
        else
          yprt = .false.
        end if

C******Set total print switch
        if(pftemp.eq.1) then
          tprt = .true.
        else
          tprt = .false.
        end if 
C*****End beginning of run inits
      end if

C*****Turn unit switch on for daily values when the user
C     selects daily totals but not unit totals

      if(dprt.and.timestep.ge.24) uprt=.true.


C*****start beginning of restarted run inits
C******Set header print switch
        header_prt = 0
        if(yprt.and.(((uprt.or.dprt).and..not.mprt).or.
     +    ((mprt.and..not.(uprt.or.dprt))))) header_prt = 1
        if((uprt.or.dprt).and.mprt) header_prt = 2
        if(uprt.and.timestep.lt.24.) then
              if(dprt)header_prt = 3
        end if
              
c
c Initial storage volumes are available from WEBMOD_RESMOD variables
c
c
c Print header with time, date, names of parameter and data files and
c if print_explanation is set to 1, print field descriptors
c

c Get the date/time and names of the data file and the parameter file 

      call DATE_AND_TIME(DATE=yyyymmdd, TIME=hhmmssdsss)
      read(yyyymmdd,'(I4,I2,I2)')Exe_YR, Exe_mo, Exe_dy
      read(hhmmssdsss,'(I2,I2,F6.0)')Exe_hr, Exe_min, Exe_sec

c      call idate(today)      ! today(1)=day, (2)=month, (3)=year
c      call itime(now)        ! now(1)=hour, (2)=minute, (3)=second

      ret = getdataname(datafile,'')
      if(ret.eq.-1) then
         print*,'Problems retrieving data file name',
     $           ' Run has been stopped.'
        return
      end if
      lengdat = index(datafile,achar(0))-1
      
      IF (control_string(paramfile, 'param_file')
     &       .NE.0 ) then
        print*,'Problems retrieving parameter file name',
     $           ' Run has been stopped.'
      end if
      lengpar = index(paramfile,achar(0))-1
c
c Print header for chemical output file. Acknowledge no content if chem_sim = 0.
c
      if(chem_sim.eq.1) then
         span1   = '==== == =='
         spansol = '  ================'
         buffercspan = span1
         if(print_type.lt.2) then  ! print properties of dischage at outlet
            write (chemout%lun, 3010)char(13),Exe_mo,Exe_dy,Exe_yr,
     $       Exe_hr,Exe_min, Exe_sec,char(13), datafile(1:lengdat),
     $       char(13),paramfile(1:lengpar),char(13),trim(cheminit_file),
     $       (char(13),i=1,8)
            write(bufferclabels, 5000)(trim(sol_name(i)),i=1,nsolute)
            write(chemout%lun,'(A)')trim(bufferclabels)
            do i=1,nsolute+3
                buffercspan=trim(buffercspan)//spansol
            enddo
            write(chemout%lun,'(A)')trim(buffercspan)
         else ! print mass balance for watershed
             write (chemout%lun, 3020)char(13),Exe_mo,Exe_dy,Exe_yr,
     $       Exe_hr,Exe_min, Exe_sec,char(13), datafile(1:lengdat),
     $       char(13),paramfile(1:lengpar),char(13),trim(cheminit_file),
     $        (char(13),i=1,9)
             write(bufferclabels, 5020)
             do i=1,nsolute
               if(sol_id(i)%iso) then
                 write(bufferclabels, 5021)trim(bufferclabels),
     $            trim(sol_name(i))//'_in',
     $            trim(sol_name(i))//'_out',trim(sol_name(i))//'_ET',
     $            trim(sol_name(i))//'_fin'
               else
                 write(bufferclabels, 5021)trim(bufferclabels),
     $            trim(sol_name(i))//'_in',
     $            trim(sol_name(i))//'_out',trim(sol_name(i))//'_rxn',
     $            trim(sol_name(i))//'_fin'
               endif
             enddo
            write(chemout%lun,'(A)')trim(bufferclabels)
            do i=1,4*nsolute+8
                buffercspan=trim(buffercspan)//spansol
            enddo
            write(chemout%lun,'(A)')trim(buffercspan)
            write(buffercdat,5022)'Init',c_chem_basin%vol(1)
            do i=1,nsolute
                if(sol_id(i)%iso) then
                   write(buffercdat,'(A,56x,F16.1)')trim(buffercdat),
     $               c_chem_basin%delta(i,1)/c_chem_basin%vol(1)
                else
                   write(buffercdat,'(A,56X,F16.1)')trim(buffercdat),
     $               c_chem_basin%M(i,1)
                endif
            enddo           
            write(chemout%lun,'(A)')trim(buffercdat)
         endif
      else
         write (chemout%lun, *)'The parameter chem_sim(one)=0 '//
     $                  ' so no geochemisty has been simulated'
      end if
c
c Print header for hydrology output files. Describe fields if print_explanation = 1.
c
      if(print_explanation.eq.1) then
         write (buffer_txt, 10)char(13),Exe_mo,Exe_dy,Exe_yr,
     $       Exe_hr,Exe_min, Exe_sec,char(13), datafile(1:lengdat),
     $       char(13),paramfile(1:lengpar),(char(13),i=1,71)
         call opstr(trim(buffer_txt))
      else
         write (buffer_txt, 20)char(13),Exe_mo,Exe_dy,Exe_yr,
     $       Exe_hr,Exe_min, Exe_sec,char(13), datafile(1:lengdat),
     $       char(13),paramfile(1:lengpar),(char(13),i=1,3)
         call opstr(trim(buffer_txt))
      end if
c
c Print field labels
c
        if(print_type.eq.0) then
         write(buffer,2000)
         call opstr(trim(buffer))
         write(buffer,2001)
         call opstr(trim(buffer))
         write(buffer,2003)
         call opstr(trim(buffer))
        elseif(print_type.eq.1) then
         write(bufferw,2010)
         call opstr(trim(bufferw))
         write(bufferw,2011)
         call opstr(trim(bufferw))
         write(bufferw,2013)
         call opstr(trim(bufferw))
         write(bufferw,2012) basin_sto_cm
         call opstr(trim(bufferw))
        else
         write(bufferd,2020)
         call opstr(trim(bufferd))
         write(bufferd,2021)
         call opstr(trim(bufferd))
         write(bufferd,2023)
         call opstr(trim(bufferd))
         write(bufferd,2022) basin_sto_cm,basin_intcp_sto_cm,
     $        basin_pweqv_cm,basin_soil_moist_cm,
     $        basin_sssto_cm,basin_gw_sto_cm,basin_chan_sto_cm
         call opstr(trim(bufferd))
        end if

      end if
c
c formats for hydro.out file
c
 10   format('WEBMOD version 1.0',A,
     $ 'Model Run: ', i2.2, '/', i2.2, '/', i4.4, '; ',
     &   i2.2, ':', i2.2, ':', i2.2,A,
     $ 'Data File ',2A ,
     $ 'Parameter File ',3A,
     $ 'Explanation:',A,
     $ 'srad      Incoming shortwave radiation',A,
     $ 'tmax      Maximum mean temperature for basin',A,
     $ 'tmin      Minimum mean temperature for basin',A,
     $ 'Pot_ET    Potential evapotranspiration',A,
     $ 'Act_ET    Actual evapotranspiration',A,
     $ 'Basin_sto Total basin storage',A,
     $ 'Basin_in  Total Basin inputs',A,
     $ 'Basin_out Total Basin outputs',A,
     $ 'Delta_sto Change in total basin storage',A,
     $ 'Precip    Atmospheric precipitation',A,
     $ 'Irr_ext   Irrigation applied from external source(i.e.canal)',A,
     $ 'Irr_sat   Irrigation pumped from shallow wells',A,
     $ 'Irr_hyd   Irrigation pumped or diverted from streams',A,
     $ 'Intcp_sto Storage on canopy and vegetation',A,
     $ 'Intcp     Amount of precipitation intercepted by canopy',A,
     $ 'Can_evap  Canopy evaporation',A,
     $ 'Transp    Transpiration to canopy on dry days',A,
     $ 'Thrufall  Throughfall from canopy',A,
     $ 'Net_dep   Net deposition below canopy',A,
     $ 'Net_rain  Net rain',A,
     $ 'Net_snow  Net snow',A,
     $ 'SWE       Snow-water equivalence (Storage in snowpack)',A,
     $ 'Snowmelt  Snowpack melt reaching ground',A,
     $ 'Snowevap  Evaporation of snowpack',A,
     $ 'Surf_dep  Rain on bare ground plus snow melt',A,
     $ 'OF_Hort   Hortonian overland flow (Infiltration excess)',A,
     $ 'OF_Dunn   Dunnian overland flow (Precip on saturated areas)',A,
     $ 'OF_tot    Total overland flow',A,
     $ 'Infil     Infiltration into soil',A,
     $ 'SoilMoist Soil moisture in root zone',A,
     $ 'UZ_Sto    Storage in unsaturated zone',A,
     $ 'Qdir_flow Lateral macropore flow',A,
     $ 'Soil_evap Evaporation from soil',A,
     $ 'GW_sto    Storage in the saturated zone',A,
     $ 'GW_in1    Leakage from Irrigation canal',A,
     $ 'GW_in2    Upgradient groundwater entering basin',A,
     $ 'Qwet      Wetting of root zone by groundwater',A,
     $ 'Qvpref    Surface deposition delivered directly to ',
     $            'saturated zone',A,
     $ 'UZ2Sat    UZ porewater engulfed by saturated zone as water ',
     $           'table rises',A,
     $ 'Recharge  Matrix recharge to saturated zone from UZ. Does ',
     $           'not include Qvpref or UZ2Sat',A,
     $ 'Sat2UZ    Groundwater stranded in UZ as water table ',
     $           'lowers',A,
     $ 'Qwell     Groundwater pumped to the surface for irrigation',A,
     $ 'Qsatpref  Preferential flow in the saturated zone ',
     $           '(i.e. tile drains)',A,
     $ 'GW_loss   Groundwater loss through aquitard',A,
     $ 'Baseflow  Discharge of groundwater directly to stream',A,
     $ 'Q_Hill    Combine hillslope inputs into stream',A,
     $ 'Chan_sto  Total water stored in channel',A,
     $ 'Diversion Water diverted from stream for irrigation',A,
     $ 'Chan_loss Seepage through streambed',A,
     $ 'Q_sim     Simulated discharge at basin outlet',A,
     $ 'Q_obs     Observed discharge at basin outlet',A,
     $ 'Wat_bal   Error in water balance\n',A,
     $ 'Summary values are total centimeters with the following ',
     $        'exceptions:',A,
     $ '  State variables indicate storage at the beginning of ',
     $        'the run or',A,
     $ '  the end of the time step:',A,
     $ '    Basin_sto, Intcp_sto, SWE,SoilMoist, UZ_sto, GW_sto, ',
     $        'Chan_sto',A,
     $ '  Mean values are reported for solar radiation and ',
     $        'discharge in m3s.',A,
     $ '  Temperature maxima and minima are for the period.',2A,
     $ 'Headline symbols',A,
     $ '== State variables',A,
     $ '++ Basin inputs',A,
     $ '%% Evapotranspiration',A,
     $ '-- Basin outputs',A,
     $ '** Internal basin fluxes',A,
     $ '(This header can be turned off by setting ',
     $        'print_explanation=0)',2A1)

20    format('WEBMOD version 1.0',A,
     $ 'Model Run: ', i2.2, '/', i2.2, '/', i4.4, '; ',
     &   i2.2, ':', i2.2, ':', F4.1,A,
     $ 'Data File ',2A ,
     $ 'Parameter File ',2A,
     $ 'Field descriptions are available if ',
     $ 'print_explanation=1',2A)

 2000 format('   Year Month Day H  M  S    Simulated   Observed')
 2001 format('                              (m^3/s)    (m^3/s) ')
 2003 format(' ------------------------------------------------')


 2010 format(
     $     '    Date  /  Time   Basin_sto  Precip   ',
     $     ' Irr_ext    GW_in     Act_ET   GW_loss  ',
     $     'Chan_loss   Q_sim     Q_obs    Wat_bal  ')
 2011 format(
     $     'Year Mo Dy Hr Mi Se   (cm)      (cm)    ',
     $     '   (cm)      (cm)      (cm)     (cm)    ',
     $     '   (cm)      (cm)      (cm)     (cm)    ')
 2012 format(20x,f9.4)
 2013 format(
     $     '=== == == == == == ========= +++++++++ ',
     $     '+++++++++ +++++++++ %%%%%%%%% --------- ',
     $     '--------- --------- --------- ========= ')
 2020 format('    Date  /  Time   srad tmax tmin  Pot_',
     $ 'ET Act_ET Basin_sto Basin_in  Basin_out ',
     $ 'Delta_sto  Precip    Irr_ext   Irr_sat  ',
     $ ' Irr_hyd  Intcp_sto   Intcp    Int_evap ',
     $ '  Transp   Thrufall  Net_dep   Net_rain ',
     $ ' Net_snow    SWE     Snowmelt  Snowevap ',
     $ ' Surf_dep  OF_Hort   OF_Dunn    OF_tot  ',
     $ '  Infil   SoilMoist   UZ_sto  Qdir_flow ',
     $ 'Soil_evap   GW_sto    GW_in1    GW_in2  ',
     $ '  Qwet      Qvpref    UZ2Sat   Recharge ',
     $ '  Sat2UZ    Qwell    Qsatpref  GW_loss  ',
     $ ' Baseflow   Q_Hill  ',
     $ ' Chan_sto Diversion Chan_loss   Q_sim   ',
     $ '  Q_obs     Q_sim     Q_obs    Wat_bal ')
 2021 format('Year Mo Dy Hr Mi Se W/m2  (C)  (C)   (cm',
     $ ')   (cm)    (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)       m3/s      m3/s      (cm)  ')
 2023 format(
     $ '==== == == == == == ++++ ==== ====  %%%%',
     $ '%% %%%%%% ========= +++++++++ --------- ',
     $ '========= +++++++++ +++++++++ ********* ',
     $ '********* ========= ********* %%%%%%%%% ',
     $ '********* ********* ********* ********* ',
     $ '********* ========= ********* %%%%%%%%% ',
     $ '********* ********* ********* ********* ',
     $ '********* ========= ========= ********* ',
     $ '%%%%%%%%% ========= +++++++++ +++++++++ ',
     $ '********* ********* ********* ********* ',
     $ '********* ********* ********* --------- ',
     $ '********* ********* ',
     $ '========= ********* --------- --------- ',
     $ '--------- --------- --------- =========')


 2022 format(50x,f9.4,71x,f9.4,71x,f9.4,70x,2(1x,f9.4),
     $     21x,f9.4,121x,f9.4)

c
c Formats for chem.out file
c
 3010 format('WEBMOD version 1.0',A,
     $ 'Model Run: ', i2.2, '/', i2.2, '/', i4.4, '; ',
     &   i2.2, ':', i2.2, ':', F4.1,A,
     $ 'Data File ',2A ,
     $ 'Parameter File ',2A,
     $ 'PHREEQC Initialization File ',3A,
     $ 'Stream discharge and chemistry at outlet:',A,
     $ 'Volume      total, in cubic meters',A,
     $ 'TempC       volume-weighted average, in degrees celsius',A,
     $ 'pH          volume-weighted average, in standard units',A,
     $ 'Solute      total, in moles, or ',A,
     $ 'D or [18O]  volume-weighted average, in permil',4A1)

 3020 format('WEBMOD version 1.0',A,
     $ 'Model Run: ', i2.2, '/', i2.2, '/', i4.4, '; ',
     &   i2.2, ':', i2.2, ':', F4.1,A,
     $ 'Data File ',2A ,
     $ 'Parameter File ',2A,
     $ 'PHREEQC Initialization File ',3A,
     $ 'Watershed Mass Balance (_in; _out; _ET, evapotranspiration;',
     $ ' _rxn, (+) added to or (-) removed from solution; _fin,',
     $ ' volume or mass in watershed at end of time step',A,
     $ 'Vol         total volume, in cubic meters',A,
     $ 'Temp        volume-weighted average of stream temperature at ',
     $              'the outlet, in degrees celsius',A,
     $ 'pH          volume-weighted average of stream pH at the ',
     $              'outlet, in standard units',A,
     $ 'Solute      total, in moles, or ',A,
     $ 'D or [18O]  volume-weighted average, in permil. _ET ',
     $              'reflects average delta of evaporated water.',A
     $ 'NA          Not Available. Value assigned to volume-weighted ',
     $             'average on days with no precipitation, ',
     $             'irrigation, or upgradient inputs',4A1)

 5000 format('Year Mo Dy               Vol',
     $ '              Temp                pH', 50(2x,A16))
 5020 format('Year Mo Dy           Temp_in          Temp_out',
     $ '             pH_in            pH_out            Vol_in',
     $ '           Vol_out            Vol_ET           Vol_Fin')
 5021 format(A,4(2x,A16))    
 5022 format(A4,134x,F16.1,100(56x,F16.1))

      sumbinit = 0

      return
      end

c***********************************************************************
c
c     sumbrun - Computes summary values

      integer function sumbrun()

      USE WEBMOD_SUM
      USE WEBMOD_IO, ONLY: nowtime, chemout
      USE WEBMOD_OBSCHEM, ONLY :phq_lut,sol_id,sol_name,n_iso,iso_list
      USE WEBMOD_RESMOD, ONLY : basin_sto_cm, basin_in_cm, basin_out_cm,
     + basin_ppt_cm, basin_irr_ext_cm, basin_irr_sat_cm, 
     + basin_irr_hyd_cm, basin_intcp_sto_cm, basin_intcp_cm, 
     + basin_intcp_evap_cm, basin_thruf_cm, basin_net_dep_cm, 
     + basin_net_rain_cm, basin_net_snow_cm, basin_pweqv_cm, 
     + basin_snowmelt_cm, basin_snowevap_cm, basin_surfdep_cm, 
     + basin_soil_moist_cm, basin_sssto_cm, basin_ofhort_cm, 
     + basin_ofdunn_cm, basin_sroff_cm, basin_infil_cm, basin_qdf_cm, 
     + basin_uz_et_cm, basin_gw_sto_cm, basin_gw1_in_cm, 
     + basin_gw2_in_cm, basin_qwet_cm, basin_vpref_cm, basin_uz2sat_cm,
     + basin_recharge_cm, basin_sat2uz_cm, basin_qwell_cm, 
     + basin_qpref_cm, basin_gwloss_cm, basin_gwflow_cm, 
     + basin_chan_sto_cm, basin_stflow_cm, 
     + basin_chan_div_cm, basin_chan_loss_cm, basin_qsim_cm, 
     + basin_qsim_m3s, basin_qobs_cm, basin_qobs_m3s, basin_et_cm, 
     + basin_transp_cm
      USE WEBMOD_PHREEQ_MMS, ONLY : 
     +  ch_outlet_m3,ch_outlet_tempC,ch_outlet_pH, ch_outlet_M,
     +  ch_outlet_permil, ch_basin_m3_in, ch_basin_m3_out, 
     +  ch_basin_m3_ET, ch_basin_m3_final, ch_basin_tempC_in, 
     +  ch_basin_tempC_out, ch_basin_pH_in, ch_basin_pH_out,
     +  ch_basin_M_in,ch_basin_M_out,ch_basin_M_rxn,ch_basin_M_final,
     $  ch_basin_permil_in, ch_basin_permil_out,ch_basin_permil_ET,
     $  ch_basin_permil_final, chem_sim, nsolute
      
      IMPLICIT NONE

      integer endper

c*****local variables

      logical end_run, end_yr, end_mo, end_dy, end_storm, NoInputs

      character(80) buffer
      character(120) bufferw
      character(560) bufferd
      character(16) :: na='              NA'
      integer year, mo, day
      integer hr, min, sec
      integer nstep
      integer i,j,jj

      real diffop, oflgo, oflgp, diflg


c$$$      double precision timestep,stepcheck,wyday,cday, jday
      double precision timestep,jday

      integer tmax_i,tmin_i,swrad_i
      real basin_tmax_c, basin_tmin_c

      sumbrun = 1

      if(print_freq.gt.0) then

c        call dattim('now', nowtime)
c$$$      wyday = julian('now', 'water')
c$$$      cday = julian('now','calendar')
        jday = djulian('now', 'absolute')
        year = nowtime(1)
        mo = nowtime(2)
        day = nowtime(3)
        hr = nowtime(4)
        min = nowtime(5)
        sec = nowtime(6)
c      jday = jday+rhr/24.+rmin/1440.+rsec/86400.

        nstep = getstep()

        timestep = deltim()


c     Zero accumlator variables of last time step was end of period
c     First step zeros all since lastper was set to all 1s in init

      do 5  i = 1,5
         if(lastper(i).eq.1) then
            tot_steps(i) = 0
            tot_swrad_W(i) = 0.
            basin_tmaxc_max(i) = -999.
            basin_tminc_min(i) = 999.
            basin_potet_cm_sum(i) = 0.
            basin_in_cm_sum(i) = 0.
            basin_out_cm_sum(i) = 0.
            basin_et_cm_sum(i) = 0.
            basin_ppt_cm_sum(i) = 0.
            basin_irr_ext_cm_sum(i) = 0.
            basin_irr_sat_cm_sum(i) = 0.
            basin_irr_hyd_cm_sum(i) = 0.
            basin_intcp_cm_sum(i) = 0.
            basin_intcp_evap_cm_sum(i) = 0.
            basin_transp_cm_sum(i) = 0.
            basin_thruf_cm_sum(i) = 0.
            basin_net_dep_cm_sum(i) = 0.
            basin_net_rain_cm_sum(i) = 0.
            basin_net_snow_cm_sum(i) = 0.
            basin_snowmelt_cm_sum(i) = 0.
            basin_snowevap_cm_sum(i) = 0.
            basin_surfdep_cm_sum(i) = 0.
            basin_ofhort_cm_sum(i) = 0.
            basin_ofdunn_cm_sum(i) = 0.
            basin_sroff_cm_sum(i) = 0.
            basin_infil_cm_sum(i) = 0.
            basin_qdf_cm_sum(i) = 0.
            basin_uz_et_cm_sum(i) = 0.
            basin_gw1_in_cm_sum(i) = 0.
            basin_gw2_in_cm_sum(i) = 0.
            basin_qwet_cm_sum(i) = 0.
            basin_vpref_cm_sum(i) = 0.
            basin_uz2sat_cm_sum(i) = 0.
            basin_recharge_cm_sum(i) = 0.
            basin_sat2uz_cm_sum(i) = 0.
            basin_qwell_cm_sum(i) = 0.
            basin_qpref_cm_sum(i) = 0.
            basin_gwloss_cm_sum(i) = 0.
            basin_gwflow_cm_sum(i) = 0.
            basin_stflow_cm_sum(i) = 0.
            basin_chan_div_cm_sum(i) = 0.
            basin_chan_loss_cm_sum(i) = 0.
            basin_qsim_cm_sum(i) = 0.
            basin_qobs_cm_sum(i) = 0.
            tot_basin_qsim_m3s(i) = 0.
            tot_basin_qobs_m3s(i) = 0.
            wat_bal_cm_sum(i) = 0.
            do 3 j= 1,5
               obj_func_sum(j,i) = 0
 3          continue
c
c Accumulators for geochemistry
c
            ch_outlet_m3_sum(i) = 0
            ch_outlet_tempC_avg(i) = 0
            ch_outlet_pH_avg(i) = 0
            ch_basin_m3_in_sum(i) = 0
            ch_basin_m3_out_sum(i) = 0
            ch_basin_m3_ET_sum(i) = 0
            ch_basin_tempC_in_avg(i) = 0
            ch_basin_tempC_out_avg(i) = 0
            ch_basin_pH_in_avg(i) = 0
            ch_basin_pH_out_avg(i) = 0
            do 4 j=1,nsolute
              ch_outlet_permil_avg(j,i) = 0
              ch_outlet_M_sum(j,i) = 0
              ch_basin_permil_in_avg(j,i) = 0
              ch_basin_permil_out_avg(j,i) = 0
              ch_basin_permil_ET_avg(j,i) = 0
              ch_basin_M_in_sum(j,i) = 0
              ch_basin_M_out_sum(j,i) = 0
              ch_basin_M_rxn_sum(j,i) = 0
 4          continue            
            lastper(i) = 0
         end if
 5    continue
      NoInputs=.false.  ! test to assign NA value to average Temp, pH, and deltas when no input
c
c Count time steps for computing average discharge values for radiation
c   and discharge. Volume weights will be used for temperatures and pH
c   entering and leaving the basin
c
      do 8 i = 1,5
         tot_steps(i) = tot_steps(i) + 1
 8    continue


      if(getvar('io', 'endper', 1, 'integer', endper)
     +   .ne.0) return

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
c Also keep track of days in the month and the year so that average values
c are computed correctly for runs beginning or ending in partial months
c or years.
c
      if (endper.ne.0.and.mod(endper,2).ne.0) then
         endper = endper - 1
         end_storm = .true.
c         lastper(1) = 1
      end if
      
      if (endper.ne.0) then
         end_dy = .true.
c Until a storm summary section is implemented just let storm
c totals (index 1) equal daily totals
         lastper(1) = 1
         lastper(2) = 1
         endper = endper - 2
         if (endper.ne.0) then
            end_mo = .true.
            lastper(3) = 1
            endper = endper - 4
            if (endper.ne.0) then
               end_yr = .true.
               lastper(4) = 1
               endper = endper - 8
               if (endper.ne.0) end_run = .true.
            end if
         end if
      end if

c
c radiation and temperature retrieved directly
c

      if(getvar('solrad', 'orad', 1, 'real',
     +   swrad_W).ne.0) return

      if(getvar('temp', 'basin_tmax_c', 1, 'real', basin_tmax_c)
     +   .ne.0) return

      if(getvar('temp', 'basin_tmin_c', 1, 'real', basin_tmin_c)
     +   .ne.0) return

c
c     Convert radiation in Langleys per day to average Watts per
c     square meter
c
      swrad_W = swrad_W * Ld2Wm2
c
c Get the basin potet and convert to cm
c
      if(getvar('potet', 'basin_potet', 1, 'real', basin_potet_cm)
     +   .ne.0) return

      basin_potet_cm = basin_potet_cm * 2.54
c
c Get unit basin variables from webmod_res through USE WEBMOD_RES
c

c*****Compute error in water balance

      wat_bal_cm = basin_ppt_cm + basin_irr_ext_cm - basin_et_cm -
     $     basin_stflow_cm + basin_gw1_in_cm + basin_gw2_in_cm -
     $     basin_gwloss_cm - basin_chan_loss_cm +
     $     last_basin_stor - basin_sto_cm +
     $     basin_chan_sto_cm - last_chan_stor

      if(nstep.eq.1) wat_bal_cm = 0

      last_basin_stor = basin_sto_cm
      last_chan_stor = basin_chan_sto_cm

c$$$      last_basin_stor = basin_storage_cm

c******Compute Objective Function

      diffop = basin_qobs_m3s - basin_qsim_m3s
      obj_func(1) = abs(diffop)
      obj_func(2) = diffop*diffop
      oflgo = alog(basin_qobs_m3s + 1.)
      oflgp = alog(basin_qsim_m3s + 1.)
      diflg = oflgo-oflgp
      obj_func(3) = abs(diflg)
      obj_func(4) = diflg*diflg
      obj_func(5) = diffop

c****** Accumulate unit values into sums for storm, day,
c       month, year, and run totals

      do 130 i = 1,5
c  Hydrology
       tot_swrad_W(i) = tot_swrad_W(i) + swrad_W
       swrad_W_avg(i) = tot_swrad_W(i) / tot_steps(i)
       if(basin_tmax_c.gt.basin_tmaxc_max(i))
     +    basin_tmaxc_max(i) = basin_tmax_c
       if(basin_tmin_c.lt.basin_tminc_min(i))
     +    basin_tminc_min(i) = basin_tmin_c
       basin_potet_cm_sum(i) = basin_potet_cm_sum(i) + basin_potet_cm
       basin_et_cm_sum(i) = basin_et_cm_sum(i) + basin_et_cm
       basin_in_cm_sum(i) = basin_in_cm_sum(i) + basin_in_cm
       basin_out_cm_sum(i) = basin_out_cm_sum(i) + basin_out_cm
       basin_ppt_cm_sum(i) = basin_ppt_cm_sum(i) + basin_ppt_cm
       basin_irr_ext_cm_sum(i) = basin_irr_ext_cm_sum(i) +
     $     basin_irr_ext_cm
       basin_irr_sat_cm_sum(i) = basin_irr_sat_cm_sum(i) +
     $     basin_irr_sat_cm
       basin_irr_hyd_cm_sum(i) = basin_irr_hyd_cm_sum(i) +
     $     basin_irr_hyd_cm
       basin_intcp_cm_sum(i) = basin_intcp_cm_sum(i) +
     $     basin_intcp_cm
       basin_intcp_evap_cm_sum(i) = basin_intcp_evap_cm_sum(i) +
     $     basin_intcp_evap_cm
       basin_transp_cm_sum(i) = basin_transp_cm_sum(i) + basin_transp_cm
       basin_thruf_cm_sum(i) = basin_thruf_cm_sum(i) + basin_thruf_cm
       basin_net_dep_cm_sum(i) = basin_net_dep_cm_sum(i) +
     $     basin_net_dep_cm
       basin_net_rain_cm_sum(i) = basin_net_rain_cm_sum(i) +
     $     basin_net_rain_cm
       basin_net_snow_cm_sum(i) = basin_net_snow_cm_sum(i) +
     $     basin_net_snow_cm
       basin_snowmelt_cm_sum(i) = basin_snowmelt_cm_sum(i) +
     $     basin_snowmelt_cm
      basin_snowevap_cm_sum(i) = basin_snowevap_cm_sum(i) +
     $     basin_snowevap_cm
      basin_surfdep_cm_sum(i) = basin_surfdep_cm_sum(i) +
     $     basin_surfdep_cm
      basin_ofhort_cm_sum(i) = basin_ofhort_cm_sum(i) + basin_ofhort_cm
      basin_ofdunn_cm_sum(i) = basin_ofdunn_cm_sum(i) + basin_ofdunn_cm
      basin_sroff_cm_sum(i) = basin_sroff_cm_sum(i) + basin_sroff_cm
      basin_infil_cm_sum(i) = basin_infil_cm_sum(i) + basin_infil_cm
      basin_qdf_cm_sum(i) = basin_qdf_cm_sum(i) + basin_qdf_cm
      basin_uz_et_cm_sum(i) = basin_uz_et_cm_sum(i) + basin_uz_et_cm
      basin_gw1_in_cm_sum(i) = basin_gw1_in_cm_sum(i) + basin_gw1_in_cm
      basin_gw2_in_cm_sum(i) = basin_gw2_in_cm_sum(i) + basin_gw2_in_cm
      basin_qwet_cm_sum(i) = basin_qwet_cm_sum(i) + basin_qwet_cm
      basin_vpref_cm_sum(i) = basin_vpref_cm_sum(i) + basin_vpref_cm
      basin_uz2sat_cm_sum(i) = basin_uz2sat_cm_sum(i) + basin_uz2sat_cm
      basin_recharge_cm_sum(i) = basin_recharge_cm_sum(i) +
     $     basin_recharge_cm
      basin_sat2uz_cm_sum(i) = basin_sat2uz_cm_sum(i) + basin_sat2uz_cm
      basin_qwell_cm_sum(i) = basin_qwell_cm_sum(i) + basin_qwell_cm
      basin_qpref_cm_sum(i) = basin_qpref_cm_sum(i) + basin_qpref_cm
      basin_gwloss_cm_sum(i) = basin_gwloss_cm_sum(i) + basin_gwloss_cm
      basin_gwflow_cm_sum(i) = basin_gwflow_cm_sum(i) + basin_gwflow_cm
      basin_stflow_cm_sum(i) = basin_stflow_cm_sum(i) + basin_stflow_cm
      basin_chan_div_cm_sum(i) = basin_chan_div_cm_sum(i) +
     $     basin_chan_div_cm
      basin_chan_loss_cm_sum(i) = basin_chan_loss_cm_sum(i) +
     $     basin_chan_loss_cm
      basin_qsim_cm_sum(i) = basin_qsim_cm_sum(i) + basin_qsim_cm
      basin_qobs_cm_sum(i) = basin_qobs_cm_sum(i) + basin_qobs_cm
      tot_basin_qsim_m3s(i) = tot_basin_qsim_m3s(i) + basin_qsim_m3s
      tot_basin_qobs_m3s(i) = tot_basin_qobs_m3s(i) +
     $     basin_qobs_m3s
      basin_qsim_m3s_avg(i) = tot_basin_qsim_m3s(i) / tot_steps(i)
      basin_qobs_m3s_avg(i) = tot_basin_qobs_m3s(i) / tot_steps(i)
      wat_bal_cm_sum(i) = wat_bal_cm_sum(i) + wat_bal_cm
         do 125 jj = 1,5
            obj_func_sum(jj,i) = obj_func_sum(jj,i) + obj_func(jj)
 125  end do
c  Geochemisty
      if(chem_sim.eq.1) then
         ch_outlet_m3_sum(i) = ch_outlet_m3_sum(i)+ch_outlet_m3 
         ch_outlet_tempC_avg(i) = ch_outlet_tempC_avg(i)
     +     +ch_outlet_tempC*ch_outlet_m3
         ch_outlet_pH_avg(i) = ch_outlet_pH_avg(i)
     +     +ch_outlet_pH*ch_outlet_m3
         ch_basin_m3_in_sum(i)= ch_basin_m3_in_sum(i)+ ch_basin_m3_in
         ch_basin_m3_out_sum(i)= ch_basin_m3_out_sum(i)+ ch_basin_m3_out
         ch_basin_m3_ET_sum(i)= ch_basin_m3_ET_sum(i)+ ch_basin_m3_ET
         ch_basin_tempC_in_avg(i) = ch_basin_tempC_in_avg(i)
     +     +ch_basin_tempC_in*ch_basin_m3_in
         ch_basin_tempC_out_avg(i) = ch_basin_tempC_out_avg(i)
     +     +ch_basin_tempC_out*ch_basin_m3_out
         ch_basin_pH_in_avg(i) = ch_basin_pH_in_avg(i)
     +     +ch_basin_pH_in*ch_basin_m3_in
         ch_basin_pH_out_avg(i) = ch_basin_pH_out_avg(i)
     +     +ch_basin_pH_out*ch_basin_m3_out
         do 129 jj=1,nsolute
          ch_outlet_permil_avg(jj,i) = ch_outlet_permil_avg(jj,i)+
     +           ch_outlet_permil(jj)*ch_outlet_m3
          ch_outlet_M_sum(jj,i) = ch_outlet_M_sum(jj,i) +
     +           ch_outlet_M(jj)
          ch_basin_permil_in_avg(jj,i) = ch_basin_permil_in_avg(jj,i)+
     +           ch_basin_permil_in(jj)*ch_basin_m3_in
          ch_basin_permil_out_avg(jj,i) = ch_basin_permil_out_avg(jj,i)+
     +           ch_basin_permil_out(jj)*ch_basin_m3_out
          ch_basin_permil_ET_avg(jj,i) = ch_basin_permil_ET_avg(jj,i)+
     +           ch_basin_permil_ET(jj)*ch_basin_m3_ET
          ch_basin_M_in_sum(jj,i) = ch_basin_M_in_sum(jj,i) +
     +           ch_basin_M_in(jj)
          ch_basin_M_out_sum(jj,i) = ch_basin_M_out_sum(jj,i)+ 
     +           ch_basin_M_out(jj)
          ch_basin_M_rxn_sum(jj,i) = ch_basin_M_rxn_sum(jj,i)+ 
     +           ch_basin_M_rxn(jj)
 129    end do    
       end if ! chem_sim.eq.1
 130  end do

     
c
c
c                   write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
c     $               ch_basin_M_in_sum(i,4), ch_basin_M_out_sum(i,4),
c     $               ch_basin_M_rxn_sum(i,4),ch_basin_M_final(i)
c       buffercdat='Heres the starting field '
c      do i=1,nsolute
c                write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
c     $              ch_basin_M_in_sum(i,4), ch_basin_M_out_sum(i,3),
c     $               ch_basin_M_rxn_sum(i,3),ch_basin_M_final(i)
c
c      end do
c      write(buffercdat,'(A)')trim(buffercdat)
c      write(buffercdat,'(A,F16.1)')trim(buffercdat),
c     $    ch_basin_M_out_sum(i,4)
c      write(chemout%lun,'(A)')trim(buffercdat)
c      pause
c                   write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
c     $               ch_basin_M_in_sum(i,4), ch_basin_M_out_sum(i,4),
c     $               ch_basin_M_rxn_sum(i,4),ch_basin_M_final(i)
c
c

c Begin unit value reporting

      if(uprt) then
         if(print_type.eq.0) then
            write(buffer,2002)year, mo, day, hr, min, sec,
     +           basin_qsim_m3s,basin_qobs_m3s
            call opstr(trim(buffer))
         end if

        if(print_type.eq.1) then
           write(bufferw,2012)year, mo, day, hr, min, sec,
     +          basin_sto_cm, basin_ppt_cm, basin_irr_ext_cm,
     $          basin_gw1_in_cm, basin_et_cm,basin_gwloss_cm,
     $          basin_chan_loss_cm,basin_qsim_cm,basin_qobs_cm,
     $          wat_bal_cm
           call opstr(trim(bufferw))
       end if
       
        if(print_type.eq.2) then
           swrad_i = swrad_W
           tmax_i = basin_tmax_c
           tmin_i = basin_tmin_c
              write(bufferd,2022)year,mo,day,hr,min,sec,swrad_i,
     +          tmax_i, tmin_i,
     $          basin_potet_cm, basin_et_cm,  basin_sto_cm,
     $          basin_in_cm,basin_out_cm,
     $          basin_in_cm-basin_out_cm-basin_et_cm,
     $          basin_ppt_cm,basin_irr_ext_cm,
     $          basin_irr_sat_cm,basin_irr_hyd_cm,
     $          basin_intcp_sto_cm,basin_intcp_cm,
     $          basin_intcp_evap_cm,basin_transp_cm,
     $          basin_thruf_cm,basin_net_dep_cm,
     $          basin_net_rain_cm,basin_net_snow_cm,
     $          basin_pweqv_cm,basin_snowmelt_cm,
     $          basin_snowevap_cm,
     $          basin_surfdep_cm,basin_ofhort_cm,
     $          basin_ofdunn_cm,basin_sroff_cm,
     $          basin_infil_cm,basin_soil_moist_cm,
     $          basin_sssto_cm,basin_qdf_cm,
     $          basin_uz_et_cm, basin_gw_sto_cm,
     $          basin_gw1_in_cm,basin_gw2_in_cm,
     $          basin_qwet_cm,basin_vpref_cm,
     $          basin_uz2sat_cm,
     $          basin_recharge_cm,basin_sat2uz_cm,
     $          basin_qwell_cm,basin_qpref_cm,
     $          basin_gwloss_cm,basin_gwflow_cm,
     $          basin_stflow_cm,basin_chan_sto_cm,
     $          basin_chan_div_cm,basin_chan_loss_cm,
     $          basin_qsim_cm,basin_qobs_cm,
     $          basin_qsim_m3s,basin_qobs_m3s,
     $          wat_bal_cm
          call opstr(trim(bufferd))
        end if
c Geochem          
        if(chem_sim.eq.1) then
          if(print_type.lt.2) then ! print outlet chem
            write(buffercdat,5012)year,mo,day,ch_outlet_m3,
     +        ch_outlet_tempC,ch_outlet_pH
              do i=1,nsolute
                if(sol_id(i)%iso) then
                   write(buffercdat,'(A,2x,F16.1)')trim(buffercdat),
     $               ch_outlet_permil(i)
                else
                   write(buffercdat,'(A,56X,F16.1)')trim(buffercdat),
     $               ch_outlet_M(i)
                endif
              enddo
            write(chemout%lun,'(A)')trim(buffercdat)
          else ! print basin mass balance
            if(ch_basin_m3_in.eq.0) NoInputs=.TRUE.
            if(NoInputs) then ! NA values for average input temp, pH, and deltas
              write(buffercdat,5020)year,mo,day,NA,
     +          ch_basin_tempC_out,NA, ch_basin_pH_out,
     +          ch_basin_m3_in, ch_basin_m3_out, ch_basin_m3_ET,
     +          ch_basin_m3_final
            else
              write(buffercdat,5022)year,mo,day,ch_basin_tempC_in,
     +          ch_basin_tempC_out,ch_basin_pH_in, ch_basin_pH_out,
     +          ch_basin_m3_in, ch_basin_m3_out, ch_basin_m3_ET,
     +          ch_basin_m3_final
            endif
            do i=1,nsolute
              if(sol_id(i)%iso) then
                if(NoInputs) then
                  write(buffercdat,'(A,2x,A16,3(2x,F16.1))')
     $               trim(buffercdat),NA, ch_basin_permil_out(i),
     $               ch_basin_permil_ET(i),ch_basin_permil_final(i)
                else
                  write(buffercdat,'(A,4(2x,F16.1))')
     $               trim(buffercdat),ch_basin_permil_in(i),
     $                ch_basin_permil_out(i),ch_basin_permil_ET(i),
     $                ch_basin_permil_final(i)
                endif
              else
               write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
     $                   ch_basin_M_in(i), ch_basin_M_out(i),
     $                   ch_basin_M_rxn(i),ch_basin_M_final(i)
              endif
            enddo
            write(chemout%lun,'(A)')trim(buffercdat)
          end if
        end if ! chem_sim
      end if
      NoInputs=.FALSE.

c
c Storm summary section (index 1 of 5) should be included here once
c implemented
c

c Begin daily summaries

      if(end_dy.and.timestep.lt.24.) then
         if(dprt) then
               if(print_type.eq.0) then
                  if(dprt.eqv..true.) then
                     write(buffer,903)
                     call opstr(trim(buffer))
                  end if
                  write(buffer,904)year,mo,day,basin_qsim_m3s_avg(2),
     $                 basin_qobs_m3s_avg(2)
                  call opstr(trim(buffer))
               end if

               if(print_type.eq.1) then
                  if(dprt.eqv..true.) then
                     write(bufferw,913)
                     call opstr(trim(bufferw))
                  end if
                  write(bufferw,914)year, mo, day,
     +          basin_sto_cm, basin_ppt_cm_sum(2),
     $          basin_irr_ext_cm_sum(2),basin_gw1_in_cm_sum(2),
     $          basin_et_cm_sum(2),basin_gwloss_cm_sum(2),
     $          basin_chan_loss_cm_sum(2),basin_qsim_cm_sum(2),
     $          basin_qobs_cm_sum(2),wat_bal_cm_sum(2)
                  call opstr(trim(bufferw))
               end if

               if(print_type.eq.2) then
                  if(dprt.eqv..true.) then
                     write(bufferd,923)
                     call opstr(trim(bufferd))
                  end if

                  swrad_i = swrad_W_avg(2)
                  tmax_i = basin_tmaxc_max(2)
                  tmin_i = basin_tminc_min(2)

                  write(bufferd,924)year,mo,day,
     $       swrad_i,tmax_i,tmin_i,
     $       basin_potet_cm_sum(2), basin_et_cm_sum(2),
     $       basin_sto_cm,
     $       basin_in_cm_sum(2),basin_out_cm_sum(2),
     $       basin_in_cm_sum(2)-basin_out_cm_sum(2)-basin_et_cm_sum(2),
     $       basin_ppt_cm_sum(2),basin_irr_ext_cm_sum(2),
     $       basin_irr_sat_cm_sum(2),basin_irr_hyd_cm_sum(2),
     $       basin_intcp_sto_cm,basin_intcp_cm_sum(2),
     $       basin_intcp_evap_cm_sum(2),basin_transp_cm_sum(2),
     $       basin_thruf_cm_sum(2),basin_net_dep_cm_sum(2),
     $       basin_net_rain_cm_sum(2),basin_net_snow_cm_sum(2),
     $       basin_pweqv_cm,basin_snowmelt_cm_sum(2),
     $       basin_snowevap_cm_sum(2),
     $       basin_surfdep_cm_sum(2),basin_ofhort_cm_sum(2),
     $       basin_ofdunn_cm_sum(2),basin_sroff_cm_sum(2),
     $       basin_infil_cm_sum(2),basin_soil_moist_cm,
     $       basin_sssto_cm,basin_qdf_cm_sum(2),
     $       basin_uz_et_cm_sum(2), basin_gw_sto_cm,
     $       basin_gw1_in_cm_sum(2),basin_gw2_in_cm_sum(2),
     $       basin_qwet_cm_sum(2),basin_vpref_cm_sum(2),
     $       basin_uz2sat_cm_sum(2),
     $       basin_recharge_cm_sum(2),basin_sat2uz_cm_sum(2),
     $       basin_qwell_cm_sum(2),basin_qpref_cm_sum(2),
     $       basin_gwloss_cm_sum(2),basin_gwflow_cm_sum(2),
     $       basin_stflow_cm_sum(2),basin_chan_sto_cm,
     $       basin_chan_div_cm_sum(2),basin_chan_loss_cm_sum(2),
     $       basin_qsim_cm_sum(2),basin_qobs_cm_sum(2),
     $       basin_qsim_m3s_avg(2),basin_qobs_m3s_avg(2),
     $       wat_bal_cm_sum(2)

                  call opstr(trim(bufferd))
               end if
            if(print_objfunc.eq.1) then
               write(bufferw,1032)
               call opstr(trim(bufferw))
               write(bufferd,1033) (obj_func_sum(jj,2), jj=1,5)
               call opstr(trim(bufferd))
            end if
         end if
      end if
      NoInputs=.FALSE.

c End Day Section

c Begin monthly summaries

      if(end_mo) then
         if(mprt) then
c
c Hydrology
c
            if(print_type.eq.0) then
               if(dprt.eqv..true.) then
                  write(buffer,2003)
                  call opstr(trim(buffer))
               end if
               write(buffer,2004)year, mo, basin_qsim_m3s_avg(3),
     $              basin_qobs_m3s_avg(3)
               call opstr(trim(buffer))
            end if

            if(print_type.eq.1) then
               if(dprt.eqv..true.) then
                  write(bufferw,2013)
                  call opstr(trim(bufferw))
               end if
               write(bufferw,2014)year, mo,
     +          basin_sto_cm, basin_ppt_cm_sum(3),
     $          basin_irr_ext_cm_sum(3),basin_gw1_in_cm_sum(3),
     $          basin_et_cm_sum(3),basin_gwloss_cm_sum(3),
     $          basin_chan_loss_cm_sum(3),basin_qsim_cm_sum(3),
     $          basin_qobs_cm_sum(3),wat_bal_cm_sum(3)
               call opstr(trim(bufferw))
            end if

            if(print_type.eq.2) then
               if(dprt.eqv..true.) then
                  write(bufferd,2023)
                  call opstr(trim(bufferd))
               end if
                  swrad_i = swrad_W_avg(3)
                  tmax_i = basin_tmaxc_max(3)
                  tmin_i = basin_tminc_min(3)
               write(bufferd,2024)year,mo,
     $       swrad_i,tmax_i,tmin_i, 
     $       basin_potet_cm_sum(3), basin_et_cm_sum(3),
     $       basin_sto_cm,
     $       basin_in_cm_sum(3),basin_out_cm_sum(3),
     $       basin_in_cm_sum(3)-basin_out_cm_sum(3)-basin_et_cm_sum(3),
     $       basin_ppt_cm_sum(3),basin_irr_ext_cm_sum(3),
     $       basin_irr_sat_cm_sum(3),basin_irr_hyd_cm_sum(3),
     $       basin_intcp_sto_cm,basin_intcp_cm_sum(3),
     $       basin_intcp_evap_cm_sum(3),basin_transp_cm_sum(3),
     $       basin_thruf_cm_sum(3),basin_net_dep_cm_sum(3),
     $       basin_net_rain_cm_sum(3),basin_net_snow_cm_sum(3),
     $       basin_pweqv_cm,basin_snowmelt_cm_sum(3),
     $       basin_snowevap_cm_sum(3),
     $       basin_surfdep_cm_sum(3),basin_ofhort_cm_sum(3),
     $       basin_ofdunn_cm_sum(3),basin_sroff_cm_sum(3),
     $       basin_infil_cm_sum(3),basin_soil_moist_cm,
     $       basin_sssto_cm,basin_qdf_cm_sum(3),
     $       basin_uz_et_cm_sum(3), basin_gw_sto_cm,
     $       basin_gw1_in_cm_sum(3),basin_gw2_in_cm_sum(3),
     $       basin_qwet_cm_sum(3),basin_vpref_cm_sum(3),
     $       basin_uz2sat_cm_sum(3),
     $       basin_recharge_cm_sum(3),basin_sat2uz_cm_sum(3),
     $       basin_qwell_cm_sum(3),basin_qpref_cm_sum(3),
     $       basin_gwloss_cm_sum(3),basin_gwflow_cm_sum(3),
     $       basin_stflow_cm_sum(3),basin_chan_sto_cm,
     $       basin_chan_div_cm_sum(3),basin_chan_loss_cm_sum(3),
     $       basin_qsim_cm_sum(3),basin_qobs_cm_sum(3),
     $       basin_qsim_m3s_avg(3),basin_qobs_m3s_avg(3),
     $       wat_bal_cm_sum(3)
               call opstr(trim(bufferd))
            endif   

            if(print_objfunc.eq.1) then
               write(bufferw,1034)
               call opstr(trim(bufferw))
               write(bufferd,1035) (obj_func_sum(jj,3), jj=1,5)
               call opstr(trim(bufferd))
            end if
c
c Geochem
c
          if(chem_sim.eq.1) then
          write(chemout%lun,'(A)')trim(buffercspan)
              if(print_type.lt.2) then ! print outlet chem
                write(buffercdat,5014)year,mo,ch_outlet_m3_sum(3),
     +            ch_outlet_tempC_avg(3)/ch_outlet_m3_sum(3),
     +            ch_outlet_pH_avg(3)/ch_outlet_m3_sum(3)
                  do i=1,nsolute
                    if(sol_id(i)%iso) then
                     write(buffercdat,'(A,2x,F16.1)')trim(buffercdat),
     $                  ch_outlet_permil_avg(i,3)/ch_outlet_m3_sum(3)
                    else
                     write(buffercdat,'(A,2X,F16.1)')trim(buffercdat),
     $                  ch_outlet_M_sum(i,3)
                    endif
                  enddo
                write(chemout%lun,'(A)')trim(buffercdat)
              else ! print basin mass balance
                if(ch_basin_m3_in_sum(3).eq.0) NoInputs=.TRUE.
                if(NoInputs) then ! NA values for average input temp, pH, and deltas
                  write(buffercdat,5024)year,mo,NA,
     +             ch_basin_tempC_out_avg(3)/ch_basin_m3_out_sum(3),NA,
     +             ch_basin_pH_out_avg(3)/ch_basin_m3_out_sum(3),
     +             ch_basin_m3_in_sum(3), ch_basin_m3_out_sum(3),
     +             ch_basin_m3_ET_sum(3),ch_basin_m3_final
                else
                  write(buffercdat,5026)year,mo,
     +             ch_basin_tempC_in_avg(3)/ch_basin_m3_in_sum(3),
     +             ch_basin_tempC_out_avg(3)/ch_basin_m3_out_sum(3),
     +             ch_basin_pH_in_avg(3)/ch_basin_m3_in_sum(3),
     +             ch_basin_pH_out_avg(3)/ch_basin_m3_out_sum(3),
     +             ch_basin_m3_in_sum(3), ch_basin_m3_out_sum(3),
     +             ch_basin_m3_ET_sum(3),ch_basin_m3_final
                endif
                do i=1,nsolute
                  if(sol_id(i)%iso) then
                   if(NoInputs) then
                    write(buffercdat,'(A,2x,A16,3(2x,F16.1))')
     $               trim(buffercdat),NA, ch_basin_permil_out_avg(i,3)/
     $               ch_basin_m3_out_sum(3),ch_basin_permil_ET_avg(i,3)/
     $               ch_basin_m3_ET_sum(3),ch_basin_permil_final(i)
                   else
                    write(buffercdat,'(A,4(2x,F16.1))')
     $              trim(buffercdat),ch_basin_permil_in_avg(i,3)/
     $              ch_basin_m3_in_sum(3),ch_basin_permil_out_avg(i,3)/
     $              ch_basin_m3_out_sum(3),ch_basin_permil_ET_avg(i,3)/
     $              ch_basin_m3_ET_sum(3),ch_basin_permil_final(i)
                   endif
                  else
                   write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
     $               ch_basin_M_in_sum(i,3), ch_basin_M_out_sum(i,3),
     $               ch_basin_M_rxn_sum(i,3),ch_basin_M_final(i)
                  endif
                enddo
                write(chemout%lun,'(A)')trim(buffercdat)
              end if
              if(end_yr.and.yprt.or.end_run.and.tprt) then
                write(chemout%lun,'(A)')trim(buffercspan)
              else
                write(chemout%lun,'(2A)')char(13),
     $           trim(bufferclabels)
                write(chemout%lun,'(A)')trim(buffercspan)
              endif
            end if  ! chem_sim
           end if
      end if
      NoInputs=.FALSE.

c****** Begin yearly summaries

      if(yprt) then
        if(end_yr) then
          if(print_type.eq.0) then
            if(mprt.eqv..true..or.dprt.eqv..true.) then
              write(buffer,2005)
              call opstr(trim(buffer))
            end if
            write(buffer,2006)year,  basin_qsim_m3s_avg(4),
     $          basin_qobs_m3s_avg(4)
            call opstr(trim(buffer))
          end if

          if(print_type.eq.1) then
            if(mprt.eqv..true..or.dprt.eqv..true.) then
              write(bufferw,2013)
              call opstr(trim(bufferw))
            end if
            write(bufferw,2016)year,
     +          basin_sto_cm, basin_ppt_cm_sum(4),
     $          basin_irr_ext_cm_sum(4),basin_gw1_in_cm_sum(4),
     $          basin_et_cm_sum(4),basin_gwloss_cm_sum(4),
     $          basin_chan_loss_cm_sum(4),basin_qsim_cm_sum(4),
     $          basin_qobs_cm_sum(4),wat_bal_cm_sum(4)
            call opstr(trim(bufferw))
         end if

          if(print_type.eq.2) then
            if(mprt.eqv..true..or.dprt.eqv..true.) then
              write(bufferd,2023)
              call opstr(trim(bufferd))
           end if

                  swrad_i = swrad_W_avg(4)
                  tmax_i = basin_tmaxc_max(4)
                  tmin_i = basin_tminc_min(4)

            write(bufferd,2026)year,
     $       swrad_i,tmax_i,tmin_i,
     $       basin_potet_cm_sum(4), basin_et_cm_sum(4),
     $       basin_sto_cm,
     $       basin_in_cm_sum(4),basin_out_cm_sum(4),
     $       basin_in_cm_sum(4)-basin_out_cm_sum(4)-basin_et_cm_sum(4),
     $       basin_ppt_cm_sum(4),basin_irr_ext_cm_sum(4),
     $       basin_irr_sat_cm_sum(4),basin_irr_hyd_cm_sum(4),
     $       basin_intcp_sto_cm,basin_intcp_cm_sum(4),
     $       basin_intcp_evap_cm_sum(4),basin_transp_cm_sum(4),
     $       basin_thruf_cm_sum(4),basin_net_dep_cm_sum(4),
     $       basin_net_rain_cm_sum(4),basin_net_snow_cm_sum(4),
     $       basin_pweqv_cm,basin_snowmelt_cm_sum(4),
     $       basin_snowevap_cm_sum(4),
     $       basin_surfdep_cm_sum(4),basin_ofhort_cm_sum(4),
     $       basin_ofdunn_cm_sum(4),basin_sroff_cm_sum(4),
     $       basin_infil_cm_sum(4),basin_soil_moist_cm,
     $       basin_sssto_cm,basin_qdf_cm_sum(4),
     $       basin_uz_et_cm_sum(4), basin_gw_sto_cm,
     $       basin_gw1_in_cm_sum(4),basin_gw2_in_cm_sum(4),
     $       basin_qwet_cm_sum(4),basin_vpref_cm_sum(4),
     $       basin_uz2sat_cm_sum(4),
     $       basin_recharge_cm_sum(4),basin_sat2uz_cm_sum(4),
     $       basin_qwell_cm_sum(4),basin_qpref_cm_sum(4),
     $       basin_gwloss_cm_sum(4),basin_gwflow_cm_sum(4),
     $       basin_stflow_cm_sum(4),basin_chan_sto_cm,
     $       basin_chan_div_cm_sum(4),basin_chan_loss_cm_sum(4),
     $       basin_qsim_cm_sum(4),basin_qobs_cm_sum(4),
     $       basin_qsim_m3s_avg(4),basin_qobs_m3s_avg(4),
     $       wat_bal_cm_sum(4)
            call opstr(trim(bufferd))
         end if

           if(print_objfunc.eq.1) then
             write(bufferw,1036)
             call opstr(trim(bufferw))
             write(bufferd,1037) (obj_func_sum(i,4), i=1,5)
             call opstr(trim(bufferd))
           end if
c
c Geochem
c
           if(chem_sim.eq.1) then
c             write(chemout%lun,'(A)')trim(buffercspan)
               if(print_type.lt.2) then ! print outlet chem
                 write(buffercdat,5016)year,ch_outlet_m3_sum(4),
     +            ch_outlet_tempC_avg(4)/ch_outlet_m3_sum(4),
     +            ch_outlet_pH_avg(4)/ch_outlet_m3_sum(4)
                  do i=1,nsolute
                   if(sol_id(i)%iso) then
                     write(buffercdat,'(A,2x,F16.1)')trim(buffercdat),
     $                  ch_outlet_permil_avg(i,4)/ch_outlet_m3_sum(4)
                    else
                     write(buffercdat,'(A,2X,F16.1)')trim(buffercdat),
     $                  ch_outlet_M_sum(i,4)
                    endif
                  enddo
                write(chemout%lun,'(A)')trim(buffercdat)
              else ! print basin mass balance
                if(ch_basin_m3_in_sum(4).eq.0) NoInputs=.TRUE.
                if(NoInputs) then ! NA values for average input temp, pH, and deltas
                  write(buffercdat,5028)year,NA,
     +             ch_basin_tempC_out_avg(4)/ch_basin_m3_out_sum(4),NA,
     +             ch_basin_pH_out_avg(4)/ch_basin_m3_out_sum(4),
     +             ch_basin_m3_in_sum(4), ch_basin_m3_out_sum(4),
     +             ch_basin_m3_ET_sum(4),ch_basin_m3_final
                else
                  write(buffercdat,5030)year,
     +             ch_basin_tempC_in_avg(4)/ch_basin_m3_in_sum(4),
     +             ch_basin_tempC_out_avg(4)/ch_basin_m3_out_sum(4),
     +             ch_basin_pH_in_avg(4)/ch_basin_m3_in_sum(4),
     +             ch_basin_pH_out_avg(4)/ch_basin_m3_out_sum(4),
     +             ch_basin_m3_in_sum(4), ch_basin_m3_out_sum(4),
     +             ch_basin_m3_ET_sum(4),ch_basin_m3_final
                endif
                do i=1,nsolute
                 if(sol_id(i)%iso) then
                  if(NoInputs) then
                   write(buffercdat,'(A,2x,A16,3(2x,F16.1))')
     $              trim(buffercdat),NA, ch_basin_permil_out_avg(i,4)/
     $              ch_basin_m3_out_sum(4),ch_basin_permil_ET_avg(i,4)/
     $              ch_basin_m3_ET_sum(4),ch_basin_permil_final(i)
                  else
                   write(buffercdat,'(A,4(2x,F16.1))')
     $             trim(buffercdat),ch_basin_permil_in_avg(i,4)/
     $             ch_basin_m3_in_sum(4),ch_basin_permil_out_avg(i,4)/
     $             ch_basin_m3_out_sum(4),ch_basin_permil_ET_avg(i,4)/
     $             ch_basin_m3_ET_sum(4),ch_basin_permil_final(i)
                  endif
                else
                   write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
     $               ch_basin_M_in_sum(i,4), ch_basin_M_out_sum(i,4),
     $               ch_basin_M_rxn_sum(i,4),ch_basin_M_final(i)
                endif
               enddo
               write(chemout%lun,'(A)')trim(buffercdat)
              end if
              if(end_run.and.tprt) then
                write(chemout%lun,'(A)')trim(buffercspan)
              else
                write(chemout%lun,'(2A)')char(13),
     $           trim(bufferclabels)
                write(chemout%lun,'(A)')trim(buffercspan)
              endif
            end if  ! chem_sim
        end if
      end if

c******Print heading if needed

      if(end_run.eqv..false.) then
         if((header_prt.eq.2.and.end_mo).or.
     $        (header_prt.eq.1.and.end_yr).or.
     $        (header_prt.eq.3.and.end_dy))
     +        then
            write(buffer,999)
            call opstr(trim(buffer))
            if(print_type.eq.0) then
               write(buffer,2000)
               call opstr(trim(buffer))
               write(buffer,2001)
               call opstr(trim(buffer))
               write(buffer,2003)
               call opstr(trim(buffer))
            else if(print_type.eq.1) then
               write(bufferw,2010)
               call opstr(trim(bufferw))
               write(bufferw,2011)
               call opstr(trim(bufferw))
               write(bufferw,2013)
               call opstr(trim(bufferw))
            else
               write(bufferd,2020)
               call opstr(trim(bufferd))
               write(bufferd,2021)
               call opstr(trim(bufferd))
               write(bufferd,2023)
               call opstr(trim(bufferd))
            end if
         end if
      end if
      NoInputs=.FALSE.
c****** Print run summary
c
      if(end_run) then
         if(tprt) then
c
c Hydrology
c
            if(print_type.eq.0) then
               write(buffer,2007)
               call opstr(trim(buffer))
               write(buffer,2008) basin_qsim_m3s_avg(5),
     $              basin_qobs_m3s_avg(5)
               call opstr(trim(buffer))
            end if

            if(print_type.eq.1) then
               write(bufferw,2013)
               call opstr(trim(bufferw))
               write(bufferw,2018)
     +          basin_sto_cm, basin_ppt_cm_sum(5),
     $          basin_irr_ext_cm_sum(5),basin_gw1_in_cm_sum(5),
     $          basin_et_cm_sum(5),basin_gwloss_cm_sum(5),
     $          basin_chan_loss_cm_sum(5),basin_qsim_cm_sum(5),
     $          basin_qobs_cm_sum(5),wat_bal_cm_sum(5)
               call opstr(trim(bufferw))
            end if

            if(print_type.eq.2) then
               write(bufferd,2023)
               call opstr(trim(bufferd))

               swrad_i = swrad_W_avg(5)
               tmax_i = basin_tmaxc_max(5)
               tmin_i = basin_tminc_min(5)

               write(bufferd,2028)swrad_i,tmax_i,tmin_i,
     $       basin_potet_cm_sum(5), basin_et_cm_sum(5),
     $       basin_sto_cm,
     $       basin_in_cm_sum(5),basin_out_cm_sum(5),
     $       basin_in_cm_sum(5)-basin_out_cm_sum(5)-basin_et_cm_sum(5),
     $       basin_ppt_cm_sum(5),basin_irr_ext_cm_sum(5),
     $       basin_irr_sat_cm_sum(5),basin_irr_hyd_cm_sum(5),
     $       basin_intcp_sto_cm,basin_intcp_cm_sum(5),
     $       basin_intcp_evap_cm_sum(5),basin_transp_cm_sum(5),
     $       basin_thruf_cm_sum(5),basin_net_dep_cm_sum(5),
     $       basin_net_rain_cm_sum(5),basin_net_snow_cm_sum(5),
     $       basin_pweqv_cm,basin_snowmelt_cm_sum(5),
     $       basin_snowevap_cm_sum(5),
     $       basin_surfdep_cm_sum(5),basin_ofhort_cm_sum(5),
     $       basin_ofdunn_cm_sum(5),basin_sroff_cm_sum(5),
     $       basin_infil_cm_sum(5),basin_soil_moist_cm,
     $       basin_sssto_cm,basin_qdf_cm_sum(5),
     $       basin_uz_et_cm_sum(5), basin_gw_sto_cm,
     $       basin_gw1_in_cm_sum(5),basin_gw2_in_cm_sum(5),
     $       basin_qwet_cm_sum(5),basin_vpref_cm_sum(5),
     $       basin_uz2sat_cm_sum(5),
     $       basin_recharge_cm_sum(5),basin_sat2uz_cm_sum(5),
     $       basin_qwell_cm_sum(5),basin_qpref_cm_sum(5),
     $       basin_gwloss_cm_sum(5),basin_gwflow_cm_sum(5),
     $       basin_stflow_cm_sum(5),basin_chan_sto_cm,
     $       basin_chan_div_cm_sum(5),basin_chan_loss_cm_sum(5),
     $       basin_qsim_cm_sum(5),basin_qobs_cm_sum(5),
     $       basin_qsim_m3s_avg(5),basin_qobs_m3s_avg(5),
     $       wat_bal_cm_sum(5)
               call opstr(trim(bufferd))
            end if

            if(print_objfunc.eq.1) then
               write(bufferw,1038)
               call opstr(trim(bufferw))
               write(bufferd,1039) (obj_func_sum(jj,5), jj=1,5)
               call opstr(trim(bufferd))
            end if  
c
c Geochem
c
           if(chem_sim.eq.1) then
c             write(chemout%lun,'(A)')trim(buffercspan)
               if(print_type.lt.2) then ! print outlet chem
                 write(buffercdat,5018)ch_outlet_m3_sum(5),
     +            ch_outlet_tempC_avg(5)/ch_outlet_m3_sum(5),
     +            ch_outlet_pH_avg(5)/ch_outlet_m3_sum(5)
                  do i=1,nsolute
                   if(sol_id(i)%iso) then
                     write(buffercdat,'(A,2x,F16.1)')trim(buffercdat),
     $                  ch_outlet_permil_avg(i,5)/ch_outlet_m3_sum(5)
                    else
                     write(buffercdat,'(A,2X,F16.1)')trim(buffercdat),
     $                  ch_outlet_M_sum(i,5)
                    endif
                  enddo
                write(chemout%lun,'(A)')trim(buffercdat)
              else ! print basin mass balance
                if(ch_basin_m3_in_sum(5).eq.0) NoInputs=.TRUE.
                if(NoInputs) then ! NA values for average input temp, pH, and deltas
                  write(buffercdat,5032)NA,
     +             ch_basin_tempC_out_avg(5)/ch_basin_m3_out_sum(5),NA,
     +             ch_basin_pH_out_avg(5)/ch_basin_m3_out_sum(5),
     +             ch_basin_m3_in_sum(5), ch_basin_m3_out_sum(5),
     +             ch_basin_m3_ET_sum(5),ch_basin_m3_final
                else
                  write(buffercdat,5034)
     +             ch_basin_tempC_in_avg(5)/ch_basin_m3_in_sum(5),
     +             ch_basin_tempC_out_avg(5)/ch_basin_m3_out_sum(5),
     +             ch_basin_pH_in_avg(5)/ch_basin_m3_in_sum(5),
     +             ch_basin_pH_out_avg(5)/ch_basin_m3_out_sum(5),
     +             ch_basin_m3_in_sum(5), ch_basin_m3_out_sum(5),
     +             ch_basin_m3_ET_sum(5),ch_basin_m3_final
                endif
                do i=1,nsolute
                 if(sol_id(i)%iso) then
                  if(NoInputs) then
                   write(buffercdat,'(A,2x,A16,3(2x,F16.1))')
     $              trim(buffercdat),NA, ch_basin_permil_out_avg(i,5)/
     $              ch_basin_m3_out_sum(5),ch_basin_permil_ET_avg(i,5)/
     $              ch_basin_m3_ET_sum(5),ch_basin_permil_final(i)
                  else
                   write(buffercdat,'(A,4(2x,F16.1))')
     $             trim(buffercdat),ch_basin_permil_in_avg(i,5)/
     $             ch_basin_m3_in_sum(5),ch_basin_permil_out_avg(i,5)/
     $             ch_basin_m3_out_sum(5),ch_basin_permil_ET_avg(i,5)/
     $             ch_basin_m3_ET_sum(5),ch_basin_permil_final(i)
                  endif
                else
                   write(buffercdat,'(A,4(2x,F16.1))')trim(buffercdat),
     $               ch_basin_M_in_sum(i,5), ch_basin_M_out_sum(i,5),
     $               ch_basin_M_rxn_sum(i,5),ch_basin_M_final(i)
                endif
               enddo
               write(chemout%lun,'(A)')trim(buffercdat)
              end if
            end if  ! chem_sim
         end if
      end if
      end if
c      end if

  999 format('                                        ')
 1000 format('   Year Month Day   Predicted  Observed ')
 1001 format('                     (m^3/s)    (m^3/s) ')
 1002 format(2x,f5.0,2x,f3.0,2x,f3.0,f11.4,1x,f11.4)
 1003 format(' ---------------------------------------')
 1004 format(2x,f5.0,2x,f3.0,5x,f11.4,1x,f11.4)
 1005 format(' =======================================')
 1006 format(2x,f5.0,10x,f11.4,1x,f11.4)
 1007 format(' ***************************************')
 1008 format(' Total for run',3x,2f12.4)

 1010 format('   Year Month Day  Precip     ET    Storage  P-Runoff' 
     +,' O-Runoff  Wat_Bal ')
 1011 format('                    (cm)     (cm)     (cm)     (cm)   '
     +,'  (cm)     (cm)  ')
 1012 format(2x,f5.0,2x,f3.0,2x,f3.0,6(1x,f8.3))
 1013 format(' -----------------------------------------------------'
     +,'-----------------')
 1014 format(2x,f5.0,2x,f3.0,5x,6(1x,f8.3))
 1015 format(' ====================================================='
     +,'=================')
 1016 format(2x,f5.0,10x,6(1x,f8.3))
 1017 format(' *****************************************************'
     +,'*****************')
 1018 format(' Total for run',3x,6(1x,f8.3))

 1032 format('Daily Objective Functions ',
     +'(Observed-Simulated Discharge)')
 1033 format(' Abs Dif= ',e9.2,' Dif Sq= ',e9.2,
     +' Abs Diflg= ',e9.2,' Diflg Sq= ',e9.2,
     +' Day Sum Unit Difs = ', e9.2)
 1034 format('Monthly Objective Functions ',
     +'(Observed-Simulated Discharge)')
 1035 format(' Abs Dif= ',e9.2,' Dif Sq= ',e9.2,
     +' Abs Diflg= ',e9.2,' Diflg Sq= ',e9.2,
     +' Month Sum Unit Difs = ', e9.2)
 1036 format('Yearly Objective Functions ',
     +'(Observed-Simulated Discharge)')
 1037 format(' Abs Dif= ',e9.2,' Dif Sq= ',e9.2,
     +' Abs Diflg= ',e9.2,' Diflg Sq= ',e9.2,
     +' Year Sum Unit Difs = ', e9.2)
 1038 format('Total Run Objective Functions ',
     +'(Observed-Simulated Discharge)')
 1039 format(' Abs Dif= ',e9.2,' Dif Sq= ',e9.2,
     +' Abs Diflg= ',e9.2,' Diflg Sq= ',e9.2,
     +' Total Run Sum Unit Difs = ', e9.2)
C  1040 format('Pred Vol= ', f9.2, '  Pred Peak= ', f9.2)
C New formats for unit value runs - RW

 903  format(' ------------------------------------------------')
 904  format(2x,f5.0,2x,f3.0,2x,f3.0,9x,f11.4,1x,f11.4)
 913  format(' -----------------------------------------------------'
     +,'--------------------------')
 914  format(i4.4,1x,2(i2.2,1x),9x,10(f9.4,1x))
 923  format(' --------------------------------------------------',
     +'--------------------------------------------------------------',
     +'-----------------------------------------------------------')
 2000 format('   Year Month Day H  M  S   Simulated   Observed')
 2001 format('                              (m^3/s)    (m^3/s) ')
 2002 format(2x,f5.0,2x,f3.0,2x,4f3.0,f11.4,1x,f11.4)
 2003 format(' ------------------------------------------------')
 2004 format(2x,f5.0,2x,f3.0,14x,f11.4,1x,f11.4)
 2005 format(' ================================================')
 2006 format(2x,f5.0,19x,f11.4,1x,f11.4)
 2007 format(' ************************************************')
 2008 format(' Total for run',12x,2f12.2)

 2010 format(
     $     '    Date  /  Time   Basin_sto  Precip   ',
     $     ' Irr_ext    GW_in     Act_ET   GW_loss  ',
     $     'Chan_loss   Q_sim     Q_obs    Wat_bal  ')
 2011 format(
     $     'Year Mo Dy Hr Mi Se   (cm)      (cm)    ',
     $     '   (cm)      (cm)      (cm)     (cm)    ',
     $     '   (cm)      (cm)      (cm)     (cm)    ')
 2012 format(i4.4,1x,5(i2.2,1x),10(f9.4,1x))
 2013 format(
     $     '==== == == == == == ========= +++++++++ ',
     $     '+++++++++ +++++++++ %%%%%%%%% --------- ',
     $     '--------- --------- --------- ========= ')

 2014 format(i4.4,1x,i2.2,1x,12x,10(f9.4,1x))
 2016 format(i4.4,16x,10(f9.4,1x))
 2018 format(' Total for run',6x,10(f9.2,1x))

 2020 format(
     $ '    Date  /  Time   srad tmax tmin  Pot_',
     $ 'ET Act_ET Basin_sto Basin_in  Basin_out ',
     $ 'Delta_sto  Precip    Irr_ext   Irr_sat  ',
     $ ' Irr_hyd  Intcp_sto   Intcp    Int_evap ',
     $ '  Transp   Thrufall  Net_dep   Net_rain ',
     $ ' Net_snow    SWE     Snowmelt  Snowevap ',
     $ ' Surf_dep  OF_Hort   OF_Dunn    OF_tot  ',
     $ '  Infil   SoilMoist   UZ_sto  Qdir_flow ',
     $ 'Soil_evap   GW_sto    GW_in1    GW_in2  ',
     $ '  Qwet      Qvpref    UZ2Sat   Recharge ',
     $ '  Sat2UZ    Qwell    Qsatpref  GW_loss  ',
     $ ' Baseflow   Q_Hill  ',
     $ ' Chan_sto Diversion Chan_loss   Q_sim   ',
     $ '  Q_obs     Q_sim     Q_obs    Wat_bal ')
2021  format(
     $ 'Year Mo Dy Hr Mi Se W/m2  (C)  (C)   (cm',
     $ ')   (cm)    (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)      (cm)    ',
     $ '  (cm)      (cm)      (cm)      (cm)    ',
     $ '  (cm)       m3/s      m3/s      (cm)  ')
 2023 format(
     $ '==== == == == == == ++++ ==== ====  %%%%',
     $ '%% %%%%%% ========= +++++++++ --------- ',
     $ '========= +++++++++ +++++++++ ********* ',
     $ '********* ========= ********* %%%%%%%%% ',
     $ '********* ********* ********* ********* ',
     $ '********* ========= ********* %%%%%%%%% ',
     $ '********* ********* ********* ********* ',
     $ '********* ========= ========= ********* ',
     $ '%%%%%%%%% ========= +++++++++ +++++++++ ',
     $ '********* ********* ********* ********* ',
     $ '********* ********* ********* --------- ',
     $ '********* ********* ',
     $ '========= ********* --------- --------- ',
     $ '--------- --------- --------- =========')
 2022 format(i4.4,1x,5(i2.2,1x),3(i4.3,1x),2f7.2,49(1x,f9.4))
 924  format(i4.4,1x,2(i2.2,1x),9x,3(i4.3,1x),2f7.2,49(1x,f9.4))
 2024 format(i4.4,1x,i2.2,13x,3(i4.3,1x),2f7.2,49(1x,f9.4))
 2026 format(i4.4,16x,3(i4.3,1x),2f7.2,49(1x,f9.4))
 2028 format(' Total for Run ',5x,3(i4.3,1x),2f7.1,49(1x,f9.2))
c
c Formats for geochem output
c print_type = 0 or 1
c
c print_type = 0 or 1
c
 5012 format(i4.4,2(1x,i2.2),3(2x,F16.1))  ! Year month day Vol tempc pH
 5014 format(i4.4,1x,i2.2,3x,3(2x,F16.1))  ! Year month     Vol tempc pH
 5016 format(i4.4,6x,3(2x,F16.1))          ! Year           Vol tempc pH
 5018 format('Run Totals',3(2x,F16.1))     ! Run Totals
c
c print_type = 2
c
 5020 format(i4.4,2(1x,i2.2),2(2x,A16,2x,F16.1),4(2x,F16.1))  ! Year month day 2tempcs 2pHs 4Vols, NA for input temp and pH
 5022 format(i4.4,1x,i2.2,1x,i2.2,8(2x,F16.1))                ! Year month day 2tempcs 2pHs 4Vols
 5024 format(i4.4,1x,i2.2,3x,2(2x,A16,2x,F16.1),4(2x,F16.1))  ! Year month     2tempcs 2pHs 4Vols, NA for input temp and pH
 5026 format(i4.4,1x,i2.2,3x,8(2x,F16.1))                     ! Year month     2tempcs 2pHs 4Vols
 5028 format(i4.4,6x,2(2x,A16,2x,F16.1),4(2x,F16.1))          ! Year           2tempcs 2pHs 4Vols, NA for input temp and pH
 5030 format(i4.4,6x,8(2x,F16.1))                             ! Year           2tempcs 2pHs 4Vols
 5032 format('Run Totals',2(2x,A16,2x,F16.1),4(2x,F16.1))     ! Run Totals     2tempcs 2pHs 4Vols, NA for input temp and pH
 5034 format('Run Totals',8(2x,F16.1))                        ! Run Totals     2tempcs 2pHs 4Vols

      sumbrun = 0
      if(end_run) print*, "Finished"
      

      return
      end


