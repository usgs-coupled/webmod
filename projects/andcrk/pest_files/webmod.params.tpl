ptf %
Andrews Creek
Version: 1.7
** Dimensions **
####
nchan
1
####
nresinp
21
####
ndays
366
####
ngw_ext
1
####
nrain
2
####
ntemp
1
####
one
1
####
nac
11
####
nchem_sets
10
####
nirrig_int
1
####
nexlag
2
####
nac_nmru_nresinp
231
####
nchemobs
1
####
ndeplval
22
####
nobs
3
####
nmru
1
####
nchemvar
10
####
nhydro
4
####
nlapse
3
####
nmonths
12
####
nconvert
3
####
five
5
####
nform
0
####
nmru_res
9
Canopy
Snowpack
Imperm
Ohoriz
UZ_Wet
UZ_Dry
UZPref
Sat
SatPref
####
nsolute
11
Ca
@Calcium
Mg
@Magnesium
Na
@Sodium
K
@Potassium
Amm
@Ammonia
Alkalinity
@ Alkalinity
Cl
@Chloride
S
@Sulfate
N(5)
@Nitrate
Si
@Silica
[18O]
@ Oxygen Isotope
####
nsnow
0
####
nirrig_ext
1
####
nhum
1
####
nxkbin
9
####
ntopchan
5
####
ndepl
2
####
nsol
1
####
nphq_lut
75
####
nevap
0
####
nchem_ext
1
** Parameters **
####
trxn_ohoriz_stat 9
1
nmru
1
1
1
####
trxn_uz_stat 9
1
nmru
1
1
1
####
trxn_sat_stat 9
1
nmru
1
1
1
####
trxn_ohoriz_days 9
1
one
1
2
% todays       %
####
trxn_uz_days 9
1
one
1
2
% tudays       %
####
trxn_sat_days 9
1
one
1
2
% tsdays       %
####
trxn_ohoriz_c_adj 9
1
nmru
1
2
% tocadj       %
####
trxn_uz_c_adj 9
1
nmru
1
2
% tucadj       %
####
trxn_sat_c_adj 9
1
nmru
1
2
% tscadj       %
####
radj_wppt 14
1
one
1
2
% radjwp       %
####
adjmix_rain 12
1
nmonths
12
2
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
% adjmix       %
####
transp_beg 11
1
nmru
1
1
5
####
SUBRATE 8
1
nmru
1
2
0.01
####
convfactor 11
1
nconvert
3
3
1.0
1.0
1.0
####
rain_adj 9
2
nmru
nmonths
12
2
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
% rainad       %
####
UADJ 8
1
nmru
1
2
% windad       %
####
snow_adj 9
2
nmru
nmonths
12
2
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
% sn_adj       %
####
init_surf_mru 14
1
nmru
1
1
1
####
pmac_sat 9
1
nmru
1
2
% pmacst       %
####
init_eq_ph_mru 15
1
nmru
1
1
1
####
irrig_int_src 14
1
nmru
1
1
1
####
nacsc 8
1
nmru
1
1
11
####
ALAT 5
1
one
1
2
40.28738
####
iout 5
1
one
1
1
2
####
sbar0 6
1
nmru
1
2
0.0
####
qdffrac 8
1
nmru
1
2
% qdffrac      %
####
c_mru 6
1
nchemvar
10
1
1
1
1
1
1
1
1
1
1
1
####
q0 8
1
one
1
2
0.06
####
covden_win 11
1
nmru
1
2
% cdenwn       %
####
init_soln_ppt 14
1
one
1
1
0
####
atmos_eq_ph 14
1
one
1
1
0
####
s_porosity 11
1
nmru
1
2
% s_thpo       %
####
exchset_table 14
2
nmru_res
nchem_sets
90
1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
####
mru2chan 12
1
nmru
1
1
1
####
solnset_table 14
2
nmru_res
nchem_sets
90
1
2
2
2
2
3
3
3
4
4
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
####
xdebug_start 4
1
one
1
1
0
####
tmax_allsnow_c 14
1
one
1
2
% tmxals       %
####
pmacro 7
1
nmru
1
2
% pmacro       %
####
s_rock_depth 13
1
nmru
1
2
% s_rock       %
####
srain_intcp 12
1
nmru
1
2
% sraini       %
####
init_exch_mru 14
1
nmru
1
1
1
####
sched_gw1 10
1
nmru
1
1
0
####
sched_gw2 10
1
nmru
1
1
0
####
hf 8
1
nmru
1
2
0.01
####
epan_coef 10
1
nmonths
12
2
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
% pancof       %
####
ccov_intcp 11
1
nmonths
12
2
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
% ccvint       %
####
strain_adj 11
2
nmru
nmonths
12
2
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
####
c_can_depth 12
1
one
1
2
% candep       %
####
init_eq_ph_hydro 17
1
nhydro
4
1
2
2
2
2
####
gwbnd_len2 11
1
nmru
1
2
0.0
####
xk_cv 8
1
nmru
1
2
% xkcvar       %
####
temp_units 14
1
one
1
1
1
####
mru_aspect 11
1
nmru
1
2
60.0
####
gwbnd_len1 11
1
nmru
1
2
0.0
####
ccov_slope 11
1
nmonths
12
2
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
% ccvslp       %
####
mru_area_frac 14
1
nmru
1
2
1.0
####
TIPM 8
1
nmru
1
2
% sn_flx       %
####
NEGHSI 8
1
nmru
1
2
0.0
####
WEI 8
1
nmru
1
2
15.0
####
irrig_int_init 15
1
nmru
1
2
0.0
####
kinset_table 13
2
nmru_res
nchem_sets
90
1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
####
snow_intcp 11
1
nmru
1
2
% sn_int       %
####
dtinit 7
1
one
1
2
24.0
####
MFMIN 8
1
nmru
1
2
% meltmn       %
####
chan_loss_rate 4
1
one
1
2
% chloss       %
####
mru_tsta 9
1
nmru
1
1
1
####
print_objfunc 14
1
one
1
1
0
####
C_units 8
1
nchemvar
10
1
13
13
13
13
7
7
7
7
7
7
####
src_ext_irrig 14
1
nmru
1
1
1
####
transp_tmax_c 12
1
nmru
1
2
% trdegc       %
####
irrig_sched_int 16
1
nmru
1
1
0
####
hamon_coef 11
1
nmonths
12
2
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
% hamcof       %
####
mru_slope 10
1
nmru
1
2
0.8
####
d 12
2
ntopchan
nchan
5
2
0.0
1000.0
2000.0
3000.0
4000.0
####
crad_coef 10
1
one
1
2
% crdcof       %
####
chem_ext 9
1
one
1
1
1
####
MRUDEPL 8
1
nmru
1
1
2
####
szm 8
1
nmru
1
2
% recess       %
####
nchan_d 8
1
nchan
1
1
5
####
rxnset_table 13
2
nmru_res
nchem_sets
90
1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
####
tsta_elev 10
1
ntemp
1
2
3505.0
####
snowmelt_18O_depl 14
1
nmru
1
2
% sO_dpl       %
####
covden_sum 11
1
nmru
1
2
% cdensm       %
####
tmax_lapse 11
1
nmonths
12
2
10.4
8.2
8.9
10.1
9.3
10.0
9.8
8.5
8.6
8.0
9.3
10.6
####
crad_exp 9
1
one
1
2
% crddexp      %
####
DAYGM 8
1
nmru
1
2
% daymlt       %
####
init_rxn_hydro 15
1
nhydro
4
1
-1
-1
-1
-1
####
mru_lat 8
1
nmru
1
2
40.29
####
C_ires 7
1
nchemvar
10
1
2
2
4
5
1
2
4
5
8
99
####
wrain_intcp 12
1
nmru
1
2
% wraini       %
####
init_soln_ext 14
1
nchem_ext
1
1
1
####
qobsta 7
1
one
1
1
1
####
surfset_table 14
2
nmru_res
nchem_sets
90
1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
####
C_metric 9
1
nchemvar
10
1
3
5
5
5
5
3
5
5
5
5
####
mru_percent_imperv 19
1
nmru
1
2
0.0
####
xdebug_stop 4
1
one
1
1
0
####
tmin_lapse 11
1
nmonths
12
2
5.0
10.2
8.4
7.2
6.8
5.2
5.4
4.9
6.3
7.6
8.6
2.2
####
td 8
1
nmru
1
2
% td           %
####
ppt_chem 9
1
one
1
1
1
####
init_kin_hydro 15
1
nhydro
4
1
-1
-1
-1
-1
####
print_freq 11
1
one
1
1
9
####
print_explanation 4
1
one
1
1
0
####
s_satpref_zmin 15
1
nmru
1
2
% s_zmin       %
####
transp_end 11
1
nmru
1
1
9
####
s_ohoriz_depth 15
1
nmru
1
2
% s_ohor       %
####
mru_elev 9
1
nmru
1
2
3510.0
####
dth 8
1
nmru
1
2
1.0
####
st 8
2
nac
nmru
11
2
8050.947754
98.354874
55.918774
41.412315
32.021423
25.06893
19.79916
15.07565
10.469971
4.462064
1.463311
####
T_decay
1
nmru
1
1
1
####
infex 8
1
one
1
1
1
####
init_exch_hydro 16
1
nhydro
4
1
-1
-1
-1
-1
####
eq_phset_table 15
2
nmru_res
nchem_sets
90
1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
####
chv 8
1
one
1
2
% chvelo       %
####
s_root_depth 13
1
nmru
1
2
% s_root       %
####
ACUMX 8
1
nmru
1
2
0.0
####
iso_n 8
1
nmru
1
3
0.5
####
iso_theta 8
1
nmru
1
3
1.0
####
iso_fac 8
1
nmru
1
3
% isofac       %
####
snow_ion_factor 14
1
nmru
1
2
% sn_ion       %
####
s_theta_fc 11
1
nmru
1
2
% s_thfc       %
####
s_theta_0 11
1
nmru
1
2
% s_thfc       %
####
init_soln_hydro 16
1
nhydro
4
1
3
3
3
3
####
s_satpref_zmax 15
1
nmru
1
2
% s_zmax       %
####
MFMAX 8
1
nmru
1
2
% meltmx       %
####
mru_area 9
1
nmru
1
2
1.74
####
cov_type 9
1
nmru
1
1
3
####
s_theta_wp 11
1
nmru
1
2
% s_thwp       %
####
riparian_thresh 6
1
nmru
1
2
60,0
####
basin_area 14
1
one
1
2
1.74
####
ADC 4
1
ndeplval
22
2
0.05
0.24
0.4
0.52
0.65
0.75
0.82
0.88
0.93
0.99
1.0
0.02
0.047
0.114
0.212
0.318
0.416
0.498
0.56
0.592
0.604
0.604
####
c_rip 6
1
nchemvar
10
1
0
0
0
0
0
0
0
0
0
0
####
To 8
1
nmru
1
2
% transm       %
####
tl 8
1
nmru
1
2
71.5705387
####
c_stindx 9
1
nchemvar
10
1
0
0
0
5
0
0
0
5
0
0
####
irrig_int_max 14
1
nmru
1
2
100.0
####
c_obs_indx 11
1
nchemvar
10
1
0
0
0
0
0
0
0
0
0
0
####
init_rxn_mru 13
1
nmru
1
1
1
####
mru_psta 9
1
nmru
1
1
1
####
gw_loss_k 4
1
nmru
1
2
% gwloss       %
####
NMF 8
1
nmru
1
2
% nmeltf       %
####
chem_sim 9
1
one
1
1
1
####
tmin_adj 9
1
nmru
1
2
% tmnadj       %
####
snowmelt_D_depl 14
1
nmru
1
2
% sD_dpl       %
####
radj_sppt 14
1
one
1
2
% radjsp       %
####
print_type 11
1
one
1
1
2
####
print_vse 11
1
one
1
1
0
####
irrig_sched_ext 16
1
nmru
1
1
0
####
PLWHC 8
1
nmru
1
2
% pklwhc       %
####
c_hyd_indx 11
1
nchemvar
10
1
0
0
0
0
0
0
0
0
0
0
####
MBASE 8
1
nmru
1
2
% meltba       %
####
init_surf_hydro 16
1
nhydro
4
1
-1
-1
-1
-1
####
ppt_rad_adj 4
1
nmonths
12
2
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
% pptrad       %
####
xk0 8
1
nmru
1
2
% xkvert       %
####
SI 8
1
nmru
1
2
% sn_thr       %
####
ach 12
2
ntopchan
nchan
5
2
0.0
0.274
0.742
0.982
1.0
####
init_soln_mru 14
1
nmru
1
1
1
####
src_gw2 7
1
nmru
1
1
1
####
src_gw1 7
1
nmru
1
1
1
####
radmax 14
1
one
1
2
% radmax       %
####
tmax_adj 9
1
nmru
1
2
% tmxadj       %
####
ac 8
2
nac
nmru
11
2
0.0
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
####
init_kin_mru 13
1
nmru
1
1
1
####
tmax_allrain_c 13
1
nmonths
12
2
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
% tmxalr       %
####
s_satpref_k 12
1
nmru
1
2
% s_satk       %
