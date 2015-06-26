ptf %
Rio Icacos
Version: 1.7
** Dimensions **
####
nrain
1
####
nresinp
21
####
nxkbin
9
####
ntemp
1
####
nchemvar
10
####
nirrig_ext
1
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
nphq_lut
75
####
nchem_sets
10
####
nac
11
####
nsol
0
####
ndepl
2
####
nchem_ext
1
####
nmonths
12
####
five
5
####
nchemobs
1
####
nobs
1
####
nhum
1
####
nac_nmru_nresinp
231
####
nhydro
1
####
one
1
####
ngw_ext
1
####
nlapse
3
####
nexlag
2
####
ntopchan
5
####
ndeplval
22
####
nmru
1
####
nevap
0
####
nsnow
0
####
ndays
366
####
nirrig_int
1
####
nchan
1
####
nconvert
3
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
** Parameters **
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
tmax_adj 9
1
nmru
1
2
% tmxadj       %
####
irrig_int_src 14
1
nmru
1
1
1
####
xdebug_start 4
1
one
1
1
2568
####
C_metric 9
1
nchemvar
10
1
5
5
5
5
5
5
5
5
5
5
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
pmacro 7
1
nmru
1
2
% pmacro       %
####
snow_ion_factor 14
1
nmru
1
2
% sn_ion       %
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
TIPM 8
1
nmru
1
2
% sn_flx       %
####
transp_beg 11
1
nmru
1
1
1
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
30
30
30
30
4
4
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
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
s_satpref_k 12
1
nmru
1
2
% s_satk       %
####
mru_aspect 11
1
nmru
1
2
150.0
####
ppt_chem 9
1
one
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
trxn_ohoriz_stat 9
1
nmru
1
1
1
####
radj_sppt 14
1
one
1
2
% radjsp       %
####
init_kin_mru 13
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
init_surf_mru 14
1
nmru
1
1
1
####
NEGHSI 8
1
nmru
1
2
0.0
####
transp_tmax_c 12
1
nmru
1
2
% trdegc       %
####
covden_sum 11
1
nmru
1
2
% cdensm       %
####
tmin_adj 9
1
nmru
1
2
% tmnadj       %
####
qdffrac 8
1
nmru
1
2
% qdffrac      %
####
init_rxn_mru 13
1
nmru
1
1
1
####
UADJ 8
1
nmru
1
2
% windad       %
####
pmac_sat 9
1
nmru
1
2
% pmacst       %
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
WEI 8
1
nmru
1
2
0.0
####
gw_loss_k 4
1
nmru
1
2
% gwloss       %
####
irrig_sched_ext 16
1
nmru
1
1
0
####
chem_sim 9
1
one
1
1
1
####
MFMAX 8
1
nmru
1
2
% meltmx       %
####
SI 8
1
nmru
1
2
% sn_thr       %
####
T_decay 10
1
nmru
1
1
0
####
mru_slope 10
1
nmru
1
2
0.4
####
nacsc 8
1
nmru
1
1
11
####
srain_intcp 12
1
nmru
1
2
% sraini       %
####
mru_percent_imperv 19
1
nmru
1
2
0.0
####
crad_coef 10
1
one
1
2
% crdcof       %
####
snowmelt_18O_depl 14
1
nmru
1
2
% sO_dpl       %
####
print_vse 11
1
one
1
1
0
####
nchan_d 8
1
nchan
1
1
5
####
st 8
2
nac
nmru
11
2
20.9
8.8
7.8
7.21
6.76
6.39
6.03
5.65
5.24
4.71
2.62
####
mru_psta 9
1
nmru
1
1
1
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
trxn_sat_c_adj 9
1
nmru
1
2
% tscadj       %
####
xdebug_stop 4
1
one
1
1
2572
####
NMF 8
1
nmru
1
2
% nmeltf       %
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
basin_area 14
1
one
1
2
3.26
####
init_rxn_hydro 15
1
nhydro
1
1
-1
####
td 8
1
nmru
1
2
% td           %
####
transp_end 11
1
nmru
1
1
13
####
mru_elev 9
1
nmru
1
2
615.0
####
C_units 8
1
nchemvar
10
1
13
13
7
13
7
13
7
13
7
7
####
iout 5
1
one
1
1
2
####
print_explanation 4
1
one
1
1
0
####
MRUDEPL 8
1
nmru
1
1
2
####
tl 8
1
nmru
1
2
6.68
####
SUBRATE 8
1
nmru
1
2
0.01
####
gwbnd_len2 11
1
nmru
1
2
0.0
####
tmax_lapse 11
1
nmonths
12
2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
####
cov_type 9
1
nmru
1
1
3
####
init_exch_mru 14
1
nmru
1
1
1
####
gwbnd_len1 11
1
nmru
1
2
0.0
####
mru_lat 8
1
nmru
1
2
18.3
####
ALAT 5
1
one
1
2
18.3
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
covden_win 11
1
nmru
1
2
% cdenwn       %
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
irrig_int_max 14
1
nmru
1
2
100.0
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
riparian_thresh 6
1
nmru
1
2
% ripthr       %
####
irrig_int_init 15
1
nmru
1
2
0.0
####
print_freq 11
1
one
1
1
9
####
mru_area 9
1
nmru
1
2
3.26
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
c_can_depth 12
1
one
1
2
% candep       %
####
chem_ext 9
1
one
1
1
1
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
trxn_sat_stat 9
1
nmru
1
1
1
####
snowmelt_D_depl 14
1
nmru
1
2
% sD_dpl       %
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
s_theta_fc 11
1
nmru
1
2
% s_thfc       %
####
dtinit 7
1
one
1
2
24.0
####
trxn_sat_days 9
1
one
1
2
% tsdays       %
####
To 8
1
nmru
1
2
% transm       %
####
trxn_uz_c_adj 9
1
nmru
1
2
% tucadj       %
####
xk0 8
1
nmru
1
2
% xkvert       %
####
s_rock_depth 13
1
nmru
1
2
% s_rock       %
####
PLWHC 8
1
nmru
1
2
% pklwhc       %
####
s_ohoriz_depth 15
1
nmru
1
2
% s_ohor       %
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
MFMIN 8
1
nmru
1
2
% meltmn       %
####
tsta_elev 10
1
ntemp
1
2
482.0
####
init_kin_hydro 15
1
nhydro
1
1
-1
####
radj_wppt 14
1
one
1
2
% radjwp       %
####
src_ext_irrig 14
1
nmru
1
1
1
####
init_surf_hydro 16
1
nhydro
1
1
-1
####
infex 8
1
one
1
1
1
####
chv 8
1
one
1
2
% chvelo       %
####
s_porosity 11
1
nmru
1
2
% s_thpo       %
####
wrain_intcp 12
1
nmru
1
2
% wraini       %
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
s_theta_wp 11
1
nmru
1
2
% s_thwp       %
####
s_root_depth 13
1
nmru
1
2
% s_root       %
####
chan_loss_rate 4
1
one
1
2
% chloss       %
####
DAYGM 8
1
nmru
1
2
% daymlt       %
####
MBASE 8
1
nmru
1
2
% meltba       %
####
snow_intcp 11
1
nmru
1
2
% sn_int       %
####
init_exch_hydro 16
1
nhydro
1
1
-1
####
iso_n 8
1
nmru
1
3
0.5
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
print_type 11
1
one
1
1
2
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
kinset_table 13
2
nmru_res
nchem_sets
90
1
-1
-1
-1
1
1
1
1
2
2
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
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
s_satpref_zmin 15
1
nmru
1
2
% s_zmin       %
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
C_ires 7
1
nchemvar
10
1
1
4
4
5
5
7
7
8
8
99
####
init_soln_mru 14
1
nmru
1
1
1
####
init_soln_hydro 16
1
nhydro
1
1
3
####
sched_gw1 10
1
nmru
1
1
0
####
print_objfunc 14
1
one
1
1
0
####
dth 8
1
nmru
1
2
1.0
####
sched_gw2 10
1
nmru
1
1
0
####
s_theta_0 11
1
nmru
1
2
% s_thfc       %
####
q0 8
1
one
1
2
0.06
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
tmax_allsnow_c 14
1
one
1
2
% tmxals       %
####
tmin_lapse 11
1
nmonths
12
2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
8.2
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
5
0
0
5
0
0
####
temp_units 14
1
one
1
1
0
####
mru_tsta 9
1
nmru
1
1
1
####
szm 8
1
nmru
1
2
% recess       %
####
trxn_ohoriz_c_adj 9
1
nmru
1
2
% tocadj       %
####
s_satpref_zmax 15
1
nmru
1
2
% s_zmax       %
####
mru2chan 12
1
nmru
1
1
1
####
trxn_uz_days 9
1
one
1
2
% tudays       %
####
atmos_eq_ph 14
1
one
1
1
0
####
ACUMX 8
1
nmru
1
2
0.0
####
init_eq_ph_hydro 17
1
nhydro
1
1
2
####
init_soln_ppt 14
1
one
1
1
0
####
init_eq_ph_mru 15
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
irrig_sched_int 16
1
nmru
1
1
0
####
sbar0 6
1
nmru
1
2
0.3
####
xk_cv 8
1
nmru
1
2
% xkcvar       %
####
src_gw1 7
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
mru_area_frac 14
1
nmru
1
2
1.0
####
hf 8
1
nmru
1
2
0.01
####
crad_exp 9
1
one
1
2
% crddexp      %
