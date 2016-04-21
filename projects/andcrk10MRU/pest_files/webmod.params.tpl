ptf %
Andrews Creek
Version: 1.7
** Dimensions **
####
nrain
2
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
nhcs
10
####
nac
11
####
nsol
1
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
3
####
nhum
1
####
nac_nmru_nresinp
2310
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
2
####
ndeplval
22
####
nmru
10
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
10
2
% tmxadj01     %
% tmxadj02     %
% tmxadj03     %
% tmxadj04     %
% tmxadj05     %
% tmxadj06     %
% tmxadj07     %
% tmxadj08     %
% tmxadj09     %
% tmxadj10     %
####
irrig_int_src 14
1
nmru
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
xdebug_start 4
1
one
1
1
0
####
C_metric 9
1
nchemvar
10
1
3
5
3
3
3
4
4
4
4
2
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
10
10
5
5
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
10
2
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
% pmacro       %
####
snow_ion_factor 14
1
nmru
10
2
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
% sn_ion       %
####
solnset_table 14
2
nmru_res
nhcs
90
1
2
2
2
3
30
40
40
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
10
2
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
% sn_flx       %
####
transp_beg 11
1
nmru
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
eq_phset_table 15
2
nmru_res
nhcs
90
1
-1
-1
-1
3
3
3
3
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
10
2
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
% s_satk       %
####
mru_aspect 11
1
nmru
10
2
135.0
0.0
135.0
45.0
155.0
335.0
75.0
330.0
67.0
330.0
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
110
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
radj_sppt 14
1
one
1
2
% radjsp       %
####
init_kinset_mru 13
1
nmru
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
trxn_ohoriz_days 9
1
one
1
2
% todays       %
####
init_surfset_mru 14
1
nmru
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
NEGHSI 8
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
####
transp_tmax_c 12
1
nmru
10
2
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
% trdegc       %
####
covden_sum 11
1
nmru
10
2
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
% cdensm       %
####
tmin_adj 9
1
nmru
10
2
% tmnadj01     %
% tmnadj02     %
% tmnadj03     %
% tmnadj04     %
% tmnadj05     %
% tmnadj06     %
% tmnadj07     %
% tmnadj08     %
% tmnadj09     %
% tmnadj10     %
####
qdffrac 8
1
nmru
10
2
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
% qdffrac      %
####
init_rxnset_mru 13
1
nmru
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
UADJ 8
1
nmru
10
2
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
% windad       %
####
pmac_sat 9
1
nmru
10
2
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
% pmacst       %
####
WEI 8
1
nmru
10
2
15.0
15.0
15.0
15.0
15.0
15.0
15.0
15.0
15.0
15.0
####
gw_loss_k 4
1
nmru
10
2
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
% gwloss       %
####
irrig_sched_ext 16
1
nmru
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
10
2
% meltmx01     %
% meltmx02     %
% meltmx03     %
% meltmx04     %
% meltmx05     %
% meltmx06     %
% meltmx07     %
% meltmx08     %
% meltmx09     %
% meltmx10     %
####
SI 8
1
nmru
10
2
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
% sn_thr       %
####
T_decay 10
1
nmru
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
mru_slope 10
1
nmru
10
2
1.0
0.686
1.022
0.65
0.717
0.548
0.445
0.747
1.025
0.825
####
nacsc 8
1
nmru
10
1
11
11
11
11
11
11
11
11
11
11
####
srain_intcp 12
1
nmru
10
2
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
% sraini       %
####
mru_percent_imperv 19
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
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
10
2
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
% sO_dpl       %
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
2
####
st 8
2
nac
nmru
110
2
7604.02
82.49
45.99
31.12
27.43
23.74
20.02
16.29
12.59
8.89
1.47
7600.64
90.81
61.33
53.93
46.45
39.08
31.63
24.11
16.73
9.41
1.99
1469.55
66.81
37.35
28.75
23.03
18.74
15.89
11.6
8.7
4.43
1.57
2377.58
136.88
57.74
43.97
34.82
27.74
23.14
18.46
13.97
7.0
2.38
4912.55
97.48
59.11
44.75
35.16
30.28
25.53
15.96
11.18
6.37
1.59
8050.95
158.87
96.03
72.34
56.64
40.85
33.16
25.27
17.44
9.54
1.74
3694.25
113.43
71.07
48.89
38.72
31.61
24.26
20.83
13.49
5.8
2.85
4014.77
68.26
41.14
33.33
29.4
25.49
17.63
13.74
9.81
5.8
1.98
2628.35
67.89
47.53
37.35
29.67
24.55
19.42
16.82
11.72
6.59
1.46
3855.45
84.08
57.97
46.54
39.22
31.59
24.14
20.41
12.9
5.37
1.64
####
mru_psta 9
1
nmru
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
d 12
2
ntopchan
nchan
2
2
0.0
600.0
####
trxn_sat_c_adj 9
1
nmru
10
2
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
% tscadj       %
####
xdebug_stop 4
1
one
1
1
0
####
NMF 8
1
nmru
10
2
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
% nmeltf       %
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
1.74
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
10
2
% td           %
% td           %
% td           %
% td           %
% td           %
% td           %
% td           %
% td           %
% td           %
% td           %
####
transp_end 11
1
nmru
10
1
9
9
9
9
9
9
9
9
9
9
####
mru_elev 9
1
nmru
10
2
3479.0
3318.0
3545.0
3437.0
3645.0
3667.0
3348.0
3429.0
3683.0
3549.0
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
13
4
4
4
4
1
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
10
1
2
2
2
2
2
2
2
2
2
2
####
tl 8
1
nmru
10
2
78.0
111.7
36.5
125.9
65.8
84.0
90.5
102.6
43.6
54.3
####
SUBRATE 8
1
nmru
10
2
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
####
gwbnd_len2 11
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
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
cov_type 9
1
nmru
10
1
3
3
3
3
3
3
3
3
3
3
####
init_exchset_mru 14
1
nmru
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
gwbnd_len1 11
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
####
mru_lat 8
1
nmru
10
2
40.29174
40.28814
40.29062
40.2881
40.2896
40.28573
40.28739
40.28624
40.2835
40.28335
####
ALAT 5
1
one
1
2
40.28738
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
10
2
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
% cdenwn       %
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
10
2
100.0
100.0
100.0
100.0
100.0
100.0
100.0
100.0
100.0
100.0
####
ach 12
2
ntopchan
nchan
2
2
0.0
1.0
####
riparian_thresh 6
1
nmru
10
2
65.0
65.0
65.0
65.0
65.0
100.0
65.0
65.0
65.0
65.0
####
irrig_int_init 15
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
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
10
2
0.269
0.132
0.181
0.069
0.203
0.297
0.051
0.09
0.276
0.172
####
exchset_table 14
2
nmru_res
nhcs
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
snowmelt_D_depl 14
1
nmru
10
2
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
% sD_dpl       %
####
rain_adj 9
2
nmru
nmonths
120
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
10
2
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
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
10
2
% transm01     %
% transm02     %
% transm03     %
% transm04     %
% transm05     %
% transm06     %
% transm07     %
% transm08     %
% transm09     %
% transm10     %
####
trxn_uz_c_adj 9
1
nmru
10
2
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
% tucadj       %
####
xk0 8
1
nmru
10
2
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
% xkvert       %
####
s_rock_depth 13
1
nmru
10
2
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
% s_rock       %
####
PLWHC 8
1
nmru
10
2
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
% pklwhc       %
####
s_ohoriz_depth 15
1
nmru
10
2
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
% s_ohor       %
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
10
2
% meltmn01     %
% meltmn02     %
% meltmn03     %
% meltmn04     %
% meltmn05     %
% meltmn06     %
% meltmn07     %
% meltmn08     %
% meltmn09     %
% meltmn10     %
####
tsta_elev 10
1
ntemp
1
2
3150.0
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
10
2
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
% s_thpo       %
####
wrain_intcp 12
1
nmru
10
2
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
% wraini       %
####
iso_fac 8
1
nmru
10
3
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
% isofac       %
####
rxnset_table 13
2
nmru_res
nhcs
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
10
2
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
% s_thwp       %
####
s_root_depth 13
1
nmru
10
2
% s_root       %
% s_root       %
% s_root       %
% s_root       %
% s_root       %
% s_root       %
% s_root       %
% s_root       %
% s_root       %
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
10
2
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
% daymlt       %
####
MBASE 8
1
nmru
10
2
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
% meltba       %
####
snow_intcp 11
1
nmru
10
2
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
% sn_int       %
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
10
3
0.5
0.5
0.5
0.5
0.5
0.5
0.5
0.5
0.5
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
0.015853
0.018364
0.019156
0.017939
0.015147
0.012033
0.010571
0.011373
0.014426
0.018272
0.020824
0.018511
####
kinset_table 13
2
nmru_res
nhcs
90
1
-1
-1
-1
1
1
1
1
1
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
####
s_satpref_zmin 15
1
nmru
10
2
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
% s_zmin       %
####
snow_adj 9
2
nmru
nmonths
120
2
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
% sn_5_6      %
% sn_5_6      %
% snadja      %
% snadja      %
% snadja      %
% snadja      %
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
8
5
8
5
8
99
####
init_solnset_mru 14
1
nmru
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
init_soln_hydro 16
1
nhydro
1
1
3
####
iso_theta 8
1
nmru
10
3
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
sched_gw1 10
1
nmru
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
10
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
####
sched_gw2 10
1
nmru
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
s_theta_0 11
1
nmru
10
2
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
% s_thfc       %
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
nhcs
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
5
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
1
####
mru_tsta 9
1
nmru
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
szm 8
1
nmru
10
2
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
% recess       %
####
trxn_ohoriz_c_adj 9
1
nmru
10
2
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
% tocadj       %
####
s_satpref_zmax 15
1
nmru
10
2
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
% s_zmax       %
####
mru2chan 12
1
nmru
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
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
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
init_eq_phset_mru 15
1
nmru
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
trxn_uz_stat 9
1
nmru
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
sbar0 6
1
nmru
10
2
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
####
xk_cv 8
1
nmru
10
2
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
% xkcvar       %
####
src_gw1 7
1
nmru
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
src_gw2 7
1
nmru
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
mru_area_frac 14
1
nmru
10
2
0.155
0.076
0.104
0.04
0.117
0.17
0.029
0.052
0.158
0.099
####
hf 8
1
nmru
10
2
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
####
crad_exp 9
1
one
1
2
% crddexp      %
