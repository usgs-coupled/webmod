ptf %
# Par2par file for PRMS parameter estimation
* parameter data
k_decr = %    k_decr    %
k_decr_pyr = %  k_decr_pyr  %
k_decr_nit = %  k_decr_nit  %
CO2uz =  %    CO2uz     %
CO2str = %    CO2str    %
O2uz =   %    O2uz      %
O2str =  %    O2str     %
koligu = %    koligu    %
kbiotu = %    kbiotu    %
kchlru = %    kchlru    %
kcalcu = %    kcalcu    %
kpyr_u = %    kpyr_u    %
knitru = %    knitru    %
kK_sum = %    kK_sum    %
kK_win = %    kK_win    %
kpyrxt = %    kpyrxt    %
tpyrxt = %    tpyrxt    %
qdffrac = %    qdffrac   %
td =      %    td        %
adjmix = %    adjmix    %
candep = %    candep    %
ccvint = %    ccvint    %
ccvslp = %    ccvslp    %
chloss = %    chloss    %
chvelo = %    chvelo    %
cdensm = %    cdensm    %
cdenwn = %    cdenwn    %
crdcof = %    crdcof    %
crddexp = %    crddexp   %
daymlt = %    daymlt    %
pancof = %    pancof    %
gwloss = %    gwloss    %
hamcof = %    hamcof    %
meltba = %    meltba    %
meltmx = %    meltmx    %
meltmn = %    meltmn    %
nmeltf = %    nmeltf    %
pklwhc = %    pklwhc    %
pmacst = %    pmacst    %
pmacro = %    pmacro    %
pptrad = %    pptrad    %
radjsp = %    radjsp    %
radjwp = %    radjwp    %
radmax = %    radmax    %
rainad = %    rainad    %
ripthr = %    ripthr    %
s_ohor = %    s_ohor    %
s_thpo = %    s_thpo    %
s_rock = %    s_rock    %
s_root = %    s_root    %
s_satk = %    s_satk    %
s_zmax = %    s_zmax    %
s_zmin = %    s_zmin    %
s_thfc = %    s_thfc    %
s_thwp = %    s_thwp    %
sn_thr = %    sn_thr    %
sn_adj = %    sn_adj    %
sn_int = %    sn_int    %
sn_ion = %    sn_ion    %
sO_dpl = %    sO_dpl    %
sD_dpl = %    sD_dpl    %
sraini = %    sraini    %
wraini = %    wraini    %
recess = %    recess    %
transm = %    transm    %
sn_flx = %    sn_flx    %
tmxadj = %    tmxadj    %
tmxalr = %    tmxalr    %
tmxals = %    tmxals    %
tmnadj = %    tmnadj    %
trdegc = %    trdegc    %
windad = %    windad    %
xkcvar = %    xkcvar    %
xkvert = %    xkvert    %

koligs = ( koligu + k_decr )
kbiots = ( kbiotu + k_decr )
kchlrs = ( kchlru + k_decr )
kcalcs = ( kcalcu + k_decr )
kpyr_s = ( kpyr_u + k_decr_pyr )
knitrs = ( knitru + k_decr_nit)


* template files
/home/dlpark/WEB_PEST/params_andcrk.tpl andcrk.params
/home/dlpark/WEB_PEST/pqi_andcrk.tpl andcrk.dat.pqi
