ptf %
PHASES
Oligoclase
Ca0.27Na0.73Al1.27Si2.73O8 + 8H2O = 0.27Ca+2 + 0.73Na+ + 1.27Al(OH)4- + 2.73H4SiO4
	log_k          -19.0  # Ab -18; An -19.714
 
Biotite
	#KMg3AlSi3O10(OH)2 + 6H3O+ + 4H2O = K+ + 3Mg+2 + Al(OH)4- + 3H4SiO4
     	K0.98Mg1.0Fe1.33Si0.18Al0.33Al1.35Si2.65O10(OH)2 + 10.68H3O+ = 0.98K+ + Mg+2 + 1.33Fe+2 + 1.68Al+3 + 2.83H4SiO4 + 0.68H2O +10.68 H2O
	#K0.98Mg1.0Fe1.33Al0.33Al1.35Si2.83O10(OH)2 + 3.96H3O+ + 6.04H2O = 0.98K+ + Mg+2 + 1.33Fe+2 + 1.68Al(OH)4- + 2.83H4SiO4
	#log_k 12.7 # Kmica phreeqc.dat
        log_k 20.0 # Kmica phreeqc.dat
Chlorite
	#Mg5Al2Si3O10(OH)8 + 16H3O+ = 5Mg+2 + 2Al+3 + 3H4SiO4 + 6H2O
	Mg1.81Fe2.72Al1.39Al1.23Si2.77O10(OH)8 + 16.92H3O+ = 1.81Mg+2 + 2.72Fe+2+ 1.39Al+3 + 1.23Al+3 + 2.77H4SiO4 + 6.92H2O + 16.92 H2O
        log_k           68.38 # phreeqc.dat Chlorite14A
        delta_h -151.494 kcal
 
Smectite-Illite              45
	#K0.6Mg0.25Al2.3Si3.5O10(OH)2 + 11.2H2O = 0.6K+ +0.25Mg+2 + 2.3Al(OH)4- + 3.5H4SiO4 + 1.2H3O+
	#K0.32Ca0.1Fe0.25Mg0.39Al1.47Al0.46Si3.54O10(OH)2 + 11.2H2O = 0.32K+ +0.1Ca+2 + 0.25Fe+2 + 0.39Mg+2 + 1.47Al(OH)4- \
	#		+ 0.46Al(OH)4- + 3.54H4SiO4 + 1.2H3O+
	# adjusted Ca to balance O10(OH)2
	K0.32Ca0.225Fe0.25Mg0.39Al1.47Al0.46Si3.54O10(OH)2 + 9.88H2O + 0.12H3O+ = 0.32K+ +0.225Ca+2 + 0.25Fe+2 + 0.39Mg+2 + 1.47Al(OH)4- \
			+ 0.46Al(OH)4- + 3.54H4SiO4 + 0.12 H2O
        log_k           -40.267   # Illite, phreeqc.dat
        delta_h 54.684 kcal
END
SOLUTION_SPECIES
   H2O + 0.01e- = H2O-0.01
   log_k     -9.0
 
NO3- + 2 H3O+ + 2 e- = NO2- + 3H2O
     #log_k           28.570
     log_k           -100
 
 
2 NO3- + 12 H3O+ + 10e- = N2 + 18 H2O
     #log_k           207.080
     log_k           -100
 
SOLUTION 0 Precipitation
	-units	mg/kgw
	pH      5.3
	Ca	0.22
	Mg	0.028
	Na	0.08
	K	0.03
	Amm	0.14
	Cl	0.1
	S	0.655
	N(5)	0.21
	Si	0.001 
	Doc	0.001
	Al  	0.001 kaolinite 0.0
	Fe	0.001
	O(0)	1 O2(g) -1.8
	[18O]   -8.0
SOLUTION 1 E.T. or DI
SOLUTION 2 Canopy,  Snow
	-units	mg/kgw
	pH      5.3
	temp	2
	Ca	0.22
	Mg	0.028
	Na	0.08
	K	0.03
	Amm	0.14
	Cl	0.1
	S	0.655
	N(5)	0.21
	Si	0.001
	Doc	0.001
	Al  	0.001
	Fe	0.001
	O(0)	1 O2(g) -1.8
	[18O]	-15.6
SOLUTION 3 Unsaturated Zone and Stream Segments
    temp      10.00
    pH        5.59
    pe        16.15
    redox     pe
    units     mol/kgw
    density   1
	Al               3.578e-008
	Amm              3.129e-007
	C                3.840e-004
	Ca               3.310e-005
	Cl               1.567e-006
	Doc              1.352e-009
	Fe               4.598e-013
	K                3.491e-006
	Mg               1.026e-005
	N(5)             2.155e-005
	Na               2.279e-005
	S                1.365e-005
	Si               4.466e-005
#    Ca        1.92E+00
#    Mg        3.45E-01
#    Na        9.38E-01
#    K         2.36E-01
#    Cl        9.32E-02
#    Amm       1.90E-05
#    S         2.61
#    N(5)      0.263
#    Si        3.68
#    Doc       0.001
#    Al        1.14E-04
#    Fe        7.73E-09
#    Alkalinity 4.3
    O(0)      1 O2(g)      -1.55
    [18O]     -15.15
    -water    1 # kg
SOLUTION 4 Saturated Zone
    temp      1.52
    pH        5.97
    pe        16.39
    redox     pe
    units     mol/kgw
    density   1
	Al               4.704e-009
	Amm              1.296e-009
	C                4.165e-004
	Ca               5.321e-005
	Cl               2.923e-006
	Doc              3.646e-008
	Fe               1.679e-013
	K                5.355e-006
	Mg               1.476e-005
	N(5)             2.442e-005
	Na               4.868e-005
	S                3.054e-005
	Si               8.033e-005
#    Ca        1.92E+00
#    Mg        3.45E-01
#    Na        9.38E-01
#    K         2.36E-01
#    Cl        9.32E-02
#    Amm       1.90E-05
#    S         2.61
#    N(5)      0.263
#    Si        3.68 gibbsite 0.0
#    Doc       0.001
#    Al        1.14E-04 kaolinite  0
#    Fe        7.73E-09 Goethite   0
#    Alkalinity 4.3
    O(0)      1 O2(g)      -2.41
    [18O]     -18.34
    -water    1 # kg
END
RATES
Oligoclase
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Oligoclase", -9.67, 65.0, 0.457, 1, 1,   -11.84, 69.8, 0, 1, 1,   -99, 0, 0, 1, 1
30  RESTORE 20
40  READ name$
50  DIM p(3,5)
60  FOR i = 1 to 3
70    FOR j = 1 to 5
80      READ p(i,j)
100   NEXT j
110 NEXT i
120 omega = SR(name$)
130 GOSUB 2000  # calculate rates
140 moles = 10^PARM(1) * rate * time
150 SAVE moles
160 PUT (rate, 1,10)
180 END
 
2000 REM Palandri and Kharaka rate subroutine
2010 R = 0.00831470	# kJ/deg-mol
2020 aH = act("H3O+")
2030 FOR i = 1 to 3
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) # * (1 - omega^p(i,4))^p(i,5)
#2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) * (1 - omega^p(i,4))^p(i,5)
2050   rate = rate + rate_i
2060 NEXT i
2070 RETURN
		-end
 
Biotite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Biotite", -9.84, 22.0, 0.525, 1, 1,   -12.55, 22.0, 0, 1, 1,   -99, 0, 0, 1, 1
30  RESTORE 20
40  READ name$
50  DIM p(3,5)
60  FOR i = 1 to 3
70    FOR j = 1 to 5
80      READ p(i,j)
100   NEXT j
110 NEXT i
120 omega = SR(name$)
130 GOSUB 2000  # calculate rates
140 moles = 10^PARM(1) * rate * time
150 SAVE moles
160 PUT (rate, 1,10)
180 END
 
2000 REM Palandri and Kharaka rate subroutine
2010 R = 0.00831470	# kJ/deg-mol
2020 aH = act("H3O+")
2030 FOR i = 1 to 3
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) # * (1 - omega^p(i,4))^p(i,5)
#2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3)  * (1 - omega^p(i,4))^p(i,5)
2050   rate = rate + rate_i
2060 NEXT i
2070 RETURN
	-end
 
Chlorite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Chlorite", -11.11, 88.0, 0.5, 1, 1,   -12.52, 88.0, 0, 1, 1,   -99, 0, 0, 1, 1
30  RESTORE 20
40  READ name$
50  DIM p(3,5)
60  FOR i = 1 to 3
70    FOR j = 1 to 5
80      READ p(i,j)
100   NEXT j
110 NEXT i
120 omega = SR(name$)
130 GOSUB 2000  # calculate rates
140 moles = 10^PARM(1) * rate * time
150 SAVE moles
160 PUT (rate, 1,10)
180 END
 
2000 REM Palandri and Kharaka rate subroutine
2010 R = 0.00831470	# kJ/deg-mol
2020 aH = act("H3O+")
2030 FOR i = 1 to 3
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) # * (1 - omega^p(i,4))^p(i,5)
#2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) * (1 - omega^p(i,4))^p(i,5)
2050   rate = rate + rate_i
2060 NEXT i
2070 RETURN
	-end
 
Calcite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Calcite", -0.30, 14.4, 1.0, 1, 1,  -5.81, 23.5, 0, 1, 1,   -3.48, 35.4, 1.0, 1, 1,
30  RESTORE 20
40  READ name$
50  DIM p(3,5)
60  FOR i = 1 to 3
70    FOR j = 1 to 5
80      READ p(i,j)
100   NEXT j
110 NEXT i
120 omega = SR(name$)
130 GOSUB 2000  # calculate rates
140 moles = 10^PARM(1) * rate * time
150 SAVE moles
160 PUT (rate, 1,10)
180 END
 
2000 REM Palandri and Kharaka rate subroutine
2010 R = 0.00831470	# kJ/deg-mol
2020 aH = act("H3O+")
2030 FOR i = 1 to 3
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) * (1 - omega^p(i,4))^p(i,5)
2050   rate = rate + rate_i
2060 NEXT i
2070 RETURN
	-end
 
Smectite-Illite # Secondary phase, rate not used
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Smectite-Illite", -10.98, 23.6, 0.34, 1, 1,   -12.78, 35.0, 0, 1, 1,   -16.52, 58.9, -0.4, 1, 1
30  RESTORE 20
40  READ name$
50  DIM p(3,5)
60  FOR i = 1 to 3
70    FOR j = 1 to 5
80      READ p(i,j)
100   NEXT j
110 NEXT i
120 omega = SR(name$)
130 GOSUB 2000  # calculate rates
140 moles = PARM(1) * rate * time
150 SAVE moles
160 PUT (rate, 1,10)
180 END
 
2000 REM Palandri and Kharaka rate subroutine
2010 R = 0.00831470	# kJ/deg-mol
2020 aH = act("H3O+")
2030 FOR i = 1 to 3
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) * (1 - omega^p(i,4))^p(i,5)
2050   rate = rate + rate_i
2060 NEXT i
2070 RETURN
	-end
 
Nitrification
 -start
   1 REM parm(1) is log factor to speed or slow reaction
   10 rate = 10^parm(1) * TOT("Amm") * MOL("O2") / (MOL("O2") + 1e-5)
   20 moles = rate * TIME
   100 SAVE moles
 -end
 
 Amm_assimilation
  -start
    1 REM parm(1) is log factor to speed or slow reaction
    10 rate = 10^parm(1) * TOT("Amm")
    20 moles = rate * TIME
    100 SAVE moles
 -end
 
 Potassium_cycling
  -start
    1 REM parm(1) is summer uptake rate constant, units per day
    2 REM parm(2) is winter release rate constant, units per day
    4 REM Winter
    5 transp = callback(cell_no,dummy,"transp_on")
    10 If transp > 0.1 then goto 100
    20 k = parm(2)/86400
    30 rate = k * M
    40 moles = rate * TIME
    50 goto 200
    100 REM Summer
    120 k = parm(1)/86400
    130 rate = k * TOT("K")
    140 moles = rate * TIME
    200 SAVE moles
 -end
 
 NO3_assimilation
  -start
    1 REM parm(1) is log factor to speed or slow reaction
    10 rate = 10^parm(1) * TOT("N(5)")
    20 moles = rate * TIME
    100 SAVE moles
 -end
 
Pyrite_O2
  -start
   1 rem        Handcart Gulch model numbers, round 1, Add Pyrite
   2 rem        parm(1) = log10(A/V, m^2/L water)      parm(2) = exp for (m/m0)
   3 rem        parm(3) = exp for O2            parm(4) = exp for H+
#   10 if (m <= 0) then goto 200
   20 if (si("Pyrite") >= 0) then goto 200
   25  rate = -9.14 + parm(1) + parm(3)*lm("O2") + parm(4)*lm("H3O+") + parm(2)*log10(m/m0)
   26  rate = 10^rate * (1 - SR("pyrite"))
#   26  rate = 10^rate
   30  moles = rate * time
#   40 if (moles > m) then moles = m
   200 save moles
  -end
 
#Pyrite_O2_xtreme
#  -start
#   1 rem        Handcart Gulch model numbers, round 1, Linear increase to add sulfate as water table drops
#   2 rem        parm(1) = log10(A/V, m^2/L water)      parm(2) = exp for (m/m0)
#   3 rem        parm(3) = exp for O2            parm(4) = exp for H+
##   10 if (m <= 0) then goto 200
#   20 if (si("Pyrite") >= 0) then goto 200
## log pO2 of -1.8 converted to log molal = -4.476
##   25  rate = -9.14 + parm(1) + parm(3)*-4.476 + parm(4)*lm("H3O+") + parm(2)*log10(m/m0)
#   25  rate = -9.14 + parm(1) + parm(3)*lm("O2")+ parm(4)*lm("H3O+") + parm(2)*log10(m/m0)
#   26  rate = 10^rate * (1 - SR("pyrite"))
##   26  rate = 10^rate
#   27 gwsto = callback(cell_no,dummy,"gwsto")
## linear enhancement when gwsto below parm(5)
#   28 if(gwsto < parm(5)) then rate = rate * (parm(5)-gwsto) else rate = 0
#   30  moles = rate * time
##   40 if (moles > m) then moles = m
#   200 save moles
#  -end
 
END
 
 
KINETICS 1 O-Horizon and UZ (normal and preferential)
 
Oligoclase
	-parm % koligu       %                # fit the first parameter as log10 Surface Area
        -m    1000

Biotite
	-parm % kbiotu       %              # fit the first parameter as log10 Surface Area
        -m    1000

Chlorite
	-parm % kchlru       %  # fit the first parameter as log10 Surface Area
        -m    1000

Calcite
        -parm % kcalcu       % # fit the first parameter as log10 Surface Area
        -m    1000
        
Pyrite_O2
    -formula  Pyrite  1
    -m        1000
    -m0       1000
    -parms    % kpyr_u       % 0 0.5 -0.11
    -tol      1e-008	

Nitrification
    -formula  Amm -1 NH3 +1
    -m        1000
    -parms    % knitru       %

Potassium_Cycling
    -formula KHCO3 1  Mg(HCO3)2 1 # Ca(HCO3)2 1 
    -m       0
    -parms   % kK_sum       % % kK_win       %
 
#Amm_assimilation
#	-formula AmmHHCO3 -1
#	-m       1000
#	-parms   -5.5
	
#NO3_assimilation
#	-formula NO3 -1 HCO3 1
#	-m       1000
#	-parms   -7.5	
 
-step 86400
#        -cvode
 
END
 
 
KINETICS 2 Sat and Preferential Sat 
Oligoclase
	-parm % koligs       %              # fit the first parameter as log10 Surface Area
        -m    1000

Biotite
	-parm % kbiots       %              # fit the first parameter as log10 Surface Area
        -m    1000

Chlorite
	-parm % kchlrs       %              # fit the first parameter as log10 Surface Area
        -m    1000

Calcite
        -parm % kcalcs       %              # fit the first parameter as log10 Surface Area
        -m    1000
        
Pyrite_O2
    -formula  Pyrite 1
    -m        1000
    -m0       1000
    -parms    % kpyr_s       % 0 0.5 -0.11
    -tol      1e-008	

Nitrification
    -formula Amm -1 NH3 +1
    -m       1000
    -parms   % knitrs       %


-step 86400

END

EQUILIBRIUM_PHASES 0 Equilibrate precipitation and irrigation with atmospheric pO2 and pCO2 @ 3500 meters above mean sea level, in log ppm
	O2(g)         -8.6551E-01 100      #
	CO2(g)        -3.643659    10      # ~ PCO2
END
 
EQUILIBRIUM_PHASES 2 Equilibrate stream pO2 and pCO2, in log ppm
# should be greater than atmosphere 
	O2(g)       % O2str        % 100      # 
	CO2(g)      % CO2str       %  10      # ~ PCO2 
END
EQUILIBRIUM_PHASES 30 Unsaturated zone
	Kaolinite   0 1e-2 # precipitate
	Goethite    0 1e-2 # precipitate
	Gibbsite    0 0 # precipitate
#	Calcite    -.5 0
	Smectite-Illite % ksmecu       % 0 # precipitate  # fit the first 0, equilibrium constant
#	Calcite     0 0
#	O2(g)       -0.9 100      # atmospheric, could adjust Use this for gradient of O2, this being the most oxygen rich
	O2(g)       % O2uz         % 100      # Less oxygen and more CO2 as a result of root respiration and pyrite oxidation
       CO2(g)       % CO2uz        %  10      # UZ PCO2, greater than -3.65 which is atmospheric at 3500 m amsl
END
EQUILIBRIUM_PHASES 4 Saturated zone
	Kaolinite   0 1e-2 # precipitate
	Goethite    0 1e-2 # precipitate
	Gibbsite    0 0 # precipitate
#	Calcite    0 0
	Smectite-Illite % ksmecs       % 0 # precipitate  # fit the first 0, equilibrium constant
END
EXCHANGE 1
X	.001
-eq solution 2
SURFACE 1
Hfo_w .001 1 50
-eq solution 2
REACTION_TEMPERATURE 1
0
END
#USE SOLUTION 3
#USE EQUILIBRIUM_PHASES 30
#USE EXCHANGE 1
#USE KINETICS 1
#SAVE SOLUTION 3
#END
 
 
PRINT
# -RESET FALSE
-STATUS FALSE

USER_PRINT
10 PRINT "Oligoclase      ", KIN("Oligoclase") - 1e3
20 PRINT "Biotite         ", KIN("Biotite") - 1e3
30 PRINT "Chlorite        ", KIN("Chlorite") - 1e3
40 PRINT "Calcite         ", KIN("Calcite") - 1e3
50 PRINT "Pyrite_O2       ", KIN("Pyrite_O2") - 1e3
60 PRINT "Nitrification   ", KIN("Nitrification") - 1e3
