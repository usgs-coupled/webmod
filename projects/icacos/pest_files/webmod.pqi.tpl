ptf %
TITLE Reactions for Icacos
PHASES
Andesine   # Kinetics
Ca0.40Na0.60Al1.4Si2.6O8 + 8H2O = 0.4Ca+2 + 0.6Na+ + 1.4Al(OH)4- + 2.6H4SiO4
	log_k          -19.0  # Ab -18; An -19.714

Albite
        NaAlSi3O8 + 8 H2O = Na+ + Al(OH)4- + 3 H4SiO4
        log_k           -18.002
        delta_h 25.896 kcal
        
Anorthite
        CaAl2Si2O8 + 8 H2O = Ca+2 + 2 Al(OH)4- + 2 H4SiO4
        log_k           -19.714
        delta_h 11.580 kcal

Kaolinite
        Al2Si2O5(OH)4 + 6 H3O+ = 7H2O + 2 H4SiO4 + 2 Al+3
        log_k           7.435
        delta_h -35.300 kcal
       
K-spar     # Kinetics
Na0.09K0.91AlSi3O8 + 8H2O = 0.09Na+ + 0.91K+ + Al(OH)4- + 3H4SiO4
        log_k           -20.573
        delta_h 30.820  kcal

Biotite    # Kinetics
K0.85(Al0.10Si0.2Fe1.3Fe0.05Mg1.25)Si2.8Al1.2O10(OH)2 + 10.0H3O+ = 0.85K+ + 0.1Al+3 + 0.2H4SiO4 + 1.3Fe+2 + 0.05Fe+3 + 1.25Mg+2 + 2.8H4SiO4 + 1.2Al+3 + 10H2O 
log_k 20 # 12.7 phreeqc.dat

Hornblende # Kinetics
(Na0.34K0.05)(Ca1.71Mg2.84Fe2.06Al0.89)Al1.0Si6.68O22(OH)2 + 19.28H3O+ = \
  0.34Na+ + 0.05K+ + 1.71Ca+2 + 2.84Mg+2 + 2.06Fe+2 + 0.00Fe+3 + \
  0.89Al+3 + 1.0Al+3 + 6.68H4SiO4 + 16.56H2O
  log_k 0
  
Regolith_biotite    # Equilibrium phase
K0.65(Al1.10Si0.15Fe0.35Fe0.15Mg0.55)Si3.2Al0.8O10(OH)2 + 8.6H3O+ = \
  0.65K+ + 1.1Al+3 + 0.15H4SiO4 + 0.35Fe+2 + 0.15Fe+3 + 0.55Mg+2 + \
  3.2H4SiO4 + 0.8Al+3 + 7.2H2O 
       # log_k           -40.267   Illite, Al(OH)4-, phreeqc.dat
       # delta_h 54.684 kcal
       # Al+3 + 8 H2O = Al(OH)4- + 4 H3O+
       #        log_k           -22.7
       # delta_h 42.30           kcal
       # log_k = -40.267 - 1.9x-22.7 = 2.683
       # delta_h = 50.685 - 1.9x42.30 =  -29.685
       log_k 2.683
       delta_h -29.685

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
	Al  	0.001
	Fe	0.001
	O(0)	1 O2(g) -1.8
	[18O]   -1.5
SOLUTION 1 E.T. or DI
SOLUTION 2 Canopy,  Snow
	-units	mg/kgw
	pH      6.5
	temp	25
	Ca	3.0
	Mg	1.0
	Na	4.2
	K	0.6
	Amm	0.06
	Cl	6.0
	S	1.2
	N(5)	0.12
	Si	15.0
	Doc	0.001
	Al  	0.001
	Fe	0.001
	O(0)	1 O2(g) -1.8
	[18O]	-2.5
SOLUTION 3 Unsaturated Zone and Stream Segments
    temp      19.2
    pH        6.18
    -units     mol/kgw
    density   1
	Al               1.343e-009
	Amm              1.085e-008
	C                5.389e-004
	Ca               4.489e-005
	Cl               7.784e-005
	Doc              9.187e-011
	Fe               1.522e-013
	K                1.387e-005
	Mg               3.424e-005
	N(5)             6.074e-006
	Na               1.373e-004
	S                1.176e-005
	Si               1.217e-004
    O(0)	1 O2(g) -1.8
    [18O]	-2.59
    -water    1 # kg
SOLUTION 4 Saturated Zone
    temp      19.84
    pH        8.89
    -units     mol/kgw
    density   1
	Al               1.463e-007
	Amm              1.006e-010
	C                4.040e-004
	Ca               1.035e-004
	Cl               9.562e-005
	Doc              2.282e-010
	Fe               3.261e-014
	K                5.824e-023
	Mg               4.802e-005
	N                7.158e-006
	Na               2.705e-004
	S                1.202e-005
	Si               3.289e-004
    O(0)	1 O2(g) -2.8  # less oxygen in saturated zone to start
    [18O]	-2.43
    -water    1 # kg
END
RATES
Andesine
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Andesine", -8.88, 53.5, 0.541, 1, 1,   -11.47, 57.4, 0, 1, 1,   -99, 0, 0, 1, 1
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
Albite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Albite", -10.16, 65.0, 0.457, 1, 1,   -12.56, 69.8, 0, 1, 1,   -15.6, 71.0, -0.572, 1, 1
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
Anorthite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Anorthite", -3.50, 16.6, 1.411, 1, 1,   -9.12, 17.8, 0, 1, 1,   -99, 0, 0, 1, 1
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
		
K-spar
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "K-spar", -10.06, 51.7, 0.500, 1, 1,   -12.41, 38.0, 0, 1, 1,   -21.20, 94.1, -0.823, 1, 1
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

Hornblende
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Hornblende", -7.00, 75.5, 0.600, 1, 1,   -10.3, 94.4, 0, 1, 1,   -99, 0, 0, 1, 1
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
Kaolinite
	-start
10  REM PARM(1) is log10 surface area in m^2
20  DATA "Kaolinite", -11.31, 65.9, 0.777, 1, 1,   -13.18, 22.2, 0, 1, 1,   -17.05, 17.9, -0.472, 1, 1
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
#2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) # * (1 - omega^p(i,4))^p(i,5)
2040   rate_i = 10^p(i,1) * exp(-p(i,2)/R*(1/TK - 1/298.15)) * aH^p(i,3) * (1 - omega^p(i,4))^p(i,5)
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

END


KINETICS 1 O-Horizon and UZ (normal and preferential)
#-cvode

#Andesine
#	-parm % kandiu       %                # fit the first parameter as log10 Surface Area
#        -m    1000
        
#Albite
#	-parm  % kalbiu       %             # fit the first parameter as log10 Surface Area
#        -m    1000
        
K-spar
	-parm % kkspru       %                # fit the first parameter as log10 Surface Area
        -m    1000
        

#Biotite
#	-parm % kbiotu       %              # fit the first parameter as log10 Surface Area
#        -m    1000

#Hornblende
#	-parm % khnbdu       %              # fit the first parameter as log10 Surface Area
#        -m    1000

#Kaolinite
#	-parm  % kkaolu       %        # fit the first parameter as log10 Surface Area
#        -m    1000
        
#Pyrite_O2
#    -formula  Pyrite  1
#    -m        1000
#    -m0       1000
#    -parms    % kpyr_u       % 0 0.5 -0.11
#    -tol      1e-008	

Nitrification
    -formula  Amm -1 NH3 +1
    -m        1000
    -parms    % knitru       %


-step 86400
#        -cvode

END


KINETICS 2 Sat and Preferential Sat 
#-cvode
Andesine
	-parm % kandis       %                # fit the first parameter as log10 Surface Area
        -m    1000
        
#Albite
#	-parm  % kalbis       %             # fit the first parameter as log10 Surface Area
#        -m    1000
        
K-spar
	-parm % kksprs       %                # fit the first parameter as log10 Surface Area
        -m    1000
        
#Biotite
#	-parm % kbiots       %              # fit the first parameter as log10 Surface Area
#        -m    1000

Hornblende
	-parm % khnbds       %              # fit the first parameter as log10 Surface Area
        -m    1000

Kaolinite
	-parm  % kkaols       %        # fit the first parameter as log10 Surface Area
        -m    1000
        
#Pyrite_O2
#    -formula  Pyrite 1
#    -m        1000
#    -m0       1000
#    -parms    % kpyr_s       % 0 0.5 -0.11
#    -tol      1e-008	

Nitrification
    -formula Amm -1 NH3 +1
    -m       1000
    -parms   % knitrs       %


-step 86400

END

EQUILIBRIUM_PHASES 0 Equilibrate precipitation and irrigation with atmospheric pO2 and pCO2 @ 686 meters above mean sea level, in log ppm
	O2(g)         -0.713  1000      # log(ppO2)  ~ ppO2 = 0.193 @ 686 m above sea level; = 0.21 at sea level
	CO2(g)        -3.43   1000      # log(ppCO2) ~ ppCO2 = 369 ppm @ 686 m above sea level; ppCO2 400ppm at sea level


EQUILIBRIUM_PHASES 2 Equilibrate stream pO2 and pCO2, in log ppm
# CO2 should be greater than atmosphere. Stream O2 not sensitive in this model.
	O2(g)       % O2str        % 1000      # 
	CO2(g)      % CO2str       % 1000      # ~ PCO2 
END
EQUILIBRIUM_PHASES 30 Unsaturated zone
#	Kaolinite   0 0 #precipitate
	Goethite    0 0 #precipitate
	Gibbsite    0 0 #precipitate
#	Regolith_biotite %   ksBiou   % 0 #precipitate  # fit the first 0, equilibrium constant
	O2(g)       % O2uz         % 1000      # Less oxygen and more CO2 as a result of root respiration and pyrite oxidation
       CO2(g)       % CO2uz        % 1000      # UZ PCO2, greater than -3.65 which is atmospheric at 3500 m amsl
END
EQUILIBRIUM_PHASES 4 Saturated zone
#	Kaolinite   0 1e-2 # precipitate
	Goethite    0 1e-2 # precipitate
	Gibbsite    0 0 # precipitate
#	Regolith_biotite %   ksBiou   % 0 #precipitate  # fit the first 0, equilibrium constant
END
REACTION 1 Sea salt
C       2.180677103e-03 
Ca      1.066093459e-02 
Cl      5.665211193e-01 
K       1.057772844e-02 
Mg      5.506609125e-02 
Na      4.854108848e-01 
S       2.925759851e-02 
Si      7.382312412e-05 
H       2.255504347e-03 
O       1.243079188e-01 
%   dryfal       %   # 1 mol makes a liter of seawater 
END
RATES
Dryfall
-start
10 rate = parm(1) / 86400
20 moles = rate * time
30 save moles
-end
END
KINETICS 10
dryfall
-m 1000
-parm %   dryfal       % 
-formula \
C       2.180677103e-03 \
Ca      1.066093459e-02 \
Cl      5.665211193e-01 \
K       1.057772844e-02 \
Mg      5.506609125e-02 \
Na      4.854108848e-01 \
S       2.925759851e-02 \
Si      7.382312412e-05 \
H       2.255504347e-03 \
O       1.243079188e-01
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
