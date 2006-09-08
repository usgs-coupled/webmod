run_test1 = 
function()
{
    sel = "SELECTED_OUTPUT\
        -temperature true\
        -totals Ca\
        -molalities Ca+2\
        -activities Ca+2\
        -saturation_indices Calcite\
        -gases CO2(g)"

    solution.25 = "SOLUTION 25 Test solution number 25\
        temp      25\
        pH        7.0     charge\
        pe        4.5\
        redox     O(-2)/O(0)\
        units     ppm\
        density   1.02\
        Ca        80.\
        S(6)      96.     as SO4\
        S(-2)     1.      as S\
        N(5) N(3) 14.     as N\
        O(0)      8.0\
        C         61.0    as HCO3      CO2(g)     -3.5\
        Fe        55.     ug/kgs as Fe S(6)/S(-2) Pyrite\
        -isotope  13C     -12.   1.  # permil PDB\
        -isotope  34S     15.    1.5 # permil CDT\
        -water    0.5     # kg\
        solution 26"

    cat("Reading phreeqc.dat from current working dir ...\n")
    phrReadDB("phreeqc.dat")

    cat("Reading input instrings with phrReadString ...\n")
    phrReadString(sel)
    phrReadString(solution.25)

    cat("Running with phrRun ...\n")
    phrRun()

    cat("Returning SELECTED_OUTPUT ...\n")
    return(phrGetSelectedOutput())
}
