! **********************************************************************
! $Id: fmodules.inc 491 2008-10-24 23:13:50Z rsregan $
! system include file for all fortran modules
! **********************************************************************
!!      INTEGER, EXTERNAL :: decldim, declfix, declparam, declvar, declpri
!!      INTEGER, EXTERNAL :: declfix, declparam, declvar, declpri
      INTEGER, EXTERNAL :: declvar, declpri
      DOUBLE PRECISION, EXTERNAL :: deltim, delnex, djulian
      INTEGER, EXTERNAL :: getdim, getparam, getstep
      INTEGER, EXTERNAL :: getvar, readvar, putvar, declparam_u
      INTEGER, EXTERNAL :: julian, isleap
      INTEGER, EXTERNAL :: getdataname, getoutname, getoutdirfile
      INTEGER, EXTERNAL :: declmodule, control_string_, control_integer

!      EXTERNAL dattim, opstr, dpint4, dpstr, dpdble, dpreal
      EXTERNAL decl_control
! fmodules.inc
!      MODULE WEBMOD_MAX
! Dimensions
        INTEGER, PARAMETER :: MAXRAIN =      25
        INTEGER, PARAMETER :: MAXTEMP =      25
        INTEGER, PARAMETER :: MAXGND =       25
        INTEGER, PARAMETER :: MAXHUM =       25
        INTEGER, PARAMETER :: MAXPRES =      25
        INTEGER, PARAMETER :: MAXOBS =       25
        INTEGER, PARAMETER :: MAXSOL =       10
        INTEGER, PARAMETER :: MAXEVAP =      10
        INTEGER, PARAMETER :: MAXWIND =      10
        INTEGER, PARAMETER :: MAXSNODPL =    10
        INTEGER, PARAMETER :: MAXCHAN =     200
! Specific to WEBMOD
        INTEGER, PARAMETER :: MAXMRU =      500
        INTEGER, PARAMETER :: MAXSNOPIL =    10
        INTEGER, PARAMETER :: MAXTOPCHAN =   50
        INTEGER, PARAMETER :: MAXNAC =       50
        INTEGER, PARAMETER :: MAXXKBIN =      9
        INTEGER, PARAMETER :: MAXSOLUTE =    75 ! reduced to the number of species in phreeq_lut using nphq_lut
        INTEGER, PARAMETER :: MAXSOLOBS =    100
        INTEGER, PARAMETER :: MAXCONVERT =    3
        INTEGER, PARAMETER :: MAXRESINP =    21
        INTEGER, PARAMETER :: MAXIRRIG =    100
        INTEGER, PARAMETER :: MAXHYDRO =    100
        INTEGER, PARAMETER :: MAXSOLNSETS = MAXMRU
        INTEGER, PARAMETER :: MAXMRURES = 9
        INTEGER, PARAMETER :: MAXSOLN = 2*(((MAXMRURES+MAXNAC)*MAXMRU)+MAXHYDRO)
        INTEGER, PARAMETER :: MAXNMR_3D = MAXNAC*MAXMRU*MAXRESINP
        INTEGER, PARAMETER :: MAXCHEMOBS = 20
        INTEGER, PARAMETER :: MAXCHEMVAR = 10
        INTEGER, PARAMETER :: MAXENTITY = 11
!        INTEGER, PARAMETER :: MAXTALLYCOL = 100
        INTEGER, PARAMETER :: MAXRTDELAY =  100
        INTEGER, PARAMETER :: MAXLAG =        5
! do not change the following PARAMETERS
        INTEGER, PARAMETER :: MAXOBJFUNC = 5
        INTEGER, PARAMETER :: MAXFORM = 1
        INTEGER, PARAMETER :: MAXDAY = 366
        INTEGER, PARAMETER :: MAXMO = 12
        INTEGER, PARAMETER :: MAXLAPSE = 3
        INTEGER, PARAMETER :: MAXDEPLVAL = MAXSNODPL*11
! Constants
!       REAL, PARAMETER :: NEARZERO = 1.0E-15, SMALL = 0.000001
!       REAL, PARAMETER :: ONETHIRD = 1.0/3.0, TWOTHIRDS = 2.0/3.0
!       REAL, PARAMETER :: FOURTHIRDS = 4.0/3.0, FIVETHIRDS = 5.0/3.0
!       REAL, PARAMETER :: FIVE_NINTHS = 5.0/9.0, NINE_FIFTHS = 9.0/5.0
!       REAL, PARAMETER :: FEET2METERS = 0.304800, CFS2CMS = 0.028316844
!       REAL, PARAMETER :: INCH2CM = 2.54, INCH2MM = 25.4
!       REAL, PARAMETER :: TSTEPCHK = 23.999, CONVRG = 0.001
!       DOUBLE PRECISION, PARAMETER :: DTSTEPCHK = 23.999D0
!      END MODULE WEBMOD_MAX
