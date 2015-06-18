#include "defines.h"
c **********************************************************************
c     obs_chem.f: read concentrations for atmospheric deposition,
c                 external chemical sources, and observed chemistry.
c    27 apr  2010 - Port to Fortran 90 with module and dynamic memory
c **********************************************************************
      MODULE WEBMOD_OBSCHEM
      IMPLICIT NONE
#include "fmodules.inc"
c   Dimensions and Local Variables
c nsolobs is the number of sites or irrigation sources with defined solute
c concentrations returned from parsechemobs that reads the content of the chemdat file.
c
c nsolute is the number of solutes to be modeled
c
      integer, save :: Nsolobs, Nconvert, Nsolute, Nchem_ext, Nphq_lut
      integer, save :: Nchemdep_obs, Nchemdat_obs, Nchemobs, n_iso
!   Declared Parameters
      integer, save :: chem_sim, chem_ext, ppt_chem
      double precision, save, allocatable:: convfactor(:)
!   Declared Variables
c
c chemdat contains the daily values of concentrations, pH, and Temp. Concentrations
c   are converted to Molal and stored temporarily in cconc before being assigned
c   to cconc_precipM, cconc_extM, or cconc_obsM.
c
      double precision, allocatable :: chemdat(:,:)
      double precision, allocatable :: cconc(:)
c
c       cconc_precipM, and cconc_extM will be populated
c       in the order of the dimensions names for nsolute
c       declared in the phreeq_mms module. Units are converted
c       to moles per liter
      double precision, save ::c_precip_pH,c_precipT
      double precision, save, allocatable:: cconc_precipM(:)
      double precision, save, allocatable:: cconc_extM(:,:)
      double precision, save, allocatable:: cconc_obsM(:,:)
      double precision, save, allocatable:: c_ext_pH(:), c_extT(:)
      double precision, save, allocatable:: c_obs_pH(:), c_obsT(:)
c internal or private variables
      logical, save:: step1
      integer, save :: chemdat_flag
!
! TYPE structures
!
!  phq_lut is read from the phreeq_lut file that is transcribed with
!  minor modifications from the phreeq.dat and iso.dat master species list.

      TYPE :: lut   ! phreeq lookup table with conversions for solutes and isotopes
         character(20) :: Spec_ID   ! ID, Species, 'As'
         character(20) :: Spec   ! ID, Species, 'As'
         character(20) :: Spec_as   ! ID, Species, 'As'
         double precision :: M2mg   ! conversions from Molal to mg/L
         double precision :: M2meq   ! conversions from Molal to meq/L
         character (6) :: isounit    ! units for isotopes
         double precision :: isoratio  ! ratio of isotope in standard
         character(40) :: isofrac_gl   ! expression for common log of alpha for equilibrium fractionation for gas from liquid at temperature of reaction
         character(40) :: isofrac_gs   ! expression for common log of alpha for equilibrium fractionation for gas from solid at temperature of reaction
         double precision :: isofrac_max   ! maximum diffusive fractionation in permil (0% rel hum, n = 1, theta=1, Craig and Gordon (1965) )
      END TYPE lut

     
      TYPE(lut), save, allocatable :: phq_lut(:)
!
! sol_name is the name of each of nsolutes. They must match one of those
! listed in the phreeq.lut file
!
      character(len=12), save, allocatable :: sol_name(:)
!
!  for each of nolutes sol_id describes the row number in phq_lut,
!  sol_id(i)%phq, the phreeq tally table created for each mix,
!  sol_id(i)%tally, and whether it is an isotope, sol_id(i)%iso
!  (phq_lut is read from the phreeq_lut file that is transcribed with
!  minor modifications from the phreeq.dat master species list).
!
!
      TYPE :: sols   ! row number in phq_lut and tally table for each of nsolutes
         integer :: phq   ! Row number of solute in phreeq_lut file
         integer :: phq_e ! Row number of uncharged element or species with no valence state (i.e. Ca) in phreeq_lut file
                          ! (i.e. S when solute is S(6), or N when solute is N(6) - populated in phreeq_mms).
         integer :: tally ! Row number of solute in tally table (populated in phreeq_mms)
         integer :: tally_e ! Row number for uncharged element in in tally table
         logical :: iso   ! Is this species an isotope (True|False)
      END TYPE sols
      
      TYPE(sols), save, allocatable :: sol_id(:)
c
c interface for unit_ID
c
      INTERFACE 
         integer function unit_id(STRING, unit_lut)
           IMPLICIT NONE
           CHARACTER(len=12), INTENT(IN) :: STRING, unit_lut(:)
         END function unit_id
       END INTERFACE
c
c The obs_sol_lut matrix is allocated in parsechemobs and made available elsewhere
c
      INTERFACE 
         subroutine parsechemobs(buf,unit_lut, nsolute,
     &                nsolobs, validobs, tph)
           CHARACTER(len=1000), INTENT(IN) :: buf
           CHARACTER(len=12), INTENT(IN) :: unit_lut(:)
           INTEGER, INTENT(IN) :: nsolute
           INTEGER, INTENT(OUT) :: nsolobs
           INTEGER, INTENT(OUT) :: validobs(:)
           LOGICAL, INTENT(OUT) :: tph
         END subroutine parsechemobs
       END INTERFACE
c     subroutine parsechemobs(buf,unit_lut, nsolute,
c                              nsolobs, validobs, obs_sol_lut,tph)
c
c     Inputs
c     buf - the line in the .chemdat file that describes
c           the species and concentration units
c     phq_lut - the lookup list of species used in phreeqc 
c               (shared via webmod_phreeq_mms module)
c     unit_lut - the lookup list of acceptable concentration units
c     nsolute - the number of solutes in the current model run
c     sol_id - the row number in phq_lut (%phq) and the tally table (%tally)
c              of the nsolutes(shared via webmod_phreeq_mms module)
c
c     Ouputs
c     nsolobs - the number of solutes in the .chemdat file
c     validobs - index in order of nsolute indicating which
c           observation field corresponds to a valid model solute.
c           Set to zero if no corresponding observation field in the
c           chemdat file.
c     obs_sol_lut - a matrix with dimension (nsolobs,2).
c           The two dimensions include integer codes for the species ID,
c           and the concentration unit.
c     tph - logical test to confirm that tempc (temperature in celsius)
c           and pH are the first data fields for each set of observations
c **********************************************************************
c
c for each observed solute, obs_sol_lut contains the element ID (from 
c phq_lut)and an integer describing the concentration units, and
c the dimension index of the solute for each column of observed solute
c concentrations.
c
c valid obs maps the solutes to be modeled to their position in
c  the observed data (chemdat) file
c
c iso_list lists which if any of the nsolutes is an isotope
c
      integer, save, allocatable:: obs_sol_lut(:,:), validobs(:),
     $  iso_list(:)
c
c The unit codes:
c 1) mg/L;
c 2) meq/L;
c 3-5) uconc1, uconc2, and uconc3, user-defined conversion factors to
c    be read as parameter convfactor(nconvert(=3)).
c    If a uconc is indicated for a species in the chemdat file, then
c    the associated parameter value (convfactor 1, 2, or 3) should have
c    the value to convert the user units to moles per liter.
c 6) mmol/L; 
c 7) mol/L; or
c 8) permil
c
c unit_lut describes the eight acceptable units described below
      character(LEN=12), parameter :: unit_lut(8) = 
     $   (/'mg/L        ','meq/L       ','uconc1      ',
     $     'uconc2      ','uconc3      ','mmol/L      ',
     $     'mol/L       ','permil      '/)
c
c M2mM is equal to 1000 to convert moles to millimoles.
      double precision, parameter:: M2mM = 1000.0
c
      END MODULE WEBMOD_OBSCHEM

c *************************************************************************
c
c     main obs_chem routine
c
      integer function obs_chem(arg)
      IMPLICIT NONE
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
      CHARACTER(len=256) SVN_ID
c
! Functions
      INTEGER, EXTERNAL :: obsc_decl, obsc_init, obsc_run, obsc_clean
      character(len=12), external :: parse
c
      SVN_ID = 
     $     '$Id: obs_chem.f 35 2010-03-30 17:34:40Z rmwebb $ '

      obs_chem = 0

      if (arg.eq.'declare') then
         obs_chem = obsc_decl()
      else if (arg.eq.'initialize') then
         obs_chem = obsc_init()
      else if (arg.eq.'run') then
         obs_chem = obsc_run()
      else if (arg.eq.'cleanup') then
         obs_chem = obsc_clean()
      end if

      if(obs_chem.eq.1) close (unit = 15)

      return
      end
c **********************************************************************
c
c     obsc_decl - makes public variable declarations for the
c                     obs_chem module
c
 
      integer function obsc_decl()

      USE WEBMOD_OBSCHEM
      implicit none
      obsc_decl = 1
!
! Get dimensions
!
      nsolute = getdim('nsolute')
      nphq_lut = getdim('nphq_lut')
      nchem_ext = getdim('nchem_ext')
      nchemobs = getdim('nchemobs')
      nconvert = getdim('nconvert')

c
c chemdat_flag indicates whether a chemdat file could be found.
c If so then concentration values for precip, external sources,
c and quality of water samples are read from there every time step.
c If the file cannot be found then constant concentrations
c described in the pqi file and assigned in the phreeq_mms
c module will be used.
c

      if(declvar('obs_chem', 'chemdat_flag', 'one', 1,
     +     'integer','Equals 1 if the chemdat file was found '//
     $     'and 0 if not','none', chemdat_flag).ne.0) return
c
c  The following six variables list the temperature, in degrees celsius,
c  and the pH for the precipitation, the external sources,
c  and the water quality samples.
c
      if(declvar('obs_chem', 'c_precipT', 'one', 1,
     +     'double','Temperature of NADP precipitation',
     +     'deg celsius', c_precipT).ne.0) return

      if(declvar('obs_chem', 'c_precip_pH', 'one', 1,
     +     'double','pH of NADP precipitation',
     +     'pH units', c_precip_pH).ne.0) return

      ALLOCATE (c_extT(Nchem_ext))
      if(declvar('obs_chem', 'c_extT', 'nchem_ext', Nchem_ext,
     +     'double','Temperature of external source',
     +     'deg celsius', c_extT).ne.0) return

      ALLOCATE (c_ext_pH(Nchem_ext))
      if(declvar('obs_chem', 'c_ext_pH', 'nchem_ext', Nchem_ext,
     +     'double','pH of external source',
     +     'pH units', c_ext_pH).ne.0) return

      ALLOCATE (c_obsT(Nchemobs))
      if(declvar('obs_chem', 'c_obsT', 'nchemobs', Nchemobs,
     +     'double','Temperature of water sample',
     +     'deg celsius', c_obsT).ne.0) return

      ALLOCATE (c_obs_pH(Nchemobs))
      if(declvar('obs_chem', 'c_obs_pH', 'nchemobs', Nchemobs,
     +     'double','pH of water sample',
     +     'pH units', c_obs_pH).ne.0) return
c
c
c  cconc_precipM is an nsolute dimensioned matrix of observed solute
c  concentrations in moles/L (assumed equivalent to moles/kg of water).
c  The matrix is indexed by the nsolute dimension names. The cconc_precipM
c  values on any time step will be used to describe concentrations for any
c  and all mrus even if they are indexed to different meteorological
c  observation stations.
c
      ALLOCATE (cconc_precipM(Nsolute))
      if(declvar('obs_chem', 'cconc_precipM', 'nsolute', Nsolute,
     +     'double','Solute concentrations in precipitation',
     +     'moles/liter', cconc_precipM).ne.0) return
c
c  cconc_extM is an similar matrix with the additional dimension of
c  nchem_ext, or the number of unique external sources.
c
      ALLOCATE (cconc_extM(nchem_ext,nsolute))
      if(declvar('obs_chem', 'cconc_extM', 'nchem_ext,nsolute',
     $     nchem_ext*nsolute,'double','Solute concentrations '//
     $     'describing an external source',
     +     'moles/liter', cconc_extM).ne.0) return
c
c  cconc_obsM is an similar matrix with the dimension of
c  nobs, or the number of unique chemical sampling sites.
c
      ALLOCATE (cconc_obsM(nchemobs,nsolute))
      if(declvar('obs_chem', 'cconc_obsM', 'nchemobs,nsolute',
     $     nchemobs*nsolute,'double','Solute concentrations '//
     $     'for a water sample',
     +     'moles/liter', cconc_obsM).ne.0) return

c
c  nsolobs is the number of solutes included in the input chemical
c  deposition file (.chemdat).
c
      if(declpri('obs_chem_nsolobs', 1, 'integer',
     +     nsolobs) .ne. 0) return 
c
c convfactor is available for up to 3 user-specified conversion factors
c to be used to convert input chemical concentration identified in
c the chemdat file as uconc1, uconc2, or uconc3.
c
      ALLOCATE (convfactor(nconvert))
      if(declparam('obs_chem', 'convfactor', 'nconvert', 'double',
     +     '1', '0.000001', '1000000',
     +     'User-defined conversion factor for solute '//
     $     'inputs or outputs',
     +     'User-defined conversion factors. Units are converted as '//
     $     'follows: '//
     $     'uconc(in)/convfactor = mol/L; or '//
     $     'mol/L*convfactor = uconc(out).',
     +     'user-defined') .ne.0) return
c
c A value of '0' for ppt_chem indicates that no time series of precipitation
c chemistry is available and therefore constant concentrations as defined by 
c init_soln_ppt will be used. If ppt_chem = 1 then a time series of precip
c chemistry will be read from the chemdat file. Whether varying or constant,
c all deposition and irrigation will be at equilibrium with atmospheric oxygen
c and carbon dioxide as described by the parameter atmos_eq_ph which points
c to a EQUILIBRIUM_PHASES block in the pqi file.
c
      if(declparam('obs_chem', 'ppt_chem', 'one', 'integer',
     +     '0', '0', '1',
     +     '0 - Precipitation chemistry is constant as defined by '//
     $     'init_soln_ppt and the pqi file; 1 - Precipitation '//
     $     'chemistry is read every time step from the chemdep file',
     +     '0 - Precipitation chemistry is constant as defined by '//
     $     'init_soln_ppt and the pqi file; 1 - Precipitation '//
     $     'chemistry is read every time step from the chemdep file',
     +     'none') .ne.0) return
c
c A value of '0' for chem_ext indicates that no time series of external
c chemistry is available for the nchem_ext unique chemical sources. Therefore
c constant concentrations as defined by init_soln_ext will be used.
c If chem_ext = 1 then a time series of solution chemistry for each of
c the nchem_ext sources will be read from the chemdat file.
c
      if(declparam('obs_chem', 'chem_ext', 'one', 'integer',
     +     '0', '0', '1',
     +     '0 - External chemistry is constant as '//
     $     'defined by init_soln_ext and the pqi file; '//
     $     '1 - External chemistry is provided in '//
     +     'the chemdep file for every time.',
     +     '0 - External chemistry is constant as '//
     $     'defined by init_soln_ext and the pqi file; '//
     $     '1 - External chemistry is provided in '//
     +     'the chemdep file for every time.',
     +     'none') .ne.0) return

      if(declparam('phreeqmms', 'chem_sim', 'one', 'integer',
     +   '1', '0', '1',
     +   'Simulate solute fluxes (0=no; 1=yes)',
     +   'Simulate solute fluxes (0=no; 1=yes)',
     +   'integer').ne.0) return

      allocate (iso_list(nsolute))
      allocate (validobs(nsolute))
      allocate (sol_name(nsolute))
      allocate (sol_id(nsolute))
      allocate (phq_lut(nphq_lut))

      obsc_decl = 0

      return
      end
c****************************************************************************
c
c     obsc_init - initializes obs_chem module. The species and units listed
c                 in the input chemistry file are indexed by lookup table
c                 phq_lut, and identified by their nsolute index as defined
c                 by the active parameters in the model.
c

      integer function obsc_init()

c$$$#if defined(_WIN32)
c$$$      USE IFPORT
c$$$#endif
      USE WEBMOD_OBSCHEM
      IMPLICIT NONE
C      intrinsic getenv
      logical filflg,tph
      integer, external :: phreeq_id
c  internal or private variables
      character(len=12) solname
      character(len=12), external :: parse
      integer i, j, ret, io
      integer, external :: length
      character(len=80) test
      !character(len=135) data_path, phreeq_lut_file, ptxt
!      character(len=135) data_path, phreeq_lut_file
      character(len=135) chemdat_file, phreeq_lut_file
      character(len=1000) buf

      obsc_init = 1


      if(getparam('phreeqmms', 'chem_sim', 1,
     $     'integer',chem_sim) .ne.0) return

c
c Skip initialization (and later run) if no chemical simulations desired
c
      if(chem_sim.eq.1) then
c

      step1 = .true.
      chemdat_flag = 1

      if(getparam('obs_chem', 'convfactor', nconvert, 'double',
     $     convfactor) .ne.0) return

      if(getparam('obs_chem', 'ppt_chem', 1, 'integer',
     $     ppt_chem) .ne.0) return

      if(getparam('obs_chem', 'chem_ext', 1, 'integer',
     $     chem_ext) .ne.0) return

      if(ppt_chem.ne.0.and.ppt_chem.ne.1) then
         print*,'The parameter ppt_chem must equal 0 or 1.'
         print*,'Run terminated'
         return
      end if

      if(chem_ext.ne.0.and.chem_ext.ne.1) then
         print*,'The parameter chem_ext must equal 0 or 1.'
         print*,'The number of unique external chemical sources '//
     $        'is described by the dimension nchem_ext.'
         print*,'Run terminated'
         return
      end if

      nchemdep_obs = nchem_ext*chem_ext+ppt_chem
      nchemdat_obs = nchemdep_obs + nchemobs

! Read the look up tables for species and conversion units from the
! phreeq_lut file.  The file is used by both the obs_chem module and the
! phreeq_mms module.
!
! 10 elements in each row of phreeq_lut file:
! species_index, name, 'as', mg/mole, meq/mole, isounit, isostd (ratio in standard), 
!     PHREEQ named expression for log(alpha)(vapor-liq), log(alpha)(vapor-solid), and maximum diffusive isotope depletion (in permil)
! examples
! species_26 Alkalinity CO3-2 60000 2000 none 1 none none 1
! species_27 S SO4-2 96064 2000 none 1 none none 1
! species_51 D D2O 20000 0 permil 0.00015576 'Log_alpha_D_H2O(g)/H2O(l)' 'Log_alpha_D_H2O(g)/H2O(s)'  -25.1
! species_57 [18O] H2[18O] 20000 0 permil 0.0020052 'Log_alpha_18O_H2O(g)/H2O(l)' 'Log_alpha_18O_H2O(g)/H2O(s)' -28.5
! species_63 [14C] [14C]O2 46000 0 pmc 1.175887709E-12 none none 1
!
! The flag for an isotope is an isounit other than "none". Acceptable
!  isounits: permil deviation from standard, permil; tritium units, TU; 
!            and percent modern carbon, pmc.
! PHREEQC converts istope units to Molal units so no coversion is necessary when 
!  reading these units from the obs chem data file. However, the inverse computation
!  from Molal to isounit requires the ratio in the standard, isostd, and the species
!  in the denominator of the ratio, isoden (total O, C, or S for example)
! 
!
!      call getenv('mms_data_dir',phreeq_lut_file)
!      phreeq_lut_file= phreeq_lut_file(1:length(phreeq_lut_file))//
!     $     '/phreeq_lut'
      IF ( control_string(phreeq_lut_file, 'input_dir').NE.0 ) RETURN
      i = index(phreeq_lut_file,CHAR(0))   ! CHAR(0) is end of strings returned from control_string call
      phreeq_lut_file= phreeq_lut_file(1:i-1)//'phreeq_lut'
      inquire(file=phreeq_lut_file,exist=filflg)
      if (filflg) then
         open (14, file=phreeq_lut_file, iostat=io)
!         read (14,*)(spec_id(i),(sol_lut(i,j),j=1,2),
!     $        (conv2moles(i,k),k=1,2),i=1,nphq_lut)
         read (14,*)(phq_lut(i),i=1,nphq_lut)
      close (unit=14)
      else
         print*,'The required look up table of phreeq species ',
     $        phreeq_lut_file,' could not be found. Run terminated.'
         return
      end if
c
c  This section retrieves the nsolute dimension names, and identifies
c  the species using the phq_lut lookup table. The number and species 
c  of isotopes are also identified (isounits not equal to 'non').
c
c  Get the dimension names that describe the solute species,
c  match them with their species ID in phq_lut and save the IDs
c  in the order of their nsolute indexes in the variable sol_id(i)%phq.
c
c  According to section 6.8.9 of the MMS manual
c  http://wwwbrr.cr.usgs.gov/mms/html/chap5.frame.new.html#HDR59
c
c$$$  getdimname (name, i, index_name)
c$$$  name - dimension name.
c$$$  i - "ith" index specifier.
c$$$  index_name - returned name of the "ith" index.
c
c First pad with dimension names with blanks. Then overwrite with
c the names listed in dimension names list. Finally initialize the
c validobs array.
c
      n_iso = 0  ! number of isotopes
      do 100 i = 1, nsolute
         solname='            '
         sol_name(i)='            '
         call getdimname('nsolute',i,solname)
c solname is terminated by mms with a null so take one less character
c and overwrite into sol_name
         sol_name(i)=solname(1:length(solname)-1)
         sol_id(i)%phq = phreeq_id(sol_name(i))
! assign elemental root if valenced species
         sol_id(i)%phq_e = phreeq_id(parse(sol_name(i), '('))
         sol_id(i)%iso=.FALSE.
         if(phq_lut(sol_id(i)%phq)%isounit.eq.'permil') then
             sol_id(i)%iso=.TRUE.
             n_iso = n_iso+1
             iso_list(n_iso)=i
         end if
 100  continue

c
c Read the chemical deposition file through to the line that describes
c the species and units in the deposition file
c
C kludge line to use only one chemdat file, comment out for normal run
c      data_path = './input/loch3.dat.chemdat'
c end kludge
      IF (control_string(chemdat_file, 'chemdat_file').NE.0 ) RETURN
      i = index(chemdat_file,CHAR(0))   ! CHAR(0) is end of strings returned from control_string call
      chemdat_file= chemdat_file(1:i-1)
      inquire(file=chemdat_file,exist=filflg)
      if (filflg.and.(ppt_chem.ne.0)) then
         open (15, file=chemdat_file)
 10      read (15,1000,END=500) test
         if (test .ne. 
     $        '** Species and Units **') goto 10
 500     continue
 1000    format(A23)
         read (15,1100)buf
 1100    format(A1000)

         call parsechemobs(buf,unit_lut,nsolute,
     $        nsolobs,validobs,tph)
c     
c     Notify user and stop the run if the following errors appear:
c     - if tempc and pH are not the first two fields
c     - if a species name in the chemdat file does not match a phreeqc
c     element name as listed in the phq_lut table
c     - if one of the units listed in the chemdat file does not match
c     one of the eight acceptable units listed in unit_lut.
c     - if meq/L was listed in the chemdat file for a non-charged species
c     
         if(.not.tph) then
            print*,"Temperature, in celsius, and pH must be precede "//
     $           "concentrations in the chemdat file. "//
     $           "Run stopped."
            return
         end if
c
c allocate size of chemdat and cconc for use in the run section
c
      ALLOCATE (chemdat(nchemdat_obs,nsolobs+2))
      ALLOCATE (cconc(nsolute))
c
         do 1200 i = 1, nsolobs
            if (obs_sol_lut(i,1).eq.-99) then
               print*,"Unrecognized species name in chemdat file. "//
     $              "Run stopped."
               return
            end if
            if(obs_sol_lut(i,2).eq.-99) then
               print*,"Unrecognized concentration unit in chemdat "//
     $              "file. Run stopped"
               return
            end if
            if (phq_lut(obs_sol_lut(i,1))%M2meq.lt.1.0
     $           .and. obs_sol_lut(i,2).eq.2) then
               print*,"The chemdat file indicates equivalence units "//
     $              "for a non-charged species. Run stopped."
               return
            end if
 1200    continue
         do 1300 i = 1, nsolute
            if (sol_id(i)%phq.eq.-99) then
               print*,"One of nsolute names not recognized. Verify "//
     $           "that it is described in the phreeq.lut file. "//
     $           "Run stopped."
               return
            end if
 1300    continue
c     
c     Read 3 descriptive lines preceding data
c     
         DO 101 I=1,3
            read (15,1100)
 101     CONTINUE
c     
c     
c     
      else

         print*,'The file ',chemdat_file,' describing time series of '//
     $      'input and/or observed concentrations could not be '//
     $      'found or inputs of constant input concentrations have ' //
     $      'been selected from the pqi file.'
         chemdat_flag = 0
      end if

c
c Initialize temperature and pH to standard values. temperature is in degrees C
c
      c_precipT = 25.0
      c_precip_pH = 7.0
      
      if(nchem_ext.gt.0) then
         do 8 i=1,nchem_ext
            c_extT(i)=25.0
            c_ext_pH(i)=7.0
 8       end do
      end if
      if(nchemobs.gt.0) then
         do 9 j=1,nchemobs
            c_obsT(j)=25.0
            c_obs_pH(j)=7.0
 9       end do
      end if
c
c Initialize concentrations with zeros. The zeros will be replaced
c with observations if present in the chemobs file
c
      do 20 i=1,nsolute
         cconc_precipM(i)=0.0
         if(nchem_ext.gt.0) then
            do 11 j=1,nchem_ext
               cconc_extM(j,i)=0.0
 11         continue
         end if
         if(nchemobs.gt.0) then
            do 13 j=1,nchemobs
               cconc_obsM(j,i)=0.0
 13         continue
         end if
 20      continue

      else
         print*,'Chem_sim parameter equals zero so no chemical '//
     $        'inputs read or processed'
      end if

      obsc_init = 0
 
      return
      end

c **********************************************************************
c
c     obsc_run - runs obs_chem module
c
 
      integer function obsc_run()

      use WEBMOD_OBSCHEM
      IMPLICIT NONE

      logical day1run
      integer datetime(6), obsc_datetime(6)
      integer i, j, k, iv,iobs, io

      obsc_run = 1
c
c If no chemical simulation (chem_dep=0) do nothing
c
      if(chem_sim.eq.1) then

c
c If chemdat_flag equals zero (no chemdat file) then there is
c nothing to do in this module. Default concentrations described
c in the pqi file and assigned in phreeq_mms will be used.
c

      if(chemdat_flag.eq.1) then

         nchem_ext = nchemdep_obs-ppt_chem
c
c Reset temperature and pH to standard values. temperature is in degrees C
c
         if(ppt_chem.eq.1) then
            c_precipT = 25.0
            c_precip_pH = 7.0
         end if
         if(chem_ext.eq.1) then
            do 8 i=1,nchem_ext
               c_extT(i)=25.0
               c_ext_pH(i)=7.0
 8          end do
         end if
         if(nchemobs.gt.0) then
            do 9 j=1,nchemobs
               c_obsT(j)=25.0
               c_obs_pH(j)=7.0
 9          end do
         end if

c
c Reset concentrations to zero. The zeros will be replaced
c with observations if present in the chemobs file
c
      do 10 i=1,nsolute
         cconc_precipM(i)=0.0
         if(chem_ext.gt.0) then
            do 11 j=1,nchem_ext
               cconc_extM(j,i)=0.0
 11         continue
         end if
         if(nchemobs.gt.0) then
            do 13 j=1,nchemobs
               cconc_obsM(j,i)=0.0
 13         continue
         end if
 10   continue
c     
c Read in the date/time and solute concentrations from
c the chemdat file
c
c Scroll ahead in the data file until the first matching date
      call dattim('now', datetime)
      read(15,*,iostat=io)(obsc_datetime(i),i=1,6),
     $     ((chemdat(iobs,j),j=1,(nsolobs+2)),iobs=1,nchemdat_obs)
      day1run = .false.
      if (step1) then
         do 12 while (.not.day1run)
            k=0
            do 14 i=1,6
               if(datetime(i).eq.obsc_datetime(i)) k=k+1
 14         continue
            if (k.eq.6) then
               day1run = .true.
               step1 = .false.
            else
               read(15,*)(obsc_datetime(i),i=1,6),
     $         ((chemdat(iobs,j),j=1,nsolobs+2),iobs=1,nchemdat_obs)
            end if
 12      continue
      end if

cccccccccccccccccccc
c
c Check that the date/time in the chemdat file equal those read
c from the data file
c 
c$$$      call dattim('now', datetime)
      do 60 i=1,6
         if(datetime(i).ne.obsc_datetime(i)) then
            print*,"The date/times in the chemical observation "//
     $           "file are not synchronized with those in the "//
     $           "data file. Run aborted."
            return
         end if
 60   continue
c
c Read tempc and pH for precip, external sources and water quality samples
c
      if(ppt_chem.eq.1) then
         c_precipT = chemdat(1,1)
         c_precip_pH = chemdat(1,2)
      end if
      if(chem_ext.eq.1) then
         do 80 i=1,nchem_ext
            c_extT(i)=chemdat(i+ppt_chem,1)
            c_ext_pH(i)=chemdat(i+ppt_chem,2)
 80      end do
      end if
      if(nchemobs.gt.0) then
         do 90 j=1,nchemobs
            c_obsT(j)=chemdat(nchemdep_obs+j,1)
            c_obs_pH(j)=chemdat(nchemdep_obs+j,2)
 90      end do
      end if
      

      do 100 iobs = 1,nchemdat_obs

c
c Convert the observations to moles/liter and place in
c cconc_precipM, cconc_extM, or cconc_obsM.
c
         do 100 i=1,nsolute
            if(validobs(i).ne.0) then
               iv = validobs(i)+2 ! 1st two fields are tempc and pH
               if(chemdat(iobs,iv).ge.0) then   ! values outside of the range 0-8 are caught in the obs_sol_lut assignments
                  if (obs_sol_lut(iv-2,2).eq.1) then
                   cconc(i)=chemdat(iobs,iv)/
     $              phq_lut(obs_sol_lut(iv-2,1))%M2mg
!     $              conv2moles(obs_sol_lut(iv-2,1),obs_sol_lut(iv-2,2))
                  else if (obs_sol_lut(iv-2,2).eq.2) then
                   cconc(i)=chemdat(iobs,iv)/
     $              phq_lut(obs_sol_lut(iv-2,1))%M2meq
!     $              conv2moles(obs_sol_lut(iv-2,1),obs_sol_lut(iv-2,2))
                  else if (obs_sol_lut(iv-2,2).le.5) then
                     cconc(i)=chemdat(iobs,iv)/
     $                convfactor(obs_sol_lut(iv-2,2)-2)
                  else if (obs_sol_lut(iv-2,2).eq.6) then
                     cconc(i)=chemdat(iobs,iv)/M2mM
                  else
                     cconc(i)=chemdat(iobs,iv)    ! This include moles per liter (unit=7) and per mil for isotopes (unit = 8)
                  end if
               else
                  cconc(i) = chemdat(iobs,iv)
               end if
               if(ppt_chem.eq.1.and.iobs.eq.1) then
                  cconc_precipM(i)=cconc(i)
               else if (iobs.le.nchemdep_obs) then
                  cconc_extM(iobs-ppt_chem,i)=cconc(i)
               else
                  cconc_obsM(iobs-nchemdep_obs,i)=cconc(i)
               end if
            end if
  100     continue
 
       end if
 950   if(io.lt.0) then
          print*,'Chemdat file terminated before the end '//
     $         'of hydrologic data. Run terminated.'
          return
       end if
c
c  Landing point if no chemical simualtions desired
c
       end if
       obsc_run = 0
 
       return
       end
c***********************************************************************
c
c     obs_chemclean - Close the chemout file
c

      integer function obsc_clean()
      IMPLICIT NONE

#include "fmodules.inc"


      obsc_clean = 1

      close (unit=15)

      obsc_clean = 0

      return
      end



c **********************************************************************
c
c     subroutine parsechemobs(buf,unit_lut, nsolute,
c                              nsolobs, validobs, obs_sol_lut,tph)
c
c     Inputs
c     buf - the line in the .chemdat file that describes
c           the species and concentration units
c     phq_lut - the lookup list of species used in phreeqc 
c               (shared via webmod_phreeq_mms module)
c     unit_lut - the lookup list of acceptable concentration units
c     nsolute - the number of solutes in the current model run
c     sol_id - the row number in phq_lut (%phq) and the tally table (%tally)
c              of the nsolutes(shared via webmod_phreeq_mms module)
c
c     Ouputs
c     nsolobs - the number of solutes in the .chemdat file
c     validobs - index in order of nsolute indicating which
c           observation field corresponds to a valid model solute.
c           Set to zero if no corresponding observation field in the
c           chemdat file.
c     obs_sol_lut - a matrix with dimension (nsolobs,2).
c           The two dimensions include integer codes for the species ID,
c           and the concentration unit.
c     tph - logical test to confirm that tempc (temperature in celsius)
c           and pH are the first data fields for each set of observations
c **********************************************************************
c
      subroutine parsechemobs(buf,unit_lut,nsolute,
     $                        nsolobs,validobs,tph)
     
      USE WEBMOD_OBSCHEM, ONLY: phq_lut, sol_id, unit_id, obs_sol_lut
      IMPLICIT NONE

#include "fmodules.inc"
c
c functions
c
      integer, external :: length, phreeq_id
      character(len=12), external :: parse
c variables
      character(len=1000) buf, rembuf
      character(len=12) elem(MAXSOLOBS), elem_units(MAXSOLOBS)
      character(len=12) unit_lut(:), tempc, ph
!      character*12 solname, sol_name(MAXSOLUTE)
!      integer obs_sol_lut(:,:)
!      integer obs_sol_lut(MAXSOLUTE,2), sol_id(MAXSOLUTE)
      integer nsolute, nsolobs, validobs(:)
      integer elem_len(MAXSOLOBS), elem_ulen(MAXSOLOBS)
      integer tot_len, buf_len, i, j
      logical matched, tph
c$$$
c$$$      print*,"Starting Up. Enter the elements and units ",
c$$$     $     "separated by spaces"
      nsolobs=0
c     read(*,'(A100)') buf
      buf_len = length(buf)
      tot_len = 0
      rembuf = buf(1:buf_len)
c
c The first two fields must be tempc and pH
c
      tph = .true.
      
      tempc = parse(rembuf, ' ')
      elem_len(1) = length(tempc)
      tot_len = tot_len + elem_len(1) +1
      rembuf= buf(tot_len+1:buf_len)
      ph = parse(rembuf, ' ')
      elem_len(1) = length(ph)
      tot_len = tot_len + elem_len(1) +1
      rembuf= buf(tot_len+1:buf_len)
      if(tempc.ne.'tempc'.or.ph.ne.'pH') tph=.false.
c
c  The 60 loop below parses the input solute descriptions into
c  species, elem(nsolobs), and concentration units, elem_units(nsolobs).
c  This routine was written so that the user does not have to
c  enter the number of observed solutes either in the data file
c  or as an additional parameter in the model.
c
c  matched is a flag to ensure that each species is accompanied
c  by a concentration unit.
c
      matched = .TRUE.

       do 60 while (tot_len.lt.buf_len)
         if(matched) then
            nsolobs = nsolobs + 1
            elem(nsolobs) = parse(rembuf, ' ')
            elem_len(nsolobs) = length(elem(nsolobs))
            tot_len = tot_len + elem_len(nsolobs) +1
            matched = .false.
         else
c            elem_units(nsolobs) = '            '
            elem_units(nsolobs) =  parse(rembuf, ' ')
            elem_ulen(nsolobs) = length(elem_units(nsolobs))
            tot_len = tot_len + elem_ulen(nsolobs) +1
            matched = .true.
         end if
         rembuf= buf(tot_len+1:buf_len)
 60   continue
      if(.not.matched)print*,"Problem: There is not a unit for each ",
     $     "element."
!
! Allocate the chem obs matrix
!
c
c  obs_sol_lut is a 2-column matrix of the species number in phq_lut
c  and the units used for that species (mg/L, meq/L, permil, etc)
c  for each of the species listed in the  
c
      ALLOCATE (obs_sol_lut(nsolobs,2))
c
c zero the validobs
c
      do 100 i = 1, nsolute
         validobs(i) = 0
 100  continue
c
c Populate the obs_sol_lut matrix
c
      do 200 i = 1, nsolobs
c
c Look up the phreeq_id and the unit_id
c
         obs_sol_lut(i,1)= phreeq_id(elem(i))
         obs_sol_lut(i,2) = unit_id(elem_units(i),unit_lut)
c
c Populate validobs with the field number of the observation corresponding
c to the nsolute dimension
c
         do 250 j = 1, nsolute
            if (obs_sol_lut(i,1).eq.sol_id(j)%phq) then
               validobs(j) = i
            end if
 250     continue
 200  continue
      return
      end
c
c Lookup the PHREEQC internal ID
c
      integer function phreeq_id(STRING)
      USE WEBMOD_OBSCHEM, ONLY: phq_lut, nphq_lut
      IMPLICIT NONE
#include "fmodules.inc"
      CHARACTER(len=12) STRING 
!      character*12 sol_lut(MAXSOLUTE,2)
      integer i, len1, len2
      integer, external :: length
      !logical test1
      phreeq_id = 0
      i = 0
      do while (phreeq_id.eq.0)
         i = i+1
         if(i.gt.nphq_lut) then
            print*,"Species ", string(1:len1) ," cannot be "//
     $           "found in the species lookup table"
            print*," "
            phreeq_id = -99
            return
         end if
         len1 = length(string)
         len2 = length(phq_lut(i)%spec)
         if(STRING(1:len1).eq.(phq_lut(i)%spec(1:len2))) phreeq_id = i
      end do
      RETURN
      END
c
c Lookup the unit ID
c
      integer function unit_id(STRING, unit_lut)
      IMPLICIT NONE
      CHARACTER(len=12) STRING 
      character(len=12) unit_lut(:)
      integer i, len1, len2
      integer, external :: length
      unit_id = 0
      i = 0
      do while (unit_id.eq.0)
         i = i+1
         len1 = length(string)
         len2 = length(unit_lut(i))
         if(STRING(1:len1) .eq. unit_lut(i)(1:len2)) unit_id = i
         if(i.gt.8) then
            print*,"Concentration Unit ",string," cannot be "//
     $           "found in the unit lookup table"
            unit_id = -99
         end if
      end do
      RETURN
      END
c
c Function to return the length of the string using the standard
c intrinsic FORTRAN LEN function
c
      INTEGER FUNCTION LENGTH(STRING)
      IMPLICIT NONE
      CHARACTER(len=*) STRING 
      character(len=1) c
      INTEGER I
      DO 15, I = LEN(STRING), 1, -1
         c=string(i:i)
         IF(STRING(I:I) .NE. ' '.and.STRING(I:I)
     $        .ne.'\f'.and.string(i:i).ne.'\r') GO TO 16 
 15   CONTINUE 
 16   LENGTH = I
      RETURN
      END
c
c Function to parse the first word of a delimited string (before space or paren for elemental root)
c
      character(len=12) FUNCTION PARSE(STRING, DEL)
      IMPLICIT NONE
      CHARACTER(len=*) STRING
      CHARACTER(len=*) DEL
      INTEGER I
      DO 15, I = 2, LEN(STRING), 1
         IF(STRING(I:I) .EQ. DEL) GO TO 16 
 15   CONTINUE 
 16   PARSE = STRING(1:I-1)
      RETURN
      END
