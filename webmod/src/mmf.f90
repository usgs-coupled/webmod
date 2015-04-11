MODULE mmf
IMPLICIT NONE

    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTERFACE declvar
        MODULE PROCEDURE declvar_real, declvar_realvec, declvar_realmat,   &
                         declvar_int,  declvar_intvec,                     &
                         declvar_dp,   declvar_dpvec,   declvar_dpmat,    declvar_dp3d
    END INTERFACE declvar

    ! long declpri (char *name, long size, char *type, char *value)
    INTERFACE declpri
        MODULE PROCEDURE declpri_real, declpri_realvec, declpri_realmat,   &
                         declpri_int,  declpri_intvec,  declpri_intmat,    &
                         declpri_dp,   declpri_dpmat
    END INTERFACE declpri

    SAVE
    CONTAINS
    
    ! long decldim (char *name, long value, long max, char *descr)
    INTEGER FUNCTION decldim(name, val, max, descr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION  decldimC  &
                (name, val, max, descr)         &
                BIND(C, NAME='decldim')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR),  INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: val
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: max
                CHARACTER(KIND=C_CHAR),  INTENT(IN)        :: descr(*)
            END FUNCTION decldimC
        END INTERFACE
        CHARACTER(*),        INTENT(IN) :: name
        INTEGER,             INTENT(IN) :: val
        INTEGER,             INTENT(IN) :: max
        CHARACTER(*),        INTENT(IN) :: descr
        INTEGER(KIND=C_LONG)            :: tval
        INTEGER(KIND=C_LONG)            :: tmax
        INTEGER(KIND=C_LONG)            :: creturn
        tval = val
        tmax = max
        creturn = decldimC(trim(name)//C_NULL_CHAR,                &
                           tval, tmax, trim(descr)//C_NULL_CHAR)
        decldim = int(creturn)
    END FUNCTION decldim

    ! long declfix (char *name, long value, long max, char *descr)
    INTEGER FUNCTION declfix(name, val, max, descr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declfixC   &
                (name, val, max, descr)         &
                BIND(C, NAME='declfix')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR),  INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: val
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: max
                CHARACTER(KIND=C_CHAR),  INTENT(IN)        :: descr(*)
            END FUNCTION declfixC
        END INTERFACE
        CHARACTER(*),        INTENT(IN) :: name
        INTEGER,             INTENT(IN) :: val
        INTEGER,             INTENT(IN) :: max
        CHARACTER(*),        INTENT(IN) :: descr
        INTEGER(KIND=C_LONG)            :: tval
        INTEGER(KIND=C_LONG)            :: tmax
        INTEGER(KIND=C_LONG)            :: creturn
        tval = val
        tmax = max
        creturn = declfixC(trim(name)//C_NULL_CHAR,                &
                           tval, tmax, trim(descr)//C_NULL_CHAR)
        declfix = int(creturn)
    END FUNCTION declfix
    
    ! long declparam (char *module, char *name, char *dimen, char *type, char *value,
    !     char *minimum, char *maximum, char *descr, char *help, char *units)     
    INTEGER FUNCTION declparam(mname, name, dimen, ptype, value, &
                               minimum, maximum, descr, help, units)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declparamC                   &
                (mname, name, dimen, ptype, value,                &
                 minimum, maximum, descr, help, units)            & 
                BIND(C, NAME='declparam')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: value(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: minimum(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: maximum(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: descr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
            END FUNCTION declparamC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: value
        CHARACTER(*), INTENT(IN)        :: minimum
        CHARACTER(*), INTENT(IN)        :: maximum
        CHARACTER(*), INTENT(IN)        :: descr
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declparamC(trim(mname)   // C_NULL_CHAR,  &
                             trim(name)    // C_NULL_CHAR,  &
                             trim(dimen)   // C_NULL_CHAR,  &
                             trim(ptype)   // C_NULL_CHAR,  &
                             trim(value)   // C_NULL_CHAR,  &
                             trim(minimum) // C_NULL_CHAR,  &
                             trim(maximum) // C_NULL_CHAR,  &
                             trim(descr)   // C_NULL_CHAR,  &
                             trim(help)    // C_NULL_CHAR,  &
                             trim(units)   // C_NULL_CHAR   &
                             )
        declparam = int(creturn)
    END FUNCTION declparam



    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_real(mname, name, dimen, maxsize, ptype, &
                                  help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER  (KIND=C_LONG), INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        REAL,         INTENT(IN), TARGET:: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_real = int(creturn)
    END FUNCTION declvar_real

    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_realvec(mname, name, dimen, maxsize, ptype, &
                                     help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        REAL,         INTENT(IN), TARGET:: value(:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_realvec = int(creturn)
    END FUNCTION declvar_realvec
                             

    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_realmat(mname, name, dimen, maxsize, ptype, &
                                        help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        REAL,         INTENT(IN), TARGET:: value(:,:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_realmat = int(creturn)
    END FUNCTION declvar_realmat
                                     
    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_int(mname, name, dimen, maxsize, ptype, &
                                 help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        INTEGER,      INTENT(IN), TARGET:: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_int = int(creturn)
    END FUNCTION declvar_int

    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_intvec(mname, name, dimen, maxsize, ptype, &
                                     help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        INTEGER,      INTENT(IN), TARGET:: value(:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_intvec = int(creturn)
    END FUNCTION declvar_intvec
                                     
    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_dp(mname, name, dimen, maxsize, ptype, &
                                help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        DOUBLE PRECISION, INTENT(IN), TARGET :: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_dp = int(creturn)
    END FUNCTION declvar_dp
                                   
                                   
    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_dpvec(mname, name, dimen, maxsize, ptype, &
                                   help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        DOUBLE PRECISION, INTENT(IN), TARGET:: value(:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_dpvec = int(creturn)
    END FUNCTION declvar_dpvec

    
    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_dpmat(mname, name, dimen, maxsize, ptype, &
                                   help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        DOUBLE PRECISION, INTENT(IN), TARGET:: value(:,:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_dpmat = int(creturn)
    END FUNCTION declvar_dpmat

    
    ! long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
    !               char *help, char *units, char *value)
    INTEGER FUNCTION declvar_dp3d(mname, name, dimen, maxsize, ptype, &
                                   help, units, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declvarC          &
                (mname, name, dimen, maxsize, ptype,   &
                 help, units, value)                   & 
                BIND(C, NAME='declvar')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dimen(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: maxsize
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: help(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: units(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declvarC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: mname
        CHARACTER(*), INTENT(IN)        :: name
        CHARACTER(*), INTENT(IN)        :: dimen
        INTEGER,      INTENT(IN)        :: maxsize
        CHARACTER(*), INTENT(IN)        :: ptype
        CHARACTER(*), INTENT(IN)        :: help
        CHARACTER(*), INTENT(IN)        :: units
        DOUBLE PRECISION, INTENT(IN), TARGET:: value(:,:,:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declvarC(trim(mname)  // C_NULL_CHAR,  &
                           trim(name)   // C_NULL_CHAR,  &
                           trim(dimen)  // C_NULL_CHAR,  &
                           maxsize,                      &
                           trim(ptype)  // C_NULL_CHAR,  &
                           trim(help)   // C_NULL_CHAR,  &
                           trim(units)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declvar_dp3d = int(creturn)
    END FUNCTION declvar_dp3d

    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_real(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        REAL,         INTENT(IN), TARGET:: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_real = int(creturn)
    END FUNCTION declpri_real
    
    

    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_realvec(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        REAL,         INTENT(IN), TARGET:: value(:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_realvec = int(creturn)
    END FUNCTION declpri_realvec
    
    

    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_realmat(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        REAL,         INTENT(IN), TARGET:: value(:,:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_realmat = int(creturn)
    END FUNCTION declpri_realmat
    
    
    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_int(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        INTEGER,      INTENT(IN), TARGET:: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_int = int(creturn)
    END FUNCTION declpri_int
    

    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_intvec(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        INTEGER,      INTENT(IN), TARGET:: value(:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_intvec = int(creturn)
    END FUNCTION declpri_intvec
    

    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_intmat(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        INTEGER,      INTENT(IN), TARGET:: value(:,:)
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_intmat = int(creturn)
    END FUNCTION declpri_intmat
    
    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_dp(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: name(*)
                INTEGER(KIND=C_LONG),   INTENT(IN)        :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                TYPE(C_PTR), VALUE                        :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*), INTENT(IN)        :: name
        INTEGER,      INTENT(IN)        :: size
        CHARACTER(*), INTENT(IN)        :: ptype
        DOUBLE PRECISION, INTENT(IN), TARGET:: value
        INTEGER(KIND=C_LONG)            :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_dp = int(creturn)
    END FUNCTION declpri_dp
    
    
    ! long declpri (char *name, long size, char *type, char *value)
    INTEGER FUNCTION declpri_dpmat(name, size, ptype, value)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER(C_LONG) FUNCTION declpriC          &
                (name, size, ptype, value)             & 
                BIND(C, NAME='declpri')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN) :: name(*)
                INTEGER  (KIND=C_LONG), INTENT(IN) :: size
                CHARACTER(KIND=C_CHAR), INTENT(IN) :: ptype(*)
                TYPE(C_PTR), VALUE                 :: value
            END FUNCTION declpriC
        END INTERFACE
        CHARACTER(*),     INTENT(IN)         :: name
        INTEGER,          INTENT(IN)         :: size
        CHARACTER(*),     INTENT(IN)         :: ptype
        DOUBLE PRECISION, INTENT(IN), TARGET :: value(:,:)
        INTEGER(KIND=C_LONG)                 :: creturn
        creturn = declpriC(trim(name)   // C_NULL_CHAR,  &
                           size,                         &
                           trim(ptype)  // C_NULL_CHAR,  &
                           C_LOC(value)                  &
                           )
        declpri_dpmat = int(creturn)
    END FUNCTION declpri_dpmat
    
END MODULE mmf
