MODULE mmf
IMPLICIT NONE
SAVE
CONTAINS
    
    ! long decldim (char *name, long value, long max, char *descr)
    INTEGER FUNCTION decldim(name, val, max, descr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER (C_LONG) FUNCTION decldimC  &
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
        decldim = creturn
    END FUNCTION decldim

    ! long declfix (char *name, long value, long max, char *descr)
    INTEGER FUNCTION declfix(name, val, max, descr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER (C_LONG) FUNCTION declfixC  &
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
        declfix = creturn
    END FUNCTION declfix
    
    ! long declparam_ (char *mname, char *pname, char *pdimen, char *ptype,
    !     char *pvalstr, char *minstr, char *maxstr, char *dstr, char *hstr,
    !     char *ustr, ftnlen mnamelen, ftnlen pnamelen, ftnlen pdimenlen, ftnlen ptypelen,
    !     ftnlen pvallen, ftnlen minlen, ftnlen maxlen, ftnlen dlen, ftnlen hlen, ftnlen ulen)
    INTEGER(C_LONG) FUNCTION declparam(mname, pname, pdimen, ptype,                   &
                                       pvalstr, minstr, maxstr, dstr, hstr,           &
                                       ustr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER (C_LONG) FUNCTION declparam_                  &
                (mname, pname, pdimen, ptype,                     &
                 pvalstr, minstr, maxstr, dstr, hstr,             &
                 ustr, mnamelen, pnamelen, pdimenlen, ptypelen,   &
                 pvallen, minlen, maxlen, dlen, hlen, ulen)       &
                BIND(C, NAME='declparam_')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: mname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: pname(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: pdimen(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ptype(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: pvalstr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: minstr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: maxstr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dstr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: hstr(*)
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ustr(*)               
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: mnamelen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: pnamelen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: pdimenlen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: ptypelen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: pvallen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: minlen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: maxlen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: dlen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: hlen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: ulen
            END FUNCTION declparam_
        END INTERFACE
        CHARACTER(*),        INTENT(IN) :: mname
        CHARACTER(*),        INTENT(IN) :: pname
        CHARACTER(*),        INTENT(IN) :: pdimen
        CHARACTER(*),        INTENT(IN) :: ptype
        CHARACTER(*),        INTENT(IN) :: pvalstr
        CHARACTER(*),        INTENT(IN) :: minstr
        CHARACTER(*),        INTENT(IN) :: maxstr
        CHARACTER(*),        INTENT(IN) :: dstr
        CHARACTER(*),        INTENT(IN) :: hstr
        CHARACTER(*),        INTENT(IN) :: ustr
        declparam = declparam_(trim(mname)  // C_NULL_CHAR,  &
                               trim(pname)  // C_NULL_CHAR,  &
                               trim(pdimen) // C_NULL_CHAR,  &
                               trim(ptype)  // C_NULL_CHAR,  &
                               trim(pvalstr)// C_NULL_CHAR,  &
                               trim(minstr) // C_NULL_CHAR,  &
                               trim(maxstr) // C_NULL_CHAR,  &
                               trim(dstr)   // C_NULL_CHAR,  &
                               trim(hstr)   // C_NULL_CHAR,  &
                               trim(ustr)   // C_NULL_CHAR,  &
                               len(mname),                   &
                               len(pname),                   &
                               len(pdimen),                  &
                               len(ptype),                   &
                               len(pvalstr),                 &
                               len(minstr),                  &
                               len(maxstr),                  &
                               len(dstr),                    &
                               len(hstr),                    &
                               len(ustr))
    END FUNCTION declparam
END MODULE mmf
