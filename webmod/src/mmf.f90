MODULE mmf
IMPLICIT NONE
SAVE
CONTAINS
    
    ! long decldim_ (char *dname, ftnint *dval, ftnint *dmax, char *ddescr, ftnlen namelen, ftnlen descrlen)    
    INTEGER(C_LONG) FUNCTION decldim(dname, dval, dmax, ddescr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER (C_LONG) FUNCTION decldimC                    &
                !(dname, dval, dmax, ddescr, namelen, descrlen)    &
                (dname, dval, dmax, ddescr)    &
                BIND(C, NAME='decldim')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dname(*)
                !INTEGER(KIND=C_INT),    INTENT(IN)        :: dval
                !INTEGER(KIND=C_INT),    INTENT(IN)        :: dmax
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: dval
                INTEGER(KIND=C_LONG),    INTENT(IN)        :: dmax
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ddescr(*)
                !INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: namelen
                !INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: descrlen
            END FUNCTION decldimC
        END INTERFACE
        CHARACTER(*),        INTENT(IN) :: dname        
        INTEGER, INTENT(IN) :: dval
        INTEGER, INTENT(IN) :: dmax
        CHARACTER(*),        INTENT(IN) :: ddescr        
        !decldim = decldim_(trim(dname)//C_NULL_CHAR,              &
        !                   dval, dmax, trim(ddescr)//C_NULL_CHAR, &
        !                   len(dname), len(ddescr))
        decldim = decldimC(trim(dname)//C_NULL_CHAR,              &
                           dval, dmax, trim(ddescr)//C_NULL_CHAR)
    END FUNCTION decldim

    ! long declfix_ (char *dname, ftnint *dval, ftnint *dmax, char *ddescr, ftnlen namelen, ftnlen descrlen)
    INTEGER(C_LONG) FUNCTION declfix(dname, dval, dmax, ddescr)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTERFACE
            INTEGER (C_LONG) FUNCTION declfix_                    &
                (dname, dval, dmax, ddescr, namelen, descrlen)    &
                BIND(C, NAME='declfix_')
                USE ISO_C_BINDING
                IMPLICIT NONE
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: dname(*)
                INTEGER(KIND=C_INT),    INTENT(IN)        :: dval
                INTEGER(KIND=C_INT),    INTENT(IN)        :: dmax
                CHARACTER(KIND=C_CHAR), INTENT(IN)        :: ddescr(*)
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: namelen
                INTEGER(KIND=C_INT),    INTENT(IN), VALUE :: descrlen
            END FUNCTION declfix_
        END INTERFACE
        CHARACTER(*),        INTENT(IN) :: dname        
        INTEGER(KIND=C_INT), INTENT(IN) :: dval
        INTEGER(KIND=C_INT), INTENT(IN) :: dmax
        CHARACTER(*),        INTENT(IN) :: ddescr        
        declfix = declfix_(trim(dname)//C_NULL_CHAR,              &
                           dval, dmax, trim(ddescr)//C_NULL_CHAR, &
                           len(dname), len(ddescr))
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
