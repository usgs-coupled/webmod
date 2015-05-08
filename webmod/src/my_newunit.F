#include "defines.h"
      INTEGER FUNCTION my_newunit()
      IMPLICIT NONE
      INTEGER, SAVE :: next = 100
      LOGICAL :: isopen

      isopen = .TRUE.
      DO WHILE(isopen)
        my_newunit = next
        INQUIRE(UNIT=my_newunit, OPENED=isopen)
        next = next + 1
      END DO
      END FUNCTION my_newunit
