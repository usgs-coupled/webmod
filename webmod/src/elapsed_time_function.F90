#include "defines.h"
function elapsed_time(t_opt) result(e_time)
      IMPLICIT NONE
  ! ... Purpose: calculate elapsed execution time, tenths of a second.
  ! ... Assumes that  execution time is no more than one month.
  ! ... First call sets reference time; subsequent calls return elapsed 
  ! ... time relative to reference time.
  ! ... Returns:
  ! ... t_opt=0: returns time in deci-seconds
  ! ... t_opt=1: returns time in centi-seconds
  ! ... t_opt=2: returns time in milli-seconds
  ! ... t_opt=other: returns time in seconds
  ! ... argument list
  ! ... 
  integer :: t_opt
  ! ... 
  ! ... result
  ! ... 
  integer :: e_time
  ! ... 
  ! ... local variables
  ! ... 
  integer, dimension(1:8) :: t_part
  integer, save, dimension(1:12) :: no_dapmo= &
       (/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer, save :: mo_s, da_s, h_s, mi_s, sec_s, msec_s
  logical, save :: initial=.true., flip
  integer :: mo_d, da_d, h_d, mi_d, sec_d, msec_d, t_sec, t_msec
  ! .....................................................................
  ! ... year=>t_part(1); month=>t_part(2); day=>t_part(3); hour=>t_part(5)
  ! ... minute=>t_part(6); second=>t_part(7); msecond=>t_part(8)
  call date_and_time(values=t_part)
  if (initial) then
     initial=.false.
     mo_s=t_part(2); da_s=t_part(3); h_s=t_part(5)
     mi_s=t_part(6); sec_s=t_part(7); msec_s=t_part(8)
     if (mod(t_part(1),4)==0) no_dapmo(2)=29
     e_time=0
     return
  else
     flip=.false.
     t_sec=0
     mo_d=t_part(2)-mo_s-1
     da_d=t_part(3)-da_s-1
     h_d=t_part(5)-h_s-1
     mi_d=t_part(6)-mi_s-1
     sec_d=t_part(7)-sec_s-1
     msec_d=t_part(8)-msec_s
     if (mo_d/=-1) then
        t_sec=t_sec+(no_dapmo(mo_s)+da_d)*86400
        flip=.true.
     elseif (da_d>-1) then
        t_sec=t_sec+da_d*86400
        flip=.true.
     endif
     if (flip) then
        t_sec=t_sec+(24+h_d)*3600
     elseif (h_d>-1) then
        t_sec=t_sec+h_d*3600
        flip=.true.
     endif
     if (flip) then
        t_sec=t_sec+(60+mi_d)*60
     elseif (mi_d>-1) then
        t_sec=t_sec+mi_d*60
        flip=.true.
     endif
     if (flip) then
        t_sec=t_sec+60+sec_d
     elseif (sec_d>-1) then
        t_sec=t_sec+sec_d
        flip=.true.
     endif
     ! ... milliseconds
     if (flip) then
        t_msec=1000+msec_d
     else
        t_msec=msec_d
     endif
     ! ... add second and milliseconds to form time measure
     if (t_opt==0) then
        e_time=t_sec*10+t_msec/100
        if (mod(t_msec,100)>50) e_time=e_time+1
     elseif (t_opt==1) then
        e_time=t_sec*100+t_msec/10
        if (mod(t_msec,10)>5) e_time=e_time+1
     elseif (t_opt==2) then
        e_time=t_sec*1000+t_msec
     else
        e_time=t_sec+t_msec/1000
        if (mod(t_msec,1000)>500) e_time=e_time+1
     endif
     ! ... 
  endif
  ! ... 
end function elapsed_time
