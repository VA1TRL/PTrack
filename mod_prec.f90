!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
module mod_prec
  implicit none

  !--Double Precision Coding------------------------------------------------------!
  integer, parameter :: DP = selected_real_kind(15, 307)

end module mod_prec

