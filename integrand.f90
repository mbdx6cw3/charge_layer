     REAL FUNCTION integrand(P1)
     USE header_file, only: zchge, bulkc, nspec 
!    ***************************************************************************
!    Evaluates the integrand for the calculation of the concentration profile.
!    Called from qag.f90 routine.
!    ***************************************************************************
     implicit none
     real :: sum
     real, intent(in) :: P1
     integer :: i

     sum = 0.0d0
     do i = 1, nspec
       sum = sum + bulkc(i) * (P1**zchge(i) - 1.0)
     end do

     if(sum.le.0.0) then
       integrand = 1.0d12/P1
     else
       integrand = 1.0/P1/sqrt(sum)
     end if
!    ***************************************************************************
     END FUNCTION integrand
!    ---------------------------------------------------------------------------
