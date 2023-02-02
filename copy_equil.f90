     MODULE equil_module
!    ---------------------------------------------------------------------------

     CONTAINS 
     SUBROUTINE equilibria(nspec1,trialc,func,derivc,nspec2,iflag)
     USE header_file, only: zchge,bulkcO,double,abstrial,coeff,tc_coeff,s_coeff,&
   & fe_coeff,switch
!    ***************************************************************************
     implicit none
     integer(kind=4), intent(in) :: nspec1, nspec2, iflag
     real(kind=double), intent(inout) :: func(nspec1)
     real(kind=double), intent(in) :: trialc(nspec1)
     real(kind=double), intent(inout) :: derivc(nspec2,nspec1)
     integer :: i

     ! if iflag = 2 then calculate the Jacobian, otherwise calculate the function
     if (iflag.eq.2) then

       call jacobian(nspec2,trialc,derivc)

     end if
     
     ! loop over all ionic species
     func(1) = 0.0d0
     do i = 1, nspec1

       func(1) = func(1) + zchge(i) * trialc(i)

     end do

     ! charge neutrality
     func(1) = func(1) / abstrial                                                       !F1  = ensures charge neutrality

     if (switch.eq.1) then

       ! fixed concentrations of FeCl3, KTcO4, K2SO4 and pH
       func(2)  = trialc(1) * trialc(2) / coeff - 1.0                                   !F2  = [H+][OH-]/Kw - 1
       func(3)  = trialc(1) * trialc(4) / (trialc(5) * tc_coeff) - 1.0                  !F3  = [H+][TcO4-]/[HTcO4]KTc - 1
       func(4)  = trialc(1) * trialc(6) / (trialc(7) * s_coeff(1)) - 1.0                !F4  = [H+][SO4--]/[HSO4-]KS1 - 1
       func(5)  = trialc(1) * trialc(7) / (trialc(8) * s_coeff(2)) - 1.0                !F5  = [H+][HSO4-]/[H2SO4]KS2 - 1
       func(6)  = trialc(1) * trialc(11) / (trialc(9) * fe_coeff(1)) - 1.0              !F6  = [H+][FeOH++]/[Fe3+]KFe1 - 1
       func(7)  = trialc(1)**2 * trialc(12) / (trialc(9) * fe_coeff(2)) - 1.0           !F7  = [H+]**2[Fe(OH)2+]/[Fe3+]KFe2 - 1
       func(8)  = trialc(1)**3 * trialc(13) / (trialc(9) * fe_coeff(3)) - 1.0           !F8  = [H+]**3[Fe(OH)3]/[Fe3+]KFe3 - 1
       func(9)  = trialc(1)**4 * trialc(14) / (trialc(9) * fe_coeff(4)) - 1.0           !F9  = [H+]**4[Fe(OH)4-]/[Fe3+]KFe4 - 1
       func(10) = trialc(1) / bulkcO(1) - 1.0    ! fixed pH                             !F10 = [H+]/[H+]o - 1
       func(11) = trialc(3) / bulkcO(3) - 1.0    ! fixed [K+]                           !F11 = [K+]/[K+]o - 1
       func(12) = (trialc(4) + trialc(5))/ bulkcO(4) - 1.0    ! fixed [Tc]              !F12  = [TcO4-]+[HTcO4]/[TcO4-]o - 1
       func(13) = (trialc(6)+trialc(7)+trialc(8))/ bulkcO(6) - 1.0 ! fixed [SO42-]      !F13  = [SO4--]+[HSO4-]+[H2SO4]/[SO4--]o - 1
       func(14) = trialc(10) / bulkcO(10) - 1.0  ! fixed [Cl-]                          !F14 = [Cl-]/[Cl-]o - 1
  
     else if (switch.eq.2) then

       func(2)  = trialc(1) * trialc(2) / coeff - 1.0                                   !F2  = [H+][OH-]/Kw - 1
       func(3)  = trialc(1) * trialc(4) / (trialc(5) * tc_coeff) - 1.0                  !F3  = [H+][TcO4-]/[HTcO4]KTc - 1
       func(4)  = trialc(1) * trialc(6) / (trialc(7) * s_coeff(1)) - 1.0                !F4  = [H+][SO4--]/[HSO4-]KS1 - 1
       func(5)  = trialc(1) * trialc(7) / (trialc(8) * s_coeff(2)) - 1.0                !F5  = [H+][HSO4-]/[H2SO4]KS2 - 1
       func(6)  = trialc(1) / bulkcO(1) - 1.0    ! fixed pH                             !F6  = [H+]/[H+]o - 1
       func(7)  = trialc(3) / bulkcO(3) - 1.0    ! fixed [K+]                           !F7  = [K+]/[K+]o - 1
       func(8)  = (trialc(4) + trialc(5))/ bulkcO(4) - 1.0    ! fixed [Tc]              !F8  = [TcO4-]+[HTcO4]/[TcO4-]o - 1
       func(9)  = (trialc(6)+trialc(7)+trialc(8))/ bulkcO(6) - 1.0 ! fixed [SO42-]      !F9  = [SO4--]+[HSO4-]+[H2SO4]/[SO4--]o - 1
       func(10) = trialc(10) / bulkcO(10) - 1.0  ! fixed [Cl-]                          !F10 = [Cl-]/[Cl-]o - 1

     else if (switch.eq.3) then

       ! fixed concentrations of FeCl3, KTcO4, K2SO4 and pH
       func(2)  = trialc(1) * trialc(2) / coeff - 1.0                                   !F2  = [H+][OH-]/Kw - 1
       func(3)  = trialc(1) * trialc(4) / (trialc(5) * tc_coeff) - 1.0                  !F3  = [H+][TcO4-]/[HTcO4]KTc - 1
       func(4)  = trialc(1) * trialc(6) / (trialc(7) * s_coeff(1)) - 1.0                !F4  = [H+][SO4--]/[HSO4-]KS1 - 1
       func(5)  = trialc(1) * trialc(7) / (trialc(8) * s_coeff(2)) - 1.0                !F5  = [H+][HSO4-]/[H2SO4]KS2 - 1
       func(6)  = trialc(1) / bulkcO(1) - 1.0                                           !F6  = [H+]/[H+]o - 1
       func(7)  = (trialc(4) + trialc(5))/ bulkcO(4) - 1.0    ! fixed [Tc]              !F7  = [TcO4-]+[HTcO4]/[TcO4-]o - 1
       func(8)  = (trialc(6)+trialc(7)+trialc(8))/ bulkcO(6) - 1.0 ! fixed [SO42-]      !F8  = [SO4--]+[HSO4-]+[H2SO4]/[SO4--]o - 1
 
     end if

!    ***************************************************************************
     END SUBROUTINE equilibria
!    ---------------------------------------------------------------------------
     
     SUBROUTINE jacobian(nspec,trialc,derivc)
     USE header_file, only: zchge,bulkcO,double,abstrial,coeff,tc_coeff,s_coeff,&
   & fe_coeff,coeff,switch
!    ***************************************************************************
     implicit none
     real(kind=double), intent(inout) :: derivc(nspec,nspec)
     real(kind=double), intent(in) :: trialc(nspec)
     integer(kind=4), intent(in) :: nspec
     integer(kind=4) :: i,j
   
     ! derivc is the jacobian partial derivative matrix containing nspec
     ! functions and nspec variables.
     derivc(:,:) = 0.0d0

     ! charge neutrality
     do i = 1, nspec
     
       derivc(1,i) = zchge(i) / abstrial

     end do

     if (switch.eq.1) then

       ! pH is controlled directly and fixed concentrations of FeCl3, KTcO4, K2SO4
       derivc(2,1)  =   trialc(2) / coeff
       derivc(3,1)  =   trialc(4) / trialc(5) / tc_coeff
       derivc(4,1)  =   trialc(6) / trialc(7) / s_coeff(1)
       derivc(5,1)  =   trialc(7) / trialc(8) / s_coeff(2)
       derivc(6,1)  =   trialc(11) / trialc(9) / fe_coeff(1)
       derivc(7,1)  =   2.0d0 * trialc(1) * trialc(12) / trialc(9) / fe_coeff(2)
       derivc(8,1)  =   3.0d0 * trialc(1)**2 * trialc(13) / trialc(9) / fe_coeff(3)
       derivc(9,1)  =   4.0d0 * trialc(1)**3 * trialc(14) / trialc(9) / fe_coeff(4)
       derivc(10,1) =   1 / bulkcO(1)
       derivc(2,2)  =   trialc(1) / coeff
       derivc(11,3)  =   1 / bulkcO(3)
       derivc(3,4)   =   trialc(1) / trialc(5) / tc_coeff
       derivc(12,4)  =   trialc(5) / bulkcO(4)
       derivc(12,5)  =   trialc(4) / bulkcO(4)
       derivc(3,5)  =  -trialc(1) * trialc(4) / trialc(5)**2 / tc_coeff
       derivc(4,6)  =   trialc(1) / trialc(7) / s_coeff(1)
       derivc(13,6)  =   (trialc(7) + trialc(8)) / bulkcO(6)
       derivc(13,7)  =   (trialc(6) + trialc(8)) / bulkcO(6)
       derivc(13,8)  =   (trialc(6) + trialc(7)) / bulkcO(6)
       derivc(4,7)  =  -trialc(1) * trialc(6) / trialc(7)**2 / s_coeff(1)
       derivc(5,7)  =   trialc(1) / trialc(8) / s_coeff(2)
       derivc(5,8)  =  -trialc(1) * trialc(7) / trialc(8)**2 / s_coeff(2)
       derivc(6,9)  =  -trialc(1) * trialc(11) / trialc(9)**2 / fe_coeff(1)
       derivc(7,9)  =  -trialc(1)**2 * trialc(12) / trialc(9)**2 / fe_coeff(2)
       derivc(8,9)  =  -trialc(1)**3 * trialc(13) / trialc(9)**2 / fe_coeff(3)
       derivc(9,9)  =  -trialc(1)**4 * trialc(14) / trialc(9)**2 / fe_coeff(4)
       derivc(14,10)=   1 / bulkcO(10)
       derivc(6,11) =   trialc(1) / trialc(9) / fe_coeff(1)
       derivc(7,12) =   trialc(1)**2 / trialc(9) / fe_coeff(2)
       derivc(8,13) =   trialc(1)**3 / trialc(9) / fe_coeff(3)
       derivc(9,14) =   trialc(1)**4 / trialc(9) / fe_coeff(4)

     else if (switch.eq.2) then

       derivc(2,1)  =   trialc(2) / coeff
       derivc(3,1)  =   trialc(4) / trialc(5) / tc_coeff
       derivc(4,1)  =   trialc(6) / trialc(7) / s_coeff(1)
       derivc(5,1)  =   trialc(7) / trialc(8) / s_coeff(2)
       derivc(6,1)  =   1 / bulkcO(1)
       derivc(2,2)  =   trialc(1) / coeff
       derivc(7,3)  =   1 / bulkcO(3)
       derivc(3,4)  =   trialc(1) / trialc(5) / tc_coeff
       derivc(8,4)  =   trialc(5) / bulkcO(4)
       derivc(8,5)  =   trialc(4) / bulkcO(4)
       derivc(3,5)  =  -trialc(1) * trialc(4) / trialc(5)**2 / tc_coeff
       derivc(4,6)  =   trialc(1) / trialc(7) / s_coeff(1)
       derivc(9,6)  =   (trialc(7) + trialc(8)) / bulkcO(6)
       derivc(9,7)  =   (trialc(6) + trialc(8)) / bulkcO(6)
       derivc(9,8)  =   (trialc(6) + trialc(7)) / bulkcO(6)
       derivc(9,7)  =  -trialc(1) * trialc(6) / trialc(7)**2 / s_coeff(1)
       derivc(5,7)  =   trialc(1) / trialc(8) / s_coeff(2)
       derivc(5,8)  =  -trialc(1) * trialc(7) / trialc(8)**2 / s_coeff(2)
       derivc(10,10)=   1 / bulkcO(10)

     else if (switch.eq.3) then

       derivc(2,1)  =   trialc(2) / coeff
       derivc(3,1)  =   trialc(4) / trialc(5) / tc_coeff
       derivc(4,1)  =   trialc(6) / trialc(7) / s_coeff(1)
       derivc(5,1)  =   trialc(7) / trialc(8) / s_coeff(2)
       derivc(6,1)  =   1 / bulkcO(1)
       derivc(2,2)  =   trialc(1) / coeff
       derivc(7,4)  =   trialc(5) / bulkcO(4)
       derivc(3,4)  =   trialc(1) / trialc(5) / tc_coeff
       derivc(3,5)  =  -trialc(1) * trialc(4) / trialc(5)**2 / tc_coeff
       derivc(7,5)  =   trialc(4) / bulkcO(4)
       derivc(4,6)  =   trialc(1) / trialc(7) / s_coeff(1)
       derivc(4,7)  =  -trialc(1) * trialc(6) / trialc(7)**2 / s_coeff(1)
       derivc(5,7)  =   trialc(1) / trialc(8) / s_coeff(2)
       derivc(5,8)  =  -trialc(1) * trialc(7) / trialc(8)**2 / s_coeff(2)
       derivc(8,6)  =   (trialc(7) + trialc(8)) / bulkcO(6)
       derivc(8,7)  =   (trialc(6) + trialc(8)) / bulkcO(6)
       derivc(8,8)  =   (trialc(6) + trialc(7)) / bulkcO(6)
     end if

!    *****************************************************************************
     END SUBROUTINE jacobian
     END MODULE equil_module
