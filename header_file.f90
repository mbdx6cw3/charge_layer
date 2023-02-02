     MODULE header_file
!    **************************************************************************
!    This module provides the global variables for the charge_layer.f90 
!    program. Includes a namelist for the input parameters. 
!    **************************************************************************

     implicit none
     integer, parameter :: double = SELECTED_REAL_KIND(15,99)

!    global variables
!    --------------------------------------------------------------------------
     integer(kind=4) :: i, j, k, switch, nspec,nspec1,nspec2,npoint,info
     real(kind=double), allocatable, dimension (:) :: conc,bulkcO,bulkc,func,trialc
     real(kind=double), allocatable, dimension (:) :: zchge,atrues,btrues
     character(len=10), allocatable, dimension (:) :: ion_type
     character(len=20) :: keyword1, keyword2
     real(kind=double), allocatable, dimension (:,:)::  derivc
     real(kind=double) :: errrel, abstrial, charge_dens, ion_strength, b0, area
     real(kind=double) :: Ka, KFe, KK, pH
     real :: P0, P
     real(kind=double) :: coeff, fe_coeff(4), tc_coeff, s_coeff(2)
     real(kind=double) :: bj_len, sum_debye, dh_len, dist0

!    Grahame equation
!    --------------------------------------------------------------------------
     real(kind=double), allocatable, dimension (:,:) :: zconc, cscoeff
     real, allocatable, dimension (:) :: charge, distance
     real, allocatable, dimension (:) :: err,  splwork 
     real(kind=double), allocatable, dimension (:) :: splderiv, potential
     complex(8), allocatable, dimension (:) :: polycoeff, roots
     real(kind=double), allocatable, dimension (:) :: realroot
     real(kind=double) :: fract, fractH, fractK, fractFe
     real(kind=double) :: realpart
     real(kind=double) :: sigma_surf, sigma_sol
     logical :: lpolish = .true., lspline = .true.
     logical :: lactive = .false., lallocate = .true., wrong_root = .false.
     integer :: froot, nroot, neval, ier
     real, external :: pchqa, integrand
     
!    activity
!    --------------------------------------------------------------------------
     integer(kind=4),parameter :: npitz = 10
     logical,allocatable,dimension(:) :: convergence
     real(kind=double),allocatable,dimension(:) :: loggam,gamma,gamma_old
     real(kind=double),allocatable,dimension (:,:) :: pitzer_b0, pitzer_b1
     real(kind=double),allocatable,dimension (:,:) :: pitzer_lamda, pitzer_phi
     real(kind=double),allocatable,dimension (:,:,:) :: pitzer_psi
     real(kind=double),allocatable,dimension (:,:) :: pitzer_c

!    parameters
!    --------------------------------------------------------------------------
     integer(kind=4), parameter   :: irule=2 
     real, parameter :: relerr = 1.0d-3, errabs = 0.0d0
     real(kind=double), parameter :: e      = 1.602d-19    ! C
     real(kind=double), parameter :: e2     = e**2         ! C**2
     real(kind=double), parameter :: eps    = 80.0
     real(kind=double), parameter :: eps0   = 8.854d-12    ! C V**-1 m**-1
     real(kind=double), parameter :: kb     = 1.3806d-23   ! J K**-1
     real(kind=double), parameter :: T      = 298          ! K
     real(kind=double), parameter :: pi     = 3.1415926d0 
     real(kind=double), parameter :: fourpi = 4*pi
     real(kind=double), parameter :: avogad = 6.022d23     ! mol**-1

!    input deck
!    --------------------------------------------------------------------------
     namelist /input_deck/ nspec,npoint,charge_dens,pH,Ka,KFe,KK,bulkcO,coeff,&
                        &  fe_coeff,tc_coeff,s_coeff,atrues,btrues,zchge
!    **************************************************************************
     END MODULE header_file
