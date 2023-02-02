     PROGRAM charge_layer
     USE header_file
     USE equil_module
     USE grahame_module
!    *************************************************************************
!    *************************************************************************
!    Description:
!
!      Gouy-Chapman calculation to determine ionic concentration profiles
!      near a charged surface of known charge density. The Grahame equation
!      is used to find the surface potential.
!
!      Three cases are considered:
!
!        1. Electrolyte solution of FeCl3, KTcO4 and K2SO4. Fe+++ is allowed
!           to precipitate out, forming FeOH++, Fe(OH)2+. Fe(OH)3 and
!           Fe(OH)4- species. This is analagous to a bare surface with no
!           tethers to complex Fe+++ (keyword1 = 'iron_ppt')
!   
!        2. Solution as above but Fe+++ is prevented from forming any
!           hydroxides. This is analagous to a surface which has surface
!           tethers to complex Fe+++. This assumes the tethers complex all
!           of the Fe+++ but does not have any contribution to the surface
!           potential. (keyword = 'iron_tether')
!
!        3. Electrolyte solution containing the component ions of KTcO4,
!           and K2SO4. Fe+++ is not involved. (keyword = 'zero_iron')
!
!    Notes:
!
!      a) the surface charge must balance the total charge of all ionic
!         species in solution    
!
!      b) The number of ionic species, nspec, is dependent on the type of
!         calculation. For arrays with nspec elements, array elements 
!         correspond to: 
!
!       case =           1.          2.          3. 
! 
!            1.          H+          H+          H+
!            2.          OH-         OH-         OH-
!            3.          K+          K+          K+
!            4.          TcO4-       TcO4-       TcO4-
!            5.          HTcO4       HTcO4       HTcO4
!            6.          SO4--       SO4--       SO4--
!            7.          HSO4-       HSO4-       HSO4-
!            8.          H2SO4       H2SO4       H2SO4   
!            9.          Fe+++       Fe+++
!           10.          Cl-         Cl-
!           11.          FeOH++
!           12.          Fe(OH)2+
!           13.          Fe(OH)3
!           14.          Fe(OH)4-
!
!      c) equilibrium constants:
!           coeff          K (H+    + OH-   ->  H2O            )
!           tc_coeff       K (TcO4- + H+    ->  HTcO4          )
!           s_coeff(1)     K (SO4-- + H+    ->  HSO4-          )
!           s_coeff(2)     K (HSO4- + H+    ->  H2SO4          )
!           fe_coeff(1)    K (Fe+++ + H2O   ->  Fe(OH)++ + H+  )
!           fe_coeff(2)    K (Fe+++ + 2H2O  ->  Fe(OH)2+ + 2H+ )
!           fe_coeff(3)    K (Fe+++ + 3H2O  ->  Fe(OH)3  + 3H+ )
!           fe_coeff(4)    K (Fe+++ + 4H2O  ->  Fe(OH)4- + 4H+ )
!
!      d) equilibrium constants for cations and the surface are defined 
!         as the formation constants instead of the normal pKa definition
!         (dissociation constants), ie. Kf = 1/Ka = [ML]/[M][L]
!
!      e) pH is fixed in each case
!
!      f) activities can be recomputed until they converge, using the
!         Truesdell-Jones equation (keyword2 = 'truesdell-jones') or the
!         Pitzer model (keyword2 = 'Pitzer'). This loop can be missed 
!         out entirely (keyword2 = 'debye-huckel')
!
!      g) coeffients for the Truesdell-Jones equation, atrues and btrues,
!         for calculating activity coefficients are from 'Aqueous
!         Environmental Chemistry', D Langmuir, Prentice Hall.
!
!      h) global variables, parameters and input deck defined in header_file
!
!      i) input data read in from input.dat file
!
!    Requires modules:
!
!      a) equil.mod   - contains subroutines required for the evaluation
!                       of the initial bulk concentration of ions:
!
!        	       i) equilibria
!                     ii) jacobian
!
!      b) grahame.mod - contains subroutines required for solving the
!                       the Grahame equation
!
!                      i) polynom
!                     ii) zroots - calls laguer
!                    iii) sol_charge
!                     iv) surf_charge
!                     v)  activity
!
!    External functions:
!
!      a) integrand.f90
!
!    Subroutines from common libraries required:
!
!      a) hybrj1 - solves a set of nonlinear equations 
!      b) qag    - approximates an integral over a finite interval
!      c) pchez  - carries out easy to use spline or Hermite interpolation
!
!    Christopher D. Williams (Sept. 2013)
!
!    Adapted from Space Charge program (JH Harding)
!    ***************************************************************************
!    ***************************************************************************  
     implicit none

!    --------------------- OPEN AND READ INPUT FILES ----------------------------
 
     open(unit = 1, file = 'input.dat', status = 'old')
     read(1,*) keyword1
     
     select case(keyword1)
     case ('iron_ppt')
     switch = 1
     nspec = 14

     case ('iron_tether')
     switch = 2
     nspec = 10

     case('zero_iron')
     switch = 3
     nspec = 8
     end select     

     ! allocate arrays
     allocate(bulkcO(nspec),zchge(nspec),atrues(nspec),btrues(nspec), &
    & ion_type(nspec),bulkc(nspec),gamma(nspec),loggam(nspec),gamma_old(nspec) &
    & ,convergence(nspec))
     atrues(:)    = 0.0d0
     btrues(:)    = 0.0d0
     bulkcO(:)    = 0.0d0
     zchge(:)     = 0.0d0
     bulkc(:)     = 0.0d0
     gamma(:)     = 0.0d0
     loggam(:)    = 0.0d0
     gamma_old(:) = 0.0d0
     convergence(:) = .false.

     read(1,*) keyword2
     select case (keyword2)

     ! switch Truesdell-Jones activity loop on
     case('truesdell-jones')
     lactive = .true.

     ! switch Pitzer model activity loop on
     case('pitzer')
     lactive = .true.

     ! switch activity loop off
     case('debye-huckel')
     lactive = .false.

     end select

     read(unit = 1, nml = input_deck)
     close(unit = 1)

     ! allocate pitzer parameter arrays
     allocate(pitzer_b0(npitz,npitz),pitzer_b1(npitz,npitz), &
    & pitzer_c(npitz,npitz),pitzer_phi(npitz,npitz),pitzer_lamda(npitz,npitz), &
    & pitzer_psi(npitz,npitz,npitz))
    
     ! read in virial coefficients for pitzer model
     open(unit=2,file='pitzer.dat',status='old')

     do i = 1, npitz
       read(2,*) (pitzer_b0(i,j),j=1,npitz)
     end do
  
     do i = 1, npitz
       read(2,*) (pitzer_b1(i,j),j=1,npitz)
     end do

     do i = 1, npitz
       read(2,*) (pitzer_c(i,j),j=1,npitz)
     end do
  
     do i = 1, npitz
       read(2,*) (pitzer_phi(i,j),j=1,npitz)
     end do

     do i = 1, npitz
       read(2,*) (pitzer_lamda(i,j),j=1,npitz)
     end do

     pitzer_psi(:,:,:) = 0.0d0

     read(2,*) pitzer_psi(1,3,6)
     read(2,*) pitzer_psi(1,3,7)
     read(2,*) pitzer_psi(3,2,6)
     read(2,*) pitzer_psi(3,6,7)
     read(2,*) pitzer_psi(3,4,6)

     ! set up pitzer psi arrays - copy array elements
     do i = 1, npitz
       do j = 1, npitz
         do k = 1, npitz
           pitzer_psi(i,k,j) = pitzer_psi(i,j,k)
           pitzer_psi(j,i,k) = pitzer_psi(i,j,k)
           pitzer_psi(j,k,i) = pitzer_psi(i,j,k)
           pitzer_psi(k,j,i) = pitzer_psi(i,j,k)
           pitzer_psi(k,i,j) = pitzer_psi(i,j,k)
         end do
       end do
     end do
     close(unit=2)
!    ---------------------------------------------------------------------------


!    ------------------------------ SETUP --------------------------------------

     open(unit = 5, file = 'output.dat', status = 'unknown')
     write(5,*) '-------------- OUTPUT DATA --------------'
     
     write(5,*)
     write(5,*) '--------- SIMULATION PARAMETERS ---------'
     write(5,*) 'Calculation type:', keyword1
     write(5,*) 'Number of species:', nspec
     write(5,*) 'Number of points:', npoint
     write(5,*) 'Charge density of surface (OH/nm**2):', charge_dens 
     write(5,*) 

     if (switch.eq.1) then 
       write(5,*)
       write(5,*) '--------- EQUILIBRIUM CONSTANTS ---------'
       write(5,*) 'H+', Ka
       write(5,*) 'K+', KK
       write(5,*) 'Fe+++', KFe
       write(5,*) 

       ! write initial bulk concentrations of FeCl3, KTcO4, K2SO42 and pH
       write(5,*)
       write(5,*) '--- INITIAL CONCENTRATIONS (mol/l) ---'
       write(5,*) 'FeCl3:', bulkcO(9)
       write(5,*) 'KTcO4:', bulkcO(4)
       write(5,*) 'K2SO4:', bulkcO(6)
       write(5,*) 'Approximate pH:', pH 
       write(5,*) 

       ! calculate rough estimate of initial bulk concentrations of all species   
       bulkcO(1)  = 10.0d0**(-pH)
       bulkcO(2)  = coeff/bulkcO(1)                                                 ! [OH-] = K/[H+]
       bulkcO(3)  = bulkcO(4)+2*bulkcO(6)                                           ! [K+] = [KTcO4] + [K2SO4]
       bulkcO(4)  = bulkcO(4)                                                       ! [TcO4-] = [KTcO4]
       bulkcO(5)  = bulkcO(1)*bulkcO(4)/tc_coeff                                    ! [HTcO4] = [H+][TcO4-]/K
       bulkcO(6)  = bulkcO(6)                                                       ! [SO4--] = [K2SO4]
       bulkcO(7)  = bulkcO(1)*bulkcO(6)/s_coeff(1)                                  ! [HSO4-] = [H+][SO4--]/K
       bulkcO(8)  = bulkcO(1)*bulkcO(7)/s_coeff(2)                                  ! [H2SO4] = [H+][HSO4-]/K
       bulkcO(9)  = bulkcO(9)                                                       ! [Fe+++] = [FeCl3]
       bulkcO(10) = 3.0d0*bulkcO(9)                                                 ! [Cl-] = 3*[FeCl3]
       bulkcO(11) = fe_coeff(1)*bulkcO(9)/bulkcO(1)                                 ! [FeOH++] = K[Fe3+]/[H+]
       bulkcO(12) = fe_coeff(2)*bulkcO(9)/(bulkcO(1)**2)                            ! [Fe(OH)2+] = K[Fe3+]/[H+]**2
       bulkcO(13) = fe_coeff(3)*bulkcO(9)/(bulkcO(1)**3)                            ! [Fe(OH)3] = K[Fe3+]/[H+]**3
       bulkcO(14) = fe_coeff(4)*bulkcO(9)/(bulkcO(1)**4)                            ! [Fe(OH)4-] = K[Fe3+]/[H+]**4

       ! assign each ion type names
       ion_type(1)  = 'H+' 
       ion_type(2)  = 'OH-'
       ion_type(3)  = 'K+'
       ion_type(4)  = 'TcO4-'
       ion_type(5)  = 'HTcO4'
       ion_type(6)  = 'SO4--'
       ion_type(7)  = 'HSO4-'
       ion_type(8)  = 'H2SO4'
       ion_type(9)  = 'Fe+++'
       ion_type(10) = 'Cl-'
       ion_type(11) = 'FeOH++'
       ion_type(12) = 'Fe(OH)2+'
       ion_type(13) = 'Fe(OH)3'
       ion_type(14) = 'Fe(OH)4-'

     else if (switch.eq.2) then
    
       write(5,*)
       write(5,*) '--------- EQUILIBRIUM CONSTANTS ---------'
       write(5,*) 'H+', Ka
       write(5,*) 'K+', KK
       write(5,*) 'Fe+++', KFe
       write(5,*)

       ! write initial bulk concentrations of FeCl3, KTcO4, K2SO42 and pH
       write(5,*)
       write(5,*) '--- INITIAL CONCENTRATIONS (mol/l) ---'
       write(5,*) 'FeCl3:', bulkcO(9)
       write(5,*) 'KTcO4:', bulkcO(4)
       write(5,*) 'K2SO4:', bulkcO(6)
       write(5,*) 'Approximate pH:', pH
       write(5,*) 
       
       ! calculate rough estimate of initial bulk concentrations of all species
       bulkcO(1)  = 10.0d0**(-pH)
       bulkcO(2)  = coeff/bulkcO(1)                                                 ! [OH-] = K/[H+]
       bulkcO(3)  = bulkcO(4)+2*bulkcO(6)                                           ! [K+] = [KTcO4] + [K2SO4]
       bulkcO(4)  = bulkcO(4)                                                       ! [TcO4-] = [KTcO4]
       bulkcO(5)  = bulkcO(1)*bulkcO(4)/tc_coeff                                    ! [HTcO4] = [H+][TcO4-]/K
       bulkcO(6)  = bulkcO(6)                                                       ! [SO4--] = [K2SO4]
       bulkcO(7)  = bulkcO(1)*bulkcO(6)/s_coeff(1)                                  ! [HSO4-] = [H+][SO4--]/K
       bulkcO(8)  = bulkcO(1)*bulkcO(7)/s_coeff(2)                                  ! [H2SO4] = [H+][HSO4-]/K
       bulkcO(9)  = bulkcO(9)                                                       ! [Fe+++] = [FeCl3]
       bulkcO(10) = 3.0d0*bulkcO(9)                                                 ! [Cl-] = 3*[FeCl3]

       ! assign each ion type names
       ion_type(1)  = 'H+' 
       ion_type(2)  = 'OH-'
       ion_type(3)  = 'K+'
       ion_type(4)  = 'TcO4-'
       ion_type(5)  = 'HTcO4'
       ion_type(6)  = 'SO4--'
       ion_type(7)  = 'HSO4-'
       ion_type(8)  = 'H2SO4'
       ion_type(9)  = 'Fe+++'
       ion_type(10) = 'Cl-'

     ! control pH, ionic strength and [K+] 
     else if (switch.eq.3) then

       write(5,*)
       write(5,*) '--------- SURFACE COMPLEXATION CONSTANTS ---------'
       write(5,*) 'H+:', Ka
       write(5,*) 'K+:', KK
       write(5,*) 'Fe+++:', KFe
       write(5,*)
       
       ! write initial bulk concentrations of KTcO4, K2SO42 and pH
       write(5,*)
       write(5,*) '--------- INITIAL CONCENTRATIONS (mol/l) ----------'
       write(5,*) 'FeCl3:', '       0.000'
       write(5,*) 'KTcO4:', bulkcO(4)
       write(5,*) 'K2SO4:', bulkcO(6)
       write(5,*) 'Approximate pH:', pH
       write(5,*)

       ! calculate rough estimate of initial bulk concentrations of all species
       bulkcO(1)  = 10.0**(-pH)
       bulkcO(2)  = coeff/bulkcO(1)                                                 ! [OH-] = K/[H+]
       bulkcO(3)  = bulkcO(4) + 2 * bulkcO(6)                                       ! [K+] = [KTcO4] + [K2SO4]
       bulkcO(4)  = bulkcO(4)                                                       ! [TcO4-] = [KTcO4]
       bulkcO(5)  = bulkcO(1)*bulkcO(4)/tc_coeff                                    ! [HTcO4] = [H+][TcO4-]/K
       bulkcO(6)  = bulkcO(6)                                                       ! [SO4--] = [K2SO4]
       bulkcO(7)  = bulkcO(1)*bulkcO(6)/s_coeff(1)                                  ! [HSO4-] = [H+][SO4--]/K
       bulkcO(8)  = bulkcO(1)*bulkcO(7)/s_coeff(2)                                  ! [H2SO4] = [H+][HSO4-]/K

       ! assign each ion type names
       ion_type(1)  = 'H+'
       ion_type(2)  = 'OH-'
       ion_type(3)  = 'K+'
       ion_type(4)  = 'TcO4-'
       ion_type(5)  = 'HTcO4'
       ion_type(6)  = 'SO4--'
       ion_type(7)  = 'HSO4-'
       ion_type(8)  = 'H2SO4'
     
     end if

     ! abstrial is used to ensure charge neutrality
     abstrial = 0.0d0
     do i = 1, nspec
       if (zchge(i).ne.0.0d0) then
         if (bulkcO(i).gt.0.1) then
           bulkcO(i) = 0.001
         end if
         abstrial = abstrial + abs(zchge(i) * bulkcO(i))
       end if
     end do
!    ---------------------------------------------------------------------------

!    ------------------- CALCULATE BULK CONCENTRATIONS -------------------------

     errrel = relerr
     allocate(derivc(nspec,nspec),trialc(nspec),func(nspec))

     derivc(:,:) = 0.0d0
     func(:)     = 0.0d0
     trialc(:)   = 0.0d0
 
 100 continue

     ! Begin Broydon secant method with analytic Jacobian.
     call hybrj1(equilibria,nspec,bulkcO,func,derivc,nspec,errrel,info)

     ! -- DESCRIPTION : hybrj1 -------------------------------------------------
     ! Minpack subroutine finds a zero of the system of nspec non-linear
     ! functions in nspec variables by a modification of the Powell hybrid
     ! method, using a user supplied subroutine which calculates the functions
     ! and the Jacobian.      
     !
     ! -- VARIABLES ------------------------------------------------------------
     ! equilibria - user supplied subroutine which calculates the function
     !              and the Jacobian
     !      nspec - the number of functions and variables
     !     bulkcO - on input, contains the initial estimate of solution
     !              on output, contains the final estimate of the solution
     !       func - the functions evaluated at the output bulkcO
     !     derivc - an nspec by nspec array containing the orthonormal
     !              matrix Q produced by the QR factorization of the final
     !              approximate Jacobian
     !      nspec - leading dimension of derivc
     !     errrel - termination occurs when the algorithm estimates that
     !              the relative error between bulkcO and the solution is
     !              at most relerr
     !       info - output integer error flag: 
     !              0 - improper input parameters
     !              1 - relative error is less than relerr 
     !              2 - maximum number of iterations exceeded
     !              3 - relerr is too small
     !              4 - iteration is making slow progress, consult MINPACK 
     !                  manual
     ! ------------------------------------------------------------------------

     ! flag messages from hybrj1 subroutine
     if (info.eq.0) then

       write(5,*) 'Input Error - check code'
       stop

     else if (info.eq.1) then 

       write(5,*)
       write(5,*) '------- EQUILIBRIUM CONCENTRATIONS ------'
       write(5,*) 'SPECIES           ', 'CONCENTRATION (mol/l)'
       do i = 1, nspec

         ! bulk concentrations after equilibration
         bulkc(i) = bulkcO(i)

         write(5,*) ion_type(i), bulkc(i)

       end do

     else if (info.eq.2) then

       write(5,*) 'Maximum number of iterations exceeded'
       stop

     else if (info.eq.3) then !recompute bulk concentrations with new relerr
      
       write(5,*) 'Toleration set too small - value is doubled'
       errrel = errrel * 2.0d0
       goto 100

     else if (info.eq.4) then

       write(5,*) 'Error - consult MINPACK manual'
       stop

     end if

!    ---------------------------------------------------------------------------


!    --------------------- CALCULATE REQUIRED CONSTANTS ------------------------    
     
     write(5,*)
     write(5,*) '--------------- CONSTANTS --------------'

     ! set activity coefficients to 1 for first pass
     gamma(:) = 1.0d0

     ! calculate the Bjerrum length
     bj_len = e2/eps0/eps/kb/T                                   ! m
     write(5,*) 'Bjerrum length (nm) =', bj_len * 1.0d9 / fourpi ! nm

     ! concentration of surface charge sites
     charge_dens = charge_dens * 1.0d18                          ! m**-2
    
     ! area per unit charge
     area = 1 / charge_dens                                      ! m**2

     ! surface density of charges
     b0 = bj_len / 2.0d0 / avogad / area**2                      ! mol m**-3
     write(5,*) 'Volume density of charges on surface (m**-3)', b0

     bj_len = bj_len *  avogad                                   ! mol**-1 m

     ! calculate Debye Huckel screening length
     sum_debye = 0.0d0
     do i = 1, nspec
        sum_debye = sum_debye + bulkc(i) * zchge(i) * zchge(i)   ! mol dm**-3
     end do
     dh_len = 1 / sqrt(bj_len * sum_debye)                       
     write(5,*) 'Debye-Huckel screening length (nm) =', dh_len * 1.0d9   ! nm

     ! calculate distance interval between points
     dist0 = (kb*T/e)*sqrt(eps*eps0/2.0/kb/T/avogad)             ! ??     

     ! calculate ionic strength
     ion_strength = 0.0d0
     do i = 1, nspec
       ion_strength =  ion_strength + bulkc(i) * zchge(i) * zchge(i)
     end do
     ion_strength = ion_strength / 2.0d0                         ! mol dm**-3
     write(5,*) 'Ionic strength of electrolyte (mol/l) =', ion_strength
     write(5,*)
 
!    ---------------------------------------------------------------------------      


!    ---------------------- SOLVE THE GRAHAME EQUATION -------------------------
!    Get the roots of the nth order polynomial in P0 = exp (-e*psi/kbT).
!    Solve the set of equations self-consistently.
!    On first pass assume activity coefficients are set to unity.

     ! allocate arrays
     allocate (zconc(npoint,nspec),charge(npoint),distance(npoint), &
   & potential(npoint),err(npoint),splwork(2*npoint),splderiv(npoint)) 

 500 continue
     print*, gamma
!    1. Obtain the coefficients for the polynomial
     call polynom

!    2. Find the zeros of the polynomial
     call zroots(polycoeff,nroot,roots,lpolish)

     ! -- DESCRIPTION : zroots -------------------------------------------------
     ! Subroutine to find the zero of a polynomial with real coefficients using
     ! Laguerre's method.
     !
     ! Uses laguer subroutine
     !
     ! -- VARIABLES ------------------------------------------------------------
     ! polycoeff - vector containing nroot + 1 coefficients of the polynomial
     !             in increasing order by degree
     !     nroot - the number of roots 
     !     roots - complex vector containing the polynomial zeros
     !   lpolish - set lpolish =  .true. if polishing is required
     ! ------------------------------------------------------------------------

!    3. Choose the positive real roots 
     froot = 0
     do i = 1, nroot 

       ! if imaginary part is much larger than zero not a positive real root
       if (abs(aimag(roots(i))).gt.1.0d-4) then
         goto 200
       end if
       
       realpart = roots(i)

       ! if realpart is negative then not a positive real root
       if (realpart.lt.0.0d0) then 
         goto 200
       end if

       ! number of positive real roots
       froot = froot + 1
       realroot(froot) = realpart

200    continue
 
     end do

!    Choose between the real roots by requiring that the charge on the
!    surface is equal to the charge in solution.

!    4. Calculate surface and solution charges for each root (in units of e)
     do i = 1, froot
        
       ! potential at surface
       P0 = realroot(i)

       ! calculate charge at surface
       call surf_charge

       ! calculate charge in solution
       call sol_charge

       ! go on to next root if distances are negative
       if (wrong_root) then
         goto 700
       end if

       ! surface and solution charges are balanced
       if (abs(sigma_surf + sigma_sol).lt.0.01) then
         write(5,*) '-------------------- CHARGE-BALANCING ---------------------'   
         write(5,*)
         write(5,*) 'Charge on the surface = ', sigma_surf
         write(5,*) 'Charge in solution = ', sigma_sol
         write(5,*) 'Difference between charges = ', sigma_sol + sigma_surf
         write(5,*)
         goto 300

       else
         
         write(5,*) i, sigma_surf, sigma_sol

       end if
              
700    continue

     end do

     write(5,*) 'ERROR - no consistent solution found.'
     stop
     
300  continue

     if (lactive) then
!      5. Redetermine activity coefficients using TJ equation
       call activity
       lallocate = .false.

       ! if all activity coefficients have converged then proceed
       do i = 1, nspec
         if (convergence(i)) then 
           continue
         ! if not recalculate the Grahame equation using the new coefficients
         else
           goto 500
         end if
       end do
     end if

!    ---------------------------------------------------------------------------

!    -------------- CALCULATE & PRINT CHARGE-DISTANCE PROFILES -----------------

     write(5,*)
     write(5,*) '---------------------- RESULTS -----------------------'
     write(5,*)
     write(5,*) 'FINAL CONVERGED ACTIVITY COEFFICIENTS:'
     write(5,*)
     do i = 1, nspec
       write(5,*) ion_type(i),gamma(i)
     end do
     write(5,*)

     ! calculate solution charge at the surface
     charge(1) = 0.0
     do i = 1, nspec

       zconc(1,i) = bulkc(i) * P0**zchge(i)
       charge(1) = charge(1) + zchge(i) * zconc(1,i)

     end do

     ! potential at surface
     potential(1) = -kb * T * log(P0) / e                           ! V

     write(5,*)
     write(5,*) 'Surface potential (mV) = ', potential(1) * 1.0d3
     write(5,*) 'pH at surface = ', -log(bulkc(1) * P0)
     write(5,*) 'K+ surface conc. (mol/l) = ', zconc(1,3)

     if (switch.ne.3) then
 
       write(5,*) 'Fe+++ surface conc. (mol/l) = ', zconc(1,9)
 
     end if

     write(5,*) 'Fraction of monolayer sites covered with H+ = ', fractH
     write(5,*) 'Fraction of monolayer sites covered with K+ = ', fractK
 
     if (switch.ne.3) then
 
       write(5,*) 'Fraction of monolayer sites covered with Fe+++ = ', fractFe
 
     end if
 
     write(5,*) 'Fraction of ionised monolayer sites = ', fract
    
     ! calculate concentration profiles by numerical integration of the field
     write(5,*)
     write(5,*) '------------ CHARGE-DISTANCE PROFILES ------------'
     write(5,*) 
     write(5,*) 'distance   ', 'charge density   ', 'potential   ', (ion_type(i),i=1,nspec)
     write(5,*) 'units: ','  nm  ','  mol/l    ','    mV    ','    mol/l    '
     write(5,*) distance(1)*1.0d9,  charge(1), potential(1) * 1.0d3, (zconc(1,i),i=1,nspec)

     ! calculate distance, charge, potential and ionic concentrations at npoints
     do i = 2, npoint

       P = P0 - float(i - 1) * (P0 - 1.0) / float(npoint)

       call qag(integrand,P,P0,errabs,1.0e-4,irule,distance(i),err(i),neval,ier)

       if ((ier.eq.1).or.(ier.eq.3)) then
         write(5,*) 'ERROR - consult qag routine header'
       else if (ier.eq.2) then
         write(5,*) 'ERROR - roundoff error prevents requested accurcay'
       else if (ier.eq.6) then
         write(5,*) 'ERROR - invalid input, errabs =', errabs
       end if  

       distance(i) = distance(i) * dist0                         ! m
       charge(i) = 0.0

       do j = 1, nspec
         
         zconc(i,j)   = bulkc(j) * P**zchge(j)                   ! mol dm**-3
         charge(i)    = charge(i) + zchge(j) * zconc(i,j)        ! mol dm**-3    
         potential(i) = -kb * T * 1.0d3 * log(P) / e             ! mV

       end do

       write(5,*) distance(i)*1.0d9, charge(i), potential(i), (zconc(i,j),j=1,nspec)     

     end do

     write(5,*)
     write(5,*) '------------------- END -----------------------'

!    ***************************************************************************
!    ***************************************************************************
     END PROGRAM

