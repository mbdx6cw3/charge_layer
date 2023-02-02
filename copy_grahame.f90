     MODULE grahame_module
!    ---------------------------------------------------------------------------
!    Module containing subroutines called from the charge_layer.f90 program
!    used to solve the Grahame equation.

     CONTAINS
     SUBROUTINE polynom
     USE header_file, only: Ka,KK,KFe,gamma,bulkc,polycoeff,&
   & roots,realroot,nroot,switch,double,b0,lallocate
!    ***************************************************************************
     implicit none
     real(kind=double) :: k1, k2, k3
     real(kind=double) :: a1, a2, a3
     real(kind=double) :: d0, d1, d2, d3, d4
     real(kind=double) :: b2, b3
     integer i

     k1 = Ka * gamma(1) * bulkc(1)
     k1 = k1 + KK * gamma(3) * bulkc(3)

     if (switch.eq.1) then
      !  nroot = 1
      ! k2 = KFe * gammaFe * bulkc(9)
      ! k3 = KK * gammaK * bulkc(3)

        ! determine polynomial coefficients
       ! polycoeff(1) = cmplx(d0                      , 0.0d0)
       ! polycoeff(2) = cmplx(d1 + d0*a1              , 0.0d0)
       ! polycoeff(3) = cmplx(d2 + d1*a1 + d0*a2 - b2 , 0.0d0)
       ! polycoeff(4) = cmplx(d3 + d2*a1 + d1*a2      , 0.0d0)
       ! polycoeff(5) = cmplx(d4     d3*a1 + d2*a2      , 0.0d0)
       ! polycoeff(6) = cmplx(d5             d3*a2      , 0.0d0)

     else if (switch.eq.2) then
     !   nroot = 1
     !  k2 = KFe * gammeFe * bulkc(3)
     !  k3 = KK * gammaK * * bulkc(5)
     else if (switch.eq.3) then

        nroot = 5
        if (lallocate) then
          allocate(polycoeff(nroot+1),roots(nroot),realroot(nroot))
          polycoeff(:) = 0.0d0
          roots(:)     = 0.0d0
          realroot(:)  = 0.0d0
        end if

        a1 = 2.0d0 * k1
        a2 = k1 * k1
        b2 = b0

        d0 = bulkc(6)
        d1 = bulkc(2) + bulkc(4) + bulkc(7)
        d3 = bulkc(1) + bulkc(3)
        d2 = -(d0 + d1 + d3)
        
        ! determine polynomial coefficients
        polycoeff(1) = cmplx(d0                      , 0.0d0)
        polycoeff(2) = cmplx(d1 + d0*a1              , 0.0d0)
        polycoeff(3) = cmplx(d2 + d1*a1 + d0*a2 - b2 , 0.0d0)
        polycoeff(4) = cmplx(d3 + d2*a1 + d1*a2      , 0.0d0)
        polycoeff(5) = cmplx(     d3*a1 + d2*a2      , 0.0d0)
        polycoeff(6) = cmplx(             d3*a2      , 0.0d0)

      end if
!    ***************************************************************************
     END SUBROUTINE polynom
!    ---------------------------------------------------------------------------
     
     SUBROUTINE zroots(a,m,roots,polish)
!    ***************************************************************************
      IMPLICIT none
      INTEGER m, maxm
      DOUBLE PRECISION eps
      COMPLEX*16 a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (eps=1.0D-6, maxm=101)
!
!   Uses laguer
!
!   Given the degree m and the complex coefficients a(1:m+1) of a polynomial
!   calls laguer to find all m complex routes. 
!   Set polish = .true. if polishing is desired
!
      INTEGER i,j,jj,its
      COMPLEX*16 ad(maxm), x,b,c
!
      DO j=1, m+1
         ad(j) = a(j)
      ENDDO
      DO j=m,1,-1
         x = DCMPLX(0.0D0,0.0D0)
         CALL laguer(ad,j,x,its)
         IF(ABS(AIMAG(x)).LE.2.0D0*eps**2*ABS(REAL(x))) x=DCMPLX(REAL(x),0.0D0)
         roots(j) = x
         b = ad(j+1)
         DO jj=j,1,-1
            c = ad(jj)
            ad(jj) = b
            b = x*b+c
         ENDDO
      ENDDO
      IF (polish) THEN
         DO j=1,m
            CALL laguer(a,m,roots(j),its)
         ENDDO
      ENDIF
      DO j=2,m
         x = roots(j)
         DO i=j-1,1,-1
            IF(REAl(roots(i)).LE.REAL(x)) GOTO 10
            roots(i+1) = roots(i)
         ENDDO
         i = 0
10       roots(i+1) = x
      ENDDO
      RETURN
      END SUBROUTINE zroots
!    -----------------------------------------------------------------------------

      SUBROUTINE laguer(a,m,x,its)
!     ****************************************************************************
      IMPLICIT none
      INTEGER m,its,maxit,mr,mt
      DOUBLE PRECISION epss
      COMPLEX*16 a(m+1), x
      PARAMETER (epss=2.0D-7, mr=8, mt=10, maxit=mr*mt)
!
!       Given the degree m and the complex coefficients a of a polynomial, 
!       this routine improves the solution x by Laguerre's method until it
!       converges within achievable roundoff error
!
      INTEGER iter, j
      DOUBLE PRECISION abx,abp,abm,err,frac(mr)
      COMPLEX*16 dx, x1, b, d, f, g, h, sq, gp, gm, g2
      SAVE frac
      DATA frac /0.5D0, 0.25D0, 0.75D0, 0.13D0, 0.38D0, 0.62D0, 0.88D0, 1.0D0/
!
      DO iter=1,maxit
         its = iter
         b = a(m+1)
         err = ABS(b)
         d = DCMPLX(0.0D0,0.0D0)
         f = DCMPLX(0.0D0,0.0D0)
         abx = ABS(x)
         DO j=m,1,-1
            f = x*f+d
            d = x*d+b
            b = x*b+a(j)
            err = ABS(b)+abx*err
         ENDDO
         err = epss*err
         IF (ABS(b).LE.err) THEN
            RETURN
         ELSE
            g = d/b
            g2 = g*g
            h = g2 - 2.0D0*f/b
            sq = SQRT((m-1)*(m*h-g2))
            gp = g+sq
            gm = g-sq
            abp = ABS(gp)
            abm = ABS(gm)
            IF (abp.LT.abm) gp=gm
            IF (MAX(abp,abm).GT.0.0D0) THEN
               dx = m/gp
            ELSE
               dx = EXP(DCMPLX(LOG(1.0D0+abx),DFLOAT(iter)))
            ENDIF
         ENDIF
         x1 = x-dx
         IF (x.EQ.x1) RETURN
         IF (MOD(iter,mt).NE.0.0D0) THEN
            x = x1
         ELSE
            x = x -dx*frac(iter/mt)
         ENDIF
      ENDDO
      WRITE (6,1000)
1000  FORMAT (//6X,'Too many iterations in laguer - try different starting point')
      STOP
      END SUBROUTINE laguer
!    ---------------------------------------------------------------------------


     SUBROUTINE surf_charge
     USE header_file, only: Ka,KK,KFe,gamma,bulkc,zchge,P0,&
   & fractH,fractK,fractFe,fract,sigma_surf,switch,double
!    ***************************************************************************
     implicit none
     real(kind=double) :: surfH, surfK, surfFe
     real(kind=double) :: totion     
     
     ! surface concentration of cations
     surfH = Ka * gamma(1) * P0**zchge(1) * bulkc(1)
     surfK = KK * gamma(3) * P0**zchge(3) * bulkc(3)

     ! sigma_surf is the charge (units of e) on the surface
     ! fract is the fraction of surface sites which are ionised
     if (switch.ne.3) then
       surfFe = KFe * gamma(9) * P0**zchge(9) * bulkc(9)
       totion = 1.0 + surfH + surfK + surfFe
       fractH = surfH / totion
       fractK = surfK / totion
       fractFe = surfFe / totion
       fract = 1 / totion
       sigma_surf = 2.0d0 * fractFe - fract     !due to 3+ charge of Fe
     else
       totion = 1.0 + surfH + surfK
       fractH = surfH / totion
       fractK = surfK / totion
       fract = 1 / totion
       sigma_surf = -fract  
     end if
!    ***************************************************************************
     END SUBROUTINE surf_charge
!    ---------------------------------------------------------------------------

     SUBROUTINE sol_charge
     USE header_file, only: charge,bulkc,zchge,errrel,sigma_sol,distance,P0,&
   & roots,realroot,nroot,integrand,pchqa,P,errabs,irule,err,neval,ier,&
   & wrong_root,dist0,nspec,npoint,splderiv,lspline,splwork,avogad,sum_debye,&
   & eps,eps0,area,kb,T,e2,dh_len
!    ***************************************************************************
!    Calculate charge in solution by first calculating the charge profile and
!    then integrate to obtain total charge.
     implicit none
     integer :: i,j
     errrel      = 1.0e-6
     sigma_sol   = 0.0d0
     charge(1)   = 0.0d0
     distance(1) = 0.0d0

     ! solution concentrations at the surface
     do i = 1, nspec
       charge(1) = charge(1) + zchge(i) * bulkc(i) * P0**zchge(i)
     end do
      
     do i = 2, npoint

       P = P0 - float(i - 1) * (P0 - 1.0) / float (npoint)

       ! approximate distances
       call qag(integrand,P,P0,errabs,1.0e-4,irule,distance(i),err(i),neval,ier)
 
       ! -- DESCRIPTION : qag ----------------------------------------------------
       ! Quadpack subroutine to approximate an integral (distance) over a finite 
       ! definite integral. I = integral of F over (P,P0).
       !
       ! -- VARIABLES ------------------------------------------------------------
       ! integrand - the name of the function routine, provided by user
       !     P, P0 - the limits of the integration
       !    errabs - absolute accuracy required
       !    errrel - relative accuracy required
       !     irule - chooses the order of the local integration rule
       !  distance - the estimated value of the integral
       !    err(i) - an estimate of the |integral - distance| 
       !     neval - the number of times the integral was evaluated
       !       ier - return flag:
       !               0 - normal and reliable termination, accurcay achieved
       !               1 - maximum number of subdivisions achieved
       !               2 - occurence of roundoff error is detected, preventing
       !                   requested tolerance from being achieved
       !               3 - extremely bad integrand behaviour occurs
       !               6 - the input is invalid (errabs<0 and errrel<0)
       ! -------------------------------------------------------------------------

       ! if the distance is negative then we have the wrong root 
       if (distance(i).lt.0.0d0) then
         wrong_root = .true.
         ! must exit here as pchez will not accept decreasing distance
         go to 600
       else
         wrong_root = .false.
       end if

       if ((ier.eq.1).or.(ier.eq.3)) then
         write(5,*) 'ERROR - consult qag routine header'
       else if (ier.eq.2) then
         write(5,*) 'ERROR - roundoff error prevents requested accurcay'
       else if (ier.eq.6) then
         write(5,*) 'ERROR - invalid input, errabs =', errabs 
       end if  

       distance(i) = distance(i) * dist0
       charge(i) = 0.0
       do j = 1, nspec 
         charge(i) = charge(i) + zchge(j) * bulkc(j) * P**zchge(j)
       end do
     end do

     ! fit the result to a set of cubic splines
     call pchez(npoint,distance,charge,splderiv,lspline,splwork,2*npoint,ier)

     ! -- DESCRIPTION : pchez --------------------------------------------------
     ! NMS subroutine that carries out spline or cubic Hermite interpolation.
     !
     ! -- VARIABLES ------------------------------------------------------------ 
     !    npoint - number of data points (>2)
     !  distance - strictly increasing independent variable values
     !    charge - function values, corresponding to the independent variable
     !  splderiv - the derivative values at the data points
     !   lspline - logical variable which specifies whether the interpolant is
     !             to be a spline with two continuous derivatives (.true.), or
     !             a Hermite cubic interpolant with one continuous derivative
     !             (.false.)
     !   splwork - array with dimension 2*npoint, only required if lspline is
     !             true
     !       ier - error flag:
     !               0 - no errors
     !               1 - N<2
     !               3 - if the distance array is not strictly increasing
     !               7 - if 2*npoint is less than 2*npoint. This won't happen!
     ! -------------------------------------------------------------------------

     sigma_sol = pchqa(npoint,distance,charge,splderiv,distance(1),distance(npoint))
     sigma_sol = sigma_sol * area * avogad
     sigma_sol = sigma_sol + eps*eps0*area*kb*T*log(P)/e2/dh_len

 600 continue

!    ***************************************************************************
     END SUBROUTINE sol_charge
!    ---------------------------------------------------------------------------


     SUBROUTINE activity
     USE header_file, only: bulkc,zchge,P0,ion_strength,convergence,double,&
    & nspec,gamma,gamma_old,loggam,keyword2
!    ***************************************************************************
!    Determine activity coefficients, gamma, set to unity on first pass,
!    assuming a dilute solution.
! 
!    After the first pass use the Truesdell-Jones equation to determine the
!    activity coefficients, which is used to renormalise pKa values, recalculate
!    the Grahame equation until the pKa values are stable to within a factor of
!    two.
!    ***************************************************************************
     implicit none
     integer :: i

     ! save old activity coefficients
     gamma_old(:) = gamma(:)

     ! calculate ionic strength at the surface
     ion_strength = 0.0d0
     do i = 1, nspec
       ion_strength = ion_strength + 0.5 * bulkc(i) * zchge(i)**2 * P0**zchge(i)
     end do

     ! if ionic strength is greater than something use the Pitzer model
     if (keyword2.eq.'pitzer') then

       ! get activity coefficients
       call pitzer

     ! otherwise use the Truesdell-Jones equation
     else

       ! get activity coefficients
       call truesdell

     end if

     ! loop over all ions to check their activity coefficients have converged
     do i = 1, nspec

       gamma(i) = 10.0d0 ** (loggam(i))      

       if(abs(gamma_old(i)-gamma(i)).lt.1.0d-8) then
         convergence(i) = .true.
       else
         convergence(i) = .false.
       end if

     end do
!    ***************************************************************************
     END SUBROUTINE activity
!    ---------------------------------------------------------------------------


     SUBROUTINE truesdell
     USE header_file, only: ion_strength,double,loggam,zchge,atrues,btrues,nspec
!    ***************************************************************************
     implicit none
     real(kind=double), parameter :: a_dh = 0.5085, b_dh = 0.3281 
     integer :: i    

     do i = 1, nspec
       loggam(i) = -a_dh * zchge(i) * zchge(i) * sqrt(ion_strength) &
      &          /(1.0d0 + b_dh * atrues(i) * sqrt(ion_strength)) &
      &          + btrues(i) * ion_strength
     end do

!    ***************************************************************************
     END SUBROUTINE truesdell
!    ---------------------------------------------------------------------------


     SUBROUTINE pitzer
     USE header_file, only: ion_strength,double,loggam,zchge,pitzer_b0, &
    &  pitzer_b1,pitzer_C,pitzer_phi,pitzer_lamda,pitzer_psi,nspec,P0,bulkc
!    ***************************************************************************
     implicit none
     integer :: i,j,k
     real(kind=double) :: alpha, alpha1, gam_sum, zchge2, z_val, dh_lim
     real(kind=double) :: pitzer_ba, pitzer_bb, gfunc

     z_val = 0.0d0
     do i = 1, nspec
       z_val = z_val + bulkc(i) * zchge(i) * P0 * zchge(i)
     end do

     ! factor of 2 since all interacting species are 1:1-4 electrolyte
     alpha = 2.0d0 * sqrt(ion_strength)

     do i = 1, nspec

       gam_sum = 0.0d0

       zchge2 = zchge(i) * zchge(i)

       dh_lim = -0.392*(sqrt(ion_strength)/(1+1.2*sqrt(ion_strength)) &
      &           + 1.667*dlog(1+1.2*sqrt(ion_strength)))
       gam_sum = zchge2 * dh_lim

       ! compute activity coefficient of cations
       if (zchge(i).gt.0.0d0) then

         ! loop over all ions
         do j = 1, nspec

           ! anionic species
           if (zchge(j).lt.0.0d0) then

             gfunc = 2.0d0*(1.0-(1.0+alpha)*exp(-alpha))/alpha/alpha
             pitzer_ba = pitzer_b0(i,j) + pitzer_b1(i,j) * gfunc
             gam_sum = gam_sum + bulkc(j) * (2.0d0 * pitzer_ba + z_val * pitzer_c(i,j))

             ! inner cation loop
             do k = 1, nspec

               if (zchge(k).gt.0.0d0) then
                 gfunc = exp(-alpha) - gfunc
                 pitzer_bb = pitzer_b1(j,k) * gfunc
                 gam_sum = gam_sum + abs(zchge(i)) * bulkc(j) * bulkc(k) * pitzer_c(j,k)
                 gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * &
               &                (pitzer_bb + z_val * pitzer_c(j,k) / 2.0)
                 gam_sum = gam_sum + bulkc(k) * bulkc(j) * pitzer_psi(i,j,k)
               ! inner anion loop
               else if (zchge(k).lt.0.0d0) then
                   gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * pitzer_phi(j,k)
                   gam_sum = gam_sum + bulkc(j) * bulkc(k) * pitzer_psi(i,j,k)
               end if

             end do

           ! cationic species
           else if (zchge(j).gt.0.0d0) then

             gam_sum = gam_sum + 2.0d0 * bulkc(j) * pitzer_phi(i,j)
       
             ! inner cation loop
             do k = 1, nspec
               if (zchge(k).gt.0.0d0) then
                 gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * pitzer_phi(j,k)
               end if
             end do

           ! neutral species
           else

             gam_sum = gam_sum + 2.0d0 * bulkc(j) * pitzer_lamda(i,j)
       
           end if

         end do

       ! computing the activity coefficient of anions
       else if (zchge(i).lt.0.0d0) then

         ! loop over all ions
         do j = 1, nspec

           ! cationic species
           if (zchge(j).gt.0.0d0) then

             gfunc = 2.0d0*(1.0-(1.0+alpha)*exp(-alpha))/alpha/alpha
             pitzer_ba = pitzer_b0(i,j) + pitzer_b1(i,j) * gfunc
             gam_sum = gam_sum + bulkc(j) * (2.0d0 * pitzer_ba + z_val * pitzer_c(i,j))

             ! inner anion loop
             do k = 1, nspec

               if (zchge(k).lt.0.0d0) then
                 gfunc = exp(-alpha) - gfunc
                 pitzer_bb = pitzer_b1(j,k) * gfunc
                 gam_sum = gam_sum + abs(zchge(i)) * bulkc(j) * bulkc(k) * pitzer_c(j,k)
                 gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * &
               &                (pitzer_bb + z_val * pitzer_c(j,k) / 2.0)
                 gam_sum = gam_sum + bulkc(k) * bulkc(j) * pitzer_psi(i,j,k)
               ! inner cation loop
               else if (zchge(k).gt.0.0d0) then
                 gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * pitzer_phi(j,k)
                 gam_sum = gam_sum + bulkc(j) * bulkc(k) * pitzer_psi(i,j,k)
               end if

             end do

           ! anionic species
           else if (zchge(j).lt.0.0d0) then

             gam_sum = gam_sum + 2.0d0 * bulkc(j) * pitzer_phi(i,j)
             ! inner anion loop
             do k = 1, nspec

               if (zchge(k).lt.0.0d0) then
                 gam_sum = gam_sum + zchge2 * bulkc(j) * bulkc(k) * pitzer_phi(j,k)
               end if

             end do
           ! neutral species
           else
             gam_sum = gam_sum + 2.0d0 * bulkc(j) * pitzer_lamda(i,j)
           end if

         end do

       ! computing the activity coefficient of a neutral species
       else

         do j = 1, nspec
           gam_sum = gam_sum * bulkc(j) * pitzer_lamda(i,j)
         end do
         gam_sum = gam_sum * 2.0d0

       end if

       loggam(i) = gam_sum

     end do 
!    *************************************************************************** 
     END SUBROUTINE pitzer
!    ---------------------------------------------------------------------------


     END MODULE grahame_module
