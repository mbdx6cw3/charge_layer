subroutine pchez ( n, x, f, d, spline, wk, lwk, ierr )
!
!*******************************************************************************
!
!! PCHEZ carries out easy to use spline or cubic Hermite interpolation.
!
!
!  Discussion:
!
!    This routine sets derivatives for spline (two continuous derivatives)
!    or Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data. 
!
!    This routine is an easy to use driver for the PCHIP routines.
!    Various boundary conditions are set to default values by PCHEZ.  
!    Many other choices are available in the subroutines PCHIC, 
!    PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!    If SPLINE is TRUE, the interpolating spline satisfies the default 
!    "not-a-knot" boundary condition, with a continuous third derivative 
!    at X(2) and X(N-1). 
!
!    If SPLINE is FALSE, the interpolating Hermite cubic will be monotone 
!    if the input data is monotone.  Boundary conditions are computed from 
!    the derivative of a local quadratic unless this alters monotonicity.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    Fred Fritsch and R Carlson, 
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch and J Butland,
!    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
!    LLNL Preprint UCRL-87559, April 1982.
!
!    Carl de Boor,
!    A Practical Guide to Splines, Chapter IV,
!    Springer-Verlag,
!    New York, 1978.
!
!    Fred Fritsch, 
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications, 
!    Lawrence Livermore National Laboratory, 
!    Computer Documentation UCID-30194, August 1982.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(N), the function values.  F(I) is the value corresponding 
!    to X(I).
!
!    Output, real D(N), the derivative values at the data points.
!
!    Input, logical SPLINE, specifies if the interpolant is to be a spline 
!    with two continuous derivaties (SPLINE is TRUE), or a Hermite cubic 
!    interpolant with one continuous derivative (SPLINE is FALSE).
!
!    Workspace, real WK(LWK), required only if SPLINE is TRUE.
!
!    Input, integer LWK, the length of the work array WK, which must
!    be at least 2*N.  However, WK is not needed if SPLINE is FALSE,
!    and in this case LWK is arbitrary.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, can only occur when SPLINE is FALSE,  means that there were
!      IERR switches in the direction of monotonicity.  When SPLINE is 
!      FALSE, PCHEZ guarantees that if the input data is monotone, the 
!      interpolant will be too.  This warning is to alert you to the fact 
!      that the input data was not monotone.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -7, if LWK is less than 2*N and SPLINE is TRUE.
!
  implicit none
!
  integer lwk
  integer n
!
  real d(n)
  real f(n)
  integer, save, dimension ( 2 ) :: ic = (/ 0, 0 /)
  integer ierr
  integer, parameter :: incfd = 1
  logical spline
  real vc(2)
  real wk(lwk)
  real x(n)
  integer i
!
  if ( spline ) then
    call  pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    call  pchim ( n, x, f, d, incfd, ierr )
  end if

  return
end


subroutine pchsp ( ic, vc, n, x, f, d, incfd, wk, nwk, ierr )
!
!*******************************************************************************
!
!! PCHSP sets derivatives for Hermite representation of cubic spline interpolant.
!
!
!  Description:
!
!    PCHSP sets derivatives needed to determine the Hermite representation 
!    of the cubic spline interpolant to given data, with specified boundary 
!    conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of Carl de Boor's cubic spline routine CUBSPL.
!
!  Reference:  
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978).
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    0, to set D(1) so that the third derivative is continuous at X(2).  
!      This is the "not a knot" condition provided by de Boor's cubic spline 
!      routine CUBSPL, and is the default boundary condition here.
!    1, if first derivative at X(1) is given in VC(1).
!    2, if second derivative at X(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 4.
!    For the "natural" boundary condition, use ibeg=2 and vc(1)=0.
!    IC(2) = IEND, desired condition at end of data.
!    IEND may take on the same values as IBEG, but applied to derivative at
!    X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    Input, real VC(2), specifies desired boundary values, as indicated above.
!    VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(INCFD,N), the dependent values to be interpolated.  
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real D(INCFD,N), the derivative values at the data points.
!    These values will determine the cubic spline interpolant with the 
!    requested boundary conditions.  The value corresponding to X(I) is 
!    stored in D(1+(I-1)*INCFD).
!
!    Input, integer INCFD, increment between successive values in F and D.
!
!    Workspace, real WK(NWK).
!
!    Input, integer NWK, the size of WK, which must be at least 2 * N.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IBEG < 0 or IBEG > 4.
!    -5, if IEND < 0 or IEND > 4.
!    -6, if both of the above are true.
!    -7, if NWK is too small.
!    -8, in case of trouble solving the linear system
!        for the interior derivative values.
!
  implicit none
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real f(incfd,n)
  real g
  integer ibeg
  integer ic(2)
  integer iend
  integer ierr
  integer index
  integer j
  integer nwk
  real pchdf
  real stemp(3)
  real vc(2)
  real wk(2,n)
  real x(n)
  real xtemp(4)
!
  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchsp -- number of data points less than two', 44, ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchsp -- increment less than one', 32, ierr, 1)
    return
  end if

  do j = 2, n
    if ( x(j) <= x(j-1) ) then
      ierr = -3
      call xerror ('pchsp -- x-array not strictly increasing', 40, ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
  if ( ibeg < 0 .or. ibeg > 4 )  ierr = ierr - 1
  if ( iend < 0 .or. iend > 4 )  ierr = ierr - 2
  if ( ierr < 0 ) go to 5004
!
!  Function definition is ok -- go on.
!
  if ( nwk < 2 * n )  go to 5007
!
!  Compute first differences of X sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j = 2, n
     wk(1,j) = x(j) - x(j-1)
     wk(2,j) = ( f(1,j) - f(1,j-1) ) / wk(1,j)
  end do
!
!  Set to default boundary conditions if N is too small.
!
  if ( ibeg > n )  ibeg = 0
  if ( iend > n )  iend = 0
!
!  Set up for boundary conditions.
!
  if ( ibeg == 1 .or. ibeg == 2 )  then
     d(1,1) = vc(1)
  else if ( ibeg > 2 )  then
!
!  Pick up first IBEG points, in reverse order.
!
     do j = 1, ibeg
       index = ibeg-j+1
       xtemp(j) = x(index)
       if ( j < ibeg ) then
         stemp(j) = wk(2,index)
       end if
     end do

     d(1,1) = pchdf ( ibeg, xtemp, stemp, ierr )
     if ( ierr /= 0 ) then
       go to 5009
     end if

     ibeg = 1
  end if

  if ( iend == 1 .or. iend == 2 )  then
     d(1,n) = vc(2)
  else if ( iend > 2 )  then
!
!  Pick up last IEND points.
!
     do j = 1, iend
       index = n - iend + j
       xtemp(j) = x(index)
       if ( j < iend ) then
         stemp(j) = wk(2,index+1)
       end if
     end do

     d(1,n) = pchdf ( iend, xtemp, stemp, ierr )

     if ( ierr /= 0 ) then
       go to 5009
     end if

     iend = 1

  end if
!
!  Begin coding from cubspl
!
!  A tridiagonal linear system for the unknown slopes S(1:N) of
!  F at X(1:N) is generated and then solved by Gauss elimination, 
!  with s(j) ending up in d(1,j), all j.
!  wk(1,.) and wk(2,.) are used for temporary storage.
!
!  Construct first equation from first boundary condition, of the form
!    wk(2,1) * s(1) + wk(1,1) * s(2) = D(1,1)
!
  if ( ibeg == 0 )  then
     if ( n == 2 )  then
!
!  No condition at left end and N = 2.
!
        wk(2,1) = 1.0E+00
        wk(1,1) = 1.0E+00
        d(1,1) = 2.0E+00 * wk(2,2)
     else
!
!  Not-a-knot condition at left end and N > 2.
!
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =(( wk(1,2) + 2.0E+00 * wk(1,1) ) * wk(2,2) * wk(1,3) &
                             + wk(1,2)**2 * wk(2,3)) / wk(1,1)
     end if
  else if ( ibeg == 1 ) then
!
!  Slope prescribed at left end.
!
     wk(2,1) = 1.0E+00
     wk(1,1) = 0.0E+00
  else
!
!  Second derivative prescribed at left end.
!
     wk(2,1) = 2.0E+00
     wk(1,1) = 1.0E+00
     d(1,1) = 3.0E+00 * wk(2,2) - 0.5E+00 * wk(1,2) * d(1,1)
  end if
!
!  If there are interior knots, generate the corresponding equations and
!  carry out the forward pass of Gauss elimination, after which the J-th
!  equation reads    
!
!    wk(2,j) * s(j) + wk(1,j) * s(j+1) = d(1,j).
!
  if ( n-1 > 1 )  then
    do j = 2, n-1
        if ( wk(2,j-1) == 0.0E+00 )  go to 5008
        g = -wk(1,j+1) / wk(2,j-1)
        d(1,j) = g * d(1,j-1) + 3.0E+00 * ( wk(1,j) * wk(2,j+1) + wk(1,j+1) * wk(2,j) )
        wk(2,j) = g * wk(1,j-1) + 2.0E+00 * ( wk(1,j) + wk(1,j+1) )
    end do
  end if
!
!  Construct last equation from second boundary condition, of the form
!
!    (-g * wk(2,n-1)) * s(n-1) + wk(2,n) * s(n) = d(1,n)
!
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since arrays happen to be set up just right for it
!  at this point.
!
  if ( iend == 1 ) then
    go to 30
  end if

  if ( iend == 0 )  then
     if ( n == 2 .and. ibeg == 0 )  then
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
        d(1,2) = wk(2,2)
        go to 30
     else if ( n == 2 .or. ( n == 3 .and. ibeg == 0 ) )  then
!
!  Either ( N = 3 and not-a-knot also at left) or (N=2 and *not*
!  not-a-knot at left end point).
!
        d(1,n) = 2.0E+00 * wk(2,n)
        wk(2,n) = 1.0E+00
        if ( wk(2,n-1) == 0.0E+00 ) then
          go to 5008
        end if
        g = -1.0E+00 / wk(2,n-1)
     else
!
!  Not-a-knot and N >= 3, and either N > 3 or also not-a-
!  knot at left end point.
!
        g = wk(1,n-1) + wk(1,n)
!
!  Do not need to check following denominators (x-differences).
!
        d(1,n) = ( ( wk(1,n) + 2.0E+00 * g ) * wk(2,n) * wk(1,n-1) &
          + wk(1,n)**2 * ( f(1,n-1) - f(1,n-2) ) / wk(1,n-1) ) / g
        if ( wk(2,n-1) == 0.0E+00 )  go to 5008
        g = -g / wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!
!  Second derivative prescribed at right endpoint.
!
     d(1,n) = 3.0E+00 *wk(2,n) + 0.5E+00 * wk(1,n) * d(1,n)
     wk(2,n) = 2.0E+00
     if ( wk(2,n-1) == 0.0E+00 ) then
       go to 5008
     end if
     g = -1.0E+00 / wk(2,n-1)
  end if
!
!  Complete forward pass of Gauss elimination.
!
  wk(2,n) = g * wk(1,n-1) + wk(2,n)

  if ( wk(2,n) == 0.0E+00 ) then
    go to 5008
  end if

  d(1,n) = ( g * d(1,n-1) + d(1,n) ) / wk(2,n)
!
!  Carry out back substitution.
!
   30 continue

  do j = n-1, 1, -1
     if ( wk(2,j) == 0.0E+00 ) then
       go to 5008
     end if
     d(1,j) = ( d(1,j) - wk(1,j) * d(1,j+1) ) / wk(2,j)
  end do

  return
!
!  error returns.
!
 5004 continue
!
!  ic out of range return.
!
  ierr = ierr - 3
  call xerror ('pchsp -- ic out of range', 24, ierr, 1)
  return

 5007 continue
!
!  nwk too small return.
!
  ierr = -7
  call xerror ('pchsp -- work array too small', 29, ierr, 1)
  return

 5008 continue
!  singular system.
!  theoretically, this can only occur if successive x-values
!  are equal, which should already have been caught (ierr=-3).
  ierr = -8
  call xerror ('pchsp -- singular linear system', 31, ierr, 1)
  return
!
 5009 continue
!  error return from pchdf.
!  this case should never occur.
  ierr = -9
  call xerror ('pchsp -- error return from pchdf', 32, ierr, 1)
  return
end

subroutine xerror ( messg, nmessg, nerr, level )
!
!*******************************************************************************
!
!! XERROR processes an error (diagnostic) message.
!
!
!  Discussion:
!
!    XERROR processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    See subroutine xsetf for details.
!
!  Example:
!
!    call xerror('smooth -- num was zero.',23,1,2)
!
!    call xerror('integ  -- less than full accuracy achieved.',43,2,1)
!
!    call xerror('rooter -- actual zero of f found before interval fully collapsed.',65,3,0)
!
!    call xerror('exp    -- underflows being set to zero.',39,1,-1)
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed, 
!    containing no more than 72 characters.
!
!    Input, integer NMESSG, the actual number of characters in MESSG.
!
!    Input, integer NERR, the error number associated with this message.
!    NERR must not be zero.
!
!    Input, integer LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if XSETF has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most once, 
!      regardless of how many times this call is executed.
!
  implicit none
!
  integer level
  character ( len = * ) messg
  integer nerr
  integer nmessg
!
  call xerrwv ( messg, nmessg, nerr, level, 0, 0, 0, 0, 0.0E+00, 0.0E+00 )

  return
end


subroutine xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )
!
!*******************************************************************************
!
!! XERRWV processes an error message that includes numeric information.
!
!
!  Discussion:
!
!    XERRWV processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    (see subroutine xsetf for details.)
!    in addition, up to two integer values and two real
!    values may be printed along with the message.
!
!  Example:
!
!    call xerrwv ( 'smooth -- num (=i1) was zero.',29,1,2,1,num,0,0,0.,0.)
!
!    call xerrwv ( 'quadxy -- requested error (r1) less than minimum (r2).,54,77,1,0,0,0,2,errreq,errmin)
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer NMESSG, the number of characters in MESSG.
!
!    Input, integer NERR, the error number associated with this message.
!    NERR must not be zero.
!
!    Input, integer LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if xsetf has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most 
!      once, regardless of how many times this call is executed.
!
!    Input, integer NI, the number of integer values to be printed. (0 to 2)
!
!    Input, integer I1, I2, the first and second integer values.
!
!    Input, integer NR, the number of real values to be printed. (0 to 2)
!
!    Input, real R1, R2, the first and second real values.
!
  implicit none
!
  character ( len = 37 ) form
  integer i
  integer i1
  integer i1mach
  integer i2
  integer ifatal
  integer isizei
  integer isizef
  integer iunit
  integer j4save
  integer junk
  integer kdummy
  integer kount
  integer kunit
  integer lerr
  integer level
  character ( len = 20 ) lfirst
  integer lkntrl
  integer llevel
  integer lmessg
  integer lun(5)
  integer maxmes
  character ( len = * ) messg
  integer mkntrl
  integer nerr
  integer ni
  integer nmessg
  integer nr
  integer nunit
  real r1
  real r2
!
!  Get flags
!
  lkntrl = j4save(2,0,.false.)
  maxmes = j4save(4,0,.false.)
!
!  Check for valid input
!
  if ( (nmessg>0) .and. (nerr/=0) .and. (level>=(-1)) .and. (level<=2) ) then
    go to 10
  end if

    if ( lkntrl > 0 ) then
      call xerprt('fatal error in...',17)
    end if

    call xerprt('xerror -- invalid input',23)

    if ( lkntrl > 0 ) then
      call xerprt('job abort due to fatal error.',29)
    end if

    if ( lkntrl > 0 ) then
      call xersav ( ' ', 0, 0, 0, kdummy )
    end if

    call xerabt('xerror -- invalid input',23)
    return

   10 continue
!
!  Record the message.
!
  junk = j4save(1,nerr,.true.)
  call xersav ( messg, nmessg, nerr, level, kount )
!
!  Let user override
!
  lfirst = messg
  lmessg = nmessg
  lerr = nerr
  llevel = level
  call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
!
!  Reset to original values.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level
  lkntrl = max ( -2, min ( 2, lkntrl ) )
  mkntrl = abs ( lkntrl )
!
!  Decide whether to print message
!
  if ( llevel < 2 .and. lkntrl == 0 ) go to 100

  if (((llevel == (-1)) .and. ( kount > min ( 1, maxmes ) ) ) &
    .or.((llevel == 0)    .and. ( kount>maxmes)) &
    .or.((llevel == 1)    .and. ( kount>maxmes).and.(mkntrl==1) ) &
    .or.((llevel == 2)    .and. ( kount > max ( 1, maxmes ) ) ) ) then
    go to 100
  end if

  if ( lkntrl > 0 ) then

    call xerprt(' ',1)

    if ( llevel == -1 ) then

      call xerprt &
      ( 'warning message...this message will only be printed once.',57)

    end if

    if ( llevel == 0 ) then
      call xerprt ( 'warning in...', 13 ) 
    else if ( llevel == 1 ) then
      call xerprt ( 'recoverable error in...', 23 )
    else if ( llevel == 2 ) then
      call xerprt ( 'fatal error in...', 17 )
    end if

  end if
!
!  Message
!
     call xerprt(messg,lmessg)
     call xgetua(lun,nunit)
     isizei = log10 ( real ( i1mach(9) ) ) + 1.0E+00
     isizef = log10 ( real ( i1mach(10) )**i1mach(11) ) + 1.0E+00

     do kunit = 1, nunit

        iunit = lun(kunit)

        do i = 1, min ( ni, 2 )
           write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
           if ( iunit == 0 ) then
             if (i == 1) write (*,form) i1
             if (i == 2) write (*,form) i2
           else
             if (i == 1) write (iunit,form) i1
             if (i == 2) write (iunit,form) i2
           end if
        end do

        do i = 1, min ( nr, 2 )
           write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
           if ( iunit == 0 ) then
             if ( i == 1 ) write (*,form) r1
             if ( i == 2 ) write (*,form) r2
           else
             if (i == 1) write (iunit,form) r1
             if (i == 2) write (iunit,form) r2
           end if
        end do

        if (lkntrl<=0) go to 40
!
!  error number
!
           if ( iunit == 0 ) then
             write(*,30) lerr
           else
             write (iunit,30) lerr
           end if
   30          format (15h error number =,i10)
   40       continue

     end do
!
!  trace-back
!
  100 continue
  ifatal = 0
  if ((llevel == 2).or.((llevel==1) .and. (mkntrl==2))) then
    ifatal = 1
  end if
!
!  quit here if message is not fatal
!
  if ( ifatal <= 0 ) then
    return
  end if

  if ( lkntrl <= 0 .or. kount > max ( 1, maxmes ) ) go to 120
!
!  Print reason for abort
!
     if ( llevel == 1 ) then
       call xerprt ('job abort due to unrecovered error.',35)
     end if

     if ( llevel == 2 ) then
       call xerprt('job abort due to fatal error.',29)
     end if
!
!  Print error summary
!
     call xersav ( ' ', -1, 0, 0, kdummy )

  120 continue
!
!  Abort
!
  if ( llevel == 2 .and. kount > max ( 1, maxmes ) ) then
    lmessg = 0
  end if

  call xerabt(messg,lmessg)

  return
end



subroutine xerprt ( messg, nmessg )
!
!*******************************************************************************
!
!! XERPRT prints a message on each file indicated by xgetua.
!
!
!  Reference:
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be printed.
!
!    Input, integer NMESSG, the actual number of characters in MESSG.
!
  implicit none
!
  integer ichar
  integer iunit
  integer kunit
  integer last
  integer lenmes
  integer lun(5)
  character ( len = * ) messg
  integer nmessg
  integer nunit
!
!  Obtain unit numbers and write line to each unit
!
  call xgetua ( lun, nunit )

  lenmes = len ( messg )

  do kunit = 1, nunit

     iunit = lun(kunit)

     do ichar = 1, lenmes, 72
        last = min ( ichar+71 , lenmes )
        if ( iunit == 0 ) then
          write (*,'(1x,a)') messg(ichar:last)
        else
          write (iunit,'(1x,a)') messg(ichar:last)
        end if
    end do

  end do

  return
end

subroutine xgetua ( iunita, n )
!
!*******************************************************************************
!
!! XGETUA returns the unit number(s) to which error messages are being sent.
!
!
!  Discussion:
!
!    XGETUA may be called to determine the unit number or numbers to which 
!    error messages are being sent.  These unit numbers may have been set 
!    by a call to XSETUN, or a call to XSETUA, or may be a default value.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!    Output, integer IUNITA(N),  an array unit numbers, 
!    A value of zero refers to the default unit, as defined by the 
!    I1MACH machine constant routine.  Only IUNITA(1),..., IUNITA(N) are
!    defined by XGETUA.  The values of IUNITA(N+1),..., IUNITA(5) are 
!    not defined (for N < 5) or altered in any way by XGETUA.
!
!    Output, integer N, the number of units to which copies of the
!    error messages are being sent.  N will be in the range from 1 to 5.
!
  implicit none
!
  integer i
  integer index
  integer iunita(5)
  integer j4save
  integer n
!
  n = j4save ( 5, 0, .false. )

  do i = 1, n

    index = i+4
    if ( i == 1 ) then
      index = 3
    end if

    iunita(i) = j4save(index,0,.false.)

  end do

  return
end

function j4save ( iwhich, ivalue, iset )
!
!*******************************************************************************
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!    Input, integer IWHICH, the index of the item desired.
!    1, the current error number.
!    2, the current error control flag.
!    3, the current unit number to which error messages are sent.
!       (0 means use standard.)
!    4, the maximum times any message is printed (as set by xermax).
!    5, the number of units to which each error message is written.
!    6, the 2nd unit for error messages.
!    7, the 3rd unit for error messages.
!    8, the 4th unit for error messages.
!    9, the 5th unit for error messages.
!
!    Input, integer IVALUE, the value to be set for the IWHICH-th parameter,
!    if ISET is TRUE.
!
!    Input, logical ISET.
!    TRUE: the IWHICH-th parameter will be given the value, IVALUE.
!
!    Output, integer J4SAVE, the old value of the IWHICH-th parameter.
!
  implicit none
!
  integer, save, dimension ( 9 ) :: iparam = (/ 0, 2, 0, 10, 1, 0, 0, 0, 0 /)
  logical iset
  integer ivalue
  integer iwhich
  integer j4save
!
  j4save = iparam(iwhich)

  if ( iset ) then
    iparam(iwhich) = ivalue
  end if

  return
end

subroutine xerabt ( messg, nmessg )
!
!*******************************************************************************
!
!! XERABT aborts program execution and prints an error message.
!
!
!  Discussion:
!
!    XERABT aborts the execution of the program.  The error message causing 
!    the abort is given in the calling sequence, in case one needs it for 
!    printing on a dayfile, for example.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed, 
!    containing no more than 72 characters.
!
!    Input, integer NMESSG, the actual number of characters in MESSG.
!    If NMESSG is 0, no message is being supplied.
!
  implicit none
!
  character ( len = * ) messg
  integer nmessg
!
  if ( nmessg > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XERABT - Termination after fatal error!'
    write ( *, '(a)' ) trim ( messg )
  end if

  stop
end


subroutine xerctl ( messg1, nmessg, nerr, level, kontrl )
!
!*******************************************************************************
!
!! XERCTL allows user control over handling of individual errors.
!
!
!  Discussion:
!
!    Allows user control over handling of individual errors.
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCTL.
!    If the user has provided his own version of XERCTL, he
!    can then override the value of KONTRL used in processing
!    this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Parameters:
!
!    Input, character ( len = 20 ) MESSG1, the first word (only) of the error 
!    message.
!
!    Input, integer NMESSG, same as in the call to xerror or xerrwv.
!
!    Input, integer NERR, same as in the call to xerror or xerrwv.
!
!    Input, integer LEVEL, same as in the call to xerror or xerrwv.
!
!    Input/output, integer KONTRL.  On input, the current value of the control 
!    flag as set by a call to XSETF.  On output, the new value of kontrl.  
!    If KONTRL is not defined, it will remain at its original value.
!    This changed value affects only the current occurrence of the current 
!    message.
!
  implicit none
!
  integer kontrl
  integer level
  character ( len = 20 ) messg1
  integer nerr
  integer nmessg
!
  return
end


subroutine xersav ( messg, nmessg, nerr, level, icount )
!
!*******************************************************************************
!
!! XERSAV records that an error occurred.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    R E Jones, D K Kahaner, 
!    XERROR, The SLATEC Error Handling Package, 
!    SAND82-0800, Sandia Laboratories, 1982.
!
!  Author:
!
!    Ron Jones
!
!  Parameters:
!
!     --input--
!
!       messg, nmessg, nerr, level are as in xerror,
!       except that when nmessg=0 the tables will be
!       dumped and cleared, and when nmessg is less than zero the
!       tables will be dumped and not cleared.
!
!    Output, integer ICOUNT, the number of times this message has
!    been seen, or zero if the table has overflowed and
!    does not contain this message specifically.
!    when nmessg=0, icount will not be altered.
!
  implicit none
!
  integer i
  integer i1mach
  integer icount
  integer ii
  integer iunit
  integer, save, dimension ( 10 ) :: kount = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer, save :: kountx = 0
  integer kunit
  integer level
  integer, save, dimension ( 10 ) :: levtab
  integer lun(5)
  character ( len = 20 ) mes
  character ( len = * ) messg
  character ( len = 20 ), save, dimension ( 10 ) :: mestab
  integer nerr
  integer, save, dimension ( 10 ) :: nertab
  integer nmessg
  integer nunit
!
!  Dump the table
!
  if ( nmessg <= 0 ) then

     if ( kount(1) == 0 ) then
       return
     end if
!
!  Print to each unit
!
     call xgetua ( lun, nunit )

     do kunit = 1, nunit

       iunit = lun(kunit)

       if ( iunit == 0 ) then
         iunit = i1mach(4)
       end if
!
!  Print table header
!
       write ( iunit, 10 )
   10  format ('0          error message summary'/ &
       ' message start             nerr     level     count')
!
!  print body of table
!
        do i = 1, 10
          if ( kount(i) == 0 ) then
            exit
          end if
          write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15     format (1x,a20,3i10)
        end do
!
!  Print number of other errors
!
        if ( kountx /= 0 ) then
          write (iunit,40) kountx
        end if

   40       format (41h0other errors not individually tabulated=,i10)
        write ( iunit, '(a)' ) ' '
     end do
!
!  Clear the error tables
!
    if ( nmessg == 0 ) then
      kount(1:10) = 0
      kountx = 0
    end if

    return

  end if
!
!  process a message...
!  search for this message, or else an empty slot for this messg,
!  or else determine that the error table is full.
!
  mes = messg

  do i = 1, 10

    ii = i

    if ( kount(i) == 0 ) then
      mestab(ii) = mes
      nertab(ii) = nerr
      levtab(ii) = level
      kount(ii)  = 1
      icount = 1
      return
    end if

    if ( mes /= mestab(i) ) go to 90
    if (nerr/=nertab(i)) go to 90
    if (level/=levtab(i)) go to 90
    go to 100

90  continue

  end do
!
!  table is full
!
  kountx = kountx+1
  icount = 1
  return
!
!  message found in table
!
  100    continue

     kount(ii) = kount(ii) + 1
     icount = kount(ii)

  return
end



subroutine pchfd ( n, x, f, d, incfd, skip, ne, xe, fe, de, ierr )
!
!*******************************************************************************
!
!! PCHFD evaluates a piecewise cubic Hermite function.
!
!
!  Discsussion:
!
!    PCHFD evaluates a piecewise cubic Hermite function and its first
!    derivative at an array of points.  PCHFD may be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.  
!
!    PCHFD evaluates the cubic Hermite function and its first derivative
!    at the points XE.
!
!    If only function values are required, use PCHFE instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Programming notes:
!
!    Most of the coding between the call to CHFDV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFDV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of PCHFD that assumes a strictly 
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points <
!    X(N-1), followed by points > X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are 
!    increasing relative to X; that is, for all K >= J,
!      XE(J) >= X(I)
!    implies   
!      XE(K) >= X(I).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(INCFD,N), the function values.  F(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, real D(INCFD,N), the derivative values.  D(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, integer INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with IERR >= 0, SKIP will be set to TRUE.
!
!    Input, integer NE, the number of evaluation points.
!
!    Input, real XE(NE), points at which the function is to be evaluated.
!
!    Output, real FE(NE), the values of the cubic Hermite function at XE.
!
!    Output, real DE(NE), the derivative of the cubic Hermite function at XE.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, if an error has occurred in the lower-level routine CHFDV.  
!
  implicit none
!
  integer incfd
  integer n
  integer ne
!
  real d(incfd,n)
  real de(ne)
  real f(incfd,n)
  real fe(ne)
  integer i
  integer ierc
  integer ierr
  integer ir
  integer j
  integer jfirst
  integer next(2)
  integer nj
  logical skip
  real x(n)
  real xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchfd -- number of data points less than two', 44, ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchfd -- increment less than one', 32, ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchfd -- x-array not strictly increasing', 40, ierr, 1)
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    call xerror ('pchfd -- number of evaluation points less than one', 50, &
      ierr, 1)
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  ( interval index is  il = ir-1  . )
!  ( interval is x(il)<=x<x(ir) . )
!
  jfirst = 1
  ir = 2

10 continue
!
!  Skip out of loop if have processed all evaluation points.
!
     if ( jfirst > ne ) then
       return
     end if
!
!  Locate all points in interval.
!
     do j = jfirst, ne
       if ( xe(j) >= x(ir) )  go to 30
     end do

     j = ne + 1
     go to 40
!
!  Have located first point beyond interval.
!
30    continue

     if ( ir == n ) then
       j = ne + 1
     end if

40    continue

     nj = j - jfirst
!
!  Skip evaluation if no points in interval.
!
     if ( nj == 0 )  go to 50
!
!  Evaluate cubic at xe(i),  i = jfirst (1) j-1 .
!
    call chfdv ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
      nj, xe(jfirst), fe(jfirst), de(jfirst), next, ierc)

     if ( ierc < 0 )  go to 5005

     if ( next(2) == 0 )  go to 42
!
!  In the current set of XE-points, there are next(2) to the right of x(ir).
!
        if ( ir < n )  go to 41
!
!  These are actually extrapolation points.
!
           ierr = ierr + next(2)
           go to 42
41       continue
!
!  We should never have gotten here.
!
           go to 5005

42    continue

     if ( next(1) == 0 )  go to 49
!
!  In the current set of xe-points, there are next(1) to the left of x(ir-1).
!
        if ( ir > 2 )  go to 43
!
!  These are actually extrapolation points.
!
           ierr = ierr + next(1)
           go to 49
43       continue
!
!  XE is not ordered relative to x, so must adjust evaluation interval.
!
!  First, locate first point to left of X(IR-1).
!
           do i = jfirst, j-1
             if ( xe(i) < x(ir-1) )  go to 45
           end do
!
!  Cannot drop through here unless there is an error in chfdv.
!
           go to 5005

   45          continue
!
!  Reset J.  This will be the new JFIRST.
!
           j = i
!
!  Now find out how far to back up in the x-array.
!
           do i = 1, ir-1
             if ( xe(j) < x(i) ) then
               exit
             end if
           end do
!
!  Can never drop through here, since xe(j)<x(ir-1).
!
   47          continue
!
!  At this point, either xe(j) < x(1) or x(i-1) <= xe(j) < x(i) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
           ir = max ( 1, i-1 )
   49    continue

     jfirst = j
!
!  End of ir-loop.
!
   50 continue
  ir = ir + 1
  if ( ir <= n )  go to 10

  return
!
!  error returns.
!


 5005 continue
!
!     error return from chfdv.
!  this case should never occur.
  ierr = -5
  call xerror ('pchfd -- error return from chfdv -- fatal', 41, ierr, 2)
  return
end

subroutine pchim ( n, x, f, d, incfd, ierr )
!
!*******************************************************************************
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.  The interpolant will have 
!    an extremum at each point where monotonicity switches direction.  
!    See PCHIC if user control is desired over boundary or switch conditions.
!
!    If the data are only piecewise monotonic, the interpolant will
!    have an extremum at each point where monotonicity switches direction.  
!    See PCHIC if user control is desired in such cases.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  References:
!
!    Fred Fritsch and R Carlson, 
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch and J Butland,
!    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
!    LLNL Preprint UCRL-87559, April 1982.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(INCFD,N), dependent variable values to be interpolated.  
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).  PCHIM is designed 
!    for monotonic data, but it will work for any F-array.  It will force 
!    extrema at points where monotonicity switches direction.  If some other 
!    treatment of switch points is desired, PCHIC should be used instead.
!
!    Output, real D(INCFD,N), the derivative values at the data points.
!    If the data are monotonic, these values will determine a monotone 
!    cubic Hermite function.  The value corresponding to X(I) is stored in
!    D(1+(I-1)*INCFD).
!
!    Input, integer INCFD, increment between successive values in F and D.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, means that IERR switches in the direction of monotonicity 
!      were detected.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real del1
  real del2
  real dmax
  real dmin
  real drat1
  real drat2
  real dsave
  real f(incfd,n)
  real h1
  real h2
  real hsum
  real hsumt3
  integer i
  integer ierr
  integer nless1
  real pchst
  real w1
  real w2
  real x(n)
!
!  Check arguments.
!
  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchim -- number of data points less than two', 44, ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchim -- increment less than one', 32, ierr, 1)
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      call xerror ('pchim -- x-array not strictly increasing', 40, ierr, 1)
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(1,2) - f(1,1) ) / h1
  dsave = del1
!
!  Special case N=2 -- use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, N >= 3.
!
  h2 = x(3) - x(2)
  del2 = ( f(1,3) - f(1,2) ) / h2
!
!  Set D(1) via non-centered three-point formula, adjusted to be
!  shape-preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0E+00 )  then
     d(1,1) = 0.0E+00
  else if ( pchst ( del1, del2 ) < 0.0E+00 )  then
!
!  Need do this check only if monotonicity switches.
!
     dmax = 3.0E+00 * del1

     if ( abs ( d(1,1) ) > abs ( dmax ) ) then
       d(1,1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( i > 2 ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(1,i+1) - f(1,i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0E+00
    if ( pchst ( del1, del2 ) )  42, 41, 45
!
!  Count number of changes in direction of monotonicity.
!
   41    continue
    if ( del2 /= 0.0E+00 ) then
      if ( pchst ( dsave, del2 ) < 0.0E+00 ) then
        ierr = ierr + 1
      end if
      dsave = del2
    end if
    go to 50

42  continue

    ierr = ierr + 1
    dsave = del2
    go to 50
!
!  Use Brodlie modification of Butland formula.
!
45  continue

    hsumt3 = 3.0E+00 * hsum
    w1 = ( hsum + h1 ) / hsumt3
    w2 = ( hsum + h2 ) / hsumt3
    dmax = max ( abs ( del1 ), abs ( del2 ) )
    dmin = min ( abs ( del1 ), abs ( del2 ) )
    drat1 = del1 / dmax
    drat2 = del2 / dmax
    d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

   50 continue

  end do
!
!  Set D(N) via non-centered three-point formula, adjusted to be
!  shape-preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0E+00 )  then
    d(1,n) = 0.0E+00
  else if ( pchst(del1,del2) < 0.0E+00 )  then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0E+00 * del2

    if ( abs ( d(1,n) ) > abs ( dmax ) ) then
      d(1,n) = dmax
    end if

  end if

  return
end


function pchdf ( k, x, s, ierr )
!
!*******************************************************************************
!
!! PCHDF approximates a derivative using divided differences of data.
!
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3 and 4 point boundary
!    derivative approximations.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978), pp. 10-16.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer K, is the order of the desired derivative approximation.
!    K must be at least 3.
!
!    Input, real X(K), contains the K values of the independent variable.
!    X need not be ordered, but the values must be distinct.
!
!    Input/output, real S(K-1).  On input, the associated slope values:
!      S(I) = ( F(I+1)-F(I))/(X(I+1)-X(I))
!    On output, S is overwritten.
!
!    Output, integer IERR, error flag.
!    0, no error.
!    -1, if K < 2.
!
!    Output, real PCHDF, the desired derivative approximation if
!    IERR=0 or to zero if IERR=-1.
!
  implicit none
!
  integer k
!
  integer i
  integer ierr
  integer j
  real pchdf
  real s(k-1)
  real value
  real x(k)
!
!  Check for legal value of K.
!
  if ( k < 3 ) then
    ierr = -1
    call xerror ( 'pchdf -- k less than three', 26, ierr, 1 )
    pchdf = 0.0E+00
    return
  end if
!
!  Compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = ( s(i+1) - s(i) ) / ( x(i+j) - x(i) )
    end do
  end do
!
!  Evaluate the derivative at X(K).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do

  ierr = 0
  pchdf = value

  return
end


function pchst ( arg1, arg2 )
!
!*******************************************************************************
!
!! PCHST: PCHIP sign-testing routine.
!
!
!  Discussion:
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, real ARG1, ARG2, two values to check.
!
!    Output, real PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  implicit none
!
  real arg1
  real arg2
  real pchst
!
  pchst = sign ( 1.0E+00, arg1 ) * sign ( 1.0E+00, arg2 )

  if ( arg1 == 0.0E+00 .or. arg2 == 0.0E+00 ) then
    pchst = 0.0E+00
  end if

  return
end

subroutine chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )
!
!*******************************************************************************
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial given in Hermite form and its
!    first derivative at an array of points.  While designed for
!    use by PCHFD, it may be useful directly as an evaluator for
!    a piecewise cubic Hermite function in applications, such as
!    graphing, where the interval is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Parameters:
!
!    Input, real X1, X2, the endpoints of the interval of definition of 
!    the cubic.  X1 and X2 must be distinct.
!
!    Input, real F1, F2, the values of the function at X1 and X2, respectively.
!
!    Input, real D1, D2, the derivative values at the ends of the interval.
!
!    Input, integer NE, the number of evaluation points.
!
!    Input, real XE(NE), the points at which the functions are to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in next.
!
!    Output, real FE(NE), DE(NE), the values of the cubic function and 
!    its derivative at the points XE(*).
!
!    Output, integer NEXT(2), indicates the number of extrapolation points:
!    NEXT(1) = number of evaluation points to left of interval.
!    NEXT(2) = number of evaluation points to right of interval.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none
!
  integer ne
!
  real c2
  real c2t2
  real c3
  real c3t3
  real d1
  real d2
  real de(ne)
  real del1
  real del2
  real delta
  real f1
  real f2
  real fe(ne)
  real h
  integer i
  integer ierr
  integer next(2)
  real x
  real x1
  real x2
  real xe(ne)
  real xma
  real xmi
!
!  Check arguments.
!
  if ( ne < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1

  if ( h == 0.0E+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The interval endpoints are equal.'
    return
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0E+00, h )
  xma = max ( 0.0E+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h

  c2 = -( del1 + del1 + del2 )
  c2t2 = c2 + c2
  c3 = ( del1 + del2 ) / h
  c3t3 = c3 + c3 + c3
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
    de(i) = d1 + x * ( c2t2 + x * c3t3 )
!
!  Count extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( x > xma ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end


function pchqa ( n, x, f, d, a, b, ierr )
!
!*******************************************************************************
!
!! PCHQA: easy to use cubic Hermite or spline integration.
!
!
!  Discussion:
!
!    PCHQA evaluates the definite integral of a cubic Hermite or spline
!    function over the interval [A, B].  This is an easy to use driver
!    for the routine PCHIA.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    Fred Fritsch and R Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real D(N), the derivative values.  D(I) is the value
!    corresponding to X(I).
!
!    Input, real A, B, the limits of integration.  The interval [A,B] is
!    normally contained in [X(1),X(N)], but this is not required.
!
!    Output, integer IERR, error flag.
!    0, no errors).
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2 .
!    -3, if the X array is not strictly increasing.
!
!    Output, real PCHQA, the value of the requested integral.
!
  implicit none
!
  integer n
!
  real a
  real b
  real d(n)
  real f(n)
  integer ierr
  integer, save :: incfd = 1
  real pchia
  real pchqa
  logical, save :: skip = .true.
  real x(n)
!
  pchqa  =  pchia ( n, x, f, d, incfd, skip, a, b, ierr )

  return
end

function pchia ( n, x, f, d, incfd, skip, a, b, ierr )
!
!*******************************************************************************
!
!! PCHIA evaluates the integral of a piecewise cubic Hermite function.
!
!
!  Description:
!
!    PCHIA evaluates the definite integral of a cubic Hermite function
!    over the interval [A, B].
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Output, real VALUE, the value of the requested integral.
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(INCFD,N), the function values.  F(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, real D(INCFD,N), the derivative values.  D(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, integer INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with IERR >= 0, SKIP will be set to TRUE.
!
!    Input, real A, B, the limits of integration.  The integration interval
!    is normally contained within [X(1),X(N)], but this is not required.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none
!
  integer incfd
  integer n
!
  real a
  real b
  real chfiv
  real d(incfd,n)
  real f(incfd,n)
  integer i
  integer ia
  integer ib
  integer ierd
  integer ierr
  integer ierv
  integer il
  integer ir
  real pchia
  real pchid
  logical skip
  real value
  real x(n)
  real xa
  real xb
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchia -- number of data points less than two', 44, ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchia -- increment less than one', 32, ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchia -- x-array not strictly increasing', 40, ierr, 1)
        return
      end if
    end do

    skip = .true.

  end if

  ierr = 0

  if ( a < x(1) .or. a > x(n) ) then
    ierr = ierr + 1
  end if

  if ( b < x(1) .or. b > x(n) ) then
    ierr = ierr + 2
  end if
!
!  Compute integral value.
!
  if ( a == b )  then
    value = 0.0E+00
  else
     xa = min (a, b)
     xb = max (a, b)
     if ( xb <= x(2) )  then
!
!  Interval is to left of X(2), so use first cubic.
!
        value = chfiv ( x(1), x(2), f(1,1),f(1,2), &
          d(1,1),d(1,2), a, b, ierv)

        if ( ierv < 0 ) then
          ierr = -4
          call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
          return
        end if

     else if ( xa >= x(n-1) )  then
!
!  Interval is to right of x(n-1), so use last cubic.
!
        value = chfiv ( x(n-1), x(n), f(1,n-1), f(1,n), &
          d(1,n-1), d(1,n), a, b, ierv )

        if ( ierv < 0 ) then
          ierr = -4
          call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
          return
        end if

     else
!
!  Normal case -- xa<xb, xa<x(n-1), xb>x(2).
!
!  Locate ia and ib such that
!  x(ia-1) < xa <= x(ia) <= x(ib) <= xb <= x(ib+1)
!
        ia = 1
        do i = 1, n-1
           if ( xa > x(i) )  ia = i + 1
        end do
!
!  IA = 1 implies xa<x(1) .  Otherwise,
!  ia is largest index such that x(ia-1)<xa,.
!
        ib = n
        do i = n, ia, -1
          if ( xb < x(i) ) then
            ib = i - 1
          end if
        end do
!
!  IB = N implies XB > X(N).  Otherwise,
!  ib is smallest index such that xb<x(ib+1) .
!
!  Compute the integral.
!
        ierv = 0
        if ( ib < ia )  then
!
!  This means IB = IA-1 and (A,B) is a subset of (x(ib),x(ia)).
!
           value = chfiv (x(ib), x(ia), f(1,ib), f(1,ia), &
             d(1,ib), d(1,ia), a, b, ierv )

           if ( ierv < 0 ) then
             ierr = -4
             call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
             return
           end if

        else
!
!  First compute integral over (x(ia),x(ib)).
!
           if ( ib == ia )  then
              value = 0.0E+00
           else

              value = pchid ( n, x, f, d, incfd, skip, ia, ib, ierd )

              if ( ierd < 0 )  go to 5005
           end if
!
!  Then add on integral over ( XA, X(IA) ).
!
           if ( xa < x(ia) )  then
              il = max ( 1, ia-1 )
              ir = il + 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), xa, x(ia), ierv )

              if ( ierv < 0 ) then
                ierr = -4
                call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
                return
              end if

           end if
!
!  Then add on integral over ( X(IB), XB ).
!
           if ( xb > x(ib) )  then
              ir = min ( ib+1, n )
              il = ir - 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), x(ib), xb, ierv )

              if ( ierv < 0 ) then
                ierr = -4
                call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
                return
              end if

           end if
!
!  Adjust sign if necessary.
!
           if ( a > b ) then
             value = -value
           end if

        end if
     end if
  end if

  pchia = value
  return
!
!  error returns.
!

 5005 continue
!     trouble in pchid.  (should never occur.)
  ierr = -5
  call xerror ('pchia -- trouble in pchid', 25, ierr, 1)
  return
end

function chfiv ( x1, x2, f1, f2, d1, d2, a, b, ierr )
!
!*******************************************************************************
!
!! CHFIV evaluates the integral of a cubic polynomial in Hermite form.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Discussion:
!
!    CHFIV is called by PCHIA to evaluate the integral of a single cubic (in
!    Hermite form) over an arbitrary interval (A,B).
!
!  Parameters:
!
!    Output, real VALUE, the value of the requested integral.
!
!    Input, real X1, X2, the endpoints of the interval of definition of 
!    the cubic.  X1 and X2 must be distinct.
!
!    Input, real F1, F2, the values of the function at X1 and X2, respectively.
!
!    Input, real D1, D2, the derivative values at the ends of the interval.
!
!    Input, real A, B, the endpoints of interval of integration.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, X1 == X2.
!
  implicit none
!
  real a
  real b
  real chfiv
  real d1
  real d2
  real dterm
  real f1
  real f2
  real fterm
  real h
  integer ierr
  real phia1
  real phia2
  real phib1
  real phib2
  real psia1
  real psia2
  real psib1
  real psib2
  real ta1
  real ta2
  real tb1
  real tb2
  real ua1
  real ua2
  real ub1
  real ub2
  real x1
  real x2
!
!  Check input.
!
  if ( x1 == x2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  X1 = X2.'
    stop
  end if

  ierr = 0
!
!  Compute integral.
!
  h = x2 - x1
  ta1 = ( a - x1 ) / h
  ta2 = ( x2 - a ) / h
  tb1 = ( b - x1 ) / h
  tb2 = ( x2 - b ) / h

  ua1 = ta1**3
  phia1 = ua1 * ( 2.0E+00 - ta1 )
  psia1 = ua1 * ( 3.0E+00 * ta1 - 4.0E+00 )
  ua2 = ta2**3
  phia2 =  ua2 * ( 2.0E+00 - ta2)
  psia2 = -ua2 * ( 3.0E+00 * ta2 - 4.0E+00 )

  ub1 = tb1**3
  phib1 = ub1 * ( 2.0E+00 - tb1 )
  psib1 = ub1 * ( 3.0E+00 * tb1 - 4.0E+00 )
  ub2 = tb2**3
  phib2 =  ub2 * ( 2.0E+00 - tb2 )
  psib2 = -ub2 * ( 3.0E+00 * tb2 - 4.0E+00 )

  fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 )
  dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) * ( h / 6.0E+00 )

  chfiv = 0.5E+00 * h * ( fterm + dterm )

  return
end

function pchid ( n, x, f, d, incfd, skip, ia, ib, ierr )
!
!*******************************************************************************
!
!! PCHID evaluates the definite integral of a piecewise cubic Hermite function.
!
!
!  Description:
!
!    PCHID evaluates the definite integral of a cubic Hermite function
!    over the interval [X(IA), X(IB)].  The endpoints of the integration 
!    interval must be data points.
!
!    To provide compatibility with PCHIM and pchic, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:  
!
!    Fred Fritsch,  
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent variable values.
!
!    Input, real F(INCFD,N), the function values.  F(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, real F(INCFD,N), the derivative values.  D(1+(I-1)*INCFD) is
!    the value corresponding to X(I).
!
!    Input, integer INCFD, increment between successive values in F and D.
!
!    Input/output, logical SKIP, should be set to TRUE if the user wishes to 
!    skip checks for validity of preceding parameters, or to FALSE otherwise.
!    This will save time in case these checks have already been performed 
!    say, in PCHIM or PCHIC.  SKIP will be set to TRUE on return with 
!    IERR = 0 or -4.
!
!    Input, integer IA, IB, the indices in the X array for the limits of 
!    integration.  Both must be in the range [1,N].
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IA or IB is out of range.
!
!    Output, real PCHID, the value of the requested integral.
!
  implicit none
!
  integer incfd
  integer  n
!
  real d(incfd,n)
  real f(incfd,n)
  real h
  integer i
  integer ia
  integer ib
  integer ierr
  integer iup
  integer low
  real pchid
  logical skip
  real sum2
  real value
  real x(n)
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      call xerror ('pchid -- number of data points less than two', 44, ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      call xerror ('pchid -- increment less than one', 32, ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        call xerror ('pchid -- x-array not strictly increasing', 40, ierr, 1)
        return
      end if
    end do

  end if

  skip = .true.

  if ( ia < 1 .or. ia > n ) then
    go to 5004
  end if

  if ( ib < 1 .or. ib > n ) then
    go to 5004
  end if

  ierr = 0
!
!  Compute integral value.
!
  if ( ia == ib )  then

    value = 0.0E+00

  else

    low = min ( ia, ib )
    iup = max ( ia, ib ) - 1
    sum2 = 0.0E+00

    do i = low, iup
      h = x(i+1) - x(i)
      sum2 = sum2 + h * ( ( f(1,i) + f(1,i+1) ) + ( d(1,i) - d(1,i+1) ) * ( h / 6.0E+00 ) )
    end do

    value = 0.5E+00 * sum2

    if ( ia > ib ) then
      value = -value
    end if

  end if

  pchid = value

  return
!
!  error returns.
!
 5004 continue
!     ia or ib out of range return.
  ierr = -4
  call xerror ('pchid -- ia or ib out of range', 30, ierr, 1)
  return
end

function i1mach ( i )
!
!*******************************************************************************
!
!! I1MACH returns integer machine constants.
!
!
!  I/O unit numbers.
!
!    I1MACH(1) = the standard input unit.
!    I1MACH(2) = the standard output unit.
!    I1MACH(3) = the standard punch unit.
!    I1MACH(4) = the standard error message unit.
!
!  Words.
!
!    I1MACH(5) = the number of bits per integer storage unit.
!    I1MACH(6) = the number of characters per integer storage unit.
!
!  Integers.
!
!  Assume integers are represented in the S digit base A form:
!
!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!  where 0<=X(I)<A for I=0 to S-1.
!
!    I1MACH(7) = A, the base.
!    I1MACH(8) = S, the number of base A digits.
!    I1MACH(9) = A**S-1, the largest integer.
!
!  Floating point numbers
!
!  Assume floating point numbers are represented in the T digit base B form:
!
!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!
!  where 0<=X(I)<B for I = 1 to T, 0<X(1) and EMIN<=E<=EMAX
!
!    I1MACH(10) = B, the base.
!
!  Single precision
!
!    I1MACH(11) = T, the number of base B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!  Double precision
!
!    I1MACH(14) = T, the number of base B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment, the desired set of DATA
!  statements should be activated by removing the C from column 1.  On rare
!  machines, a STATIC statement may need to be added, but probably more systems
!  prohibit than require it.
!
!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!  consistency with the local operating system.  For FORTRAN 77, you may wish
!  to adjust the data statement so imach(6) is set to 1, and then to comment
!  out the executable test on I.EQ.6 below.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!
  implicit none
!
  integer i
  integer i1mach
  integer imach(16)
  integer output
!
  equivalence (imach(4),output)
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!  8087 based micros such asthe IBM PC and ATT 6300.
!
   data imach( 1) /    5 /
   data imach( 2) /    6 /
   data imach( 3) /    7 /
   data imach( 4) /    6 /
   data imach( 5) /   32 /
   data imach( 6) /    4 /
   data imach( 7) /    2 /
   data imach( 8) /   31 /
   data imach( 9) / 2147483647 /
   data imach(10) /    2 /
   data imach(11) /   24 /
   data imach(12) / -125 /
   data imach(13) /  128 /
   data imach(14) /   53 /
   data imach(15) / -1021 /
   data imach(16) /  1024 /
!
!  ALLIANT FX/8 UNIX FORTRAN compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  AMDAHL machines.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  BURROUGHS 1700 system.
!
!      data imach( 1) /    7 /
!      data imach( 2) /    2 /
!      data imach( 3) /    2 /
!      data imach( 4) /    2 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   33 /
!      data imach( 9) / Z1FFFFFFFF /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -256 /
!      data imach(13) /  255 /
!      data imach(14) /   60 /
!      data imach(15) / -256 /
!      data imach(16) /  255 /
!
!  BURROUGHS 5700 system.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -50 /
!      data imach(16) /  76 /
!
!  BURROUGHS 6700/7700 systems.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -32754 /
!      data imach(16) /  32780 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   48 /
!      data imach(12) / -974 /
!      data imach(13) / 1070 /
!      data imach(14) /   96 /
!      data imach(15) / -927 /
!      data imach(16) / 1070 /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     7 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) / 9223372036854775807 /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -4095 /
!      data imach(13) /  4094 /
!      data imach(14) /    94 /
!      data imach(15) / -4095 /
!      data imach(16) /  4094 /
!
!  CDC CYBER 200 series
!
!      data imach( 1) /      5 /
!      data imach( 2) /      6 /
!      data imach( 3) /      7 /
!      data imach( 4) /      6 /
!      data imach( 5) /     64 /
!      data imach( 6) /      8 /
!      data imach( 7) /      2 /
!      data imach( 8) /     47 /
!      data imach( 9) / X'00007FFFFFFFFFFF' /
!      data imach(10) /      2 /
!      data imach(11) /     47 /
!      data imach(12) / -28625 /
!      data imach(13) /  28718 /
!      data imach(14) /     94 /
!      data imach(15) / -28625 /
!      data imach(16) /  28718 /
!
!  CDC 6000/7000 series using FTN4.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / 00007777777777777777B /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CDC 6000/7000 series using FTN5.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CONVEX C-1.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1023 /
!      data imach(13) /  1023 /
!      data imach(14) /    53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1021 /
!      data imach(13) /  1024 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /   102 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) /  777777777777777777777B /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -8189 /
!      data imach(13) /  8190 /
!      data imach(14) /    94 /
!      data imach(15) / -8099 /
!      data imach(16) /  8190 /
!
!  DATA GENERAL ECLIPSE S/200.
!
!      data imach( 1) /   11 /
!      data imach( 2) /   12 /
!      data imach( 3) /    8 /
!      data imach( 4) /   10 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) /32767 /
!      data imach(10) /   16 /
!      data imach(11) /    6 /
!      data imach(12) /  -64 /
!      data imach(13) /   63 /
!      data imach(14) /   14 /
!      data imach(15) /  -64 /
!      data imach(16) /   63 /
!
!  ELXSI 6400
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  HARRIS 220
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /   43 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   63 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  HP 2100, 3 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   39 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 2100, 4 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   55 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 9000
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     7 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1015 /
!      data imach(16) /  1017 /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z7FFFFFFF /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data imach( 1) /     4 /
!      data imach( 2) /     7 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!  constants with Y's.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   6 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z'7FFFFFFF' /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  62 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  62 /
!
!  PDP-10 (KA processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   54 /
!      data imach(15) / -101 /
!      data imach(16) /  127 /
!
!  PDP-10 (KI processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   62 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 32-bit integer arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 16-bit integer arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!  supplied by Igor Bray.
!
!      data imach( 1) /            1 /
!      data imach( 2) /            1 /
!      data imach( 3) /            2 /
!      data imach( 4) /            1 /
!      data imach( 5) /           32 /
!      data imach( 6) /            4 /
!      data imach( 7) /            2 /
!      data imach( 8) /           31 /
!      data imach( 9) / :17777777777 /
!      data imach(10) /            2 /
!      data imach(11) /           23 /
!      data imach(12) /         -127 /
!      data imach(13) /         +127 /
!      data imach(14) /           47 /
!      data imach(15) /       -32895 /
!      data imach(16) /       +32637 /
!
!  SEQUENT BALANCE 8000.
!
!      data imach( 1) /     0 /
!      data imach( 2) /     0 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     1 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) /  2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -125 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  SUN 3 (68881 or FPA)
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    6 /
!      data imach( 4) /    0 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  UNIVAC 1100 series.
!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!  instead.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    6 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   60 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  VAX.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  Z80 microprocessor.
!
!      data imach( 1) /    1 /
!      data imach( 2) /    1 /
!      data imach( 3) /    0 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
  if ( i < 1 .or. i > 16 )then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a,i6)' ) '  I is out of bounds:', i
    i1mach = 0
    stop
  else
    i1mach = imach(i)
  end if

  return
end


