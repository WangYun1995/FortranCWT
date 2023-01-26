Module Tools
  Use Parameters
  
  Implicit None
  Contains
  !--------------------------------------------------------------!
  ! The Gaussian-derived wavelet
  !--------------------------------------------------------------!
  Function GDW( var )
    ! Calling variables
    Real(8), Intent(in) :: var
    Real(8) :: GDW

    ! Local variables
    Real(8) :: x, cN, p1, p2

    x   = var
    cN  = 1.0d0/Sqrt( Sqrt(18.0d0*Pi) )
    p1  = Exp(-0.25d0*x**2)
    p2  = (2.0d0-x**2)
    GDW = cN*p1*p2
  End Function GDW

  !--------------------------------------------------------------!
  ! The cosine weighted Gaussian-derived wavelet
  !--------------------------------------------------------------!
  Function CWGDW( var )
    ! Calling variables
    Real(8), Intent(in) :: var
    Real(8) :: CWGDW

    ! Local variables
    Real(8) :: x, cN, p1, p2

    x   = var
    cN  = Sqrt( 8.0d0/((1.0d0+5.0d0*Exp(1.0d0))*Sqrt(Pi)) )
    p1  = Exp(0.5d0*(1.0d0-x**2))
    p2  = (1.0d0-x**2)*Cos(x)-x*Sin(x)
    CWGDW = cN*p1*p2
  End Function CWGDW

  !--------------------------------------------------------------!
  ! The real part of the Morlet wavelet 
  !--------------------------------------------------------------!
  Function rMW( var )
    ! Calling variables
    Real(8), Intent(in) :: var
    Real(8) :: rMW

    ! Local variables
    Real(8) :: x, cN, p1, p2

    x   = var
    cN  = Exp(8.0d0)/Sqrt( (1.0d0-2.0d0*Exp(4.0d0)+Exp(16.0d0)*Sqrt(Pi)) )
    p1  = Exp(-0.5d0*x**2)
    p2  = Cos(4.0d0*x)-Exp(-8.0d0)
    rMW = cN*p1*p2
  End Function rMW

  !--------------------------------------------------------------!
  ! The imaginary part of the Morlet wavelet 
  !--------------------------------------------------------------!
  Function iMW( var )
    ! Calling variables
    Real(8), Intent(in) :: var
    Real(8) :: iMW

    ! Local variables
    Real(8) :: x, cN, p1, p2

    x   = var
    cN  = Exp(8.0d0)/Sqrt( (1.0d0-2.0d0*Exp(4.0d0)+Exp(16.0d0)*Sqrt(Pi)) )
    p1  = Exp(-0.5d0*x**2)
    p2  = -Sin(4.0d0*x)
    iMW = cN*p1*p2
  End Function iMW

  !--------------------------------------------------------------!
  ! Cumulative summation of the input array
  !--------------------------------------------------------------!
  Function CumSum( arr_in )
    ! Calling variables
    Real(8), Dimension(:), Intent(in) :: arr_in
    Real(8), Dimension(Size(arr_in)) :: CumSum
    ! Local variables
    Integer :: i, N
    
    N         = Size(arr_in)
    CumSum(1) = arr_in(1)
    Do i = 2, N
      CumSum(i) = CumSum(i-1)+arr_in(i) 
    End Do
  End Function CumSum

  !--------------------------------------------------------------!
  ! Cumulative summation of the array Four times
  !--------------------------------------------------------------!
  Function CumSum4( arr_in )
    ! Calling variables
    Real(8), Dimension(:), Intent(in) :: arr_in
    Real(8), Dimension(Size(arr_in)) :: CumSum4
    ! Local variables
    Integer :: i, N

    N          = Size(arr_in)
    CumSum4(:) = arr_in(:)
    Do i = 1, 4
      CumSum4(:) = CumSum( CumSum4(:) )
    End Do
  End Function CumSum4

  !--------------------------------------------------------------!
  ! Factorial function
  !--------------------------------------------------------------!
  Integer recursive function fact(n)
    Integer, Intent(in) :: n
    Integer :: i

    fact = 1
    Do i = 1, n
      fact = fact * i
    End do
    return
  End function fact

  !--------------------------------------------------------------!
  ! The 7th order B-spline function
  !--------------------------------------------------------------!
  Function BS7th( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, x2, x3, x4, x5, x6, x7, BS7th

    x  = Abs(var)
    x2 = x**2
    x3 = x*x2
    x4 = x*x3
    x5 = x*x4
    x6 = x*x5
    x7 = x*x6
    If ( (x-1.0d0)<=Error ) Then
      BS7th = 151.0d0/315.0d0 - x2/3.0d0+x4/9.0d0-x6/36.0d0+x7/144.0d0
    Else if ( ((x-1.0d0)>Error).And.((x-2.0d0)<=Error) ) Then
      BS7th = (2472.0d0-392.0d0*x-504.0d0*x2-1960.0d0*x3+2520.0d0*x4                         &
              -1176.0d0*x5+252.0d0*x6-21.0d0*x7)/5040.0d0
    Else if ( ((x-2.0d0)>Error).And.((x-3.0d0)<=Error) ) Then
      BS7th = (-1112.0d0+12152.0d0*x-19320.0d0*x2+13720.0d0*x3-5320.0d0*x4                   &
              +1176.0d0*x5-140.0d0*x6+7.0d0*x7)/5040.0d0
    Else if ( ((x-3.0d0)>Error).And.((x-4.0d0)<=Error) ) Then
      BS7th = -(x-4.0d0)**7/5040.0d0
    Else
      BS7th = 0.0d0
    End If
  End Function BS7th

  !--------------------------------------------------------------!
  ! Coefficients of the cubic B-spline approximation of the wavelet
  !--------------------------------------------------------------!
  Function dcoeffs_BS3rd( wavelet_samples, Nd )
    ! Calling variables
    Integer, Intent(in) :: Nd
    Real(8), Dimension(-Nd-2:Nd+2), Intent(in) :: wavelet_samples
    
    ! Local variables
    Integer :: k
    Real(8), Allocatable, Dimension(:) :: arr, c_plus, c_minus
    Real(8), Dimension(-Nd:Nd) :: dcoeffs_BS3rd

    ! Intialize
    Allocate( arr(-(Nd+Nz1+2):Nd+Nz1+2),                                     &
              c_plus(-(Nd+Nz1+2):Nd+Nz1+2),                                  &
              c_minus(-(Nd+Nz1+2):Nd+Nz1+2) )
    arr(:-Nd-3)             = 0.0d0
    arr(-Nd-2:Nd+2) = wavelet_samples
    arr(Nd+3:)              = 0.0d0
  
    ! Compute c_plus
    c_plus(-Nd-Nz1-2) = 0.0d0
    Do k = -Nd-Nz1-1, Nd+Nz1+2
      c_plus(k) = arr(k) + z1*c_plus(k-1)
    End Do
    Deallocate( arr )

    ! Compute c_minus
    c_minus(Nd+Nz1+2) = 0.0d0
    Do k = Nd+Nz1+1, -Nd-Nz1-2, -1
      c_minus(k) = z1*( c_minus(k+1) - c_plus(k) )
    End Do
    ! Deallocate c_plus
    Deallocate( c_plus )

    ! The dcoeffs
    dcoeffs_BS3rd(:) = 6.0d0*c_minus(-Nd:Nd)
  End Function dcoeffs_BS3rd

  !--------------------------------------------------------------!
  ! p coefficients
  !--------------------------------------------------------------!
  Function pcoeffs( wavelet_name, chi, Nd )
    ! Calling variables
    Integer, Intent(in) :: Nd
    Real(8), Intent(in) :: chi
    Character(len=*), Intent(in) :: wavelet_name
    
    ! Local variables
    Integer :: k, l
    Real(8) :: cN, hh
    Real(8), Dimension(-Nd-2:Nd+2) :: wavelet_samples
    Real(8), Dimension(-Nd:Nd) :: dk
    Real(8), Dimension(0:4) :: ql
    Real(8), Dimension(-Nd-4:Nd+4) :: dk_tem
    Real(8), Dimension(-Nd:Nd+4) :: pcoeffs
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set, realMW_set, imagMW_set

    ! Check if there is a wavelet in this set
    CBSW_set   = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set    = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set  = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    realMW_set = "rmw, rMW, real part of Morlet wavelet"
    imagMW_set = "imw, iMW, imaginary part of Morlet wavelet"
    If ( index(CBSW_set, wavelet_name)/= 0 ) Then
      cN    = Sqrt(30.0d0/31.0d0)
      dk(:) =  (/-1.0d0*cN, 2.0d0*cN,-1.0d0*cN /)
    Else if ( index(GDW_set, wavelet_name)/= 0 ) Then
      hh = chi/(Nd+2)
      Do k = -Nd-2, Nd+2
        wavelet_samples(k) = GDW( k*hh )
      End Do
      dk = dcoeffs_BS3rd(wavelet_samples, Nd)
    Else if ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      hh = chi/(Nd+2)
      Do k = -Nd-2, Nd+2
        wavelet_samples(k) = CWGDW( k*hh )
      End Do
      dk = dcoeffs_BS3rd(wavelet_samples, Nd)
    Else if ( index(realMW_set, wavelet_name)/= 0 ) Then
      hh = chi/(Nd+2)
      Do k = -Nd-2, Nd+2
        wavelet_samples(k) = rMW( k*hh )
      End Do
      dk = dcoeffs_BS3rd(wavelet_samples, Nd)
    Else if ( index(imagMW_set, wavelet_name)/= 0 ) Then
      hh = chi/(Nd+2)
      Do k = -Nd-2, Nd+2
        wavelet_samples(k) = iMW( k*hh )
      End Do
      dk = dcoeffs_BS3rd(wavelet_samples, Nd)
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If

    ! Compute ql coefficients
    Do l = 0, 4
      ql(l) = ((-1.0d0)**l)*Real( fact(4)/( fact(l)*fact(4-l) ), 8)
    End Do

    ! Compute pcoeffs
    pcoeffs(:)                    = 0.0d0
    dk_tem(:-Nd-1)         = 0.0d0
    dk_tem(-Nd:Nd) = dk(:)
    dk_tem(Nd+1:)          = 0.0d0
    Do k = -Nd, Nd+4
      Do l = 0, 4
        pcoeffs(k) = pcoeffs(k) + dk_tem(k-l)*ql(l)
      End Do
    End Do
  End Function pcoeffs

  !--------------------------------------------------------------!
  ! Coefficients of the cubic B-spline approximation of the signal,
  ! which satisfies the zero boundary condition
  !--------------------------------------------------------------!
  Subroutine InterpCoeffs_zerobc( arr_in, Nmesh, coeffs )
    Integer, Intent(in) :: Nmesh
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: arr_in
    Real(8), Dimension(-6:Nmesh+1), Intent(out) :: coeffs

    ! Local variables
    Integer :: i, k
    Real(8) :: z1_ps
    Real(8), Allocatable, Dimension(:) :: c_plus
    Real(8), Allocatable, Dimension(:) :: c_minus
    
    ! Allocate arrays
    Allocate( c_plus(0:Nmesh-1), c_minus(-6:Nmesh+1) )

    ! Compute c_plus
    c_plus(0) = arr_in(0)
    Do k = 1, Nmesh-1
      c_plus(k) = arr_in(k) + z1*c_plus(k-1)
    End Do

    ! Compute c_minus(Nmesh-1)
    z1_ps = 0.0d0
    Do i = 1, Nz1
      z1_ps = z1_ps + z1_2**i
    End Do
    c_minus(Nmesh-1) = - c_plus(Nmesh-1)*z1_ps/z1
    ! Compute c_minus(0:Nmesh-2)
    Do k = Nmesh-2, 0, -1
      c_minus(k) = z1*( c_minus(k+1) - c_plus(k) )
    End Do
    ! Compute c_minus(-6:-1)
    Do k = -6, -1
      c_minus(k) = c_minus(0)*(z1**(-k))
    End Do
    ! Compute c_minus(Nmesh:Nmesh+1)
    c_minus(Nmesh)   = - c_plus(Nmesh-1)*z1_ps
    c_minus(Nmesh+1) = - c_plus(Nmesh-1)*z1_ps*z1
    ! Deallocate c_plus
    Deallocate( c_plus )
    ! The final c_minus
    coeffs(:)  = 6.0d0*c_minus(:)
    ! Deallocate c_minus
    Deallocate( c_minus )
  End Subroutine InterpCoeffs_zerobc

  !--------------------------------------------------------------!
  ! Coefficients of the cubic B-spline approximation of the signal,
  ! which satisfies the periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine InterpCoeffs_periodbc( arr_in, Nmesh, coeffs )
    Integer, Intent(in) :: Nmesh
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: arr_in
    Real(8), Dimension(0:Nmesh-1), Intent(out) :: coeffs

    ! Local variables
    Integer :: i, k
    Real(8), Allocatable, Dimension(:) :: c_plus
    Real(8), Allocatable, Dimension(:) :: c_minus
    
    ! Allocate arrays
    Allocate( c_plus(0:Nmesh-1), c_minus(0:Nmesh-1) )

    ! Compute c_plus(0)
    c_plus(0) = arr_in(0)
    Do i = 1, Nz1
      c_plus(0) = c_plus(0) + arr_in(Nmesh-i)*z1**i  
    End Do
    ! Compute c_plus(1:Nmesh-1)
    Do k = 1, Nmesh-1
      c_plus(k) = arr_in(k) + z1*c_plus(k-1)
    End Do

    ! Compute c_minus(Nmesh-1)
    c_minus(Nmesh-1) = c_plus(Nmesh-1)
    Do i = 1, Nz1
      c_minus(Nmesh-1) = c_minus(Nmesh-1) + c_plus(i-1)*z1**i 
    End Do
    c_minus(Nmesh-1) = - z1*c_minus(Nmesh-1)
    ! Compute c_minus(0:Nmesh-2)
    Do k = Nmesh-2, 0, -1
      c_minus(k) = z1*( c_minus(k+1) - c_plus(k) )
    End Do
    
    ! Deallocate c_plus
    Deallocate( c_plus )
    ! The final c_minus
    coeffs(:)  = 6.0d0*c_minus(:)
    ! Deallocate c_minus
    Deallocate( c_minus )
  End Subroutine InterpCoeffs_periodbc
   
End Module Tools 
  