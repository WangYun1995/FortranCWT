Module Tools
  Use Parameters
  
  Implicit None
  Contains
  !-------------------------------------------------------------!
  ! The 2nd-order derivative of the CBSW
  !-------------------------------------------------------------!
  Function d2CBSWdx2( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, CN
    Real(8) :: d2CBSWdx2

    CN = Sqrt(30.0d0/31.0d0)
    x  = var
    If ( x<=-3 ) Then
      d2CBSWdx2 = 0.0d0
    Else if ( (x>-3).And.(x<=-2) ) Then 
      d2CBSWdx2 = -CN*(x+3d0)
    Else if ( (x>-2).And.(x<=-1) ) Then 
      d2CBSWdx2 = CN*(9d0+5d0*x)
    Else if ( (x>-1).And.(x<=0) ) Then 
      d2CBSWdx2 = -2*CN*(3d0+5d0*x)
    Else if ( (x>0).And.(x<=1) ) Then 
      d2CBSWdx2 = 2*CN*(-3d0+5d0*x)
    Else if ( (x>1).And.(x<=2) ) Then 
      d2CBSWdx2 = CN*(9d0-5d0*x)
    Else if ( (x>2).And.(x<=3) ) Then 
      d2CBSWdx2 = CN*(-3d0+x)
    Else if ( x>3 ) Then 
      d2CBSWdx2 = 0.0d0
    End If
  End Function d2CBSWdx2
  
  !-------------------------------------------------------------!
  ! The 2nd-order derivative of the GDW
  !-------------------------------------------------------------!
  Function d2GDWdx2( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, CN
    Real(8) :: d2GDWdx2

    CN = -1d0/(4d0*Sqrt(3d0*Sqrt(2d0*Pi)))
    x  = var

    d2GDWdx2 = CN*(12d0-12d0*x**2+x**4)*Exp(-0.25d0*x**2)
  End Function d2GDWdx2

  !-------------------------------------------------------------!
  ! The 2nd-order derivative of the CWGDW
  !-------------------------------------------------------------!
  Function d2CWGDWdx2( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, CN, p1, p2
    Real(8) :: d2CWGDWdx2

    CN = -2d0*Sqrt(2d0/(1d0+5d0*Exp(1d0)))/Sqrt( Sqrt(Pi) )
    x  = var
    p1 = CN*Exp(0.5d0-0.5d0*x**2)
    p2 = (6d0-9d0*x**2+x**4)*Cos(x)+x*(3d0*x**2-10d0)*Sin(x)
    d2CWGDWdx2 = p1*p2
  End Function d2CWGDWdx2

  !-------------------------------------------------------------!
  ! The 2nd-order derivative of the real part of the MW
  !-------------------------------------------------------------!
  Function d2rMWdx2( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, CN, p1, p2
    Real(8) :: d2rMWdx2

    CN = 1d0/Sqrt( (1d0-2d0*Exp(4d0)+Exp(16d0))*Sqrt(Pi) )
    x  = var
    p1 = CN*Exp(-0.5d0*x**2)
    p2 = 1d0-x**2+Exp(8d0)*((x**2-17d0)*Cos(4d0*x)+8d0*x*Sin(4d0*x))
    d2rMWdx2 = p1*p2
  End Function d2rMWdx2

  !-------------------------------------------------------------!
  ! The 2nd-order derivative of the imaginary part of the MW
  !-------------------------------------------------------------!
  Function d2iMWdx2( var )
    Real(8), Intent(in) :: var
    Real(8) :: x, CN, p1, p2
    Real(8) :: d2iMWdx2

    CN = 1d0/Sqrt( (1d0-2d0*Exp(4d0)+Exp(16d0))*Sqrt(Pi) )
    x  = var
    p1 = CN*Exp(8d0-0.5d0*x**2)
    p2 = 8d0*x*Cos(4d0*x)-(x**2-17d0)*Sin(4d0*x)
    d2iMWdx2 = p1*p2
  End Function d2iMWdx2

  !-------------------------------------------------------------!
  ! The B-spline of degree 7
  !-------------------------------------------------------------!
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
      BS7th = (2472.0d0-392.0d0*x-504.0d0*x2-1960.0d0*x3+2520.0d0*x4                  &
              -1176.0d0*x5+252.0d0*x6-21.0d0*x7)/5040.0d0
    Else if ( ((x-2.0d0)>Error).And.((x-3.0d0)<=Error) ) Then
      BS7th = (-1112.0d0+12152.0d0*x-19320.0d0*x2+13720.0d0*x3-5320.0d0*x4            &
              +1176.0d0*x5-140.0d0*x6+7.0d0*x7)/5040.0d0
    Else if ( ((x-3.0d0)>Error).And.((x-4.0d0)<=Error) ) Then
      BS7th = -(x-4.0d0)**7/5040.0d0
    Else
      BS7th = 0.0d0
    End If
  End Function BS7th
  
  !--------------------------------------------------------------!
  ! Coefficients of the cubic B-spline approximation of the signal,
  ! which satisfies the zero bouNchiary coNchiition
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
  ! which satisfies the periodic bouNchiary coNchiition
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
  ! B coefficients
  !--------------------------------------------------------------!
  Subroutine Bcoeffs( wavelet_name, chi, Nchi, chik, Bk )
    ! Calling variables
    Integer, Intent(in) :: Nchi
    Real(8), Intent(in) :: chi
    Real(8), Dimension(-Nchi:Nchi), Intent(out) :: chik 
    Real(8), Dimension(-Nchi:Nchi), Intent(out) :: Bk
    Character(len=*), Intent(in) :: wavelet_name

    ! Local variables
    Integer :: i, k
    Real(8) :: Da
    Real(8), Dimension(-Nchi:Nchi-1) :: alpha3
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set, realMW_set, imagMW_set

    Da  = chi/Real(Nchi,8)
    Do i = -Nchi, Nchi
      chik(i) = Da*Real(i,8)
    End Do
    
    ! Check if there is a wavelet in this set
    CBSW_set   = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set    = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set  = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    realMW_set = "rmw, rMW, real part of Morlet wavelet"
    imagMW_set = "imw, iMW, imaginary part of Morlet wavelet"
    If ( index(CBSW_set, wavelet_name)/= 0 ) Then
      Do k = -Nchi, Nchi-1
        alpha3(k) = (d2CBSWdx2(chik(k+1))-d2CBSWdx2( chik(k) ))/Da
      End Do
    Else if ( index(GDW_set, wavelet_name)/= 0 ) Then
      Do k = -Nchi, Nchi-1
        alpha3(k) = (d2GDWdx2(chik(k+1))-d2GDWdx2( chik(k) ))/Da
      End Do
    Else if ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      Do k = -Nchi, Nchi-1
        alpha3(k) = (d2CWGDWdx2(chik(k+1))-d2CWGDWdx2( chik(k) ))/Da
      End Do
    Else if ( index(realMW_set, wavelet_name)/= 0 ) Then
      Do k = -Nchi, Nchi-1
        alpha3(k) = (d2rMWdx2(chik(k+1))-d2rMWdx2( chik(k) ))/Da
      End Do
    Else if ( index(imagMW_set, wavelet_name)/= 0 ) Then
      Do k = -Nchi, Nchi-1
        alpha3(k) = (d2iMWdx2(chik(k+1))-d2iMWdx2( chik(k) ))/Da
      End Do
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    Bk(-Nchi) = alpha3(-Nchi)
    Do k = 1-Nchi, Nchi-1
      Bk(k) = (alpha3(k) - alpha3(k-1))
    End Do
    Bk(Nchi) = -alpha3(Nchi-1)
  End Subroutine Bcoeffs
   
End Module Tools 
  