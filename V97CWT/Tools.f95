Module Tools
  Use Parameters
  
  Implicit None
  Contains
    !--------------------------------------------------------------!
    ! f0_initial is equal to the convolution between
    ! the signal and the cubic B-spline
    !--------------------------------------------------------------!
    Function f0_initial( signal, Nmesh )
      ! Calling variables
      Integer, Intent(in) :: Nmesh
      Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
  
      ! Local variables
      Real(8), Dimension(0:Nmesh-1) :: f0_initial
      Real(8) :: inv6
  
      inv6 = 1.0d0/6.0d0
  
      f0_initial(0)         = signal(Nmesh-1) + 4.0d0*signal(0) + signal(1)
      f0_initial(1:Nmesh-2) = signal(0:Nmesh-3) + 4.0d0*signal(1:Nmesh-2) + signal(2:Nmesh-1)
      f0_initial(Nmesh-1)   = signal(Nmesh-2) + 4.0d0*signal(Nmesh-1) + signal(0)
      f0_initial(:)         = inv6*f0_initial(:)
    End Function f0_initial
  
    !--------------------------------------------------------------!
    ! h fir filtering
    !--------------------------------------------------------------!
    Function h_fir_filtering(fi, Nmesh, power)
      ! Calling variables
      Integer, Intent(in) :: Nmesh
      Integer, Intent(in) :: power
      Real(8), Dimension(0:Nmesh-1), Intent(in) :: fi
  
      ! Local variables
      Integer :: i, i1, i2, i3, i4, p, twop, twop_1
      Real(8) :: fi23, fi14
      Real(8), Dimension(0:Nmesh-1) :: h_fir_filtering
      
      p                  = power
      twop               = 2**p
      twop_1             = 2**(p+1)
      Do i = 0, Nmesh-1
        i1                 = Modulo(i-twop_1, Nmesh)
        i2                 = Modulo(i-twop, Nmesh)
        i3                 = Modulo(i+twop, Nmesh)
        i4                 = Modulo(i+twop_1, Nmesh)
        fi23               = fi(i2)+fi(i3)
        fi14               = fi(i1)+fi(i4)
        h_fir_filtering(i) = ( 6.0d0*fi(i)+ 4.0d0*fi23 + fi14 ) 
      End Do
      h_fir_filtering(:) = 0.125d0*h_fir_filtering(:)
    End Function h_fir_filtering
  
    !--------------------------------------------------------------!
    ! Recursive IIR filter with periodic boundary condition
    !--------------------------------------------------------------!
    Function q12_iir_filtering( fi, Nmesh, power, z1)
      ! Calling variables
      Integer, Intent(in) :: Nmesh
      Integer, Intent(in) :: power
      Real(8), Intent(in) :: z1
      Real(8), Dimension(0:Nmesh-1), Intent(in) :: fi
  
      ! Local varibales
      Integer :: i, j, p, twop, Nz1
      Integer, Allocatable, Dimension(:) :: jj, jj_2p, indx
      Real(8), Allocatable, Dimension(:) :: c_plus, z1_pj
      Real(8), Dimension(0:Nmesh-1) :: q12_iir_filtering
  
      ! Allocate arrays
      Allocate( c_plus(0:Nmesh-1) )
  
      ! Initialize c_plus
      p         = power 
      twop      = 2**p
      Nz1       = -Nint( 16.0d0*Log(10.0d0) / Log(Abs(z1)) )
      Allocate( jj(0:Nz1), jj_2p(0:Nz1), indx(0:Nz1), z1_pj(0:Nz1) )
      c_plus(:) = 0.0d0
      jj        = [ (j, j = 0, Nz1) ]
      jj_2p     = jj*twop
      z1_pj     = z1**jj
      Do i = 0, twop-1
        indx      = Modulo(i-jj_2p, Nmesh)
        c_plus(i) = Sum( z1_pj*fi(indx) )
      End Do
      ! Compute c_plus
      Do i = twop, Nmesh-1
        c_plus(i) = z1*c_plus(i-twop) + fi(i)
      End Do
      
      Deallocate( jj, jj_2p, indx, z1_pj )
      Allocate( jj(0:Nz1-1), jj_2p(0:Nz1-1), indx(0:Nz1-1), z1_pj(0:Nz1-1) )
      ! Initialize c_minus
      q12_iir_filtering(:) = 0.0d0
      jj        = [ (j, j = 0, Nz1-1) ]
      jj_2p     = jj*twop
      z1_pj     = z1**(jj+1)
      Do i = Nmesh-1, Nmesh-twop, -1
        indx       = Modulo( i + jj_2p, Nmesh)
        q12_iir_filtering(i) = -Sum( z1_pj*c_plus(indx) )
      End Do
      Deallocate( jj, jj_2p, indx, z1_pj )
      ! Compute c_minus
      Do i = Nmesh-twop-1, 0, -1
        q12_iir_filtering(i) = z1*( q12_iir_filtering(i+twop) - c_plus(i) )
      End Do
      Deallocate( c_plus )
    End Function q12_iir_filtering
    
    !--------------------------------------------------------------!
    ! qj filter coefficients for the CBSW
    !--------------------------------------------------------------!
    Function qj_cbsw(scal_dimls0, Nsubs, j, k)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: j
      Integer, Intent(in) :: k
      Real(8), Intent(in) :: scal_dimls0
      
  
      ! Local variables
      Real(8) :: ss, invss, kk, cN, aa, bb, p1, p2
      Real(8) :: qj_cbsw
  
      ss    = scal_dimls0/( 2.0**( Real(j,8)/Real(Nsubs,8) ) )
      invss = 1.0d0/ss
      kk    = Real(k,8)
      cN    = Sqrt(30.0d0/31.0d0)
  
      aa = ss*(kk+0.5d0)
      bb = ss*(kk-0.5d0)
      p1 = BS4th( aa-0.5d0 ) - BS4th( bb-0.5d0 )
      p2 = BS4th( bb+0.5d0 ) - BS4th( aa+0.5d0 )
  
      qj_cbsw = cN*invss*(p1+p2)
    End Function qj_cbsw

    !--------------------------------------------------------------!
    ! qj filter coefficients for the GDW
    !--------------------------------------------------------------!
    Function qj_gdw(scal_dimls0, Nsubs, j, k)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: j
      Integer, Intent(in) :: k
      Real(8), Intent(in) :: scal_dimls0
      
  
      ! Local variables
      Real(8) :: ss, ss2, kk, cN, inv16, p1, p2, p12, p22 
      Real(8) :: qj_gdw
  
      ss    = scal_dimls0/( 2.0**( Real(j,8)/Real(Nsubs,8) ) )
      ss2   = ss**2
      kk    = Real(k,8)
      inv16 = 1.0d0/16.0d0
  
      cN     = 1.0d0/Sqrt( Sqrt(18.0*Pi) )
      p1     = 1.0d0+2.0d0*kk
      p2     = 1.0d0-2.0d0*kk
      p12    = p1**2
      p22    = p2**2
  
      qj_gdw = cN*( p1*Exp(-inv16*p12*ss2) + p2*Exp(-inv16*p22*ss2) )
    End Function qj_gdw
    
    !--------------------------------------------------------------!
    ! qj filter coefficients for the CWGDW
    !--------------------------------------------------------------!
    Function qj_cwgdw(scal_dimls0, Nsubs, j, k)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: j
      Integer, Intent(in) :: k
      Real(8), Intent(in) :: scal_dimls0
  
      ! Local variables
      Real(8) :: ss, kk, cN, p1, p1s, p1s2, p2, p2s, p2s2
      Real(8) :: qj_cwgdw
  
      ss    = scal_dimls0/( 2.0**( Real(j,8)/Real(Nsubs,8) ) )
      kk    = Real(k, 8)
      cN    = Exp(0.5d0)*Sqrt( 8.0d0/(1.0d0+5.0d0*Exp(1.0d0)) )/Sqrt( Sqrt(Pi) )
      p1    = kk+0.5d0
      p2    = kk-0.5d0 
      p1s   = p1*ss
      p2s   = p2*ss
      p1s2  = p1s**2
      p2s2  = p2s**2
  
      qj_cwgdw = cN*( p1*Cos(p1s)*Exp(-0.5d0*p1s2) - p2*Cos(p2s)*Exp(-0.5*p2s2) )
    End Function qj_cwgdw
    
    !--------------------------------------------------------------!
    ! qj filter coefficients for the MW
    !--------------------------------------------------------------!
    Function qj_mw(scal_dimls0, Nsubs, j, k)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: j
      Integer, Intent(in) :: k
      Real(8), Intent(in) :: scal_dimls0
  
      ! Local variables
      Real(8) :: cN, ss, invss, p1, p2, p3, kk, invr2
      Complex(8) :: p4 
      Complex(8) :: qj_mw
  
      cN    = Sqrt( 0.5d0*Sqrt(Pi)/(1.0d0-2.0d0*Exp(4.0d0)+Exp(16.0d0)) )
      ss    = scal_dimls0/( 2.0**( Real(j,8)/Real(Nsubs,8) ) )
      invss = 1.0d0/ss
      invr2 = 1.0/Sqrt(2.0d0)
      kk    = Real(k, 8)
      p1    = (kk-0.5d0)*ss
      p2    = (kk+0.5d0)*ss
      p3    = Erf(p1*invr2)+Erfc(p2*invr2)
      p4    = cErf(p2*invr2, 4.0d0*invr2) + cErfc(p1*invr2, 4.0d0*invr2)
  
      qj_mw = cN*invss*(p3+p4-2.0d0)
    End Function qj_mw
  
    !--------------------------------------------------------------!
    ! If the qj coefficients are even, 
    ! then we use this qj FIR filter.
    !--------------------------------------------------------------!
    Function qj_fir_filtering_even(scalesperLevs, qj, Nq, Nsubs, ffi, Nmesh, power)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: Nq
      Integer, Intent(in) :: Nmesh
      Integer, Intent(in) :: power
      Real(8), Dimension(1:Nsubs), Intent(in) :: scalesperLevs
      Real(8), Dimension(0:Nq,1:Nsubs), Intent(in) :: qj
      Real(8), Dimension(0:Nmesh-1), Intent(in) :: ffi
  
      ! Local variables
      Integer :: i, j, k, p, twop
      Integer, Dimension(1:Nq) :: kk, indx1, indx2, kk_2p
      Real(8), Dimension(1:Nsubs) :: sqrt_scals
      Real(8), Dimension(0:Nmesh-1,1:Nsubs) :: qj_fir_filtering_even
  
      p                   = power
      twop                = 2**p
      kk                  = [ (k, k = 1, Nq) ]
      kk_2p               = kk*twop
      
      sqrt_scals = Sqrt( scalesperLevs )
      Do j = 1, Nsubs
        Do i = 0, Nmesh-1
          indx1                      = Modulo(i+kk_2p, Nmesh)
          indx2                      = Modulo(i-kk_2p, Nmesh)
          qj_fir_filtering_even(i,j) = qj(0,j)*ffi(i)+Sum( qj(1:Nq,j)*( ffi(indx1)+ffi(indx2) ) )
        End Do
        qj_fir_filtering_even(:,j) = sqrt_scals(j)*qj_fir_filtering_even(:,j)
      End Do
    End Function qj_fir_filtering_even
    
    !--------------------------------------------------------------!
    ! If the qj coefficients are odd, 
    ! then we use this qj FIR filter.
    !--------------------------------------------------------------!
    Function qj_fir_filtering_odd(scalesperLevs, qj, Nq, Nsubs, ffi, Nmesh, power)
      ! Calling variables
      Integer, Intent(in) :: Nsubs
      Integer, Intent(in) :: Nq
      Integer, Intent(in) :: Nmesh
      Integer, Intent(in) :: power
      Real(8), Dimension(1:Nsubs), Intent(in) :: scalesperLevs
      Real(8), Dimension(0:Nq,1:Nsubs), Intent(in) :: qj
      Real(8), Dimension(0:Nmesh-1), Intent(in) :: ffi
  
      ! Local variables
      Integer :: i, j, k, p, twop
      Integer, Dimension(1:Nq) :: kk, indx1, indx2, kk_2p
      Real(8), Dimension(1:Nsubs) :: sqrt_scals
      Real(8), Dimension(0:Nmesh-1,1:Nsubs) :: qj_fir_filtering_odd
  
      p                   = power
      twop                = 2**p
      kk                  = [ (k, k = 1, Nq) ]
      kk_2p               = kk*twop

      sqrt_scals = Sqrt( scalesperLevs )
      Do j = 1, Nsubs
        Do i = 0, Nmesh-1
          indx1                      = Modulo(i+kk_2p, Nmesh)
          indx2                      = Modulo(i-kk_2p, Nmesh)
          qj_fir_filtering_odd(i,j)  = Sum( qj(1:Nq,j)*( ffi(indx2)-ffi(indx1) ) )
        End Do
        qj_fir_filtering_odd(:,j) = sqrt_scals(j)*qj_fir_filtering_odd(:,j)
      End Do
    End Function qj_fir_filtering_odd
    
    
    !--------------------------------------------------------------! 
    ! The B-spline of degree 4
    !--------------------------------------------------------------!
    Function BS4th( var )
      Real(8), Intent(in) :: var
      Real(8) :: x, BS4th
  
      x = Abs(var)
      If ( (x-0.5d0)<=Error ) Then
        BS4th = 115.0d0/192.0d0-5.0d0*x**2/8.0d0+x**4/4.0d0
      Else if ( ((x-0.5d0)>Error).And.((x-1.5d0)<=Error) ) Then
        BS4th = (55.0d0+20.0d0*x-120.0d0*x**2+80.0d0*x**3-16.0d0*x**4)/96.0d0
      Else if ( ((x-1.5d0)>Error).And.((x-2.5d0)<=Error) ) Then
        BS4th = (5.0d0-2.0d0*x)**4/384.0d0
      Else
        BS4th = 0.0d0
      End If
    End Function BS4th
    
    !--------------------------------------------------------------!
    ! The complex Erfc function
    !--------------------------------------------------------------!
    Function cErfc(xi, yi)
      Real(8), Intent(in) :: xi
      Real(8), Intent(in) :: yi
  
      Real(8) :: u, v, cosxy, sinxy, rerfc, ierfc
      Complex(8) ::  cErfc
      Logical :: flag 
      
      Call wofz(-yi, xi, u, v, flag)
  
      cosxy = cos(2.0*xi*yi)
      sinxy = sin(2.0*xi*yi)
      rerfc = (u*cosxy+v*sinxy)*exp(yi**2-xi**2)
      ierfc = -(u*sinxy-v*cosxy)*exp(yi**2-xi**2)
  
      cErfc = Complex(rerfc,ierfc)
    
    End Function cErfc
    
    !--------------------------------------------------------------!
    ! The complex Erf function
    !--------------------------------------------------------------!
    Function cErf(xi, yi)
      Real(8), Intent(in) :: xi
      Real(8), Intent(in) :: yi
      Complex(8) ::  cErf
  
      cErf = 1.0d0 - cErfc(xi,yi)
    End Function cErf
  
    !--------------------------------------------------------------!
    ! the subroutine computeing 'EXP(-Z**2)*ERFC(-I*Z)', which is
    ! drawn from https://netlib.org/toms/680.gz
    !--------------------------------------------------------------!
    Subroutine wofz(xi, yi, u, v, flag)
      !
      !  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
      !  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
      !  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
      !  MEANS SQRT(-1).
      !  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
      !  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
      !  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
      !  OF THE FUNCTION.
      !  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
      !
      !
      !  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
      !     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
      !                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
      !                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
      !                FLOATING-POINT ARITHMETIC
      !     RMAXEXP  = LN(RMAX) - LN(2)
      !     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
      !                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
      !  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
      !  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
      !
      !
      !  PARAMETER LIST
      !     XI     = REAL      PART OF Z
      !     YI     = IMAGINARY PART OF Z
      !     U      = REAL      PART OF W(Z)
      !     V      = IMAGINARY PART OF W(Z)
      !     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
      !              OCCUR OR NOT; TYPE LOGICAL;
      !              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
      !              MEANING :
      !              FLAG=.FALSE. : NO ERROR CONDITION
      !              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
      !                             BECOMES INACTIVE
      !  XI, YI      ARE THE INPUT-PARAMETERS
      !  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
      !
      !  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
      !
      !  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
      !  PUT TO 0 UPON UNDERFLOW;
      !
      !  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
      !  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
      !
        Implicit None
      !
        Real(8), Parameter :: factor=1.12837916709551257388D0
        Real(8), Parameter :: rmaxreal=0.5D+154
        Real(8), Parameter :: rmaxexp=708.503061461606D0
        Real(8), Parameter :: rmaxgoni=3.53711887601422D+15
        
        Real(8), Intent(in) :: xi
        Real(8), Intent(in) :: yi
        Real(8), Intent(out) :: u 
        Real(8), Intent(out) :: v 
        Logical, Intent(out) :: flag
        
        Logical :: a, b
        Real(8) :: c, daux, h, h2, qlambda, qrho, rx, ry, sx, sy, tx, ty
        Real(8) :: u1, u2, v1, v2, w1, x, y, xabs, xabsq, xaux, xquad, xsum
        Real(8) :: yabs, yquad, ysum
        Integer :: i, j, n, kapn, np1, nu
        
        
      !
        flag = .False.
      !
        xabs = dabs(xi)
        yabs = dabs(yi)
        x = xabs/6.3
        y = yabs/4.4
      !
      !
      !     THE FOLLOWING IF-STATEMENT PROTECTS
      !     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
      !
        If ((xabs>rmaxreal) .Or. (yabs>rmaxreal)) Goto 100
      !
        qrho = x**2 + y**2
      !
        xabsq = xabs**2
        xquad = xabsq - yabs**2
        yquad = 2*xabs*yabs
      !
        a = qrho < 0.085264D0
      !
        If (a) Then
      !
      !  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
      !  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
      !  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
      !  ACCURACY
      !
          qrho = (1-0.85*y)*dsqrt(qrho)
          n = idnint(6+72*qrho)
          j = 2*n + 1
          xsum = 1.0/j
          ysum = 0.0D0
          Do i = n, 1, -1
            j = j - 2
            xaux = (xsum*xquad-ysum*yquad)/i
            ysum = (xsum*yquad+ysum*xquad)/i
            xsum = xaux + 1.0/j
          End Do
          u1 = -factor*(xsum*yabs+ysum*xabs) + 1.0
          v1 = factor*(xsum*xabs-ysum*yabs)
          daux = dexp(-xquad)
          u2 = daux*dcos(yquad)
          v2 = -daux*dsin(yquad)
      !
          u = u1*u2 - v1*v2
          v = u1*v2 + v1*u2
      !
        Else
      !
      !  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
      !  CONTINUED FRACTION
      !  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
      !  ACCURACY
      !
      !  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
      !  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
      !  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
      !  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
      !  TO OBTAIN THE REQUIRED ACCURACY
      !  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
      !  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
      !
          If (qrho>1.0) Then
            h = 0.0D0
            kapn = 0
            qrho = dsqrt(qrho)
            nu = idint(3+(1442/(26*qrho+77)))
          Else
            qrho = (1-y)*dsqrt(1-qrho)
            h = 1.88*qrho
            h2 = 2*h
            kapn = idnint(7+34*qrho)
            nu = idnint(16+26*qrho)
          End If
      !
          b = (h>0.0)
      !
          If (b) qlambda = h2**kapn
      !
          rx = 0.0
          ry = 0.0
          sx = 0.0
          sy = 0.0
      !
          Do n = nu, 0, -1
            np1 = n + 1
            tx = yabs + h + np1*rx
            ty = xabs - np1*ry
            c = 0.5/(tx**2+ty**2)
            rx = c*tx
            ry = c*ty
            If ((b) .And. (n<=kapn)) Then
              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda/h2
            End If
          End Do
      !
          If (h==0.0d0) Then
            u = factor*rx
            v = factor*ry
          Else
            u = factor*sx
            v = factor*sy
          End If
      !
          If (yabs==0.0d0) u = dexp(-xabs**2)
      !
        End If
      !
      !  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
      !
        If (yi<0.0) Then
      !
          If (a) Then
            u2 = 2*u2
            v2 = 2*v2
          Else
            xquad = -xquad
      !
      !         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
      !         AGAINST OVERFLOW
      !
            If ((yquad>rmaxgoni) .Or. (xquad>rmaxexp)) Goto 100
      !
            w1 = 2*dexp(xquad)
            u2 = w1*dcos(yquad)
            v2 = -w1*dsin(yquad)
          End If
      !
          u = u2 - u
          v = v2 - v
          If (xi>0.0) v = -v
        Else
          If (xi<0.0) v = -v
        End If
      !
        Return
      !
        100 flag = .True.
        Return
      !
    End Subroutine wofz
   
End Module Tools 
  