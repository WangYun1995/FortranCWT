Module Tools
  Use Parameters
  
  Implicit None
  Contains
  !--------------------------------------------------------------!
  ! Fourier transform of the CBSW
  !--------------------------------------------------------------!
  Function FCBSW( scale, wavenum )
    ! Calling variables
    Real(8), Intent(in) :: scale
    Real(8), Intent(in) :: wavenum

    ! Local variables
    Real(8) :: k, cN, p1, p2
    Real(8) :: FCBSW

    k     = wavenum/scale
    If ( (Abs(k)-0d0)<Error ) Then
      FCBSW = 0d0
    Else 
      cN    = 64d0*Sqrt(30d0/31d0)
      p1    = Sin(0.5d0*k)/k
      p2    = cN*(p1**4)*(Sin(0.5d0*k))**2
      FCBSW = p2/Sqrt(scale)
    End If
  End Function FCBSW

  !--------------------------------------------------------------!
  ! Fourier transform of the GDW
  !--------------------------------------------------------------!
  Function FGDW( scale, wavenum )
    ! Calling variables
    Real(8), Intent(in) :: scale
    Real(8), Intent(in) :: wavenum

    ! Local variables
    Real(8) :: k, cN
    Real(8) :: FGDW

    k    = wavenum/scale
    cN   = 8d0*Sqrt(Pi)/Sqrt( Sqrt(18d0*Pi) )
    FGDW = cN*Exp(-k**2)*(k**2)/Sqrt(scale)
  End Function FGDW

  !--------------------------------------------------------------!
  ! Fourier transform of the CWGDW
  !--------------------------------------------------------------!
  Function FCWGDW( scale, wavenum )
    ! Calling variables
    Real(8), Intent(in) :: scale
    Real(8), Intent(in) :: wavenum

    ! Local variables
    Real(8) :: k, cN, p1
    Real(8) :: FCWGDW

    k      = wavenum/scale
    cN     = Sqrt(2d0*Pi)
    cN     = cN*Sqrt( 8d0/( Sqrt(Pi)*(1d0+5d0*Exp(1d0)) ) )
    p1     = 0.5d0*(1d0+k)*Exp(-0.5d0*k**2-k)+0.5d0*(k-1d0)*Exp(k-0.5d0*k**2) 
    FCWGDW = cN*k*p1/Sqrt(scale)
  End Function FCWGDW

  !--------------------------------------------------------------!
  ! Fourier transform of the MW
  !--------------------------------------------------------------!
  Function FMW( scale, wavenum )
    ! Calling variables
    Real(8), Intent(in) :: scale
    Real(8), Intent(in) :: wavenum

    ! Local variables
    Real(8) :: k, cN
    Real(8) :: FMW

    k      = wavenum/scale
    cN     = Sqrt(2d0*Pi)
    cN     = cN/Sqrt( Sqrt(Pi)*(1d0-2d0*Exp(4d0)+Exp(16d0)) )
    FMW = cN*(Exp(4d0*k-0.5d0*k**2)-Exp(-0.5d0*k**2))/Sqrt(scale)
  End Function FMW
   
End Module Tools 
  