!--------------------------------------------------------------!
! 'A19CWT' Module contains
! (1) 'rV97CWT_periodbc' subroutine
!      which use the real-valued wavelet (CBSW, GDW or CWGDW) to
!      perform CWT of the signal with periodic boundary condition.
! (2) 'cV97CWT_periodbc' subroutine
!      which use the complex-valued wavelet (MW) to perform 
!      CWT of the signal with periodic boundary condition.
!--------------------------------------------------------------!
Module V97CWT
  Use Parameters
  Use Tools
  Implicit None
  Contains

  !--------------------------------------------------------------!
  ! The real Continuous Wavelet Transform
  ! with periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine rV97CWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet_name, scales, Nscales, CWT )
    ! Input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name

    ! Local variables
    Integer :: i, j, k, j1, j2, Nlevs, Nq
    Real(8) :: Dx, s0, scal0, cw
    Real(8), Dimension(0:Nmesh-1) :: f0, f1, f1_
    Real(8), Dimension(1:Nsubs) :: scalesperLevs
    Real(8), Allocatable, Dimension(:,:) :: qj_coeffs
    Real(8), Allocatable, Dimension(:,:) :: CCWT
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set

    ! Output
    Integer, Intent(out) :: Nscales
    Real(8), Allocatable, Dimension(:),Intent(out) :: scales
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: CWT
    
    !Sampling interval of signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    CBSW_set  = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set   = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    If ( index(GDW_set, wavelet_name)/= 0 ) Then
      cw      = 2.0d0/Sqrt(5.0d0)     ! The relationship between Fourier freq and wavelet scale
      scal0   = 1.34d0*cw*Real(Nmesh,8)/Lbox
      s0      = scal0*Dx
      Nlevs   = Nint( Log( s0*Real(Nmesh,8)/(cw*Pi) )/Log(2.0d0) )
      Nscales = Nsubs*Nlevs
      Nq      = Nq_gdw
      Allocate( scales(0:Nscales-1), qj_coeffs(0:Nq,0:Nsubs-1) )
      Do i = Nscales-1, 0, -1
        scales(Nscales-1-i) = scal0/(2.0d0**(Real(i,8)/Nsubs))
      End Do
      Do j = 0, Nsubs-1   
        Do k = 0, Nq
          qj_coeffs(k,j) = qj_gdw(s0, Nsubs, j, k)
        End Do
      End Do

    Else if ( index(CBSW_set, wavelet_name)/= 0 ) Then
      cw      = 0.466094761079290d0
      scal0   = 1.34d0*cw*Real(Nmesh,8)/Lbox
      s0      = scal0*Dx
      Nlevs   = Nint( Log( s0*Real(Nmesh,8)/(cw*Pi) )/Log(2.0d0) )
      Nscales = Nsubs*Nlevs
      Nq      = Nq_cbsw
      Allocate( scales(0:Nscales-1), qj_coeffs(0:Nq,0:Nsubs-1) )
      Do i = Nscales-1, 0, -1
        scales(Nscales-1-i) = scal0/(2.0d0**(Real(i,8)/Nsubs))
      End Do
      Do j = 0, Nsubs-1   
        Do k = 0, Nq
          qj_coeffs(k,j) = qj_cbsw(s0, Nsubs, j, k)
        End Do
      End Do
    
    Else if ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      cw      = 0.428218886729052d0
      scal0   = 1.34d0*cw*Real(Nmesh,8)/Lbox
      s0      = scal0*Dx
      Nlevs   = Nint( Log( s0*Real(Nmesh,8)/(cw*Pi) )/Log(2.0d0) )
      Nscales = Nsubs*Nlevs
      Nq      = Nq_cwgdw
      Allocate( scales(0:Nscales-1), qj_coeffs(0:Nq,0:Nsubs-1) )
      Do i = Nscales-1, 0, -1
        scales(Nscales-1-i) = scal0/(2.0d0**(Real(i,8)/Nsubs))
      End Do
      Do j = 0, Nsubs-1   
        Do k = 0, Nq
          qj_coeffs(k,j) = qj_cwgdw(s0, Nsubs, j, k)
        End Do
      End Do

    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Initialization: i = 0
    Allocate( CCWT(0:Nmesh-1,0:Nscales-1) )
    f0(:)            = f0_initial( signal(:), Nmesh )
    f1_(:)           = q12_iir_filtering( f0(:), Nmesh, 0, z11)
    f1_(:)           = 384.0d0*q12_iir_filtering(f1_(:), Nmesh, 0, z12)
    j1               = Nscales-1
    j2               = Nscales-1-(Nsubs-1)
    scalesperLevs    = scales(j1:j2:-1)
    CCWT(:,j1:j2:-1) = qj_fir_filtering_even(scalesperLevs, qj_coeffs(:,:), Nq, Nsubs, f1_(:), Nmesh, 0)

    ! perform CWT at i>=1
    Do i = 1, Nlevs-1
      j1               = Nscales-1-i*Nsubs
      j2               = Nscales-1-((i+1)*Nsubs-1)
      scalesperLevs    = scales(j1:j2:-1)
      f1(:)            = h_fir_filtering(f0(:), Nmesh, i-1)
      f1_(:)           = q12_iir_filtering( f1(:), Nmesh, i, z11)
      f1_(:)           = 384.0d0*q12_iir_filtering(f1_(:), Nmesh, i, z12)
      CCWT(:,j1:j2:-1) = qj_fir_filtering_even(scalesperLevs, qj_coeffs(:,:), Nq, Nsubs, f1_(:), Nmesh, i)
      f0(:)            = f1(:)
    End Do
    Allocate( CWT(0:Nscales-1,0:Nmesh-1) )
    CWT(:,:) = Dx*Transpose(CCWT(:,:))
    Deallocate( CCWT )
  
  End Subroutine rV97CWT_periodbc

  !--------------------------------------------------------------!
  ! The complex Continuous Wavelet Transform
  ! with periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine cV97CWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet_name, scales, Nscales, realCWT, imagCWT )
    ! Input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name

    ! Local variables
    Integer :: i, j, k, j1, j2, Nlevs, Nq
    Real(8) :: Dx, s0, scal0, cw
    Real(8), Dimension(0:Nmesh-1) :: f0, f1, f1_
    Real(8), Dimension(1:Nsubs) :: scalesperLevs
    Real(8), Allocatable, Dimension(:,:) :: qj_coeffs_r, qj_coeffs_i
    Real(8), Allocatable, Dimension(:,:) :: CCWT_r, CCWT_i
    Complex(8) :: qj_coeffs_c
    Character(len=80) :: MW_set

    ! Output
    Integer, Intent(out) :: Nscales
    Real(8), Allocatable, Dimension(:),Intent(out) :: scales
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: realCWT
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: imagCWT

    !Sampling interval of signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    MW_set  = "mw, MW, MorletW, Morlet wavelet"

    If ( index(MW_set, wavelet_name)/= 0 ) Then
      cw      = 0.242640671273266d0         ! The relationship between Fourier freq and wavelet scale
      scal0   = 1.34d0*cw*Real(Nmesh,8)/Lbox
      s0      = scal0*Dx
      Nlevs   = Nint( Log( s0*Real(Nmesh,8)/(cw*Pi) )/Log(2.0d0) )
      Nscales = Nsubs*Nlevs
      Nq      = Nq_mw
      Allocate( scales(0:Nscales-1), qj_coeffs_r(0:Nq,0:Nsubs-1), qj_coeffs_i(0:Nq,0:Nsubs-1) )
      Do i = Nscales-1, 0, -1
        scales(Nscales-1-i) = scal0/(2.0d0**(Real(i,8)/Nsubs))
      End Do
      Do j = 0, Nsubs-1   
        Do k = 0, Nq
          qj_coeffs_c      = qj_mw(s0, Nsubs, j, k)
          qj_coeffs_r(k,j) = Dble(qj_coeffs_c)
          qj_coeffs_i(k,j) = Dimag(qj_coeffs_c)
        End Do
      End Do

    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Compute the realCWT and imagCWT
    ! Initialization: i = 0
    Allocate( CCWT_r(0:Nmesh-1,0:Nscales-1), CCWT_i(0:Nmesh-1,0:Nscales-1) )
    f0(:)              = f0_initial( signal(:), Nmesh )
    f1_(:)             =  q12_iir_filtering( f0(:), Nmesh, 0, z11)
    f1_(:)             = 384.0d0*q12_iir_filtering(f1_(:), Nmesh, 0, z12)
    j1                 = Nscales-1
    j2                 = Nscales-1-(Nsubs-1)
    scalesperLevs      = scales(j1:j2:-1)
    CCWT_r(:,j1:j2:-1) = qj_fir_filtering_even(scalesperLevs, qj_coeffs_r(:,:), Nq, Nsubs, f1_(:), Nmesh, 0)
    CCWT_i(:,j1:j2:-1) = qj_fir_filtering_odd(scalesperLevs, qj_coeffs_i(:,:), Nq, Nsubs, f1_(:), Nmesh, 0)

    ! perform CWT at i>=1
    Do i = 1, Nlevs-1
      j1                 = Nscales-1-i*Nsubs
      j2                 = Nscales-1-((i+1)*Nsubs-1)
      scalesperLevs      = scales(j1:j2:-1)
      f1(:)              = h_fir_filtering(f0(:), Nmesh, i-1)
      f1_(:)             = q12_iir_filtering( f1(:), Nmesh, i, z11)
      f1_(:)             = 384.0d0*q12_iir_filtering(f1_(:), Nmesh, i, z12)
      CCWT_r(:,j1:j2:-1) = qj_fir_filtering_even(scalesperLevs, qj_coeffs_r(:,:), Nq, Nsubs, f1_(:), Nmesh, i)
      CCWT_i(:,j1:j2:-1) = qj_fir_filtering_odd(scalesperLevs, qj_coeffs_i(:,:), Nq, Nsubs, f1_(:), Nmesh, i)
      f0(:)              = f1(:)
    End Do
    Allocate( realCWT(0:Nscales-1,0:Nmesh-1) )
    realCWT(:,:) = Dx*Transpose(CCWT_r(:,:))
    Deallocate( CCWT_r )

    Allocate( imagCWT(0:Nscales-1,0:Nmesh-1) )
    imagCWT(:,:) = Dx*Transpose(CCWT_i(:,:))
    Deallocate( CCWT_i )

  End Subroutine cV97CWT_periodbc

End Module V97CWT