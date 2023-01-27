!--------------------------------------------------------------!
! 'FFTCWT' Module contains
! (1) 'rFFTCWT_periodbc' subroutine
!      which use the real-valued wavelet (CBSW, GDW or CWGDW) to
!      perform CWT of the signal with periodic boundary condition.
! (2) 'cFFTCWT_periodbc' subroutine
!      which use the complex-valued wavelet (MW) to perform 
!      CWT of the signal with periodic boundary condition.
!--------------------------------------------------------------!
Module FFTCWT
  Use, intrinsic :: iso_c_binding
  Use Parameters
  Use Tools
  Implicit None
  include 'fftw3.f03'
  Contains
  !--------------------------------------------------------------!
  ! The real Continuous Wavelet Transform
  ! with periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine rFFTCWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet_name, scales, Nscales, CWT )
    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    Real(8), Dimension(1:Nmesh), Intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name

    ! output variables
    Integer, Intent(out) :: Nscales
    Real(8), Allocatable, Dimension(:), Intent(out) :: scales
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: CWT

    ! Local variables
    Type(c_ptr) :: plan_fft, plan_ifft
    Integer :: i, s, k, Nlevs
    Real(8) :: Dx, scal, Dk, cw, scal0
    Real(8), Allocatable, Dimension(:) :: Fwavelet
    Complex(c_double_complex), Allocatable, Dimension(:) :: ssignal
    Complex(c_double_complex), allocatable, Dimension(:) :: Fssignal,FCWT_perscale
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)
    ! The frequency interval
    Dk = 2.0d0*Pi/Lbox
    
    Allocate( Fssignal(0:Nmesh-1), ssignal(0:Nmesh-1), FCWT_perscale(0:Nmesh-1), Fwavelet(0:Nmesh-1) )
    ssignal(:) = signal(:)

    ! Fourier Transform of the signal
    plan_fft = fftw_plan_dft(1, [Nmesh], ssignal, Fssignal, FFTW_FORWARD, FFTW_ESTIMATE)
    Call fftw_execute_dft(plan_fft, ssignal, Fssignal)
    Fssignal = Dx*Fssignal

    ! Choose a wavelet
    CBSW_set  = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set   = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    If ( index(CBSW_set, wavelet_name)/= 0 ) Then
      ! Set scales
      cw       = 0.466094761079290d0    ! The relationship between fourier freqs and wavelet scales
      scal0    = cw*Pi/Lbox
      Nlevs    = Nint( Log( Real(Nmesh,8) )/Log(2.0d0) )
      Nscales  = Nsubs*Nlevs
      Allocate( scales(0:Nscales-1) )
      Do i = 0, Nscales-1
        scales(i) = scal0*(2.0d0**(Real(i,8)/Nsubs))
      End Do

      ! Perform CWT 
      Allocate( CWT(0:Nscales-1,0:Nmesh-1) )
      Do s = 0, Nscales-1
        scal = scales(s)
        ! Values of the Fourier transform of the wavelet at each scale
        Do k = 0, Nmesh/2
          Fwavelet(k) = FCBSW(scal, k*Dk)
        End Do
        Do k = Nmesh/2+1, Nmesh-1
          Fwavelet(k) = FCBSW(scal, (k-Nmesh)*Dk)
        End Do
        ! Convolution is just the multiplication in Fourier space
        FCWT_perscale = Fssignal*Real(Fwavelet,c_double)
        ! Inverse FFT
        plan_ifft = fftw_plan_dft(1, [Nmesh], FCWT_perscale, ssignal, FFTW_BACKWARD, FFTW_ESTIMATE)
        Call fftw_execute_dft(plan_ifft, FCWT_perscale, ssignal)
        CWT(s,:) = Real(ssignal/Lbox,8)
        Call fftw_destroy_plan(plan_ifft)
      End Do
    Else If ( index(GDW_set, wavelet_name)/= 0 ) Then
      ! Set scales
      cw       = 2.0d0/Sqrt(5.0d0)    ! The relationship between fourier freqs and wavelet scales
      scal0    = cw*Pi/Lbox
      Nlevs    = Nint( Log( Real(Nmesh,8) )/Log(2.0d0) )
      Nscales  = Nsubs*Nlevs
      Allocate( scales(0:Nscales-1) )
      Do i = 0, Nscales-1
        scales(i) = scal0*(2.0d0**(Real(i,8)/Nsubs))
      End Do

      ! Perform CWT 
      Allocate( CWT(0:Nscales-1,0:Nmesh-1) )
      Do s = 0, Nscales-1
        scal = scales(s)
        ! Values of the Fourier transform of the wavelet at each scale
        Do k = 0, Nmesh/2
          Fwavelet(k) = FGDW(scal, k*Dk)
        End Do
        Do k = Nmesh/2+1, Nmesh-1
          Fwavelet(k) = FGDW(scal, (k-Nmesh)*Dk)
        End Do
        ! Convolution is just the multiplication in Fourier space
        FCWT_perscale = Fssignal*Real(Fwavelet,c_double)
        ! Inverse FFT
        plan_ifft = fftw_plan_dft(1, [Nmesh], FCWT_perscale, ssignal, FFTW_BACKWARD, FFTW_ESTIMATE)
        Call fftw_execute_dft(plan_ifft, FCWT_perscale, ssignal)
        CWT(s,:) = Real(ssignal/Lbox,8)
        Call fftw_destroy_plan(plan_ifft)
      End Do
    Else If ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      ! Set scales
      cw       = 0.428218886729052d0    ! The relationship between fourier freqs and wavelet scales
      scal0    = cw*Pi/Lbox
      Nlevs    = Nint( Log( Real(Nmesh,8) )/Log(2.0d0) )
      Nscales  = Nsubs*Nlevs
      Allocate( scales(0:Nscales-1) )
      Do i = 0, Nscales-1
        scales(i) = scal0*(2.0d0**(Real(i,8)/Nsubs))
      End Do

      ! Perform CWT 
      Allocate( CWT(0:Nscales-1,0:Nmesh-1) )
      Do s = 0, Nscales-1
        scal = scales(s)
        ! Values of the Fourier transform of the wavelet at each scale
        Do k = 0, Nmesh/2
          Fwavelet(k) = FCWGDW(scal, k*Dk)
        End Do
        Do k = Nmesh/2+1, Nmesh-1
          Fwavelet(k) = FCWGDW(scal, (k-Nmesh)*Dk)
        End Do
        ! Convolution is just the multiplication in Fourier space
        FCWT_perscale = Fssignal*Real(Fwavelet,c_double)
        ! Inverse FFT
        plan_ifft = fftw_plan_dft(1, [Nmesh], FCWT_perscale, ssignal, FFTW_BACKWARD, FFTW_ESTIMATE)
        Call fftw_execute_dft(plan_ifft, FCWT_perscale, ssignal)
        CWT(s,:) = Real(ssignal/Lbox,8)
        Call fftw_destroy_plan(plan_ifft)
      End Do
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    ! clean up
    Call fftw_destroy_plan(plan_fft)
    Deallocate( Fssignal, ssignal, FCWT_perscale, Fwavelet )

  End Subroutine rFFTCWT_periodbc

  !--------------------------------------------------------------!
  ! The complex Continuous Wavelet Transform
  ! with periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine cFFTCWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet_name, scales, Nscales, realCWT, imagCWT )
    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    Real(8), Dimension(1:Nmesh), Intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name

    ! output variables
    Integer, Intent(out) :: Nscales
    Real(8), Allocatable, Dimension(:), Intent(out) :: scales
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: realCWT
    Real(8), Allocatable, Dimension(:,:), Intent(out) :: imagCWT


    ! Local variables
    Type(c_ptr) :: plan_fft, plan_ifft
    Integer :: i, s, k, Nlevs
    Real(8) :: Dx, scal, Dk, cw, scal0
    Real(8), Allocatable, Dimension(:) :: Fwavelet
    Complex(c_double_complex), Allocatable, Dimension(:) :: ssignal
    Complex(c_double_complex), allocatable, Dimension(:) :: Fssignal,FCWT_perscale
    Character(len=60) :: MW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)
    ! The frequency interval
    Dk = 2.0d0*Pi/Lbox
    
    Allocate( Fssignal(0:Nmesh-1), ssignal(0:Nmesh-1), FCWT_perscale(0:Nmesh-1), Fwavelet(0:Nmesh-1) )
    ssignal(:) = signal(:)

    ! Fourier Transform of the signal
    plan_fft = fftw_plan_dft(1, [Nmesh], ssignal, Fssignal, FFTW_FORWARD, FFTW_ESTIMATE)
    Call fftw_execute_dft(plan_fft, ssignal, Fssignal)
    Fssignal = Dx*Fssignal

    ! Choose a wavelet
    MW_set  = "mw, MW, MorletW, Morlet wavelet"
    If ( index(MW_set, wavelet_name)/= 0 ) Then
      ! Set scales
      cw       = 0.242640671273266d0    ! The relationship between fourier freqs and wavelet scales
      scal0    = cw*Pi/Lbox
      Nlevs    = Nint( Log( Real(Nmesh,8) )/Log(2.0d0) )
      Nscales  = Nsubs*Nlevs
      Allocate( scales(0:Nscales-1) )
      Do i = 0, Nscales-1
        scales(i) = scal0*(2.0d0**(Real(i,8)/Nsubs))
      End Do

      ! Perform CWT 
      Allocate( realCWT(0:Nscales-1,0:Nmesh-1), imagCWT(0:Nscales-1,0:Nmesh-1) )
      Do s = 0, Nscales-1
        scal = scales(s)
        ! Values of the Fourier transform of the wavelet at each scale
        Do k = 0, Nmesh/2
          Fwavelet(k) = FMW(scal, k*Dk)
        End Do
        Do k = Nmesh/2+1, Nmesh-1
          Fwavelet(k) = FMW(scal, (k-Nmesh)*Dk)
        End Do
        ! Convolution is just the multiplication in Fourier space
        FCWT_perscale = Fssignal*Real(Fwavelet,c_double)
        ! Inverse FFT
        plan_ifft = fftw_plan_dft(1, [Nmesh], FCWT_perscale, ssignal, FFTW_BACKWARD, FFTW_ESTIMATE)
        Call fftw_execute_dft(plan_ifft, FCWT_perscale, ssignal)
        realCWT(s,:) = Dble(ssignal/Lbox)
        imagCWT(s,:) = -Dimag(ssignal/Lbox)
        Call fftw_destroy_plan(plan_ifft)
      End Do
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    ! clean up
    Call fftw_destroy_plan(plan_fft)
    Deallocate( Fssignal, ssignal, FCWT_perscale, Fwavelet )

  End Subroutine cFFTCWT_periodbc

End Module FFTCWT
