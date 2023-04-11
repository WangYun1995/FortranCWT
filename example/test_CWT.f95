Program test
  Use FFTCWT
  Implicit None

  ! Declaring Parameters
  Integer, Parameter :: Nsubs = 12              ! Number of sub-levels of scales
  Integer, Parameter :: Nmesh = 4096            ! Data points of the signal
  Real(8), Parameter :: Lbox  = 1d0             ! Spatial length of the signal

  ! Declaring variables
  Integer :: i, Nscales 
  Real(8) :: t1, t2, Dx
  Real(8), Dimension(0:Nmesh-1) :: signal
  Real(8), Allocatable, Dimension(:) :: scales
  Real(8), Allocatable, Dimension(:,:) :: cwt   ! cwt, realcwt, imagcwt
  Character(len=:), Allocatable :: wavelet

  ! Choose a real valued wavelet, e.g. cbsw, gdw, and cwgdw 
  wavelet = "cwgdw"

  ! Load the test signal 
  Dx = Lbox/Real(Nmesh,8)        
  Open(101, file="path/delta_100.dat")        ! the 1D signal with periodic boundary condition        
  ! Open(101, file="path/delta_1.dat")        ! the 1D signal with periodic boundary condition                                                                                         
  Do i = 0, Nmesh-1                                                
    Read(101,*) signal(i)                                           
  End Do
  Close(101)  


  Call CPU_time(t1)

  ! Perform the CWT
  Call rFFTCWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet, scales, Nscales, cwt )
  
  Call CPU_time(t2)
  
  ! Output the CPU time
  Write (*,*) (t2-t1)

  ! Output the CWT of the signal
  Open(102, file="cwt.dat")
  Do i = 0, Nscales-1
    Write (102,*) cwt(i,:)
  End Do
  Close(102)

  ! Output the scales 
  Open(103, file="scales.dat")
  Do i = 0, Nscales-1
    Write (103,*) scales(i)
  End Do
  Close(103)

End Program test