!----------------------------------------------------------------------------------
Program test
  Use M02CWT
  Implicit None

  ! Declaring Parameters
  Integer, Parameter :: Nsubs = 12       ! Number of sub levels of scales
  Integer, Parameter :: Nmesh = 512      ! Sampling number of the signal
  Real(8), Parameter :: Lbox  = 2*Pi     ! Physical length of the signal
  ! Real(8), Parameter :: Lbox  = 12d0     ! Physical length of the signal
    
  ! Declaring variables
  Integer :: i, Nscales 
  Real(8) :: t1, t2, Dx
  Real(8), Dimension(0:Nmesh-1) :: signal
  Real(8), Allocatable, Dimension(:) :: scales
  Real(8), Allocatable, Dimension(:,:) :: realcwt, imagcwt
  Character(len=:), Allocatable :: wavelet

  ! Choose a complex valued wavelet, e.g. mw
  wavelet = "mw"

  ! Load the test signal 
  Dx = Lbox/Real(Nmesh,8)        
  Open(101, file="path/periodbc_signal_L2pi_N512.dat")   ! the signal with periodic boundary condition
  ! Open(101, file="path/zerobc_signal_L12_N512.dat")    ! the signal with zero boundary condition                                                                                             
  Do i = 0, Nmesh-1                                                
    Read(101,*) signal(i)                                           
  End Do
  Close(101)  

  Call CPU_time(t1)
  ! Perform the CWT
  Call cM02CWT_periodbc( signal, Nmesh, Lbox, Nsubs, wavelet, scales, Nscales, realcwt, imagcwt )
  ! Perform the CWT
  ! Call cM02CWT_zerobc( signal, Nmesh, Lbox, Nsubs, wavelet, scales, Nscales, realcwt, imagcwt )
  Call CPU_time(t2)

  ! Output the CPU time
  Write (*,*) (t2-t1)
  
  ! Output the real part of the CWT
  Open(102, file="realcwt.dat")
  Do i = 0, Nscales-1
    Write (102,*) realcwt(i,:)
  End Do
  Close(102)

  ! Output the imaginary part of the CWT
  Open(103, file="imagcwt.dat")
  Do i = 0, Nscales-1
    Write (103,*) imagcwt(i,:)
  End Do
  Close(103)
  
  ! Output the scales 
  Open(104, file="scales.dat")
  Do i = 0, Nscales-1
    Write (104,*) scales(i)
  End Do
  Close(104)

End Program test