# FortranCWT
**The Fortran 95 codes for fast implementation of the Continuous Wavelet Transform (CWT) of the one-dimensional signals.**

If you use our codes or your research is related to our paper, please kindly cite the following paper, in which the fast CWT algorithms are reviewed and compared.

``` bib
@article{Wang2023,
    author = {{Wang}, Yun and {He}, Ping},
    title = "{Comparisons between fast algorithms for the continuous wavelet transform and applications in cosmology: The one-dimensional case}",
    journal = {RAS Techniques and Instruments},
    year = {2023},
    month = {06},
    doi = {10.1093/rasti/rzad020},
    url = {https://doi.org/10.1093/rasti/rzad020},
    note = {rzad020},
    eprint = {https://academic.oup.com/rasti/advance-article-pdf/doi/10.1093/rasti/rzad020/50599483/rzad020.pdf},
}
```


## Tutorial:

Our codes can be used to process 1D signals with periodic boundary conditions as well as zero boundary conditions. In our codes, four types of wavelets are available, namely the cubic B-spline wavelet (CBSW), the Gaussian-derived wavelet (GDW), the cosine-weighted Gaussian-derived wavelet (CW-GDW), and the Morlet wavelet (MW). The former three wavelets are real valued and the last one is complex valued, which are shown in the figure below.
<div align=left><img width="800" height="390" src="https://github.com/WangYun1995/FortranCWT/blob/main/figures/wavelets.png"/>

- The module ```FFTCWT.f95``` contains two subroutines, the ```rFFTCWT_periodbc```, and the ```cFFTCWT_periodbc```.
- The module ```V97CWT.f95``` contains two subroutines, the ```rV97CWT_periodbc```, and the ```cV97CWT_periodbc```.
- The module ```M02CWT.f95``` contains four subroutines, the ```rM02CWT_zerobc```, the ```rM02CWT_periodbc```, the ```cM02CWT_zerobc```, and the ```cM02CWT_periodbc```.
- The module ```A19CWT.f95``` contains four subroutines, the ```rA19CWT_zerobc```, the ```rA19CWT_periodbc```, the ```cA19CWT_zerobc```, and the ```cA19CWT_periodbc```.

The prefix "```r```" ("```c```") indicates that the subroutine uses real-valued (complex-valued) wavelets. The suffix "```periodbc```" ("```zerobc```") indicates that the the signals are assumed to satisfy periodic boundary conditions (zero boundary conditions).

All these subroutines receive five input arguments, which are
- ```signal```, type: real(kind=8), the signal to be analyzed, which is a 1D array.
- ```Nmesh```, type: integer, the number of data points in the array ```signal```.
- ```Lbox```, type: real(kind=8), the spatial length or time duration of ```signal```.
- ```Nsubs```, type: integer, in our codes, the scales are set as $w\propto w_02^{i+j/Nsubs}$, which divides the scales into ```Nlevs``` levels, numbered by $i$ and determined automatically by subroutines; then each level is divided into ```Nsubs```, numbered by $j$. Thus, there is a total of ```Nscales```=```Nlevs*Nsubs``` scales, and ```Nsubs``` controls the scale resolution.
- ```wavelet_name```, type: character, ```cbsw```, ```gdw```, and ```cwgdw``` are available for the subroutines with the prefix ```r```, while ```mw``` is available for the subroutines with the prefix ```c```.

Subroutines with the prefix ```r``` will output three arguments, which are
- ```scales```, type: real(kind=8), the 1D array contains all scale parameters.
- ```Nscales```, type: integer, the number of scale parameters, which satisfies ```Nscales```=```Nlevs*Nsubs```.
- ```CWT```, type: real(kind=8), dimension: (Nscales,Nmesh), the CWT of the signal, which is a 2D array.

Subroutines with the prefix ```c``` will output four arguments, which are
- ```scales```, type: real(kind=8), the 1D array contains all scale parameters.
- ```Nscales```, type: integer, the number of scale parameters, which satisfies ```Nscales```=```Nlevs*Nsubs```.
- ```realCWT```, type: real(kind=8), dimension: (Nscales,Nmesh), the real part of the CWT, which is a 2D array.
- ```imagCWT```, type: real(kind=8), dimension: (Nscales,Nmesh), the imaginary part of the CWT, which is a 2D array.

## Example:

As an example, we use ```rFFTCWT_periodbc``` to perform the CWT of 1D density fields (see Wang & He 2023 ) with periodic boundary conditions, and the program is shown below. Note that the use of the ```FFTCWT``` module requires that ```FFTW3``` be installed beforehand, while other modules do not.
  ``` fortran
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
  ```
As shown in the figure below, the outcome of the program is
<div align=left><img width="800" height="390" src="https://github.com/WangYun1995/FortranCWT/blob/main/figures/dens_cwts.png"/>
  
 
