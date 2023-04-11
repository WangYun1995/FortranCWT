# FortranCWT
**The Fortran 95 codes for fast implementation of the Continuous Wavelet Transform (CWT) of the one-dimensional signals.**

If you use our codes or your research is related to our paper, please kindly cite the following paper, in which the fast CWT algorithms are reviewed and compared.

*Wang, Y., & He, P. (2023). Comparisons between fast algorithms for the continuous wavelet transform and applications in cosmology: the one-dimensional case. arXiv preprint arXiv:2302.03909.*

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

## Examples:
