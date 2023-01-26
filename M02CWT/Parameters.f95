Module Parameters 
  Implicit None 

  Real(8), Parameter :: chi_cbsw  = 3d0
  Integer, Parameter :: Nd_cbsw   = 1

  Real(8), Parameter :: chi_gdw   = 12d0 
  Integer, Parameter :: Nd_gdw    = 21

  Real(8), Parameter :: chi_cwgdw = 8d0 
  Integer, Parameter :: Nd_cwgdw  = 27

  Real(8), Parameter :: chi_mw    = 7.5d0 
  Integer, Parameter :: Nd_mw     = 43

  Real(8), Parameter :: Error     = Epsilon(1.0d-15)
  Real(8), Parameter :: Pi        = 3.1415926535897932384626433d0
  Real(8), Parameter :: z1        = -2.0d0 + Sqrt(3.0d0)
  Real(8), Parameter :: z1_2      = z1**2
  Integer, Parameter :: Nz1       = -Nint( 16.0d0*Log(10.0d0) / Log(Abs(z1)) )
End Module Parameters