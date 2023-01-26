Module Parameters 
  Implicit None 

  Integer, Parameter :: Nq_cbsw  = 9
  Integer, Parameter :: Nq_gdw   = 20 
  Integer, Parameter :: Nq_cwgdw = 27 
  Integer, Parameter :: Nq_mw    = 46
  Real(8), Parameter :: Error    = Epsilon(1.0d-15)
  Real(8), Parameter :: Pi       = 3.141592653589793238d0
  Real(8), Parameter :: z11      = -0.36056667770681444d0
  Real(8), Parameter :: z12      = -0.013725466379618684d0
End Module Parameters