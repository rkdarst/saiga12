# Richard Darst, March 2011

import saiga12

S0 = saiga12.Sys(N=10) ; S0.N = 10
S1 = saiga12.Sys(N=20) ; S1.N = 20
S2 = saiga12.Sys(N=30) ; S2.N = 30
ARRAY = saiga12.SimData_p * 3
array = ARRAY()
array[0] = S0.SD_p ; array[1] = S1.SD_p ; array[2] = S2.SD_p
S0.C.ctest_array(3, array)
