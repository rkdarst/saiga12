# Richard Darst, May 2009


import saiga12
#import saiga12.geom
from saiga12.geom.grid import *
from saiga12.geom.util import checkConnDistances, checkConnReversibility



ccr = checkConnReversibility
ccd = lambda S: checkConnDistances(S, 1., doAssert=True)


S = Grid1d()    ; S.makegrid(6)     ; ccd(S) ; ccr(S)
S = Grid2d()    ; S.makegrid(6,6)   ; ccd(S) ; ccr(S)
S = Grid3d()    ; S.makegrid(6,6,6) ; ccd(S) ; ccr(S)
S = GridHex2d() ; S.makegrid(6,6)   ; ccd(S) ; ccr(S)
S = Grid3dHCP() ; S.makegrid(6,6,6) ; ccd(S) ; ccr(S)
