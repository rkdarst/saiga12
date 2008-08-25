# Richard Darst 2008

import math, numpy
import time

import saiga12
from saiga12.geom.grid import GridHex2d, GridHex3d

from saiga12.geom.util import findAsymmetry, drawConn, checkConnDistances

a = b = c = 4
density=.999
type_ = 12

def getS(density=density):
    S = GridHex3d()
    S.makegrid(a, b, c)
    S.hardness = 1
    S.addParticleRandomDensity(density, type_=type_)
    return S

S = getS(density=density)

# set the first two particles on the x-axis to be a different color:
S.atomtype[S.lattsite[0]] = 3
#S.atomtype[S.lattsite[36]] = 3
S.atomtype[S.lattsite[16]] = 3

print S.coords(numpy.arange(64))

# Print out distance of every atom from every other atom it is
# connected to.
# All of these distances should be 1:
checkConnDistances(S, 1, setType=2)

# See if any connections are not symmetric:
findAsymmetry(S)

import saiga12.viz
V = saiga12.viz.VizSystem(S)
V.vizMakeBox()
V.vizDisplay()

#drawConn(V, S, 1, whichi=(8, 10))
drawConn(V, S, pos=1, )

from rkddp.interact import interact ; interact()
