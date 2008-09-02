# Richard Darst 2008

import math, numpy
import time

import saiga12
import saiga12.viz
from saiga12.geom.grid import GridHex2d, Grid3dHCP

from saiga12.geom.util import findAsymmetry, drawConn, checkConnDistances

try:
    from rkddp.interact import interact
except ImportError:
    from code import interact

a = b = c = 6
density=.9999
type_ = 12


# 
# First, test the 3d grid (hexagonal close-packed)
#
S = Grid3dHCP()
S.makegrid(a, b, c)
S.addParticleRandomDensity(density, type_=type_)

# set the first two particles on the x-axis to be a different color:
S.atomtype[S.lattsite[0]] = 3
S.atomtype[S.lattsite[36]] = 3

#print S.coords(numpy.arange(64))

# Print out distance of every atom from every other atom it is
# connected to.
# All of these distances should be 1:
checkConnDistances(S, 1, setType=2, doAssert=True)

# See if any connections are not symmetric:
findAsymmetry(S, doAssert=True)

if not globals().has_key("noviz"):
    V = saiga12.viz.VizSystem(S)
    V.vizMakeBox()
    V.vizDisplay()

    # draw connections:
    drawConn(V, S, pos=1, )
    #drawConn(V, S, 1, whichi=(8, 10))

if __name__ == "__main__":
    interact(local=locals(), banner="")



# 
# Now test the 2d grids:
#
S = GridHex2d()
S.makegrid(a, b)
S.addParticleRandomDensity(density, type_=type_)

# set the first two particles on the x-axis to be a different color:
S.atomtype[S.lattsite[0]] = 3
S.atomtype[S.lattsite[6]] = 3

checkConnDistances(S, 1, setType=2, doAssert=True)

findAsymmetry(S, doAssert=True)

if not globals().has_key("noviz"):
    V = saiga12.viz.VizSystem(S)
    V.vizMakeBox()
    V.vizDisplay()

    # draw connections:
    drawConn(V, S, pos=1, )
    #drawConn(V, S, 1, whichi=(8, 10))

if __name__ == "__main__":
    interact(local=locals(), banner="")
