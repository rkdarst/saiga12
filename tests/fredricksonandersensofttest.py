# Richard Darst, August 2008

import numpy

import saiga12
from saiga12.geom.grid import Grid2d, Grid3d
import saiga12.util

dim = 30, 30
temp = .5
density = saiga12.util.TtoC(temp)
type_ = 1

S = Grid2d()
S.makegrid(*dim)
S.inserttype = type_
S.setCycleMode('fredricksonandersen')
S.addParticles({type_: density})
S.setCycleMoves(1)
S.eddEnable()
S.hardness = 1
S.beta = 1/.05
#S.beta = -10.


for i in range(100):
    S.cycle(100)
    assert not S.eddConsistencyCheck()

for i in range(100):
    S.cycle(1)
    assert not S.eddConsistencyCheck()
    S.printLattice()
    print
    #raw_input('> ')


# Start with everything down and
S.delAll()
S.resetTime()
S.naccept = 0
S.beta = 1/1.
S.eddEnable()
S.printLattice()
print
assert not S.eddConsistencyCheck()
S.cycle(100)
assert not S.eddConsistencyCheck()
S.printLattice()
assert S.naccept > 0
