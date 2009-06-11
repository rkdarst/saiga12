# Richard Darst, August 2008

import numpy
import os
import resource

from saiga12.geom.grid import Grid2d, Grid3d

def within(a, b, ratio=.01):
    return abs(   (a-b) / ((a+b)/2.)   )  < ratio

dim = 5, 5, #15 #15, 15
particles = {1:.4 }
def makeS():
    S = Grid2d()
    S.makegrid(*dim)
    S.setCycleMode('ctcc')
    S.addParticles(particles)
    S.setCycleMoves()
    return S
S1 = makeS()

S2 = makeS()
S2.setCycleMoves()
S2.eddEnable()

assert not S2.eddConsistencyCheck()

cycleTime = 1000
if globals().has_key('short'): cycleTime = 20
#S2.eddFindBestMode()
for i in xrange(100):
    S1.cycle(cycleTime)
    S2.cycle(cycleTime)

    assert S2.eddConsistencyCheck() == 0
    print "-->", i, "  ", S1.mctime, S1.naccept, "  ", S2.mctime, S2.naccept

    #raw_input('> ')
    if i > 5:
        assert within(S1.naccept, S2.naccept, ratio=.05), \
               "CTCC test, naccept is getting too misaligned"
        assert S1.naccept > 0, \
               "CTCC test not running!"
        assert S2.naccept > 0, \
               "CTCC test not running!"
