# Richard Darst, May 2008

import numpy

from saiga12.geom.grid import Grid2d

a, b = 4, 4
S = Grid2d()
S.makegrid(a, b)

S.addParticle(pos=4, type_=3)
S.addParticle(pos=5, type_=3)
S.addParticle(pos=6, type_=3)
S.addParticle(pos=9, type_=3)
#S.addParticleRandomDensity(.7, type_=3)
#S.setCycleMoves()

assert S.energy() == 0
S.eddEnable()
S.lattsite.shape = S.lattShape
#print S.MLL
#print S.MLLr
#print S.lattsite
#print S.conn

print
for i in range(10):
    S.cycle(10)
    #print S.MLL
    #print S.MLLr
    #print S.lattsite
    S.eddConsistencyCheck()
    print i, S.mctime, S.naccept, numpy.sum(S.lattsite != -1 ), S.energy(), S.hash()

S1 = Grid2d()
S1.makegrid(25, 25)
S1.addParticleRandomDensity(.5, type_=3)
S1.setCycleMoves()

S2 = Grid2d()
S2.makegrid(25, 25)
S2.addParticleRandomDensity(.5, type_=3)
S2.setCycleMoves()
S2.eddEnable()

cycleTime = 100
if globals().has_key('short'): cycleTime = 20
for i in xrange(15):
    S1.cycle(cycleTime)
    S2.cycle(cycleTime)
    print i, S1.mctime, S1.naccept, S2.mctime, S2.naccept
    if i > 10:
        assert abs(S1.naccept - S2.naccept) / float(S1.naccept) < .01, \
               "BM Event-driven dyn test, naccept is getting too misaligned"
assert S2.eddConsistencyCheck() == 0
