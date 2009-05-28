# Richard Darst, May 2009

import numpy
import time

import saiga12
from saiga12.geom.grid import Grid2d, Grid3d

expectedChemPotential = None
if False:
    S = Grid2d()
    S.makegrid(6, 6)
    expectedChemPotential = 1.812
    particles = {1:.5}
if True:
    S = Grid3d()
    S.makegrid(15, 15, 15)
    particles = {1:.5}
    expectedChemPotential = 1.8102

S.setCycleMode('ctcc')
#S.uVTchempotential = .5

S.addParticle(4, type_=1)

#
# Test 
#
print "Detecting bad (non-directional) moves:"
for i in xrange(10000):
    pos = S.atompos[0]
    old = pos, S.conn[pos,S.orient[pos]]
    S.cycle()
    pos = S.atompos[0]
    if pos != old[0]:
        assert old[0] == S.conn[pos,S.orient[pos]], "Non-allowed transition"
        assert old[1] == pos, "Non-allowed transition"

S.hardness = 1.
S.addParticles(particles)
S.anneal()
S.setCycleMoves()
averager = saiga12.util.Averager()
print "running with many particles"
numOrient = sum(S.orient != -1)
for i in xrange(100):
    assert S.energy() == 0., "Disallowed overlap found"
    S.cycle()
    mu = S.chempotential(1, store=False)
    if mu != float('inf'): averager.add(mu)
    #print averager.mean, averager.var
    assert sum(S.orient != -1) == numOrient

print averager.mean
if expectedChemPotential:
    assert abs(averager.mean - expectedChemPotential) < .05, \
           "%s vs %s"%(averager.mean, expectedChemPotential)

# Now test for insertion:
averager = saiga12.util.Averager()
S.setInsertType({1: (1, expectedChemPotential)})
S.setCycleMoves('grandcanonical')
print "running with many particles"
for i in xrange(100):
    assert S.energy() == 0., "Disallowed overlap found"
    S.cycle()
    averager.add(S.chempotential(1, store=False))
    #print averager.mean, averager.var

S.coords()

print '='*15
import pickle
x = pickle.dumps(S)
S2 = pickle.loads(x)
print S.hash(), S2.hash()
assert S.hash() == S2.hash()
assert tuple(S.orient) == tuple(S2.orient)


#from fitz.interact import interact ; interact()
