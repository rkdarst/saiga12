# Richard Darst, August 2008

import numpy
import os
import resource

from saiga12.geom.grid import Grid2d, Grid3d

def within(a, b, ratio=.01):
    return abs(   (a-b) / ((a+b)/2.)   )  < ratio

dim = 25, 25 #15, 15
density = .7
type_ = 2

S1 = Grid2d()
S1.setCycleMode('kobandersen')
S1.makegrid(*dim)
S1.addParticleRandomDensity(density, type_=type_)
S1.setCycleMoves()


S2 = Grid2d()
S2.setCycleMode('kobandersen')
S2.makegrid(*dim)
S2.addParticleRandomDensity(density, type_=type_)
S2.setCycleMoves()
S2.eddEnable()

#t1 = resource.getrusage(resource.RUSAGE_SELF).ru_utime

cycleTime = 1000
if globals().has_key('short'): cycleTime = 20
#S2.eddFindBestMode()
for i in xrange(10):
    S1.cycle(cycleTime)
    S2.cycle(cycleTime)
    #print i, S1.mctime, S1.naccept
    #S1.printLattice()

    #print i, S2.mctime, S2.naccept, S2.MLLextraTime
    #assert S2.eddConsistencyCheck() == 0
    #S2.printLattice()
    
    print "-->", i, "  ", S1.mctime, S1.naccept, "  ", S2.mctime, S2.naccept, S2.eddConsistencyCheck()
    #raw_input('> ')
    if i > 5:
        assert within(S1.naccept, S2.naccept, ratio=.05), \
               "Kob-Andersen test, naccept is getting too misaligned"
        assert S1.naccept > 0, \
               "Kob-Andersen test not running!"
        assert S2.naccept > 0, \
               "Kob-Andersen test not running!"

#print "time:", resource.getrusage(resource.RUSAGE_SELF).ru_utime-t1
assert S2.eddConsistencyCheck() == 0


#
# test persistence function stuff
#
S1 = Grid2d()
S1.setCycleMode('kobandersen')
S1.makegrid(3, 3)
S1.addParticleRandom(n=1, type_=3)

print S1.lattsite, S1.persist
S1.cycle()
print S1.lattsite, S1.persist

S1._allocPersistArray()

print S1.lattsite, S1.persist
S1.cycle()
print S1.lattsite, S1.persist

print

import saiga12 ; saiga12.randomSeed(174)
S1 = Grid2d()
S1.setCycleMode('kobandersen')
S1.makegrid(3, 3)
S1.addParticleRandom(n=1, type_=3)
S1.eddEnable()

print S1.lattsite, S1.persist
S1.cycle()
print S1.lattsite, S1.persist

S1._allocPersistArray()

print S1.lattsite, S1.persist
S1.cycle()
print S1.lattsite, S1.persist


#from rkddp.interact import interact ; interact()




#### Test Soft KA: ####
S = Grid2d()
S.setCycleMode('kobandersen')
S.makegrid(5, 5)
# Make a grid that has two open rows:
# x  xx
# x  xx
# x  xx
# x  xx
# x  xx
# Everything (with type 2 particles) should be immobile
for pos in range(25):
    if pos%25 in (1,2):  continue
    S.addParticle(pos, type_=2)
S.setCycleMoves()
assert S.eddCheckAllowed()

# Regular dynamics.  Nothing should be able to move:
origContents = S.atompos.copy()
S.cycle(1000)
assert (origContents == S.atompos).all()
assert S.naccept == 0

# Now make it soft, and things should be able to move:
S.hardness = 1
S.beta = 1/.5
S.flags |= saiga12.S12_FLAG_KA_SOFT
assert not S.eddCheckAllowed()

S.cycle(1000)
assert not (origContents == S.atompos).all()
assert S.naccept != 0


#### Test grand canonical: ###
# This test is like the one above - it should be completly blocked -
# but we use grand canonical dynamics to allow sampling.
S = Grid2d()
S.setCycleMode('kobandersen')
S.makegrid(5, 5)
# Everything (with type 2 particles) should be immobile
for pos in range(25):
    if pos%25 in (1,2):  continue
    S.addParticle(pos, type_=2)
S.setCycleMoves()

# Regular dynamics.  Nothing should be able to move:
origContents = S.atompos.copy()
S.cycle(1000)
assert (origContents == S.atompos).all()
assert S.naccept == 0

# Now enable grand-canonical dynamics
S.setCycleMoves(shift=.5*S.N, insertdel=.5*S.N)
S.setInsertType({2: (1., .5)})
# And insist that there is some motion once it moves
S.cycle(1000)
assert not (origContents == S.atompos).all()
assert S.naccept != 0



#### Test that type3 particles with the grid this can move ####
S = Grid2d()
S.setCycleMode('kobandersen')
S.makegrid(5, 5)
# Same grid with two open rows, but type 3 particles so they _should_ be modile
for pos in range(25):
    if pos%25 in (1,2):  continue
    S.addParticle(pos, type_=3)
S.setCycleMoves()

# Regular dynamics.  Nothing should be able to move:
origContents = S.atompos.copy()
S.cycle(1000)
assert not (origContents == S.atompos).all()
assert S.naccept != 0


