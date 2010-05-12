# Richard Darst, August 2008

import numpy
import os
import resource
import sys

from saiga12.geom.grid import Grid2d, Grid3d

S = Grid2d()
S.makegrid(10, 10)
S.inserttype = 1
S.setCycleMode('spiral')

if False:
    S = Grid2d()
    S.makegrid(26, 26)
    S.inserttype = 1
    S.setCycleMode('spiral')
    #for density in (0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .99):
    for density in numpy.arange(0, 1+.02, .02):

        # Test printing of empty and full lattices
        print "== density: %s =="%density
        S.delAll()
        if density > 0:
            S.addParticles({1:density})
        S.eddEnable()
        S.printLattice()
        print
    S.delAll()
    sys.exit()

# Random sites
S.addParticles({1:.5})
S.eddEnable()
S.setCycleMoves(1)
S.printLattice()
print


assert not S.eddConsistencyCheck()
S.cycle(1)
S.printLattice()
print
assert not S.eddConsistencyCheck()
sys.exit()


dim = 30, 30
density = .5
type_ = 1

#S1 = Grid2d()
#S1.makegrid(*dim)
#S1.inserttype = type_
#S1.setCycleMode('fredricksonandersen')
#S1.addParticles({type_: density})
#S1.setCycleMoves(1)
#S1.eddEnable()

S2 = Grid2d()
S2.makegrid(*dim)
S2.inserttype = type_
S2.setCycleMode('spiral')
S2.addParticles({type_: density})
S2.setCycleMoves(1)
S2.eddEnable()

#t1 = resource.getrusage(resource.RUSAGE_SELF).ru_utime

#S2.eddConsistencyCheck()
##import sys ; sys.exit()
#for i in range(10):
#    S2.printLattice()
#    S2.cycle(1)
#    print S2.naccept
#    #print S.MLL
#    #print S.MLLr
#    #print S.lattsite
#    S2.eddConsistencyCheck()
#    print i, S2.mctime, S2.naccept, S2.N, S2.energy(), S2.hash()
#import sys ; sys.exit()


cycleTime = 100
if globals().has_key('short'): cycleTime = 5
S2.resetTime()
#S2.eddFindBestMode()
for i in xrange(10):
    #S1.printLattice()
    #S1.cycle(cycleTime)
    S2.cycle(cycleTime)
    #print i, S1.mctime, S1.naccept

    #print i, S2.mctime, S2.naccept, S2.MLLextraTime
    #assert S2.eddConsistencyCheck() == 0
    #S2.printLattice()
    
    print "-->", i, "  ", \
          S2.mctime, S2.N, S2.naccept, "cc:", S2.eddConsistencyCheck()
          #S1.mctime, S1.N, S1.naccept, "  ", \
    assert S2.naccept > 0, \
               "Fredrickson-Andersen test not running!"
    #raw_input('> ')

#print "time:", resource.getrusage(resource.RUSAGE_SELF).ru_utime-t1

assert S2.eddConsistencyCheck() == 0
