# Richard Darst, August 2008

import numpy
import os
import resource

from saiga12.geom.grid import Grid2d, Grid3d

dim = 3, 3
density = .5
type_ = 1

S1 = Grid2d(cycleMode='fredricksonandersen', energyMode='fredricksonandersen')
S1.makegrid(*dim)
S1.addParticleRandomDensity(density, type_=type_)
S1.setCycleMoves(1)
S1.inserttype = type_

S2 = Grid2d(cycleMode='fredricksonandersen', energyMode='fredricksonandersen')
S2.makegrid(*dim)
S2.addParticleRandomDensity(density, type_=type_)
S2.setCycleMoves(1)
S2.eddEnable()
S2.inserttype = type_


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
#    print i, S2.mctime, S2.naccept, numpy.sum(S2.lattsite != -1 ), S2.energy(), S2.hash()
#import sys ; sys.exit()


cycleTime = 10000
#S2.eddFindBestMode()
for i in xrange(10):
    S1.printLattice()
    S1.cycle(cycleTime)
    S2.cycle(cycleTime)
    #print i, S1.mctime, S1.naccept

    #print i, S2.mctime, S2.naccept, S2.MLLextraTime
    #assert S2.eddConsistencyCheck() == 0
    #S2.printLattice()
    
    print "-->", i, "  ", \
          S1.mctime, S1.N, S1.naccept, "  ", \
          S2.mctime, S2.N, S2.naccept, "cc:", S2.eddConsistencyCheck()
    #raw_input('> ')

#print "time:", resource.getrusage(resource.RUSAGE_SELF).ru_utime-t1

S2.eddConsistencyCheck()
