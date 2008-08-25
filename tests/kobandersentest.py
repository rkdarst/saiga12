# Richard Darst, August 2008

import numpy
import os
import resource

from saiga12.geom.grid import Grid2d, Grid3d

dim = 25, 25 #15, 15
density = .7
type_ = 2

S1 = Grid2d(cycleMode='kobandersen', energyMode='kobandersen')
S1.makegrid(*dim)
S1.addParticleRandomDensity(density, type_=type_)
S1.setCycleMoves()


S2 = Grid2d(cycleMode='kobandersen', energyMode='kobandersen')
S2.makegrid(*dim)
S2.addParticleRandomDensity(density, type_=type_)
S2.setCycleMoves()
S2.eddEnable()

#t1 = resource.getrusage(resource.RUSAGE_SELF).ru_utime

cycleTime = 1000
S2.eddFindBestMode()
for i in xrange(100):
    S1.cycle(cycleTime)
    S2.cycle(cycleTime)
    #print i, S1.mctime, S1.naccept
    #S1.printLattice()

    #print i, S2.mctime, S2.naccept, S2.MLLextraTime
    #assert S2.eddConsistencyCheck() == 0
    #S2.printLattice()
    
    print "-->", i, "  ", S1.mctime, S1.naccept, "  ", S2.mctime, S2.naccept, S2.eddConsistencyCheck(),
    #raw_input('> ')

#print "time:", resource.getrusage(resource.RUSAGE_SELF).ru_utime-t1

S2.eddConsistencyCheck()
