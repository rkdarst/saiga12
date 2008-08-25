# Richard Darst, August 2008

import numpy
import os
import resource

from saiga12.geom.grid import Grid3d


S1 = Grid3d()
S1.makegrid(15, 15, 15)
S1.addParticleRandomDensity(.5, type_=3)
S1.setCycleMoves()

#S2 = Grid3d()
#S2.makegrid(25, 25)
#S2.addParticleRandomDensity(.7, type_=3)
#S2.setCycleMoves()
#S2.eddEnable()

t1 = resource.getrusage(resource.RUSAGE_SELF).ru_utime

for i in xrange(20):
    S1.cycle(100)
    #S2.cycle(100)
    print i, S1.mctime, S1.naccept
    #print i, S1.mctime, S1.naccept, S2.mctime, S2.naccept

print "time:", resource.getrusage(resource.RUSAGE_SELF).ru_utime-t1

