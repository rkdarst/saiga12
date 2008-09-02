
import sys
import time

import saiga12
from saiga12.geom.grid import Grid2d, Grid3d

#class Sc(saiga12.Sys, Grid3d):
#    pass

a = b = c = 10
N = a*b*c
density = .58
type_ = 3

MaxTime = 1000
skip = 5
iterations = 100

S = Grid3d(N=N)
S.beta = 1.
S.makegrid(a, b, c)
S.hardness = 1
S.addParticleRandomDensity(density, type_=type_)
S.anneal()
S.hardness = saiga12.inf
S.inserttype = type_
S.setCycleMoves()   # defaults to translating (canonical), shift=self.N

S.cycle(100) #  equilibrate

corrfunc = { }

for ii in range(iterations):
    print "iteration:", ii
    S.resetTime()
    latticeStart = S.lattsite.copy()
    while S.mctime < MaxTime:
        corr = saiga12.correlation(latticeStart,S.lattsite) - N*density*density
        l = corrfunc.setdefault(S.mctime, [])
        l.append(corr)
    
        S.cycle(skip)
        print "\r",
        print S.mctime,
        sys.stdout.flush()

    # do it every time, so see it build up slowly.
    #logfile = file("logfile.txt", "w")
    #for mctime in sorted(corrfunc.keys()):
    #    print >> logfile, mctime, sum(corrfunc[mctime])/len(corrfunc[mctime])
    #logfile.flush()
    
