# Richard Darst, 2008

import math
import time

import saiga12
from saiga12.geom.grid import Grid2d

a = b = 100
beta = 1.
type_ = 3
density = .5
#density = .1
#uVTchempotential = 1.0
uVTchempotential = .3290

udpairs = (#(-6.8114, .001),
           (-5.2735, .005),
           (-4.5852, .01),
           (-2.9424, .05),
           (-2.1957, .1),
           (-1.3778, .2),
           (-1.0791, .25),
           ( -.8080, .3),
           ( -.2813, .4),
           (  .3222, .5),
           ( 1.1914, .6),
           ( 2.8902, .7),)

#class Sc(saiga12.Sys, Grid2d):
#    pass
def getS(density=density):
    S = Grid2d(N=a*b)
    S.beta = beta
    S.makegrid(a, b)
    S.hardness = 1
    S.addParticleRandomDensity(density, type_=type_)
    S.anneal(False)
    S.hardness = saiga12.inf
    S.inserttype = type_
    return S

#for density in (.001, .005, .01, .05, .1, .15, .2, .25, .3, .4, .5, .6, .7):
#for density in (density, ):
#for density in ( ):
if False:
    S = getS()
    #print "E:", S.energy(), S.N
    #print "chempot:", uVTchempotential, "density:", density
    
    S.setMoveProb(S.N, 0)
    for i in range(1000):
        S.cycle(1)
        S.storeAvg("density", S.density)
        mu = S.chempotential()
    print i, S.energy(), S.N, \
          "  %.4f %.4f %.4f %.4f "%(S.density,
                                    S.avg("density"),
                                    S.stddev("density"),
                                    S.secondmoment("density")), \
          " %1.4f %1.4f %1.4f "%(mu, S.avg("chempotential"),
                                 S.stddev("chempotential"))
    S.resetAvgs()

#for uVTchempotential in (uVTchempotential, ):
for i in range(len(udpairs)):
#if False:
    uVTchempotential = udpairs[i][0]
    density = udpairs[i][1]
    S = getS(density=density)
    # now do grand canonical (!)
    S.uVTchempotential = uVTchempotential
    #S.setMoveProb(S.N, 10)
    S.setCycleMoves(shift=0, insertdel=S.N)
    for i in range(1000):
        if i in (5000, ):
            S.resetAvgs()
        S.cycle(1)
        S.avgStore("density", S.density)
        # warning: calling this function automatically stores it in an
        # averaging function.  If you do the chempotential for
        # different types, you want to have it not storet it.  See the
        # method code.
        mu = S.chempotential(inserttype=type_)
        #print S.uVTchempotential
        if False:
            print i, S.energy(), S.N, \
                  "  %.4f %.4f %.4f %.4f "%(S.density, S.avg("density"),
                                            S.avgStddev("density"),
                                            S.avgSecondmoment("density")), \
                  " %1.4f %1.4f "%(mu, S.avg("chempotential"))
            S.consistencyCheck()
    print i, S.energy(), S.N, \
          "  %.4f %.4f %.4f %.4f "%(S.density, S.avg("density"),
                                    S.avgStddev("density"),
                                    S.avgSecondmoment("density")), \
              " %1.4f %1.4f %1.4f"%(mu, S.avg("chempotential"),
                                    S.avgStddev("chempotential"))

