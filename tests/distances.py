# Richard Darst, April 2008

import math
import numpy

import saiga12
from saiga12.geom.grid import Grid2d

a = 100
b = 50
S = Grid2d()
S.makegrid(a, b)

S.hardness = 1
#S.addParticleRandomDensity(.75, 3)
S.addParticle(10, type_=1)
S.addParticle(550, type_=1)
S.addParticle(1075, type_=1)
S.anneal()

#print S.lattsite
#print S.atompos
S.setCycleMoves(1)

iterations = 1000
if globals().has_key('short'): iterations = 10
for i in range(iterations):
    #print
    #print
    #print
    
    startpos = S.getPos()
    S.cycle()
    #S.consistencyCheck()
    endpos = S.getPos()

    d = S.distance(endpos, startpos)
    #print d
    x = sum(d) / S.N
    #print S.atompos[0], x
    print S.mctime, x
    #if x != 1.0:
    #    raise
    #print endpos
    
#print S.lattsite
#print S.atompos



S = saiga12.Grid2d()
S.makegrid(5, 5)
sqrt2 = math.sqrt(2)

assert S.distance(6, 12) == sqrt2
assert (S.distance((6, 11), (12, 17)) == sqrt2).all()
assert (S.distance((6, 11), 12) == (sqrt2, 1)).all()

assert S.distanceToPoint(6, (2,2)) == sqrt2
assert (S.distanceToPoint((6, 11), ((2,2), (3,2))) == sqrt2).all()
assert (S.distanceCoords(numpy.asarray(((1.,1), (2,1))), ((2,2), (3,2))) ==
                                                        (sqrt2, sqrt2)).all()

