
import saiga12
from saiga12.geom.grid2d import Grid2d

class Sc(saiga12.Sys, Grid2d):
    pass

a, b = 25, 25
S = Sc(N=a*b)
S.makeconn_2Dgrid(a, b)

density = .5
S.hardness = 1
S.addParticleRandomDensity(density, type_=3)
S.anneal()
S.hardness = saiga12.inf

S.inserttype = 3

print "E:", S.energy(), S.N
S.setMoveProb(S.N, 0)

for i in range(1000):
    S.cycle(1)
    mu = S.chempotentialEx()
    print S.energy(), S.N, "%1.3f"%mu


S.chempotentialEx = .75
# now do grand canonical (!)
S.setMoveProb(S.N, 10)
for i in range(100000):
    if i in (1000, ):
        S.resetAvgs()
    S.cycle(1)
    mu = S.chempotentialEx()
    print i, S.energy(), S.N, mu, S.avg("chempotentialEx")


