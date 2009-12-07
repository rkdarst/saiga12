# Richard Darst, December 2009

import numpy

import saiga12
from saiga12.geom.grid import Grid2d

S = Grid2d()
S.makegrid(10,10)

S.addParticleRandomDensity(.1, 1)
S.addParticleRandomDensity(.25, 2)
S.cycle(1000)

from saiga12.corrfunc import SpinGlass, SpinGlassList
S0 = S.copy()

# Assertion test for all particles
SG = SpinGlass(S0)
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
#print SG._siteCorrelation
assert numpy.sum(SG._siteCorrelation) == S.N*S.N*3

# Test diagonal:
SG = SpinGlass(S0)
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
#print SG._siteCorrelation
assert numpy.sum(SG._siteCorrelation.diagonal()) == S.N*3


# Assertion test for mixed-type particles
types = 1, 2
SG = SpinGlass(S0, types=types)
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
S.cycle(10) ; SG += S, S
#print SG._siteCorrelation
assert numpy.sum(SG._siteCorrelation) == S.ntype[types[0]]*S.ntype[types[1]]*3


SGList = SpinGlassList(S)
S.cycle(100)
SGList += S, S
print SGList.resultMatrix()
print SGList.resultList()
print SGList.typePairs()
