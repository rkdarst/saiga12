
import cPickle as pickle

import saiga12
from saiga12.geom.grid import Grid3d

a = b = c = 15


S = Grid3d()

S.makegrid(a, b, c)
fracA = .5
initialdensity = .25
S.hardness = 1.
S.addParticleRandomDensity(initialdensity * fracA , 1)  # A
S.addParticleRandomDensity(initialdensity * 1.    , 3)  # B
S.anneal()
S.cycle(100)

import saiga12.io
#S.ioSaveVersion = 1


print S.hash(), hash((tuple(S.nneighbors)))
print S.ntype
s = pickle.dumps(S, -1)
print "save size:", len(s)

S2 = pickle.loads(s)
print S2.hash(), hash((tuple(S2.nneighbors)))
print S2.ntype
S2.consistencyCheck()

