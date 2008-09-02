
import pickle

import saiga12
import saiga12.geom.grid
from saiga12.geom.grid import Grid3d

a = b = c = 15

import saiga12.io
#S.ioSaveVersion = 1

def dotest(S):
    firsthash = S.hash(), hash((tuple(S.nneighbors)))
    print firsthash[0], firsthash[1]
    s = pickle.dumps(S, -1)
    print "save size:", len(s)
    
    S2 = pickle.loads(s)
    #S2.eddEnable()
    #S2 = S
    #print S2.hash(), hash((tuple(S2.nneighbors)))
    S2.consistencyCheck()
    S2.eddConsistencyCheck()
    #
    #assert firsthash == (S2.hash(), hash((tuple(S2.nneighbors))))
    #
    saiga12.randomSeed(9999)
    S.cycle(10)
    saiga12.randomSeed(9999)
    S2.cycle(10)
    #printstuff(S, S2)
    #
    #assert S.hash() == S2.hash()
    print S.persist, S2.persist

def printstuff(S, S2):
    for k in ("lattsite", "conn", "connN", "atomtype", "atompos",
        "MLL", "MLL_down", "MLLlen", "MLLlen_down", "MLLr", "MLL"):
        print k, getattr(S.SD, k), getattr(S2.SD, k)

#S = Grid3d()
#S.makegrid(a, b, c)
#fracA = .5
#initialdensity = .25
#S.hardness = 1.
#S.addParticleRandomDensity(initialdensity * fracA , 1)  # A
#S.addParticleRandomDensity(initialdensity * 1.    , 3)  # B
#S.anneal()
#S.cycle(100)
#
#dotest(S)

print "testing fredrickson-andersen dynamics:"
S = saiga12.geom.grid.Grid1d()
S.makegrid(10)
S.setCycleMode("fredricksonandersen")
print S.persist
S.inserttype = 1
S.addParticles({1:.25})
S.eddEnable()
dotest(S)

print "done with iotest.py"
