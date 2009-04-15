# Richard Darst, April 2009

from saiga12 import Grid3d

def getS():
    S = Grid3d()
    S.beta = 1.
    S.makegrid(16, 16, 16)
    S.hardness = 1
    S.addParticles({1:.15, 3:.35})
    S.setCycleMoves()
    S.anneal(False)
    S.hardness = 1
    return S

S = getS()

S.vib_init()
S.vib_offset()
S.cycle(10)
S.vib_offset()
S.cycle(10)
S.vib_offset()

