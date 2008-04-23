

import saiga12
from saiga12.geom.grid import Grid2d


b = 100
S = Grid2d()
S.makegrid(b, b)

S.hardness = 1
S.addParticleRandomDensity(.75, 3)
S.anneal()

#print S.lattsite
#print S.atompos

def getpos(S):
    return S.atompos[S.atompos != saiga12.S12_EMPTYSITE]

startpos = getpos(S)
#print startpos
print "."

S.setCycleMoves(10)
for i in range(10000):
    endpos = getpos(S)
    d =S.distance(endpos, startpos)
    print sum(d) / S.N
    #print endpos
    S.cycle()
    
#print S.lattsite
#print S.atompos

