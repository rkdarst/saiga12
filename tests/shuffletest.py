# Richard Darst, May 2011

import random

import saiga12
from saiga12.geom.grid import Grid3d


S = Grid3d()
S.makegrid(15, 15, 15)
S.hardness = 1
S.addParticles({1:.01, 2:.3, 3:.2})
S.anneal()
S.hardness = saiga12.inf
S.setCycleMoves()


# make a list of sites that are in a 5x5x5 cube
L = 5
o = 3
cavitysites = [ ]
for i in range(15**3):
    if (i)%15<L and (i)%15**2<15*L and (i)%15**3<15**2*L:
        cavitysites.append(i)

assert len(cavitysites) == len(set(cavitysites))

for i in range(100):
    saiga12.util.shuffleSites(S, cavitysites)
    S = S.copy()

    cavitysites = random.sample(range(S.N), random.randint(100, 1000))
