# Richard Darst, January 2010

import saiga12

S = saiga12.Grid2d()
length = 15
S.makegrid(length, length)
S.addParticles({2:.5})


print S.lattsite
def frozen_bubble(S, length, radius):
    center = length/2., length/2.
    frozen = [ ]
    dists = S.distanceToPoint(range(S.lattSize), center)
    allPos = range(S.lattSize)
    for pos, coord, dist in zip(allPos, S.coords(), dists):
        #print pos, coord, dist
        if dist > radius:
            frozen.append(pos)
    return frozen

frozenSites = frozen_bubble(S, length, radius=4)

# Test freezing
origContents = S.lattsite[frozenSites]
S.setFrozenSites(frozenSites)
# original atom numbers at these sites
S.cycle(1000)
assert (origContents == S.lattsite[frozenSites]).all()

# Test unfreezing

S.setFrozenSites(None)
S.cycle(1000)
assert not (origContents == S.lattsite[frozenSites]).all()

