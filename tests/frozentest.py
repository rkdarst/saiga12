# Richard Darst, January 2010

import saiga12

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

for cycleMode in ('montecarlo', 'ctcc'):
    print "Mode:", cycleMode
    S = saiga12.Grid2d()
    length = 15
    S.makegrid(length, length)
    S.setCycleMode(cycleMode)
    if cycleMode == 'montecarlo':
        S.addParticles({2:.5})
    elif cycleMode == 'ctcc':
        S.addParticles({1:.3})
    S.setCycleMoves()

    print S.lattsite

    frozenSites = frozen_bubble(S, length, radius=5)

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


    # event driven vs regular compriasons.
    origContents = S.lattsite[frozenSites]
    S.setFrozenSites(frozenSites)
    nacceptRegular = S.cycle(10000)
    S.eddEnable()
    print S.MLLlen
    nacceptEDD = 0
    for i in range(10):
        nacceptEDD += S.cycle(1000)
        assert not S.eddConsistencyCheck()
    print nacceptRegular, nacceptEDD
    assert (origContents == S.lattsite[frozenSites]).all()
    assert abs(nacceptRegular-nacceptEDD)/(.5*(nacceptRegular+nacceptEDD)) < .05
