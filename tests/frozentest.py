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

for cycleMode in ('montecarlo', 'ctcc', 'kobandersen',):
    print
    print "Mode:", cycleMode
    S = saiga12.Grid2d()
    length = 15
    S.makegrid(length, length)
    S.setCycleMode(cycleMode)
    if cycleMode == 'montecarlo':
        S.addParticles({2:.5})
    elif cycleMode == 'ctcc':
        S.addParticles({1:.3})
    elif cycleMode == 'kobandersen':
        S.addParticles({3:.5})
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
    assert not S.eddConsistencyCheck()
    print S.MLLlen
    nacceptEDD = 0
    for i in range(10):
        nacceptEDD += S.cycle(1000)
        assert not S.eddConsistencyCheck()
    print nacceptRegular, nacceptEDD
    assert (origContents == S.lattsite[frozenSites]).all()
    assert abs(nacceptRegular-nacceptEDD)/(.5*(nacceptRegular+nacceptEDD)) < .05


for cycleMode in ('fredricksonandersen', 'east',):
    print
    print "Mode:", cycleMode
    S = saiga12.Grid2d()
    #length = 5
    length = 10
    S.makegrid(length, length)
    S.setCycleMode(cycleMode)
    if cycleMode == 'fredricksonandersen':
        S.addParticles({1:.4})
        S.inserttype = 1
        S.beta = 1/.5
    elif cycleMode == 'east':
        S.addParticles({1:.4})
        S.inserttype = 1
        S.beta = 1/.5
    #S.setCycleMoves()

    #frozenSites = frozen_bubble(S, length, radius=5)
    #frozenSites = [0,1,2,3,4, 5,9, 10,14, 15,19, 20,21,22,23,24]
    frozenSites = [0,1,2,3,4,5,6,7,8,9,
                   10,19, 20,29, 30,39, 40,49, 50,59, 60,69, 70,79, 80,89,
                   90,91,92,93,94,95,96,97,98,99, ]
    S.setFrozenSites(frozenSites)
    origContents = S.lattsite[frozenSites] != saiga12.S12_EMPTYSITE

    S.eddEnable()
    assert not S.eddConsistencyCheck()
    S.printLattice()
    for i in range(100):
        S.cycle(100)
        assert not S.eddConsistencyCheck()
        assert ((S.lattsite[frozenSites] != saiga12.S12_EMPTYSITE)
                == origContents).all()


        #print S.printLattice().reshape(length, length)
        #raw_input('> ')



#### Soft KA model ###
if True:
    print "\nMode: kob-andersen soft"
    S = saiga12.Grid2d()
    length = 15
    S.makegrid(length, length)
    S.setCycleMode('kobandersen')
    S.addParticles({3:.5})
    S.setCycleMoves()
    frozenSites = frozen_bubble(S, length, radius=5)

    S.hardness = 1
    S.beta = 1/.5
    S.flags |= saiga12.S12_FLAG_KA_SOFT

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
