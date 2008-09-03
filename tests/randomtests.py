
import saiga12
import saiga12.geom.grid
import cPickle as pickle

import random

rand = 0
numOptions = 17
dontCycle = False
S_FA = None
for count1 in xrange(1000):
    if count1 > 0:
        rand = random.choice(range(numOptions))

    if rand == 0:
        print count1, "making a regular system"
        #if count1 != 0: continue
        S = saiga12.geom.grid.Grid2d()
        S.makegrid(5, 5)
        S.addParticles( {1: .1, 3:.3})
        S.setCycleMoves(10)

    elif rand == 1:
        print count1, "making a FA system"
        if S_FA is None:
            S_FA = saiga12.geom.grid.Grid2d()
            #_FAprint "."
            S_FA.makegrid(5, 5)
            S_FA.setCycleMode('fredricksonandersen')
            S_FA.addParticles( {2: .5})
            S_FA.inserttype = 2
            S_FA.eddEnable()
        S = S_FA
        #print S.printLattice()
        assert S.eddConsistencyCheck() == 0
    elif rand == 2:
        print count1, "making a KA system"
        #print "."
        S = saiga12.geom.grid.Grid2d()
        #print "."
        S.makegrid(5, 5)
        S.setCycleMode('kobandersen')
        #print S.lattsite
        S.addParticles( {1: .25, 2: .25})

    elif rand == 3:
        #continue
        print count1, "enabling edd"
        S.eddEnable()

    elif rand == 4:
        if S.cycleModeStr != "fredricksonandersen":
            print count1, "disabling edd"
            S.eddDisable()
        else:
            print count1, "not disabling edd (is FA)"
    elif rand == 5:
        print count1, "doing I/O test"
        s = pickle.dumps(S, -1)
        S2 = pickle.loads(s)
        assert S.hash() == S2.hash()
        S = S2

    elif rand == 6:
        execfile("tests/decay.py", {'short':1})
    elif rand == 7:
        execfile("tests/diffusion.py", {'short':1})
    elif rand == 8:
        execfile("tests/distances.py", {'short':1})
    elif rand == 9:
        execfile("tests/fredricksonandersentest.py", {'short':1})
    elif rand == 10:
        execfile("tests/hextest.py", {'noviz':True})
    elif rand == 11:
        execfile("tests/iotest.py", {'short':1})
    elif rand == 12:
        execfile("tests/kobandersentest.py", {'short':1})
    elif rand == 13:
        execfile("tests/eventdrivendynamics.py", {'short':1})
    elif rand == 14:
        execfile("tests/timing.py", {'short':1})


    else:
        if dontCycle:
            continue
        print count1, "cycling 100"
        S.cycle(100)
        #for i in range(100):
        #    S.cycle(1)
        if S._eddEnabled:
            assert S.eddConsistencyCheck() == 0
        S.consistencyCheck()

