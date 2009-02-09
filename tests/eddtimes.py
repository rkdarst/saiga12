# Richard Darst, Jan 2009

import numpy
import sys
import saiga12

S = saiga12.Grid3d()
S.makegrid(15, 15, 15)
if False:  # BM dynamics
    S.hardness = 1
    S.addParticles({3:.5})
    S.anneal()
    cycleMovesList = (1000, 100, 10, 1)
if False:
    S.setCycleMode('kobandersen')
    S.addParticles({3:.7})
    cycleMovesList = (10000, 1000, 100, 10, 1)
if True:
    S.setCycleMode('fredricksonandersen')
    S.inserttype = 1
    S.addParticles({1:.1})
    S.beta = 1./.1
    cycleMovesList = (10000, 1000, 100, 10, 1)
    S.eddEnable()
    S.cycle(1000000000)

    

print "N:", S.N
S.cycle(10000)
totalMoves = 100000
S.setCycleMoves(100)


print "regular moves"
for cycleMoves in cycleMovesList:
    Ratios = [ ]
    print "running: cycleMoves = %6d"%cycleMoves,
    sys.stdout.flush()
    for i in range(totalMoves // cycleMoves):
        S.resetTime() ; S.naccept = 0
        S.cycle(cycleMoves)
        Ratios.append(S.naccept / float(cycleMoves))
    print "    move ratio:", numpy.mean(Ratios), \
          numpy.std(Ratios) / numpy.sqrt(len(Ratios))

print "EDD, initilizing every start"
for cycleMoves in cycleMovesList:
    Ratios = [ ]
    print "running: cycleMoves = %6d"%cycleMoves,
    sys.stdout.flush()
    for i in range(totalMoves // cycleMoves):
        S.resetTime() ; S.naccept = 0
        S.eddEnable()
        S.cycle(cycleMoves)
        Ratios.append(S.naccept / float(cycleMoves))
        #S.eddConsistencyCheck()
        S.eddDisable()
    print "    move ratio:", numpy.mean(Ratios), \
          numpy.std(Ratios) / numpy.sqrt(len(Ratios))

print "EDD, initilizing only once"
S.eddEnable()
for cycleMoves in cycleMovesList:
    Ratios = [ ]
    print "running: cycleMoves = %6d"%cycleMoves,
    sys.stdout.flush()
    for i in range(totalMoves // cycleMoves):
        S.resetTime() ; S.naccept = 0
        S.cycle(cycleMoves)
        Ratios.append(S.naccept / float(cycleMoves))
    print "    move ratio:", numpy.mean(Ratios), \
          numpy.std(Ratios) / numpy.sqrt(len(Ratios))
S.eddDisable()


