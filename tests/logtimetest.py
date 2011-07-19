from saiga12.geom.grid import Grid2d
from saiga12.logtime import logTimeGenerator

totalTime = 50

cbAlist = [ ]
cbAtime = 10
def callbackA(S, deltaT):
    print "A", S.mctime, deltaT
    cbAlist.append(deltaT)
cbBlist = [ ]
cbBtime = 50
def callbackB(S, deltaT):
    print "B", S.mctime, deltaT
    cbBlist.append(deltaT)


S = Grid2d()
S.makegrid(4, 4)
S.addParticles({4:.1})

for S, deltaT in logTimeGenerator(S, maxTime=totalTime-1,
                                  callbacks=[(cbAtime, callbackA),
                                             (cbBtime, callbackB)]):
    print deltaT

assert len(cbAlist) == len(cbAlist)
for t in cbAlist:
    assert (t//cbAtime) * cbAtime == t
for t in cbBlist:
    assert (t//cbBtime) * cbBtime == t
assert len(cbAlist) == totalTime / cbAtime
assert len(cbBlist) == totalTime / cbBtime

assert S.mctime == totalTime
