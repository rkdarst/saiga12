# Richard Darst, 2008

import math
import numpy

from saiga12.geom.grid import Grid3d
import saiga12.io
import saiga12.corrfunc

S = Grid3d()
S.makegrid(15, 15, 15)
S.addParticles({1:.25})
S.setCycleMoves('grandcanonical')
S.setInsertType({1:(1, 5)})
print S.N
for i in range(5):
    S.cycle(100)
    print S.N
S.setCycleMoves('canonical')
print S.N
type_ = 1

def assertwithin(a, b, ratio=.01):
    assert abs(   (a-b) / ((a+b)/2)   )  < ratio

SkList = saiga12.corrfunc.StructCorrList(
    S=S, kmag2s=range(1, 50), type_=type_)

FsList = saiga12.corrfunc.StructCorrList(
    S=S, kmag2s=range(1, 50), type_=type_)


for i in range(25):
    Sold = S.copy()
    S.cycle(100)
    print i
    
    SkList.calcSk(S)
    FsList.calcFs(S, Sold)



v_kmag = SkList.kmags()
v_ssf = SkList.SkAverages()
v_bykvec = SkList.SkArraysByKvec()
# test that the different ways of calculating are within tolerances of
# each other
for a, c in zip(v_ssf, v_bykvec):
    c = numpy.average(c)
    assertwithin(a, c)
if not globals().has_key('noviz'):
    from rpy import r
    SkList.plotSk()
    SkList.plotVertical()
    SkList.plotSkByKvec()
    import code ; code.interact(local=locals(), banner="" )


v_kmag = FsList.kmags()
v_ssf = FsList.SkAverages()
v_byatom = FsList.SkArraysByAtom()
v_bykvec = FsList.SkArraysByKvec()
# test that the different ways of calculating are within tolerances of
# each other
for a, b, c in zip(v_ssf, v_byatom, v_bykvec):
    b = numpy.average(b)
    c = numpy.average(c)
    assertwithin(a, b)
    assertwithin(b, c)
    assertwithin(a, c)

if not globals().has_key('noviz'):
    from rpy import r
    FsList.plotSk()
    FsList.plotVertical()
    FsList.plotSkByKvec()
    FsList.plotSkByAtom()
    import code ; code.interact(local=locals(), banner="")

