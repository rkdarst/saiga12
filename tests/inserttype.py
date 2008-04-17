# Richard Darst, 2008

import saiga12
from saiga12.geom.grid import Grid2d

S = Grid2d(N=100)
S.makegrid(10, 10)
S.uVTchempotential = .5


S.addParticleRandomDensity(.25, 3)
S.addParticleRandomDensity(.5, 5)
print S.lattsite
print "should be 2:", float(sum(S.lattsite.flat)) / len(S.lattsite.flat)



# we will insert a 50/50 mixture
S.setInsertType( ((3, .2), (5, .8)))
print S.inserttypes_prob, S.inserttypes_type


results = [ ]
for i in range(10):
    results.append(S.getInsertType())
print results


print
S.setCycleMoves(0, insertdel=S.N)
for i in range(10000):
    S.cycle()

print S.lattsite
print "number of 3s:", len(S.lattsite[S.lattsite==3])
print "number of 5s:", len(S.lattsite[S.lattsite==5])
print "number empty:", len(S.lattsite[S.lattsite==0].flat)
#print "should be 2:", float(sum(S.lattsite.flat)) / len(S.lattsite.flat)


