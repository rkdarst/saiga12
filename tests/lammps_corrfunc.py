# Richard Darst, 2009

import saiga12.logtime
import saiga12.lammps.lmp
#saiga12.lammps.lmp.mpi_autoinit(procs=2)
import saiga12.lammps.liblmp as liblmp
from lammpstest import make_system
S = make_system()
S.makecoords()

# equilibrate
S.cycle(10000)

from saiga12 import corrfunc
type_ = 1


SkList = corrfunc.StructCorrList(
    S=S, kmag2s=range(1, 50), type_=type_, orthogonal=False)
def make_FsList():
    FsList = corrfunc.StructCorrList(
        S=S, kmag2s=range(25, 26), type_=type_, orthogonal=True)
    return FsList

import collections
FsList_t = collections.defaultdict(make_FsList)
f = file('blah', 'w')

for i in range(10):
    S0 = S.copy()
    S0.makecoords()
    cur_t = 0

    FsList_t[cur_t].calcFs(S0, S)
    for next_t in saiga12.logtime.logTimeTimesGenerator():
        S.cycle(next_t-cur_t)
        cur_t = next_t
        print "*"*10, next_t, cur_t, next_t-cur_t
        print S0.realcoords(0), S.realcoords(0)
        print [(x, l.SsfDict[5].Sk()) for x, l in sorted(FsList_t.items())]
        #from fitz import interactnow
        FsList_t[cur_t].calcFs(S0, S)
        print       cur_t, FsList_t[cur_t].SsfDict[5].Sk()
        print >> f, cur_t, FsList_t[cur_t].SsfDict[5].Sk()
        f.flush()
        if FsList_t[cur_t].SsfDict[5].Sk() < .1: # or cur_t > 500:
            break
    for t in sorted(FsList_t):
        print t, FsList_t[t].SsfDict[5].Sk()
    
