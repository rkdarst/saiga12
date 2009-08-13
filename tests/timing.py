# Richard Darst, August 2008

import numpy
import os
import resource
import timeit

import saiga12
from saiga12.geom.grid import Grid1d, Grid2d, Grid3d, GridHex2d, Grid3dHCP

testRuns = [
    {'name':"overhead",
     'n': 500,
     'setup':'from saiga12 import Grid3d \nGrid3d().makegrid(15,15,15)',
     'command':'Grid3d().makegrid(15,15,15)',},
    {'name':"BM 1d",
     'grid':'Grid1d', 'size':'150', 'particles':{2:.45}, 'n':10000,},
    {'name':"BM 2d",
     'grid':'Grid2d', 'size':'15,15', 'particles':{4:.45}, 'n':10000,},
    {'name':"BM 3d",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{6:.45}, 'n':500,},
    {'name':"BM Hex2d",
     'grid':'GridHex2d', 'size':'16,16', 'particles':{6:.45}, 'n':10000,},
    {'name':"BM 3dHCP",
     'grid':'Grid3dHCP', 'size':'16,16,16', 'particles':{12:.45}, 'n':200,},
    #{'name':"BM 3d",
    # 'grid':'Grid3d', 'size':'15,15,15', 'particles':{3:.45}, 'n':1000,},
    # Grand Canonical stuff
    {'name':"BM 3d GC",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{3:.45}, 'n':500,
     'cycleMoves':{'mode':'grandCanonical'}},
    # CTCC Dynamics
    {'name':"CTCC 3d",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{1:.6}, 'n':500,
     'cycleMode':'ctcc', 'hardness':'anneal'},
    {'name':"CTCC 3d classic",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{1:.6}, 'n':500,
     'cycleMode':saiga12.S12_CYCLE_CTCCclassic, 'hardness':'anneal'},
    {'name':"CTCC 3d GC",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{1:.6}, 'n':500,
     'cycleMode':'ctcc', 'hardness':'anneal',
     'cycleMoves':{'mode':'grandCanonical'}},
    {'name':"CTCC coords",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{1:.6}, 'n':100,
     'cycleMode':'ctcc', 'hardness':'anneal',
     'command':'S.coords()'},
    ]
ToDo = None
#ToDo = set(('CTCC coords', ))
if __name__ == "__main__":
    import sys
    ToDo = set(sys.argv[1:])

def run_test(kwargs):
    stmt = [ ]
    a = stmt.append
    a('from saiga12.geom.grid import %s'%kwargs['grid'])
    a("S = %s()"%kwargs['grid'])
    a("S.makegrid(%s)"%kwargs['size'])
    if 'cycleMode' in kwargs:
        a('S.setCycleMode(%s)'%repr(kwargs['cycleMode']))
    # deal with hardness stuff
    hardness = kwargs.get('hardness', None)
    if hardness:
        if hardness == 'anneal':
            a('S.hardness = 1')
        else:
            a('S.hardness = %s'%hardness)
    # add particles and anneal if needed
    a("S.addParticles(%s)"%repr(kwargs['particles']))
    if hardness == 'anneal':
        a('S.anneal(verbose=False)')
    # set cycle moves
    cycleMoves = kwargs.get('cycleMoves', {})
    a('S.setCycleMoves(**%s)'%repr(cycleMoves))
    return '\n'.join(stmt)

for run in testRuns:
    if ToDo and run['name'] not in ToDo: continue
    # Setup part: either custom or default cycle testing
    if 'setup' in run:
        setup = run['setup']
    else:
        setup = run_test(run)
    nMoves = run['n']
    if 'fast' in globals(): nMoves /= 10
    # Run custom command or default cycle testing?
    if 'command' in run:
        command = "for i in xrange(%s): %s"%(nMoves, run['command'])
    else:
        command = 'S.cycle(%s)'%nMoves

    # Actually do it
    T = timeit.Timer(command, setup)
    number = 3
    try:
        times = T.repeat(3, number=number)
    except:
        T.print_exc()
        import sys ; sys.exit()
    time = min(times)
    print ("%-15s %9.2f / real s  (%1.1e) (%3.1fs here)"%
           (run['name'], run['n']/time, time/run['n'], sum(times)))


"""
Latest:
BM 1d            60078.09 / real s  (1.7e-05) (0.5s here)
BM 2d            27763.21 / real s  (3.6e-05) (1.1s here)
BM 3d             1195.59 / real s  (8.4e-04) (1.3s here)
BM Hex2d         17720.04 / real s  (5.6e-05) (1.7s here)
BM 3dHCP           641.03 / real s  (1.6e-03) (0.9s here)
BM 3d GC           959.46 / real s  (1.0e-03) (1.6s here)
CTCC 3d           1441.41 / real s  (6.9e-04) (1.0s here)
CTCC 3d classic   1175.62 / real s  (8.5e-04) (1.3s here)
CTCC 3d GC        1442.59 / real s  (6.9e-04) (1.0s here)
CTCC coords        263.93 / real s  (3.8e-03) (1.1s here)
"""
