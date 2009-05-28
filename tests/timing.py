# Richard Darst, August 2008

import numpy
import os
import resource
import timeit

from saiga12.geom.grid import Grid1d, Grid2d, Grid3d, GridHex2d, Grid3dHCP

testRuns = [
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
    {'name':"CTCC 3d GC",
     'grid':'Grid3d', 'size':'15,15,15', 'particles':{1:.6}, 'n':500,
     'cycleMode':'ctcc', 'hardness':'anneal',
     'cycleMoves':{'mode':'grandCanonical'}},
    ]

def run_test(kwargs):
    stmt = [ ]
    a = stmt.append
    a('from saiga12.geom.grid import %s'%kwargs['grid'])
    a("S = %s()"%kwargs['grid'])
    a("S.makegrid(%s)"%kwargs['size'])
    if 'cycleMode' in kwargs:
        a('S.setCycleMode("%s")'%kwargs['cycleMode'])

    hardness = kwargs.get('hardness', None)
    if hardness:
        if hardness == 'anneal':
            a('S.hardness = 1')
        else:
            a('S.hardness = %s'%hardness)
    a("S.addParticles(%s)"%repr(kwargs['particles']))
    if hardness == 'anneal':
        a('S.anneal(verbose=False)')

    cycleMoves = kwargs.get('cycleMoves', {})
    a('S.setCycleMoves(**%s)'%repr(cycleMoves))
    return '\n'.join(stmt)

for run in testRuns:
    setup = run_test(run)
    nMoves = run['n']
    if 'fast' in globals(): nMoves /= 10
    T = timeit.Timer('S.cycle(%s)'%nMoves, setup)
    number = 3
    try:
        times = T.repeat(3, number=number)
    except:
        T.print_exc()
        import sys ; sys.exit()
    time = min(times)
    print ("%-15s %9.2f mctime per real second (%3.1fs here)"%
           (run['name'], run['n']/time, sum(times)))


"""
Latest:
BM 1d            59684.75 mctime per real second (0.5s here)
BM 2d            28165.29 mctime per real second (1.1s here)
BM 3d             1203.85 mctime per real second (1.2s here)
BM Hex2d         17865.05 mctime per real second (1.7s here)
BM 3dHCP           647.55 mctime per real second (0.9s here)
BM 3d GC           967.99 mctime per real second (1.6s here)
CTCC 3d           1369.39 mctime per real second (1.1s here)
CTCC 3d GC        1368.52 mctime per real second (1.1s here)
"""
