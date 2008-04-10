# Richard Darst, April 2008
from __future__ import division

import ctypes
import numpy
import os
import sys
import time
import random
RandomSeed = 1361
random.seed(RandomSeed+165)

from common import *

class SimData(ctypes.Structure):
    _fields_ = [
        ("beta", ctypes.c_double),
        ("N", ctypes.c_int),
        ("NMax", ctypes.c_int),
        ("hardness", ctypes.c_double),
        ("chempotentialEx", ctypes.c_double),
        ("inserttype", ctypes.c_int),
        ("lattSize", ctypes.c_int),
        #("partpos", ctypes.c_void_p),
        ("lattsite", ctypes.c_void_p),    # position (occupancy)
        ("conn", ctypes.c_void_p),   # connections
        ("connN", ctypes.c_void_p),  # num of connections for each
        ("connMax", ctypes.c_int),

        ("cumProbAdd", ctypes.c_double),
        ("cumProbDel", ctypes.c_double),
        ]
SimData_p = ctypes.POINTER(SimData)

def getClib():
    filename = "saiga12c.so"
    C = numpy.ctypeslib.load_library(filename,
                                     os.path.dirname(__file__))
    c_int = ctypes.c_int
    c_double = ctypes.c_double

    C.neighbors_pos.restype = c_int
    C.neighbors_pos.argtypes = SimData_p, c_int

    C.energy_pos.restype = c_double
    C.energy_pos.argtypes = SimData_p, c_int

    C.energy_posNeighborhood.restype = c_double
    C.energy_posNeighborhood.argtypes = SimData_p, c_int

    C.energy.restype = c_double
    C.energy.argtypes = SimData_p,

    C.chempotentialEx.restype = c_double
    C.chempotentialEx.argtypes = SimData_p,

    C.cycle.restype = c_int
    C.cycle.argtypes = SimData_p, c_int
    
    C.init_gen_rand.restype = None
    C.init_gen_rand(RandomSeed+641)
    return C


class Sys(object):
    def __init__(self, N):

        SD = SimData()
        self.__dict__["SD"] = SD
        self.C = getClib()
        self.SD_p = ctypes.pointer(SD)

        self.beta = 1
        self.hardness = float("inf")
        self.cumProbAdd = 0
        self.cumProbDel = 0
        self._avgs = {}

        self.N = SD.N = 0

        #self.partpos = numpy.zeros(shape=(self.NMax), dtype=numpy.int_)
        #SD.partpos   = self.partpos.ctypes.data
        #self.partpos[:] = S12_EMPTYSITE

    def initboard(self):
        # particle 1 at site 0
        # particle 2 at site 7
        pass
    def addParticle(self, pos, type_=1):
        """Add particle i at lattice site `pos`"""
        if self.N >= self.NMax:
            print "ERROR: No more room to insert particles"
        if pos < 0 or pos >= self.lattSize:
            print "ERROR: Lattice position is out of bounds"
        if self.lattsite[pos] != S12_EMPTYSITE:
            print "ERROR: lattice site already occupied"
        self.lattsite[pos] = type_
        #self.partpos[self.N] = pos
        self.N += 1
        #self.SD = self.N
    def addParticleRandom(self, n, type_=1):
        """Randomly add n particles to the system.

        If we run out of spots which can produce a non-infinite
        energy, then we return not having added the right number.

        In order to add more particles, try setting hardness
        (invhardness) to make it not hard, add the particles (since it
        isn't hard anymore, the energy won't become infinite), and
        then use the `anneal` method."""
        # make a list of all unoccupied lattice sites.
        spots = \
              list(numpy.arange(self.lattSize)[self.lattsite == S12_EMPTYSITE])
        #print len(spots)
        inserted = 0  # number of particles successfully inserted.
        # repeat until we have inserted the requested number.
        while inserted < n :
            # if no more spots will hold a particle, quit.
            if len(spots) == 0:
                raise Exception("run out of spots to insert particles")
                break
            # pick a random lattice site.  Try inserting one.  If
            # energy becomes infinite, abort that insertion and
            # repeat.
            i = random.randrange(len(spots))
            pos = spots.pop(i)
            #print "inserting at site:", pos
            self.lattsite[pos] = type_
            if self.energy_posNeighborhood(pos) == float("inf"):
                self.lattsite[pos] = S12_EMPTYSITE
                continue
            inserted += 1
            self.N += 1
    def addParticleRandomDensity(self, density, type_=1):
        """Add particles to give system requested density.

        This method produces a total system density equal to the
        target, if there are already particles in the system, we only
        add the number of particles required to get up to that
        density.

        It's just a wrapper around addParticleRandom.
        """
        target = int(round(self.lattSize * density))
        extra = target - self.N
        if extra <= 0:
            print "We already have too many particles in the system."
            return
        self.addParticleRandom(extra, type_=type_)
    def cycle(self, n=1):
        """Run MC trial moves for once cycle.

        Cycles defined by setMoveProb()
        """
        moves = int(n * self.movesPerCycle)
        self.C.cycle(self.SD_p, moves)
    def setMoveProb(self, shift=None, insertdel=0):
        """Sets moves executed each cycle.

        To run GCE, do:
        - set this method
        - set self.inserttype
        - set self.chempotentialEx
        """
        if shift is None:
            shift = self.N
        self.movesPerCycle = shift + insertdel

        self.cumProbAdd = (insertdel/2.) / self.movesPerCycle
        self.cumProbDel = (insertdel) / self.movesPerCycle
        print self.cumProbAdd, self.cumProbDel
        
        
        
    def removeParticle(self, pos):
        """not inplemented yet"""
        pass
    def neighbors_pos(self, i):
        """Number of neighbors of a lattice site"""
        return self.C.neighbors_pos(self.SD_p, 0)
    def energy(self):
        """Total energy of the system."""
        return self.C.energy(self.SD_p)
    def energy_pos(self, pos):
        """Energy due to a single lattice site."""
        return self.C.energy_pos(self.SD_p, pos)
    def energy_posNeighborhood(self, pos):
        """Energy due to a single position and all neighboring sites"""
        return self.C.energy_posNeighborhood(self.SD_p, pos)
    def findInfiniteEnergy(self):
        """Utility function to find the lattice sites which have inf energy.
        """
        for pos in range(self.lattSize):
            if self.energy_posNeighborhood(pos) == float("inf"):
                print pos, self.energy_pos(pos), \
                      self.energy_posNeighborhood(pos)
                self.printLatticeLocal(pos, width=2)
    def hash(self):
        """Utility function to print checksum of state."""
        x = ( tuple(self.lattsite.flat),
              tuple(self.conn.flat), tuple(self.connN.flat),
              self.hardness, self.lattSize, 
              )
        return hash(x)
    
    def anneal(self):
        """Slowly increase the hardness of the system until energy is zero.

        Used to go from an initial configuration which has overlaps,
        to a hard configuration with no overlaps.
        """
        hardness = 1.
        while self.energy() > 0:
            self.hardness = hardness
            self.C.cycle(self.SD_p, self.N * 2)
            print "hardness: %5d  energy: %8.2f"%(hardness, self.energy()), \
                  self.N
            hardness += 1
        self.hardness == float("inf")
    def chempotentialEx(self):
        """Chemical potential of the system, test inserting at every site.
        """
        if self.inserttype == S12_EMPTYSITE:
            raise Exception("Must set self.inserttype to the type of particle to insert")
        
        mu = self.C.chempotentialEx(self.SD_p)
        if mu != inf:
            self.storeAvg("chempotentialEx", mu)
        return mu

    def storeAvg(self, name,  value):
        x = self._avgs.setdefault(name, [0., 0. ])
        x[0] += 1
        x[1] += value
    def avg(self, name):
        x = self._avgs[name]
        return x[1] / x[0]
    def resetAvgs(self):
        self._avgs = { }

    def _density_get(self):
        """Density of the system, N/lattSize"""
        return self.N / self.lattSize
    density = property(fget=_density_get)



    def __getattr__(self, attrname):
        if hasattr(self.SD, attrname):
            return getattr(self.SD, attrname)
        raise AttributeError("unknown ")
    def __setattr__(self, name, value):
        if hasattr(self.SD, name):
            setattr(self.SD, name, value)
        else:
            self.__dict__[name] = value
    

def correlation(lattice0, lattice1, S):
    avg = S.density * S.density * S.lattSize * 3 * 3
    #print avg
    return float(sum(lattice0 * lattice1)) - avg
    

if __name__ == "__main__":
    from saiga12.geom.grid2d import Grid2d
    class Sc(Sys, Grid2d):
        pass

    a, b = 500, 500
    S = Sc(N=a*b)
    print "x"
    S.makeconn_2Dgrid(a, b)
    print "y"
    S.hardness = 10.
    S.addParticleRandomDensity(.75, type_=3)
    
    print "E:", S.energy()
    #sys.exit()
    latticeOrig = S.lattsite.copy()

    print S.energy()
    print

    S.anneal()
    #for i in range(1000):
    #    S.C.cycle(S.SD_p, S.N)
    #    print S.energy(), correlation(latticeOrig, S.lattsite, S), S.N
    #    #S.printLattice()
    #    #time.sleep(.005)
