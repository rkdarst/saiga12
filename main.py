# Richard Darst, April 2008
from __future__ import division

import ctypes
import math
import numpy
import os
import sys
import time
import random
RandomSeed = 1361
random.seed(RandomSeed+165)

from saiga12.common import *

class SimData(ctypes.Structure):
    _fields_ = [
        ("beta", ctypes.c_double),
        ("N", ctypes.c_int),
        #("NMax", ctypes.c_int),  # lattSize is NMax
        ("hardness", ctypes.c_double),
        ("uVTchempotential", ctypes.c_double),
        ("inserttype", ctypes.c_int),
        ("widominserttype", ctypes.c_int),
        ("inserttypes_n", ctypes.c_int),
        ("inserttypes_prob", ctypes.c_void_p),
        ("inserttypes_type", ctypes.c_void_p),
        ("inserttypes_plookup", ctypes.c_void_p),
        ("inserttypes_mulookup", ctypes.c_void_p),
        ("lattSize", ctypes.c_int),
        #("partpos", ctypes.c_void_p),
        ("lattsite", ctypes.c_void_p),    # position (occupancy)
        ("conn", ctypes.c_void_p),   # connections
        ("connN", ctypes.c_void_p),  # num of connections for each
        ("connMax", ctypes.c_int),
        ("nneighbors", ctypes.c_void_p),
        ("atomtype", ctypes.c_void_p),
        ("atompos", ctypes.c_void_p),
        ("ntype", ctypes.c_void_p),

        ("cumProbAdd", ctypes.c_double),
        ("cumProbDel", ctypes.c_double),
        ]
SimData_p = ctypes.POINTER(SimData)

neighlist = 1

def getClib():
    filename = "saiga12c.so"
    C = numpy.ctypeslib.load_library(filename,
                                     os.path.dirname(__file__))
    c_int = ctypes.c_int
    c_double = ctypes.c_double

    C.neighbors_pos.restype = c_int
    C.neighbors_pos.argtypes = SimData_p, c_int

    C.getInsertType.restype = c_int
    C.getInsertType.argtypes = SimData_p, 
    if neighlist:
        C.addParticle.restype = None
        C.addParticle.argtypes = SimData_p, c_int, c_int

        C.delParticle.restype = None
        C.delParticle.argtypes = SimData_p, c_int

    C.energy_pos.restype = c_double
    C.energy_pos.argtypes = SimData_p, c_int

    C.energy_posNeighborhood.restype = c_double
    C.energy_posNeighborhood.argtypes = SimData_p, c_int

    C.energy.restype = c_double
    C.energy.argtypes = SimData_p,

    C.chempotential.restype = c_double
    C.chempotential.argtypes = SimData_p, c_int

    C.cycle.restype = c_int
    C.cycle.argtypes = SimData_p, c_int
    
    C.init_gen_rand.restype = None
    C.init_gen_rand(RandomSeed+641)
    return C


class Sys(object):
    def __init__(self, N=None):

        SD = SimData()
        self.__dict__["SD"] = SD
        self.C = getClib()
        self.SD_p = ctypes.pointer(SD)

        self.beta = 1.
        self.hardness = float("inf")
        self.cumProbAdd = 0
        self.cumProbDel = 0
        self.inserttype = S12_EMPTYSITE
        self.widominserttype = S12_EMPTYSITE

        self.N = 0
        self.resetTime()
        self.avgReset()
        self.setCycleMoves(shift=1)  # this must be reset once N is known.

        #self.partpos = numpy.zeros(shape=(self.NMax), dtype=numpy.int_)
        #SD.partpos   = self.partpos.ctypes.data
        #self.partpos[:] = S12_EMPTYSITE

    def _initArrays(self, lattSize, connMax):
        """Function to initilize various data structures.

        This is used with the various grid initialization methods.
        """
        self.lattSize = lattSize
        self.connMax = connMax
        maxTypes = connMax + 1  # we don't use zero  [1, N+1)

        # nneighbors   (the only reason these comments are here is to
                      # make it easier to pick out the separations by eye)
        self.__dict__["nneighbors"]=numpy.zeros(shape=(lattSize),
                                          dtype=numpy.int_)
        self.SD.nneighbors = self.nneighbors.ctypes.data

        # conn
        self.__dict__["conn"]=numpy.zeros(shape=(lattSize, self.connMax),
                                          dtype=numpy.int_)
        self.SD.conn = self.conn.ctypes.data
        self.conn.shape = lattSize, self.connMax

        # connN
        self.__dict__["connN"] = numpy.zeros(shape=(lattSize),
                                             dtype=numpy.int_)
        self.SD.connN = self.connN.ctypes.data

        # lattsite
        self.__dict__["lattsite"] = numpy.zeros(shape=(self.lattSize),
                                                dtype=numpy.int_)
        self.SD.lattsite = self.lattsite.ctypes.data
        self.lattsite[:] = S12_EMPTYSITE

        # atomtype                                        # lattSize is NMax
        self.__dict__["atomtype"] = numpy.zeros(shape=(self.lattSize),
                                                dtype=numpy.int_)
        self.SD.atomtype = self.atomtype.ctypes.data
        self.atomtype[:] = S12_EMPTYSITE

        # atomtype                                       # lattSize is NMax
        self.__dict__["atompos"] = numpy.zeros(shape=(self.lattSize),
                                               dtype=numpy.int_)
        self.SD.atompos = self.atompos.ctypes.data
        self.atompos[:] = S12_EMPTYSITE

        # atomtype                                        # lattSize is NMax
        self.__dict__["ntype"] = numpy.zeros(shape=(maxTypes),
                                             dtype=numpy.int_)
        self.SD.ntype = self.ntype.ctypes.data

        
        

    def addParticle(self, pos, type_=1):
        """Add particle at lattice site `pos`.

        This does minimal checks, and then calls the proper C
        function.  You really should be doing this from C, not python."""
        if self.N >= self.lattSize:
            print "ERROR: No more room to insert particles"
        if pos < 0 or pos >= self.lattSize:
            print "ERROR: Lattice position is out of bounds"
        if self.lattsite[pos] != S12_EMPTYSITE:
            print "ERROR: lattice site already occupied"
        #if not neighlist:   self.lattsite[pos] = type_
        #else:               self.C.addParticle(self.SD_p, pos, type_)
        self.C.addParticle(self.SD_p, pos, type_)
    def delParticle(self, pos):
        """Delete particle at lattice site `pos`"""
        #if not neighlist:  self.lattsite[pos] = S12_EMPTYSITE
        #else:              self.delParticle(pos)
        self.C.delParticle(self.SD_p, pos)
    def consistencyCheck(self, type_=3):
        """Check that all stored neighbor numbers are correct.
        """
        if False:
            print self.lattsite
            print self.atomtype
            print self.atompos
            print self.ntype
        # check total number of atoms
        # via lattsite
        Ncalc = len(self.lattsite[self.lattsite != S12_EMPTYSITE].flat)
        if Ncalc != self.N:
            raise Exception("Inconsistent N (lattsite): %s %s"%(Ncalc, self.N))
        # via atomtype
        Ncalc = len(self.atomtype[self.atomtype != S12_EMPTYSITE].flat)
        if Ncalc != self.N:
            raise Exception("Inconsistent N (atomtype): %s %s"%(Ncalc, self.N))
        # via atompos
        Ncalc = len(self.atompos[self.atompos != S12_EMPTYSITE].flat)
        if Ncalc != self.N:
            raise Exception("Inconsistent N (atompos): %s %s"%(Ncalc, self.N))
        # number of atoms of each type:
        for type_ in range(len(self.ntype)):
            if self.ntype[type_] != \
                   len(self.atomtype[self.atomtype==type_].flat):
                raise Exception("wrong atom type count: type %d"%type_)
    
        # lookup for each lattsite
        for pos in range(len(self.lattsite.flat)):
            if self.lattsite[pos] != S12_EMPTYSITE:
                if pos != self.atompos[self.lattsite[pos]]:
                    raise Exception(
                      "lookup error: %d != self.atompos[self.lattsite[%d]]"%
                      (self.atompos[self.lattsite[pos]], pos))
        # lookup for each atompos
        for i in range(self.N):
            if i != self.lattsite[self.atompos[i]]:
                raise Exception(
                    "lookup error: %d != self.lattsite[self.atompos[%d]]"%
                                (self.lattsite[self.atompos[i]], i))
        # neighborlists
        for pos in range(self.lattSize):
            if self.neighbors_pos(pos) != self.nneighbors[pos]:
                raise Exception("wrong number of neighbors cached at pos %s",
                                pos)

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
            #if not neighlist:      self.lattsite[pos] = type_
            #else:                  self.addParticle(pos, type_)
            self.addParticle(pos, type_)
            if self.energy_posNeighborhood(pos) == float("inf"):
                #if not neighlist:  self.lattsite[pos] = S12_EMPTYSITE
                #else:              self.delParticle(pos)
                self.delParticle(pos)
                continue
            inserted += 1
            #self.N += 1
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
        self.mctime += n
    def resetTime(self):
        """Sets time back to zero.

        Time is measured in cycles, and stored in self.mctime.  Cycles
        are set via self.setCycleMoves().
        """
        self.mctime = 0
    def setCycleMoves(self, shift=None, insertdel=0):
        """Sets moves executed each cycle.

        To run GCE, do:
        - set this method
        - set self.inserttype
        - set self.uVTchempotential
        """
        if shift is None:
            shift = self.N
        self.movesPerCycle = shift + insertdel

        self.cumProbAdd = (insertdel/2.) / self.movesPerCycle
        self.cumProbDel = (insertdel) / self.movesPerCycle
        #print self.cumProbAdd, self.cumProbDel
        
        
        
    def neighbors_pos(self, i):
        """Number of neighbors of a lattice site"""
        return self.C.neighbors_pos(self.SD_p, i)
    def setInsertType(self, probs):
        """Set the types of atoms to be inserted
        
        
        type_prob_array = ((t0, p0),
                           (t1, p1),
                           (t2, p2))
        """
        self.inserttype = S12_EMPTYSITE
        n = len(probs)
        maxtype = max(probs.keys())
        inserttypes_prob = numpy.empty(shape=(n), dtype=numpy.double)
        self.__dict__["inserttypes_prob"] = inserttypes_prob
        self.SD.inserttypes_prob = inserttypes_prob.ctypes.data

        inserttypes_type = numpy.empty(shape=(n), dtype=numpy.int_)
        self.__dict__["inserttypes_type"] = inserttypes_type
        self.SD.inserttypes_type = inserttypes_type.ctypes.data

        inserttypes_plookup = numpy.zeros(shape=(maxtype+1),
                                          dtype=numpy.double)
        self.__dict__["inserttypes_plookup"] = inserttypes_plookup
        self.SD.inserttypes_plookup = inserttypes_plookup.ctypes.data

        inserttypes_mulookup = numpy.zeros(shape=(maxtype+1),
                                          dtype=numpy.double)
        self.__dict__["inserttypes_mulookup"] = inserttypes_mulookup
        self.SD.inserttypes_mulookup = inserttypes_mulookup.ctypes.data

        cumulProb = 0.
        for i, (type_, prob) in enumerate(probs.iteritems()):
            mu = prob[1]    # hack twiddling with variables
            prob = prob[0]
            cumulProb += prob
            inserttypes_prob[i] = cumulProb
            inserttypes_type[i] = type_
            inserttypes_plookup[type_] = prob
            inserttypes_mulookup[type_] = mu
        if inserttypes_prob[-1] != 1:
            raise Exception, "sum of insert types must be one!"
    def getInsertType(self):
        """Get the type of particle to be inserted-- fixed or random."""
        return self.C.getInsertType(self.SD_p)
    def getPos(self):
        """Return positions of all existant particles."""
        return self.atompos[self.atompos != S12_EMPTYSITE]
    def energy(self):
        """Total energy of the system."""
        return self.C.energy(self.SD_p)
    def energy_pos(self, pos):
        """Energy due to a single lattice site."""
        return self.C.energy_pos(self.SD_p, pos)
    def energy_posNeighborhood(self, pos):
        """Energy due to a single position and all neighboring sites."""
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
        """Utility function to print checksum of current state."""
        x = ( tuple(self.lattsite.flat),
              tuple(self.conn.flat), tuple(self.connN.flat),
              self.hardness, self.lattSize, 
              )
        return hash(x)
    def numberOfType(self, type_):
        """Return the number of lattice sites with type type_
        """
        #return len(self.lattsite[self.lattsite==type_].flat)
        return len(self.atomtype[self.atomtype==type_].flat)
        #return self.ntype[type_]
    
    def anneal(self, verbose=True):
        """Slowly increase the hardness of the system until energy is zero.

        Used to go from an initial configuration which has overlaps,
        to a hard configuration with no overlaps.

        Set self.hardness = inf at the end.
        """
        hardness = 1.
        while self.energy() > 0:
            self.hardness = hardness
            self.C.cycle(self.SD_p, self.N * 2)
            if verbose:
                print "\rannealing: hardness: %5d  energy: %8.2f"% \
                      (hardness, self.energy()), \
                  self.N
            hardness += 1
            print "\033[2A\r"
            #self.printLattice()
        print
        self.hardness == float("inf")
    def chempotential(self, inserttype, store=True):
        """Chemical potential of the system, test inserting at every site.
        """
        #if self.widominserttype == S12_EMPTYSITE:
        #    raise Exception(
        #        "Must set self.inserttype to the type of particle to insert")
        
        mu = self.C.chempotential(self.SD_p, inserttype)
        if store and mu != inf:
            self.avgStore("chempotential", mu)
        return mu

    def avgStore(self, name,  value):
        """Add an average to lists, to easily compute avgs and stddevs later.

        `name` is the identifier to use when storing `value`.  Each
        time a certain quantity is computed, it should be stored with
        this function, with the same name.

        Related methods:
          - avg()
          - secondmoment()
          - stddev()
          - resetAvgs()
        """
        x = self._avgs.setdefault(name, [0., 0., 0. ])
        x[0] += 1
        x[1] += value
        x[2] += value*value
    def avg(self, name):
        """Return < x >
        """
        x = self._avgs[name]
        return x[1] / x[0]
    def avgSecondmoment(self, name):
        """Return < x**2 >.
        """
        x = self._avgs[name]
        return x[2]/x[0]
    def avgStddev(self, name, population=True):
        """Return the *population* std dev of item `name`.

        The population std dev *does* divide by N-1.  To get the
        sample std dev, set population=False.
        """
        x = self._avgs[name]
        if population:
            try:
                return math.sqrt((x[2]/x[0] - (x[1]/x[0])**2) / (x[0]-1))
            except (ValueError, ZeroDivisionError):
                return -1.
        else:
            try:
                return math.sqrt((x[2]/x[0] - (x[1]/x[0])**2))
            except (ValueError, ZeroDivisionError):
                return -1.
            
    def avgReset(self):
        """Reset all averages.
        """
        self._avgs = { }

    def _density_get(self):
        """Density of the system, N/lattSize"""
        return self.N / self.lattSize
    density = property(fget=_density_get)



    def __getattr__(self, attrname):
        """Wrapper to proxy attribute gets to the C SimData struct
        """
        if hasattr(self.SD, attrname):
            return getattr(self.SD, attrname)
        raise AttributeError("unknown ")
    def __setattr__(self, name, value):
        """Wrapper to proxy attribute sets to the C SimData struct
        """
        if hasattr(self.SD, name):
            setattr(self.SD, name, value)
        else:
            self.__dict__[name] = value
    

def correlation(lattice0, lattice1):
    lattice0 = lattice0.copy()
    lattice1 = lattice1.copy()
    
    lattice0[lattice0 != S12_EMPTYSITE] = 1
    lattice1[lattice1 != S12_EMPTYSITE] = 1

    #avg = S.density * S.density * S.lattSize * 3 * 3
    #print avg
    return float(numpy.sum(lattice0 * lattice1))# - avg
    

if __name__ == "__main__":
    from saiga12.geom.grid import Grid2d

    a, b = 4, 4
    S = Grid2d(N=a*b)
    print "x"
    S.makegrid(a, b)
    print "y"
    S.hardness = 10.
    S.addParticle(pos=1, type_=3)
    S.addParticle(pos=4, type_=3)
    S.addParticle(pos=5, type_=3)
    S.addParticle(pos=6, type_=3)
    S.addParticle(pos=9, type_=3)
    #S.addParticleRandomDensity(.7, type_=3)
    
    print "E:", S.energy()
    #sys.exit()
    latticeOrig = S.lattsite.copy()
    S.printLattice()

    print S.energy()
    print
    print S.consistencyCheck()
    S.anneal()
    #for i in range(1000):
    #    S.C.cycle(S.SD_p, S.N)
    #    print S.energy(), correlation(latticeOrig, S.lattsite, S), S.N
    #    #S.printLattice()
    #    #time.sleep(.005)
