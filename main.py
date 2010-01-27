# Richard Darst, April 2008
from __future__ import division

import copy
import ctypes
import math
import numpy
import os
import sys
import time
import random
RandomSeed = 1361
random.seed(RandomSeed+165)

from saiga12 import io
from saiga12 import vibration
from saiga12 import ctccdynamics
from saiga12.common import *


c_int = ctypes.c_int
c_int_p = ctypes.POINTER(c_int)
c_double = ctypes.c_double
c_double_p = ctypes.POINTER(c_double)
c_void_p = ctypes.c_void_p
numpy_int = ctypes.c_int
numpy_double = ctypes.c_double
def get_cdoublep(array):
    return array.ctypes.data_as(c_double_p)
def get_cintp(array):
    return array.ctypes.data_as(c_int_p)

# Shared state between python and C.
class SimData(ctypes.Structure):
    _fields_ = [
        ("beta", c_double),
        ("N", c_int),             # total number of particles
        ("ntype", c_void_p),      # number of each atomtype.
        ("ntypeMax", c_int),      # Maximum allowable particle type.
        #("NMax", c_int),         # lattSize is NMax
        ("hardness", c_double),   # hardness of the hard spheres
        ("cycleMode", c_int),     # 1=MC, 2=kob-andersen
        ("energyMode", c_int),    # 1=BM, 2=all zero
        ("flags", c_int),         # flags for certain sim parameters

        ("uVTchempotential", c_double),
        ("inserttype", c_int),
        ("widominserttype", c_int),
        ("inserttypes_n", c_int),
        ("inserttypes_prob", c_void_p),
        ("inserttypes_type", c_void_p),
        ("inserttypes_plookup", c_void_p),
        ("inserttypes_mulookup", c_void_p),

        ("MLL", c_void_p),        # Move Lookup List, for event-driven dynamics
        ("MLLr",c_void_p),        # Move Lookup List, reverse
        ("MLLlen", c_int),     # MLL length
        ("MLL_down", c_void_p),     # Move Lookup List, down spins (FA model)
        ("MLLlen_down", c_int),     # MLL length
        ("MLLextraTime", c_double), # time-surplus
        
        

        ("lattSize", c_int),      # integer, length of lattice
        ("lattsite", c_void_p),   # pos -> atomnumber
        ("conn", c_void_p),       # (connMax*pos + conni) -> neighboring pos
        ("connN", c_void_p),      # pos -> N connections of that site
        ("connMax", c_int),       # int, maximum number of connections
        ("nneighbors", c_void_p), # pos -> num of neighbers
        ("atomtype", c_void_p),   # atomnumber -> atomtype
        ("atompos", c_void_p),    # atomnumber -> pos (latt position)
        ("persist", c_void_p),    # pos -> persistence func (has it moved?)
        ("orient", c_void_p),     # pos -> direction particle points.
        ("frozen", c_void_p),     # pos -> is site frozen (cannot move)?
        ("selected", c_void_p),   # pos -> is site selected for analysis?


        ("cumProbAdd", c_double),
        ("cumProbDel", c_double),
        ]
SimData_p = ctypes.POINTER(SimData)
SimDataFields = { }
for name, t in SimData._fields_:
    SimDataFields[name] = True


class Saiga12Exception(Exception):
    pass

_clibCache = { }
def getClib():
    if _clibCache.has_key("C"):
        C = _clibCache['C']
        return C
    filename = 'saiga12c.so'
    C = ctypes.cdll[os.path.join(os.path.dirname(__file__), filename)]
    C.addParticle.restype = None
    C.addParticle.argtypes = SimData_p, c_int, c_int

    C.delParticle.restype = None
    C.delParticle.argtypes = SimData_p, c_int

    cfuncs = (
        ("neighbors_pos",          c_int,    (SimData_p, c_int)),
        ("addParticle",            None,    (SimData_p, c_int, c_int)),
        ("delParticle",            None,    (SimData_p, c_int)),
        ("moveParticle",           None,    (SimData_p, c_int, c_int)),
        ("randomizeSystem",        None,    (SimData_p, c_int, c_double)),
        ("getInsertType",          c_int,    (SimData_p,)),
        ("loadStateFromSave",      None,     (SimData_p,)),
        ("energy_posLocal",             c_double, (SimData_p, c_int)),
        ("energy_pos", c_double, (SimData_p, c_int) ),
        ("energy",                 c_double, (SimData_p,)),
        ("chempotential",          c_double, (SimData_p, c_int)),
        ("chempotential_innersum", c_double, (SimData_p, c_int)),
        ("cycle",                  c_int,    (SimData_p, c_double)),
        ("calc_structfact",        c_double, (SimData_p, SimData_p, # SD1, SD2
                                        c_void_p, c_int, # *atomlist, N
                                        c_void_p, c_int, # *kvecs, Nk
                                        c_void_p, c_void_p, #  *cords *cords2
                                        c_void_p, c_int, # shape, nDim
                                        c_void_p, c_void_p,#*result,SkByAtom
                                        c_int )),        # flags
        ("fourpoint",              c_int, (SimData_p, SimData_p, # SD1, SD2
                                        c_void_p, c_void_p, #  *cords *cords2
                                        c_void_p,c_int,c_int, #shape,nDim,type
                                        c_void_p, c_int, # *kvecs, Nk
                                        c_void_p,        # q
                                        c_void_p,c_void_p,c_void_p, #*A,*B,*C
                                        c_void_p, c_void_p,#*result,SkByAtom
                                        c_int )),        # flags
        ("istructure",             c_int,    (SimData_p, )),
        ("spinGlass",              None,    (SimData_p, SimData_p, # SD0, SD1
                                             c_int, c_int, # type0, type1
                                             c_int_p,    # *siteCorrelation
                                             c_int)),    # flags
        ("spinGlass_sumArray",     c_double,(c_int_p, # *data
                                             c_int, c_int, # *size, *nPoints
                                             c_double, c_double)),# mean, pref
        ("fourpointDensity",       None,    (SimData_p, SimData_p, # SD0, SD1
                                             c_int,      # type
                                             c_int_p, c_int_p, # sc4_p, sc2_p
                                             c_int)),    # flags
        ("Q",                      c_int,   (SimData_p, SimData_p, # SD0, SD1
                                             c_int,      # type
                                             c_int)),    # flags

        ("addToMLL",               None,     (SimData_p, c_int, c_int)),
        ("removeFromMLL",          None,     (SimData_p, c_int, c_int)),

        ("EddBM_updateLatPos",     None,     (SimData_p, c_int, )),
        ("EddBM_init",             None,     (SimData_p, )),
        ("EddBM_consistencyCheck", c_int,    (SimData_p, )),
        ("EddBM_cycle",            c_int,    (SimData_p, c_double)),

        ("EddKA_updateLatPos",     None,     (SimData_p, c_int, )),
        ("EddKA_init",             None,     (SimData_p, )),
        ("EddKA_consistencyCheck", c_int,    (SimData_p, )),
        ("EddKA_cycle",            c_int,    (SimData_p, c_double)),

        ("EddFA_updateLatPos",     None,     (SimData_p, c_int, )),
        ("EddFA_init",             None,     (SimData_p, )),
        ("EddFA_consistencyCheck", c_int,    (SimData_p, )),
        ("EddFA_cycle",            c_int,    (SimData_p, c_double)),

        ("EddCTCC_updateLatPos",   None,     (SimData_p, c_int, )),
        ("EddCTCC_init",           None,     (SimData_p, )),
        ("EddCTCC_consistencyCheck",c_int,   (SimData_p, )),
        ("EddCTCC_cycle",          c_int,    (SimData_p, c_double)),

        ("EddEast_updateLatPos",     None,     (SimData_p, c_int, )),
        ("EddEast_init",             None,     (SimData_p, )),
        ("EddEast_consistencyCheck", c_int,    (SimData_p, )),
        ("EddEast_cycle",            c_int,    (SimData_p, c_double)),
        )
    for name, restype, argtypes in cfuncs:
        getattr(C, name).restype  = restype
        getattr(C, name).argtypes = argtypes
    
    C.init_gen_rand.restype = None
    C.init_gen_rand(RandomSeed+641)
    _clibCache['C'] = C
    return C

def randomSeed(data):
    C = getClib()
    C.init_gen_rand(int(hash(data)) + 165)
    random.seed(int(hash(data)) + 641)
    

class Sys(io.IOSys, vibration.SystemVibrations, ctccdynamics.CTCCDynamics,
          object):
    def __init__(self, N=None):

        SD = SimData()
        self.__dict__["SD"] = SD
        self.C = getClib()
        self.SD_p = ctypes.pointer(SD)

        # Start off lattSize as negative -- if it is used before it is
        # set, then hopefully it will raise some errors.
        self.lattSize = -1
        self.beta = 1.
        self.hardness = float("inf")
        self.cumProbAdd = 0
        self.cumProbDel = 0
        self.inserttype = S12_EMPTYSITE
        self.widominserttype = S12_EMPTYSITE
        self.naccept = 0
        self._eddEnabled = False

        self.resetTime()
        self.avgReset()
        self.setCycleMoves(shift=1)  # this must be reset once N is known.

        self.setCycleMode('montecarlo')
        self.setEnergyMode('birolimezard')

    def setCycleMode(self, cycleMode):
        """Set the dynamics cycle mode.

        Options are:
        'montecarlo'  -- Normal monte carlo dynamics.
                         Has translate and grand canonical modes.
                         Event driven dynamics is *only* for biroli-mezard
                           type systems.
        'kobandersen' -- Kob-Andersen kinetically constrained glass dynamics.
                         This automatically sets energymode to zero.
                         There *is* event-driven dynamics in this mode.
        'fredricksonandersen' -- Fredrickson-Andersen kinetically
                         constrained glass dynamics.  This automatically
                         sets energymode to zero.  There is *only*
                         event-driven dynamics in this mode.  You must make
                         the grid (to set lattSize) before enabling FA
                         dynamics.
        'ctcc'        -- On a lattice, with one directional degree of
                         freedom.
        'east'        -- East model (one sided Fredrickson-Anderson).
        """
        # self.cycleModeStr is what is used when save/reloading.
        self.cycleModeStr = cycleMode
        # self.cycleMode is what is used in actual dynamics
        if isinstance(cycleMode, int):
            if cycleMode in (S12_CYCLE_CTCCclassic, ):
                self.setCycleMode('ctcc')
                self.cycleModeStr = cycleMode # above call resets it
            self.cycleMode    = cycleMode
            return
        if cycleMode.lower() == 'montecarlo':
            self.cycleMode = S12_CYCLE_MC
            # event driven dynamics only work for biroli-mezard dynamics!
            self._eddInit = self.C.EddBM_init
            self._eddUupdateLatPos = self.C.EddBM_updateLatPos
            self._eddConsistencyCheck = self.C.EddBM_consistencyCheck
            self._eddCycle = self.C.EddBM_cycle
            
        elif cycleMode.lower() == 'kobandersen':
            self.cycleMode = S12_CYCLE_KA
            self._eddInit = self.C.EddKA_init
            self._eddUpdateLatPos = self.C.EddKA_updateLatPos
            self._eddConsistencyCheck = self.C.EddKA_consistencyCheck
            self._eddCycle = self.C.EddKA_cycle
            self.setEnergyMode('zero')
        elif cycleMode.lower() == 'fredricksonandersen':
            self.cycleMode = S12_CYCLE_FA
            self._eddInit = self.C.EddFA_init
            self._eddUpdateLatPos = self.C.EddFA_updateLatPos
            self._eddConsistencyCheck = self.C.EddFA_consistencyCheck
            self._eddCycle = self.C.EddFA_cycle
            self.setEnergyMode('zero')
            if self.lattSize <= 0:
                raise Exception, ("lattSize is not set-- you must set up "+
                                 "your arrays before enabling F-A dynamics")
            self._allocPersistArray() # automatically set to zeros
        elif cycleMode.lower() == 'east':
            self.cycleMode = S12_CYCLE_EAST
            self._eddInit = self.C.EddEast_init
            self._eddUpdateLatPos = self.C.EddEast_updateLatPos
            self._eddConsistencyCheck = self.C.EddEast_consistencyCheck
            self._eddCycle = self.C.EddEast_cycle
            self.setEnergyMode('zero')
            if self.lattSize <= 0:
                raise Exception, ("lattSize is not set-- you must set up "+
                                 "your arrays before enabling East dynamics")
            self._allocPersistArray() # automatically set to zeros
        elif cycleMode.lower() == 'ctcc':
            self.cycleMode = S12_CYCLE_CTCC
            self._eddInit =             self.C.EddCTCC_init
            self._eddUpdateLatPos =     self.C.EddCTCC_updateLatPos
            self._eddConsistencyCheck = self.C.EddCTCC_consistencyCheck
            self._eddCycle =            self.C.EddCTCC_cycle
            if self.N != 0:
                raise Exception("You must set ctcc cycle mode before you "+
                                "have any particles.")
            self._allocArray("orient", shape=(self.lattSize),
                             dtype=numpy_int)
            self.orient[:] = S12_EMPTYSITE
            self.setEnergyMode('ctcc')

        else:
            raise Exception("Unknown cycle mode: %s", cycleMode)
    def _allocPersistArray(self):
        """Allocate a persistence function array. (public method)

        Allocate just the persistence function array--- models can
        then update it.

        FA dynamics - automatically created.
        KA dynamics - implemented
        MC dynamics (BM) - implemented
        """
        self._allocArray("persist", shape=(self.lattSize),
                         dtype=numpy_int)
        # automatically set to zeros

    def setEnergyMode(self, energyMode):
        """Set the energy calculation mode.

        Options are:
        'birolimezard' -- Biroli-Mezard lattice glass model dynamics
        'zero'         -- energy is always zero (but no overlaps)
        'ctcc'        -- On a lattice, with one directional degree of
                         freedom, no overlaps.
        """
        self.energyModeStr = energyMode
        if isinstance(energyMode, int):
            if energyMode not in S12_ENERGY_AVAIL:
                raise Exception("Unknown energy mode: %s", energyMode)
            self.energyMode = energyMode
        elif energyMode.lower() == 'birolimezard':
            self.energyMode = S12_ENERGY_BM
        elif energyMode.lower() in ('zero', 'kobandersen',
                                  'fredricksonandersen'):
            self.energyMode = S12_ENERGY_ZERO
        elif energyMode.lower() == 'ctcc':
            self.energyMode = S12_ENERGY_CTCC
        else:
            raise Exception("Unknown energy mode: %s", energyMode)
    

    def _allocArray(self, name, **args):
        """Allocate an array in a way usable by both python and C.

        `name` is the name of the array to allocate: it will be
        self.`name` and self.SD.`name`.  `**args` are arguments passed
        to `numpy.zeros` to allocate the arrays.

        If you want to initilize the array to something, you have to
        do it afterwards, like this:
        self.lattsite[:] = S12_EMPTYSITE
        """
        # This code can be hard to understand, since we can't use the
        # name, but it is equivalent to this for the name `lattsite`:
        ##  self.__dict__["lattsite"] = numpy.zeros(**args)
        ##  self.SD.lattsite = self.lattsite.ctypes.data
        self.__dict__[name] = array = numpy.zeros(**args)
        setattr(self.SD, name, array.ctypes.data)

    def _initArrays(self, lattSize, connMax):
        """Function to initilize various data structures.

        This is used with the various grid initialization methods.
        """
        self.lattSize = lattSize
        self.connMax = connMax
        self.ntypeMax = connMax  # valid types are [0, connMax]

        # nneighbors
        self._allocArray("nneighbors", shape=(lattSize),
                                       dtype=numpy_int)
        # conn
        self._allocArray("conn", shape=(lattSize, self.connMax),
                                 dtype=numpy_int)
        # connN
        self._allocArray("connN", shape=(lattSize),
                                  dtype=numpy_int)
        # lattsite
        self._allocArray("lattsite", shape=(self.lattSize),
                                     dtype=numpy_int)
        self.lattsite[:] = S12_EMPTYSITE
        # atomtype                          # lattSize is NMax
        self._allocArray("atomtype", shape=(self.lattSize),
                                     dtype=numpy_int)
        self.atomtype[:] = S12_EMPTYSITE
        # atompos
        self._allocArray("atompos", shape=(self.lattSize),
                                    dtype=numpy_int)
        self.atompos[:] = S12_EMPTYSITE
        # ntype
        self._allocArray("ntype", shape=(self.ntypeMax+1),
                                  dtype=numpy_int)
    def setFrozenSites(self, sites):
        """Set or unselect certain sites as immobile

        `sites` is a list of the sites to freeze.  Note that freezing
        is by lattice site, not by particle number.  Not all
        cycle modes are guarenteed to support this.

        If `sites` is None, then unset frozen sites, so that all
        particles are mobile again.
        """
        if sites is None:
            self.flags &= ~S12_FLAG_FROZEN
            del self.frozen
            self.SD.frozen = None
            return
        if self.frozen is None:
            self._allocArray("frozen", shape=(self.lattSize), dtype=numpy_int)
        self.flags |= S12_FLAG_FROZEN
        self.frozen[:] = 0
        self.frozen[sites] = 1
    def setSelectSites(self, sites):
        """Set or unselect certain sites as ``selected''.

        `sites` is a list of the sites to select.  Selected simply
        means that certain analyzes methods can operate on only a
        subset of all particles.  Respecting this is up to each
        individual analysis calculation.

        If `sites` is None, then unset the list of selected sites.
        """
        if sites is None:
            del self.selected
            self.SD.selected = None
            return
        if self.selected is None:
            self._allocArray("selected",shape=(self.lattSize), dtype=numpy_int)
        self.selected[:] = 0
        self.selected[sites] = 1
        

    def addParticle(self, pos, type_=1):
        """Add particle at lattice site `pos`.

        This does minimal checks, and then calls the proper C
        function.  You really should be doing this from C, not python.

        The maximum type of particles is self.ntypesMax, which is currently
        0 to connMax (inclusive)."""
        if self.N >= self.lattSize:
            print "ERROR: No more room to insert particles"
        if pos < 0 or pos >= self.lattSize:
            print "ERROR: Lattice position is out of bounds"
        if self.lattsite[pos] != S12_EMPTYSITE:
            print "ERROR: lattice site already occupied"
        self.C.addParticle(self.SD_p, pos, type_)
    def delParticle(self, pos):
        """Delete particle at lattice site `pos`"""
        self.C.delParticle(self.SD_p, pos)
    def delAll(self):
        """Remove all particles from the system"""
        for pos in self.atompos[self.atompos != S12_EMPTYSITE]:
            self.delParticle(pos)
    def consistencyCheck(self, type_=3):
        """Check of internal data structures.

        Shouldn't be needed in production, only during testing of code
        changes.

        type_ is unused
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
        if type(n) not in (int, float, long):
            raise Exception("n must be of a numeric type, try .addParticles?")
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
            self.addParticle(pos, type_)
            if self.energy_pos(pos) == float("inf"):
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
            raise Exception("We already have too many particles in the system.")
            return
        self.addParticleRandom(extra, type_=type_)
    def addParticles(self, densities):
        """Add particles of various types and densities to the system.

        densities must be a dictionary, with keys being types and
        values being densities of that type:
        {type1: density1,
         type2: density2,
         ... }
        """
        cumulDensity = 0.
        for type_ in sorted(densities.keys()):
            cumulDensity += densities[type_]
            self.addParticleRandomDensity(cumulDensity, type_=type_)


    def cycle(self, n=1):
        """Run MC trial moves for once cycle.

        Cycles defined by setMoveProb()

        Timing: 45000 function calls takes one second (in an single
        unscientific test).  This excludes the time taken to do the
        simulation.
        """
        moves = int(n * self.movesPerCycle)
        if self._eddEnabled:
            self.naccept += self._eddCycle(self.SD_p, moves)
        else:
            self.naccept += self.C.cycle(self.SD_p, moves)
        self.mctime += n
    def resetTime(self):
        """Sets time back to zero.

        Time is measured in cycles, and stored in self.mctime.  Cycles
        are set via self.setCycleMoves().
        """
        self.mctime = 0
    def setCycleMoves(self, mode=None,
                      shift=0, insertdel=0):
        """Sets moves executed each cycle.

        To run GCE, do:
        - set this method
        - do one of these:
          single particle type:
            - set self.inserttype
            - set self.uVTchempotential
          multiple particle types:
            - use self.setInsertTypes()
        """
        # This is for backwards compatibility
        if type(mode) in (int, float):
            shift = mode
            mode = None
        if mode == 'canonical':
            shift = self.N
        elif mode == 'grandcanonical':
            insertdel = self.N
        if hasattr(self, "_dontSetCycleMoves"):
            # the F-A model is defined intensively already
            shift = 1
            insertdel = 0
        # If nothing is set, default to canonical.
        if shift + insertdel == 0:
            shift = self.N

        # Now set the cumulative probabilities.
        self.movesPerCycle = shift + insertdel

        self.cumProbAdd = (insertdel/2.) / self.movesPerCycle
        self.cumProbDel = (insertdel) / self.movesPerCycle
        #print self.cumProbAdd, self.cumProbDel
        
        
        
    def neighbors_pos(self, i):
        """Number of neighbors of a lattice site"""
        return self.C.neighbors_pos(self.SD_p, i)
    def setInsertType(self, probs):
        """Set the types of atoms to be inserted

        probs = {type0: (prob0, mu0 ),
                 type1: (prob1, mu1 ),
                 ... }
        """
        self._currentInsertTypes = probs
        if len(probs) == 1:
            inserttype = probs.keys()[0]
            self.inserttype = inserttype
            self.uVTchempotential = probs[inserttype][1]
            self.inserttypes_prob = None
            self.inserttypes_type = None
            self.inserttypes_plookup = None
            self.inserttypes_mulookup = None
            return

        self.inserttype = S12_EMPTYSITE
        n = len(probs)
        maxtype = max(probs.keys())
        inserttypes_prob = numpy.empty(shape=(n), dtype=numpy_double)
        self.__dict__["inserttypes_prob"] = inserttypes_prob
        self.SD.inserttypes_prob = inserttypes_prob.ctypes.data

        inserttypes_type = numpy.empty(shape=(n), dtype=numpy_int)
        self.__dict__["inserttypes_type"] = inserttypes_type
        self.SD.inserttypes_type = inserttypes_type.ctypes.data

        inserttypes_plookup = numpy.zeros(shape=(maxtype+1),
                                          dtype=numpy_double)
        self.__dict__["inserttypes_plookup"] = inserttypes_plookup
        self.SD.inserttypes_plookup = inserttypes_plookup.ctypes.data

        inserttypes_mulookup = numpy.zeros(shape=(maxtype+1),
                                          dtype=numpy_double)
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
    def getPos(self, type_=None):
        """Return positions of all existant particles.

        If type_ is None or S12_TYPE_ANY (default), return positions
        of all particles.  Otherwise, return all positions
        of that particular type."""
        if type_ == None or type_ == S12_TYPE_ANY:
            return self.atompos[self.atompos != S12_EMPTYSITE].copy()
        else:
            return self.atompos[self.atomtype == type_].copy()
    def atomIndexOfType(self, type_=None):
        """Return atom indexes (numbers) of a certain type of particle.

        type_ can be None or S12_TYPE_ANY (default) to return indexes
        of all particles, otherwise return indexes of that particular
        type."""
        if type_ == None or type_ == S12_TYPE_ANY:
            return numpy.where(self.atomtype != S12_EMPTYSITE)[0]
        else:
            return numpy.where(self.atomtype == type_)[0]
    @property
    def L(self):
        """Return the length of the system"""
        L = self.lattShape[0]
        # Check to ensure it is square
        for L2 in self.lattShape[1:]:
            assert self.lattShape[0] == L2, "Lattice is not square"
        return L
    def energy(self):
        """Total energy of the system."""
        return self.C.energy(self.SD_p)
    def energy_posLocal(self, pos):
        """Energy due to a single lattice site."""
        return self.C.energy_posLocal(self.SD_p, pos)
    def energy_pos(self, pos):
        """Energy due to a single position.

        This function is what is used in monte-carlo and other functions."""
        return self.C.energy_pos(self.SD_p, pos)
    def findInfiniteEnergy(self):
        """Utility function to find the lattice sites which have inf energy.

        Print information about infinite areas.
        """
        for pos in range(self.lattSize):
            if self.energy_pos(pos) == float("inf"):
                print pos, self.energy_posLocal(pos), \
                      self.energy_pos(pos)
                self.printLatticeLocal(pos, width=2)
    def hash(self):
        """Utility function to find checksum of current state.

        This can be used to see if two systems are in exactly the same
        state after some changes."""
        if self.persist is not None:  persist = tuple(self.persist)
        else:                         persist = None
        if self.orient is not None:   orient = tuple(self.orient)
        else:                         orient = None
        x = ( tuple(self.lattsite.flat),
              tuple(self.conn.flat), tuple(self.connN.flat),
              persist, orient,
              self.hardness, self.lattSize, tuple(self.ntype),
              self.beta, self.cumProbAdd, self.cumProbDel, self.inserttype,
              self.movesPerCycle,
              self.cycleMode,    self.energyMode,
              self.cycleModeStr, self.energyModeStr
              )
        return hash(x)
    def copy(self):
        """Return a independent copy of the System object.

        Essentially, this is the equivalent of pickling and unpickling
        the System object (I think).  This does make independent data
        arrays, so that you can propogate one in time without the
        other one changing.  To see what is preserved and what isn't,
        look at the io.py methods.

        I have seen this from some tests I did, however, I make *no
        guarentees* about how well this method works.  You should
        explicitely test yourself if the systems are connected in some
        way.
        """
        return copy.copy(self)
                                                
    def numberOfType(self, type_):
        """Return the number of lattice sites with type type_

        If type_ is None or S12_TYPE_ANY (default), return number of
        all particles.  Otherwise, return number of that particular
        type.
        """
        #return len(self.lattsite[self.lattsite==type_].flat)
        #return len(self.atomtype[self.atomtype==type_].flat)
        if type_ == None or type_ == S12_TYPE_ANY:
            return self.N
        if type_ == S12_EMPTYSITE:
            return len(self.atomtype[self.atomtype==type_].flat)
        return self.ntype[type_]
    n = numberOfType
    def densityOf(self, type_):
        """Density of a certin atomtype.

        Uses same atom selection logic as numberOfType()"""
        return float(self.numberOfType(type_)) / self.lattSize
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
            if verbose: print "\033[2A\r"
            #self.printLattice()
        if verbose: print
        self.hardness = inf
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
    def chempotential_innersum(self, inserttype):
        """Chemical potential of the system, test inserting at every site.

        mu = -log(average_of_returned_values)/S.beta
        """
        sum_ = self.C.chempotential_innersum(self.SD_p, inserttype)
        return sum_
    def persistFunction(self):
        return (self.lattSize - numpy.sum(self.persist)) / float(self.lattSize)


    def eddCheckAllowed(self, raiseException=False):
        """Check if event-driven dynamics is allowed with our current setup.

        Returns True if EDD is allowed, False if it's not.  So far,
        EDD is allowed with:
        - FA model
        - KA model
        - BM model with infinite hardness (T=0/beta=inf or hardness=inf)
        """
        if (
            self.cycleModeStr == 'fredricksonandersen' or
            self.cycleModeStr == 'east' or
            self.cycleModeStr == 'kobandersen' or
            ( self.cycleModeStr == 'montecarlo' and
              self.energyModeStr == 'birolimezard' and
              ( self.beta == inf or self.hardness == inf )
            ) or
            ( self.cycleModeStr == 'ctcc' and
              self.energyModeStr == 'ctcc' and
              ( self.beta == inf or self.hardness == inf )
            )
          ):
            return True
        if raiseException:
            raise Exception("EDD is not allowed with our paremeters")
        return False
    def eddEnable(self):
        if not self.eddCheckAllowed():
            print "can not run event-driven dynamics: ignoring"
            return
        assert self.lattSize > 0, "lattSize <= zero... is grid initialized?"
        MLLsize = self.lattSize*self.connMax
        if self.cycleModeStr in ('fredricksonandersen', 'east'):
            MLLsize = self.lattSize
            self._allocArray("MLL_down", shape=MLLsize, dtype=numpy_int)
            self.MLL_down[:] = -1
            assert self.inserttype != S12_EMPTYSITE, "inserttype is still S12_EMPTYSITE... likely you shouldn't initialize until you specify the type of particle to insert."
        self._allocArray("MLL", shape=MLLsize, dtype=numpy_int)
        self._allocArray("MLLr", shape=MLLsize, dtype=numpy_int)

        self.MLL [:] = -1
        self.MLLr[:] = -1
        self.MLLlen = 0
        self.MLLlen_down = 0
        self.MLLextraTime = -1.
        self._eddInit(self.SD_p)
        self._eddEnabled = True
    def eddDisable(self):
        if self._eddEnabled:
            del self.MLL, self.MLLr
            self._eddEnabled = False
    def eddCycle(self, nraw):
        return self._eddCycle(self.SD_p, nraw)
    def eddConsistencyCheck(self):
        if not self._eddEnabled:
            print "event-driven dynamics is not enabled! (can't c. check)"
            return 0
        x = self._eddConsistencyCheck(self.SD_p)
        #print "consistency check:", x
        return x
    def eddFindBestMode(self, n=None, nomodify=False, quiet=False):
        """Enable Event Driven Dynamics, if a test shows it is more efficient.

        `n` is the number of test cycles of a sample to take,
        defaulting to 100.

        If `copy` is false (default), then the system will be modified
        (propogated by 2n cycles total).  This is fine if you are in
        the process of equilibrating the system.  However, if you set
        `copy` to true, it will make a copy before it 
        """
        if not self.eddCheckAllowed():
            return
        if self.cycleModeStr == 'montecarlo':
            n = 100
        elif self.cycleModeStr == 'kobandersen':
            n = 100
        elif self.cycleModeStr in ('fredricksonandersen', 'east'):
            return # it should always be enabled for FA.
        if copy:
            origHash = self.hash()
            origSelf = self
            self = copy.copy(self)
        
        self.eddDisable()
        t = time.time()
        self.cycle(n)
        t1 = time.time() - t   # without EDD
        if not quiet: print 'regular moves:', t1

        self.eddEnable()
        t = time.time()
        self.cycle(n)
        t2 = time.time() - t   # with EDD
        if not quiet: print 'event driven dynamics:', t2
        if t1 < t2:
            if not quiet: print "disabling event driven dynamics"
            self.eddDisable()

        if copy:
            # enable EDD for the original one, if we need
            if self._eddEnabled: origSelf.eddEnable()
            else:                origSelf.eddDisable()
            if origSelf.hash() != origHash:
                raise Exception, "Hashes don't match after enabling even driven dynamics - the original system has changed."


    def _corrfunc_makeCoordLookup(self, StructCorr):
        # DANGER -- only works for integer coordinates (so far!
        c = self.coords()
        c = numpy.asarray(c, dtype=numpy_double)
        if not c.flags.carray:
            c = c.copy()
        StructCorr.coordLookup = c
        StructCorr.coordLookup_p = StructCorr.coordLookup.ctypes.data
    def _corrfunc_calcFs(self, StructCorr, S0):
        flags = 0
        # self is S
        if getattr(S0, 'vibEnabled', False) or S0.orient is not None:
            # Vibrations enabled: if we have vibrations enabled, we
            # have to use different coordinates for every timestep.
            flags = saiga12.S12_FLAG_VIB_ENABLED
            c1 = S0  .getCCords(returnPointer=True)
            c2 = self.getCCords(returnPointer=True)
        else:
            # No vibrations: use the same coordinates at all times
            c1 = c2 = StructCorr.coordLookup_p

        totalsum = S0.C.calc_structfact(
            S0.SD_p, self.SD_p, StructCorr.atomlist_p, StructCorr.N,
            StructCorr.kvecs_p, len(StructCorr.kvecs),
            c1, c2, StructCorr.physicalShape_p, len(StructCorr.physicalShape),
            StructCorr._SkArrayByKvecTMP_p, StructCorr._SkArrayByAtomTMP_p,
            flags)
        return totalsum


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
        if name in SimDataFields:
            return getattr(self.SD, attrname)
        raise AttributeError("No Such Attribute: %s"%attrname)
    def __setattr__(self, name, value):
        """Wrapper to proxy attribute sets to the C SimData struct
        """
        if name in SimDataFields:
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

def diff(S1, S2):
    #if type(S1) != numpy.ndarray:
    #    S1 = S1.lattsite
    #if type(S2) != numpy.ndarray:
    #    S2 = S2.lattsite
    if False:
        diffs = (S1.lattsite != S2.lattsite)
        if not numpy.any(diffs):
            print "no differences"
            return
        diffIndexes = numpy.where(diffs)
        #print diffIndexes[0]
        for i, index in enumerate(diffIndexes[0]):
            #print index, S1[index], S2[index]
            print "(at %4d) %4d -> %4d   "%(index, S1.lattsite[index],
                                                   S2.lattsite[index]),
            if i%3 == 2: print
        print
    # diff by atom numbers
    if True:
        diffs = (S1.atompos != S2.atompos)
        if not numpy.any(diffs):
            print "no differences"
            return
        diffIs = numpy.where(diffs)  # get indexes of moved atoms
        for i, index in enumerate(diffIs[0]):
            print "(#%5d)%5d ->%5d   "%(index, S1.atompos[index],
                                                  S2.atompos[index]),
            if i%3 == 2: print
        print
        
    from rkddp.interact import interact ; interact()
def distance(coord1, coord2):
    """Naive distance algorithm - does not consider periodic boundries
    """
    if len(coord1) != len(coord2):
        raise Exception("coord1 different dimensionality from coord2")
    return math.sqrt(sum([(x0-x1)**2 for (x0,x1) in zip(coord1, coord2) ]))


if __name__ == "__main__":
    # Some test scripts, but the main tests are in the tests/
    # directory.
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
