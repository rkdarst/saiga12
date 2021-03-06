# Richard Darst, 2008


import ctypes
import math
import numpy

import saiga12
from saiga12.geom.grid import cartesianproduct

try:
    from rkddp.interact import interact
except ImportError:
    pass

# at least on Debian systems, scipy is linked to FFTW and numpy isn't.
try:
    from scipy.fftpack import fftn, ifftn
except ImportError:
    from numpy.fft import fftn, ifftn


_kvecCache = { }

class Averager(object):
    # Averaging methods
    def avgStore(self, name,  value):
        """Add an average to lists, to easily compute avgs and stddevs later.
        """
        #if not hasattr(self, "_avgs"):
        #    self._avgs = { }
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


class ConvolvingCorrelation(Averager, object):

    def __init__(self):
        self._singleFrameItems = { }
        self.items = { }
        self.frameList = [ ]

    corr_len = 1
    def doFrames(self, frames):
        frameList = self.frameList
        for f in frames:
            frameList.append(f)
            if len(frameList) > self.corr_len:
                del self.frameList[0]
            self.runItems()

    def runItems(self):
        for k, v in self._singleFrameItems.iteritems():
            self.avgStore(k, v(self.frameList[-1]))
    def addItemSingleFrame(self, name, func):
        self._singleFrameItems[name] = func



def getLattice(S, type_):
    lattice = S.lattsite.copy()
    lattice[:] = 0
    lattice[S.getPos(type_)] = 1
    lattice.shape = S.lattShape
    return lattice
def getFromCache(S, type_, function):
    """function should be fftn or ifftn
    """
    cachepar = (type_, function)
    cache = S.__dict__.setdefault("_fftcache", { })
    val = None
    if cache.has_key(cachepar):
        val, mctime = cache[cachepar]
        # if our mctime has changed (we have advanced in time),
        # then we need to regenerate the FFTn.
        if mctime != S.mctime:
            val = None
    if val is None:
        lattice = getLattice(S, type_)
        val = function(lattice)
        S._fftcache[cachepar] = val, S.mctime
    return val


def getFftArrays(S, type_, SOverlap=None):
    # This is getting the FFTs from the cache
    cachepar = (type_, SOverlap is None)
    cache = S.__dict__.setdefault("_fftcache", { })
    val = None
    if cache.has_key(cachepar):
        val, mctime = cache[cachepar]
        # if our mctime has changed (we have advanced in time),
        # then we need to regenerate the FFTn.
        if mctime != S.mctime:
            val = None

    # Do the actual regeneration of the functions:
    if val is None:
        #lattice = getLattice(S, type_)
        lattice = S.lattsite.copy()
        lattice[:] = 0
        if SOverlap is None:
            lattice[S.getPos(type_)] = 1
            norm = N
        else:
            lattice[S       .getPos(type_)] += 1
            lattice[SOverlap.getPos(type_)] += 1
            lattice[lattice != 2] = 0
            lattice[lattice == 2] = 1
            norm = N * S.densityOf(type_)

        #from rkddp import interact ; interact.interact()
        lattice.shape = S.lattShape
        val = fftn(lattice), ifftn(lattice), norm
        #val = function(lattice)
        S._fftcache[cachepar] = val, S.mctime
    return val

class NoKvecsException(saiga12.Saiga12Exception):
    pass

_kvecCache = { }

def makeKvecs(n, orthogonal, dim):
    n2 = n**2
    if orthogonal:
        if dim == 2:
            kvecs = numpy.asarray([(n , 0.  ),
                                   (0., n),
                                   ],
                                  dtype=saiga12.c_double)
            return kvecs
        elif dim == 3:
            kvecs = numpy.asarray([(n , 0., 0.),
                                   (0., n , 0.),
                                   (0., 0., n ),
                                   ],
                                  dtype=saiga12.c_double)
            return kvecs
        else:
            raise Exception("dimensions must be 2 or 3.")

    else:
        import shelve
        cachepar = str((n2, dim))
        nmax = int(math.ceil(math.sqrt(n2)))

        useCache = False
        if _kvecCache.has_key(cachepar):  # memory cache always used
            kvecs = _kvecCache[cachepar]
        else:
            cache = shelve.open("kvecCache")
            if cache.get('version', 0) < 1:
                for k in cache.keys(): del cache[k]
                cache['version'] = 1
            if cache.has_key(cachepar) and useCache:
                kvecs = cache[cachepar]
                _kvecCache[cachepar] = kvecs
            else:
                # kvecs is a list of all k-vectors consistent with our
                # magnitude.
                full = [ x for x in range(-nmax+1, nmax+1) ]
                half = [ x for x in range(0, nmax+1) ]

                if dim == 1:
                    kvecs = cartesianproduct(half, )
                elif dim == 2:
                    kvecs = cartesianproduct(full, half)
                elif dim == 3:
                    kvecs = cartesianproduct(full, full, half)
                kvecs = numpy.asarray([ _ for _ in kvecs
                                    # the following lines limits to only
                                    # orthogonal vectors, which is handled
                                    # at the very first of this method
                                    #if numpy.max(numpy.abs(_))**2 == kmag2
                                        ],
                                      dtype=saiga12.c_double)
                magnitudes2 = numpy.sum(kvecs * kvecs, axis=1)
                kvecs = kvecs[magnitudes2 == n2]
                if len(kvecs) == 0:
                    raise NoKvecsException

                if useCache:
                    cache[cachepar] = kvecs
                _kvecCache[cachepar] = kvecs   # memory cache always
            del cache
        return kvecs

def makeKvecs2D(num):
    kvecs = numpy.zeros(shape=(num, 3), dtype=saiga12.c_double)
    thetas = (numpy.arange(num) / num) * 2*math.pi
    for i, t in enumerate(thetas):
        kx = math.cos(i)
        ky = math.sin(i)
        kvecs[i, 0] = kx
        kvecs[i, 1] = ky
    return kvecs

class StructCorr(object):
    def __init__(self, S, type_, L, kmag2, orthogonal=True):
        self._type_ = type_

        # the only importance of S is the shape.
        self.L = L
        self.kmag = math.sqrt(kmag2)
        self.makeKvecs(kmag2, S, orthogonal=orthogonal)
        self.kvecsOrig = self.kvecs
        self.kvecs = self.kvecs.copy()
        self.kvecs *= (2*math.pi / L)
        self.kvecs_p = saiga12.get_cdoublep(self.kvecs)
        self.atomlist = numpy.asarray(S.atomIndexOfType(type_),
                                      dtype=ctypes.c_int)
        self.atomlist_p = saiga12.get_cintp(self.atomlist)
        self.N = len(self.atomlist) # Number of _included_ atoms...


        #self.makeCoordLookup(S)
        S._corrfunc_makeCoordLookup(self)
        self.physicalShape = numpy.asarray(S.physicalShape,
                                           dtype=saiga12.c_double)
        self.physicalShape_p = saiga12.get_cdoublep(self.physicalShape)

        #self.SkArrayAvgs    = numpy.zeros(shape=len(kvecs),
        #                                  dtype=saiga12.c_double)
        self._SkArrayByKvec    = numpy.zeros(shape=len(self.kvecs),
                                             dtype=saiga12.c_double)
        self._SkArrayByKvecTMP = numpy.zeros(shape=len(self.kvecs),
                                             dtype=saiga12.c_double)
        self._SkArrayByKvecTMP_p = saiga12.get_cdoublep(self._SkArrayByKvecTMP)
        self._SkArrayByAtom    = numpy.zeros(shape=self.N,
                                             dtype=saiga12.c_double)
        self._SkArrayByAtomTMP = numpy.zeros(shape=self.N,
                                             dtype=saiga12.c_double)
        self._SkArrayByAtomTMP_p = saiga12.get_cdoublep(self._SkArrayByAtomTMP)
        self.reset()
        # we rely on people to not go non-contiguouizing this later.
        assert self.kvecs.flags.c_contiguous, "kvecs is not c_contiguous"

    def reset(self):
        #self.avgReset()
        self._SkTotal = 0.
        self._niterSk = 0
        #self.SkArrayAvgs[:] = 0
        self._SkArrayByKvec[:] = 0
        self._SkArrayByAtom[:] = 0

    def makeKvecs(self, kmag2, S, orthogonal=True):
        """Make k-vectors needed for space fourier transform.

        `orthogonal` is a a keyword indicating that we will only use
        (k,0,0), (0,k,0), (0,0,k) as our vectors, instead of that in
        addition to all diagonal vectors consistent with a certain
        angle.  This is because in a cubic lattice, diagonals pointing
        in different directions (0,0,3) vs (1,2,2) do not have the
        same symmetry properties.  Default to only using orthogonal
        vectors (orthogonal=True)
        """
        if orthogonal:
            kmag = int(math.sqrt(kmag2))
            # If we didn't get an exact integer used with orthogonal,
            # is is incorrect.  We could just round it, but likely the
            # user did not want that, since most likely they are using
            # all values of integer $k^2$ instead of all values of
            # integer $k$.
            if kmag != math.sqrt(kmag2):
                raise NoKvecsException
            if len(S.physicalShape) == 1:
                kvecs = numpy.asarray([(kmag, ),
                                       ],
                                      dtype=saiga12.c_double)
            elif len(S.physicalShape) == 2:
                kvecs = numpy.asarray([(kmag, 0.  ),
                                       (0.  , kmag),
                                       ],
                                      dtype=saiga12.c_double)
            elif len(S.physicalShape) == 3:
                kvecs = numpy.asarray([(kmag, 0.  , 0.  ),
                                       (0.  , kmag, 0.  ),
                                       (0.  , 0.  , kmag),
                                       ],
                                      dtype=saiga12.c_double)
        else:
            import shelve
            cachepar = str((kmag2, len(S.physicalShape)))
            kmax = int(math.ceil(math.sqrt(kmag2)))

            useCache = False
            if _kvecCache.has_key(cachepar):  # memory cache always used
                kvecs = _kvecCache[cachepar]
            else:
                cache = shelve.open("kvecCache")
                if cache.get('version', 0) < 1:
                    for k in cache.keys(): del cache[k]
                    cache['version'] = 1
                if cache.has_key(cachepar) and useCache:

                    kvecs = cache[cachepar]
                    _kvecCache[cachepar] = kvecs
                else:
                    # kvecs is a list of all k-vectors consistent with our
                    # magnitude.
                    full = [ x for x in range(-kmax+1, kmax+1) ]
                    half = [ x for x in range(0, kmax+1) ]

                    if len(S.physicalShape) == 1:
                        kvecs = cartesianproduct(half, )
                    elif len(S.physicalShape) == 2:
                        kvecs = cartesianproduct(full, half)
                    elif len(S.physicalShape) == 3:
                        kvecs = cartesianproduct(full, full, half)
                    kvecs = numpy.asarray([ _ for _ in kvecs
                                        # the following lines limits to only
                                        # orthogonal vectors, which is handled
                                        # at the very first of this method
                                        #if numpy.max(numpy.abs(_))**2 == kmag2
                                            ],
                                          dtype=saiga12.c_double)
                    magnitudes2 = numpy.sum(kvecs * kvecs, axis=1)
                    kvecs = kvecs[magnitudes2 == kmag2]
                    if len(kvecs) == 0:
                        raise NoKvecsException

                    if useCache:
                        cache[cachepar] = kvecs
                    _kvecCache[cachepar] = kvecs   # memory cache always
                del cache
        self.kvecs = kvecs

    #def makeCoordLookup(self, S):
    #    # DANGER -- only works for integer coordinates (so far!
    #    c = S.coords()
    #    c = numpy.asarray(c, dtype=saiga12.numpy_double)
    #    if not c.flags.carray:
    #        c = c.copy()
    #    self.coordLookup = c
    #    self.coordLookup_p = self.coordLookup.ctypes.data


    def calcSk(self, S, SOverlap=None):
        """Static Structure Factor

        SOverlap is used for only taking a subset where there is an
        overlap, used for the Q analysis.
        """
        # this used to be method 1
        #Sk takes MUCH longer than Fs

        type_ = self._type_
        self._niterSk += 1
        N = S.numberOfType(type_)
        lattSize = S.lattSize

        # This is getting the FFTs from the cache
        cachepar = (type_, SOverlap is None)
        cache = S.__dict__.setdefault("_fftcache", { })
        val = None
        if cache.has_key(cachepar):
            val, mctime = cache[cachepar]
            # if our mctime has changed (we have advanced in time),
            # then we need to regenerate the FFTn.
            if mctime != S.mctime:
                val = None

        # Do the actual regeneration of the functions:
        if val is None:
            #lattice = getLattice(S, type_)
            #lattice = S.lattsite.copy()
            #lattice[:] = 0
            lattice = numpy.zeros(S.lattsite.size)
            if SOverlap is None:
                lattice[S.getPos(type_)] = 1
                norm = N
                density = N / float(lattSize)
                #shape = lattice.shape
                #lattice.shape = lattice.size
                #import random
                #density = .2     # adjust
                #norm = int(math.floor(density * 15**3))
                #poss = range(15**3)
                #random.shuffle(poss)
                #poss = poss[:norm]
                #lattice[poss] = 1
                #lattice.shape = shape
            else:
                lattice[S       .getPos(type_)] += 1
                lattice[SOverlap.getPos(type_)] += 1
                lattice[lattice != 2] = 0
                lattice[lattice == 2] = 1
                #norm = N * S.densityOf(type_)
                density = N * S.densityOf(type_) / float(lattSize)

            #from rkddp import interact ; interact.interact()
            lattice.shape = S.lattShape
            val = fftn(lattice), ifftn(lattice), density
            #S._fftcache[cachepar] = val, S.mctime

        ForwardFFT, InverseFFT, density = val

        totalsum2 = 0.
        for i, k in enumerate(self.kvecsOrig):
            k = tuple(k)
            #x = ((InverseFFT[k]) * (ForwardFFT[k])).real * lattSize / norm
            x = ((InverseFFT[k]) * (ForwardFFT[k])).real /(density*(1-density))
            totalsum2 += x
            self._SkArrayByKvec[i] += x
        Sk = totalsum2 / len(self.kvecs)
        Sk = Sk.real

        self._SkTotal += Sk



    def calcFs(self, S0, S):
        """Intermediate scattering function
        """
        # this used to be method 0
        # old method: self.staticStructureFactor(S1=S1, S2=S2, method=0)

        self._SkArrayByKvecTMP[:] = 0
        self._SkArrayByAtomTMP[:] = 0
        self._niterSk += 1
        N = S.numberOfType(self._type_)
        if S.numberOfType(self._type_) != self._SkArrayByAtomTMP.size:
            raise Exception, "Number of atoms has changed... "\
                  "Fs assumes you aren't doing that."

        totalsum = S._corrfunc_calcFs(self, S0)

        Sk = totalsum / (N * len(self.kvecs))
        self._SkTotal += Sk

        self._SkArrayByKvecTMP /= N
        self._SkArrayByAtomTMP /= len(self.kvecs)

        self._SkArrayByKvec += self._SkArrayByKvecTMP
        self._SkArrayByAtom += self._SkArrayByAtomTMP
        return {'average':Sk,
                'byKvec':self._SkArrayByKvecTMP,
                'byAtom':self._SkArrayByAtomTMP}

    def Sk(self):
        return self._SkTotal /  self._niterSk
    def SkArraysByKvec(self):
        return self._SkArrayByKvec / self._niterSk
    def SkArraysByAtom(self):
        return self._SkArrayByAtom / self._niterSk


    def staticStructureFactor(self, S1, method, S2=None):
        print "don't use the staticStructureFactor method anymore."
        if method == 0:
            self.calcFs(S1, S2)
        elif method == 1:
            if S2 is not None and S1 is not S2:
                raise
            self.calcSk(S1)
        else:
            raise

    def fourpoint(self, S1, S2, qprime, dosin=False):
        """Density definition (old way).  Should be named fourpointDensity_0
        """
        q = 2 * qprime * math.pi / self.L
        q = numpy.asarray(((q, 0, 0),
                           (0, q, 0),
                           (0, 0, q),), dtype=saiga12.numpy_double)
        A = ctypes.c_double(0)
        B = ctypes.c_double(0)
        C = ctypes.c_double(0)

        c1 = c2 = self.coordLookup_p
        flags = 0
        if dosin:
            flags |= saiga12.S12_FLAG_DOSIN

        for qvec in q:
            q_p = q.ctypes.data
            assert q.flags.carray
            S1.C.fourpoint(S1.SD_p, S2.SD_p,
                       c1, c2,
                       self.physicalShape_p,len(S1.physicalShape), self._type_,
                       self.kvecs_p, len(self.kvecs),
                       q_p, # q
                       ctypes.byref(A), ctypes.byref(B), ctypes.byref(C),
                       self._SkArrayByKvecTMP_p, self._SkArrayByAtomTMP_p,
                       flags, # flags.
                       )
        return A.value/3., B.value/3., C.value/3.


class StructCorrList(object):
    def __init__(self, S, type_=saiga12.S12_TYPE_ANY,
                 kmag2s=None, kmags=None, L=None, orthogonal=True):
        if L is None:
            L = S.L
        if kmags is not None:
            kmag2s = [ kmag*kmag for kmag in kmags ]

        SsfList = [ ]
        SsfDict = { }

        for kmag2 in kmag2s:
            try:
                Ssf = StructCorr(kmag2=kmag2, S=S, type_=type_, L=L,
                                 orthogonal=orthogonal)
            except NoKvecsException:
                continue
            SsfList.append(Ssf)
            SsfDict[math.sqrt(kmag2)] = Ssf
        self.SsfList = SsfList
        self.SsfDict = SsfDict

    def calcSk(self, S, SOverlap=None):
        for Ssf in self.SsfList:
            Ssf.calcSk(S=S, SOverlap=SOverlap)
    def calcFs(self, S0, S=None):
        for Ssf in self.SsfList:
            Ssf.calcFs(S0=S0, S=S)
    def reset(self):
        for Ssf in self.SsfList:
            Ssf.reset()
    def kmags(self):
        return tuple( Ssf.kmag      for Ssf in self.SsfList )
    def SkAverages(self):
        return tuple( Ssf.Sk()      for Ssf in self.SsfList )
    def SkArraysByKvec(self):
        return tuple( Ssf.SkArraysByKvec() for Ssf in self.SsfList )
    def SkArraysByAtom(self):
        return tuple( Ssf.SkArraysByAtom() for Ssf in self.SsfList )


    def plotSk(self, xlab='', ylab='', type='l', col=1, ylim=(0, 15),
               **other):
        from rpy import r
        r.plot(self.kmags(), self.SkAverages(),
               xlab=xlab, ylab=ylab, type=type, col=col,
               ylim=ylim, **other
               )
    def plotVertical(self, lty=3, col='lightgray', **other):
        from rpy import r
        for kmag in self.kmags():
            r.abline(v=kmag, lty=lty, col=col, **other)
    def plotSkByKvec(self, col='blue', pch='+', function=numpy.average,
                     **other):
        from rpy import r
        for kmag, SkArray in zip(self.kmags(), self.SkArraysByKvec()):
            if function is not None:
                SkArray = function(SkArray)
            try: len_ = len(SkArray)
            except: len_ = 1
            r.points(x=[kmag]*len_, y=SkArray,
                     col=col, pch=pch, **other)
    def plotSkByAtom(self, col='green', pch='x', function=numpy.average,
                     **other):
        from rpy import r
        for kmag, SkArray in zip(self.kmags(), self.SkArraysByAtom()):
            if function is not None:
                SkArray = function(SkArray)
            try: len_ = len(SkArray)
            except: len_ = 1
            r.points(x=[kmag]*len_, y=SkArray,
                     col=col, pch=pch, **other)


def makeSsfList(S, type_, kmag2s, L):
    """Create a list of Static Structure Factors to the parameters,
    mainly being the different kmag's.
    """
    SsfList = [ ]
    for kmag2 in kmag2s:
        Ssf = StructCorr(kmag2=kmag2,
                         S=S,
                         type_=type_)
        if len(Ssf.kvecs) == 0:
            continue
        #Ssf.kvecsOrig = Ssf.kvecs
        #Ssf.kvecs = Ssf.kvecs.copy()
        #Ssf.kvecs *= (2*math.pi / L)
        SsfList.append(Ssf)
    return SsfList



class FourpointDensity(object):
    def __init__(self, S, type_=saiga12.S12_TYPE_ANY):
        #self._siteCorr4 = numpy.zeros(shape=(S.lattSize, S.lattSize),
        #                              dtype=saiga12.c_int)
        #self._siteCorr2 = numpy.zeros(shape=(S.lattSize, ),
        #                              dtype=saiga12.c_int)
        #self._siteCorr4_p = self._siteCorr4.ctypes.data_as(saiga12.c_int_p)
        #self._siteCorr2_p = self._siteCorr2.ctypes.data_as(saiga12.c_int_p)
        self._type = type_
        self.S = S
        self.n = 0
        self.Q = saiga12.util.Averager()
    def __repr__(self):
        return "<%s object (0x%x)>"%\
               (self.__class__.__name__,id(self))
    def __iadd__(self, (S0, S)):
        #S.C.fourpointDensity(S0.SD_p, S.SD_p, self._type,
        #                     self._siteCorr4_p, self._siteCorr2_p,
        #                     0) # flags
        Q = S.C.Q(S0.SD_p, S.SD_p, self._type, 0)
#        print "Chi4Density:", Q
        self.Q += Q
        self.n += 1
        return self
    def resultArray(self):
        N       = self.S.N
        density = self.S.density

        nn4 = self._siteCorr4/float(self.n)
        nn2 = self._siteCorr2/float(self.n)
        d = self.S.density
        L = self.S.lattSize
        return ((1. / N) *
               numpy.sum((nn4 - numpy.outer(nn2,nn2))))
    def resultQ(self):
        return (1./self.S.N)*self.Q.var
    result = resultQ



class FourpointData(object):
    pass
class FourpointFs(object):
    """New S4 and chi4 calculation method.

    Berthier, Biroli, Bouchaud, Kob, Miyazaki, Reichman.  J Chem Phys
    126 184503 (2007).
    """
    def __init__(self, S, qmag, type_):
        self._type = type_
        self.ntype = S.numberOfType(self._type)
        # CorrectMode = True: do cos**2 + sin**2 method.  Add them at the end
        #             = False: use only cos**2, which simplifies.
        self.CorrectMode = True

        kmag = self.kmag = 10
        qmag = self.qmag = qmag
        # Two independent modes: sin^2 and cos^2, need to store
        # their data separately and add at the end.
        def makedat(qmag=qmag):
            dat = FourpointData()
            dat.corr = saiga12.corrfunc.StructCorr(
                S=S, type_=self._type, L=S.lattShape[0],
                kmag2=kmag**2,
                orthogonal=True)
            for name in ("A", "B", "C", "AmBC2", "AA", "AB", "BB"):
                setattr(dat, name, saiga12.util.Averager())
            return dat
        def makecorr():
            return makedat(), makedat()
        if not self.CorrectMode:
            makecorr = makedat
        self.Sdat = makecorr()
        #self.makecorr = makecorr

    def __iadd__(self, (S0, S)):
        # The loop over deltaT is above this function.
        qmag = self.qmag
        Sdat_orig = self.Sdat
        for dosin in (True, False):
            if self.CorrectMode:
                Sdat = Sdat_orig[int(dosin)]
            else:
                if dosin: continue
                Sdat = Sdat_orig
            A,B,C = Sdat.corr.fourpoint(S0, S, qprime=qmag, dosin=dosin)
            #A /= self.ntype ; B /= self.ntype ; C /= self.ntype
            #if qmag == 0: print A, B, C ;
            if (not self.CorrectMode or not dosin) and qmag == 0: assert B==self.ntype
            Sdat.A.add(A) ; Sdat.B.add(B) ; Sdat.C.add(C)
            Sdat.AmBC2.add((A-B*C)**2)
            Sdat.AA.add(A**2)
            Sdat.AB.add(A*B)
            Sdat.BB.add(B**2)
            #print A,B,C
            n = self.ntype
            #print A/n - B*(C/n)/n
        return self

    def S4(self):
        def S4(Sdat):
            """S4 from one data set: note"""
            # method 1
            A,B,C = Sdat.A.mean, Sdat.B.mean, Sdat.C.mean
            Fs = C/self.ntype
            S4 = (A - B*(C))**2 / self.ntype**2
            #return S4
            # method 2
            S4 = Sdat.AA.mean - 2*Sdat.AB.mean*Fs + \
                 Sdat.BB.mean*Fs**2
            S4 /= self.ntype**2
            S4 *= self.ntype
            return S4
            # method 3
            #S4 = Sdat.AmBC2.mean  # wrong
            #return S4
        if self.CorrectMode: S4 = S4(self.Sdat[0]) + S4(self.Sdat[1])
        else:                S4 = S4(self.Sdat)
        if self.qmag == 0:
            chi4 = self.chi4()
            assert chi4==S4 or (chi4-S4)/(chi4+S4) < 1e-5, \
                   "Error: S4(q=0) != chi4  (chi4=%s S4=%s)"%(chi4, S4)
        return S4

    def accumulators(self):
        if self.CorrectMode:
            return (.5*(self.Sdat[0].A.mean+self.Sdat[1].A.mean),
                    .5*(self.Sdat[0].B.mean+self.Sdat[1].B.mean),
                    .5*(self.Sdat[0].C.mean+self.Sdat[1].C.mean),)
        else:
            return self.Sdat.A.mean, self.Sdat.B.mean, self.Sdat.C.mean

    def Fs(self):
        if self.CorrectMode:
            assert self.Sdat[0].C.mean == self.Sdat[1].C.mean or True
            return .5*(self.Sdat[0].C.mean + self.Sdat[1].C.mean)
        else:
            return self.Sdat.C.mean
    def chi4(self):
        if self.CorrectMode:
            chi4_0 = self.Sdat[0].C.var/self.ntype
            chi4_1 = self.Sdat[1].C.var/self.ntype
            n = self.ntype
            #print "*", chi4_0, chi4_1, self.Sdat[0].C.mean, self.Sdat[1].C.mean
            chi4 = .5*(chi4_0+chi4_1) + numpy.var((self.Sdat[0].C.mean/n,
                                                   self.Sdat[1].C.mean/n))
        else:
            chi4 = self.Sdat.C.var/self.ntype
        #chi4 = self.FPData.C.var/self.ntype
        return chi4
class FourpointFsList(object):
    def __init__(self, S, qmags, type_):
        self.FPList = [ ]
        self.FPDict = { }
        for q in qmags:
            FP = FourpointFs(S, q, type_)
            self.FPList.append(FP)
            self.FPDict[q] = FP
    def __iadd__(self, (S0, S)):
        for FP in self.FPList:
            FP += S0, S
        return self
    def qList(self):
        return [ FP.qmag for FP in self.FPList ]
    def S4List(self):
        return [ FP.S4() for FP in self.FPList ]
    def S4Dict(self):
        dict( [(q, FP.S4()) for (q, FP) in self.FPDict.iteritems() ])
    def FsList(self):
        return [ FP.Fs() for FP in self.FPList ]



class SpinGlass(object):
    def __init__(self, S, types=(saiga12.S12_TYPE_ANY, saiga12.S12_TYPE_ANY)):
        gridSize = S.lattSize
        self.orient = False
        if S.orient is not None:
            self._nConn = S.connMax
            gridSize = S.lattSize*self._nConn
            self.orient = True
        self._siteCorrelation = numpy.zeros(shape=(gridSize, gridSize),
                                            dtype=saiga12.c_int)
        self._siteCorrelation_p = self._siteCorrelation.ctypes.data_as(
                                                               saiga12.c_int_p)
        self.types = types
        self.S = S
        self.n = 0
    def __repr__(self):
        return "<%s object (0x%x), types (%s,%s)>"%\
               (self.__class__.__name__,id(self),self.types[0],self.types[1])

    def __iadd__(self, (S0, S)):
        S.C.spinGlass(S0.SD_p, S.SD_p, self.types[0], self.types[1],
                      self._siteCorrelation_p,
                      0) # flags
        self.n += 1
        return self
    def result(self):
        n0   = self.S.numberOfType(self.types[0])
        rho0 = self.S.densityOf(   self.types[0])
        n1   = self.S.numberOfType(self.types[1])
        rho1 = self.S.densityOf(   self.types[1])

        #print self.types, n0, rho0, n1, rho1
        #print self._siteCorrelation
        #return ((1. / math.sqrt(n0*n1)) *
        #       numpy.sum((self._siteCorrelation/float(self.n) - rho0*rho1)**2))
        #print self.n
        d = self.S.density
        L = self.S.lattSize
        mean = rho0*rho1
        prefactor = 1. / math.sqrt(n0*n1)
        if self.orient:
            mean /= self._nConn*self._nConn
            #prefactor
        #print "*", "%.2f"%d, nn.sum(), nn.diagonal().sum(), \
        #      "% -12.6g"%((((nn.diagonal()-rho0*rho1)**2).sum()/L) - (d*(1-d))**2), \
        #      "% -12.6g"%((((nn.diagonal()-rho0*rho1)**2).sum()/L)), \
        #      "% -12.6g"%((((nn-rho0*rho1)**2).sum() - ((nn.diagonal()-rho0*rho1)**2).sum())/L), \
        #      "% -12.6g"%((((nn-rho0*rho1)**2).sum()/L)) #, \
              #"% -12.6g"%((((nn-rho0*rho1)**2).sum() - ((nn.diagonal()-rho0*rho1)**2).sum())/L) , \
              #"% -12.6g"%(nn.diagonal()-rho0).mean(), \
              #"% -12.6g"%((nn.sum() - nn.diagonal().sum())/(nn.size-nn.diagonal().size)), \
              #"% -12.6g"%(((nn.sum() - nn.diagonal().sum())/(nn.size-nn.diagonal().size)
              #            - (rho0)*(rho1))) #, \
        #      "% -12.6g"%nn.mean(), \
        #      "% -12.6g"%(nn - rho0*rho1).mean(), \
        #      "% -12.6g"%nn.var()#, \
        #      "% -12.6g"%((nn - rho0*rho1)**2).mean()
        #print

        ##nn = self._siteCorrelation/float(self.n)
        ##nn = numpy.subtract(nn, mean, nn)
        ##nn = numpy.square(nn, nn)
        ##result = ( prefactor * numpy.sum(nn))
        ##print "1"*10, result

        nn = self._siteCorrelation
        result = self.S.C.spinGlass_sumArray(
            nn.ctypes.data_as(saiga12.c_int_p), nn.size,
            self.n, mean, prefactor
            )
        #print "2"*10, result

        del nn
        return result
        #return ( prefactor * numpy.sum((nn - mean)**2))

class SpinGlassList(object):
    def __init__(self, S, types=None):
        if types is None:
            types = [t for t,n in enumerate(S.ntype) if n > 0 ]
        self.types = types
        # Make a matrix to hold all of our info.  This convoluted
        # procedure is needed since otherwise we end up with shallow
        # copies.
        typeMatrix = [ ]
        for _ in range(len(types)):
            typeMatrix.append([ None ] * len(types))
        # A simple list to hold all elements.  This only works if
        # commutability of the types holds.
        typeList = [ ]
        # Go through and build up all of our lists.
        for i,type0 in enumerate(types):
            for j,type1 in list(enumerate(types))[i:]:
                SG = SpinGlass(S, types=(type0, type1))
                typeMatrix[i][j] = SG
                typeMatrix[j][i] = SG
                if i <= j:
                    typeList += [SG]
        self.typeMatrix = typeMatrix
        self.typeList = typeList

    def __iadd__(self, (S0, S)):
        for SG in self.typeList:
            SG += S0, S
        return self
    def typePairs(self):
        return [ sg.types for sg in self.typeList ]
    def resultMatrix(self):
        result = [ ]
        for row in self.typeMatrix:
            result.append( [x.result() for x in row] )
        return result
    def resultList(self):
        return [ x.result() for x in self.typeList]



if __name__ == "__main__":
    import sys

    from rpy import r

    fname = sys.argv[2]
    try:
        type_ = int(sys.argv[1])
    except ValueError:
        type_ = saiga12.S12_TYPE_ANY

    S = saiga12.io.io_open(file(fname))

    maxAllowedK =  50
    maxKmag = min(max(S.lattShape), maxAllowedK)
    kmag2s = range(1, maxKmag**2)
    kmags = [math.sqrt(x) for x in kmag2s]

    StructCorr = StructCorrList(S=S, type_=type_,
                                kmag2s=kmag2s,)
    L = S.lattShape[0]
    for S in (S, ):
        StructCorr.calcSk(S)

    v_kmag = StructCorr.kmags()
    v_sk = StructCorr.SkAverages()

    #thisRun.sort()
    print "top modes:"

    print "%6s %6s %7s %9s"%('kmag', 'Sk', 'kmag/L', 'L/kmag')
    for Sk, kmag in zip(v_sk, v_kmag)[-10:]:
        print "%6.3f %6.3f %7.5f %9.5f"%(kmag, Sk, kmag/L, L/kmag)

    r.plot(v_kmag, v_sk,
           xlab="", ylab="", type="l",
           ylim=(0., 15)
           #ylim=(0., 15000)
           )
    # plots a line at each k-vector point.
    for kmag in v_kmag:
        r.abline(v=kmag, lty=3, col="lightgray")
    # plots each individual k-vector
    for kmag, SkArray in zip(v_kmag, StructCorr.SkArrays()):
        r.points(x=[kmag]*len(SkArray), y=SkArray,
                 col="blue", pch="x")
    # plots the standard deviation of the by-kvector lists.
    for kmag, SkArray in zip(v_kmag, StructCorr.SkArrays()):
        r.points(x=kmag, y=numpy.std(SkArray),
                 col="red", pch="s")
    # plots it using the average of all k-vectors -- this should line up
    # exactly
    #for kmag, SkArray in zip(v_kmag, StructCorr.SkArrays()):
    #    r.points(x=kmag, y=numpy.average(SkArray),
    #             col="green", pch="x")


    import code ; code.interact(local=locals(), banner="" )
