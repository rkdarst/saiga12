# Richard Darst, 2008


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

class StructCorr(object):
    def __init__(self, S, type_, L, kmag2, orthogonal=True):
        self._type_ = type_
        
        # the only importance of S is the shape.
        self.kmag = math.sqrt(kmag2)
        self.makeKvecs(kmag2, S, orthogonal=orthogonal)
        self.kvecsOrig = self.kvecs
        self.kvecs = self.kvecs.copy()
        self.kvecs *= (2*math.pi / L)
        self.kvecs_p = self.kvecs.ctypes.data
        self.makeCoordLookup(S)
        self.physicalShape = numpy.asarray(S.physicalShape,
                                           dtype=saiga12.c_double)
        self.physicalShape_p = self.physicalShape.ctypes.data
        self._whichAtoms = S.atomIndexOfType(type_)

        #self.SkArrayAvgs    = numpy.zeros(shape=len(kvecs),
        #                                  dtype=saiga12.c_double)
        self._SkArrayByKvec    = numpy.zeros(shape=len(self.kvecs),
                                             dtype=saiga12.c_double)
        self._SkArrayByKvecTMP = numpy.zeros(shape=len(self.kvecs),
                                             dtype=saiga12.c_double)
        self._SkArrayByKvecTMP_p = self._SkArrayByKvecTMP.ctypes.data
        self._SkArrayByAtom    = numpy.zeros(shape=S.N,
                                             dtype=saiga12.c_double)
        self._SkArrayByAtomTMP = numpy.zeros(shape=S.N,
                                             dtype=saiga12.c_double)
        self._SkArrayByAtomTMP_p = self._SkArrayByAtomTMP.ctypes.data
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
            kvecs = numpy.asarray([(kmag, 0.  , 0.  ),
                                   (0.  , kmag, 0.  ),
                                   (0.  , 0.  , kmag),
                                   ],
                                  dtype=saiga12.c_double)
        else:
            import shelve
            cachepar = str((kmag2, len(S.lattShape)))
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
                    
                    if len(S.lattShape) == 1:
                        kvecs = cartesianproduct(half, )
                    elif len(S.lattShape) == 2:
                        kvecs = cartesianproduct(full, half)
                    elif len(S.lattShape) == 3:
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

    def makeCoordLookup(self, S):
        # DANGER -- only works for integer coordinates (so far!
        c = S.coords()
        c = numpy.asarray(c, dtype=saiga12.numpy_double)
        if not c.flags.carray:
            c = c.copy()
        self.coordLookup = c
        self.coordLookup_p = self.coordLookup.ctypes.data
        

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
            S._fftcache[cachepar] = val, S.mctime

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
        type_ = self._type_
        self._niterSk += 1
        N = S0.numberOfType(type_)
        if S.N != self._SkArrayByAtomTMP.size:
            raise Exception, "Number of atoms has changed... "\
                  "Fs assumes you aren't doing that."

        flags = 0
        if getattr(S0, 'vibEnabled', False):
            # Vibrations enabled: if we have vibrations enabled, we
            # have to use different coordinates for every timestep.
            flags = saiga12.S12_FLAG_VIB_ENABLED
            c1 = S0.getCCords(returnPointer=True)
            c2 = S .getCCords(returnPointer=True)
        else:
            # No vibrations: use the same coordinates at all times
            c1 = c2 = self.coordLookup_p

        totalsum = S0.C.calc_structfact(S0.SD_p, S.SD_p,
                                        self.kvecs_p,
                                        len(self.kvecs), type_,
                                        #self.coordLookup_p,
                                        c1, c2,
                                        self.physicalShape_p,
                                        len(S0.physicalShape),
                                        #self._SkArray.ctypes.data,
                                        #self._SkArrayByAtom.ctypes.data)
                                        self._SkArrayByKvecTMP_p,
                                        self._SkArrayByAtomTMP_p,
                                        flags)
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
        return self._SkArrayByAtom[self._whichAtoms] / self._niterSk


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


class StructCorrList(object):
    def __init__(self, S, type_=saiga12.S12_TYPE_ANY,
                 kmag2s=None, kmags=None, L=None, orthogonal=True):
        if L is None:
            L = S.lattShape[0]  # assumes square!
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
