# Richard Darst, 2008


import math
import numpy

import saiga12
from geom.grid import cartesianproduct, coords

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
            #print self._avgs
    def addItemSingleFrame(self, name, func):
        self._singleFrameItems[name] = func

class StructCorrAll(Averager, object):
    def __init__(self, S, ):
        pass


def getLattice(S, type_):
    lattice = S.lattsite.copy()
    lattice[:] = 0
    lattice[S.atompos[S.atomtype == type_]] = 1
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
        #print "calculating", function
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
        #print "calculating", function
        #lattice = getLattice(S, type_)
        lattice = S.lattsite.copy()
        lattice[:] = 0
        if SOverlap is None:
            lattice[S.atompos[S.atomtype == type_]] = 1
            norm = N
        else:
            lattice[S       .atompos[S       .atomtype == type_]] += 1
            lattice[SOverlap.atompos[SOverlap.atomtype == type_]] += 1
            lattice[lattice != 2] = 0
            lattice[lattice == 2] = 1
            norm = N * S.densityOf(type_)

        #from rkddp import interact ; interact.interact()
        lattice.shape = S.lattShape
        val = fftn(lattice), ifftn(lattice), norm
        #val = function(lattice)
        S._fftcache[cachepar] = val, S.mctime
    return val

class StructCorr(Averager, object):
    def __init__(self, S, kmag=None, kmag2=None, type_=None):
        if type_ == None:
            type_ = saiga12.S12_TYPE_ANY
        self._type_ = type_

        # the only importance of S is the shape.
        if kmag2 is not None:
            kmag = math.sqrt(kmag2)
        self.kmag = kmag
        self.makeKvecs(kmag, S)
        self.makeCoordLookup(S)
        self.reset()

    def reset(self):
        self.avgReset()
        self._niterSk = 0
        self.SkArrayAvgs[:] = 0
        
    def makeKvecs(self, kmag, S):
        import shelve
        cachepar = str((kmag, len(S.lattShape)))
        kmax = int(math.ceil(kmag))

        if _kvecCache.has_key(cachepar):
            kvecs = _kvecCache[cachepar]
        else:
            cache = shelve.open("kvecCache")
            if cache.has_key(cachepar):
            
                kvecs = cache[cachepar]
                _kvecCache[cachepar] = kvecs
            else:
                # kvecs is a list of all k-vectors consistent with our
                # magnitude.
                full = [ x for x in range(-kmax, kmax+1) ]
                half = [ x for x in range(0, kmax+1) ]
                
                if len(S.lattShape) == 2:
                    kvecs = cartesianproduct(full, half)
                elif len(S.lattShape) == 3:
                    kvecs = cartesianproduct(full, full, half)
                kvecs = numpy.asarray([ _ for _ in kvecs ],
                                      dtype=saiga12.c_double)
                magnitudes2 = numpy.sum(kvecs * kvecs, axis=1)
                kvecs = kvecs[magnitudes2 == kmag*kmag]
            
                cache[cachepar] = kvecs
                _kvecCache[cachepar] = kvecs
            del cache
        self.kvecs = kvecs
        #print kvecs
        self.SkArray_ = numpy.zeros(shape=len(kvecs),
                                    dtype=saiga12.c_double)
        self.SkArrayAvgs = numpy.zeros(shape=len(kvecs),
                                       dtype=saiga12.c_double)

    def makeCoordLookup(self, S):
        c = S.coords(numpy.arange(S.lattSize))
        c = numpy.asarray(c, dtype=saiga12.c_int)
        #print c.dtype
        if not c.flags.carray:
            c = c.copy()
        #print c, c.flags
        self.coordLookup = c
        

    def calcSk(self, S, SOverlap=None):
        #Sk takes MUCH longer than Fs
        #self.staticStructureFactor(S1=S1, S2=S2, method=1)

        type_ = self._type_
        self._niterSk += 1
        N = S.numberOfType(type_)
        self.SkArray_[:] = 0

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
            #print "calculating", function
            #lattice = getLattice(S, type_)
            lattice = S.lattsite.copy()
            lattice[:] = 0
            if SOverlap is None:
                lattice[S.atompos[S.atomtype == type_]] = 1
                norm = N
            else:
                lattice[S       .atompos[S       .atomtype == type_]] += 1
                lattice[SOverlap.atompos[SOverlap.atomtype == type_]] += 1
                lattice[lattice != 2] = 0
                lattice[lattice == 2] = 1
                norm = N * S.densityOf(type_)
        
            #from rkddp import interact ; interact.interact()
            lattice.shape = S.lattShape
            val = fftn(lattice), ifftn(lattice), norm
            #val = function(lattice)
            S._fftcache[cachepar] = val, S.mctime

        ForwardFFT, InverseFFT, norm = val
        #ForwardFFT, InverseFFT, norm = getFftArrays(S, type_, SOverlap, )
        #ForwardFFT = getFromCache(S, type_, fftn )
        #InverseFFT = getFromCache(S, type_, ifftn)
            
        totalsum2 = 0.
        lattSize = S.lattSize
        for i, k in enumerate(self.kvecsOrig):
            k = tuple(k)
            x = ((InverseFFT[k]) * (ForwardFFT[k])).real * lattSize / norm
            totalsum2 += x
            self.SkArray_[i] += x
        Sk2 = totalsum2 / len(self.kvecs)
        Sk2 = Sk2.real
        Sk = Sk2

        self.avgStore('Sk', Sk)
        self.SkArrayAvgs += self.SkArray_



    def calcFs(self, S1, S2=None):
        self.staticStructureFactor(S1=S1, S2=S2, method=0)
    def staticStructureFactor(self, S1, method, S2=None):
        type_ = self._type_
        self.SkArray_[:] = 0
        self._niterSk += 1
        if S2 == None:
            S2 = S1
        N = S1.numberOfType(type_)

        #print self.kvecs
        #print self.kvecs.dtype
        #sys.exit()

        if method == 0:
            # old method: using my custom C code:
            physicalShape = numpy.asarray(S1.physicalShape,
                                          dtype=saiga12.c_double).ctypes.data
            totalsum = S1.C.calc_structfact(S1.SD_p, S2.SD_p,
                                            self.kvecs.ctypes.data,
                                            len(self.kvecs), type_,
                                            self.coordLookup.ctypes.data,
                                            physicalShape,
                                            len(S1.physicalShape),
                                            self.SkArray_.ctypes.data)
            #print totalsum, N, len(self.kvecs), type_
            Sk = totalsum / (N * len(self.kvecs))
            #print Sk

        if method == 1:
            # Do it using FFTn.
            # caching code:
            raise Exception('use new calcSk method instead.')
            ForwardFFT =  getFromCache(S1, type_, fftn )
            InverseFFT = getFromCache(S2, type_, ifftn)
            
            totalsum2 = 0.
            lattSize = S1.lattSize
            for i, k in enumerate(self.kvecsOrig):
                k = tuple(k)
                x = ((InverseFFT[k]) * (ForwardFFT[k])).real * lattSize / N
                totalsum2 += x
                self.SkArray_[i] += x
            Sk2 = totalsum2 / len(self.kvecs)
            Sk2 = Sk2.real
            Sk = Sk2

        self.avgStore('Sk', Sk)
        #print self._avgs

        self.SkArrayAvgs += self.SkArray_ #/ ((N * (N-1)/2.))
        #print Sk
        #print self.SkArray_
        return Sk

    def Sk(self):
        return self.avg('Sk')
    def SkArray(self):
        return self.SkArrayAvgs / self._niterSk

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
        Ssf.kvecsOrig = Ssf.kvecs
        Ssf.kvecs = Ssf.kvecs.copy()
        Ssf.kvecs *= (2*math.pi / L)
        SsfList.append(Ssf)
    return SsfList

if __name__ == "__main__":
    import sys

    from rpy import r
    
    fname = sys.argv[2]
    type_ = int(sys.argv[1])

    S = saiga12.io.io_open(file(fname))

    maxAllowedK =  50
    maxKmag = min(max(S.lattShape), maxAllowedK)
    kmag2s = range(1, maxKmag**2)
    kmags = [math.sqrt(x) for x in kmag2s]
    L = S.lattShape[0]
    SsfList = makeSsfList(S=S, type_=type_,
                          kmag2s=kmag2s,
                          L=L)
    for S in (S, ):
        for Ssf in SsfList:
            Ssf.calcSk(S)
        
    thisRun = [ ]
    for Ssf in SsfList:
        thisRun.append((Ssf.Sk(), Ssf.kmag))
    v_kmag = zip(*thisRun)[1]
    v_ssf = zip(*thisRun)[0]

    thisRun.sort()
    print "top modes:"
    
    print "%6s %6s %7s %9s"%('kmag', 'Sk', 'kmag/L', 'L/kmag')
    for Sk, kmag in thisRun[-10:]:
        print "%6.3f %6.3f %7.5f %9.5f"%(kmag, Sk, kmag/L, L/kmag)
    
    r.plot(v_kmag, v_ssf,
           xlab="", ylab="", type="l",
           ylim=(0., 15)
           #ylim=(0., 15000)
           )
    for kmag in v_kmag:
        r.abline(v=kmag, lty=3, col="lightgray")

    for Ssf in SsfList:
        r.points(x=[Ssf.kmag]*len(Ssf.kvecs), y=Ssf.SkArray(),
                 col="blue", pch="x")
    
    
    import code ; code.interact(local=locals(), banner="" )
