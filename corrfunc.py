# Richard Darst, 2008


import math
import numpy

import saiga12
from geom.grid import cartesianproduct, coords

try:
    from rkddp.interact import interact
except ImportError:
    pass

class Averager(object):
    # Averaging methods
    def avgStore(self, name,  value):
        """Add an average to lists, to easily compute avgs and stddevs later.
        """
        if not hasattr(self, "_avgs"):
            self._avgs = { }
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
    

class StructCorr(Averager, object):
    def __init__(self, S, kmag=None, kmag2=None, type_=None):
        self._niterSk = 0
        self._niterNneigh = 0
        if type_ == None:
            type_ = saiga12.S12_TYPE_ANY
        self._type_ = type_

        # the only importance of S is the shape.
        if kmag2 is not None:
            kmag = math.sqrt(kmag2)
        self.kmag = kmag
        self.makeKvecs(kmag)
        self.makeCoordLookup(S)

    def makeKvecs(self, kmag):
        kmax = int(math.ceil(kmag))

        # kvecs is a list of all k-vectors consistent with our magnitude.
        full = [ x for x in range(-kmax+1, kmax+1) ]
        half = [ x for x in range(0, kmax+1) ]

        kvecs = cartesianproduct(full, full, half)
        kvecs = numpy.asarray([ _ for _ in kvecs ], dtype=saiga12.c_double)

        magnitudes2 = numpy.sum(kvecs * kvecs, axis=1)
        kvecs = kvecs[magnitudes2 == kmag*kmag]

        self.kvecs = kvecs
        self.SkArray_ = numpy.zeros(shape=len(kvecs),
                                    dtype=saiga12.c_double)
        self.SkArrayAvgs = numpy.zeros(shape=len(kvecs),
                                       dtype=saiga12.c_double)


    def makeCoordLookup(self, S):
        c = S.coords(numpy.arange(S.lattSize))
        if not c.flags.carray:
            c = c.copy()
        #print c, c.flags
        self.coordLookup = c
        


    def staticStructureFactor(self, S1):
        type_ = self._type_
        self.SkArray_[:] = 0
        self._niterSk += 1

        lattShape = numpy.asarray(S1.lattShape,
                                  dtype=saiga12.c_double).ctypes.data

        #totalsum = S1.C.calc_structfact(S1.SD_p, self.kvecs.ctypes.data,
        #                                len(self.kvecs), type_,
        #                                self.coordLookup.ctypes.data,
        #                                lattShape,
        #                                self.SkArray_.ctypes.data)
        totalsum = 0
        N = S1.numberOfType(type_)

        Sk = totalsum / ((N * (N-1)/2.) * len(self.kvecs))


        # Do it using FFTn.

        # this is caching code.  
        # fftpar is a parameter relating to our fourier transform.
        # For example, we will need different transforms for different
        # particle typs.
        fftpar = (type_, )
        fftcache = S1.__dict__.setdefault("_fftcache", { })
        fftn = None
        # Do we have the FFT already stored in the cache?
        if fftcache.has_key(fftpar):
            fftn, ifftn, mctime = S1._fftcache[fftpar]
            # if our mctime has changed (we have advanced in time),
            # then we need to regenerate the FFTn.
            if mctime != S1.mctime:
                fftn = None
        if fftn is None:
            print "calculating fft"
            lattsite = S1.lattsite.copy()
            lattsite[:] = 0
            lattsite[S1.atompos[S1.atomtype == type_]] = 1
            lattsite.shape = S1.lattShape
            fftn = numpy.fft.fftn(lattsite)
            ifftn = numpy.fft.ifftn(lattsite)
            S1._fftcache[fftpar] = fftn, ifftn, S1.mctime

        totalsum2 = 0.
        for i, k in enumerate(self.kvecsOrig):
            k = tuple(k)
            #print k
            x = ((ifftn[k]) * (fftn[k])).real * S1.lattSize / N #(N * (N-1)/2.)
            totalsum2 += x
            self.SkArray_[i] += x
        #Sk2 = totalsum2 / ((N * (N-1)/2.) * len(self.kvecs))
        Sk2 = totalsum2 / len(self.kvecs)
        Sk2 = Sk2.real #abs(Sk2)
        #print Sk, Sk2, Sk2/Sk


        Sk = Sk2
        self.SkArrayAvgs += self.SkArray_ #/ ((N * (N-1)/2.))
        
        #interact()

        self.avgStore('Sk', Sk)
        return Sk

    def Sk(self):
        return self.avg('Sk')
    def SkArray(self):
        return self.SkArrayAvgs / self._niterSk

def makeSsfList(S, type_, kmag2s, L=15.):
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
        Ssf.kvecsOrig = Ssf.kvecs.copy()
        Ssf.kvecs *= (2*math.pi / L)
        SsfList.append(Ssf)
    return SsfList

if __name__ == "__main__":
    import sys

    from rpy import r
    
    fname = sys.argv[2]
    type_ = int(sys.argv[1])

    S = saiga12.io.io_open(file(fname))
    SsfList = [ ]

    for kmag2 in range(1, 50):
        kmag = math.sqrt(kmag2)
    
        Ssf = StructCorr(kmag2=kmag2,
                         S=S,
                         type_=type_)
        if len(Ssf.kvecs) == 0:
            continue
        Ssf.kvecsOrig = Ssf.kvecs.copy()
        Ssf.kvecs *= (2*math.pi / 15.)
        SsfList.append(Ssf)
    
    
    thisRun = [ ]
    
    for S in (S, ):
        for Ssf in SsfList:
            #s =
            Ssf.staticStructureFactor(S)
            #Ssf.avgStore('ssf', s)
            #print Ssf._avgs
        
    for Ssf in SsfList:
        thisRun.append((Ssf.kmag, Ssf.Sk()))
    v_kmag = zip(*thisRun)[0]
    v_ssf = zip(*thisRun)[1]
    
    i = 0
    #if mu > 10: pch = "x"
    #else:       pch = "+" #1
    pch = 1
    r.plot(v_kmag, v_ssf,
           xlab="", ylab="", type="l", col=i+1,
           ylim=(0., 15)
           #ylim=(0., 15000)
           )
    for kmag in v_kmag:
        r.abline(v=kmag, lty=3, col="lightgray")

    for Ssf in SsfList:
        r.points(x=[Ssf.kmag]*len(Ssf.kvecs), y=Ssf.SkArray(),
                 col=i+1, pch=pch)
    
    
    import code ; code.interact(local=locals(), banner="" )
