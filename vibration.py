# Richard Darst, 2009

import math
import numpy

import saiga12

#mode = 'orthogonal'
mode = None

class SystemVibrations(object):
    def vib_init(self, loadFrom=None):
        """Generation of the random vibration parameters.

        If loadFrom is given, load the phase and frequency parameters
        from this object.
        """
        self.vibEnabled = True
        realFrequency = .0125
        angFrequency = 2 * math.pi * realFrequency
        amplitude = .35
        #self.vib_scale = .35

        if loadFrom is None:
            #phases = lambda length:
            #    numpy.random.uniform(low=0, high=2*math.pi, size=length)
            phases = lambda length: numpy.zeros(length)
            freqs = lambda length: \
                    numpy.random.normal(loc=0,scale=angFrequency, size=length)
            ampls = lambda length: \
                numpy.random.normal(loc=0, scale=amplitude, size=length)
            if mode == "orthogonal":
                direction = numpy.random.randint(3, size=self.lattSize)
                sliceobj = numpy.arange(self.lattSize), direction

                # it is important that we use sin instead of cosine to
                # keep all non-significant ones zero.
                self.vib_phases = numpy.zeros(shape=(self.lattSize, 3))
                self.vib_freqs  = numpy.zeros(shape=(self.lattSize, 3))
                self.vib_ampls  = numpy.zeros(shape=(self.lattSize, 3))
                self.vib_phases[sliceobj] = phases(self.lattSize)
                self.vib_freqs[sliceobj]  = freqs(self.lattSize)
                self.vib_ampls[sliceobj]  = ampls(self.lattSize)
            else:
                self.vib_phases = phases(3*self.lattSize)
                self.vib_freqs  = freqs(3*self.lattSize)
                self.vib_ampls  = ampls(3*self.lattSize)
                self.vib_phases.shape = self.lattSize, 3
                self.vib_freqs.shape  = self.lattSize, 3
                self.vib_ampls.shape  = self.lattSize, 3

            #from rkddp.interact import interact ; interact()

        else:
            self.vib_phases = loadFrom.vib_phases
            self.vib_freqs = loadFrom.vib_freqs
            self.vib_ampls = loadFrom.vib_ampls

    def vib_offset(self, index=None):
        """Return an array of vibration offsets for coordinates.

        index = a slice of subset of sites to act on
        """
        if not hasattr(self, 'vib_cache'):
            # self.vib_cache is [cache_paremeter , cached_coordinates]
            self.vib_cache = [ None, self.vib_phases.copy() ]
        if self.vib_cache[0] != self.mctime:
            # do coordOffset = scale*sin(mctime*freqs + phases)  inplace
            coordOffset = self.vib_cache[1]    # load from cache
            coordOffset[:] = self.vib_freqs    # = vib_freqs
            numpy.multiply(coordOffset, self.mctime, coordOffset)
            numpy.add(     coordOffset, self.vib_phases, coordOffset)
            numpy.sin(     coordOffset, coordOffset)
            #numpy.multiply(coordOffset, self.vib_scale, coordOffset)
            numpy.multiply(coordOffset, self.vib_ampls, coordOffset)
            # update cache par:
            self.vib_cache[0] = self.mctime
        else:
            coordOffset = self.vib_cache[1]
        # return proper slice by index:
        if index is None:
            return coordOffset
        else:
            return coordOffset[index]
    def vib_adjustCoords(self, coords, index=None):
        """Take a coordinate array and add vibrations to it.
        """
        if getattr(self, 'vibEnabled', False):
            offset = self.vib_offset(index)
            return coords + offset
        else:
            return coords

    coordLookupCacheKey = None
    def getCCords(self, returnPointer=False):
        """Return a pointer to coordinate array suitable for use in ctypes.

        If `pointer` is True, return the <array>.ctypes.data object
        needed for passing to a C function.
        """
        cachepar = hash((id(self), self.mctime))

        if self.coordLookupCacheKey != cachepar:
            c = self.coords()
            # make sure we have proper shape/dtype
            if c.dtype != saiga12.numpy_double:
                c = numpy.asarray(c, dtype=saiga12.numpy_double)
            if not c.flags.carray:
                c = c.copy()
            # store it:
            self.coordLookupCacheKey = cachepar
            self.coordLookupCache = c
            self.coordLookupCache_p = c.ctypes.data

        if returnPointer:
            return self.coordLookupCache_p
        else:
            return self.coordLookupCache

