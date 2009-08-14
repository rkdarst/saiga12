# Richard Darst, May 2009

import numpy

class CTCCDynamics(object):

    def ctcc_adjustCoords(self, coords, coordsRaw, index=None):
        """Adjust coords based on CTCC dynamics orientation.
        """
        #activeParticles = self.orient != -1
        #firstCoords = coords[activeParticles]
        #neighCoords = \
        #   coords[self.conn[activeParticles, self.orient[activeParticles]]]
        #newCoords = (3*firstCoords + neighCoords) / 4.
        #coords0 = coords.copy()
        #coords0[activeParticles] = newCoords
        #return coords0

        coords = coords.copy() # we *need* to make a copy here...

        if index is None:
            activeParticles = self.orient != -1
        else:
            index = numpy.asarray(index)
            # We need to handle the case where a single integer is
            # passed in as an index.  This makes slicing harder, and
            # must be expanded to the right number of dimensions to be
            # sliced.
            if index.shape == ():
                index.shape = 1,
                coords.shape = 1,3
            activeParticles = self.orient[index] != -1
            active_FromAll = index[self.orient[index] != -1]
        firstCoords = coords[activeParticles]
        if index is None:
            neighCoords = \
            coordsRaw[self.conn[activeParticles, self.orient[activeParticles]]]
        else:
            neighCoords = \
              coordsRaw[self.conn[active_FromAll, self.orient[active_FromAll]]]

        deltas = (neighCoords - firstCoords)
        physicalShape = numpy.asarray(self.physicalShape)
        deltas += .5 * physicalShape
        numpy.mod(deltas, physicalShape, deltas)
        deltas -= .5 * physicalShape
        deltas /= 4.

        coords[activeParticles] += deltas
        # If our input `index` was an integer, then we need to be sure
        # the right dimenionality of array comes out.  See `.coords()`
        # method in geom/grid.py for a comment describing this.
        if coords.shape == (1, 3):
            coords.shape = 3,
        return coords


