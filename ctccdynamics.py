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

        if index is None:
            activeParticles = self.orient != -1
        else:
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

        coords = coords.copy() # we *need* to make a copy here...
        coords[activeParticles] += deltas
        #from fitz.interact import interact ; interact()
        return coords


