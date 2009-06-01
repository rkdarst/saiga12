# Richard Darst, May 2009

import numpy

class CTCCDynamics(object):

    def ctcc_adjustCoords(self, coords, index=None):
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


        activeParticles = self.orient != -1
        firstCoords = coords[activeParticles]
        neighCoords = \
           coords[self.conn[activeParticles, self.orient[activeParticles]]]
        #import fitz.interactnow
        deltas = neighCoords - firstCoords

        physicalShape = numpy.asarray(self.physicalShape)
        deltas += .5 * physicalShape
        numpy.mod(deltas, physicalShape, deltas)
        deltas -= .5 * physicalShape

        deltas /= 4.
        coords = coords.copy() # we *need* to make a copy here...
        coords[activeParticles] += deltas

        #coords[coords==14.25] = 10.5
        #coords[coords==-.25] = 3.5

        #if not numpy.all(coords == coords0):
        #    import fitz.interactnow

        #from fitz.interact import interact ; interact()
        return coords


