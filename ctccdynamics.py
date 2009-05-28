# Richard Darst, May 2009

import numpy

class CTCCDynamics(object):

    def ctcc_adjustCoords(self, coords, index=None):
        activeParticles = self.orient != -1
        firstCoords = coords[activeParticles]
        neighCoords = \
           coords[self.conn[activeParticles, self.orient[activeParticles]]]
        newCoords = (3*firstCoords + neighCoords) / 4.
        return newCoords

        #from fitz.interact import interact ; interact()

