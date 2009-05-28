# Richard Darst, April 2008

import math
import numpy
from saiga12.common import *
from saiga12.util import cartesianproduct
import saiga12

# the following are from Richard's dimcalc.py file.
# this is the cryptic, impossible to understand version of these...
from operator import mul as operator_mul
product = lambda l: reduce(operator_mul, l, 1)
def coords(dims, index):
    """Given array dimensions and a (linear) position index,
    return the coordinates.
    """
    return tuple( (index // product(dims[i+1:])) % dims[i]
                  for i in range(len(dims)))
def index(dims, coords):
    """Given array dimensions and coordinates, return the (linear) index.
    """
    return reduce(lambda a, x:a*x[0]+x[1],
                  zip(dims, coords),
                  0)

# These globals are used to cache various dynamically-generated data.
# It is useful to do this, since you often have many copies of the
# exact same size system, so all this data is the same.  One warning:
# this means that modifying the data might modify it for ALL objects--
# right now this only applies to the coordinate data, but later be
# changed to apply to the connection data also.
grid_datacache = { }   # conn and connN
coord_datacache = { }  # coords


class GridNd(saiga12.Sys):
    
    def makegrid(self, *dimensions):
        lattSize = reduce(lambda x,y: x*y, dimensions) # product of numbers
        self.lattShape = dimensions
        self.physicalShape = self._makePhysicalShape(self.lattShape)
        # v-- used for re-initilizing the lattice.

        self._initArrays(lattSize=lattSize,
                         connMax=len(self._neighborlist))

        connCacheKey = (self.__class__.__name__, dimensions)
        if grid_datacache.has_key(connCacheKey):
            self.conn[:] = grid_datacache[connCacheKey][0]
            self.connN[:] = grid_datacache[connCacheKey][1]
            return
        # x is our 2D array of maps of shapes
        x = numpy.arange(lattSize)
        x.shape = dimensions
        #
        dimranges = [ range(i) for i in dimensions ]
        celllist = [ tuple(zz) for zz in cartesianproduct(*dimranges)]

        celllist = numpy.asarray(celllist)
        self._makegrid_convolve(celllist, x, dimensions)
        
        #self.printLattice(x)
        #print self.conn
        grid_datacache[connCacheKey] = self.conn, self.connN
    def _makegrid_convolve(self, celllist, x, dimensions):
        """Convolve the neighborlist across the lattice sites.

        This has been split into a separate function since some grid
        objects need to override it (for example, due to different
        handling of alternate rows).
        """
        for c in celllist:
            for n in self._neighborlist:
                cur = x[tuple(c)]
                neighbor = x[tuple(numpy.mod(c+n, dimensions))]
                i = self.connN[cur]
                self.conn[cur, i] = neighbor
                self.connN[cur] += 1
        
    def _makePhysicalShape(self, lattShape):
        """Return the physical size, used for periodic boundary conditions.

        lattShape happens to be equal to the periodicities in every
        dimension for rectangular cartesion grids, however, that isn't
        true for other grid arrangements.  This function should
        tranform a lattShape to the physical periodicity size in
        cartesian 3d space.

        For square rectangular grids, that is easy-- it is just the
        same lattShape.
        """
        return lattShape
    def distance(self, index0, index1):
        """Distance between any two lattice points.

        Depends on distance2 to work.
        """
        return numpy.sqrt(self.distance2(index0, index1))
    def distance2(self, index0, index1):
        """Distance-squared between any two lattice points.

        This works for arbitrary dimensions, as well as arrays!
        Depends on a functioning coords() method to work.
        """
        coords0 = self.coords(index0)
        coords1 = self.coords(index1)
        lindistances = coords0 - coords1
        # v-- this 
        delta = (numpy.floor(lindistances/self.physicalShape + .5)) * \
                self.physicalShape
        # Use the ufuncs so that all operations are done in-place.
        numpy.subtract(lindistances, delta, lindistances)
        numpy.multiply(lindistances, lindistances, lindistances)
        dists2 = numpy.sum(lindistances, axis=-1)
        return dists2
    def latticeReInitData(self):
        """Get state data needed to re-create our grid.

        Used in pickling.
        """
        return self.lattShape
    def latticeReInit(self, latticeReInitData):
        """Recreate our grid arrays using the data from above.

        Used when un-pickling.
        """
        self.makegrid(*latticeReInitData)



class SquareGrid(GridNd):
    """Class to make rectangular cartesion grid in any dimension.

    Must be subclassed and have connection data added to it.
    """
    def coords(self, index=None):
        """Mapping from lattice indexa to coordinates in real space.

        'index' is sort of misnamed here-- it should be 'pos' to be
        consistent with the rest of the code.  However, 'index'
        reminds you that it's simply an integer labeling, not really a
        *position*.

        This method is a replacement for grid_coords()

        if index is None, return all coordinates.
        """
        coordCacheKey = (self.__class__.__name__, self.lattShape)
        if not coord_datacache.has_key(coordCacheKey):
            c = numpy.asarray(coords(self.lattShape,
                                     numpy.arange(self.lattSize)),
                              dtype=saiga12.numpy_double)
            c = c.transpose()
            coord_datacache[coordCacheKey] = c.copy()
        if index is None:
            c = coord_datacache[coordCacheKey]
        else:
            index = numpy.asarray(index)
            c = coord_datacache[coordCacheKey][index]
        if getattr(self, 'vibEnabled', False):
            c = self.vib_adjustCoords(c)
        if self.orient is not None:
            c = self.ctcc_adjustCoords(c)
        return c

    def gridIndex(self, coords):
        """Mapping from coordinates in real space to lattice index
        """
        return index(self.lattShape, coords)

class Grid1d(SquareGrid):
    _neighborlist = numpy.asarray(
        ((-1, ),
         ( 1, ),
         ))
    def printLattice(self, lattice=None):
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        for e in lattice:
            if e == S12_EMPTYSITE:
                e = "__"
            print "%2s"%e,
    def printLatticeDots(self, lattice=None):
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        result = [ ]
        for e in lattice:
            if e == S12_EMPTYSITE:
                result.append("_")
            else:
                result.append("%s"%self.atomtype[e])
        print "".join(result)

class Grid2d(SquareGrid):
    _neighborlist = numpy.asarray(
        (( 0, 1 ),
         ( 1, 0 ),
         ( 0,-1 ),
         (-1, 0 ),
         ))

    def printLattice(self, lattice=None):
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        for row in lattice:
            for e in row:
                if e == S12_EMPTYSITE:
                    e = "__"
                print "%2s"%e,
            print
    def printLatticeDots(self, lattice=None):
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        result = [ ]
        for row in lattice:
            for e in row:
                if e == S12_EMPTYSITE:
                    result.append("_")
                else:
                    result.append("%s"%self.atomtype[e])
            result.append("\n")
        print "".join(result[:-1]) # exclude trailing newline

    def printLatticeLocal(self, pos, width=2, lattice=None):
        center = pos//self.lattShape[1], pos%self.lattShape[1]
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        # iterate over rows
        #lattice = lattice[center[0]-width:center[0]+width+1,
        #                   center[1]-width:center[1]+width+1]
        print "lattice centered on pos %s %s:"%(pos, center)
        for rown in range(center[0]-width, center[0]+width+1):
            for coln in range(center[1]-width, center[1]+width+1):
                #print "xxx", rown, coln
                e = lattice[rown%self.lattShape[0],
                            coln%self.lattShape[1]]
                if e == S12_EMPTYSITE:
                    e = "__"
                print "%2s"%e,
            print


class Grid3d(SquareGrid):
    _neighborlist = numpy.asarray(
        (( 0, 1, 0 ),
         ( 1, 0, 0 ),
         ( 0, 0, 1 ),
         ( 0,-1, 0 ),
         (-1, 0, 0 ),
         ( 0, 0,-1 ),
         ))

class GridHex2d(GridNd):
    """Hexagonal two-dimensional grid.

    Y dimension must be multiples of two.  To make the grid, call
    .makegrid(x, y).
    """
    _neighborlist = [ None ] * 6   
    
    def _makePhysicalShape(self, lattShape):
        return numpy.multiply(lattShape,
                              (1, math.sqrt(3)/2))
    def _makegrid_convolve(self, celllist, x, dimensions):
        if dimensions[1]%2 != 0:
            raise Exception("Second dimension must be multiple of two.")
        neighborlist_Orig = ( ( 1,  0 ), (-1,  0 ),
                              ( 0,  1 ), ( 0, -1 ), )
        neighborlist_evenY = ( (-1,  1 ),
                               (-1, -1 ), )
        neighborlist_oddY =  ( ( 1,  1 ),
                               ( 1, -1 ), )

        for c in celllist:
            y = c[1]
            neighborlist = list(neighborlist_Orig) # make a copy
            # We need every other y-row to be different:
            if y%2 == 0:
                neighborlist += neighborlist_evenY
            else:
                neighborlist += neighborlist_oddY
            neighborlist = neighborlist[0::2] + neighborlist[1::2]
            for n in neighborlist:
                cur = x[tuple(c)]
                neighbor = x[tuple(numpy.mod(c+n, dimensions))]
                i = self.connN[cur]
                self.conn[cur, i] = neighbor
                self.connN[cur] += 1
    def coords(self, index=None):
        """Mapping from lattice index to coordinates in real space.

        'index' is sort of misnamed here-- it should be 'pos' to be
        consistent with the rest of the code.  However, 'index'
        reminds you that it's simply an integer labeling, not really a
        *position*.

        This method is a replacement for grid_coords()
        """
        coordCacheKey = (self.__class__.__name__, self.lattShape)
        if not coord_datacache.has_key(coordCacheKey):
            c = numpy.asarray(coords(self.lattShape,
                                     numpy.arange(self.lattSize)),
                              dtype=saiga12.numpy_double)
            hexmat = numpy.asarray(((  1,        0         ),
                                    (  0,  math.sqrt(3)/2) )  )
            c = numpy.dot(hexmat, c)
            c = c.transpose()
            c[1::2,0] +=  .5
            coord_datacache[coordCacheKey] = c.copy()
        if index is None:
            c = coord_datacache[coordCacheKey]
        else:
            index = numpy.asarray(index)
            c = coord_datacache[coordCacheKey][index]
        if getattr(self, 'vibEnabled', False):
            c = self.vib_adjustCoords(c)
        return c

    # The print lattice methods here are the same as for the 2d grid--
    # they should probably be unified sometime.
    def printLattice(self, lattice=None):
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        for row in lattice:
            for e in row:
                if e == S12_EMPTYSITE:
                    e = "__"
                print "%2s"%e,
            print

    def printLatticeLocal(self, pos, width=2, lattice=None):
        center = pos//self.lattShape[1], pos%self.lattShape[1]
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        # iterate over rows
        #lattice = lattice[center[0]-width:center[0]+width+1,
        #                   center[1]-width:center[1]+width+1]
        print "lattice centered on pos %s %s:"%(pos, center)
        for rown in range(center[0]-width, center[0]+width+1):
            for coln in range(center[1]-width, center[1]+width+1):
                #print "xxx", rown, coln
                e = lattice[rown%self.lattShape[0],
                            coln%self.lattShape[1]]
                if e == S12_EMPTYSITE:
                    e = "__"
                print "%2s"%e,
            print


class Grid3dHCP(GridNd):
    """Hexagonal close-packed grid.

    Y and Z dimensions must be multiples of two.  To make the grid,
    call .makegrid(x, y, z).
    """

    # this must have the right length, but isn't used here:
    _neighborlist = [ None ] * 12

    def _makePhysicalShape(self, lattShape):
        """Set the periodic box size.
        """
        return numpy.multiply(lattShape,
                              (1, math.sqrt(3)/2, math.sqrt(6)/3))

    def _makegrid_convolve(self, celllist, x, dimensions):
        if dimensions[1]%2 != 0 or dimensions[2]%2 != 0:
            raise Exception("Second and third dims must be multiples of two.")
        neighborlist_Orig = ( ( 1,  0,  0), (-1,  0,  0), # a
                              ( 0,  1,  0), ( 0, -1,  0), # b
                              ( 0,  0,  1), ( 0,  0, -1), # d
                              )

        neighborlist_evenY = ( (-1,  1,  0),  # c
                               (-1, -1,  0),)
        neighborlist_oddY =  ( (+1,  1,  0),  # c
                               (+1, -1,  0),)

        # These neighborlists can be found by realizing that every
        # even/odd y/z layer is different, drawing out the actual
        # configuration (since you did so much work checking the
        # coordinates function, you now know where any site `i` is
        # placed).  Then, go through and very carefully figure out
        # what is next to every site.  Use pictures.
        neighborlist_evenZevenY = ( (-1, +1, +1), (-1,  0, +1), # e
                                    (-1, +1, -1), (-1,  0, -1),)# f
        neighborlist_oddZevenY  = ( (+1,  0, -1), ( 0, -1, -1), # e
                                    (+1,  0, +1), ( 0, -1, +1),)# f
        neighborlist_evenZoddY  = ( ( 0, +1, +1), (-1,  0, +1), # e
                                    ( 0, +1, -1), (-1,  0, -1),)# f
        neighborlist_oddZoddY   = ( (+1,  0, -1), (+1, -1, -1), # e
                                    (+1,  0, +1), (+1, -1, +1),)# f
        for c in celllist:
            y = c[1]
            z = c[2]
            neighborlist = list(neighborlist_Orig) # make a copy
            if y%2 == 0:
                neighborlist += neighborlist_evenY
            else:
                neighborlist += neighborlist_oddY

            if z%2 == 1 and y%2==0:   # oddZ evenY
                neighborlist += neighborlist_oddZevenY
            elif z%2 == 0 and y%2==1: # evenZ oddY
                neighborlist += neighborlist_evenZoddY
            elif z%2 == 0:            # evenZ evenY
                neighborlist += neighborlist_evenZevenY
            else:                     # oddZ oddY
                neighborlist += neighborlist_oddZoddY
            neighborlist = neighborlist[0::2] + neighborlist[1::2]
            for n in neighborlist:
                cur = x[tuple(c)]
                neighbor = x[tuple(numpy.mod(c+n, dimensions))]
                i = self.connN[cur]
                self.conn[cur, i] = neighbor
                self.connN[cur] += 1
    def coords(self, index=None):
        """Mapping from lattice indexa to coordinates in real space.

        'index' is sort of misnamed here-- it should be 'pos' to be
        consistent with the rest of the code.  However, 'index'
        reminds you that it's simply an integer labeling, not really a
        *position*.

        This method is a replacement for grid_coords()
        """
        coordCacheKey = (self.__class__.__name__, self.lattShape)
        if not coord_datacache.has_key(coordCacheKey):
            c = numpy.asarray(coords(self.lattShape,
                                     numpy.arange(self.lattSize)),
                              dtype=saiga12.numpy_double)

            # Linear transform -- just shrink the coordinates down.
            # Offsets still need to be adjusted below.
            hexmat = numpy.asarray(((1,              0  ,  0                ),
                                    (0,  math.sqrt(3)/2.,  0                ),
                                    (0,                0,  math.sqrt(6)/3 )))
            c = numpy.dot(hexmat, c)
            c = c.transpose()
            lattShape = self.lattShape
            c.shape = lattShape[0], lattShape[1], lattShape[2], 3

            # adjust x offset by .5 for every other y layer
            c[:, 1::2, :, 0] +=  .5
            # adjust z offset by the proper amount for ever other z layer
            c[:, :, 1::2, :2] +=  (.5, -math.sqrt(3)/6)

            c.shape = product(lattShape), 3
            self._coordCache = c.copy()
            coord_datacache[coordCacheKey] = c.copy()
        if index is None:
            c = coord_datacache[coordCacheKey]
        else:
            index = numpy.asarray(index)
            c = coord_datacache[coordCacheKey][index]
        if getattr(self, 'vibEnabled', False):
            c = self.vib_adjustCoords(c)
        return c

    def printLattice(self, lattice=None):
        # expanded for 3d.
        if lattice is None:
            lattice = self.lattsite.reshape(self.lattShape)
        for plane in lattice:
            for row in plane:
                for e in row:
                    if e == S12_EMPTYSITE:
                        e = "__"
                    print "%2s"%e,
                print
            print "----------------------------------------------"
