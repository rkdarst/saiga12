# Richard Darst, April 2008

import numpy
from saiga12.common import *
import saiga12

# the below is from
#  http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/302478
# it basically defines a cartesian product.
# I think the activestate sample recipees are freely-usable and distributable
def xcombine(*seqin):
    '''returns a generator which returns combinations of argument sequences
    for example xcombine((1,2),(3,4)) returns a generator; calling the next()
    method on the generator will return [1,3], [1,4], [2,3], [2,4] and
    StopIteration exception.  This will not create the whole list of
    combinations in memory at once.'''
    def rloop(seqin,comb):
        '''recursive looping function'''
        if seqin:                   # any more sequences to process?
            for item in seqin[0]:
                newcomb=comb+[item]     # add next item to current combination
                # call rloop w/ remaining seqs, newcomb
                for item in rloop(seqin[1:],newcomb):
                    yield item          # seqs and newcomb
        else:                           # processing last sequence
            yield comb                  # comb finished, add to list
    return rloop(seqin,[])
cartesianproduct = xcombine

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

grid_datacache = { }


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
        for c in celllist:
            for n in self._neighborlist:
                cur = x[tuple(c)]
                neighbor = x[tuple(numpy.mod(c+n, dimensions))]
                i = self.connN[cur]
                self.conn[cur, i] = neighbor
                self.connN[cur] += 1
        #self.printLattice(x)
        #print self.conn
        grid_datacache[connCacheKey] = self.conn, self.connN
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
        """
        return numpy.sqrt(self.distance2(index0, index1))
    def distance2(self, index0, index1):
        """Distance-squared between any two lattice points.

        This works for arbitrary dimensions, as well as arrays!
        """
        coords0 = self.coords(index0)
        coords1 = self.coords(index1)
        lindistances = coords0 - coords1
        # v-- this 
        delta = (numpy.floor(lindistances/self.physicalShape + .5)) * \
                self.physicalShape
        lindistances = lindistances - delta
        dists2 = numpy.sum(lindistances * lindistances, axis=-1)
        return dists2
    def latticeReInitData(self):
        """Get state data needed to re-create our grid.
        """
        return self.lattShape
    def latticeReInit(self, latticeReInitData):
        """Recreate our grid arrays using the data from above.
        """
        self.makegrid(*latticeReInitData)



class SquareGrid(GridNd):
    """Class to make rectangular cartesion grid in any dimension.

    Must be subclassed and have connection data added to it.
    """
    def coords(self, index):
        """Mapping from lattice indexa to coordinates in real space.

        'index' is sort of misnamed here-- it should be 'pos' to be
        consistent with the rest of the code.  However, 'index'
        reminds you that it's simply an integer labeling, not really a
        *position*.

        This method is a replacement for grid_coords()
        """
        c = numpy.asarray(coords(self.lattShape, index),
                             dtype=saiga12.numpy_double)
        c = c.transpose()
        # should I really copy this ?
        if not c.flags.c_contiguous:
            c = c.copy()
        return c
    def gridIndex(self, coords):
        """Mapping from coordinates in real space to lattice index
        """
        return index(self.lattShape, coords)

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
         ( 0,-1, 0 ),
         (-1, 0, 0 ),
         ( 0, 0, 1 ),
         ( 0, 0,-1 ),
         ))
