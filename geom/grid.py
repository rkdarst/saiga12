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
def cords(dims, index):
    """Given array dimensions and a (linear) position index,
    return the coordinates.
    """
    return tuple( (index // product(dims[i+1:])) % dims[i]
                  for i in range(len(dims)))
def index(dims, cords):
    """Given array dimensions and coordinates, return the (linear) index.
    """
    return reduce(lambda a, x:a*x[0]+x[1],
                  zip(dims, cords),
                  0)


class GridNd(saiga12.Sys):
    
    def makegrid(self, *dimensions):
        lattSize = reduce(lambda x,y: x*y, dimensions) # product of numbers
        self.lattShape = dimensions
        # v-- used for re-initilizing the lattice.
        self.latticeReInitData = dimensions

        self._initArrays(lattSize=lattSize,
                         connMax=len(self._neighborlist))

        # x is our 2D array of maps of shapes
        x = numpy.arange(lattSize)
        x.shape = dimensions
        # celllist is a list of all possible (x,y) coordinates.
        #celllist = [ ]
        #self.fillcelllist(celllist, dimensions)
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
    latticeReInit = makegrid

class Grid2d(GridNd):
    _neighborlist = numpy.asarray(
        (( 0, 1 ),
         ( 1, 0 ),
         ( 0,-1 ),
         (-1, 0 ),
         ))
    #def fillcelllist(self, celllist, dimensions):
    #    a = dimensions[0]
    #    b = dimensions[1]
    #    for ai in range(a):
    #        for bi in range(b):
    #            celllist.append((ai, bi))

    def distance(self, index0, index1):
        """Distance between any two lattice points

        This works for arbitrary dimensions, as well as arrays!
        """
        cords0 = numpy.asarray(cords(self.lattShape, index0), dtype=float)
        cords1 = numpy.asarray(cords(self.lattShape, index1), dtype=float)
        lindistances = cords0 - cords1
        dists = numpy.sqrt(numpy.sum(lindistances * lindistances, axis=0))
        return dists


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


class Grid3d(GridNd):
    _neighborlist = numpy.asarray(
        (( 0, 1, 0 ),
         ( 1, 0, 0 ),
         ( 0,-1, 0 ),
         (-1, 0, 0 ),
         ( 0, 0, 1 ),
         ( 0, 0,-1 ),
         ))
    #def fillcelllist(self, celllist, dimensions):
    #    a = dimensions[0]
    #    b = dimensions[1]
    #    c = dimensions[1]
    #    for ai in range(a):
    #        for bi in range(b):
    #            for ci in range(c):
    #                celllist.append((ai, bi, ci))
