# Richard Darst, April 2008

import numpy
from common import *

class Grid2d(object):
    def makeconn_2Dgrid(self, a, b):
        self.lattShape = a, b
        lattSize = a*b
        self.lattSize = lattSize
        self.connMax = 4



        #self.conn=numpy.zeros(shape=(lattSize, self.connMax), dtype=numpy.int_)
        self.__dict__["conn"]=numpy.zeros(shape=(lattSize, self.connMax), dtype=numpy.int_)
        self.SD.conn = self.conn.ctypes.data
        self.conn.shape = lattSize, self.connMax

        #self.connN = numpy.zeros(shape=(lattSize), dtype=numpy.int_)
        self.__dict__["connN"] = numpy.zeros(shape=(lattSize), dtype=numpy.int_)
        self.SD.connN = self.connN.ctypes.data

        #self.lattsite = numpy.zeros(shape=(self.lattSize), dtype=numpy.int_)
        self.__dict__["lattsite"] = numpy.zeros(shape=(self.lattSize), dtype=numpy.int_)
        self.SD.lattsite = self.lattsite.ctypes.data
        self.lattsite[:] = S12_EMPTYSITE
        

        # x is our 2D array of maps of shapes
        x = numpy.arange(a*b)
        x.shape = a, b
        # celllist is a list of all possible (x,y) coordinates.
        celllist = [ ]
        for ai in range(a):
            for bi in range(b):
                celllist.append((ai, bi))
        celllist = numpy.asarray(celllist)
        neighborlist = numpy.asarray(
            (( 0, 1 ),
             ( 1, 0 ),
             ( 0,-1 ),
             (-1, 0 ),
             ))
        for c in celllist:
            for n in neighborlist:
                cur = x[tuple(c)]
                neighbor = x[tuple(numpy.mod(c+n, (a, b)))]
                i = self.connN[cur]
                self.conn[cur, i] = neighbor
                self.connN[cur] += 1
        #self.printLattice(x)
        #print self.conn
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


