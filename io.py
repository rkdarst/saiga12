# Richard Darst, April 2008

#import shelve
import cPickle as pickle

class IOSys(object):
    def io_state(self):
        state = {
            "N": self.N,
            "beta": self.beta,
            "hardness": self.hardness,
            "lattSize": self.lattSize,
            "lattShape": self.lattShape,
            "latticeReInitData": self.latticeReInitData,

            # arrays
            "lattsite": self.lattsite,
            #"conn": self.conn,
            #"connN": self.connN,
            #"connMax": self.connMax,
            "nneighbors": self.nneighbors,
            "atomtype": self.atomtype,
            "atompos": self.atompos,
            "ntype": self.ntype,
            
            }
        return state
    def io_writeToFile(self, filename):
        state = self.io_state()
        fo = file(filename, "wb")
        pickle.dump(state, fo, protocol=-1)
