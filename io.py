# Richard Darst, April 2008

#import shelve
import cPickle as pickle

classvars = ("N", "beta", "hardness", "mctime")
            # "lattSize", "lattShape", "latticeReInitData"
arrays = ("lattsite", "nneighbors", "atomtype", "atompos", "ntype")
            # self.conn, self.connN, self.connMax (not array),

class IOSys(object):
    def io_state(self, otherData=None):
        """Return a dict containing enough state to reconstruct the system.
        """
        state = { }
        for key in classvars:
            state[key] = getattr(self, key)
        for key in arrays:
            state[key] = getattr(self, key)
        state["latticeReInitData"] = self.latticeReInitData()
        state["stateSaveVersion"] = 1
        state["otherData"] = otherData
        return state
    def io_writeToFile(self, filename, otherData=None):
        """Write state to a file.
        """
        state = self.io_state(otherData)
        state.update(otherData)
        pickle.dump(state, file(filename, "wb"), protocol=-1)
    def io_loadFromFile(self, filename):
        """Load state from a file and load onto self
        """
        state = pickle.load(file(filename))
        self.io_loadState(state)
    def io_loadState(self, state):
        """Set self's state based on passed dict of state.
        """
        self.latticeReInit(state["latticeReInitData"])
        for key in classvars:
            if not state.has_key(key):
                continue
            setattr(self, key, state[key])
        for key in arrays:
            getattr(self, key)[:] = state[key]
        self.consistencyCheck()
        
        
