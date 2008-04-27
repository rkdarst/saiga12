# Richard Darst, April 2008

#import shelve
import cPickle as pickle

classvars = ("N", "beta", "hardness", "mctime",
             "otherData")
            # "lattSize", "lattShape", "latticeReInitData"
arrays = ("lattsite", "nneighbors", "atomtype", "atompos", "ntype")
            # self.conn, self.connN, self.connMax (not array),

consistencyCheck = False

class IOSys(object):
    def io_state(self, otherData=None):
        """Return a dict containing enough state to reconstruct the system.
        """
        state = { }
        for key in classvars:
            if hasattr(self, key):
                state[key] = getattr(self, key)
        for key in arrays:
            state[key] = getattr(self, key)
        state["latticeReInitData"] = self.latticeReInitData()
        state["stateSaveVersion"] = 1
        if otherData:
            state["otherData"] = otherData
        return state
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
        if state.has_key("otherData"):
            self.otherData = state["otherData"]
        if consistencyCheck:
            self.consistencyCheck()
    def __getstate__(self):
        return self.io_state()
    def __setstate__(self, state):
        self.__init__()
        self.io_loadState(state)


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


def io_open(state):
    """Open state from a pickle and return it.

    If it's a file-object, the next pickle object from the stream is
    read and processed as below.

    Return the Sys object.
    """
    if type(state) == file:
        state = pickle.load(state)

    if not type(state) == dict:
        return state
    # we have a dictionary, manually reconstruct it...
    from saiga12.geom.grid import Grid3d
    S = Grid3d()
    S.io_loadState(state)
    return S

if __name__ == "__main__":
    import code
    import readline
    import sys

    S = io_open(file(sys.argv[1]))
    code.interact(local=locals(), banner="")
    
