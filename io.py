# Richard Darst, April 2008

#import shelve
import cPickle as pickle

classvars = ("N", "beta", "mctime",
             "otherData")
            # "lattSize", "lattShape", "latticeReInitData"
arrays = ("lattsite", "nneighbors", "atomtype", "atompos", "ntype")
            # self.conn, self.connN, self.connMax (not array),
specialvars = ("hardness", "latticeReInitData", "stateSaveVersion",
               "otherData")

consistencyCheck = False

class IOSys(object):
    def io_state(self):
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

        # python2.4 can't pickle float("inf"), but 2.5 can.
        if self.hardness == float("inf"):
            state["hardness"] = "inf"
        else:
            state["hardness"] = self.hardness

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

        if state["hardness"] == "inf":
            self.hardness == float("inf")
        else:
            self.hardness == state["hardness"]

        if consistencyCheck:
            self.consistencyCheck()
    def __getstate__(self):
        return self.io_state()
    def __setstate__(self, state):
        self.__init__()
        self.io_loadState(state)


    #def io_writeToFile(self, filename):
    #    """Write state to a file.
    #
    #    This function has problems, 
    #    """
    #    pickle.dump(self, file(filename, "wb"), protocol=-1)
    #def io_loadFromFile(self, filename):
    #    """Load state from a file and load onto self
    #    """
    #    state = pickle.load(file(filename))
    #    self.io_loadState(state)


def io_open(state):
    """Open state from a pickle and return it.

    If it's a file-object, the next pickle object from the stream is
    read and processed as below.

    Return the Sys object.
    """
    # if it is a file-object, load from there.
    if hasattr(state, "read") and hasattr(state, "readline"):
        state = pickle.load(state)
    # if we we just loaded a dictionary, manually reconstruct it...
    if type(state) == dict:
        from saiga12.geom.grid import Grid3d
        S = Grid3d()
        S.io_loadState(state)
        return S

    return state


if __name__ == "__main__":
    import code
    import sys
    try: import readline
    except ImportError: pass
    else:
        import rlcompleter
        readline.parse_and_bind("tab: complete")

    try:
        from urllib import urlopen as open
    except ImportError:
        pass
    
    S = io_open(open(sys.argv[1]))
    code.interact(local=locals(), banner="")
    
