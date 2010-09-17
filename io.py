# Richard Darst, April 2008

#import shelve
import cPickle as pickle
import zlib
import numpy

classvars = ("N", "beta", "mctime",
             "movesPerCycle", "cumProbAdd", "cumProbDel",
             "otherData", "cycleModeStr", "energyModeStr",
             "inserttype", "_dontSetCycleMoves",
             "vibEnabled", "flags")
            # "lattSize", "lattShape", "latticeReInitData"
arrays = ("lattsite", "nneighbors", "atomtype", "atompos", "ntype")
            # self.conn, self.connN, self.connMax (not array),
specialvars = ("hardness", "latticeReInitData", "stateSaveVersion",
               "otherData")

consistencyCheck = False

class IOSys(object):
    ioSaveVersion = 2
    def io_state(self):
        """Return a dict containing enough state to reconstruct the system.
        """
        state = { }
        for key in classvars:
            if hasattr(self, key):
                state[key] = getattr(self, key)
        state["latticeReInitData"] = self.latticeReInitData()
        # python2.4 can't pickle float("inf"), but 2.5 can.
        if self.hardness == float("inf"):
            state["hardness"] = "inf"
        else:
            state["hardness"] = self.hardness
        # record history of event-driven dynamics (this won't
        # automatically enable it, though
        if hasattr(self, "_eddEnabled"):
            state['_eddWasEnabled'] = self._eddEnabled
        # Saving arrays is different for different methods:
        state["stateSaveVersion"] = self.ioSaveVersion
        if self.ioSaveVersion == 1:
            # save full arrays
            for key in arrays:
                state[key] = getattr(self, key)
        elif self.ioSaveVersion == 2:
            #print "saving with version 2"
            # be more clever about the arrays we save
            import numpy
            atompos = self.atompos[:self.N]
            if len(atompos) > 0 and numpy.max(atompos) < 32768:
                atompos = numpy.asarray(atompos, dtype=numpy.int16)
            state["atompos"] = atompos
            state["atomtype"] = numpy.asarray(self.atomtype[:self.N],
                                              dtype=numpy.int8)
        else:
            raise Exception("Invalid save version when saving: %s",
                            self.ioSaveVersion)
        # way of saving persistence arrays
        if self.persist is not None:
            persist = (0,  # save version
                       zlib.compress(pickle.dumps(    # compressed array
                         numpy.asarray(self.persist, dtype=numpy.bool_),
                         -1), 6)
                       )
            state["persist"] = persist
        if self.orient is not None:
            orientversion = 0
            orient = numpy.asarray(self.orient[self.atompos], dtype=numpy.int8)
            state["orient"] = orientversion, orient
        if hasattr(self, "_currentInsertTypes"):
            state['insertTypes'] = self._currentInsertTypes


        return state
    def ioSave(self, fname):
        """Equivalent to pickle.dump(S, file(fname, 'wb'), -1)
        """
        pickle.dump(self, file(fname, "wb"), pickle.HIGHEST_PROTOCOL)
    def io_loadState(self, state):
        """Set self's state based on passed dict of state.
        """
        self.latticeReInit(state["latticeReInitData"])
        if state.has_key("cycleModeStr"):
            self.setCycleMode(state["cycleModeStr"])
            del state["cycleModeStr"]
        if state.has_key("energyModeStr"):
            self.setEnergyMode(state["energyModeStr"])
            del state["energyModeStr"]
        for key in classvars:
            if not state.has_key(key):
                continue
            setattr(self, key, state[key])
        if state["hardness"] == "inf":
            self.hardness = float("inf")
        else:
            self.hardness = state["hardness"]
        if 'insertTypes' in state:
            self.setInsertType(state['insertTypes'])
            del state['insertTypes']
        # was event-driven dynamics enabled before?  (note: you *must*
        # regenerate the move-lists, you can't just set the variable
        # and call it enabled)
        if state.has_key('_eddWasEnabled'):
            self._eddWasEnabled = state['_eddWasEnabled']
        if state['stateSaveVersion'] == 1:
            for key in arrays:
                getattr(self, key)[:] = state[key]
        elif state['stateSaveVersion'] == 2:
            self.atompos[:self.N]  = state['atompos']
            self.atomtype[:self.N] = state['atomtype']
            self.C.loadStateFromSave(self.SD_p)
        else:
            raise Exception("Invalid save version when loading: %s"%
                            state['stateSaveVersion'])
        if self.cycleModeStr in ("fredricksonandersen", "east", "spiral"):
            self.eddEnable()
        if state.has_key("persist"):
            version, persist = state["persist"]
            if version == 0:
                if self.persist is None: self._allocPersistArray()
                persist = zlib.decompress(persist)
                persist = pickle.loads(persist)
                self.persist[:] = persist
            else:  # unknown version
                if not getattr(self, "_ignorePersistLoading"):
                    raise Exception("Invalid persist array version: %s"%
                                    version)
        if "orient" in state:
            orientversion, orient = state["orient"]
            if orientversion != 0:
                raise Exception("Unknown orient save version")
            self.orient[self.atompos] = orient
        if state.has_key("vibEnabled"):
            self.vib_init()
        

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
    if isinstance(state, str):
        state = file(state, "rb")
    # if it is a file-object, load from there.
    if hasattr(state, "read") and hasattr(state, "readline"):
        # This code is because of a symlink 'saiga12' -> '.' that
        # existed in my source directories.  It caused some things to
        # be saved with the class name
        # saiga12.saiga12.geom.grid.Grid3d.
        # This class name can't be imported without that symlink... 
        p = pickle.Unpickler(state)
        def find_global(module, className):
            if 'saiga12.geom.grid' in module:
                import saiga12.geom.grid
                return getattr(saiga12.geom.grid, className)
            #exec "import %s"%module in {}, {}
            mod = __import__(module, globals(), locals(), 'x')
            return getattr(mod, className)
            raise
        p.find_global = find_global
        state = p.load()
        #state = pickle.load(state)
    # if we we just unpickled a dictionary, manually reconstruct it...
    if type(state) == dict:
        from saiga12.geom.grid import Grid3d
        S = Grid3d()
        S.io_loadState(state)
        return S

    return state
open = io_open # too bad this is a namespace collision, but hopefully
                # it's ok here...

if __name__ == "__main__":
    import code
    import os
    import sys
    try: import readline
    except ImportError: pass
    else:
        import rlcompleter
        readline.parse_and_bind("tab: complete")

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-e', '--exec', dest='exec_', default=None)
    options, args = parser.parse_args()

    try:
        from urllib import urlopen as open
    except ImportError:
        pass
    
    if options.exec_:
        for fname in args:
            if not os.access(fname, os.F_OK):
                print "Not accessible:", fname
                continue
            S = io_open(open(fname))
            exec(options.exec_)
    else:
        S = io_open(open(args[0]))
        code.interact(local=locals(), banner="")
    
