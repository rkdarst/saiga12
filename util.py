# Richard Darst, September 2008

import math
import numpy
import sys

try:
    import visual
    import saiga12.viz as viz
except ImportError:
    pass


import saiga12
from saiga12.io import io_open

def getch():
    import tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch
def vector(x):
    if len(x) == 1:
        return x[0], 0, 0
    return x


def getNewFrameIndex(frame_index, nFrames, V=None, otherObjects=()):
    #ch = getch()
    inLoop = True
    while inLoop:
        inLoop = False
        ch = '~'
        while ch not in '>.<,09xqcCtnvbI':
            ch = visual.scene.kb.getkey()
            
        if ch in '>.': frame_index += 1
        if ch in '<,': frame_index -= 1
        if ch == '0': frame_index = 0
        if ch == '9': frame_index = nFrames-1
        if frame_index < 0: frame_index = 0
        if frame_index > nFrames-1: frame_index = nFrames-1
        if ch in 'xq': return None
        if ch == 'c': 
            if visual.scene.mouse.pick != None:
                visual.scene.center = visual.scene.mouse.pick.pos
            inLoop = True
        if ch == 'C':
            if hasattr(visual.scene, "originalCenter"):
                visual.scene.center = visual.scene.originalCenter
            inLoop = True
        if ch == 't': viz.tagToggle(); inLoop = True
        if ch == 'n': viz.tagToggle2(otherObjects); inLoop = True
        if ch == 'v': viz.toggleViz(V,otherObjects=otherObjects); inLoop = True
        if ch == 'b': viz.toggleBG(); inLoop = True
        if ch == 'I': from code import interact ; interact(local=locals())
    
    return frame_index

class Averager(object):
    """Numerically Stable Averager

    Calculate averages and standard deviations in a numerically stable way.

    From the 'On-Line Algorithm' from
    http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    The optional initialization argument datatype= should be a
    generator function which returns zeros of the type being averaged.
    For example, to average numpy arrays of ten values, use:
      Averager(datatype=lambda: numpy.zeros(10))
    """
    def __init__(self, datatype=float):
        self.n       = 0
        self._mean   = datatype()   # mean
        self._M2     = datatype()   # variance accumulator
        #self._var_n  = 0.
        #self._var_n1 = 0.
    def add(self, value):
        """Add a new number to the dataset.
        """
        n = self.n + 1
        delta = value - self._mean
        mean = self._mean + delta/n
        M2 = self._M2 + delta*(value - mean)
        self.n = n
        self._mean = mean
        self._M2 = M2
        return self
    __iadd__ = add

    @property
    def mean(self):
        """Mean"""
        return self._mean
    @property
    def std(self):
        """Population Variance"""
        return math.sqrt(self._M2 / self.n)
    @property
    def stdsample(self):
        """Sample Variance"""
        return math.sqrt(self._M2 / (self.n-1))
    @property
    def var(self):
        """Population Standard Deviation"""
        return self._M2 / self.n
    @property
    def varsample(self):
        """Sample Standard Deviation"""
        return self._M2 / (self.n-1)

def cartesianproduct(*args):
    """Cartesion product of iterable arguments.

    This implementation is thanks to MrGreen.
    """
    if len(args) == 0:
        yield ()
    else:
        for x in args[0]:
            for xs in cartesianproduct(*args[1:]):
                yield (x,) + xs

def diff(fname0, *fileNames):
    frame0 = io_open(fname0)

    if frame0.cycleModeStr in ('fredricksonandersen', ):
        difftype = 'spin'
    else:
        difftype = 'particle'
    frame_index = 0
    V = viz.VizSystem(frame0)
    V.vizMakeBox()
    objs = [ ]
    print frame0.lattSize
    while True:
        print '\r',
        print fileNames[frame_index][-20:],
        frame = io_open(fileNames[frame_index])
        #print frame, frame0
        if frame0.lattShape != frame.lattShape:
            raise "different lattice shapes!"
        for i in range(len(objs)):
            objs[0].visible = 0
            del objs[0]
        moveTypes = { }
        nMove = 0
        nMoved = 0
        moveDistance = 0.
        moveDistances = [ ]
        usePersist = False


        if difftype == 'particle':
            if frame.persist is not None  and  usePersist:
                for i in range(frame0.lattSize):
                    if frame.persist[i] != 0 and \
                           frame.lattsite[i] == saiga12.S12_EMPTYSITE:
                        objs.append(visual.sphere(
                            pos=vector(frame0.coords(i)),
                            radius=.25,
                            color = visual.color.red
                            ))
            for n in range(frame0.N):
                if frame0.atompos[n] == frame.atompos[n]:
                    continue
                type_ = frame0.atomtype[n]
                moveTypes[type_] = moveTypes.get(type_, 0) + 1
                # old position - red
                startCoords = frame0.coords(frame0.atompos[n])
                if frame.lattsite[frame0.atompos[n]] == saiga12.S12_EMPTYSITE:
                    objs.append(visual.sphere(
                        pos=vector(startCoords),
                        radius=.25,
                        color = visual.color.red
                        ))
                # new position: whatever
                endCoords = frame0.coords(frame.atompos[n])
                objs.append(visual.sphere(
                    pos=endCoords,
                    radius=.25,
                    color = V.vizColors.get(frame0.atomtype[n],
                                            visual.color.white)
                    ))
                d = frame0.distance(frame0.atompos[n], frame.atompos[n])
                d_raw = endCoords - startCoords
                d_raw = math.sqrt(numpy.sum(d_raw*d_raw))
                if d_raw < 7.5:
                 objs.append(
                    visual.arrow(pos=startCoords, axis=endCoords-startCoords,
                                 shaftwidth=.1, headlength=1, fixedwidth=1
                                 ))
                moveDistance += d
                moveDistances.append(d)
            print moveTypes, "%.4f"%numpy.mean(moveDistances), \
                  "%.4f"%numpy.std(moveDistances), \
                  "%.4f"%(len(moveDistances)/frame0.lattSize),
        if difftype == 'spin':
            for i in range(frame0.lattSize):
                if frame.persist[i] != 0:
                    objs.append(visual.sphere(
                        pos=vector(frame0.coords(i)),
                         radius=1,
                        color = visual.color.red
                        ))
                    nMoved += 1
                if frame.lattsite[i] != saiga12.S12_EMPTYSITE:
                    objs.append(visual.sphere(
                        pos=vector(frame0.coords(i)),
                        radius=1,
                        color = V.vizColors.get(
                            frame.atomtype[frame.lattsite[i]],
                            visual.color.white)
                        ))
                    nMove += 1
            
            print "%.4f"%(float(nMove )/frame.lattSize), \
                  "%.4f"%(float(nMoved)/frame.lattSize),
        sys.stdout.flush()
        frame_index = getNewFrameIndex(frame_index, len(fileNames),
                                       V=V, otherObjects=objs)
        if frame_index == None: break
        
    for i in range(len(objs)):
        objs[0].visible = 0
        del objs[0]
    visual.scene.visible = 0
    print            

def moves(fileNames):
    S = io_open(fileNames[0])
    V = viz.VizSystem(S)
    V.vizMakeBox()

    objs = [ ]
    connMax = S.connMax
    frame_index = 0



    while True:
        print fileNames[frame_index][-20:],
        S = io_open(fileNames[frame_index])
        S.eddEnable()
        nMovesPerParticle = [ ]
        
        for i in range(len(objs)):
            objs[-1].visible = 0
            del objs[-1]
        for i in range(S.N):
            nMovesPerParticle.append(0)
            pos = S.atompos[i]
            startCoords = S.coords(pos)
            for conni in range(connMax):
                if S.MLLr[pos*S.connMax + conni] == -1:
                    continue
                nMovesPerParticle[-1] += 1
                end = S.conn[pos, conni]
                endCoords = S.coords(end)
                    
                axis = endCoords-startCoords
                if math.sqrt(numpy.sum(axis*axis)) > 5:
                    continue
                objs.append(
                    visual.arrow(pos=startCoords, axis=axis,
                                 shaftwidth=.1, headlength=1, fixedwidth=1
                                 ))
        print "all: %.3f %.3f    mobile: %.3f %.3f "%(
            numpy.mean(nMovesPerParticle),
            numpy.std(nMovesPerParticle),
            numpy.mean([ x for x in nMovesPerParticle if x > 0]),
            numpy.std( [ x for x in nMovesPerParticle if x > 0])),
        sys.stdout.flush()
        
        frame_index = getNewFrameIndex(frame_index, len(fileNames))
        if frame_index == None: break
    for i in range(len(objs)):
        objs[-1].visible = 0
        del objs[-1]
    visual.scene.visible = 0
    print        

def visualizeKvectors(fileNames, mode='avg'):
    import saiga12.corrfunc as corrfunc
    scene = visual.scene
    scene.exit = 0
    scene.fov = .5 * scene.fov
    scene.range = 7
    scene.width, scene.height = 800, 800
    fileNames.sort()

    if mode == 'each':
        frame_index = 0
        while True:
            print '\r',
            print fileNames[frame_index][-20:],
            frame = io_open(fileNames[frame_index])

            # XXX we still need a way to set the type of this.
            type_ = list(frame.ntype).index(max(frame.ntype))
            #type_ = 2
            SsfList = corrfunc.StructCorrList(
                frame, kmags=range(1, 8), type_=type_,
                orthogonal=False)
            SsfList.calcSk(frame)
            scene = viz.visualizeKvectors(SsfList)
            sys.stdout.flush()
        
            frame_index = getNewFrameIndex(frame_index, len(fileNames))
            if frame_index == None: break
            #scene.visible = 0
            for i in range(len(scene.objects)):
                scene.objects[0].visible = 0
        
        #for i in range(len(scene.objects)):
        #    scene.objects[0].visible = 0
        visual.scene.visible = 0
        del visual.scene
    if mode == 'avg':
        SsfList = None
        for fname in fileNames:
            frame = io_open(fname)
            type_ = list(frame.ntype).index(max(frame.ntype))
            if SsfList==None: 
                SsfList = corrfunc.StructCorrList(
                    frame, kmag2s=range(1, 8**2), type_=type_,
                    orthogonal=False)
            SsfList.calcSk(frame)
        scene = viz.visualizeKvectors(SsfList)
        while True:
            frame_index = getNewFrameIndex(0, 1)
            if frame_index == None: break
        # clean up
        scene.visible = 0
        for i in range(len(scene.objects)):
            scene.objects[0].visible = 0
    

def msdPosition(S0, S, type_=saiga12.S12_TYPE_ANY, otherS=None):
    """MSD displacement of all particles of a certain type between two times.
    """
    startpos = S0.getPos(type_)
    endpos =   S.getPos(type_)
    d2 = S0.distance2(startpos, endpos, otherS)  # distance squared array
    #print d2
    maxd = max(numpy.sqrt(d2))
    msd = sum(d2) / len(d2)
    return msd


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] in ('--help', '-h', 'help'):
        print """usage: util.py <command> <arguments>...

        commands are:
        hash -- print command line hashes of all filename arguments
        """
    
    elif sys.argv[1] == 'hash':
        for name in sys.argv[2:]:
            S = io_open(name)
            print name, S.hash()

    elif sys.argv[1] == 'diff':
        diff(sys.argv[2], *sys.argv[3:])

    elif sys.argv[1] == 'moves':
        moves(sys.argv[2:])

    elif sys.argv[1] == 'kvecs':
        visualizeKvectors(sys.argv[2:])
    elif sys.argv[1] == 'kvecs-each':
        visualizeKvectors(sys.argv[2:], mode='each')
    elif sys.argv[1] == 'kvecs-avg':
        visualizeKvectors(sys.argv[2:], mode='avg')

    elif sys.argv[1] == 'info':
        for fname in sys.argv[2:]:
            print fname
            S = io_open(fname)
            print "density:", S.density
            for t,n in enumerate(S.ntype):
                if n == 0: continue
                print ("  type %d, density: %0.4f(%0.4f)(%d)"%
                      (t, S.densityOf(t), S.densityOf(t)/S.density,
                       S.numberOfType(t)))
            print "energy:", S.energy()
            print "coords:", S.__class__.__name__, S.lattShape

    else:
        print "command not found: %s"%sys.argv[1]
