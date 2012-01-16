# Richard Darst, April 2008

"""Visualization utilities.

This entire module requires python-visual, also known as vpython,
debian package name 'python-visual', (import name `visual`).

This module is primarily used via the `VizSystem` object.
"""

import itertools

import saiga12
from saiga12 import S12_EMPTYSITE
import numpy
import visual
import math
from math import log, exp, sqrt

class VizSystem(object):
    """Visualize the system using python-visual.

    Here is how it works:
      V = saiga12.viz.VizSystem(S)  # note, initialize with System object
      V.vizMakeBox()                # from then on, you just use the V obj
      V.vizDisplay()

    then, to update atom locations, call display again:
      V.vizDisplay()

    There isn't yet an interface on System objects to display things--
    maybe someday?

    RADIUS:
    To configure default particle radius, set
       V.radius = radius
    Default particle radious is .25
    To set a radius for a specific type of particle, use
       V.vizRadius[particleType] = radius

    COLOR:
    To configure particle colors, use:
       V.vizColors[particleType] = color
    To see what colors are available, do `pydoc visual.color` for
    pre-defined colors and visual docs.
    Default colors are 0:red, 1:white, 2:cyan, 3:green,
                       4:darkblue, 5:darkgreen, 6:blue, otherwise white.

    """
    tagColor = (.5, 0, .5)
    def __init__(self, S, mode="atoms", S0=None, limit=None):
        self.S = S
        self.S0 = S0
        self._mode = mode
        self._limit = limit
        self.scene = visual.scene
        self._atoms = [ ]
        self._otherObjects = [ ]
        self.vizColors = {
            0: visual.color.red,    # (1,0,0)
            1: visual.color.white,
            2: visual.color.cyan,
            3: visual.color.green,  # (0,1,0)
            4: (0.,  .5, 0. ),
            5: (0., 0. ,  .5),
            6: visual.color.blue,   # (0,0,1)
            }
        self.radius = .25
        self.vizRadius = { }
        self._vizMode = 'full'  # 'full', 'show', 'hide'

    def vizMakeBox(self):
        """Put the simulation box on the visual display.

        This method should only be called once (but won't die if you
        call it more than once, but it might leak memory long-term).
        Visual keeps references to objects around...
        """
        if visual is None:
            print "visual module was not imported, visual siletly deactivated."
            return
        radius = .02
        if len(self.S.physicalShape) == 1:
            x, = self.S.physicalShape
            y, z = 0, 0
        elif len(self.S.physicalShape) == 2:
            x, y = self.S.physicalShape
            z = 0
        else:
            x,y,z = self.S.physicalShape
        self.scene.center = self.scene.originalCenter = (x/2., y/2., z/2.)
        c = visual.color.blue
        self._otherObjects.extend((
        visual.cylinder(pos=(0,0,0), axis=(x, 0, 0), radius=radius, color=c),
        visual.cylinder(pos=(0,y,0), axis=(x, 0, 0), radius=radius, color=c),
        visual.cylinder(pos=(0,0,z), axis=(x, 0, 0), radius=radius, color=c),
        visual.cylinder(pos=(0,y,z), axis=(x, 0, 0), radius=radius, color=c),
        
        visual.cylinder(pos=(0,0,0), axis=(0, y, 0), radius=radius, color=c),
        visual.cylinder(pos=(x,0,0), axis=(0, y, 0), radius=radius, color=c),
        visual.cylinder(pos=(0,0,z), axis=(0, y, 0), radius=radius, color=c),
        visual.cylinder(pos=(x,0,z), axis=(0, y, 0), radius=radius, color=c),
        
        visual.cylinder(pos=(0,0,0), axis=(0, 0, z), radius=radius, color=c),
        visual.cylinder(pos=(x,0,0), axis=(0, 0, z), radius=radius, color=c),
        visual.cylinder(pos=(0,y,0), axis=(0, 0, z), radius=radius, color=c),
        visual.cylinder(pos=(x,y,0), axis=(0, 0, z), radius=radius, color=c),
        ))

    def vizDisplay(self):
        if   self._mode == 'atoms':    self._displayAtoms()
        elif self._mode == 'persist':  self._displayPersist()
        elif self._mode == 'overlap':  self._displayOverlap()
        elif self._mode == 'diff':     self._displayDiff()
        elif self._mode == 'moves':    self._displayMoves()
        elif self._mode == 'empty':    self._displayEmpty()
        else: print "Unknown display mode: '%s'"%self._mode
    __call__ = vizDisplay
    def cycleMode(self):
        """Cycle through the available visualization modes
        """
        if self.S0 is not None:
            modes = ['atoms', 'persist', 'overlap', 'diff', 'moves', 'empty']
        else:
            modes = ['atoms', 'persist', 'empty']
        modes = ['atoms', 'empty']
        mode = self._mode
        if mode not in modes:
            return
        # Find new mode
        i = modes.index(mode)
        i = (i+1) % len(modes)
        newmode = modes[i]
        print "New mode %s"%newmode
        # Remove old atoms, set new mode and re-display everything
        for i in range(len(self._atoms)):
            self._atoms[0].visible = 0
            del self._atoms[0]
        self._mode = newmode
        self()
    
    def _displayAtoms(self):
        """Display the position of atoms in the system.

        Call this the first time to display the atoms on the visual.
        Call it again to update all the positions.  You don't need to
        re-make the box at every step, it behaves smartly and just
        moves the particles around.
        """
        if visual is None:
            return
        vizColors = self.vizColors
        vizRadius = self.vizRadius
        atoms = self._atoms
        S = self.S
        # Now go add/update all atom positions, etc.
        for i in range(self.S.N):
            pos =    S.atompos[i]
            coords = S.coords(pos, raw=True)
            type_ =  S.atomtype[i]
            radius = vizRadius.get(type_, self.radius)
            color = vizColors.get(type_, visual.color.white)
            # create the particle if not existing yet:
            if len(atoms) <= i:
                atoms.append(visual.sphere(pos=coords, radius=radius))
                atoms[i].opacity = .2
            # update it if it's already there (yes, it re-sets pos...)
            atoms[i].visible = 1
            atoms[i].pos = coords
            if not hasattr(atoms[i], 's12viz'):
                atoms[i].color = color
            atoms[i].radius = radius
        # hide all higher particle numbers:
        for i in range(self.S.N, len(atoms)):
            atoms[i].visible = 0
    def _displayPersist(self):
        """Display the position of atoms in the system.

        Call this the first time to display the atoms on the visual.
        Call it again to update all the positions.  You don't need to
        re-make the box at every step, it behaves smartly and just
        moves the particles around.
        """
        if visual is None:  return
        if self.S.persist is None: return
        vizColors = self.vizColors
        vizRadius = self.vizRadius
        atoms = self._atoms
        S = self.S
        # Now go add/update all atom positions, etc.
        for pos in range(S.lattSize):
            coords = S.coords(pos, raw=True)
            if S.lattsite[pos] != S12_EMPTYSITE:
                type_ =  S.atomtype[S.lattsite[pos]]
                radius = vizRadius.get(type_, self.radius)
                color = vizColors.get(type_, visual.color.white)
            else:
                type_ = 0
                radius = self.radius
                color = visual.color.red
            # non-persisted atoms are invisible
            if S.persist[pos]:
                visible = 1
            else:
                visible = 0
            if self._limit and pos not in self._limit:
                visible = 0
            # Show particle:
            self._setatom(i=pos, coords=coords, radius=radius, visible=visible,
                          color=color)
    def _displayOverlap(self):
        """Display the position of atoms in the system.

        Call this the first time to display the atoms on the visual.
        Call it again to update all the positions.  You don't need to
        re-make the box at every step, it behaves smartly and just
        moves the particles around.
        """
        if visual is None:  return
        if self.S0 is None: return
        atoms = self._atoms
        S = self.S
        S0 = self.S0

        type_ = 0
        radius = self.radius
        color = visual.color.white
        # Now go add/update all atom positions, etc.
        for pos in range(S.lattSize):
            coords = S.coords(pos, raw=True)
            if S.lattsite[pos]  != S12_EMPTYSITE and \
               S0.lattsite[pos] != S12_EMPTYSITE:
                visible = 1
            else:
                visible = 0
            if self._limit and pos not in self._limit:
                visible = 0
            # Show particle:
            self._setatom(i=pos, coords=coords, radius=radius, visible=visible)

    def _displayDiff(self):
        """Display arrows for where particles have moved."""
        if self.S0 is None: return
        objs = self._atoms
        S = self.S
        S0 = self.S0
        # Delete all objects:
        for i in range(len(objs)):
            objs[-1].visible = 0
            del objs[-1]
        for i in range(S.N):
            startPos = S.atompos[i]
            endPos = S0.atompos[i]
            if startPos not in self._limit:
                continue
            if startPos == endPos:
                objs.append(visual.sphere(pos=S.coords(startPos), radius=.15))
            else:
                startCoords = S.coords(startPos)
                endCoords   = S0.coords(endPos)
                dist = math.sqrt(sum((endCoords-startCoords)**2))
                if dist > S.L/2:
                    continue
                objs.append(
                    visual.arrow(pos=startCoords, axis=endCoords-startCoords,
                                 shaftwidth=.1, headlength=1, fixedwidth=1
                                 ))

    def _displayDiffParticle(self):
        pass
    def _displayDiffSpin(self):
        pass
    def _displayMoves(self):
        pass

    def _displayEmpty(self):
        """Display a sphere at every EMPTY lattice site.
        """
        if visual is None:  return
        atoms = self._atoms
        S = self.S

        type_ = 0
        radius = self.radius
        color = visual.color.white
        # Now go add/update all atom positions, etc.
        for pos in range(S.lattSize):
            coords = S.coords(pos, raw=True)
            if S.lattsite[pos]  == S12_EMPTYSITE:
                visible = 1
            else:
                visible = 0
            if self._limit and pos not in self._limit:
                visible = 0
            # Show particle:
            self._setatom(i=pos, coords=coords, radius=radius, visible=visible)

    def _setatom(self, i, coords, radius, color=visual.color.white,visible=1):
        """Add atom to self._objects if needed, otherwise sets its properties.
        """
        atoms = self._atoms
        if len(atoms) <= i:
            atoms.append(visual.sphere(pos=coords, radius=radius,
                                       visible=visible))
            atoms[i].opacity = .2
        # update it if it's already there (yes, it re-sets pos...)
        atoms[i].visible = visible
        atoms[i].pos = coords
        if not hasattr(atoms[i], 's12viz'):
            atoms[i].color = color
        atoms[i].radius = radius


    def __del__(self):
        """Remove all shapes from the display.

        Deleting the object removes all shapes from the display.
        """
        atoms = self._atoms
        for i in range(len(self._atoms)):
            atoms[0].visible = 0
            del atoms[0]
        for i in range(len(self._otherObjects)):
            self._otherObjects[0].visible = 0
            del self._otherObjects[0]
    def tagToggle(self, obj=None):
        """Toggle tag on object under the pointer.
        """
        # Pick object if not passed
        if obj is None:
            obj = self.scene.mouse.pick
        # If no object found, then we don't do anything.
        if obj is None: return
        info = gets12viz(obj)
        if info.get('tag', False) == False:  # tag it
            obj.color = self.tagColor
            info['tag'] = True
        else:                                # untag it
            obj.color = info['color']
            info['tag'] = False
    def tagToggle2(self):
        # Selecting arrows, and other objects for which "pick" does not work
        closest = 1e9
        closestObj = None
        # Iterate through all objects, and select the closest arrow object.
        for obj in self._otherObjects:
            # only select arrows - other objects can be directly selected.
            if not isinstance(obj, visual.primitives.arrow): continue
            # If this object is closest to pointer, save it for later.
            displacement = visual.mag(self.scene.mouse.pickpos - obj.pos)
            if displacement < closest:
                closest = displacement
                closestObj = obj
        # Now set the tags like normal
        obj = closestObj
        if obj == None: return
        self.tagToggle(obj=obj)

    def cycleViz(self, otherObjects=()):
        """Toggle between highlighting modes in the 'atoms' mode.
        """
        if self.mode != 'atoms': return
        if   self._vizMode == 'full': mode = 'show'
        elif self._vizMode == 'show': mode = 'hide'
        elif self._vizMode == 'hide': mode = 'normal'
        elif self._vizMode == 'normal': mode = 'full'
        #print "new mode:", mode
        self._vizMode = mode

        for obj in itertools.chain(self._atoms, self._otherObjects,
                                   otherObjects):
            info = gets12viz(obj, relaxed=True)
            if mode == 'show':
                if info == None:
                    obj.visible = False
                else:
                    if info['tag'] == True:
                        obj.visible = True
                        obj.color = info['color']
                    else:  obj.visible = False
            elif mode == 'hide':
                if info == None:
                    obj.visible = True
                else:
                    if info['tag'] == True:
                        obj.visible = False
                    else:  obj.visible = True
            elif mode == 'full':
                if info == None:
                    obj.visible = True
                else:
                    if info['tag'] == True:
                        obj.visible = True
                        obj.color = tagColor
                    else:  obj.visible = True
            elif mode == 'normal':
                if info == None:
                    obj.visible = True
                else:
                    if info['tag'] == True:
                        obj.visible = True
                        obj.color = info['color']
                    else:  obj.visible = True
    def toggleBG(self):
        print self.scene.background, visual.color.white, visual.color.black
        # is white, make it black
        if self.scene.background == visual.color.white:
            self.scene.background = visual.color.black
        # is black, make it white
        elif self.scene.background == visual.color.black:
            self.scene.background = visual.color.white


def gets12viz(obj, relaxed=False):
    """Get visualization paremeter dictionary.
    """
    if not hasattr(obj, "s12viz"):
        if relaxed: return None
        else: obj.s12viz = { }
    info = obj.s12viz
    if not info.has_key('color'):
        info['color'] = obj.color
    return info


def visualizeKvectors(SsfList):
    scene = visual.scene
    #visual.scene = visual.display()
    #scene = visual.scene
    #scene.exit = 0
    #scene.fov = .5 * scene.fov
    #scene.range = 7
    #scene.width, scene.height = 800, 800
    import saiga12.util
    avg = saiga12.util.Averager()
    for kmag in SsfList.kmags():
        Ssf = SsfList.SsfDict[kmag]
        for i, (kvec, Sk) in enumerate(
              zip(Ssf.kvecsOrig, Ssf.SkArraysByKvec())):
            avg.add(Sk)
            #print kvec, Sk

            widthfactor = .2
            col = (Sk-.25)*.25
            #col = Sk
            col = min(1., max(.1,  col))
            col = (col, col, col)
            radius = max(.05, (Sk)*widthfactor)
            #visual.arrow(pos=(0,0,0), axis=kvec,
            #             shaftwidth=Sk*widthfactor,
            #             fixedwidth=1)
            visual.cylinder(pos=(0,0,0), axis=kvec, color=col,
                            radius=radius)
    print "mean:", avg.mean
    return scene

                    
if __name__ == "__main__":
    import saiga12.io
    import saiga12.util
    import sys
    try:
        from urllib import urlopen as open
    except ImportError:
        pass

    mode = "atoms"
    if sys.argv[1] == "persist":
        mode = 'persist'
        del sys.argv[1]
    if sys.argv[1] == "empty":
        mode = 'empty'
        del sys.argv[1]

    arguments = sys.argv[1:]
    # bash filename sort ignores "-" which is annoying for negative numbers
    #print arguments
    arguments.sort()
    listOfFrames, listOfNames = saiga12.util.openFiles(arguments)

    frame_index = 0
    wait = True
    V = None
    visual.scene.width, visual.scene.height = 800, 800

    while True:
        S = listOfFrames[frame_index]
        if isinstance(S, str):
            # We have delayed opening of individual objects, open them now
            S = saiga12.io.io_open(S)

        if V == None or V.S.lattSize != S.lattSize:
            V = VizSystem(S, mode=mode)
            V.vizMakeBox()
        V.S = S
        V.vizDisplay()


        print listOfNames[frame_index][-20:], \
                " d=%04.4f e=%4.4g "%(S.density, S.energy()),
        sys.stdout.flush()
        frame_index = saiga12.util.getNewFrameIndex(
                                           frame_index,len(listOfFrames), V=V)
        if frame_index == None: break
    del V
    print
    visual.scene.visible = 0


