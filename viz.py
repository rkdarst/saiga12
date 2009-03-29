# Richard Darst, April 2008

"""Visualization utilities.

This entire module requires python-visual, also known as vpython,
debian package name 'python-visual', (import name `visual`).

This module is primarily used via the `VizSystem` object.
"""

import saiga12
import numpy
import visual

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
    def __init__(self, S):
        self.S = S
        self._display = [ ]
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
        visual.scene.center = visual.scene.originalCenter = (x/2., y/2., z/2.)
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
        display = self._display
        S = self.S
        # Now go add/update all atom positions, etc.
        for i in range(self.S.N):
            pos =    S.atompos[i]
            coords = S.coords(pos)
            type_ =  S.atomtype[i]
            radius = vizRadius.get(type_, self.radius)
            color = vizColors.get(type_, visual.color.white)
            # create the particle if not existing yet:
            if len(display) <= i:
                display.append(visual.sphere(pos=coords, radius=radius))
                display[i].opacity = .2
            # update it if it's already there (yes, it re-sets pos...)
            display[i].visible = 1
            display[i].pos = coords
            if not hasattr(display[i], 's12viz'):
                display[i].color = color
            display[i].radius = radius
        # hide all higher particle numbers:
        for i in range(self.S.N, len(display)):
            display[i].visible = 0

    def __del__(self):
        """Remove all shapes from the display.

        Deleting the object removes all shapes from the display.
        """
        display = self._display
        for i in range(len(self._display)):
            display[0].visible = 0
            del display[0]
        for i in range(len(self._otherObjects)):
            self._otherObjects[0].visible = 0
            del self._otherObjects[0]
tagColor = (.5, 0, .5)
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
def tagToggle():
    """Toggle tag on object under the pointer.
    """
    obj = visual.scene.mouse.pick
    #print obj
    if obj == None: return
    info = gets12viz(obj)
    if info.get('tag', False) == False:  # tag it
        obj.color = tagColor
        info['tag'] = True
    else:                                # untag it
        obj.color = info['color']
        info['tag'] = False
def tagToggle2(otherObjects):
    # Selecting arrows, and other objects for which "pick" does not work
    closest = 1e9
    closestObj = None
    # Iterate through all objects, and select the closest arrow object.
    for obj in otherObjects:
        # only select arrows - other objects can be directly selected.
        if not isinstance(obj, visual.primitives.arrow): continue
        # If this object is closest to pointer, save it for later.
        displacement = visual.mag(visual.scene.mouse.pickpos - obj.pos)
        if displacement < closest:
            closest = displacement
            closestObj = obj
    # Now set the tags like normal
    obj = closestObj
    if obj == None: return
    info = gets12viz(obj)
    if info.get('tag', False) == False:  # tag it
        obj.color = tagColor
        info['tag'] = True
    else:                                # untag it
        obj.color = info['color']
        info['tag'] = False

def toggleViz(V, otherObjects=()):
    if   V._vizMode == 'full': mode = 'show'
    elif V._vizMode == 'show': mode = 'hide'
    elif V._vizMode == 'hide': mode = 'normal'
    elif V._vizMode == 'normal': mode = 'full'
    #print "new mode:", mode
    V._vizMode = mode

    import itertools
    for obj in itertools.chain(V._display, V._otherObjects, otherObjects):
        info = gets12viz(obj, relaxed=True)   # we could optimize by not using
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
def toggleBG():
    print visual.scene.background, visual.color.white, visual.color.black
    # is white, make it black
    if visual.scene.background == visual.color.white:
        visual.scene.background = visual.color.black
    # is black, make it white
    elif visual.scene.background == visual.color.black:
        visual.scene.background = visual.color.white

                    
if __name__ == "__main__":
    import saiga12.io
    import saiga12.util
    import sys
    try:
        from urllib import urlopen as open
    except ImportError:
        pass
                        

    fileNames = sys.argv[1:]
    frame_index = 0
    wait = True
    # bash filename sort ignores "-" which is annoying for negative numbers
    fileNames.sort()  
    #print fileNames
    V = None

    while True:
        fname = fileNames[frame_index]
        S = saiga12.io.io_open(open(fname))

        visual.scene.title = fname
        if V == None or V.S.lattSize != S.lattSize:
            V = VizSystem(S)
            V.vizMakeBox()
        V.S = S
        V.vizDisplay()
        

        print fname,
        sys.stdout.flush()
        frame_index = saiga12.util.getNewFrameIndex(
                                              frame_index,len(fileNames), V=V)
        if frame_index == None: break
    del V
    print
    visual.scene.visible = 0


