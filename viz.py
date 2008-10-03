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

    def vizMakeBox(self):
        """Put the simulation box on the visual display.

        This method should only be called once (but won't die if you
        call it more than once, but it might leak memory long-term).
        Visual keeps references to objects around...
        """
        if visual is None:
            print "visual module was not imported, visual siletly deactivated."
            return
        visual.scene.center = numpy.asarray(self.S.physicalShape) / 2.
        radius = .02
        if len(self.S.physicalShape) == 2:
            x, y = self.S.physicalShape
            z = 0
        else:
            x,y,z = self.S.physicalShape
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
            display[i].visible = 0
        for obj in self._otherObjects:
            obj.visible = 0
                    
if __name__ == "__main__":
    import saiga12.io
    import sys
    try:
        from urllib import urlopen as open
    except ImportError:
        pass
                        

    wait = True

    for i, fname in enumerate(sys.argv[1:]):
        if fname == "-w":
            wait = False
            continue
        
        S = saiga12.io.io_open(open(fname))

        visual.scene.title = fname
        V = VizSystem(S)
        V.vizMakeBox()
        V.vizDisplay()

        print fname,
        if wait:
            raw_input(", ...")


