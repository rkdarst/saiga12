# Richard Darst, April 2008

import saiga12
import numpy
import visual

class VizSystem(object):
    def __init__(self, S):
        self.S = S
        self._otherObjects = [ ]

    def vizMakeBox(self):
        """Put the simulation box on the visual display"""
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
    vizColors = {
        0: visual.color.red,    # (1,0,0)
        1: visual.color.white,
        2: visual.color.cyan,
        3: visual.color.green,  # (0,1,0)
        4: (0.,  .5, 0. ),
        5: (0., 0. ,  .5),
        6: visual.color.blue,   # (0,0,1)
        }

    def vizDisplay(self):
        if visual is None:
            return
        c = self.vizColors
        if not hasattr(self, "_display"):
            display = [ ]
            #print self.S.N
            #print self.S.lattShape
            for i in range(self.S.N):
                pos = self.S.atompos[i]
                coords = self.S.coords(pos)
                display.append(visual.sphere(pos=coords,
                                             radius=.25, #r, not d
                         color=c.get(self.S.atomtype[i], visual.color.white)))
                display[-1].opacity = .2
                self._display = display
        else:
            display = self._display
            for i in range(self.S.N):
                pos = self.S.atompos[i]
                coords = self.S.grid_coords(pos)
                display[i].pos = coords
    def __del__(self):
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


