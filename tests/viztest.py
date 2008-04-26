

# 2D

import saiga12
from saiga12.geom.grid import Grid2d
import saiga12.viz


S = Grid2d()
S.makegrid(25, 25)
S.addParticleRandomDensity(.05, 0)
S.addParticleRandomDensity(.10, 1)
S.addParticleRandomDensity(.15, 2)
S.addParticleRandomDensity(.20, 3)
S.addParticleRandomDensity(.25, 4)
S.addParticleRandomDensity(.30, 5)
S.addParticleRandomDensity(.35, 6)
S.cycle(100000)
V = saiga12.viz.VizSystem(S)


V.vizMakeBox()
V.vizDisplay()
