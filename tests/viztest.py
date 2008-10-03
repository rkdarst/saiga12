# Richard Darst, 2008

#Generic tests:
import time
import saiga12
from saiga12.geom.grid import Grid3d
import saiga12.viz

S = Grid3d()
S.makegrid(10, 10, 10)
S.addParticleRandomDensity(.05, 0)
S.addParticleRandomDensity(.10, 1)
S.addParticleRandomDensity(.15, 2)
S.addParticleRandomDensity(.20, 3)
S.addParticleRandomDensity(.25, 4)
S.addParticleRandomDensity(.30, 5)
S.addParticleRandomDensity(.35, 6)
S.cycle(100000)

V = saiga12.viz.VizSystem(S)

V.vizColors[4] = saiga12.viz.visual.color.white
V.vizRadius[4] = .1
V.vizRadius[5] = .5
# V.radius = .5  # set default radius

V.vizMakeBox()
V.vizDisplay()
V.vizDisplay()

S.addParticleRandomDensity(.40, 6)
V.vizDisplay()

# del V   # removes all objects
time.sleep(2)




#Try a grand canonical test:
S = Grid3d()
S.makegrid(10, 10, 10)
S.setCycleMode('fredricksonandersen')
S.beta = 1.
S.addParticleRandomDensity(.10, 1)
S.inserttype = 1
S.eddEnable()

V = saiga12.viz.VizSystem(S)
V.vizMakeBox()
V.vizDisplay()

maxN = 0
while True:
    S.cycle(100)
    maxN = max(maxN, S.N)
    print S.mctime, maxN, S.N, len(V._display), S.beta
    #time.sleep(.5)
    V.vizDisplay()
    if S.mctime > 2000: S.beta = 1/.1
    if S.mctime > 5e6: S.beta = 1
