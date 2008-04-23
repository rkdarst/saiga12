

import saiga12
import saiga12.statefinder
from saiga12.geom.grid import Grid3d

a = b = c = 15
S = Grid3d(N=a*b*c)
S.makegrid(a, b, c)

initialdensity = 1/2.45
S.hardness = 1.
S.addParticleRandomDensity(initialdensity * .3, 1)  # A
S.addParticleRandomDensity(initialdensity * 1., 3)  # B
S.anneal()
S.setInsertType( {1: (.3, 1.00),
                  3: (.7, 0.82)} )


logfile = file("logfile-sft.txt", "w")


SF = saiga12.statefinder.StateFinder()
SF.S = S
SF.mu3 = .82
SF.mu1 = 1.
SF.logfile = logfile

SF.begin()
SF.status()
SF.S.resetTime()
while SF.mu1 <= 8.:
    S.avgReset()
    for j in range(10):
        SF.runPass()
        if j == 5:
            SF.S.avgReset()
    S.io_writeToFile("outputs/gceeq_mu1=%s.ugh"%SF.mu1)
    SF.mu1 += .1
#SF.runPass()
