
import saiga12
from saiga12.geom.grid import Grid3d

beta = 1/3.

def getS():
    S = Grid3d()
    S.beta = beta
    S.makegrid(16, 16, 16)
    S.hardness = 1
    S.addParticles({1:.15, 3:.35})
    S.setCycleMoves()
    S.anneal(False)
    S.hardness = 1
    return S


S = getS()
for i in range(10):
    print "setting to notzero"
    S.setEnergyMode(saiga12.S12_ENERGY_BMnotzero)
    print S.energy() # the first after the change should have wild behavior
    S.cycle(100) ; print S.energy()
    #S.cycle(100) ; print S.energy()
    #S.cycle(100) ; print S.energy()
    print "setting to immobile1"
    S.setEnergyMode(saiga12.S12_ENERGY_BMimmobile1)
    print S.energy() # the first after the change should have wild behavior
    S.cycle(100) ; print S.energy()
    #S.cycle(100) ; print S.energy()
    #S.cycle(100) ; print S.energy()
