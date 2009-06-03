# Richard Darst, April 2008
# Part of Saiga-12


S12_EMPTYSITE = -1
inf = float("inf")
S12_TYPE_ANY = -2

S12_ENERGY_AVAIL  = (1, 2, 3, 10, 11, 12)
S12_ENERGY_BM     = 1             # biroli-mezard thermodynamic lattice glass
S12_ENERGY_ZERO   = 2             # all configs have zero energy
S12_ENERGY_CTCC   = 3             # an orientation and an overlap exclusion.
# somewhat experimental modes:
S12_ENERGY_BMnotzero     = 10
S12_ENERGY_BMimmobile1   = 11
S12_ENERGY_BMimmobile1b  = 12
S12_CYCLE_MC      = 1             # monte carlo
S12_CYCLE_KA      = 2             # kob-andersen
S12_CYCLE_FA      = 3             # fredrickson-andersen
S12_CYCLE_CTCC    = 4             # an orientation and an overlap exclusion.
S12_CYCLE_CTCCclassic    = 10     #
S12_FLAG_VIB_ENABLED     = 1      # vibrations are enabled
