# Richard drast, September 2008

import subprocess
import os

import saiga12

g = saiga12.Grid2d()
g.makegrid(10, 10)
g.addParticles({2:.25})
g.ioSave("tests/test-tmp-00.ugh")
subprocess.call(("python", "-m", "saiga12.util", "hash",
                 "tests/test-tmp-00.ugh"),
                stdout=subprocess.PIPE)
os.unlink("tests/test-tmp-00.ugh")
