import unittest

import glob

class TestSequenceFunctions(unittest.TestCase):
    #def testAll(self):
    #    files = glob.glob("tests/*.py")
    #    if "unittests_run.py" in files:
    #        files.remove("unittests_run.py")
    #    for x in files:
    #        print x
    #        execfile(x)
    #def testChemPotential(self):
    #    execfile("tests/chempotential.py", {})
    def testCorrFunc(self):
        execfile("tests/corrfunctest.py", {'noviz':True})
    def testDecay(self):
        # writes out to logfile.txt, don't use this test.
        execfile("tests/decay.py", {})
    def testDiffusion(self):
        execfile("tests/diffusion.py", {})
    def testDistances(self):
        execfile("tests/distances.py", {})


    def testFA(self):
        execfile("tests/fredricksonandersentest.py", {})
    def testHexGrids(self):
        execfile("tests/hextest.py", {'noviz':True})
    #def testInsertType(self):
    #    # needs to be updated
    #    execfile("tests/inserttype.py", {})
    def testIO(self):
        execfile("tests/iotest.py", {})
    def testKobAndersen(self):
        execfile("tests/kobandersentest.py", {})
    def testEventDrivenDynamics(self):
        execfile("tests/eventdrivendynamics.py", {})
    #def testStateFinder(self):
    #    execfile("tests/statefindertest.py", {})
    def testTiming(self):
        execfile("tests/timing.py", {})
    #def testviztest(self):
    #    execfile("tests/viztest.py", {})

print unittest.main()

#suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
#unittest.TextTestRunner(verbosity=2).run(suite)
