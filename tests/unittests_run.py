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
    def testCtccTest(self):
        execfile("tests/ctcctest.py", {})
    def testCtccEDDTest(self):
        execfile("tests/ctcceddtest.py", {'short':True})
    def testDecay(self):
        # writes out to logfile.txt, don't use this test.
        execfile("tests/decay.py", {})
    def testDiffusion(self):
        execfile("tests/diffusion.py", {})
    def testDistances(self):
        execfile("tests/distances.py", {})
    def testEast(self):
        execfile("tests/easttest.py", {})
    def testEnergymodes(self):
        execfile("tests/energymodes.py", {})


    def testFA(self):
        execfile("tests/fredricksonandersentest.py", {})
    def testFAsoft(self):
        execfile("tests/fredricksonandersensofttest.py", {})
    def testFrozen(self):
        execfile("tests/frozentest.py")
    def testGridUnitTests(self):
        execfile("tests/gridunittests.py", {})
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
    def testLogTime(self):
        execfile("tests/logtimetest.py", {})
    #def testStateFinder(self):
    #    execfile("tests/statefindertest.py", {})
    def testTiming(self):
        execfile("tests/timing.py", {'fast':True})
    #def testviztest(self):
    #    execfile("tests/viztest.py", {})
    def testScripts(self):
        execfile("tests/scripttest.py", {})
    def testShuffle(self):
        execfile("tests/shuffletest.py", {})
    def testSpinGlass(self):
        execfile("tests/spinglasstest.py", {})
    def testSpiral(self):
        execfile("tests/spiraltest.py", {})
    def testSquarePlaquette(self):
        execfile("tests/squareplaquettetest.py", {'short': True})
    def testVibrations(self):
        execfile("tests/vibrationtest.py", {})

import sys
sys.stdout = file('/dev/null', 'w')

print unittest.main()
sys.stdout = sys.__stdout__

#suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
#unittest.TextTestRunner(verbosity=2).run(suite)
