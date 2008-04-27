# Richard Darst, 2008


import saiga12

class StateFinder(object):
    # S = <set manually>
    # S.mu2 = <set manually at a guess>
    skip = 100
    cyclelen = 100
    logfile = None
    controlScale = 1.
    def setCE(self):
        self.S.setCycleMoves(shift=self.S.N)
    def setGCE(self):
        self.S.setCycleMoves(shift=0, insertdel=self.S.N)
    def begin(self, total, skip):
        # equilibrate in CE
        self.setCE()
        for i in range(total // skip):
            self.S.cycle(skip)
            self.status(log=False)
        self.mu3 = self.S.avg("mu3")
    def actualX(self):
        return self.S.densityOf(1) / self.S.densityOf(3)
        #return self.S.densityOf(3) / self.S.density
    def targetX(self):
        return self.fracA / (1 - self.fracA)
        
    def adjustMu(self):
        S = self.S

        # Establish how far from ideal we are
        deltaX = (self.actualX() - self.targetX()) * self.controlScale
        #deltaX = -deltaX
        #print "target: %.4f, actual: %.4f, delta: %.4f"%\
        #      (self.targetX(), self.actualX(), deltaX),

        # set us to our new proper state
        #print "old mu3: %.4f"%self.mu3,
        self.mu3 += deltaX
        #print "new mu3: %.4f"%self.mu3
        S.setInsertType( {1: (self.fracA   , self.mu1),
                          3: (1-self.fracA , self.mu3)} )
        #print


        ## run at the new conditions, GCE -- get number of particles closer.
        #self.setGCE()
        #for i in range(self.cyclelen):
        #    S.cycle(self.skip)
        #    self.status()


    def eqCE(self):
        S = self.S
        

    def status(self, log=True, printstdout=True):
        S = self.S
        logfile = self.logfile
        S.avgStore("density", S.density)
        S.avgStore("density1", S.densityOf(1))
        S.avgStore("density3", S.densityOf(3))
        mu1 = S.chempotential(1, store=False)
        mu3 = S.chempotential(3, store=False)
        S.avgStore("mu1", mu1)
        S.avgStore("mu3", mu3)

        if printstdout:
            print "\033[2A\r"
            print "\r", S.mctime, S.N, \
                  " ", (5*"%.4f ")%(S.density,
                                    S.densityOf(1),
                                    S.densityOf(3),
                                    S.avg("mu1"),
                                    S.avg("mu3"))
        if log and logfile:
            print >> logfile, \
                  S.mctime, (7*"%0.6f ")%(S.avg("density"),
                                          S.densityOf(1),
                                          S.densityOf(3),
                                          S.avg("mu1"),
                                          S.avg("mu3"),
                                          self.mu1,
                                          self.mu3)
            #logfile.flush()
        

