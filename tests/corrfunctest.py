# Richard Darst, 2008

import math
import numpy


import saiga12.io
import saiga12.corrfunc

from rpy import r

       # (mu, type_)
runs = ((1.0, 1),
        (1.0, 3),
        #(5.0, None),
        #(11.0,None),
        (11.0, 1),
        (11.0, 3),
        )
iterations = -1

for i, (mu, type_) in enumerate(runs):

    fname = "/home/richard/research/lgm/bm/crstl-1d-2-runs/"\
            "output_mu1.0/trj_mu1.0_mctime00005.ugh"
    S = saiga12.io.io_open(file(fname))
    SsfList = [ ]
    for kmag2 in range(1, 50):
        Ssf = saiga12.corrfunc.StructCorr(kmag2=kmag2,
                                          S=S,
                                          type_=type_)
        if len(Ssf.kvecs) == 0:
            continue
        Ssf.kvecsOrig = Ssf.kvecs.copy()
        Ssf.kvecs *= (2*math.pi / 15.)
        SsfList.append(Ssf)


    print mu
    thisRun = [ ]

    def generateFrame():
        import glob
        fnames = "/home/richard/research/lgm/bm/crstl-1d-2-runs/output_mu%s/trj_mu%s_mctime*.ugh"%(mu, mu)
        fnames= glob.glob(fnames)
        fnames.sort()
        for i, fname in enumerate(fnames):
            #print fname
            if i == iterations: return
            S = saiga12.io.io_open(file(fname))
            #from rkddp.interact import interact ; interact()
            yield S

    for S in generateFrame():
        for Ssf in SsfList:
            #s =
            Ssf.staticStructureFactor(S)
            #Ssf.avgStore('ssf', s)
            #print Ssf._avgs
        
    for Ssf in SsfList:
        thisRun.append((Ssf.kmag, Ssf.Sk()))
    v_kmag = zip(*thisRun)[0]
    v_ssf = zip(*thisRun)[1]

    if mu > 10: pch = "x"
    else:       pch = "+" #1
    if i == 0:
        r.plot(v_kmag, v_ssf,
               xlab="", ylab="", type="l", col=i+1,
               ylim=(0., 15)
               #ylim=(0., 15000)
               )
        for kmag in v_kmag:
            r.abline(v=kmag, lty=3, col="lightgray")
    else:
        r.lines(v_kmag, v_ssf, type="l", col=i+1)
    for Ssf in SsfList:
        r.points(x=[Ssf.kmag]*len(Ssf.kvecs), y=Ssf.SkArray(),
                 col=i+1, pch=pch)
    

import code ; code.interact(local=locals(), banner="" )
    
