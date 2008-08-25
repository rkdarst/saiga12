import numpy

def findAsymmetry(S, doAssert=False):
    """Find if position B not connected to A, if A is connected to B.

    This iterates through all lattice positions.  Look at all of it's
    connections.  If the adject position does not connect back, print
    out an error message.
    """
    for n in range(S.lattSize):
        conn = S.conn[n]
        for i, adjPos in enumerate(conn):
            if n not in S.conn[adjPos]:
                print "missing connection: %s is connected to %s (%s), but %s is not connected to %s"%(n, adjPos, i, adjPos, n)
            assert (not doAssert) or \
                   n in S.conn[adjPos], "Asymmetric connection: see above"



def drawConn(V, S, pos, whichi=None):
    """Debug -- draw all connections on the visual, for some position

    Using visual object V, draw all connections from lattice site
    `pos` to all of it's neighbors.  If `whichi` is given, then only
    draw these connection numbers.
    """
    import visual
    start = S.coords(pos)
    for i, adjpos in enumerate(S.conn[pos]):
        if whichi is not None and i not in whichi:
            continue
        end = S.coords(adjpos)
        V._otherObjects.append(
            visual.arrow(pos=start, axis=end - start, shaftwidth=.1,
                         headlength=.5, fixedwidth=1
                         ))


def checkConnDistances(S, distanceShouldBe=1, setType=None,
                       doAssert=False):
    """For every latt

    If `settype` is given, then change the atomtype of a particle at
    that position to the given type.  This can help to get a visual
    indication of what is connected improperly (by color).

    density=.999
    type_ = 12
    S = GridHex3d()
    S.addParticleRandomDensity(density, type_=type_)

    """
    for i in range(S.lattSize):
        x = S.distance(i, S.conn[i])
        # print connection distances
        print i, x
        # change atomtype if it's wrong
        if setType and numpy.any(numpy.abs(x - distanceShouldBe) > 1e-5):
            S.atomtype[S.lattsite[i]] = setType
            print i
        assert (not doAssert) or \
               not numpy.any(numpy.abs(x - distanceShouldBe) > 1e-5), \
               "Connection distances outside of normal parameters"
