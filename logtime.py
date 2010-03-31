# Richard Darst, October 2008

from math import log

def logTime(initialTime, timeIndex):
    """
    timeIndex = 0
    while S.mctime <= maxTime:
        timeIndex += 1
        nextTime = fautil.logTime(1, timeIndex)
        if nextTime <= S.mctime: continue
        cycleTime = nextTime - S.mctime
        S.cycle(cycleTime)
    """
    # 1.25892541179417 = exp(log(10) / 10)
    return int(initialTime * 1.25892541179417**timeIndex)

def getAllTimes(maxTime):
    """Return all times in the set less than this maxTime.
    """
    times = [ ]
    timeIndex = 0
    curTime = 0
    times.append(curTime)
    while curTime < maxTime:
        timeIndex += 1
        nextTime = logTime(1, timeIndex)
        if nextTime <= curTime: continue
        curTime = nextTime
        times.append(curTime)
    return times
def logTimeTimesGenerator():
    timeIndex = 0
    curTime = 0
    yield curTime
    while True:
        timeIndex += 1
        nextTime = logTime(1, timeIndex)
        if nextTime <= curTime: continue
        curTime = nextTime
        yield curTime



def roundTime(time):
    if time == 0:
        return time
    initialTime = 1.0
    x = log(time) * 10. / (log(10.) * initialTime)
    x = int(round(x, 0))
    timeRounded = int(1.0 * 1.25892541179417**x) # 1.0 decreed elsewhere
    return timeRounded


def runLogTime(S, maxTime, callback, ):
    S0 = S.copy()

    for S, deltaT in logTimeGenerator(S, maxTime):
        callback(S0=S0, S=S, deltaT=curTime, )
    # return new state
    return S

def logTimeGenerator(S, maxTime):
    """Iterate over (S, deltaT) pairs, excluding deltaT = 0.
    """
    curTime = 0
    timeIndex = 0
    while curTime < maxTime:
        timeIndex += 1
        nextTime = logTime(1., timeIndex)
        if nextTime <= curTime: continue
        cycleTime = nextTime - curTime
        S.cycle(cycleTime)
        curTime += cycleTime

        yield S, curTime
    return  # can't return values in generators.


class LogtimeCallback(object):
    def __init__(self, maxTime=None):
        self.callbacks = { }
        self.callbacks_allTimes = [ ]
        self.maxTime = maxTime
    def maxTimeSet(self, maxTime):
        """Sets the max time to at least this long
        """


    def addCallback(self, deltaT, function):
        deltaT = roundTime(deltaT)
        self.callbacks.setdefault(deltaT, []).append(function)
        if deltaT == 'all':
            self.callbacks_allTimes.append(function)


    def runCallbacks(self, **args):
        deltaT = args['deltaT']
        callbacks = self.callbacks.get(deltaT, ())
        for c in callbacks:
            c(**args)
        for c in self.callbacks_allTimes:
            c(**args)

    def runSystem(self, S):
        newS = runLogTime(S, self.maxTime, callback=self.runCallbacks)
        return newS
