import time

class StopWatch():

    def __init__(self):
        self._elapsedtime = 0.0
        self._start = 0.0
        self._setTime(self._elapsedtime)

    def _setTime(self, elap):
        # Set the time string to Minutes:Seconds:Hundreths
        minutes = int(elap/60)
        seconds = int(elap - minutes*60.0)
        hseconds = int((elap - minutes*60.0 - seconds)*100)
        timestr = ('%02d:%02d:%02d' % (minutes, seconds, hseconds))
        print timestr

    def _update(self):
        print 'update called' 
        self._elapsedtime = time.time() - self._start
        self._setTime(self._elapsedtime)
        self._update()
        self._timer = self.after(50, self._update)

    def Start(self):
        self._start = time.time() - self._elapsedtime
        self._update()

    def Stop(self):
        self.after_cancel(self._timer)
        self._elapsedtime = time.time() - self._start
        self._setTime(self._elapsedtime)

    def Reset(self):
        self._start = time.time()
        self._elapsedtime = 0.0
