from math import log10, floor

class RoundSig:

    def __init__(self, x, sig=2):
        print 'x: ',x,' sig:', sig, ' result: ', self.round_sig(x,sig) 

    def round_sig(self,x, sig):
        if x == 0.0 or x == 0 or sig is 0:
            return 'invalid combination'
        if x < 0.0:
            return round(x,sig-int(floor(log10((abs(x)))))-1)
        return round(x, sig-int(floor(log10((x))))-1)

