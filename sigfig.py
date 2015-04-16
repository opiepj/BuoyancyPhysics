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

test = RoundSig(0.0232)
test = RoundSig(0.0232,1)
test = RoundSig(0.02323,4)
test = RoundSig(-0.0232,1)
test = RoundSig(122232.23,1)
test = RoundSig(100232,4)
test = RoundSig(50.20322,3)
test = RoundSig(50.20322,0)
test = RoundSig(0,3)
test = RoundSig(0.0,3)
test = RoundSig(-999902.1,1)
test = RoundSig(-999902.1,2)
test = RoundSig(-999902.1,3)
test = RoundSig(-99820.1,4)
