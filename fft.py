"""
    fft.py
    Patrick Opie
    Fast Fourier Transform class

"""
import math, cmath, vector
from math import log, pi, cos, sin

class cFFT:

    def __init__(self, N, T):
        c = [[0] * 2]*2
        c[0] = 0
        c[1] = 0
        log_2_N = log(N)/log(2)
        _reversed = [0] * N
        T = [[0 for x in range(log_2_N)] for x in range(log_2_N)]
        for i in xrange(len(_reversed)):
            _reversed[i] = self.reverse(i)
        pow2 = 1
        for i in xrange(len(T)):
            T[i] = complex(pow2,pow2)
            for j in xrange(len(T)):
                T[i][j] = self.t(j, pow2 * 2)
            pow2 *= 2

        c[0] = complex(N,N)
        c[1] = complex(N,N)
        self._which = 0
        self._c = c
        self._N = N
        self._log_2_N = log_2_N
        self._T = T
        self._reversed = _reversed

# TO DO BITWISE OPERATION
    def reverse(self, i):
        for i in self._log_2_N:
            return 0

    def t(self, x, N):
        return complex(cos(2*pi * x / N), sin(2*pi * x / N))

    def fft(self, inp, outp, stride, offset):
        for i in xrange(N):
            c[self._which][i] = inp[self._reversed[i] * stride + offset]

        loops = N >> 1
        size = 1 << 1
        size_over_2 = 1
        w = 0

        i = 1
        while i < self._log_2_N:
            which ^= 1
            for j in xrange(loops):
                for k in xrange(size_over_2):
                    c[which][size * j + k] = c[which^1][size * j + k] + c[which^1][size * j + size_over_2 + k] * self._T[w][k]

                k = size_over_2
                while k < size:
                    c[which][size * j + k] = c[which^1][size * j - size_over_2 + k] - c[which^1][size * j + k] * self._T[w][k - size_over_two]
                    k += 1

            loops >>= 1
            size <<= 1
            size_over_2 <<= 1
            w += 1

            i += 1

        i = 0
        for i in xrange(self._N):
            outp[i * stride + offset] = c[which][i]









        

