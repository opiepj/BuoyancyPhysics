"""
    ocean.py
    Patrick Opie
    Ocean script that will be used for our buoyancy simulation

"""

import math, cmath, vector, fft
from math import sin, cos, pi, sqrt, fabs, exp
from random import gauss

class ocean:
    
    def __init__(self, geometry, g, N, Nplus1, A, w, length, fft):
        h_tilde = complex(N*N)
        h_tilde_slopex = complex(N*N)
        h_tilde_slopez = complex(N*N)
        h_tilde_dx = complex(N*N)
        h_tilde_dz = complex(N*N)
        vertices = vertex_ocean(Nplus1*Nplus1)
        indices = abs(Nplus1*Nplus1*10)
        self._g = g
        self._w = w
        self._A = A
        self._length = length
        self._fft = fft
        self._N = N
        self._Nplus1 = Nplus1
        m_prime = 0
        while m_prime < Nplus1:
            n_prime = 0
            while n_prime < Nplus1:
                index = m_prime*Nplus1 + n_prime
                
                htilde_0 = hTilde_0(n_prime,m_prime)
                htilde_0_conj = hTilde_0(-n_prime,m_prime)

                data = vertices[index]
                data['htilde_0'].x = htilde_0.x
                data['htilde_0'].y = htilde_0.y
                data['htilde_0_conj'].x = htilde_0_conj.x
                data['htilde_0_conj'].y = htilde_0_conj.y

                data['original_position'].x = data['vertex'].x = (n_prime - N/2.0) * (length/N)
                data['original_position'].y = data['vertex'].y = 0.0
                data['original_position'].x = data['vertex'].x = (m_prime - N/2.0) * (length/N)

                data['normal'].x = 0.0
                data['normal'].y = 1.0
                data['normal'].z = 0.0

                vertices[index] = data

                n_prime+=1
            m_prime+=1

        indices_count = 0
        m_prime = 0
        while m_prime < N:
            n_prime = 0
            while n_prime < N:
                if geometry:
                    indices_count +=1
                    indices[indices_count] = index
                    indices_count +=1
                    indices[indices_count] = index + 1
                    indices_count +=1
                    indices[indices_count] = index
                    indices_count +=1
                    indices[indices_count] = index + Nplus1
                    indices_count +=1
                    indices[indices_count] = index
                    indices_count +=1
                    indices[indices_count] = index + Nplus1 + 1
                    if n_prime == N - 1:
                        indices_count +=1
                        indices[indices_count] = index + 1
                        indices_count +=1
                        indices[indices_count] = index + Nplus1 + 1
                    if m_prime == N - 1:
                        indices_count +=1
                        indices[indices_count] = index + Nplus1
                        indices_count +=1
                        indices[indices_count] = index + Nplus1 + 1
                else:
                    indices_count +=1
                    indices[indices_count] = index
                    indices_count +=1
                    indices[indices_count] = index + Nplus1
                    indices_count +=1
                    indices[indices_count] = index + Nplus1 + 1
                    indices_count +=1
                    indices[indices_count] = index
                    indices_count +=1
                    indices[indices_count] = index + Nplus1 + 1
                    indices_count +=1
                    indices[indices_count] = index + 1

                n_prime+=1
            m_prime+=1
        self._vertices = vertices


    def dispersion(self, n_prime, m_prime):
        w_0 = 2.0 * pi / 120.0
        kx = pi * (2 * n_prime - self._N) / self._length
        ky = pi * (2 * m_prime - self._N) / self._length
        return (sqrt(self._g * sqrt(pow(kx,2.0) + pow(kz,2.0))) / w_0) * w_0

    def phillips(self, n_prime, m_prime):
        k = Vector(pi * (2*n_prime - self._N) / self._length, 0.0, pi * (2 * m_prime) / self._length)
        k_length = k.magnitude()
        if k_length < 0.000001:
            return 0.0
        k_length2 = pow(k_length,2.0)
        k_length4 = pow(k_length2,2.0)

        k_dot_w = k.dot(w)
        k_dot_w2 = pow(k_dot_w,2.0)

        w_length = w.magnitude()
        L = pow(w,2.0)/self._g
        L2 = pow(L,2.0)

        damping = 0.001
        l2 = L2 * pow(damping,2.0)

        return A * exp(-1.0/(k_length2 * L2)) / k_length4 * k_dot_w2 * exp(-k_length2 * l2)

    def hTilde_0(self, n_prime, m_prime):
        r = complex(gauss(10,2.5), gauss(10,2.5))
        return r * sqrt(self.phillips(n_prime, m_prime)/2.0)

    def hTilde(self, t, n_prime, m_prime):
        index = m_prime * self._Nplus1 + n_prime
        data = self._vertices[index]

        htilde0 = complex(data['htilde_0'].x, data['htilde_0'].y)
        htilde0_conj = complex(data['htilde_0_conj'].x, data['htilde_0'].y)

        omegat = dispersion(n_prime, m_prime) * t

        cos_ = cos(omegat)
        sin_ = sin(omegat)
        c0 = complex(cos_, sin_)
        c1 = complex(cos_, -sin_)

        return htilde0 * c0 + htilde0_conj*c1

    def h_D_n(self, x, t):
        h = complex(0.0,0.0)
        D = Vector(0.0,0.0,0.0)
        n = Vector(0.0,0.0,0.0)

        m_prime = 0
        while m_prime < self._N:
            n_prime = 0
            kz = 2.0 * pi * (m_prime - self._N/2.0)/self._length
            while n_prime < self._N:
                kx = 2.0 * pi * (n_prime - self._N/2.0)/self.length
                k = Vector(kx,0.0,kz)

                k_length = k.magnitude()
                k_dot_x = k.dot(x)

                c = complex(cos(k_dot_x), sin(k_dot_x))
                htilde_c = self.hTilde(t, n_prime, m_prime)

                h = h + htilde_c
                n = n + Vector(-kx * htilde_c.imag, 0.0, -kz*htilde_c.imag)
                
                if k_length  < 0.000001: continue
                D = D + Vector(kx/k_length * htilde_c.imag, 0.0, kz * htilde_c.imag)

                n_prime+=1
            m_prime+=1

        n = Vector(0.0, 1.0, 0.0) - n

        cvn = complex_vector_normal(h, D, n)

        return cvn

    def evaluate_waves(self, t):
        l = -1.0
        m_prime = 0
        while m_prime < self._N:
            while n_prime < self._N:
                index = m_prime * self._Nplus1 + n_prime
                data = self._vertices[index]
                a = data['vertex'].x
                b = 0.0
                c = data['vertex'].z
                x = Vector(a,b,c)

                h_d_n = self.h_D_n(x,t)
                h = h_d_n['waveheight']
                n = h_d_n['normal']
                D = h_d_n['displacement']

                data['vertex'].y = h.real
                data['vertex'].x = data['original_position'].x + l*D.x
                data['vertex'].z = data['original_position'].z + l*D.y

                data['normal'].x = n.x
                data['normal'].y = n.y
                data['normal'].z = n.z

                self._vertices[index] = data

                if n_prime == 0 and m_prime == 0:
                    data_prime = self._vertices[index + self._N + self._Nplus1 * self._N]
                    data_prime['vertex'].y = h.real
                    data_prime['vertex'].x = data_prime['original_position'].x + l*D.x
                    data_prime['vertex'].z = data_prime['original_position'].z + l*D.y

                    data_prime['normal'].x = n.x
                    data_prime['normal'].y = n.y
                    data_prime['normal'].z = n.z

                    self._vertices[index + self._N + self._Nplus1 * self._N] = data_prime

                if n_prime == 0:
                    data_prime = self._vertices[index + self._N]
                    data_prime['vertex'].y = h.real
                    data_prime['vertex'].x = data_prime['original_position'].x + l*D.x
                    data_prime['vertex'].z = data_prime['original_position'].z + l*D.y

                    data_prime['normal'].x = n.x
                    data_prime['normal'].y = n.y
                    data_prime['normal'].z = n.z

                    self._vertices[index + self._N] = data_prime

                if m_prime == 0:
                    data_prime = self._vertices[index + self._N * self._Nplus1]
                    data_prime['vertex'].y = h.real
                    data_prime['vertex'].x = data_prime['original_position'].x + l*D.x
                    data_prime['vertex'].z = data_prime['original_position'].z + l*D.y

                    data_prime['normal'].x = n.x
                    data_prime['normal'].y = n.y
                    data_prime['normal'].z = n.z

                    self._vertices[index + self._N * self.Nplus1] = data_prime

                n_prime+=1

            m_prime+=1
            
class vertex_ocean:

    def __init__(self, N):
        vertices = []
        i = 0
        while i < N:
            vertex = Vector(0.0,0.0,0.0)
            normal = Vector(0.0,0.0,0.0)
            htilde_0 = Vector(0.0,0.0,0.0)
            htilde_0_conj = Vector(0.0*j,0.0*j,0.0*j)
            original_position = Vector(0.0,0.0,0.0)
            data = {
                    'vertex':vertex,
                    'normal':normal,
                    'htilde_0':htilde_0,
                    'htilde_0_conj':htilde_0_conj, 
                    'original_position':original_position
                    }
            vertices.append(data)
        self = vertices

class complex_vector_normal:

    def __init__(self, height, displacement, normal):
        h = height # complex
        D = displacement # 2D Vector
        n = normal # 3D Vector
        data = {
                'waveheight':h,
                'displacement':d,
                'normal':n
                }
        self = data

static_test = ocean(False, -9.8, 1, 2, Vector(1.0,0.0,1.0), cFFT(fft)) 
