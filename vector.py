import math
from math import sqrt

class Vector3D(namedtuple('Vector3D',('x','y','z'))):
    __slots__ = ()

    def __abs__(self):
        return type(self)(abs(self.x),abs(self.y),abs(self.z))

    def __init__(self):
        return type(self)(int(self.x), int(self.y), int(self.z))

    def __length__(self):
        return sqrt(pow(self.x,2.0) + pow(self.y,2.0) + pow(self.z,2.0))

    def __add__(self,v):
        return type(self)(self.x + v.x, self.y + v.y, self.z + v.z)

    def __sub__(self,v):
        return type(self)(self.x - v.x, self.y - v.y, self.z - v.z)

    def __mul__(self,scalar):
        return type(self)(self.x * scalar, self.y * scalar, self.z * scalar) 

    def __div__(self,scalar):
        return type(self)(self.x / scalar, self.y / scalar, self.z / scalar)

    def dot_product(self,v):
        return ((self.x * v.x) + (self.y * v.y) + (self.z * v.z)) 

    def cross_product(self,v):
        return type(self)((self.y*v.z)-(self.z*v.y),(self.z*v.x)-(self.x*v.z),(self.x*v.y)-(self.y*v.x))
    
