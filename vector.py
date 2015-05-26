""" 
    vector.py
    Patrick Opie 
    3D Vector class
"""
import operator, math

class Vector(object):
    __slots__ = ['x', 'y', 'z']

    def __init__(self, x_or_triple, y = None, z = None):
        if y == None:
            self.x = x_or_triple[0]
            self.y = x_or_triple[1]
            self.z = x_or_triple[2]
        else:
            self.x = x_or_triple
            self.y = y
            self.z = z

    def __len__(self):
        return 3

    def __getitem__(self, key):
        if key == 0:
            return self.x
        elif key == 1:
            return self.y
        elif key == 2:
            return self.z
        else:
            raise IndexError("Invalid subscript "+str(key)+" to Vector")

    def __setitem__(self, key, value):
        if key == 0:
            self.x = value
        elif key == 1:
            self.y = value
        elif key == 2:
            self.z = value
        else:
            raise IndexError("Invalid subscript "+str(key)+" to Vector")

# String representaion (for debugging)
    def __repr__(self):
        return 'Vector(%s, %s, %s)' % (self.x, self.y, self.z)
    
    # Comparison
    def __eq__(self, other):
        if hasattr(other, "__getitem__") and len(other) == 3:
            return self.x == other[0] and self.y == other[1] and self.z == other[2]
        else:
            return False
    
    def __ne__(self, other):
        if hasattr(other, "__getitem__") and len(other) == 3:
            return self.x != other[0] or self.y != other[1] or self.z != other[2]
        else:
            return True
 
    def __nonzero__(self):
        return self.x or self.y or self.z

# Generic operator handlers
    def _o2(self, other, f):
        "Any two-operator operation where the left operand is a Vector"
        if isinstance(other, Vector):
            return Vector(f(self.x, other.x),
                         f(self.y, other.y),
                         f(self.z, other.z))
        elif (hasattr(other, "__getitem__")):
            return Vector(f(self.x, other[0]),
                         f(self.y, other[1]),
                         f(self.z, other[2]))
        else:
            return Vector(f(self.x, other),
                         f(self.y, other),
                         f(self.z, other))
 
    def _r_o2(self, other, f):
        "Any two-operator operation where the right operand is a Vector"
        if (hasattr(other, "__getitem__")):
            return Vector(f(other[0], self.x),
                         f(other[1], self.y),
                         f(other[2], self.z))
        else:
            return Vector(f(other, self.x),
                         f(other, self.y),
                         f(other, self.z))
 
    def _io(self, other, f):
        "inplace operator"
        if (hasattr(other, "__getitem__")):
            self.x = f(self.x, other[0])
            self.y = f(self.y, other[1])
            self.z = f(self.z, other[2])
        else:
            self.x = f(self.x, other)
            self.y = f(self.y, other)
            self.z = f(self.z, other)
        return self

    # Multiplication
    def __mul__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x*other.x, self.y*other.y, self.z*other.z)
        if (hasattr(other, "__getitem__")):
            return Vector(self.x*other[0], self.y*other[1], self.z*other[2])
        else:
            return Vector(self.x*other, self.y*other, self.z*other)
    __rmul__ = __mul__
    
    def __imul__(self, other):
        if isinstance(other, Vector):
            self.x *= other.x
            self.y *= other.y
            self.z *= other.z
        elif (hasattr(other, "__getitem__")):
            self.x *= other[0]
            self.y *= other[1]
            self.z *= other[2]
        else:
            self.x *= other
            self.y *= other
            self.z *= other
        return self
    
        
    # Division
    def __div__(self, other):
        return self._o2(other, operator.div)
    def __rdiv__(self, other):
        return self._r_o2(other, operator.div)
    def __idiv__(self, other):
        return self._io(other, operator.div)
 
 
    def __floordiv__(self, other):
        return self._o2(other, operator.floordiv)
    def __rfloordiv__(self, other):
        return self._r_o2(other, operator.floordiv)
    def __ifloordiv__(self, other):
        return self._io(other, operator.floordiv)
 
    def __truediv__(self, other):
        return self._o2(other, operator.truediv)
    def __rtruediv__(self, other):
        return self._r_o2(other, operator.truediv)
    def __itruediv__(self, other):
        return self._io(other, operator.floordiv)
 
    
    # Addition
    def __add__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
        elif hasattr(other, "__getitem__"):
            return Vector(self.x + other[0], self.y + other[1], self.z + other[2])
        else:
            return Vector(self.x + other, self.y + other, self.z + other)
    __radd__ = __add__
    
    def __iadd__(self, other):
        if isinstance(other, Vector):
            self.x += other.x
            self.y += other.y
            self.z += other.z
        elif hasattr(other, "__getitem__"):
            self.x += other[0]
            self.y += other[1]
            self.z += other[2]
        else:
            self.x += other
            self.y += other
            self.z += other
        return self
 
    # Subtraction
    def __sub__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
        elif (hasattr(other, "__getitem__")):
            return Vector(self.x - other[0], self.y - other[1], self.z - other[2])
        else:
            return Vector(self.x - other, self.y - other, self.z - other)
    def __rsub__(self, other):
        if isinstance(other, Vector):
            return Vector(other.x - self.x, other.y - self.y, other.z - self.z)
        if (hasattr(other, "__getitem__")):
            return Vector(other[0] - self.x, other[1] - self.y, other[2] - self.z)
        else:
            return Vector(other - self.x, other - self.y, other - self.z)
    def __isub__(self, other):
        if isinstance(other, Vector):
            self.x -= other.x
            self.y -= other.y
            self.z -= other.z
        elif (hasattr(other, "__getitem__")):
            self.x -= other[0]
            self.y -= other[1]
            self.z -= other[2]
        else:
            self.x -= other
            self.y -= other
            self.z -= other
        return self

    def magnitude(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def dot(self, other):
        return float(self.x*other[0] + self.y*other[1] + self.z*other[2])

    def cross(self, other):
        if isinstance(other, Vector):
            return Vector((self.y * other.z) - (self.z * other.y), (self.z*other.x)-(self.x*other.z), (self.x*other.y)-(self.y*other.x))
        else:
            return Vector(self.y*other[2] - self.z*other[1], self.z*other[0] - self.x*other[2], self.x*other[1] - self.y*other[0])
     
if __name__ == "__main__":
 
    import unittest
    import pickle
 
    ####################################################################
    class UnitTestVector(unittest.TestCase):
        def setUp(self):
            pass
        
        def testCreationAndAccess(self):
            v = Vector(111,222,333)
            self.assert_(v.x == 111 and v.y == 222 and v.z == 333)
            v.x = 333
            v[1] = 444
            v.z = 555
            self.assert_(v[0] == 333 and v[1] == 444 and v[2] == 555)
 
        def testMath(self):
            v = Vector(111,222,333)
            self.assertEqual(v + 1, Vector(112,223,334))
            self.assert_(v - 2 == [109,220,331])
            self.assert_(v * 3 == (333,666,999))
            self.assert_(v / 2.0 == Vector(55.5, 111, 166.5))
            self.assert_(v / 2 == (55, 111, 166))

        def testVectorOps(self):
            v1 = Vector(10,20,30)
            v2 = Vector(30,20,10)
            self.assert_(v1.cross(v2) == Vector(-400,800,-400))
            self.assert_(v2.cross(v1) == Vector(400,-800,400))
            self.assert_(v1.dot(v2) == 1000)
####################################################################
    unittest.main()