from openmm.app import *
import openmm
import numpy as np

class BasePotential(openmm.CustomExternalForce):
    """
    """
    def __init__(self, force):
        super().__init__(force)
    
    def CalcForce(self, pos):
        pass

    def CalcPotential(self, pos):
        pass

class DoubleWellPotential(BasePotential):
    """
    """
    def __init__(self, *kwargs):
        force = "(x^2-1)^2 + y^2"
        super().__init__(force)
    
    def CalcForce(self, pos):
        derivative = np.zeros_like(pos)
        derivative[:,0] = 2 * 2 * pos[:,0] * (pos[:,0] ** 2 - 1)
        derivative[:,1] = 2 * pos[:,1]

        return -derivative

class MullerBrown(BasePotential):
    def __init__(self):
        self.a = [-1, -1, -6.5, 0.7]
        self.b = [0, 0, 11, 0.6]
        self.c = [-10, -10, -6.5, 0.7]
        self.A = [-200, -100, -170, 15]
        self.x_bar = [1, 0, -0.5, -1]
        self.y_bar = [0, 0.5, 1.5, 1]
        for i in range(4):
            fmt = dict(a=self.a[i], b=self.b[i], c=self.c[i], A=self.A[i], x_bar=self.x_bar[i], y_bar=self.y_bar[i])
            if i == 0:
                self.force = '''{A} * exp({a} * (x - {x_bar})^2 + {b} * (x - {x_bar}) * (y - {y_bar}) + {c} * (y - {y_bar})^2)'''.format(**fmt)
            else:
                self.force += ''' + {A} * exp({a} * (x - {x_bar})^2 + {b} * (x - {x_bar}) * (y - {y_bar}) + {c} * (y - {y_bar})^2)'''.format(**fmt)

        super().__init__(self.force)