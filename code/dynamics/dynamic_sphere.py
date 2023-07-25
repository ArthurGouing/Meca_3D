from .abstract_dynamic_system import AbstractDynamicSystem
import numpy as np
from scipy.linalg import expm
from math import atan # enlever

# Dynamics class
class DynamicSphere(AbstractDynamicSystem):

    def __init__(self, object_list, dt):
        ## Constructor
        # @param self
        # @param mesh
        super().__init__()
        # Simulation parameters
        self.object = object_list
        self.t = 0.
        self.dt = dt

        for o in self.object:
            o.set_object_list(object_list)

    def step(self):
        print("")
        print(f"temps : {self.t}")
        self.t += self.dt
        for o in self.object:
            o.in_contact(self.dt)

        # Advance : calcul the new Xnp1 and Rnp1
        for o in self.object:
            o.advance(self.dt)
