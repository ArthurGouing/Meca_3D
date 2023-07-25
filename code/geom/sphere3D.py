    #!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# This file is part of SimulationTeachingElan, a python code used for teaching at Elan Inria.
#
# Copyright 2022 Thibaut Metivet <thibaut.metivet@inria.fr> (Elan / Inria - Université Grenoble Alpes)
# SimulationTeachingElan is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SimulationTeachingElan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SimulationTeachingElan.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
from math import sqrt, pi
from .rigidBody3D import RigidBody3D

## Class defining a 3D sphere
class Sphere3D(RigidBody3D):
    def __init__(self, radius=1.0, mass=1.0,
                 center=[0, 0, 0], rotation=np.eye(3, 3),
                 velocity=[0, 0, 0], angular_velocity=[0.0, 0, 0], lar=2, lon=7):

        ## Constructor
        # @param center           1-D Numpy array with 3D position of the center
        # @param velocity         1-D Numpy array with 3D velocity of the center
        # @param rotation         2-D Numpy array for rotation matrix
        # @param angular_velocity 1-D Numpy array with 3D w in the world Reference
        # @param radius           float: radius of the sphere
        # @param mass             float: mass of the sphere

        super().__init__(center, rotation)

        # Sphere property
        self.radius = radius
        self.mass = mass
        self.inertia_moment = 2./5. * mass * radius

        # Velocity parameter
        self.w_R0 = angular_velocity
        self.w_R1 = np.dot(self.rotation.transpose(), self.w_R0)
        self.V = velocity

        # Force and moment
        fz =  -9.81*mass
        self.force  = [0.0, 0.0, fz]
        self.moment = [0.0, 0.0, 0.0]

        # Initialize matrix M for computation
        M = np.zeros([6, 6])
        for i in range(0, 6):
            M[i][i] = 1./self.mass if i < 3 else 1./self.inertia_moment
        self.M_1 = M

        # self.contact contient tous les objets avec lesquelles self est en contact
        self.contact = list()
        # self.points of contacts contient tous les points de contacts pour chaque objet, triée comme
        self.points_of_contacts = list()
        self.force_contact = np.array([0.0, 0.0, 0.0])
        self.moment_contact = np.array([0., 0., 0.])

        self.young_modulus_norm = 200*(10**9)/(1-0.3**2)
        self.deq = (3 * 9.81 * self.mass / (4 * self.young_modulus_norm * sqrt(self.radius))) ** (2/3)

        self.previous_pos = [0, 0, 0]

        self.lar = lar
        self.lon = lon


    def __str__(self):
        return str(self.center)

    def get_center_radius(self):
        return (self.center, self.radius)

    def set_object_list(self, object_list):
        obj_list = object_list.copy()
        for i, o in enumerate(obj_list):
            if o == self:
                obj_list.pop(i)
        self.object_list = obj_list

    def in_contact(self, dt):
        """ Fonction qui recalcul les points de contact de l'objet """
        contact_points   = list()
        contact_distance = list()
        contact_list     = list()

        for obj in self.object_list:
            (center_obj, radius_obj) = obj.get_center_radius()
            distance = np.linalg.norm(self.center - center_obj)
            if distance <= self.radius + radius_obj:
                q = self.center + (center_obj - self.center) * self.radius/distance
                # Update contact list
                contact_list.append(obj)
                contact_points.append(q)
                # pour chaque objet on garde le point de contact le plus proche
                contact_points.append(q)

        if self.center[2] < self.radius:
            contact_points.append(abs(self.center[2] - self.radius)*np.array([0, 0, 1.]))
            contact_list.append("sol")

        self.contact = contact_list
        self.points_of_contacts = contact_points
        self.is_contact = True
        f = np.array([0., 0., 0.])
        m = np.array([0., 0., 0.])
        for obj, q  in zip(contact_list, contact_points):
            print(self.contact_force(obj, q, dt))
            [f, m] = np.add([f, m], self.contact_force(obj, q, dt))
        self.force_contact = f
        self.moment_contact = m
        return


    def advance(self, dt):
        # Initialize dq and force
        self.w_R1 =  np.dot(self.rotation.transpose(), self.w_R0)
        dq = np.array([*self.V, *self.w_R1], np.float64) 
        f =  self.force
        f += self.force_contact
        f += self.air_friction()
        print("force dans le advance", f)

        moment= np.add(self.moment, self.moment_contact)

        # Create matrix to compute dq_n+1
        Mf_R1 = np.dot(np.transpose(self.rotation), np.array(moment, np.float64))
        F = np.array([*f, *Mf_R1], np.float64)
        E = self.exp_mat(-dt, self.w_R1)
        L = np.block([[np.eye(3)      , np.zeros([3,3])],
                      [np.zeros([3,3]), E             ]])
        H = np.dot(self.M_1, L)

        # Compute dq_n+1 (velocity and anguar velocity)
        dq = np.dot(L, dq) + dt * np.dot(H, F)

        # Actualise the Position and Rotation from the computed velocity
        self.previous_pos = self.center - dt * dq[:3]
        self.center += dt*dq[:3]
        self.rotation = np.dot(self.rotation.transpose(), self.exp_mat(dt, dq[3:])).transpose()

        # Actualise the velocity
        if (self.center[0] <= -self.lar+self.radius or self.center[0] >= self.lar-self.radius):
            self.V[0], self.V[1], self.V[2] = 0, dq[1], dq[2]
            self.w_R1[0], self.w_R1[1], self.w_R1[2] = dq[3], 0, dq[5]
        elif (self.center[1] <= -self.lon+self.radius or self.center[1] >= self.lon-self.radius):
            self.V[0], self.V[1], self.V[2] = dq[0], 0, dq[2]
            self.w_R1[0], self.w_R1[1], self.w_R1[2] = 0, dq[4], dq[5]
        else:
            self.V    = dq[:3]
            self.w_R1 = dq[3:]
        self.w_R0 = np.dot(self.rotation, self.w_R1)

        Energie = self.mass*(1./2. * np.linalg.norm(self.V)**2 + 9.81*(self.center[2] - self.radius))

    def contact_force(self, obj, q, dt): # q est le point de contact
        if obj == "sol":
            direction = np.array([0., 0.0, 1.0])
            d = abs(self.radius-self.center[2])
            maxD = self.deq + abs(1.37* self.V[2]*10**(-6))
            mu=2.0

        else: # else c'est une sphere
            direction = (self.center - obj.center) / np.linalg.norm(obj.center-self.center)
            print(obj.radius + self.radius, self.center - obj.center, np.linalg.norm(self.center - obj.center))
            d = (obj.radius + self.radius) - np.linalg.norm(self.center - obj.center)
            maxD = self.deq + 1 * np.linalg.norm(self.V) * 10**(-6)
            mu = 1.0
        d = min(d, maxD)
        f = 4./3. * self.young_modulus_norm * sqrt(self.radius) * pow(d, 3./2.)
        # force tangentielle
        ds = self.calcul_ds(obj, dt)
        ft = 4./3. * self.young_modulus_norm * sqrt(self.radius) * pow(np.linalg.norm(ds), 3./2.)
        ft = min(ft, mu*f)
        dir_ft = ds/np.linalg.norm(ds)
        m = self.calcul_moment_contact(obj, dt, ft)
        f = f*direction #+ ft*dir_ft
        return [f, m]

    def air_friction(self):
        cx = 0.07
        rho = self.mass / (4*pi*(self.radius**3)/3)
        s = pi * (self.radius**2)
        f = - np.multiply(cx * rho * s * np.linalg.norm(self.V), self.V)
        if self.center[2] < self.radius:
            return f
        else:
            return 0

    def calcul_ds(self, obj, dt):
        if obj =="sol":
            ub = np.zeros(3)
            d  = np.array([.0, .0, -1.])
            ua = self.V + self.radius * np.cross(self.w_R1, d)
        else :
            d = (self.center-obj.center)/np.linalg.norm(self.center-obj.center)
            db = -d
            ua = self.V + self.radius * np.cross(self.w_R1, d)
            ub = obj.V + obj.radius * np.cross(obj.w_R1, db)
        uc = ub - ua
        ds = dt * uc
        print(np.dot(ds, d))
        ds = ds-np.dot(ds, d)*d
        return ds

    def calcul_moment_contact(self, obj, dt, ft):
        ds = self.calcul_ds(obj, dt)
        if obj == "sol":
            d  = np.array([.0, .0, -1.])
        else:
            d = (self.center-obj.center)/np.linalg.norm(self.center-obj.center)
        print(ds)
        dir_F = np.zeros(3) if np.linalg.norm(ds)==0. else ds/np.linalg.norm(ds)
        F = ft * dir_F
        m = np.cross(self.radius*d, F)
        return m
