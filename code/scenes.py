#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# This file is part of SimulationTeachingElan, a python code used for teaching at Elan Inria.
#
# Copyright 2020 Mickael Ly <mickael.ly@inria.fr> (Elan / Inria - Université Grenoble Alpes)
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
from math import *

from graphics import *
from dynamics import *
from geom import *

from graphics.sphere3D_renderable import *
from geom.sphere3D import *
from graphics.cube3D_renderable import *
from geom.cube3D import *
from dynamics.dynamic_sphere import *

def indexedTest(viewer):
    """
    @brief Demonstration for a basic static rendering
           Renders a simple square
    """

    # Indexed square
    positions = np.array([0., 0.,   # x0, y0
                          1., 0.,   # x1, y1
                          0., 1.,   # x2, y2
                          1., 1.],  # x3, y3
                         np.float64)
    colours = np.array([1., 0., 0.,  # (r, g, b) for vertex 0
                        0., 0., 1.,  # (r, g, b) for vertex 1
                        0., 1., 0.,  # ...
                        1., 1., 1.]) # ...
    indices = np.array([0, 1, 2,   # First triangle composed by vertices 0, 1 and 2
                        1, 2, 3])  # Second triangle composed by vertices 1, 2 and 3

    # Create the object
    squareMesh = Mesh2D(positions, indices, colours)
    # Create the correspondung GPU object
    squareMeshRenderable = Mesh2DRenderable(squareMesh)
    # Add it to the list of objects to render
    viewer.addRenderable(squareMeshRenderable)

def dynamicTest(viewer):
    """
    @brief Demonstration for a basic dynamic rendering
           Renders a simple square, moved by a dummy dynamic system
    """

    # Indexed square
    positions = np.array([0., 0.,   # x0, y0
                          1., 0.,   # x1, y1
                          0., 1.,   # x2, y2
                          1., 1.],  # x3, y3
                         np.float64)
    colours = np.array([1., 0., 0.,  # (r, g, b) for vertex 0
                        0., 0., 1.,  # (r, g, b) for vertex 1
                        0., 1., 0.,  # ...
                        1., 1., 1.]) # ...
    indices = np.array([0, 1, 2,   # First triangle composed by vertices 0, 1 and 2
                        1, 2, 3])  # Second triangle composed by vertices 1, 2 and 3

    # Create the object
    squareMesh = Mesh2D(positions, indices, colours)
    # Create the correspondung GPU object
    squareMeshRenderable = Mesh2DRenderable(squareMesh)
    # Add it to the list of objects to render
    viewer.addRenderable(squareMeshRenderable)

    # Create a dynamic system
    dyn = DummyDynamicSystem(squareMesh)
    # And add it to the viewer
    # Each frame will perform a call to the 'step' method of the viewer
    viewer.addDynamicSystem(dyn)



def rodTest(viewer):

    """
    @brief Demonstration for a rendering of a rod object
           Specific case, as a rod is essentialy a line, we
           need to generate a mesh over it to git it a thickness
           + demonstration of the scaling matrix for the rendering
    """
    positions = np.array([-1., 1.,
                          -1., 0.,
                          -0.5, -0.25],
                         np.float64)
    colours = np.array([1., 0., 0.,
                        0., 1., 0.,
                        0., 0., 1.])

    rod = Rod2D(positions, colours)

    rodRenderable = Rod2DRenderable(rod, thickness = 0.005)
    viewer.addRenderable(rodRenderable)

    positionsScaled = np.array([0., 1.,
                                0., 0.,
                                0.5, -0.25],
                               np.float64)
    rodScaled = Rod2D(positionsScaled, colours)

    rodRenderableScaled = Rod2DRenderable(rodScaled, thickness = 0.005)
    rodRenderableScaled.modelMatrix[0, 0] = 2.   # scale in X
    rodRenderableScaled.modelMatrix[1, 1] = 0.75 # scale in Y
    viewer.addRenderable(rodRenderableScaled)

def sphere_test(viewer):
    pi=4*atan(1.0)
    # Ground
    lar = 2.0
    lon = 7
    positions = np.array([-lar, -lon,   # x0, y0
                          +lar, -lon,   # x1, y1
                          -lar, +lon,   # x2, y2
                          +lar, +lon],  # x3, y3
                         np.float64)
    s1, s2, s3 = 219/255, 203/255, 174/255
    colours = np.array([s1, s2, s3,
                        s1, s2, s3,
                        s1, s2, s3,
                        s1, s2, s3]) # ...
    indices = np.array([0, 1, 2,   # First triangle composed by vertices 0, 1 and 2
                        1, 2, 3])  # Second triangle composed by vertices 1, 2 and 3
    ground_mesh = Mesh2D(positions, indices, colours)
    ground_renderable = Mesh2DRenderable(ground_mesh)
    viewer.addRenderable(ground_renderable)

    m1, m2 = 110/255, 51/255
    colours = np.array(
                [m1, m2, 0.,  # (r, g, b) for vertex 0
                 m1, m2, 0,
                 m1, m2, 0,
                 m1, m2, 0,
                 m1, m2, 0,
                 m1, m2, 0,
                 m1, m2, 0,
                 m1, m2, 0])
    poutre1_mesh = Cube3D(center=[-lar-0.1, 0, 0.1], lengths=[0.2, 14, 0.2])
    poutre1_renderable = Cube3DRenderable(poutre1_mesh, colours)
    viewer.addRenderable(poutre1_renderable)
    poutre2_mesh = Cube3D(center=[0, lon+0.1, 0.1], lengths=[4.4, 0.2, 0.2])
    poutre2_renderable = Cube3DRenderable(poutre2_mesh, colours)
    viewer.addRenderable(poutre2_renderable)
    poutre3_mesh = Cube3D(center=[lar+0.1, 0, 0.1], lengths=[0.2, 14, 0.2])
    poutre3_renderable = Cube3DRenderable(poutre3_mesh, colours)
    viewer.addRenderable(poutre3_renderable)
    poutre4_mesh = Cube3D(center=[0, -lon-0.1, 0.1], lengths=[4.4, 0.2, 0.2])
    poutre4_renderable = Cube3DRenderable(poutre4_mesh, colours)
    viewer.addRenderable(poutre4_renderable)

    # Sphere
    pos = np.array([0.0, -0.7, 1.0])
    theta = 0.#pi/4. # angle de rotation autour de l'axe z
    print(theta)
    c = cos(theta)
    s = sin(theta)
    rot = np.array([[c  ,-s  , 0.0],
                    [s  , c  , 0.0],
                    [0.0, 0.0, 1.0]])

    g1, g2, g3, g4 = 110/255, 120/255, 130/255, 140/255
    colors = np.array( [g1, g1, g1,  # (r, g, b) for (+,+,+)
                        g2, g2, g2,  # (r, g, b) for (+,+,-)
                        g3, g3, g3,  # (r, g, b) for (+,-,+)
                        g4, g4, g4,
                        g1, g1, g1,
                        g2, g2, g2,
                        g3, g3, g3,
                        g4, g4, g4])

    sphere = Sphere3D(center=pos, rotation=rot, radius=0.1, mass=1.0, velocity=[0, 5, 1], angular_velocity=[0., 0., 0.], lar=lar, lon=lon)
    sphere_renderable = Sphere3DRenderable(sphere, colors)
    viewer.addRenderable(sphere_renderable)

    sphere2 = Sphere3D(center=[0.0, 0.0, 1.], radius=0.1, mass=1.0, velocity=[0,0,0], angular_velocity=[0., 6.0, 0.], lar=lar , lon=lon)
    sphere_renderable2 = Sphere3DRenderable(sphere2, colors)
    viewer.addRenderable(sphere_renderable2)

    all_object = [sphere]
    all_object = [sphere, sphere2]

    # Create a dynamic system
    v0 = np.array([1.0]) # vitesse de translation et rotation ??
    dyn = DynamicSphere(all_object, 0.005)
    viewer.addDynamicSystem(dyn)
