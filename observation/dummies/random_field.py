#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PCRaster import *
from PCRaster.Framework import *

import random
from PIL import Image

imgx = 800
imgy = 600
image = Image.new("RGB", (imgx, imgy))

maxIt = 1000000 # number of pixels to draw
# maxIt = 1000
size = 30
xa = -size
xb = size
ya = -size
yb = size

delta = float(10) # Prandtl number
r = float(28)
b = float(8) / 3
h = 1e-3 # time step

def Lorenz(x, y, z):
    dx_dt = delta * (y - x)
    dy_dt = r * x - y - x * z
    dz_dt = x * y - b * z
    x += dx_dt * h
    y += dy_dt * h
    z += dz_dt * h
    return (x, y, z)

class LorenzModel(DynamicModel):

    def __init__(self):
        DynamicModel.__init__(self)
        setclone("clone.map")

    def initial(self):
        self.x = random.random() * size * 2 - 1
        self.y = random.random() * size * 2 - 1
        self.z = random.random() * size * 2 - 1
        
    def dynamic(self):
        (self.x, self.y, self.z) = Lorenz(self.x, self.y, self.z)
        xi = int((imgx - 1) * (self.x - xa) / (xb - xa))
        yi = int((imgy - 1) * (self.y - ya) / (yb - ya))
        if xi >=0 and xi < imgx and yi >= 0 and yi < imgy:
            image.putpixel((xi, yi), (255, 255, 255))

myModel = LorenzModel()
dynamicModel = DynamicFramework(myModel, lastTimeStep=maxIt, firstTimestep=1)
dynamicModel.run()
image.save("Lorenz_Attractor.png", "PNG")

