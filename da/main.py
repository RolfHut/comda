#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster import *
from pcraster.framework import *


from observation import *
from transformation import *

from models.lorenz import *

import random
from PIL import Image


delta = float(10) # Prandtl number
r = float(28)
b = float(8) / 3
h = 1e-3 # time step

observations=IdentityObservations(1)
transformations=Lorenz2identity()

myModel = LorenzModel(delta,r,h,t,observations,transformations)
dynamicModel = DynamicFramework(myModel, lastTimeStep=180, firstTimestep=1)
mcModel = MonteCarloFramework(dynamicModel, nrSamples=10)
ekfModel = EnsKalmanFilterFramework(mcModel)
ekfModel.setFilterTimesteps([70, 100, 150])
ekfModel.run()