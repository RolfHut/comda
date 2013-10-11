#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster import *
from pcraster.framework import *

import random

class LorenzModel(DynamicModel):

    def __init__(self,delta,r,h,t,observations=None,transformations=None):
        DynamicModel.__init__(self)
        setclone("clone.map")
        self.delta
        self.r=r
        self.h=h
        self.t=t
        self.observations
        self.transformations

    def initial(self):
        self.x = random.random() * size * 2 - 1
        self.y = random.random() * size * 2 - 1
        self.z = random.random() * size * 2 - 1
        
    @staticmethod
    def lorenz(x, y, z, delta, r, h, t):
        dx_dt = delta * (y - x)
        dy_dt = r * x - y - x * z
        dz_dt = x * y - b * z
        x += dx_dt * h
        y += dy_dt * h
        z += dz_dt * h
        return (x, y, z)
        
    def dynamic(self):
        (self.x, self.y, self.z) = lorenz(self.x, self.y, self.z, self.delta, self.r, self.h, self.t)
        report(self.x,"x")

    def setState(self):
        modelledData = self.readmap("x")
        return self.transformations.model2Observation(modelledData)

    def setObservations(self):
        observations, covariance = self.observations.getObservations()
        self.setObservedMatrices(observations, covariance)

    def resume(self):
        vec = self.getStateVector(self.currentSampleNumber())
        self.x = self.transformations.observation2Model(vec)




