#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster import *
from pcraster.framework import *

import numpy

class IdentityObservations(Observations):

    def __init__(self,observation_count):
        Observation.__init__(self)
        self.observation_count=observation_count
        
    def getCount():
        return self.observation_count        

    def getObservations(self):
        values = numpy.ones(self.observation_count)
        
        # creating the observation matrix (nrObservations x nrSamples)
        # here without added noise
        observations = numpy.array([values,]*self.nrSamples()).transpose()
        
        # creating the covariance matrix (nrObservations x nrObservations)
        # here just random values
        covariance = numpy.ones((5, 5))

        return observations, covariance
