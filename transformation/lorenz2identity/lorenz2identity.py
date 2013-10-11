#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster import *
from pcraster.framework import *

from numpy import *

class Lorenz2identity(Transformation):
  def __init__(self):
    Transformation.__init__(self)

  def model2Observation(modelledData):
    values = numpy.zeros(1)
    values[0] = cellvalue(modelledData, 0, 0)[0]
    return values


  def observation2Model(observationData):
     ''' normally something like this:
     oldMap=self.readmap("x")
     oldMap(observationLocations)=observationData
     self.report(oldMap,"x")
     '''
     #this time, we have only 1 pixel
     return numpy2pcr(observationData)
      
 