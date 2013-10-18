#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script will be used for the DA's project with Niko. 

import os
import sys
import math
import gc
import datetime

from pickle import *

from pcraster.framework import *
from pcraster import *

import numpy as np

import virtualOS as vos
import iniFile
import logProcess
import spinUpTwoLayers as spinUp
import currTimeStep
import meteo
import landSurfaceTwoLayers as landSurface
import groundwater
import routingTwoLayers as routing

from nikoEnKF import *

import pcraster as pcr
 
class PCRGlobWB(DynamicModel, MonteCarloModel, EnKfModel):

    def __init__(self, iniStorage = None):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        EnKfModel.__init__(self)
        setclone(iniItems.cloneMap);

        cmd = 'cp '+str(iniItems.cloneMap)+' clone.map'
        os.system(cmd)
        os.system('pcrcalc id.map = "uniqueid(clone.map)"')
        
        ## Map with station ID
        self.ids = pcr.nominal(vos.readPCRmapClone(str(iniItems.globalOptions['observationDir'])+"ObsPoints.map", iniItems.cloneMap, "."))
        ## Number of stations
        self.ObsNum = int(cellvalue(mapmaximum(scalar(self.ids)),1,1)[0])
        ## Check for monthyl or daily observations
        self.TRes = pcr.nominal(vos.readPCRmapClone(str(iniItems.globalOptions['observationDir'])+"TResObs.map", iniItems.cloneMap, "."))
        
        Qxloc = []
        Qyloc = []
        for x in range(1,vos.getMapAttributes(iniItems.cloneMap,"rows")):
            for y in range(1,vos.getMapAttributes(iniItems.cloneMap,"cols")):
                if cellvalue(self.ids, x, y)[0] >= 1:
	                Qxloc.append(x)
	                Qyloc.append(y)
        self.Qxloc = Qxloc
        self.Qyloc = Qyloc        
        
        self.ReportTime = iniItems.dataAssimilationOptions['filterTimeSteps']
        self.ReportTime = map(int, self.ReportTime.split(','))
        self.ReportTime.append(0)

    def premcloop(self): 

        # initializing routing module
        # To get the "landmask", the routing module is initiated first.  
        self.routing = routing.Routing(iniItems)
        self.landmask = self.routing.landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(\
           iniItems.globalOptions['landmask'],
           iniItems.cloneMap,iniItems.tmpDir,iniItems.inputDir)	
        
        # initializing other modules
        self.meteo = meteo.Meteo(iniItems,\
                           self.landmask)
        self.landSurface = landSurface.LandSurface(iniItems,\
                           self.landmask)
        self.groundwater = groundwater.Groundwater(iniItems,\
                           self.landmask)
        
        # obtain the reference value(s) of the parameter(s) 
        # that we want to perturb
        self.refrecessionCoeff = self.groundwater.recessionCoeff

        # coverTypes (and their shortNames) that we will use  
        self.coverTypes = \
          iniItems.landSurfaceOptions['landCoverTypes'].split(",")
        self.shortNames = \
          iniItems.landSurfaceOptions['shortNames'].split(",")

        # Niko's way to identify coordinates of active/landmask cells: 
        self.mask = boolean(self.routing.lddMap)
        xlocs2 = []
        ylocs2 = []
        for x in range(1,vos.getMapAttributes(iniItems.cloneMap,"rows")):
            for y in range(1,vos.getMapAttributes(iniItems.cloneMap,"cols")):
                if cellvalue(self.mask, x, y)[1] == True:
	                xlocs2.append(x)
	                ylocs2.append(y)
        self.xlocs2 = xlocs2
        self.ylocs2 = ylocs2

        # set initial conditions:
        self.landSurface.getICs(iniItems)
        self.groundwater.getICs(iniItems)
        self.routing.getICs(iniItems)
        # set initial conditions
        self.minSoilDepthOrig = {}
        self.maxSoilDepthOrig = {}
        self.degreeDayOrig = {}
        for coverType in self.coverTypes:
	  self.minSoilDepthOrig[coverType] = self.landSurface.landCoverObj[coverType].minSoilDepthFrac
	  self.maxSoilDepthOrig[coverType] = self.landSurface.landCoverObj[coverType].maxSoilDepthFrac
	  self.degreeDayOrig[coverType] = self.landSurface.landCoverObj[coverType].degreeDayFactor
        self.KSat1Orig = self.landSurface.KSat1
        self.KSat2Orig = self.landSurface.KSat2 
        self.THEFF1_50org = self.landSurface.THEFF1_50
        self.THEFF2_50org = self.landSurface.THEFF2_50
        self.recessionOrig = self.groundwater.recessionCoeff
        self.routingOrig = self.routing.manningsN
        
    
        self.CatchTotal = catchmenttotal(self.routing.cellArea,self.routing.lddMap)
        
        self.Height = vos.readPCRmapClone(str(iniItems.globalOptions['inputDir'])+"elevation/dem30min.map", iniItems.cloneMap, ".")
        self.Slope = vos.readPCRmapClone(str(iniItems.globalOptions['inputDir'])+"elevation/slope30min.map", iniItems.cloneMap, ".")
	self.Aspect = scalar(aspect(self.Height))
	self.Aspect = cover(self.Aspect, 0.0)
	self.Slope = cover(self.Slope, 0.0)
        
    def initial(self): 
    
    	try:
	  cmd = 'mkdir '+str(self.currentSampleNumber())+'/stateVar'
	  os.system(cmd)
	except:
	  foo = 0
        dumpfile("time.obj", int(0),"1")
        dumpfile("WriteTime.obj", int(0),"1")
	folderName = iniItems.dataAssimilationOptions['folderName']
    
        ##perturb the parameter(s):
        if folderName == "None":
	  self.minSoilDepthAdjust = scalar(1) #mapnormal() * 0.1 + 1
	  self.maxSoilDepthAdjust = mapnormal() * 0.1 + 1
	  self.degreeDayAdjust = mapnormal() * 0.1 + 1
	  self.Theta50Adjust = scalar(1) #mapnormal() * 0.1 + 1
	  self.KSatAdjust = mapnormal() * 0.1 + 1
	  self.recessionAdjust = mapnormal() * 0.1 + 1
	  self.routingAdjust = mapnormal() * 0.03 + 1
	  self.precBiasAdjust = scalar(0) #mapnormal() * 0.03
	  self.precConvectAdjust = mapnormal() * 0.03
	  self.precHeightAdjust = mapnormal() * 0.03
	## Load parameter perturbation from previous run
        if folderName != "None":
	  self.minSoilDepthAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/wmin.map")
	  self.maxSoilDepthAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/wmax.map")
	  self.degreeDayAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/ddf.map")
	  self.Theta50Adjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/theta.map")
	  self.KSatAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/ksat.map")
	  self.recessionAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/j.map")
	  self.routingAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/n.map")
	  self.precBiasAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/alpha1.map")
	  self.precConvectAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/alpha2.map")
	  self.precHeightAdjust = readmap(folderName+"/"+str(self.currentSampleNumber())+"/alpha3.map")
    
        # Perturb maps
        for coverType in self.coverTypes:
	  self.landSurface.landCoverObj[coverType].minSoilDepthFrac = min(self.minSoilDepthOrig[coverType] * self.minSoilDepthAdjust,0.9999)
	  self.landSurface.landCoverObj[coverType].maxSoilDepthFrac = max(self.maxSoilDepthOrig[coverType]* self.maxSoilDepthAdjust,1.0001)
	  self.landSurface.landCoverObj[coverType].degreeDayFactor = self.degreeDayOrig[coverType] * self.degreeDayAdjust
        self.landSurface.KSat1 = self.KSat1Orig * self.KSatAdjust
        self.landSurface.KSat2 = self.KSat2Orig * self.KSatAdjust
        self.landSurface.THEFF1_50 = self.THEFF1_50org * self.Theta50Adjust
        self.landSurface.THEFF2_50 = self.THEFF2_50org * self.Theta50Adjust
        self.groundwater.recessionCoeff = self.recessionOrig * self.recessionAdjust
        self.routing.manningsN = self.routingOrig * self.routingAdjust
  
        self.Qaccu = scalar(0)

    def dynamic(self):

        ### Precipitation
        currTimeStep.update(currTimeStep.startTime,self.currentTimeStep())
	WriteTime = int(loadfile("WriteTime.obj", "1"))
        if currTimeStep.day == 1 or self.currentTimeStep() == 1:
	  ### Load winddistance
	  self.TravelDist = vos.netcdf2PCRobjCloneWindDist(str(iniItems.globalOptions['inputDir'])+"meteo/Climatology_winddis.nc", "Dist", currTimeStep.month, useDoy = "Yes")
	  ### Reset Qaccu
	  self.Qaccu = scalar(0)
          dumpfile("time.obj", self.currentTimeStep()-1,str(self.currentSampleNumber()))
    
        self.meteo.update(currTimeStep)

        ## Modify precipitation
        self.Cor1 = self.precHeightAdjust * self.Slope
	self.Cor2 = self.precConvectAdjust * self.TravelDist
	self.orgPrec = self.meteo.precipitation
        self.meteo.precipitation = self.meteo.precipitation * min(max((1 + self.precBiasAdjust + self.Cor1 + self.Cor2),0.01),2.0)

        
        self.landSurface.update(self.meteo,self.groundwater,self.routing,currTimeStep)
        print "Land"
        
        self.groundwater.update(self.landSurface,self.routing,currTimeStep)
        print "Groundwater"
        self.routing.update(self.landSurface,self.groundwater,currTimeStep,self.meteo)
        print "Routing"
        
 
	self.Qaccu += self.routing.discharge
                
        # writing the model states to disk 
        # - that will be re-used in the "resume" method:
        dumpfile("month.obj", currTimeStep.month, str(self.currentSampleNumber()))
        dumpfile("day.obj", currTimeStep.day, str(self.currentSampleNumber()))
        if self.currentTimeStep() == self.ReportTime[WriteTime]:
	  idx = 0
	  for coverType in self.coverTypes:
	    report(\
              self.landSurface.landCoverObj[coverType].interceptStor, 
              str(self.currentSampleNumber())+"/"+"si"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].snowCoverSWE, 
              str(self.currentSampleNumber())+"/"+"sc"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].snowFreeWater, 
              str(self.currentSampleNumber())+"/"+"sf"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].topWaterLayer, 
              str(self.currentSampleNumber())+"/"+"st"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].storUpp, 
              str(self.currentSampleNumber())+"/"+"su"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].storLow, 
              str(self.currentSampleNumber())+"/"+"sl"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].interflow, 
              str(self.currentSampleNumber())+"/"+"qi"+str(self.shortNames[idx])+".map")
            idx = idx + 1
	  report(self.groundwater.storGroundwater, str(self.currentSampleNumber())+"/"+"sg"+".map")
	  report(self.routing.channelStorage, str(self.currentSampleNumber())+"/"+"Qc"+".map")
	  report(self.routing.avgDischarge, str(self.currentSampleNumber())+"/"+"Qa"+".map")
	  
	  report(self.routing.timestepsToAvgDischarge, str(self.currentSampleNumber())+"/"+"t"+".map")
	  
	  report(self.routing.WaterBodies.waterBodyStorage, str(self.currentSampleNumber())+"/"+"rs"+".map")
	  report(self.routing.WaterBodies.avgInflow,str(self.currentSampleNumber())+"/"+"in"+".map")
	  report(self.routing.WaterBodies.avgOutflow,str(self.currentSampleNumber())+"/"+"out"+".map")
	  try:
	    report(self.routing.WaterBodies.waterBodyTyp,str(self.currentSampleNumber())+"/"+"ty"+".map")
	    report(self.routing.WaterBodies.fracWat,str(self.currentSampleNumber())+"/"+"rfr"+".map") 
	    report(self.routing.WaterBodies.waterBodyIds,str(self.currentSampleNumber())+"/"+"ri"+".map")
	    report(self.routing.WaterBodies.waterBodyArea,str(self.currentSampleNumber())+"/"+"ra"+".map")
	    report(self.routing.WaterBodies.waterBodyOut,str(self.currentSampleNumber())+"/"+"ro"+".map")
	    report(self.routing.WaterBodies.waterBodyCap,str(self.currentSampleNumber())+"/"+"rc"+".map")
	  except:
	    foo = 0
	  
	  report(self.routing.Q, str(self.currentSampleNumber())+"/"+"Q"+".map")
	  report(self.Resistance, str(self.currentSampleNumber())+"/"+"Res"+".map")
	  report(self.AvgSlope, str(self.currentSampleNumber())+"/"+"Slope"+".map")
	  # Parameter report
	  self.report(self.minSoilDepthAdjust, "wmin")
	  self.report(self.maxSoilDepthAdjust, "wmax")
	  self.report(self.degreeDayAdjust, "ddf")
	  self.report(self.KSatAdjust, "ksat")
  	  self.report(self.Theta50Adjust, "theta")
	  self.report(self.recessionAdjust, "j")
	  self.report(self.routingAdjust, "n")
	  self.report(self.precBiasAdjust, "alpha1")
	  self.report(self.precConvectAdjust, "alpha2")
	  self.report(self.precHeightAdjust, "alpha3")
	  self.report(self.routing.discharge, "q")
	  prev_time = int(loadfile("time.obj", str(self.currentSampleNumber())))
	  timestep = self.currentTimeStep()
	  dayCur = int(loadfile("day.obj", str(self.currentSampleNumber())))
	  qmod = ifthenelse(self.TRes == 1, self.Qaccu/scalar(timestep-prev_time), self.Qaccu/scalar(dayCur))
	  self.report(qmod, "q_acc")
	  self.report(maptotal(self.orgPrec), "orPrec")
	  self.report(maptotal(self.meteo.precipitation), "Prec")

    def setState(self):
        timestep = self.currentTimeStep()
        
        prev_time = int(loadfile("time.obj", str(self.currentSampleNumber())))

        qmod = self.readmap("q_acc")
        
        wmin = self.readmap("wmin")
        wmax = self.readmap("wmax")
        ddf = self.readmap("ddf")
        ksat = self.readmap("ksat")
        theta = self.readmap("theta")
        j = self.readmap("j")
        n = self.readmap("n")
        alpha1 = self.readmap("alpha1")
        alpha2 = self.readmap("alpha2")
        alpha3 = self.readmap("alpha3")


        # Principle (at this moment): If there is at least a observation, 
        #                             we get the model state. 

        # The following will be replaced by observation data:
        stats = self.ids                                     # station code
        year = int(iniItems.globalOptions['startTime'][0:4])
        month = int(iniItems.globalOptions['startTime'][5:7])
        day = int(iniItems.globalOptions['startTime'][8:10])
        dif = datetime.date(year, month, day) - datetime.date(1900, 01, 01)
        obs = scalar(0)
	dayCur = int(loadfile("day.obj", str(self.currentSampleNumber())))
        if self.currentSampleNumber() == 1:
	  ## Current timeStep
	  monthCur = int(loadfile("month.obj", str(self.currentSampleNumber())))
	  dayCur = int(loadfile("day.obj", str(self.currentSampleNumber())))
	  month_end = 0

	  
	  if (monthCur == 1 and dayCur == 31) or (monthCur == 2 and dayCur == 28) or (monthCur == 2 and dayCur == 29) or\
	    (monthCur == 3 and dayCur == 31) or (monthCur == 4 and dayCur == 30) or (monthCur == 5 and dayCur == 31) or\
	    (monthCur == 6 and dayCur == 30) or (monthCur == 7 and dayCur == 31) or (monthCur == 8 and dayCur == 31) or\
	    (monthCur == 9 and dayCur == 30) or (monthCur == 10 and dayCur == 31) or (monthCur == 11 and dayCur == 30) or (monthCur == 12 and dayCur == 31):
	      month_end = 1
	      
          dumpfile("Monthend.obj", month_end, str(self.currentSampleNumber()))
          for time in range(prev_time+1,timestep+1):
	    day_obs = vos.readPCRmapClone(generateNameT(str(iniItems.globalOptions['observationDir'])+"update/Qobs", time+dif.days), iniItems.cloneMap, ".")
	    obs = ifthenelse(day_obs > 0.0, obs+day_obs, obs)
	  print range(prev_time+1,timestep+1)
	  self.qobs = ifthenelse(self.TRes == 1, obs/scalar(timestep-prev_time), day_obs)
	  report(self.qobs, generateNameT(str(self.currentSampleNumber())+"/obs", timestep))
	      
        values = []
        Qtel = []
        Lim = 3.0
        if self.currentSampleNumber() == 1:
	  for s in range(len(self.Qxloc)):
	    val  = cellvalue(qmod, self.Qxloc[s],self.Qyloc[s])[0]
	    val2 = cellvalue(self.qobs, self.Qxloc[s],self.Qyloc[s])[0]
	    if val2 > 0.0 and val/val2 < Lim and val/val2 > 1/Lim:
	      values.append(val)
	      Qtel.append(s)
	  print Qtel
	  if len(values) == 0:
	    values.append(numpy.random.normal(0,1,1)[0])
	    Qtel.append(-999)
	  dumpfile("Qtel.obj", Qtel, str(self.currentSampleNumber()))
	else:
	  Qlist = loadfile("Qtel.obj", str(1))
	  for s in Qlist:
	    val  = cellvalue(qmod, self.Qxloc[s],self.Qyloc[s])[0]
	    values.append(val)
	    Qtel.append(s)  
        print len(Qtel)
        xlocs2 = self.xlocs2 ## All locations with data
        ylocs2 = self.ylocs2 ## All locations with data
        endx = len(xlocs2)
        self.endx = endx
        
        x = xlocs2[0]
        y = ylocs2[0]
        values.append(numpy.log10(cellvalue(wmin, x, y)[0]))
        values.append(numpy.log10(cellvalue(wmax, x, y)[0]))
        values.append(numpy.log10(cellvalue(ddf, x, y)[0]))
        values.append(numpy.log10(cellvalue(ksat, x, y)[0]))
        values.append(numpy.log10(cellvalue(theta, x, y)[0]))
        values.append(numpy.log10(cellvalue(j, x, y)[0]))
        values.append(numpy.log10(cellvalue(n, x, y)[0]))
        values.append(cellvalue(alpha1, x, y)[0])
        values.append(cellvalue(alpha2, x, y)[0])
        values.append(cellvalue(alpha3, x, y)[0])
        
        values2 = numpy.array(values)

        dumpfile("x2.obj", xlocs2, str(self.currentSampleNumber()))
        dumpfile("y2.obj", ylocs2, str(self.currentSampleNumber()))
        return values2

    def setObservations(self):
		
        timestep = self.currentTimeStep()
        values = []
        Per = 0.15
        Qerr = float(0.20)
        stats = self.ids                                   # station code
        obs = self.qobs

        Qtel = loadfile("Qtel.obj", str(1))
        Qlen = len(Qtel)
        covariance = numpy.eye(Qlen, Qlen, dtype=float)

        Qplace = 0
        for s in Qtel:
	  if s != -999:
	    val = cellvalue(obs, self.Qxloc[s],self.Qyloc[s])[0]
	    values.append(numpy.maximum(val,1.0))
	    covariance[Qplace, Qplace] = (Qerr*val)**2
	    Qplace += 1
	  else:
	    val = numpy.random.normal(0,1000,1)[0]
	    values.append(val)
	    covariance[Qplace, Qplace] = (1000)**2
	    Qplace += 1
        dumpfile("time.obj", timestep,"1")
        dumpfile("WriteTime.obj", int(loadfile("WriteTime.obj", "1"))+1, "1")
        values2 = numpy.array(values)
#        print values2
        observations = numpy.array([values2,]*self.nrSamples()).transpose()
        for i in range(Qplace):
            observations[i,:] += numpy.random.normal(0,covariance[i,i]**0.5,self.nrSamples())
            observations[i,:]  = numpy.maximum(observations[i,:], 0.01)            

        month_end = loadfile("Monthend.obj", "1")
        if month_end == 1:
	  print "End of month"
	print Qtel
            
        self.giveLenghts(Qlen, self.endx)
        self.setObservedMatrices(observations, covariance)
        return values

    def resume(self):
        beta = 0.7
        vec = self.getStateVector(self.currentSampleNumber())
        timestep = self.currentTimeStep()-1
        
        self.minSoilDepthAdjust = self.readmap("wmin")
        self.maxSoilDepthAdjust = self.readmap("wmax")
        self.degreeDayAdjust = self.readmap("ddf")
        self.KSatAdjust = self.readmap("ksat")
	self.Theta50Adjust = self.readmap("theta")
        self.recessionAdjust = self.readmap("j")
        self.routingAdjust = self.readmap("n")
        self.precBiasAdjust = self.readmap("alpha1")
        self.precConvectAdjust = self.readmap("alpha2")
        self.precHeightAdjust = self.readmap("alpha3")
    
        tel = 0
        xlocs2 = loadfile("x2.obj", str(self.currentSampleNumber()))
        ylocs2 = loadfile("y2.obj", str(self.currentSampleNumber()))
        lenx = len(xlocs2)
        Qtel = loadfile("Qtel.obj", str(1))
        
	### Some parameters are turned off for updating.
        Qlen = len(Qtel)
        tel = Qlen
        self.minSoilDepthAdjust = scalar(1) #beta * self.minSoilDepthAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.maxSoilDepthAdjust = beta * self.maxSoilDepthAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.degreeDayAdjust = beta * self.degreeDayAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.KSatAdjust = beta * self.KSatAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.Theta50Adjust = scalar(1) #beta * self.Theta50Adjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.recessionAdjust = beta * self.recessionAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.routingAdjust = beta * self.routingAdjust + (1- beta)* scalar(10**vec[tel])
        tel += 1
        self.precBiasAdjust = scalar(0) #beta * self.precBiasAdjust + (1- beta)* scalar(vec[tel])
        tel += 1
        self.precConvectAdjust = beta * self.precConvectAdjust + (1- beta)* scalar(vec[tel])
        tel += 1
        self.precHeightAdjust = beta * self.precHeightAdjust + (1- beta)* scalar(vec[tel])
        
        # Recalculate parameter maps based on update
        for coverType in self.coverTypes:
	  self.landSurface.landCoverObj[coverType].minSoilDepthFrac = min(self.minSoilDepthOrig[coverType] * self.minSoilDepthAdjust,0.9999)
	  self.landSurface.landCoverObj[coverType].maxSoilDepthFrac = max(self.maxSoilDepthOrig[coverType] * self.maxSoilDepthAdjust,1.0001)
	  self.landSurface.landCoverObj[coverType].degreeDayFactor = cover(self.degreeDayOrig[coverType] * self.degreeDayAdjust,0.0)
        self.landSurface.KSat1 = self.KSat1Orig * self.KSatAdjust
        self.landSurface.KSat2 = self.KSat2Orig * self.KSatAdjust
        self.landSurface.THEFF1_50 = self.THEFF1_50org * self.Theta50Adjust
        self.landSurface.THEFF2_50 = self.THEFF2_50org * self.Theta50Adjust
        self.groundwater.recessionCoeff = self.recessionOrig * self.recessionAdjust
        self.routing.manningsN = self.routingOrig * self.routingAdjust
        
   
        self.routing.discharge = self.readmap("q")
        self.Qaccu = self.readmap("q_acc")

        idx = 0
        for coverType in self.coverTypes:
            self.landSurface.landCoverObj[coverType].interceptStor =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"si"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].snowCoverSWE  =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sc"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].snowFreeWater =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sf"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].topWaterLayer =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"st"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].storUpp =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"su"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].storLow =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sl"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].interflow =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"qi"+str(self.shortNames[idx])+".map"), 0.0)
            idx = idx + 1
        self.groundwater.storGroundwater = readmap(str(self.currentSampleNumber())+"/"+"sg"+".map")
        self.routing.channelStorage = readmap(str(self.currentSampleNumber())+"/"+"Qc"+".map")
        self.routing.avgDischarge = readmap(str(self.currentSampleNumber())+"/"+"Qa"+".map")
        self.routing.Q = readmap(str(self.currentSampleNumber())+"/"+"Q"+".map")
        self.Resistance = readmap(str(self.currentSampleNumber())+"/"+"Res"+".map")
        self.AvgSlope = readmap(str(self.currentSampleNumber())+"/"+"Slope"+".map")
        
        self.routing.timestepsToAvgDischarge = readmap(str(self.currentSampleNumber())+"/"+"t"+".map")
        self.routing.WaterBodies.waterBodyStorage = readmap(str(self.currentSampleNumber())+"/"+"rs"+".map")
        self.routing.WaterBodies.avgInflow = readmap(str(self.currentSampleNumber())+"/"+"in"+".map")
        self.routing.WaterBodies.avgOutflow = readmap(str(self.currentSampleNumber())+"/"+"out"+".map")

	try:
  	  self.routing.WaterBodies.waterBodyTyp = readmap(str(self.currentSampleNumber())+"/"+"ty"+".map")
	  self.routing.WaterBodies.fracWat = readmap(str(self.currentSampleNumber())+"/"+"rfr"+".map")
	  self.routing.WaterBodies.waterBodyIds = readmap(str(self.currentSampleNumber())+"/"+"ri"+".map")
	  self.routing.WaterBodies.waterBodyArea = readmap(str(self.currentSampleNumber())+"/"+"ra"+".map")
	  self.routing.WaterBodies.waterBodyOut = readmap(str(self.currentSampleNumber())+"/"+"ro"+".map")
	  self.routing.WaterBodies.waterBodyCap = readmap(str(self.currentSampleNumber())+"/"+"rc"+".map")
	except:
	  foo = 0        
        
        month_end = loadfile("Monthend.obj", "1")

	self.Qaccu = ifthenelse(boolean(month_end), scalar(0), ifthenelse(self.TRes != 2, scalar(0), self.Qaccu))
                
     

# Global Objects:
iniItems = iniFile.IniList()               # content of configuration file)
logProcess = logProcess.LogProcess()       # object for logging
currTimeStep = currTimeStep.CurrTimeStep() # timeStep info: year, month, day, doy, hour, etc

def main():
    
    # reading argument and reading iniFile
    iniItems.readIniFile()    

    # start logFile
    logProcess.start(iniItems.logFileDir,iniItems.iniFileName)
    
    # Running the model 
    currTimeStep.getStartEndTimeSteps(iniItems.globalOptions['startTime'],
                                      iniItems.globalOptions['endTime'])
    logProcess.write('START.',"Yes")

    myModel = PCRGlobWB()
    dynamicModel = DynamicFramework(myModel,currTimeStep.nrOfTimeSteps)

    if iniItems.dataAssimilationOptions['method'] == "None":
        dynamicModel.run()
    else:    
        nrSamples = int(iniItems.dataAssimilationOptions['nrSamples'])
        mcModel = MonteCarloFramework(dynamicModel,nrSamples)
        mcModel.setForkSamples(True, nrCPUs=int(iniItems.dataAssimilationOptions['nrCores']))
    if iniItems.dataAssimilationOptions['method'] == "MonteCarlo":
        mcModel.run()
    if iniItems.dataAssimilationOptions['method'] == "EnKF":
        ekfModel = EnsKalmanFilterFramework(mcModel)
        filterTime = iniItems.dataAssimilationOptions['filterTimeSteps']
        filterTime = map(int, filterTime.split(','))
        ekfModel.setFilterTimesteps(filterTime)    #range(365,6900,30)
        ekfModel.run()
    
    # end of logFile
    logProcess.end()

if __name__ == '__main__':
    sys.exit(main())

