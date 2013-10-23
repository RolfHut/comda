#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import gc

from pcraster.framework import *
import pcraster as pcr

# paths such that the eWaterCycle model built by Edwin will be used
absolutePATH = os.environ["HOME"]+"/"
eWaterCycleModelPath = absolutePATH + "eWaterCycle/model"
sys.path.insert(0, eWaterCycleModelPath)

import virtualOS as vos
import iniFile
import logProcess
import spinUp
import currTimeStep
import meteo
import landSurface
import groundwater
import routing

class PCRGlobWB(DynamicModel, MonteCarloModel, EnKfModel):

    def __init__(self, iniStorage = None):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        EnKfModel.__init__(self)
        setclone(iniItems.cloneMap);

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

        # short name for every land cover type (needed for file name)
        self.shortNames = ['f','g','p','n']

    def premcloop(self):

        # set initial conditions for all MC members:
        self.landSurface.getICs(iniItems,spinUp.iniLandSurface)
        self.groundwater.getICs(iniItems,spinUp.iniGroundwater)
        self.routing.getICs(iniItems,spinUp.iniRouting)
        
    def initial(self):

	  # factor for perturbing the precipitation
	  self.precBiasAdjust = mapnormal() * 0.03

	  # factor for perturbing the initial storGroundwater
	  self.groundwater.storGroundwater = self.groundwater.storGroundwater * (mapnormal()*0.5+1)
	  self.groundwater.storGroundwater = pcr.max(0.,self.groundwater.storGroundwater)

    def dynamic(self):

        currTimeStep.update(currTimeStep.startTime,self.currentTimeStep())

        self.meteo.update(currTimeStep)                                            

        # perturb the precipitation 
        self.meteo.precipitation = self.meteo.precipitation * \
             pcr.min(pcr.max((1 + self.precBiasAdjust),0.01),2.0)
        # PS: Please also make sure that precipitation >= 0     

        # THE MAIN COMPONENTS OF THE MODEL:
        self.landSurface.update(self.meteo,self.groundwater,self.routing,currTimeStep)      
        self.groundwater.update(self.landSurface,self.routing,currTimeStep)
        self.routing.update(self.landSurface,self.groundwater,currTimeStep,self.meteo)

        #write all state-outputs that differ between Ensemble Members

        if 1: #self.currentTimeStep() == self.ReportTime[WriteTime]:
	  idx = 0
	  for coverType in self.landSurface.coverTypes:
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
              self.landSurface.landCoverObj[coverType].storUpp000005, 
              str(self.currentSampleNumber())+"/"+"su"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].storUpp005030, 
              str(self.currentSampleNumber())+"/"+"sm"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].storLow030150, 
              str(self.currentSampleNumber())+"/"+"sl"+str(self.shortNames[idx])+".map")
            report(\
              self.landSurface.landCoverObj[coverType].interflow, 
              str(self.currentSampleNumber())+"/"+"qi"+str(self.shortNames[idx])+".map")
            idx = idx + 1
          
	  report(self.groundwater.storGroundwater, str(self.currentSampleNumber())+"/"+"sg"+".map")
	  self.report(self.groundwater.storGroundwater, "sg")

	  report(self.routing.channelStorage, str(self.currentSampleNumber())+"/"+"Qc"+".map")
	  report(self.routing.avgDischarge, str(self.currentSampleNumber())+"/"+"Qa"+".map")
	  
	  report(self.routing.timestepsToAvgDischarge, str(self.currentSampleNumber())+"/"+"t"+".map")
	  
	  report(self.precBiasAdjust, str(self.currentSampleNumber())+"/"+"p"+".map")
	  
        self.report(self.routing.discharge, "q")

    def postmcloop(self):
        #~ names = ["q","sg"]
        #~ mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        #~ percentiles = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        #~ mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())
        pass

    def setState(self):

        # average groundwater storage
        grAvg = pcr.areaaverage(readmap(str(self.currentSampleNumber())+"/"+"sg"+".map"),pcr.nominal(self.landmask))

        values=numpy.zeros(3)
        
        # average groundwater storage
        
        discharge = self.readmap("q")
        #Location of station 1 in grid-coordinates # Lobith
        loc0x = 7 #373
        loc0y = 2 #77
        #Location of station 2 in grid-coordinates # Borgharen
        loc1x = 5 #372
        loc1y = 3 #78
        values[0] = cellvalue(discharge, loc0y, loc0x)[0]
        values[1] = cellvalue(discharge, loc1y, loc1x)[0]
        values[2] = cellvalue(grAvg, 5, 5)[0]
        print "model states"
        print values
        return values

    def setObservations(self):
        timestep = self.currentTimeStep()
        Qerr = float(0.001)

        values = numpy.zeros(2)
        #~ values = numpy.zeros(3)

        f = open(iniItems.globalOptions['observationDir'] + "/LOBITH_data_integrated_start_1950.txt",'r')
        b = f.readlines()
        values[0] = float(b[0+18629+timestep-1].split("\t")[1].split("\n")[0])
        f.close()
        f = open(iniItems.globalOptions['observationDir'] + "/BORGHAREN_data_integrated_start_1950.txt",'r')
        b = f.readlines()
        values[1] = float(b[0+18629+timestep-1].split("\t")[1].split("\n")[0])
        f.close()
        
        #~ values[0] = 10000000000.
        #~ values[1] = 10000000000.
        #~ values[2] = 0.
        
        ## Create covariance matrix
        Qlen = len(values)
        covariance = numpy.eye(Qlen, Qlen, dtype=float)
        for Qplace in range(Qlen):
			covariance[Qplace, Qplace] = (Qerr*values[Qplace])**2
			#~ covariance[Qplace, Qplace] = 1#(Qerr*values[Qplace])**2
        
        # creating the observation matrix (nrObservations x nrSamples)
        # here without added noise
        observations = numpy.array([values,]*self.nrSamples()).transpose()
        for Qplace in range(Qlen):
            observations[Qplace,:] += numpy.random.normal(0,covariance[Qplace,Qplace]**0.5,self.nrSamples())
            observations[Qplace,:]  = numpy.maximum(observations[Qplace,:], 0.01)     
        #~ print observations

        self.setObservedMatrices(observations, covariance)

    def resume(self):
        vec = self.getStateVector(self.currentSampleNumber())
        print vec
        
        self.groundwater.storGroundwater = readmap(str(self.currentSampleNumber())+"/"+"sg"+".map")
        self.groundwater.storGroundwater = self.groundwater.storGroundwater * (vec[2]/pcr.areaaverage(self.groundwater.storGroundwater,self.landmask))
        
        self.routing.discharge = self.readmap("q")

        self.precBiasAdjust = readmap(str(self.currentSampleNumber())+"/"+"p"+".map")

        idx = 0
        for coverType in self.landSurface.coverTypes:
            self.landSurface.landCoverObj[coverType].interceptStor =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"si"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].snowCoverSWE  =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sc"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].snowFreeWater =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sf"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].topWaterLayer =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"st"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].storUpp000005 =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"su"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].storUpp005030 =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sm"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].storLow030120 =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"sl"+str(self.shortNames[idx])+".map"), 0.0)
            self.landSurface.landCoverObj[coverType].interflow =\
              cover(readmap(str(self.currentSampleNumber())+"/"+"qi"+str(self.shortNames[idx])+".map"), 0.0)
            idx = idx + 1
        self.routing.channelStorage = readmap(str(self.currentSampleNumber())+"/"+"Qc"+".map")
        self.routing.avgDischarge = readmap(str(self.currentSampleNumber())+"/"+"Qa"+".map")
        
        self.routing.timestepsToAvgDischarge = readmap(str(self.currentSampleNumber())+"/"+"t"+".map")
        #~ self.routing.WaterBodies.waterBodyStorage = readmap(str(self.currentSampleNumber())+"/"+"rs"+".map")
        #~ self.routing.WaterBodies.avgInflow = readmap(str(self.currentSampleNumber())+"/"+"in"+".map")
        #~ self.routing.WaterBodies.avgOutflow = readmap(str(self.currentSampleNumber())+"/"+"out"+".map")

	#~ try:
  	  #~ self.routing.WaterBodies.waterBodyTyp = readmap(str(self.currentSampleNumber())+"/"+"ty"+".map")
	  #~ self.routing.WaterBodies.fracWat = readmap(str(self.currentSampleNumber())+"/"+"rfr"+".map")
	  #~ self.routing.WaterBodies.waterBodyIds = readmap(str(self.currentSampleNumber())+"/"+"ri"+".map")
	  #~ self.routing.WaterBodies.waterBodyArea = readmap(str(self.currentSampleNumber())+"/"+"ra"+".map")
	  #~ self.routing.WaterBodies.waterBodyOut = readmap(str(self.currentSampleNumber())+"/"+"ro"+".map")
	  #~ self.routing.WaterBodies.waterBodyCap = readmap(str(self.currentSampleNumber())+"/"+"rc"+".map")
	#~ except:
	  #~ foo = 0        
        
# Global Objects:
iniItems = iniFile.IniList()               # content of configuration file)
logProcess = logProcess.LogProcess()       # object for logging
spinUp = spinUp.SpinUp()                   # object for spinUp
currTimeStep = currTimeStep.CurrTimeStep() # timeStep info: year, month, day, doy, hour, etc

def main():
    
    # reading argument and reading iniFile
    iniItems.readIniFile()    

    # start logFile
    logProcess.start(iniItems.logFileDir,iniItems.iniFileName)
    
    # spinningUp
    noSpinUps = int(iniItems.globalOptions['maxSpinUpsInYears'])
    if noSpinUps > 0:
        logProcess.write('Spin-Up #Total Years: '+str(noSpinUps),"Yes")
        spinUp.setupConvergence(iniItems) 
        spinUpRun = 0 ; convCondition = False
        while spinUpRun < noSpinUps and convCondition == False:
            spinUpRun += 1
            currTimeStep.getStartEndTimeStepsForSpinUp(
                    iniItems.globalOptions['startTime'],
                    spinUpRun, noSpinUps)
            logProcess.write('Spin-Up Run No. '+str(spinUpRun),"Yes")
            myModel = PCRGlobWB()
            dynamicModel = DynamicFramework(myModel,currTimeStep.nrOfTimeSteps)
            dynamicModel.setQuiet(True)
            dynamicModel.run()
            spinUp.getEndStates(myModel,currTimeStep)
            convCondition = spinUp.checkConvergence(myModel)
            
            convRep = 'Delta SoilStorage = %.2f percent ; SpinUp No. %i of %i' \
                    %(spinUp.convSoilSto, spinUpRun, noSpinUps)
            print(convRep); logProcess.write(convRep)
            convRep = 'Delta GwatStorage = %.2f percent' \
                    %(spinUp.convGwatSto)
            print(convRep); logProcess.write(convRep)
            convRep = 'Delta ChanStorage = %.2f percent' \
                    %(spinUp.convChanSto)
            print(convRep); logProcess.write(convRep)
            convRep = 'Delta TotlStorage = %.2f percent' \
                    %(spinUp.convTotlSto)
            print(convRep); logProcess.write(convRep)
    #
    # Running the model (excluding DA scheme)
    currTimeStep.getStartEndTimeSteps(iniItems.globalOptions['startTime'],
                                      iniItems.globalOptions['endTime'])
    logProcess.write('Transient simulation run started.',"Yes")
    myModel = PCRGlobWB()
    dynamicModel = DynamicFramework(myModel,currTimeStep.nrOfTimeSteps)
    dynamicModel.setQuiet(True)

    MCModel=MonteCarloFramework(dynamicModel, nrSamples=75)
    MCModel.setForkSamples(True, nrCPUs=15)
    ekfModel = EnsKalmanFilterFramework(MCModel)
    #~ FilterTimesteps = range(1,44,2)
    #~ FilterTimesteps = [30,40]
    FilterTimesteps = range(75,2000,25)
    ekfModel.setFilterTimesteps(FilterTimesteps)
    #~ ekfModel.setFilterTimesteps([5,15, 30, 45])
    ekfModel.run()

    # end of logFile
    logProcess.end()

if __name__ == '__main__':
    sys.exit(main())

