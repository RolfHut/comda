#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import gc

from pcraster.framework import *
import pcraster as pcr


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
    def premcloop(self):

        # set initial conditions for all MC members:
        self.landSurface.getICs(iniItems,spinUp.iniLandSurface)
        self.groundwater.getICs(iniItems,spinUp.iniGroundwater)
        self.routing.getICs(iniItems,spinUp.iniRouting)

        # only needed if the SpinUp is introduced: 
        if spinUp.noSpinUps != None:
            spinUp.getIniStates(self)
        
    def initial(self):
        self.landSurface.storUpp000005=self.landSurface.storUpp000005 + (mapnormal() * 0.001)
        self.test=mapnormal()

 
 

    def dynamic(self):

        if iniItems.globalOptions['debugWaterBalance'] == "True":
            preTotStor = self.landSurface.interceptStor  +\
                         self.landSurface.snowFreeWater  +\
                         self.landSurface.snowCoverSWE   +\
                         self.landSurface.topWaterLayer  +\
                         self.landSurface.storUpp000005  +\
                         self.landSurface.storUpp005030  +\
                         self.landSurface.storLow030150  +\
                         self.groundwater.storGroundwater
            preTotStor = ifthen(self.landmask, preTotStor)

        #
        # THE MAIN COMPONENTS OF THE MODEL:
        currTimeStep.update(currTimeStep.startTime,self.currentTimeStep())
        self.meteo.update(currTimeStep)                                            
        self.landSurface.update(self.meteo,self.groundwater,self.routing,currTimeStep)      
        self.groundwater.update(self.landSurface,self.routing,currTimeStep)
        self.routing.update(self.landSurface,self.groundwater,currTimeStep,self.meteo)

        if iniItems.globalOptions['debugWaterBalance'] == "True":
            endTotStor = self.landSurface.interceptStor   +\
                         self.landSurface.snowFreeWater   +\
                         self.landSurface.snowCoverSWE    +\
                         self.landSurface.topWaterLayer   +\
                         self.landSurface.storUpp000005   +\
                         self.landSurface.storUpp005030   +\
                         self.landSurface.storLow030150   +\
                         self.groundwater.storGroundwater
            endTotStor = ifthen(self.landmask, endTotStor)             

            precipitation  = ifthen(self.landmask,\
                                    self.meteo.precipitation)
            irrGrossDemand = ifthen(self.landmask,\
                                    self.landSurface.irrGrossDemand)
            actualET       = ifthen(self.landmask,\
                                    self.landSurface.actualET)
            directRunoff   = ifthen(self.landmask,\
                                    self.landSurface.directRunoff)
            gwRecharge     = ifthen(self.landmask,\
                                    self.landSurface.gwRecharge)
            infiltration   = ifthen(self.landmask,\
                                    self.landSurface.infiltration)
            runoff         = ifthen(self.landmask,self.routing.runoff)
            gwAbstractions = ifthen(self.landmask,self.groundwater.gwAbstractions)
            waterbalance   = ifthen(self.landmask,
                             preTotStor - endTotStor +\
                             precipitation + irrGrossDemand -\
                             actualET - runoff -gwAbstractions)

            vos.waterBalanceCheck([precipitation,irrGrossDemand],\
                                  [runoff,actualET,gwAbstractions],\
                                  [preTotStor],\
                                  [endTotStor],\
                                  'all Modules except river/routing',\
                                   True,\
                                   currTimeStep.fulldate,threshold=1e-4)
                             
            if currTimeStep.day == 1 or self.currentTimeStep() == 1:
                if currTimeStep.month == 1 or self.currentTimeStep() == 1:
                    self.precipitationAcc  = pcr.scalar(0.0) 
                    self.irrGrossDemandAcc = pcr.scalar(0.0)
                    self.actualETAcc       = pcr.scalar(0.0)
                    self.directRunoffAcc   = pcr.scalar(0.0)
                    self.gwRechargeAcc     = pcr.scalar(0.0)
                    self.infiltrationAcc   = pcr.scalar(0.0)
                    self.runoffAcc         = pcr.scalar(0.0)
                    self.waterbalanceAcc   = pcr.scalar(0.0)
            self.precipitationAcc  = self.precipitationAcc  + precipitation
            self.irrGrossDemandAcc = self.irrGrossDemandAcc + irrGrossDemand
            self.actualETAcc       = self.actualETAcc + actualET
            self.directRunoffAcc   = self.directRunoffAcc + directRunoff      
            self.gwRechargeAcc     = self.gwRechargeAcc + gwRecharge
            self.infiltrationAcc   = self.infiltrationAcc + infiltration 
            self.runoffAcc         = self.runoffAcc + runoff
            self.waterbalanceAcc   = self.waterbalanceAcc + waterbalance

            if currTimeStep.month == 12:
                if currTimeStep.day == 31:
                    #~ if currTimeStep.year%5 == 0: self.getEndStates(currTimeStep,iniItems) # writing the endStates every 5 year.
                    if currTimeStep.year%1 == 0: self.getEndStates(currTimeStep,iniItems)    # writing the endStates every 1 year.
                    # reporting the endStates at the end of the Year:
                    for var in ['precipitation','irrGrossDemand',\
                                'actualET','runoff',\
                                'directRunoff',\
                                'infiltration','gwRecharge',\
                                'waterbalance']:
                        volume = vos.getMapVolume(\
                                    self.__getattribute__(var + 'Acc'),\
                                    self.routing.cellArea)
                        msg = 'Accumulated %s days 1 to %i in %i = %e km3'\
                            % (var,int(currTimeStep.doy),\
                                   int(currTimeStep.year),volume/1e9)
                        logProcess.write(msg) ; print(msg)
                    logProcess.write('')

        #write outputs for MC
        self.report(self.landSurface.storUpp000005, "sm")
        self.report(self.test, "test")
        self.report(self.routing.discharge, "d")
        

    def getEndStates(self,currTimeStep,iniItems): # THIS PART SHOULD MOVE TO the spinUp module

        self.iniLandSurface = {}
        self.iniGroundwater = {}
        self.iniRouting     = {} 
        
        for coverType in self.landSurface.coverTypes:
            self.iniLandSurface[coverType] = {}
            self.iniLandSurface[coverType]['interceptStorIni'] = \
              self.landSurface.landCoverObj[coverType].interceptStor
            self.iniLandSurface[coverType]['snowCoverSWEIni' ] = \
              self.landSurface.landCoverObj[coverType].snowCoverSWE
            self.iniLandSurface[coverType]['snowFreeWaterIni'] = \
              self.landSurface.landCoverObj[coverType].snowFreeWater
            self.iniLandSurface[coverType]['topWaterLayerIni'] = \
              self.landSurface.landCoverObj[coverType].topWaterLayer
            self.iniLandSurface[coverType]['storUpp000005Ini'] = \
              self.landSurface.landCoverObj[coverType].storUpp000005
            self.iniLandSurface[coverType]['storUpp005030Ini'] = \
              self.landSurface.landCoverObj[coverType].storUpp005030
            self.iniLandSurface[coverType]['storLow030150Ini'] = \
              self.landSurface.landCoverObj[coverType].storLow030150
            self.iniLandSurface[coverType]['interflowIni'    ] = \
              self.landSurface.landCoverObj[coverType].interflow
 
        self.iniGroundwater = {}
        self.iniGroundwater['storGroundwaterIni'] = \
           self.groundwater.storGroundwater

        self.iniRouting = {}  
        self.iniRouting['channelStorageIni'] = \
           self.routing.channelStorage
        self.iniRouting['avgDischargeIni'  ] = \
           self.routing.avgDischarge

        self.iniRouting['timestepsToAvgDischargeIni'] = \
           self.routing.timestepsToAvgDischarge

        self.iniRouting['waterBodyStorageIni'] = \
           self.routing.WaterBodies.waterBodyStorage
        self.iniRouting['avgInflowIni'] = \
           self.routing.WaterBodies.avgInflow
        self.iniRouting['avgOutflowIni'] = \
           self.routing.WaterBodies.avgOutflow

        self.writeEndStates(currTimeStep,iniItems)

    def writeEndStates(self,currTimeStep,iniItems): # THIS PART SHOULD MOVE TO the spinUp module

        for coverType in self.landSurface.coverTypes:
            for var in list(self.iniLandSurface[coverType]):
                vos.writePCRmapToDir(\
                 self.iniLandSurface[coverType][var],\
                 str(var)+"_"+coverType+"_"+
                 str(currTimeStep.fulldate)+".map",\
                 iniItems.endStateDir)

        for var in list(self.iniGroundwater):
            vos.writePCRmapToDir(\
             self.iniGroundwater[var],\
             str(var)+"_"+
             str(currTimeStep.fulldate)+".map",\
             iniItems.endStateDir)
        
        for var in list(self.iniRouting):
            vos.writePCRmapToDir(\
             self.iniRouting[var],\
             str(var)+"_"+
             str(currTimeStep.fulldate)+".map",\
             iniItems.endStateDir)

    def postmcloop(self):
        pass

    def setState(self):
        modelledData = self.readmap("sm")

        #load a map that contains the pixels with measurements

        #convert modeldata to estimates at the location of the measurements
        
        modelledAverageMap = areaaverage(modelledData, iniItems.cloneMap)
        self.report(modelledAverageMap, "modAv")
        values = numpy.zeros(1)
        values[0] = cellvalue(modelledAverageMap, 1, 1)[0]
        return values

    def setObservations(self):
        timestep = self.currentTimeStep()
        #observedData = readmap(generateNameT("obsAv", timestep))
        values = numpy.zeros(1)
        values[0] = 0.5
        #values[0] = cellvalue(observedData, 1, 1)[0]
        
        # creating the observation matrix (nrObservations x nrSamples)
        # here without added noise
        observations = numpy.array([values,]*self.nrSamples()).transpose()

        # creating the covariance matrix (nrObservations x nrObservations)
        # here just random values
        covariance = numpy.random.random((1, 1))

        self.setObservedMatrices(observations, covariance)

    def resume(self):
        vec = self.getStateVector(self.currentSampleNumber())
        modelledAverageMap = self.readmap("modAv")
        self.landSurface.storUpp000005 = pcr.scalar(vec[0])
    

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

    MCModel=MonteCarloFramework(dynamicModel, nrSamples=2)
    ekfModel = EnsKalmanFilterFramework(MCModel)
    ekfModel.setFilterTimesteps([2,20])
    ekfModel.run()

    # end of logFile
    logProcess.end()

if __name__ == '__main__':
    sys.exit(main())

