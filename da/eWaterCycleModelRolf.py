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

        self.shortNames = ['f','g','p','n']
        
    def initial(self):
        self.landSurface.storUpp000005=self.landSurface.storUpp000005 + (mapnormal() * 0.001)

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
                             
            ## if currTimeStep.day == 1 or self.currentTimeStep() == 1:
            ##     if currTimeStep.month == 1 or self.currentTimeStep() == 1:
            ##         self.precipitationAcc  = pcr.scalar(0.0) 
            ##         self.irrGrossDemandAcc = pcr.scalar(0.0)
            ##         self.actualETAcc       = pcr.scalar(0.0)
            ##         self.directRunoffAcc   = pcr.scalar(0.0)
            ##         self.gwRechargeAcc     = pcr.scalar(0.0)
            ##         self.infiltrationAcc   = pcr.scalar(0.0)
            ##         self.runoffAcc         = pcr.scalar(0.0)
            ##         self.waterbalanceAcc   = pcr.scalar(0.0)
            ## self.precipitationAcc  = self.precipitationAcc  + precipitation
            ## self.irrGrossDemandAcc = self.irrGrossDemandAcc + irrGrossDemand
            ## self.actualETAcc       = self.actualETAcc + actualET
            ## self.directRunoffAcc   = self.directRunoffAcc + directRunoff      
            ## self.gwRechargeAcc     = self.gwRechargeAcc + gwRecharge
            ## self.infiltrationAcc   = self.infiltrationAcc + infiltration 
            ## self.runoffAcc         = self.runoffAcc + runoff
            ## self.waterbalanceAcc   = self.waterbalanceAcc + waterbalance

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

        self.shortNames = ['f','g','p','n']


        #write all state-outputs that differ between Ensemble Members

        ## dumpfile("month.obj", currTimeStep.month, str(self.currentSampleNumber()))
        ## dumpfile("day.obj", currTimeStep.day, str(self.currentSampleNumber()))
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
	  
	  self.report(self.routing.discharge, "q")
          self.report(self.landSurface.storUpp000005, "su")
	  ## prev_time = int(loadfile("time.obj", str(self.currentSampleNumber())))
	  ## timestep = self.currentTimeStep()
	  ## dayCur = int(loadfile("day.obj", str(self.currentSampleNumber())))
	  ## qmod = ifthenelse(self.TRes == 1, self.Qaccu/scalar(timestep-prev_time), self.Qaccu/scalar(dayCur))
	  ## self.report(qmod, "q_acc")
	  ## self.report(maptotal(self.orgPrec), "orPrec")
	  ## self.report(maptotal(self.meteo.precipitation), "Prec")


          #We are not changing the parameters now, so we don't need to output them.
	  #report(self.routing.Q, str(self.currentSampleNumber())+"/"+"Q"+".map")
	  #report(self.Resistance, str(self.currentSampleNumber())+"/"+"Res"+".map")
	  #report(self.AvgSlope, str(self.currentSampleNumber())+"/"+"Slope"+".map")

	  # Parameter report
	  #self.report(self.minSoilDepthAdjust, "wmin")
	  #self.report(self.maxSoilDepthAdjust, "wmax")
	  #self.report(self.degreeDayAdjust, "ddf")
	  #self.report(self.KSatAdjust, "ksat")
  	  #self.report(self.Theta50Adjust, "theta")
	  #self.report(self.recessionAdjust, "j")
	  #self.report(self.routingAdjust, "n")
	  #self.report(self.precBiasAdjust, "alpha1")
	  #self.report(self.precConvectAdjust, "alpha2")
	  ## self.report(self.precHeightAdjust, "alpha3")
          

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
        discharge = self.readmap("q")
        #Location of station 1 in grid-coordinates
        loc1x=2
        loc1y=5
        #Location of station 2 in grid-coordinates
        loc2x=2
        loc2y=5
        values = cellvalue(discharge, loc1x, loc1y)[0]
        values = values.append(cellvalue(discharge, loc2x, loc2y)[0])
        return values

    def setObservations(self):
        timestep = self.currentTimeStep()
        #observedData = readmap(generateNameT("obsAv", timestep))
        values = numpy.zeros(2)
        values[0] = GETFROMFILE(station1,timestep)
        values[1] = GETFROMFILE(station2,timestep)
        #values[0] = cellvalue(observedData, 1, 1)[0]
        
        # creating the observation matrix (nrObservations x nrSamples)
        # here without added noise
        observations = numpy.array([values,]*self.nrSamples()).transpose()

        # creating the covariance matrix (nrObservations x nrObservations)
        # here just random values
        covariance = numpy.cov(observations)
        
        self.setObservedMatrices(observations, covariance)

    def resume(self):
        beta = 0.7
        vec = self.getStateVector(self.currentSampleNumber())
        timestep = self.currentTimeStep()-1
        
        ## self.minSoilDepthAdjust = self.readmap("wmin")
        ## self.maxSoilDepthAdjust = self.readmap("wmax")
        ## self.degreeDayAdjust = self.readmap("ddf")
        ## self.KSatAdjust = self.readmap("ksat")
	## self.Theta50Adjust = self.readmap("theta")
        ## self.recessionAdjust = self.readmap("j")
        ## self.routingAdjust = self.readmap("n")
        ## self.precBiasAdjust = self.readmap("alpha1")
        ## self.precConvectAdjust = self.readmap("alpha2")
        ## self.precHeightAdjust = self.readmap("alpha3")
    
        ## tel = 0
        ## xlocs2 = loadfile("x2.obj", str(self.currentSampleNumber()))
        ## ylocs2 = loadfile("y2.obj", str(self.currentSampleNumber()))
        ## lenx = len(xlocs2)
        ## Qtel = loadfile("Qtel.obj", str(1))
        
	### Some parameters are turned off for updating.
        ## Qlen = len(Qtel)
        ## tel = Qlen
        ## self.minSoilDepthAdjust = scalar(1) #beta * self.minSoilDepthAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.maxSoilDepthAdjust = beta * self.maxSoilDepthAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.degreeDayAdjust = beta * self.degreeDayAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.KSatAdjust = beta * self.KSatAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.Theta50Adjust = scalar(1) #beta * self.Theta50Adjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.recessionAdjust = beta * self.recessionAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.routingAdjust = beta * self.routingAdjust + (1- beta)* scalar(10**vec[tel])
        ## tel += 1
        ## self.precBiasAdjust = scalar(0) #beta * self.precBiasAdjust + (1- beta)* scalar(vec[tel])
        ## tel += 1
        ## self.precConvectAdjust = beta * self.precConvectAdjust + (1- beta)* scalar(vec[tel])
        ## tel += 1
        ## self.precHeightAdjust = beta * self.precHeightAdjust + (1- beta)* scalar(vec[tel])
        
        # Recalculate parameter maps based on update
        ## for coverType in self.coverTypes:
	##   self.landSurface.landCoverObj[coverType].minSoilDepthFrac = min(self.minSoilDepthOrig[coverType] * self.minSoilDepthAdjust,0.9999)
	##   self.landSurface.landCoverObj[coverType].maxSoilDepthFrac = max(self.maxSoilDepthOrig[coverType] * self.maxSoilDepthAdjust,1.0001)
	##   self.landSurface.landCoverObj[coverType].degreeDayFactor = cover(self.degreeDayOrig[coverType] * self.degreeDayAdjust,0.0)
        ## self.landSurface.KSat1 = self.KSat1Orig * self.KSatAdjust
        ## self.landSurface.KSat2 = self.KSat2Orig * self.KSatAdjust
        ## self.landSurface.THEFF1_50 = self.THEFF1_50org * self.Theta50Adjust
        ## self.landSurface.THEFF2_50 = self.THEFF2_50org * self.Theta50Adjust
        ## self.groundwater.recessionCoeff = self.recessionOrig * self.recessionAdjust
        ## self.routing.manningsN = self.routingOrig * self.routingAdjust
        
   
        self.routing.discharge = self.readmap("q")
        ## self.Qaccu = self.readmap("q_acc")

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
        self.groundwater.storGroundwater = readmap(str(self.currentSampleNumber())+"/"+"sg"+".map")
        self.routing.channelStorage = readmap(str(self.currentSampleNumber())+"/"+"Qc"+".map")
        self.routing.avgDischarge = readmap(str(self.currentSampleNumber())+"/"+"Qa"+".map")
        ## self.routing.Q = readmap(str(self.currentSampleNumber())+"/"+"Q"+".map")
        ## self.Resistance = readmap(str(self.currentSampleNumber())+"/"+"Res"+".map")
        ## self.AvgSlope = readmap(str(self.currentSampleNumber())+"/"+"Slope"+".map")
        
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
        
        ##month_end = loadfile("Monthend.obj", "1")

	##self.Qaccu = ifthenelse(boolean(month_end), scalar(0), ifthenelse(self.TRes != 2, scalar(0), self.Qaccu))
    

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

    MCModel=MonteCarloFramework(dynamicModel, nrSamples=4)
    MCModel.setForkSamples(True, nrCPUs=2)
    ekfModel = EnsKalmanFilterFramework(MCModel)
    ekfModel.setFilterTimesteps([2,4])
    ekfModel.run()

    # end of logFile
    logProcess.end()

if __name__ == '__main__':
    sys.exit(main())

