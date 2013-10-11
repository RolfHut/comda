#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster import *
from pcraster.framework import *

class Transformation():

    def __init__(self):
        pass
   
    def readmap(self, name, style=1):
        """
        Read map data from disk.
    
        `name`
        Name used as filename. Use filename with less than eight characters
        and without extension. File extension for dynamic models is ".map"
        in initial section and the 8.3 style format name in the dynamic
        section. File extensions will be appended automatically.
    
        .. todo::
    
        `style` argument is not used.
        """
        return self._readmapNew(name)
    
    def model2Observation():
        msg = "Class needs to implement 'model2Observation' method"
        raise NotImplementedError(msg)
        
    def observation2Model(self):
        msg = "Class needs to implement 'observation2Model' method"
        raise NotImplementedError(msg)
        

        

