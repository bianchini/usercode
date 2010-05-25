#!/usr/bin/env python

import commands
import re
import os

command = "cmsBatch.py -o Data10May -q 1nd 40 Demo_cfg.py -r /castor/cern.ch/user/b/bianchi/CMSSW356/root/Data05May/PFAnalysis_Data07May.root -p analyzeElectrons"

print command

os.system(command)
