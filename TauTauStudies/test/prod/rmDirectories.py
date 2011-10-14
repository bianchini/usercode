#!/usr/bin/env python

import commands
import re
import os

import sys
#sys.path.append('')

###########################################
###########################################

def remove( directory ):
    
    outSamples = re.split(r'\n',commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory))

    for sample in  outSamples:
          
        if  (re.search("Stream-A", sample )==None ):
            continue
        if  (re.search("GluGluToHToTauTau_M-120", directory )!=None or
             re.search("MuTauStream-A-05AugReReco", directory )!=None):
            continue

        rm = "rfrm -rf /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory+"/"+sample
        print rm
        os.system(rm)
 
###########################################
###########################################


remove( "DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola" )
remove( "GluGluToHToTauTau_M-100_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-105_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-110_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-115_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-120_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-125_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-130_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-135_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-140_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-145_7TeV-powheg-pythia6" )
remove( "GluGluToHToTauTau_M-160_7TeV-powheg-pythia6" )
remove( "QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6" )
remove( "TTJets_TuneZ2_7TeV-madgraph-tauola" )
remove( "TauPlusX" )
remove( "VBF_HToTauTau_M-100_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-105_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-110_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-125_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-130_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-135_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-140_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-145_7TeV-powheg-pythia6-tauola" )
remove( "VBF_HToTauTau_M-160_7TeV-powheg-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-100_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-110_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-115_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-125_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-130_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-135_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-140_7TeV-pythia6-tauola" )
remove( "WH_ZH_TTH_HToTauTau_M-145_7TeV-pythia6-tauola" )
remove( "WJetsToLNu_TuneZ2_7TeV-madgraph-tauola" )
remove( "WW_TuneZ2_7TeV_pythia6_tauola")
remove( "WZ_TuneZ2_7TeV_pythia6_tauola")
remove( "ZZ_TuneZ2_7TeV_pythia6_tauola")
