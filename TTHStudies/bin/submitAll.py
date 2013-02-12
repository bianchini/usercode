#!/usr/bin/env python


import commands
import re
import os

import sys



def processAll():

    pathI   = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"
    pathO   = "/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"
    fileN   = "DiJetPt_"
    command ="srmcp -2 "

    
    samples = [  #['WW_TuneZ2star_8TeV_pythia6_tauola', 'WW'],
                 #['WZ_TuneZ2star_8TeV_pythia6_tauola', 'WZ'],
                 #['ZZ_TuneZ2star_8TeV_pythia6_tauola', 'ZZ'],
                 #['DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball', 'DYJets'],
                 #['T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola', 'TtW'],
                 ['T_s-channel_TuneZ2star_8TeV-powheg-tauola', 'Ts'],
                 #['T_t-channel_TuneZ2star_8TeV-powheg-tauola', 'Tt'],
                 #['Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola', 'Tbars'],
                 #['Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola', 'Tbart'],
                 #['Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola', 'TbartW'],
                 #['TTJets_HadronicMGDecays_8TeV-madgraph-part', 'TTJetsFullHad'],
                 #['TTJets_FullLeptMGDecays_8TeV-madgraph-part', 'TTJetsFullLept'],
                 #['TTJets_SemiLeptMGDecays_8TeV-madgraph-part', 'TTJetsSemiLept'],
                 #['WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph',     'WJets100'],
                 #['WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph', 'WJets70100'],
                 #['ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp',  'ZH25'],
                 #['WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp', 'WH25'],
                 #['SingleElectronRun2012AJul13EdmV42b',  'SingleElectron_1'],
                 #['SingleElectronRun2012AAug06EdmV42',   'SingleElectron_2'],
                 #['SingleElectronRun2012BJul13EdmV42',   'SingleElectron_3'],
                 #['SingleElectronRun2012CAug24RerecoEdmV42', 'SingleElectron_4'],
                 #['SingleElectronRun2012CPromptv2EdmV42', 'SingleElectron_5'],
                 #['SingleElectronRun2012CPromptV2TopUpEdmV42', 'SingleElectron_6'],
                 #['SingleMuRun2012AJul13EdmV42', 'SingleMu_1'],
                 #['SingleMuRun2012AAug06EdmV42', 'SingleMu_2'],
                 #['SingleMuRun2012BJul13EdmV42', 'SingleMu_3'],
                 #['SingleMuRun2012CAug24RerecoEdmV42', 'SingleMu_4'],
                 #['SingleMuRun2012CPromptv2EdmV42','SingleMu_5 '],
                 #['SingleMuRun2012CPromptV2TopUpEdmV42','SingleMu_6']
                 ]


    for sample in samples:
        print command+pathI+fileN+sample[0]+".root" +" file:///"+pathO+fileN+sample[0]+".root"
        os.system(command+pathI+fileN+sample[0]+".root" +" file:///"+pathO+fileN+sample[0]+".root")
        print "Step3 step3.py "+sample[0]+" "+sample[1]
        os.system("Step3 step3.py "+sample[0]+" "+sample[1])
        print "Skim  skim.py "+sample[0]+" "+sample[1]
        os.system("Skim  skim.py "+sample[0]+" "+sample[1])
        #print "rm "+pathO+"*"+sample[1]+"*.root"
        print "rm "+pathO+"*"+sample[0]+"*.root"
        #os.system("rm "+pathO+"*"+sample[1]+"*.root")
        os.system("rm "+pathO+"*"+sample[0]+"*.root")

###########################################
###########################################




processAll()
