#!/usr/bin/env python

import commands
import re
import os

import sys

###########################################
###########################################

def publishSkim( stream, tasks ):

    f = open(stream+'_Skim.txt', 'w')
    
    for task in tasks:
        READ = 0
        SKIM = 0
        status = "crab -status -c "+task+"_skim/"
        print status
        os.system( status )
        getoutput = "crab -getoutput -c "+task+"_skim/"
        print getoutput
        os.system( getoutput )
        report = "crab -report -c "+task+"_skim/"
        print report
        output0 = commands.getoutput( report )
        reportLines = re.split(r'\n',output0)
        for line in reportLines:
            if  re.search("Total Events read:", line )!=None:
                words = re.split(r'\s', line)
                print words
                READ = float(words[3])
                
        publish = "crab -publish -c "+task+"_skim/"
        print publish
        #os.system( publish )
        output1 = commands.getoutput( publish )
        pusblishLines = re.split(r'\n',output1)
        for line in pusblishLines:
            if  re.search("total events:", line )!=None:
                words0 = re.split(r'/',line)
                words1 = re.split(r'\s', line)
                print words1, words0
                SKIM = float(words1[3])
                datasetpath = "/"+words0[1]+"/"+words0[2]+"/"+words0[3]
                #print words
                f.write('['+task+'_run]\n')
                f.write('#skim efficiency: '+str(SKIM)+'/'+str(READ)+' = '+str(SKIM/READ)+'\n')
                f.write('CMSSW.datasetpath='+datasetpath+'\n')
                f.write('CMSSW.total_number_of_events= -1\n')
                f.write('CMSSW.events_per_job = 5000\n')
                f.write('\n')
    f.close()

###########################################
###########################################


samples = ['DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola',
           'DYToMuMu_M-20_TuneZ2_7TeV-pythia6',
           'DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola',
           'DYToEE_M-20_TuneZ2_7TeV-pythia6',
           'GluGluToHToTauTau_M-125_7TeV-powheg-pythia6',
           'GluGluToHToTauTau_M-135_7TeV-powheg-pythia6',
           'VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola',
           'VBF_HToTauTau_M-125_7TeV-powheg-pythia6-tauola',
           'VBF_HToTauTau_M-135_7TeV-powheg-pythia6-tauola',
           'WJetsToLNu_TuneZ2_7TeV-madgraph-tauola',
           'WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola',
           'ZZTo2L2Nu_TuneZ2_7TeV_pythia6_tauola',
           'TT_TuneZ2_7TeV-pythia6-tauola']

tasksMu    = ['DYToTauTau-Mu-20-PUS3',
              'DYToMuMu-20-PUS3',
              #'VBFH115-Mu-powheg-PUS4',
              #'VBFH125-Mu-powheg-PUS4',
              #'VBFH135-Mu-powheg-PUS4',
              #'GGFH115-Mu-powheg-PUS1',
              #'GGFH135-Mu-powheg-PUS4',
              #'GGFH125-Mu-powheg-PUS4',
              'GGFH125-Mu-powheg-PUS4',
              'VBFH125-Mu-powheg-PUS4',
              'TT-Mu-pythia-PUS3',
              'DYJets-Mu-50-madgraph-PUS4',
              'WW-Mu-pythia-PUS4',
              #'WZ-Mu-pythia-PUS4',
              'ZZ-Mu-pythia-PUS4',
              'WJets-Mu-madgraph-PUS4']

tasksElec  = ['DYToEE-20-PUS3',
              'DYToTauTau-20-PUS3',
              'TT-pythia-PUS3',
              #'GGFH115-powheg-PUS1',
              'GGFH125-powheg-PUS4',
              #'GGFH135-powheg-PUS4',
              #'VBFH115-powheg-PUS4',
              'VBFH125-powheg-PUS4'
              #'VBFH135-powheg-PUS4',
              'WJets-madgraph-PUS4',
              'DYJets-50-madgraph-PUS4',
              #'WW-pythia-PUS4',
              #'WZ-pythia-PUS4',
              'ZZ-pythia-PUS4']

#publishSkim( "ElecTauStream",  tasksElec )
publishSkim( "ElecTauStream_patch1",  ['VBFH125-powheg-PUS4']   )
