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
        #status = "crab -status -c "+task+"_skim/"
        status = "crab -status -c "+task
        print status
        os.system( status )
        #getoutput = "crab -getoutput -c "+task+"_skim/"
        getoutput = "crab -getoutput -c "+task
        print getoutput
        os.system( getoutput )
        #report = "crab -report -c "+task+"_skim/"
        report = "crab -report -c "+task
        print report
        output0 = commands.getoutput( report )
        reportLines = re.split(r'\n',output0)
        for line in reportLines:
            if  re.search("Total Events read:", line )!=None:
                words = re.split(r'\s', line)
                print words
                READ = float(words[3])
                
        #publish = "crab -publish -c "+task+"_skim/"
        publish = "crab -publish -c "+task
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




tasksMuTauSummer12 = [
    #'Run2012A-MuTau-PromptReco-v1-p1-v2_skim'
    #'Run2012A-MuTau-PromptReco-v1-p2_skim',
    #'DYJets-MuTau-madgraph-Tarball-v2_skim'
    #'DYJets-MuTau-madgraph-50_skim',
    'TTJets-MuTau-madgraph-50_skim',
    #'WMuNu-MuTau-pythia-50_skim',
    #'WTauNu-MuTau-pythia-50_skim',
    #'WENu-MuTau-pythia-50_skim',
    #'WZ-MuTau-50_skim',
    #'WW-MuTau-50_skim',
    #'ZZ-MuTau-50_skim',
    #'VBFH120-MuTau-pythia-50_skim'
    ]

tasksElecTauSummer12 = [
    #'Run2012A-ElecTau-PromptReco-v1-p1-v2_skim'
    #'Run2012A-ElecTau-PromptReco-v1-p2_skim',
    #'DYJets-ElecTau-madgraph-Tarball-v2_skim'
    #'DYJets-ElecTau-madgraph-50_skim',
    'TTJets-ElecTau-madgraph-50_skim',
    #'WENu-ElecTau-pythia-50_skim',
    #'WTauNu-ElecTau-pythia-50_skim',
    #'WMuNu-ElecTau-pythia-50_skim',
    #'WW-ElecTau-50_skim',
    #'WZ-ElecTau-50_skim',
    #'ZZ-ElecTau-50_skim',
    #'VBFH120-ElecTau-pythia-50_skim',
    ]




#publishSkim( "MuTauStream_30Mar2012_patch50X_v4",  tasksMuTauSummer12)
publishSkim( "ElecTauStream_30Mar2012_patch50X_v4",tasksElecTauSummer12)

