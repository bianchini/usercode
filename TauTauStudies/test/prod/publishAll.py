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

tasksMu    = [
    #'DYToTauTau-Mu-20-PUS3',
    #'DYToMuMu-20-PUS3',
    #'VBFH115-Mu-powheg-PUS4',
    #'VBFH125-Mu-powheg-PUS4',
    #'VBFH135-Mu-powheg-PUS4',
    #'GGFH115-Mu-powheg-PUS1',
    #'GGFH135-Mu-powheg-PUS4',
    #'GGFH125-Mu-powheg-PUS4',
    #'GGFH125-Mu-powheg-PUS4',
    #'VBFH125-Mu-powheg-PUS4',
    #'TT-Mu-pythia-PUS3',
    #'DYJets-Mu-50-madgraph-PUS4',
    #'WW-Mu-pythia-PUS4',
    #'WZ-Mu-pythia-PUS4',
    #'ZZ-Mu-pythia-PUS4',
    #'WJets-Mu-madgraph-PUS4'
    #'T-Mu-t-PUS1_skim',
    #'Tbar-Mu-t-PUS1_skim',
    #'WZIncl-Mu-pythia-PUS4_skim',
    #'TTJets-Mu-madgraph-PUS4_skim',
    #'GGFH105-Mu-powheg-PUS1_skim',
    #'GGFH110-Mu-powheg-PUS1_skim',
    #'GGFH115-Mu-powheg-PUS1_skim',
    #'GGFH120-Mu-powheg-PUS1_skim',
    #'GGFH125-Mu-powheg-PUS4_skim-v2'
    #'GGFH130-Mu-powheg-PUS4_skim-v2',
    #'GGFH135-Mu-powheg-PUS4_skim-v2',
    #'GGFH140-Mu-powheg-PUS4_skim-v2',
    #'VBFH105-Mu-powheg-PUS4_skim-v2',
    #'VBFH110-Mu-powheg-PUS4_skim-v2',
    #'VBFH115-Mu-powheg-PUS4_skim',
    #'VBFH120-Mu-powheg-PUS4_skim',
    #'VBFH120-Mu-powheg-PUS4_skim',
    #'VBFH125-Mu-powheg-PUS4_skim-v2',
    #'VBFH130-Mu-powheg-PUS4_skim-v2',
    #'VBFH135-Mu-powheg-PUS4_skim-v2',
    #'VBFH140-Mu-powheg-PUS4_skim-v2',
    'ZZ-Mu-pythia-PUS4_skim',
    'T-Mu-s-PUS4_skim',
    'T-Mu-tW-PUS4_skim',
    'Tbar-Mu-tW-PUS4_skim',
    'WW-Mu-pythia-PUS4_skim',
    ]

tasksElec  = [
    #'DYToEE-20-PUS3',
    #'DYToTauTau-20-PUS3',
    #'TT-pythia-PUS3',
    #'GGFH115-powheg-PUS1',
    #'GGFH125-powheg-PUS4',
    #'GGFH135-powheg-PUS4',
    #'VBFH115-powheg-PUS4',
    #'VBFH125-powheg-PUS4'
    #'VBFH135-powheg-PUS4',
    #'WJets-madgraph-PUS4',
    #'DYJets-50-madgraph-PUS4',
    #'WW-pythia-PUS4',
    #'WZ-pythia-PUS4',
    #'ZZ-pythia-PUS4',
    #'T-t-PUS1_skim',
    #'Tbar-t-PUS1_skim',
    #'WZIncl-pythia-PUS4_skim',
    #'TTJets-madgraph-PUS4_skim',
    #'GGFH105-powheg-PUS1_skim',
    #'GGFH110-powheg-PUS1_skim',
    #'GGFH115-powheg-PUS1_skim',
    #'GGFH120-powheg-PUS1_skim',
    #'GGFH125-powheg-PUS4_skim-v4',
    #'GGFH130-powheg-PUS4_skim-v2',
    #'GGFH135-powheg-PUS4_skim-v4',
    #'GGFH140-powheg-PUS4_skim-v2',
    #'GGFH145-powheg-PUS4_skim',
    #'GGFH160-powheg-PUS4_skim',
    #'VBFH105-powheg-PUS4_skim-v2',
    #'VBFH110-powheg-PUS4_skim-v2',
    #'VBFH115-powheg-PUS4_skim',
    #'VBFH120-powheg-PUS4_skim-v2',
    #'VBFH125-powheg-PUS4_skim-v4',
    #'VBFH130-powheg-PUS4_skim-v2',
    #'VBFH135-powheg-PUS4_skim-v4',
    #'VBFH140-powheg-PUS4_skim-v2',
    #'VBFH145-powheg-PUS4_skim',
    #'VBFH160-powheg-PUS4_skim',
    #'Run2011-05AugReReco_skim',
    #'Run2011-PromptReco-v6_skim/'
    'ZZ-pythia-PUS4_skim',
    'Tbar-tW-PUS4_skim',
    'T-tW-PUS4_skim',
    'Tbar-s-PUS4_skim',
    'WW-pythia-PUS4_skim',
    ]

tasksElecMu =[
    #'DYJets-ElecMu-50-madgraph-PUS4_skim',
    #'T-ElecMu-s-PUS4_skim',
    #'T-ElecMu-t-PUS4_skim',
    #'T-ElecMu-tW-PUS4_skim',
    #'Tbar-ElecMu-s-PUS4_skim',
    #'Tbar-ElecMu-t-PUS4_skim',
    #'Tbar-ElecMu-tW-PUS4_skim',
    #'TTJets-ElecMu-madgraph-PUS4_skim',
    #'WJets-ElecMu-madgraph-PUS4_skim',
    #'WW-ElecMu-madgraph-PUS4_skim',
    #'WW-ElecMu-pythia-PUS4_skim',
    #'WZ-ElecMu-pythia-PUS4_skim',
    #'ZZ-ElecMu-pythia-PUS4_skim',
    #'GGFH105-ElecMu-powheg-PUS1_skim',
    #'GGFH110-ElecMu-powheg-PUS1_skim',
    #'GGFH115-ElecMu-powheg-PUS1_skim',
    #'GGFH120-ElecMu-powheg-PUS1_skim',
    #'GGFH125-ElecMu-powheg-PUS4_skim',
    #'GGFH130-ElecMu-powheg-PUS4_skim',
    #'GGFH135-ElecMu-powheg-PUS4_skim',
    #'GGFH140-ElecMu-powheg-PUS4_skim',
    #'VBFH105-ElecMu-powheg-PUS4_skim',
    #'VBFH110-ElecMu-powheg-PUS4_skim',
    #'VBFH115-ElecMu-powheg-PUS4_skim',
    #'VBFH120-ElecMu-powheg-PUS4_skim',
    #'VBFH125-ElecMu-powheg-PUS4_skim',
    #'VBFH130-ElecMu-powheg-PUS4_skim',
    #'VBFH135-ElecMu-powheg-PUS4_skim',
    #'VBFH140-ElecMu-powheg-PUS4_skim'    
    ]

tasksMuSUSY    = [
    'SUSYBB90-MuTau-powheg-PUS6_skim',
    'SUSYBB100-MuTau-powheg-PUS6_skim',
    'SUSYBB120-MuTau-powheg-PUS6_skim',
    'SUSYBB130-MuTau-powheg-PUS6_skim',
    'SUSYBB140-MuTau-powheg-PUS6_skim',
    'SUSYBB160-MuTau-powheg-PUS6_skim',
    'SUSYBB180-MuTau-powheg-PUS6_skim',
    'SUSYBB200-MuTau-powheg-PUS6_skim',
    'SUSYBB250-MuTau-powheg-PUS6_skim',
    'SUSYBB300-MuTau-powheg-PUS6_skim',
    'SUSYBB350-MuTau-powheg-PUS6_skim',
    'SUSYBB400-MuTau-powheg-PUS6_skim',
    #'SUSYBB450-MuTau-powheg-PUS6_skim',
    'SUSYBB500-MuTau-powheg-PUS6_skim',
    'SUSYBB600-MuTau-powheg-PUS6_skim',
    'SUSYBB700-MuTau-powheg-PUS6_skim',
    'SUSYBB800-MuTau-powheg-PUS6_skim',
    'SUSYBB900-MuTau-powheg-PUS6_skim',
    'SUSYGG90-MuTau-powheg-PUS6_skim',
    'SUSYGG100-MuTau-powheg-PUS6_skim',
    'SUSYGG120-MuTau-powheg-PUS6_skim',
    'SUSYGG130-MuTau-powheg-PUS6_skim',
    'SUSYGG140-MuTau-powheg-PUS6_skim',
    'SUSYGG160-MuTau-powheg-PUS6_skim',
    'SUSYGG180-MuTau-powheg-PUS6_skim',
    'SUSYGG200-MuTau-powheg-PUS6_skim',
    'SUSYGG250-MuTau-powheg-PUS6_skim',
    'SUSYGG300-MuTau-powheg-PUS6_skim',
    'SUSYGG350-MuTau-powheg-PUS6_skim',
    'SUSYGG400-MuTau-powheg-PUS6_skim',
    'SUSYGG450-MuTau-powheg-PUS6_skim',
    'SUSYGG500-MuTau-powheg-PUS6_skim',
    'SUSYGG600-MuTau-powheg-PUS6_skim',
    'SUSYGG700-MuTau-powheg-PUS6_skim',
    'SUSYGG800-MuTau-powheg-PUS6_skim',
    'SUSYGG900-MuTau-powheg-PUS6_skim',
    ]

tasksElecSUSY    = [
    #'SUSYBB90-ElecTau-powheg-PUS4_skim',
    #'SUSYBB100-ElecTau-powheg-PUS4_skim',
    #'SUSYBB120-ElecTau-powheg-PUS4_skim',
    #'SUSYBB130-ElecTau-powheg-PUS4_skim',
    #'SUSYBB140-ElecTau-powheg-PUS4_skim',
    #'SUSYBB160-ElecTau-powheg-PUS4_skim',
    #'SUSYBB180-ElecTau-powheg-PUS4_skim',
    #'SUSYBB200-ElecTau-powheg-PUS4_skim',
    #'SUSYBB250-ElecTau-powheg-PUS4_skim',
    #'SUSYBB300-ElecTau-powheg-PUS4_skim',
    #'SUSYBB350-ElecTau-powheg-PUS4_skim',
    #'SUSYBB400-ElecTau-powheg-PUS4_skim',
    #'SUSYBB450-ElecTau-powheg-PUS4_skim',
    #'SUSYBB500-ElecTau-powheg-PUS4_skim',
    #'SUSYBB600-ElecTau-powheg-PUS4_skim',
    #'SUSYBB700-ElecTau-powheg-PUS4_skim',
    #'SUSYBB800-ElecTau-powheg-PUS4_skim',
    #'SUSYBB900-ElecTau-powheg-PUS4_skim',
    #'SUSYGG90-ElecTau-powheg-PUS6_skim',
    #'SUSYGG100-ElecTau-powheg-PUS6_skim',
    #'SUSYGG120-ElecTau-powheg-PUS6_skim',
    #'SUSYGG130-ElecTau-powheg-PUS6_skim',
    #'SUSYGG140-ElecTau-powheg-PUS6_skim',
    #'SUSYGG160-ElecTau-powheg-PUS6_skim',
    #'SUSYGG180-ElecTau-powheg-PUS6_skim',  
    #'SUSYGG200-ElecTau-powheg-PUS6_skim',
    #'SUSYGG250-ElecTau-powheg-PUS6_skim',
    #'SUSYGG300-ElecTau-powheg-PUS6_skim',
    'SUSYGG350-ElecTau-powheg-PUS6_skim',
    #'SUSYGG400-ElecTau-powheg-PUS6_skim',
    #'SUSYGG450-ElecTau-powheg-PUS6_skim',
    #'SUSYGG500-ElecTau-powheg-PUS6_skim',
    #'SUSYGG600-ElecTau-powheg-PUS6_skim',
    #'SUSYGG700-ElecTau-powheg-PUS6_skim',
    #'SUSYGG800-ElecTau-powheg-PUS6_skim',
    #'SUSYGG900-ElecTau-powheg-PUS6_skim',
    ]

tasksMuTauA = [
    #'Run2011A-MuTau-05AugReReco_skim',
    #'Run2011A-MuTau-May10ReReco_skim',
    #'Run2011A-MuTau-PromptReco-v4_skim',
    #'Run2011A-MuTau-PromptReco-v6-p1_skim',
    #'Run2011A-MuTau-PromptReco-v6-p2_skim',
    #'DYJets-MuTau-50-madgraph-PUS6_skim',
    #'GGFH100-MuTau-powheg-PUS6_skim',
    #'GGFH105-MuTau-powheg-PUS6_skim',
    #'GGFH110-MuTau-powheg-PUS6_skim',
    #'GGFH115-MuTau-powheg-PUS6_skim',
    #'GGFH120-MuTau-powheg-PUS6_skim',
    #'GGFH125-MuTau-powheg-PUS6_skim',
    #'GGFH130-MuTau-powheg-PUS6_skim',
    #'GGFH135-MuTau-powheg-PUS6_skim',
    #'GGFH140-MuTau-powheg-PUS6_skim',
    #'GGFH145-MuTau-powheg-PUS6_skim',
    #'GGFH160-MuTau-powheg-PUS6_skim',
    #'QCDmu-MuTau-pythia-20-15-PUS6_skim',
    #'TTJets-MuTau-madgraph-PUS6_skim',
    #'VBFH100-MuTau-powheg-PUS6_skim',
    #'VBFH105-MuTau-powheg-PUS6_skim',
    #'VBFH110-MuTau-powheg-PUS6_skim',
    #'VBFH115-MuTau-powheg-PUS6_skim',
    #'VBFH120-MuTau-powheg-PUS6_skim',
    #'VBFH125-MuTau-powheg-PUS6_skim',
    #'VBFH130-MuTau-powheg-PUS6_skim',
    #'VBFH135-MuTau-powheg-PUS6_skim',
    #'VBFH140-MuTau-powheg-PUS6_skim',
    #'VBFH145-MuTau-powheg-PUS6_skim',
    #'VBFH160-MuTau-powheg-PUS6_skim',
    #'VH100-MuTau-pythia-PUS6_skim',
    #'VH105-MuTau-pythia-PUS6_skim',
    #'VH110-MuTau-pythia-PUS6_skim',
    #'VH115-MuTau-pythia-PUS6_skim',
    #'VH120-MuTau-pythia-PUS6_skim',
    #'VH125-MuTau-pythia-PUS6_skim',
    #'VH130-MuTau-pythia-PUS6_skim',
    #'VH135-MuTau-pythia-PUS6_skim',
    #'VH140-MuTau-pythia-PUS6_skim',
    #'VH145-MuTau-pythia-PUS6_skim',
    #'VH160-MuTau-pythia-PUS6_skim',
    #'WJets-MuTau-madgraph-PUS6_skim',
    'W3Jets-MuTau-madgraph-PUS6_skim'
    #'WW-MuTau-pythia-PUS6_skim',
    #'WZ-MuTau-pythia-PUS6_skim',
    #'ZZ-MuTau-pythia-PUS6_skim'
    ]

tasksMuTauB = [
    #'Run2011B-MuTau-PromptReco-v1-p1_skim',
    #'Run2011B-MuTau-PromptReco-v1-p2_skim',
    #'Run2011B-MuTau-PromptReco-v1-p3_skim'
    'Run2011B-MuTau-PromptReco-v1-p4_skim'
    ]

tasksElecTauA = [
    #'Run2011A-ElecTau-05AugReReco_skim',
    #'Run2011A-ElecTau-May10ReReco_skim',
    #'Run2011A-ElecTau-PromptReco-v4_skim',
    #'Run2011A-ElecTau-PromptReco-v6-part1_skim',
    #'Run2011A-ElecTau-PromptReco-v6-part2_skim',
    #'Run2011A-ElecTau-03OctReReco_skim',
    #'Run2011B-ElecTau-PromptReco-v1_skim',
    #'DYJets-ElecTau-50-madgraph-PUS6-v2_skim',
    #'GGFH100-ElecTau-powheg-PUS6_skim',
    #'GGFH105-ElecTau-powheg-PUS6_skim',
    #'GGFH110-ElecTau-powheg-PUS6_skim',
    #'GGFH115-ElecTau-powheg-PUS6_skim',
    #'GGFH120-ElecTau-powheg-PUS6_skim',
    #'GGFH125-ElecTau-powheg-PUS6_skim',
    #'GGFH130-ElecTau-powheg-PUS6_skim',
    #'GGFH135-ElecTau-powheg-PUS6_skim',
    #'GGFH140-ElecTau-powheg-PUS6_skim',
    #'GGFH145-ElecTau-powheg-PUS6_skim',
    #'GGFH160-ElecTau-powheg-PUS6_skim',
    #'TTJets-ElecTau-madgraph-PUS6_skim',
    #'VBFH100-ElecTau-powheg-PUS6_skim',
    #'VBFH105-ElecTau-powheg-PUS6_skim',
    #'VBFH110-ElecTau-powheg-PUS6_skim',
    #'VBFH115-ElecTau-powheg-PUS6_skim',
    #'VBFH120-ElecTau-powheg-PUS6_skim',
    #'VBFH125-ElecTau-powheg-PUS6_skim',
    #'VBFH130-ElecTau-powheg-PUS6_skim',
    #'VBFH135-ElecTau-powheg-PUS6_skim',
    #'VBFH140-ElecTau-powheg-PUS6_skim',
    #'VBFH145-ElecTau-powheg-PUS6_skim',
    #'VBFH160-ElecTau-powheg-PUS6_skim',
    #'VH100-ElecTau-pythia-PUS6_skim',
    #'VH105-ElecTau-pythia-PUS6_skim',
    #'VH110-ElecTau-pythia-PUS6_skim',
    #'VH115-ElecTau-pythia-PUS6_skim',
    #'VH120-ElecTau-pythia-PUS6_skim',
    #'VH125-ElecTau-pythia-PUS6_skim',
    #'VH130-ElecTau-pythia-PUS6_skim',
    #'VH135-ElecTau-pythia-PUS6_skim',
    #'VH140-ElecTau-pythia-PUS6_skim',
    #'VH145-ElecTau-pythia-PUS6_skim',
    #'VH160-ElecTau-pythia-PUS6_skim',
    #'WJets-ElecTau-madgraph-PUS6_skim',
    'W3Jets-ElecTau-madgraph-PUS6_skim',
    #'WW-ElecTau-pythia-PUS6_skim',
    #'WZ-ElecTau-pythia-PUS6_skim',
    #'ZZ-ElecTau-pythia-PUS6_skim'
    ]

tasksElecSync  = [
    'GGFH120-ElecTau-powheg-PUS4_skim',
    'Run2011A-ElecTau-05AugReReco_skim'
    ]

tasksElecEmbedding  = [
    #'Run2011A-ElecTau-May10ReReco-Embedded_skim',
    #'Run2011A-ElecTau-PromptReco-v4-Embedded_skim',
    #'Run2011A-ElecTau-05AugReReco-Embedded_skim',
    #'Run2011A-ElecTau-PromptReco-v6-Embedded_skim'
    'Run2011B-ElecTau-PromptReco-v1-Embedded_skim'
    ]

tasksMuEmbedding  = [
    #'Run2011A-MuTau-May10ReReco-Embedded_skim',
    #'Run2011A-MuTau-PromptReco-v4-Embedded_skim',
    #'Run2011A-MuTau-PromptReco-v6-Embedded_skim'
    'Run2011A-MuTau-05AugReReco-Embedded_skim'
    #'Run2011B-MuTau-PromptReco-v1-Embedded_skim'
    ]
#publishSkim( "MuTauStream_patch2",  tasksMu )
#publishSkim( "ElecTauStream_patchVV",  tasksElec   )
#publishSkim( "ElecMuStream",  tasksElecMu   )
#publishSkim( "MuTauStream_patchSUSY",  tasksMuSUSY   )
#publishSkim( "MuTauStream_A",    tasksMuTauA   )
#publishSkim( "ElecTauStream_A",   tasksElecSync  )

#publishSkim( "MuTauStream_13Oct2011",  tasksMuTauA   )
#publishSkim( "ElecTauStream_08Nov2011_patch1",  tasksElecTauA   )
#publishSkim( "ElecTauStream_08Nov2011_Embedded_patch1",  tasksElecEmbedding   )
#publishSkim( "MuTauStream_16Nov2011_patch5",  tasksMuTauA   )
publishSkim( "ElecTauStream_01Dec2011_patch2",  tasksElecTauA   )
#publishSkim( "ElecTauStream_01Dec2011_patchSUSY_patch1",  tasksElecSUSY   )
#publishSkim( "MuTauStream_16Nov2011_patchSUSY",  tasksMuSUSY   )
#publishSkim( "MuTauStream_13Oct2011_Embedded_patch1",  tasksMuEmbedding   )
#publishSkim( "MuTauStream_13Oct2011_patchB",  tasksMuTauB   )
#publishSkim( "ElecTauStream_08Nov2011_Embedded", tasksElecEmbedding)
