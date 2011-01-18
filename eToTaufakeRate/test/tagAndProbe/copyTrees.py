#!/usr/bin/env python

import commands
import re
import os

import sys
#sys.path.append('/home/llr/cms/lbianchini/CMSSW_3_8_2/src/PFAnalyses/Z/test/')

###########################################
###########################################

def copyAndAdd( directory ):
    
    outSamples = re.split(r'\n',commands.getoutput("lcg-ls srm://polgrid4.in2p3.fr:8446/srm/managerv2?SFN=/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory))

    for sample in  outSamples:
        samplesplit =  re.split(r'/',sample)
        #print samplesplit
        output = commands.getoutput("lcg-ls srm://polgrid4.in2p3.fr:8446/srm/managerv2?SFN=/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory+"/"+samplesplit[len(samplesplit)-1])

        ### put here the namse of the samples to be imported:
        if  re.search("PYTHIA", samplesplit[len(samplesplit)-1])==None:
            continue
        
        outFiles = re.split(r'\n',output)
        for name in outFiles:
            splitname = re.split(r'/',name)
            print "copying file "+splitname[ len(splitname) -1 ]
            cp = "rfcp /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory+"/"+samplesplit[len(samplesplit)-1]+"/"+splitname[ len(splitname) -1 ]+" ."
            print cp
            os.system(cp)
        hadd = "hadd -f  testNewWriteFromPAT_"+samplesplit[len(samplesplit)-1]+".root testNewWriteFromPAT_*_*_*.root"
        print hadd
        os.system(hadd)
        remove = "rm -r testNewWriteFromPAT_*_*_*.root"
        os.system(remove)
 
###########################################
###########################################


copyAndAdd( "treesForTauEffVsFakeRate/" )
