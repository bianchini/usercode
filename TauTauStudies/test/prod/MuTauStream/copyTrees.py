#!/usr/bin/env python

import commands
import re
import os

import sys
#sys.path.append('')

###########################################
###########################################

def copyAndAdd( directory ):
    
    outSamples = re.split(r'\n',commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory))

    for sample in  outSamples:
        samplesplit =  re.split(r'/',sample)
        output = commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory+"/"+samplesplit[len(samplesplit)-1])
        
        if (re.search("05Aug", samplesplit[len(samplesplit)-1])==None ):
            continue

        outFiles = re.split(r'\n',output)
        for name in outFiles:
            splitname = re.split(r'/',name)
            print "copying file "+splitname[ len(splitname) -1 ]
            cp = "rfcp /dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/"+directory+"/"+samplesplit[len(samplesplit)-1]+"/"+splitname[ len(splitname) -1 ]+" ."
            print cp
            os.system(cp)
        hadd = "hadd -f  treeMuTauStream_"+samplesplit[len(samplesplit)-1]+".root treeMuTauStream_*_*_*.root"
        print hadd
        os.system(hadd)
        #os.system("cp  treeMuTauStream_"+samplesplit[len(samplesplit)-1]+".root ../")
        remove = "rm -r treeMuTauStream_*_*_*.root"
        os.system(remove)
 
###########################################
###########################################


copyAndAdd( "MuTauStreamFall11_06Dec2011_LooseIso" )
