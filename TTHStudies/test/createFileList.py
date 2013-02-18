#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')


def createList(path, step):

     out = open('fileListToCopy.txt','w')

     for k in range(10):
          output = commands.getoutput("lcg-ls -o "+str(k*step)+" -c "+str(step)+" -b -D srmv2 "+path)
          outFiles = re.split(r'\n',output)
          for name in outFiles:
              if re.search(".root",name)!=None:
                  out.write(name+'\n')
     
     out.close()


#################################################
#################################################


def useDataReplica(dest):
    print "data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy.txt "+dest
    os.system("data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy.txt "+dest)








#createList("srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples", 500)
useDataReplica("/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball")
