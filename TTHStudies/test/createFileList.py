#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')


def createList(sampleName, path, step, split):

     out = open('fileListToCopy_'+sampleName+'_'+str(split[0])+'-'+str(split[1])+'.txt','w')

     counter = 1
     for k in range(10):
          output = commands.getoutput("lcg-ls -o "+str(k*step)+" -c "+str(step)+" -b -D srmv2 "+path)
          outFiles = re.split(r'\n',output)
          for name in outFiles:
               if re.search(".root",name)!=None:
                    #foo = name.replace("/pnfs/lcg.cscs.ch/cms/trivcat","")
                    foo = name.replace("/pnfs/psi.ch/cms/trivcat", "")
                    if counter>=split[0] and counter<=split[1]:
                         out.write(foo+'\n')
                    counter = counter + 1
     out.close()


#################################################
#################################################


def useDataReplica(sampleName, dest, split):
    print "data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy_"+sampleName+'_'+str(split[0])+'-'+str(split[1])+".txt "+dest
    os.system("data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy_"+sampleName+'_'+str(split[0])+'-'+str(split[1])+".txt "+dest)


# T2
#createList("srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples", 500)
#createList("srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/TTZJets_8TeV-madgraph_v2",20)

#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [   1, 500])
#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 501,1000])
#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1001,2000])


# T3
#createList("WJets", "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball",500)
createList("DYJets10to50", "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph",500, [1, 5000])



#useDataReplica("/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball")
#useDataReplica("/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/TTZJets_8TeV-madgraph_v2")

#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 1, 500])
#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 501, 1000])
#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 1001, 2000])
