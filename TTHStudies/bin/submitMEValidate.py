#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms


def submitMEValidate(script,
                     vegasPoints,
                     mode,
                     norm,
                     useME, useJac, useMET, useTF, usePDF,
                     doParton, doSmear, doMassScan, doPermutations,
                     met,
                     masses,
                     evLow,evHigh):

    print "Overload meValidator.py..."
    os.system('cp ../python/meValidator.py ./')

    from meValidator import process
    
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    process.fwliteInput.outFileName      = cms.string('./root/MEValidator_'+script+'.root')
    process.fwliteInput.vegasPoints      = cms.int32(vegasPoints)
    process.fwliteInput.mode             = cms.untracked.int32(mode)
    process.fwliteInput.norm             = cms.untracked.int32(norm)
    process.fwliteInput.useME            = cms.int32(useME)
    process.fwliteInput.useJac           = cms.int32(useJac)
    process.fwliteInput.useMET           = cms.int32(useMET)
    process.fwliteInput.useTF            = cms.int32(useTF)
    process.fwliteInput.usePDF           = cms.int32(usePDF)
    process.fwliteInput.doParton         = cms.int32(doParton)
    process.fwliteInput.doSmear          = cms.int32(doSmear)   
    process.fwliteInput.doMassScan       = cms.int32(doMassScan)
    process.fwliteInput.doPermutations   = cms.int32(doPermutations)

    process.fwliteInput.masses           = masses
    process.fwliteInput.evLimits         = cms.vint32(evLow,evHigh)

    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/bin/\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('MEValidator ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+jobName+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
    #os.system('rm testME_tmp.py')
  
###########################################
###########################################

masses = cms.vdouble(60,  65,  70,  75, 80 , 85,  90, 95, 100, 105, 110, 
                     115, 120, 125, 130, 135, 140, 145, 150, 155, 
                     160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250)

'''
counter = 0
for i in range(20):
    counter = counter + 1
    submitMEValidate('SL2wj_gen_acc_p'+str(counter),     2000, 0, 2, 1,1,1,1,1,  1,0,1,1,  125,  masses,  i*25+1, (i+1)*25 )
    submitMEValidate('SL2wj_rec_acc_p'+str(counter),     2000, 0, 2, 1,1,1,1,1,  0,1,1,1,  125,  masses,  i*25+1, (i+1)*25 )
    submitMEValidate('SL2wj_gen_unnorm_p'+str(counter),  2000, 0, 0, 1,1,1,1,1,  1,0,1,1,  125,  masses,  i*25+1, (i+1)*25 )
    submitMEValidate('SL2wj_rec_unnorm_p'+str(counter),  2000, 0, 0, 1,1,1,1,1,  0,1,1,1,  125,  masses,  i*25+1, (i+1)*25 )

counter = 0
for i in range(20):
    counter = counter + 1
    submitMEValidate('SL1wj_gen_acc_p'+str(counter),     4000, 1, 2, 1,1,1,1,1,  1,0,1,1,  125,  masses,  i*10+1, (i+1)*10 )
    submitMEValidate('SL1wj_rec_acc_p'+str(counter),     4000, 1, 2, 1,1,1,1,1,  0,1,1,1,  125,  masses,  i*10+1, (i+1)*10 )
'''

counter = 0
for i in range(20):
    counter = counter + 1
    submitMEValidate('SLNoBHad_gen_acc_p'+str(counter),   10000, 2, 2, 1,1,1,1,1,  1,0,1,1,  125,  masses,  i*10+1, (i+1)*10 )
    submitMEValidate('SLNoBHad_rec_acc_p'+str(counter),   10000, 2, 2, 1,1,1,1,1,  0,1,1,1,  125,  masses,  i*10+1, (i+1)*10 )

counter = 0
for i in range(20):
    counter = counter + 1
    submitMEValidate('SLNoBLep_gen_acc_p'+str(counter),   10000, 3, 2, 1,1,1,1,1,  1,0,1,1,  125,  masses,  i*10+1, (i+1)*10 )
    submitMEValidate('SLNoBLep_rec_acc_p'+str(counter),   10000, 3, 2, 1,1,1,1,1,  0,1,1,1,  125,  masses,  i*10+1, (i+1)*10 )
                
    
