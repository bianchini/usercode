#!/usr/bin/env python

import commands
import re
import os

import sys
sys.path.append('./')

###########################################
###########################################

def treeSkim( sample, xsection, runInSeries=False):

    os.system( 'mkdir batch/' )
    
    if runInSeries:
         print "Running in series via the command ..."
         runInSeries   = 'treeSkimmerMuTau_MVA '+sample+' '+str(xsection)
         print runInSeries
         os.system(runInSeries)
         mv     = 'mv nTuple'+sample+'_Open_MuTauStream_*.root batch/' 
         hadd   = 'hadd -f batch/nTuple'+sample+'_Open_MuTauStream.root batch/nTuple'+sample+'_Open_MuTauStream_*.root'
         remove = 'rm batch/nTuple'+sample+'_Open_MuTauStream_*.root'
         print 'Now doing hadd...'
         os.system(mv)
         os.system(hadd)
         os.system(remove)

    else:
        print "Creating the shell file for the batch..."
        f = open('batch/job'+'_'+sample+'.sh','w')    
        f.write('#!/bin/sh\n\n')
        f.write('export WORKINGDIR="/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7/src/Bianchi/Limits/bin/"\n')
        f.write('source /opt/exp_soft/cms/cmsset_default.sh\n')
        f.write('cd $WORKINGDIR\n')
        f.write('eval `scram runtime -sh`\n')
        f.write('source /opt/exp_soft/cms/cms_ui_env_crab.sh\n')
        f.write('treeSkimmerMuTau_MVA '+sample+' '+str(xsection)+' \n')
        f.write('mv nTuple'+sample+'_Open_MuTauStream_*.root batch/\n')
        f.write('hadd -f batch/nTuple'+sample+'_Open_MuTauStream.root batch/nTuple'+sample+'_Open_MuTauStream_*.root\n')
        f.write('rm batch/nTuple'+sample+'_Open_MuTauStream_*.root\n')
        f.close()
        #os.system('source /opt/exp_soft/cms/t3/t3setup')
        print "Submitting the job to the cms queue via the command ..."
        submitToQueue = '/opt/exp_soft/cms/t3/t3submit -q cms -k eo -N '+sample+' /home/llr/cms/lbianchini/CMSSW_4_2_8_patch7/src/Bianchi/Limits/bin/batch/job'+'_'+sample+'.sh > batch/'+sample+'.txt'
        print submitToQueue
        os.system(submitToQueue)
  
###########################################
###########################################

#treeSkim("GGFH120-MuTau-powheg-PUS6_run",      7.10e-02*16.63 * 1.0         * 0.076531)  


treeSkim("Run2011-MuTau-All_run",              0)                            
treeSkim("Run2011-MuTau-Embedded-All_run",     0)                           
treeSkim("DYJets-MuTau-50-madgraph-PUS6_run",  3048           * 0.009631    * 0.641987) 
treeSkim("TTJets-MuTau-madgraph-PUS6_run",     157.5          * 0.020998    * 0.823613)  
treeSkim("WJets-MuTau-madgraph-PUS6_run",      31314.0        * 0.001261    * 0.572776)  
treeSkim("W3Jets-MuTau-madgraph-PUS6_run",     304.0          * 1.0         * 0.119569)   
treeSkim("WZ-MuTau-pythia-PUS6_run",           18.2           * 0.0068962   * 0.740770)         
treeSkim("ZZ-MuTau-pythia-PUS6_run",            5.9           * 0.0060357   * 0.753123)    
treeSkim("WW-MuTau-pythia-PUS6_run",           43.0           * 1.0         * 0.065573)       
treeSkim("VBFH100-MuTau-powheg-PUS6_run",      8.36e-02*1.546 * 0.04953     * 0.809326) 
treeSkim("GGFH100-MuTau-powheg-PUS6_run",      8.36e-02*24.02 * 0.03816     * 0.818372) 
treeSkim("VBFH105-MuTau-powheg-PUS6_run",      8.25e-02*1.472 * 0.05278     * 0.804249) 
treeSkim("GGFH105-MuTau-powheg-PUS6_run",      8.25e-02*21.78 * 0.04137     * 0.814457) 
treeSkim("VBFH110-MuTau-powheg-PUS6_run",      8.02e-02*1.398 * 1.0         * 0.110006)
treeSkim("GGFH110-MuTau-powheg-PUS6_run",      8.02e-02*19.84 * 1.0         * 0.066984)    
treeSkim("VBFH115-MuTau-powheg-PUS6_run",      7.65e-02*1.332 * 0.05825     * 0.808179)
treeSkim("GGFH115-MuTau-powheg-PUS6_run",      7.65e-02*18.13 * 1.0         * 0.072958)  
treeSkim("VBFH120-MuTau-powheg-PUS6_run",      7.10e-02*1.269 * 1.0         * 0.118150)   
treeSkim("GGFH120-MuTau-powheg-PUS6_run",      7.10e-02*16.63 * 1.0         * 0.076531)  
treeSkim("VBFH125-MuTau-powheg-PUS6_run",      6.37e-02*1.211 * 0.06337     * 0.813953)   
treeSkim("GGFH125-MuTau-powheg-PUS6_run",      6.37e-02*15.31 * 1.0         * 0.079825) 
treeSkim("VBFH130-MuTau-powheg-PUS6_run",      5.48e-02*1.154 * 0.06515     * 0.811363) 
treeSkim("GGFH130-MuTau-powheg-PUS6_run",      5.48e-02*14.12 * 0.05479     * 0.828420)
treeSkim("VBFH135-MuTau-powheg-PUS6_run",      4.52e-02*1.100 * 0.06740     * 0.816084) 
treeSkim("GGFH135-MuTau-powheg-PUS6_run",      4.52e-02*13.08 * 1.0         * 0.088014)    
treeSkim("VBFH140-MuTau-powheg-PUS6_run",      3.54e-02*1.052 * 1.0         * 0.128525)  
treeSkim("GGFH140-MuTau-powheg-PUS6_run",      3.54e-02*12.13 * 0.05808     * 0.822400)  
treeSkim("VBFH145-MuTau-powheg-PUS6_run",      2.61e-02*1.004 * 0.07349     * 0.817493)  
treeSkim("GGFH145-MuTau-powheg-PUS6_run",      2.61e-02*11.27 * 1.0         * 0.095074) 
treeSkim("VBFH160-MuTau-powheg-PUS6_run",      5.32e-04*0.8787* 1.0         * 0.137951)  
treeSkim("GGFH160-MuTau-powheg-PUS6_run",      5.32e-04*9.080 * 1.0         * 0.104164)  
treeSkim("VH100-MuTau-pythia-PUS6_run",        2.61e-02*(1.186+ 0.6313+0.1638  ) * 0.082128 * 0.834434)  
treeSkim("VH110-MuTau-pythia-PUS6_run",        8.02e-02*(0.8754+0.4721+0.1257  ) * 0.088174 * 0.839732)
treeSkim("VH115-MuTau-pythia-PUS6_run",        7.65e-02*(0.7546+0.4107+0.1106  ) * 0.092032 * 0.841206)
treeSkim("VH120-MuTau-pythia-PUS6_run",        7.10e-02*(0.6561+0.3598+0.09756 ) * 1.0      * 0.1654)   
treeSkim("VH125-MuTau-pythia-PUS6_run",        6.37e-02*(0.5729+0.3158+0.08634 ) * 0.09885  * 0.842507)
treeSkim("VH130-MuTau-pythia-PUS6_run",        5.48e-02*(0.4942+0.2778+0.07658 ) * 1.0      * 0.177053)
treeSkim("VH135-MuTau-pythia-PUS6_run",        4.52e-02*(0.4390+0.2453+0.06810 ) * 0.104662 * 0.841589)
treeSkim("VH140-MuTau-pythia-PUS6_run",        3.54e-02*(0.3857+0.2172+0.06072 ) * 0.108278 * 0.842954)  
treeSkim("VH145-MuTau-pythia-PUS6_run",        2.61e-02*(0.3406+0.1930+0.05435 ) * 1.0      * 0.191843)  
treeSkim("VH160-MuTau-pythia-PUS6_run",        5.32e-04*(0.2291+0.1334+0.03942 ) * 1.0      * 0.203857)  
