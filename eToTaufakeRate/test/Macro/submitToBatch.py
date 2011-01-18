#!/bin/sh

export WORKINGDIR="/home/llr/cms/lbianchini/CMSSW_3_8_6_patch1/src/PFAnalyses/VBFHTauTau/test/Macro/"
source /opt/exp_soft/cms/cmsset_default.sh
cd $WORKINGDIR
eval `scram runtime -sh`
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/slc4_ia32_gcc34/root/bin/thisroot.sh
root -b submitMacro.C
