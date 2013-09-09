#! /bin/sh


ls MEAnalysis_$1_p*.root
hadd -f MEAnalysis_$1.root MEAnalysis_$1_p*.root 
rm  MEAnalysis_$1_p*.root