#! /bin/sh

ls -ltr

#./hadd.sh  SL_$1_$2_TTJetsSemiLept
./hadd.sh  SL_$1_$2_TTJetsFullLept
./hadd.sh  SL_$1_$2_TTH125
./hadd.sh  SL_$1_$2_TTZ
./hadd.sh  SL_$1_$2_DYJets10to50
./hadd.sh  SL_$1_$2_DYJets50
./hadd.sh  SL_$1_$2_TTW
./hadd.sh  SL_$1_$2_Tt
./hadd.sh  SL_$1_$2_WJets
./hadd.sh  SL_$1_$2_Ts
./hadd.sh  SL_$1_$2_ZZ
./hadd.sh  SL_$1_$2_Tbars
./hadd.sh  SL_$1_$2_WW
./hadd.sh  SL_$1_$2_WZ
./hadd.sh  SL_$1_$2_Tbart
./hadd.sh  SL_$1_$2_TbartW
./hadd.sh  SL_$1_$2_TtW
 

hadd -f MEAnalysis_SL_$1_$2_EWK.root  MEAnalysis_SL_$1_$2_DYJets10to50.root MEAnalysis_SL_$1_$2_DYJets50.root MEAnalysis_SL_$1_$2_WJets.root
rm   MEAnalysis_SL_$1_$2_DYJets10to50.root MEAnalysis_SL_$1_$2_DYJets50.root MEAnalysis_SL_$1_$2_WJets.root

hadd -f  MEAnalysis_SL_$1_$2_SingleT.root MEAnalysis_SL_$1_$2_Tt.root MEAnalysis_SL_$1_$2_Ts.root MEAnalysis_SL_$1_$2_Tbars.root MEAnalysis_SL_$1_$2_Tbart.root MEAnalysis_SL_$1_$2_TbartW.root MEAnalysis_SL_$1_$2_TtW.root
rm MEAnalysis_SL_$1_$2_Tt.root MEAnalysis_SL_$1_$2_Ts.root MEAnalysis_SL_$1_$2_Tbars.root MEAnalysis_SL_$1_$2_Tbart.root MEAnalysis_SL_$1_$2_TbartW.root MEAnalysis_SL_$1_$2_TtW.root

hadd -f  MEAnalysis_SL_$1_$2_DiBoson.root MEAnalysis_SL_$1_$2_ZZ.root MEAnalysis_SL_$1_$2_WW.root MEAnalysis_SL_$1_$2_WZ.root
rm  MEAnalysis_SL_$1_$2_ZZ.root MEAnalysis_SL_$1_$2_WW.root MEAnalysis_SL_$1_$2_WZ.root

hadd -f  MEAnalysis_SL_$1_$2_TTV.root MEAnalysis_SL_$1_$2_TTZ.root MEAnalysis_SL_$1_$2_TTW.root
rm MEAnalysis_SL_$1_$2_TTZ.root MEAnalysis_SL_$1_$2_TTW.root


# syst
hadd -f MEAnalysis_SL_$1_csvDown_v5_TTH125.root MEAnalysis_SL_$1_csvDown_v5_TTH125_p*.root
rm MEAnalysis_SL_$1_csvDown_v5_TTH125_p*.root
hadd -f MEAnalysis_SL_$1_csvUp_v5_TTH125.root MEAnalysis_SL_$1_csvUp_v5_TTH125_p*.root
rm MEAnalysis_SL_$1_csvUp_v5_TTH125_p*.root
hadd -f MEAnalysis_SL_$1_JECDown_v5_TTH125.root MEAnalysis_SL_$1_JECDown_v5_TTH125_p*.root
rm MEAnalysis_SL_$1_JECDown_v5_TTH125_p*.root
hadd -f MEAnalysis_SL_$1_JECUp_v5_TTH125.root MEAnalysis_SL_$1_JECUp_v5_TTH125_p*.root
rm MEAnalysis_SL_$1_JECUp_v5_TTH125_p*.root
hadd -f MEAnalysis_SL_$1_csvDown_v5_TTJetsFullLept.root MEAnalysis_SL_$1_csvDown_v5_TTJetsFullLept_p*.root
rm MEAnalysis_SL_$1_csvDown_v5_TTJetsFullLept_p*.root
hadd -f MEAnalysis_SL_$1_csvUp_v5_TTJetsFullLept.root MEAnalysis_SL_$1_csvUp_v5_TTJetsFullLept_p*.root
rm MEAnalysis_SL_$1_csvUp_v5_TTJetsFullLept_p*.root
hadd -f MEAnalysis_SL_$1_JECDown_v5_TTJetsFullLept.root MEAnalysis_SL_$1_JECDown_v5_TTJetsFullLept_p*.root
rm MEAnalysis_SL_$1_JECDown_v5_TTJetsFullLept_p*.root
hadd -f MEAnalysis_SL_$1_JECUp_v5_TTJetsFullLept.root MEAnalysis_SL_$1_JECUp_v5_TTJetsFullLept_p*.root
rm MEAnalysis_SL_$1_JECUp_v5_TTJetsFullLept_p*.root