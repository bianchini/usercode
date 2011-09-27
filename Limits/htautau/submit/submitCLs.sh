#! /bin/sh

echo "Clean the area"
rm eleTau_SM_*
rm muTau_SM_*
rm Comb_SM_*
rm *txt
rm *root

for mass in  105 110 115 120 125 130 135 140
#for mass in  105 110 115 120
do

echo ">>>>>>> Now doing eleTau"

cp ../datacards/eleTauSM_diTauNSVfitMass.root ./
cp ../datacards/eTau_SM0_mH${mass}_diTauNSVfitMass.txt ./
cp ../datacards/eTau_SM2_mH${mass}_diTauNSVfitMass.txt ./
combineCards.py Name1=eTau_SM0_mH${mass}_diTauNSVfitMass.txt Name2=eTau_SM2_mH${mass}_diTauNSVfitMass.txt  > eTau_SM_mH${mass}_diTauNSVfitMass.txt
python makeGridUsingCrab.py eTau_SM_mH${mass}_diTauNSVfitMass.txt 1 30. -n 100 -I 10 -T 400 -r -o eleTau_SM_diTauNSVfitMass_mH${mass}
crab -create -submit -cfg eleTau_SM_diTauNSVfitMass_mH${mass}.cfg -USER.ui_working_dir=crab_eleTau_SM_diTauNSVfitMass_mH${mass}

cp ../datacards/eleTauSM_diTauVisMass.root ./
cp ../datacards/eTau_SM0_mH${mass}_diTauVisMass.txt ./
cp ../datacards/eTau_SM2_mH${mass}_diTauVisMass.txt ./
combineCards.py Name1=eTau_SM0_mH${mass}_diTauVisMass.txt Name2=eTau_SM2_mH${mass}_diTauVisMass.txt  > eTau_SM_mH${mass}_diTauVisMass.txt
python makeGridUsingCrab.py eTau_SM_mH${mass}_diTauVisMass.txt   1 30. -n 100 -I 10 -T 400 -r  -o eleTau_SM_diTauVisMass_mH${mass}
crab -create -submit -cfg eleTau_SM_diTauVisMass_mH${mass}.cfg -USER.ui_working_dir=crab_eleTau_SM_diTauVisMass_mH${mass}

echo ">>>>>>> Now doing muTau"

cp ../datacards/muTauSM_diTauNSVfitMass.root ./
cp ../datacards/muTau_SM0_mH${mass}_diTauNSVfitMass.txt ./
cp ../datacards/muTau_SM2_mH${mass}_diTauNSVfitMass.txt ./
combineCards.py Name1=muTau_SM0_mH${mass}_diTauNSVfitMass.txt Name2=muTau_SM2_mH${mass}_diTauNSVfitMass.txt  > muTau_SM_mH${mass}_diTauNSVfitMass.txt
python makeGridUsingCrab.py muTau_SM_mH${mass}_diTauNSVfitMass.txt 1 30. -n 100 -I 10 -T 400 -r  -o muTau_SM_diTauNSVfitMass_mH${mass}
crab -create -submit -cfg muTau_SM_diTauNSVfitMass_mH${mass}.cfg -USER.ui_working_dir=crab_muTau_SM_diTauNSVfitMass_mH${mass}

cp ../datacards/muTauSM_diTauVisMass.root ./
cp ../datacards/muTau_SM0_mH${mass}_diTauVisMass.txt ./
cp ../datacards/muTau_SM2_mH${mass}_diTauVisMass.txt ./
combineCards.py Name1=muTau_SM0_mH${mass}_diTauVisMass.txt Name2=muTau_SM2_mH${mass}_diTauVisMass.txt  > muTau_SM_mH${mass}_diTauVisMass.txt
python makeGridUsingCrab.py muTau_SM_mH${mass}_diTauVisMass.txt   1 30. -n 100 -I 10 -T 400 -r  -o muTau_SM_diTauVisMass_mH${mass}
crab -create -submit -cfg muTau_SM_diTauVisMass_mH${mass}.cfg -USER.ui_working_dir=crab_muTau_SM_diTauVisMass_mH${mass}

done
