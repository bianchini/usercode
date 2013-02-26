#include "Bianchi/TTHStudies/interface/DiJetCandidateFinder.h"

#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/Test.h"
#include "Math/GenVector/LorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <algorithm>

using namespace std;

typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;




void DiJetCandidateFinder::run(int VType , std::vector<JetByPt> jets, std::vector<LV> leptons){


  for(unsigned int i = 0; i < jets.size(); i++){
    jMapPt_[ jets[i].pt ] = i;
  }
  for(unsigned int i = 0; i < leptons.size(); i++){
    lMapPt_[ leptons[i].pt ] = i;
  }
  
  switch (VType){
  case 0:
    this->runLL(jets,leptons);
    break;
  case 1:
    this->runLL(jets,leptons);
    break;
  case 2:
    this->runLJ(jets,leptons);
    break;
  case 3:
    this->runLJ(jets,leptons);
    break;
  default:
    cout << "VType " << VType << " unsupported for DiJetCandidateFinder" << endl;
    break;
   }

  return;

}

float jetCountMult(){
  return numJetCount_;
}

float jetRecoverMult(){
  return numJetRecov_;
}

int eventFlag(){
  return eventFlag_;
}

int errorFlag(){
  return errorFlag_;
}


void DiJetCandidateFinder::runLL(std::vector<JetByPt> jets, std::vector<LV> leptons){


}


void DiJetCandidateFinder::runLJ(std::vector<JetByPt> jets, std::vector<LV> leptons){

  if(leptons.size()!=2){
    error_ = 1;
    cout << "Two leptons are needed to run DiJetCandidateFinder::runLJ: the charged letpon and the MET... return" << endl;
    return;
  }

  std::map<unsigned int, int> flagMask;
 
  std::vector<JetByPt> jetsCount;
  std::vector<JetByPt> jetsRecover;

  for(unsigned int i = 0; i < jets.size(); i++){
    bool jetCount = jets[i].pt>ptCutCount_ && TMath::Abs(jets[i].eta)<etaCutCount_;
    bool jetRecov = jets[i].pt<ptCutCount_ && jets[i].pt>ptCutRecover_ && TMath::Abs(jets[i].eta)<etaCutRecover_;
    if( jetCount ) 
      jetsCount.push_back(jets[i]); // all high pt jets
    if( jetRecov ) 
      jetsRecover.push_back(jets[i]); // jets that can be recovered
  }

  int numJetCount =   jetsCount.size();
  int numJetRecov =   jetsRecover.size();

  if( numJetCount >= 4 ){
    if(vebose_){
      cout << "The event contains " << numJetCount << " jets, proceed with findMatch..." << endl;
    }
    eventFlag_ = 0;
    bool foundMatch = this->findMatch(eventFlag_, jetsCount, leptons, flagMask);

  }
  else if( numJetRecov>0 ){
    if(vebose_){
      cout << "The event contains " << numJetCount << " jets, try to recover additional " << numJetRecov << " jets..." << endl;
    }
    eventFlag_ = 1;
    jetCount.extend( jetRecov.begin(), jetRecov.size() );
    bool foundMatch = this->findMatch(eventFlag_, jetsCount, leptons, flagMask);
  }
  else{
    if(vebose_){
      cout << "The event contains " << numJetCount << " jets, no other jets can be found, mask those available..." << endl;
    }
    eventFlag_ = 2;
  }


}


bool findMatchLJ(int eventFlag, std::vector<JetByPt>& jetsCount, std::vector<LV> leptons, std::map<unsigned int, int>& flagMask){

  std::sort(jetsCount.begin(), jetsCount.end(), sorterByJetPt);

  int counter = 0;
  if(verbose_) std::cout << "Find match... start permutations" << endl;
  do {

    counter++;

    unsigned int indexWdau1Cand = (jMapPt_.find(jetsCount[0].pt)!=jMapPt_.end()) ? jMapPt_[jetsCount[0].pt] : 999;
    unsigned int indexWdau2Cand = (jMapPt_.find(jetsCount[1].pt)!=jMapPt_.end()) ? jMapPt_[jetsCount[1].pt] : 999;
    unsigned int indexbCand     = (jMapPt_.find(jetsCount[2].pt)!=jMapPt_.end()) ? jMapPt_[jetsCount[2].pt] : 999;
    unsigned int indexbbarCand  = (jMapPt_.find(jetsCount[3].pt)!=jMapPt_.end()) ? jMapPt_[jetsCount[3].pt] : 999;

    if( indexWdau1Cand == 999 || indexWdau2Cand == 999 || indexbCand == 999 || indexbbarCand == 999){
      error_ = 1;
      cout << "Error in permutation... return FALSE" << endl;
      return false;
    }

    if(verbose_){
      cout << "Perm. num " << counter << ": " << endl;
      cout << "W(1): jet[" << indexWdau1Cand << "]" << endl;
      cout << "W(2): jet[" << indexWdau2Cand << "]" << endl;
      cout << "b:    jet[" << indexbCand << "]" << endl;
      cout << "bbar: jet[" << indexbbarCand << "]" << endl;
    }

    flagMask[indexWdau1Cand] = -99;

    LV Wdau1Cand( jetsCount[0].pt, jetsCount[0].eta, jetsCount[0].phi, jetsCount[0].mass);
    LV Wdau2Cand( jetsCount[1].pt, jetsCount[1].eta, jetsCount[1].phi, jetsCount[1].mass);
    LV bCand(     jetsCount[2].pt, jetsCount[2].eta, jetsCount[2].phi, jetsCount[2].mass);
    LV bbarCand(  jetsCount[3].pt, jetsCount[3].eta, jetsCount[3].phi, jetsCount[3].mass);
    LV chargedLep = leptons[0];
    LV met        = leptons[1];
    std::vector<LV> neutrinos   = buildNeutrinos(chargedLep, met);

    double likelihood     = 1.0;
    double WHadLikelihood = 1.0;
    double bbarLikelihood = 1.0;

    if(func1D_.find("csvLF")!=func1D_.end()){
      likelihood *= func1D_["csvLF"]->Eval( jetsCount[0].csv ); // W->LF
      likelihood *= func1D_["csvLF"]->Eval( jetsCount[1].csv ); // W->LF
      WHadLikelihood *= func1D_["csvLF"]->Eval( jetsCount[0].csv );
      WHadLikelihood *= func1D_["csvLF"]->Eval( jetsCount[1].csv );
    }
    if(func1D_.find("csvHF")!=func1D_.end()){
      likelihood *= func1D_["csvHF"]->Eval( jetsCount[2].csv ); // b tag
      likelihood *= func1D_["csvHF"]->Eval( jetsCount[3].csv ); // b tag
      bbarLikelihood *= func1D_["csvHF"]->Eval( jetsCount[3].csv ); 
    }
    if(func1D_.find("wHadMass")!=func1D_.end()){
      WHadLikelihood *= func1D_["wHadMass"]->Eval( Wdau1Cand+Wdau2Cand).M() ); // W had mass
    }
    if(func1D_.find("topLepMass")!=func1D_.end()){
      double topLepLike = TMath::Max(func1D_["topLepMass"]->Eval( (bbarCand+neutrinos[0]).M() ), // W lep mass
				     func1D_["topLepMass"]->Eval( (bbarCand+neutrinos[1]).M() ));
      likelihood *= topLepLike;
      bbarLikelihood *= topLepLike;
    }
    if(func2D_.find("topHadMass")!=func2D_.end()){
      likelihood *= func2D_["topHadMass"]->Eval( (Wdau1Cand+Wdau2Cand).M(), (Wdau1Cand+Wdau2Cand+bCand).M() ); // top vs W mass
    }

    switch(eventFlag){
    case 0: // 
      if( WHadLikelihood> 1.0 ){ // mask these jets as W candidates
	flagMask[ indexWdau1Cand ] = 1;
	flagMask[ indexWdau2Cand ] = 1;
      }
      if( bbarLikelihood > 1.0 ){ // mask these jets as bbar candidates
	flagMask[ indexbbarCand ]  = 3;
      }
      if( likelihood > 1.0 ){
	flagMask[ indexWdau1Cand ] = 1;
	flagMask[ indexWdau2Cand ] = 1;
	flagMask[ indexbCand ]     = 2;
	flagMask[ indexbbarCand ]  = 3;
      }
      break;
    default:
      break;
    }

  } while ( std::next_permutation(jetsCount.begin(), jetsCount.end(), sorterByJetPt) );

}


std::vector<LV> DiJetCandidateFinder::buildNeutrinos(LV lep, LV met){

  vector<LV> output;
  output.push_back( met);

  return output;

}
