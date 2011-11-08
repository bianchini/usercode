#include <TFile.h>
#include <TMath.h>
#include "Bianchi/Utilities/interface/AntiElectronIDMVA.h"


AntiElectronIDMVA::AntiElectronIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader_[i] = 0;
  }
}


AntiElectronIDMVA::~AntiElectronIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader_[i]) delete fTMVAReader_[i];
  }
}


void AntiElectronIDMVA::Initialize(std::string methodName,
				   std::string oneProng0Pi0_BL,
				   std::string oneProng1pi0wGSF_BL,
				   std::string oneProng1pi0woGSF_BL,
				   std::string oneProng0Pi0_EC,
				   std::string oneProng1pi0wGSF_EC,
				   std::string oneProng1pi0woGSF_EC
				   ){

  methodName_ = methodName;

  TMVA::Tools::Instance();

  TMVA::Reader *readerX0BL = new TMVA::Reader( "!Color:!Silent" );  
  readerX0BL->AddVariable("mva",       &TauLeadPFChargedHadrMva_);
  readerX0BL->AddVariable("HoP",       &TauLeadPFChargedHadrHoP_);
  readerX0BL->AddVariable("emFraction",&TauEmFraction_);
  readerX0BL->AddVariable("hasGsf",    &TauHasGsf_);
  readerX0BL->BookMVA( methodName_, oneProng0Pi0_BL );

  TMVA::Reader *reader11BL = new TMVA::Reader( "!Color:!Silent" );   
  reader11BL->AddVariable("mva",               &TauLeadPFChargedHadrMva_);
  reader11BL->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader11BL->AddVariable("visMass",           &TauVisMass_);
  reader11BL->AddVariable("etaMom1",           &GammadEta_);
  reader11BL->AddVariable("phiMom1",           &GammadPhi_);
  reader11BL->AddVariable("gammaFrac",         &GammadPt_);
  reader11BL->BookMVA( methodName_, oneProng1pi0wGSF_BL );

  TMVA::Reader *reader01BL = new TMVA::Reader( "!Color:!Silent" );   
  reader01BL->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader01BL->AddVariable("visMass",           &TauVisMass_);
  reader01BL->AddVariable("etaMom1",           &GammadEta_);
  reader01BL->AddVariable("phiMom1",           &GammadPhi_);
  reader01BL->AddVariable("gammaFrac",         &GammadPt_);
  reader01BL->BookMVA( methodName_, oneProng1pi0woGSF_BL ); 

  //////////////////

  TMVA::Reader *readerX0EC = new TMVA::Reader( "!Color:!Silent" );
  readerX0EC->AddVariable("mva",       &TauLeadPFChargedHadrMva_);
  readerX0EC->AddVariable("HoP",       &TauLeadPFChargedHadrHoP_);
  readerX0EC->AddVariable("emFraction",&TauEmFraction_);
  readerX0EC->AddVariable("hasGsf",    &TauHasGsf_);
  readerX0EC->BookMVA( methodName_, oneProng0Pi0_EC );

  TMVA::Reader *reader11EC = new TMVA::Reader( "!Color:!Silent" );
  reader11EC->AddVariable("mva",               &TauLeadPFChargedHadrMva_);
  reader11EC->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader11EC->AddVariable("visMass",           &TauVisMass_);
  reader11EC->AddVariable("etaMom1",           &GammadEta_);
  reader11EC->AddVariable("phiMom1",           &GammadPhi_);
  reader11EC->AddVariable("gammaFrac",         &GammadPt_);
  reader11EC->BookMVA( methodName_, oneProng1pi0wGSF_EC );

  TMVA::Reader *reader01EC = new TMVA::Reader( "!Color:!Silent" );
  reader01EC->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader01EC->AddVariable("visMass",           &TauVisMass_);
  reader01EC->AddVariable("etaMom1",           &GammadEta_);
  reader01EC->AddVariable("phiMom1",           &GammadPhi_);
  reader01EC->AddVariable("gammaFrac",         &GammadPt_);
  reader01EC->BookMVA( methodName_, oneProng1pi0woGSF_EC );


  fTMVAReader_[0] = readerX0BL;
  fTMVAReader_[1] = reader11BL;
  fTMVAReader_[2] = reader01BL;
  fTMVAReader_[3] = readerX0EC;
  fTMVAReader_[4] = reader11EC;
  fTMVAReader_[5] = reader01EC;


}


double AntiElectronIDMVA::MVAValue(Float_t TauEta, Float_t TauPt,
				   Float_t TauSignalPFChargedCands, Float_t TauSignalPFGammaCands, 
				   Float_t TauLeadPFChargedHadrMva, Float_t TauLeadPFChargedHadrHoP, 
				   Float_t TauHasGsf, Float_t TauVisMass,  Float_t TauEmFraction,
				   vector<Float_t>* GammasdEta, vector<Float_t>* GammasdPhi, vector<Float_t>* GammasPt
				   ){

  double mva;

  TauSignalPFGammaCands_   = TauSignalPFGammaCands; 
  TauHasGsf_               = TauHasGsf;
  TauVisMass_              = TauVisMass; 
  TauLeadPFChargedHadrMva_ = TMath::Max(TauLeadPFChargedHadrMva,float(-1.0));
  TauLeadPFChargedHadrHoP_ = TauLeadPFChargedHadrHoP;
  TauEmFraction_           = TMath::Max(TauEmFraction,float(0.0));

  float sumPt  = 0;
  float dEta   = 0;
  float dEta2  = 0;
  float dPhi   = 0;
  float dPhi2  = 0;
  float sumPt2 = 0;

  for(unsigned int k = 0 ; k < GammasPt->size() ; k++){
    float pt_k  = (*GammasPt)[k];
    float phi_k = (*GammasdPhi)[k];
    if ((*GammasdPhi)[k] > TMath::Pi()) phi_k = (*GammasdPhi)[k] - 2*TMath::Pi();
    else if((*GammasdPhi)[k] < -TMath::Pi()) phi_k = (*GammasdPhi)[k] + 2*TMath::Pi();
    float eta_k = (*GammasdEta)[k];
    sumPt  +=  pt_k;
    sumPt2 += (pt_k*pt_k);
    dEta   += (pt_k*eta_k);
    dEta2  += (pt_k*eta_k*eta_k);
    dPhi   += (pt_k*phi_k);
    dPhi2  += (pt_k*phi_k*phi_k);  
  }

  GammadPt_ = 0.;

  if(sumPt>0){
    dEta  /= sumPt;
    dPhi  /= sumPt;
    dEta2 /= sumPt;
    dPhi2 /= sumPt;
    GammadPt_ = sumPt/TauPt;
  }

  GammadEta_ = dEta;
  GammadPhi_ = dPhi;
  
  //GammadEta_ = TMath::Sqrt(dEta2);
  //GammadPhi_ = TMath::Sqrt(dPhi2);


  if( TauSignalPFChargedCands==3 ) 
    mva = 1.0;
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands==0){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[0]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[3]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf>0.5){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[1]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[4]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf<0.5){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[2]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[5]->EvaluateMVA( methodName_ );
  }
  else{
    mva = -1.0;
  }

  return mva;

}

double AntiElectronIDMVA::MVAValue(Float_t TauEta, Float_t TauPt,
				   Float_t TauSignalPFChargedCands, Float_t TauSignalPFGammaCands, 
				   Float_t TauLeadPFChargedHadrMva, Float_t TauLeadPFChargedHadrHoP, 
				   Float_t TauHasGsf, Float_t TauVisMass,  Float_t TauEmFraction,
				   Float_t GammadEta, Float_t GammadPhi, Float_t GammadPt
				   ){

  double mva;

  TauSignalPFGammaCands_   = TauSignalPFGammaCands; 
  TauHasGsf_               = TauHasGsf;
  TauVisMass_              = TauVisMass; 
  TauLeadPFChargedHadrMva_ = TMath::Max(TauLeadPFChargedHadrMva,float(-1.0));
  TauLeadPFChargedHadrHoP_ = TauLeadPFChargedHadrHoP;
  TauEmFraction_           = TMath::Max(TauEmFraction,float(0.0));
  GammadPt_                = GammadPt;
  GammadEta_               = GammadEta;
  GammadPhi_               = GammadPhi;
  
  if( TauSignalPFChargedCands==3 ) 
    mva = 1.0;
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands==0){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[0]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[3]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf>0.5){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[1]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[4]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf<0.5){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[2]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[5]->EvaluateMVA( methodName_ );
  }
  else{
    mva = -1.0;
  }

  return mva;

}



double AntiElectronIDMVA::MVAValue(pat::Tau* myTau){

  double mva;

  TauSignalPFGammaCands_   = (myTau->signalPFGammaCands()).size(); 
  TauHasGsf_               = ((myTau->leadPFChargedHadrCand())->gsfTrackRef()).isNonnull() ? 1. : 0.;
  TauVisMass_              = myTau->mass(); 
  TauLeadPFChargedHadrMva_ = TMath::Max(myTau->electronPreIDOutput(),float(-1.0));
  TauLeadPFChargedHadrHoP_ = (myTau->leadPFChargedHadrCand())->hcalEnergy()/myTau->leadPFChargedHadrCand()->p();
  TauEmFraction_           = TMath::Max(myTau->emFraction(),float(0.0));

  vector<float> GammasdEta;
  vector<float> GammasdPhi;
  vector<float> GammasPt;

  for(unsigned int k = 0 ; k < (myTau->signalPFGammaCands()).size() ; k++){
    reco::PFCandidateRef gamma = (myTau->signalPFGammaCands()).at(k);
    if( (myTau->leadPFChargedHadrCand()).isNonnull() ){
      GammasdEta.push_back( gamma->eta() - myTau->leadPFChargedHadrCand()->eta() );
      GammasdPhi.push_back( gamma->phi() - myTau->leadPFChargedHadrCand()->phi() );
    }
    else{
      GammasdEta.push_back( gamma->eta() - myTau->eta() );
      GammasdPhi.push_back( gamma->phi() - myTau->phi() );
    }
    GammasPt.push_back(  gamma->pt() );
  }

  float sumPt  = 0;
  float dEta   = 0;
  float dEta2  = 0;
  float dPhi   = 0;
  float dPhi2  = 0;
  float sumPt2 = 0;

  for(unsigned int k = 0 ; k < GammasPt.size() ; k++){
    float pt_k  = GammasPt[k];
    float phi_k = GammasdPhi[k];
    if (GammasdPhi[k] > TMath::Pi()) phi_k = GammasdPhi[k] - 2*TMath::Pi();
    else if(GammasdPhi[k] < -TMath::Pi()) phi_k = GammasdPhi[k] + 2*TMath::Pi();
    float eta_k = GammasdEta[k];
    sumPt  +=  pt_k;
    sumPt2 += (pt_k*pt_k);
    dEta   += (pt_k*eta_k);
    dEta2  += (pt_k*eta_k*eta_k);
    dPhi   += (pt_k*phi_k);
    dPhi2  += (pt_k*phi_k*phi_k);  
  }

  GammadPt_ = 0.;

  if(sumPt>0){
    dEta  /= sumPt;
    dPhi  /= sumPt;
    dEta2 /= sumPt;
    dPhi2 /= sumPt;
    GammadPt_ = sumPt/myTau->pt();
  }

  GammadEta_ = dEta;
  GammadPhi_ = dPhi;
  
  //GammadEta_ = TMath::Sqrt(dEta2);
  //GammadPhi_ = TMath::Sqrt(dPhi2);


  if( (myTau->signalPFChargedHadrCands()).size() == 3) 
    mva = 1.0;
  else if( (myTau->signalPFChargedHadrCands()).size()==1 && TauSignalPFGammaCands_==0){
    if(TMath::Abs(myTau->eta())<1.5) 
      mva = fTMVAReader_[0]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[3]->EvaluateMVA( methodName_ );
  }
  else if( (myTau->signalPFChargedHadrCands()).size()==1 && TauSignalPFGammaCands_>0 && TauHasGsf_>0.5){
    if(TMath::Abs(myTau->eta())<1.5) 
      mva = fTMVAReader_[1]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[4]->EvaluateMVA( methodName_ );
  }
  else if( (myTau->signalPFChargedHadrCands()).size()==1 && TauSignalPFGammaCands_>0 && TauHasGsf_<0.5){
    if(TMath::Abs(myTau->eta())<1.5) 
      mva = fTMVAReader_[2]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[5]->EvaluateMVA( methodName_ );
  }
  else{
    mva = -1.0;
  }

  return mva;

}
