#include <TFile.h>
#include <TMath.h>
#include "UserCode/Bianchi/Utilities/interface/AntiElectronIDMVA.h"


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
  readerX0BL->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&TauLeadPFChargedHadrMva_);
  readerX0BL->AddVariable("leadPFChargedHadrHcalEnergy/leadPFChargedHadrTrackP",&TauLeadPFChargedHadrHoP_);
  readerX0BL->AddVariable("TMath::Max(emFraction, 0.0)",&TauEmFraction_);
  readerX0BL->AddVariable("hasGsf",&TauHasGsf_);
  readerX0BL->BookMVA( methodName_, oneProng0Pi0_BL );

  TMVA::Reader *reader11BL = new TMVA::Reader( "!Color:!Silent" );   
  reader11BL->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&TauLeadPFChargedHadrMva_);
  reader11BL->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader11BL->AddVariable("visMass",&TauVisMass_);
  reader11BL->AddVariable("abs(gammadEta[0])",&GammadEta_);
  reader11BL->AddVariable("TMath::Min(abs(gammadPhi[0]),0.3)",&GammadPhi_);
  reader11BL->AddVariable("gammaPt[0]/pt",&GammaPt_);
  reader11BL->BookMVA( methodName_, oneProng1pi0wGSF_BL );

  TMVA::Reader *reader01BL = new TMVA::Reader( "!Color:!Silent" );   
  reader01BL->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader01BL->AddVariable("visMass",&TauVisMass_);
  reader01BL->AddVariable("abs(gammadEta[0])",&GammadEta_);
  reader01BL->AddVariable("TMath::Min(abs(gammadPhi[0]),0.3)",&GammadPhi_);
  reader01BL->AddVariable("gammaPt[0]/pt",&GammaPt_);
  reader01BL->BookMVA( methodName_, oneProng1pi0woGSF_BL ); 

  //////////////////

  TMVA::Reader *readerX0EC = new TMVA::Reader( "!Color:!Silent" );
  readerX0EC->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&TauLeadPFChargedHadrMva_);
  readerX0EC->AddVariable("leadPFChargedHadrHcalEnergy/leadPFChargedHadrTrackP",&TauLeadPFChargedHadrHoP_);
  readerX0EC->AddVariable("TMath::Max(emFraction, 0.0)",&TauEmFraction_);
  readerX0EC->AddVariable("hasGsf",&TauHasGsf_);
  readerX0EC->BookMVA( methodName_, oneProng0Pi0_EC );

  TMVA::Reader *reader11EC = new TMVA::Reader( "!Color:!Silent" );
  reader11EC->AddVariable("TMath::Max(leadPFChargedHadrMva,-1.0)",&TauLeadPFChargedHadrMva_);
  reader11EC->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader11EC->AddVariable("visMass",&TauVisMass_);  
  reader11EC->AddVariable("abs(gammadEta[0])",&GammadEta_);
  reader11EC->AddVariable("TMath::Min(abs(gammadPhi[0]),0.3)",&GammadPhi_);
  reader11EC->AddVariable("gammaPt[0]/pt",&GammaPt_);
  reader11EC->BookMVA( methodName_, oneProng1pi0wGSF_EC );

  TMVA::Reader *reader01EC = new TMVA::Reader( "!Color:!Silent" );
  reader01EC->AddVariable("signalPFGammaCands",&TauSignalPFGammaCands_);
  reader01EC->AddVariable("visMass",&TauVisMass_);
  reader01EC->AddVariable("abs(gammadEta[0])",&GammadEta_);
  reader01EC->AddVariable("TMath::Min(abs(gammadPhi[0]),0.3)",&GammadPhi_);
  reader01EC->AddVariable("gammaPt[0]/pt",&GammaPt_);
  reader01EC->BookMVA( methodName_, oneProng1pi0woGSF_EC );


  fTMVAReader_[0] = readerX0BL;
  fTMVAReader_[1] = reader11BL;
  fTMVAReader_[2] = reader01BL;
  fTMVAReader_[3] = readerX0EC;
  fTMVAReader_[4] = reader11EC;
  fTMVAReader_[5] = reader01EC;


}


double AntiElectronIDMVA::MVAValue(Float_t TauEta, 
				   Float_t TauSignalPFChargedCands, Float_t TauSignalPFGammaCands, 
				   Float_t TauLeadPFChargedHadrMva, Float_t TauLeadPFChargedHadrHoP, 
				   Float_t TauHasGsf, Float_t TauVisMass,  Float_t TauEmFraction,
				   Float_t GammadEta, Float_t GammadPhi, Float_t GammaPt
				   ){

  double mva;

  TauSignalPFGammaCands_   = TauSignalPFGammaCands; 
  TauHasGsf_               = TauHasGsf;
  TauVisMass_              = TauVisMass; 
  GammadEta_               = GammadEta; 
  GammadPhi_               = GammadPhi; 
  GammaPt_                 = GammaPt;
  TauLeadPFChargedHadrMva_ = TauLeadPFChargedHadrMva;
  TauLeadPFChargedHadrHoP_ = TauLeadPFChargedHadrHoP;
  TauEmFraction_           = TauEmFraction;

  if( TauSignalPFChargedCands==3 ) 
    mva = 1.0;
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands==0){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[0]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[3]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf==1){
    if(TMath::Abs(TauEta)<1.5) 
      mva = fTMVAReader_[1]->EvaluateMVA( methodName_ );
    else  
      mva = fTMVAReader_[4]->EvaluateMVA( methodName_ );
  }
  else if( TauSignalPFChargedCands==1 && TauSignalPFGammaCands>0 && TauHasGsf==0){
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
