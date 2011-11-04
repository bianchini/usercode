//--------------------------------------------------------------------------------------------------
// $Id $
//
// AntiElectronIDMVA
//
// Helper Class for applying MVA anti-electron discrimination
//
// Authors: L.Bianchini
//--------------------------------------------------------------------------------------------------

#ifndef BIANCHI_UTILITIES_AntiElectronIDMVA_H
#define BIANCHI_UTILITIES_AntiElectronIDMVA_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class AntiElectronIDMVA {
  public:

    AntiElectronIDMVA();
    ~AntiElectronIDMVA(); 

    void   Initialize(std::string methodName,
                      std::string oneProng0Pi0_BL,
                      std::string oneProng1pi0wGSF_BL,
                      std::string oneProng1pi0woGSF_BL,
		      std::string oneProng0Pi0_EC,
                      std::string oneProng1pi0wGSF_EC,
                      std::string oneProng1pi0woGSF_EC
                      );

    double MVAValue(Float_t TauEta, 
		    Float_t TauSignalPFChargedCands, Float_t TauSignalPFGammaCands, 
		    Float_t TauLeadPFChargedHadrMva, Float_t TauLeadPFChargedHadrHoP, 
		    Float_t TauHasGsf, Float_t TauVisMass,  Float_t TauEmFraction,
		    Float_t GammadEta, Float_t GammadPhi, Float_t GammaPt
		    );

    /* 
    // proposed WP:    epsilonB ~ 15%  epsilonS ~ 90% wrt signal taus passing discr. ag. electrons Medium. 
    bool pass = 
    (abs(TauEta)<1.5 && TauSignalPFGammaCands==0 && MVAValue(...)>0.073) ||
    (abs(TauEta)<1.5 && TauSignalPFGammaCands>0  && TauHasGsf>0.5 && MVAValue(...)>0.091) ||
    (abs(TauEta)<1.5 && TauSignalPFGammaCands>0  && TauHasGsf<0.5 && MVAValue(...)>0.088) ||
    (abs(TauEta)>1.5 && TauSignalPFGammaCands==0 && MVAValue(...)>0.096) ||
    (abs(TauEta)>1.5 && TauSignalPFGammaCands>0  && TauHasGsf>0.5 && MVAValue(...)>0.055) ||
    (abs(TauEta)>1.5 && TauSignalPFGammaCands>0  && TauHasGsf<0.5 && MVAValue(...)>0.042);

    where:

    TauEta                  = (*taus)[i].eta();
    TauSignalPFChargedCands = (*taus)[i].signalPFChargedHadrCands().size();
    TauSignalPFGammaCands   = (*taus)[i].signalPFGammaCands().size();
    TauLeadPFChargedHadrMva = (*taus)[i].electronPreIDOutput();
    TauLeadPFChargedHadrHoP = (*taus)[i].leadPFChargedHadrCand()->hcalEnergy()/(*taus)[i].leadPFChargedHadrCand()->p();
    TauHasGsf               = ((*taus)[i].leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
    TauVisMass              = (*taus)[i].mass();
    TauEmFraction           = (*taus)[i].emFraction();
  
    gammadEta_     = new std::vector< float >();
    gammadPhi_     = new std::vector< float >();
    gammaPt_       = new std::vector< float >();
    
    for(unsigned int k = 0 ; k < ((*taus)[i].signalPFGammaCands()).size() ; k++){
    reco::PFCandidateRef gamma = ((*taus)[i].signalPFGammaCands()).at(k);
    if( ((*taus)[i].leadPFChargedHadrCand()).isNonnull() ){
    gammadEta_->push_back( gamma->eta() - (*taus)[i].leadPFChargedHadrCand()->eta() );
    gammadPhi_->push_back( gamma->phi() - (*taus)[i].leadPFChargedHadrCand()->phi() );
    }
    else{
    gammadEta_->push_back( gamma->eta() - (*taus)[i].eta() );
    gammadPhi_->push_back( gamma->phi() - (*taus)[i].phi() );
    }
    gammaPt_->push_back(  gamma->pt() );
    }
    
    GammadEta               = fabs(gammadEta_[0]);
    GammadPhi               = dPhiG:=TMath::Min(abs(gammadPhi_[0]),0.3);
    GammaPt                 = gammaPt_[0]/pt;

    */

 private:

    std::string methodName_;
    TMVA::Reader* fTMVAReader_[6];
    Float_t TauSignalPFGammaCands_; 
    Float_t TauHasGsf_;
    Float_t TauVisMass_; 
    Float_t GammadEta_; 
    Float_t GammadPhi_; 
    Float_t GammaPt_;
    Float_t TauLeadPFChargedHadrMva_;
    Float_t TauLeadPFChargedHadrHoP_;
    Float_t TauEmFraction_;
    
};

#endif
