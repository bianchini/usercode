#ifndef Bianchi_TauTauStudies_ElecTauStreamAnalyzer_h
#define Bianchi_TauTauStudies_ElecTauStreamAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <TRandom3.h>

#include "Bianchi/Utilities/interface/PUWeight.h"


#include <string>


class ElecTauStreamAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit ElecTauStreamAnalyzer(const edm::ParameterSet&);
  ~ElecTauStreamAnalyzer();

  unsigned int jetID( const pat::Jet* jet);

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;

  TRandom3* tRandom_;
 
  edm::InputTag diTauTag_;
  edm::InputTag jetsTag_;

  bool isMC_;
  bool verbose_;
  float minCorrPt_;
  float minJetID_;
  float deltaRLegJet_;

  std::vector< double >* jetsBtagHE_;
  std::vector< double >* jetsBtagHP_;
  std::vector< int >* tauXTriggers_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genJetsIDP4_; 
 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauVisP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauCAP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauICAP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauSVfit1P4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauSVfit2P4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauSVfit3P4_; 

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauLegsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genDiTauLegsP4_; 

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  >* METP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  >* genMETP4_;
  
  float run_,event_,lumi_;
  float sumEt_;
  float chIsoLeg1_,nhIsoLeg1_,phIsoLeg1_;
  float chIsoPULeg1_,nhIsoPULeg1_,phIsoPULeg1_;
  float chIsoLeg2_,nhIsoLeg2_,phIsoLeg2_;
  float dxy1_,dxy2_;
  float MtLeg1_;
  float numPV_;
  int numOfDiTaus_;
  int decayMode_;
  float visibleTauMass_;
  float leadPFChargedHadrCandTrackPt_;
  int tightestHPSWP_;
  bool applyTauSignalSel_;
  int isTauLegMatched_;
  int isElecLegMatched_;

  // ele specific variables
  float nBrehm_;
  float likelihood_;
  float nHits_;
  float sihih_;
  float dPhi_;
  float dEta_;
  float HoE_;
  float EoP_;
  float fbrem_;
  int isEleLikelihoodID_;

  float diTauCharge_;
  float rhoFastJet_;
  int nPUVertices_;
  int nOOTPUVertices_;

  PUWeight* fpuweight_;
  float mcPUweight_;

};


#endif
