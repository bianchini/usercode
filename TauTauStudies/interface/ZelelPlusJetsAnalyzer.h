#ifndef Bianchi_TauTauStudies_ZelelPlusJetsAnalyzer_h
#define Bianchi_TauTauStudies_ZelelPlusJetsAnalyzer_h

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
#include "TMath.h"

#include <string>
#include "Bianchi/Utilities/interface/PUWeight.h"


class ZelelPlusJetsAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit ZelelPlusJetsAnalyzer(const edm::ParameterSet&);
  ~ZelelPlusJetsAnalyzer();

  unsigned int jetID( const pat::Jet* jet);

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;
  edm::InputTag diElectronTag_;
  edm::InputTag jetsTag_;
  edm::InputTag triggerResultsTag_;


  bool isMC_;
  bool verbose_;
  float minCorrPt_;
  float minJetID_;
  double deltaRLegJet_;

  std::vector< double >* jetsBtagHE_;
  std::vector< double >* jetsBtagHP_;
  std::vector< int >* triggerBits_;

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genJetsIDP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diElectronLegsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  >* METP4_;

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* extraElectrons_; 

  float run_,event_,lumi_;
  float Zmass_,sumEt_;
  float chIsoLeg1_,nhIsoLeg1_,phIsoLeg1_;
  float chIsoPULeg1_,nhIsoPULeg1_,phIsoPULeg1_;
  float chIsoLeg2_,nhIsoLeg2_,phIsoLeg2_;
  float chIsoPULeg2_,nhIsoPULeg2_,phIsoPULeg2_;
  float dxy1_,dxy2_;
  float dz1_,dz2_;
  float MtLeg1_,MtLeg2_;

  int isLegFromTau_;
  float numPV_;
  float rhoFastJet_;
  int nPUVertices_;
  int nOOTPUVertices_;
  PUWeight* fpuweight_;
  float mcPUweight_;

};


#endif
