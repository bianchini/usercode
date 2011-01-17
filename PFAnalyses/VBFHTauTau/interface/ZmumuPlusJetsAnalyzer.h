#ifndef PFAnalyses_VBFHTauTau_ZmumuPlusJetsAnalyzer_h
#define PFAnalyses_VBFHTauTau_ZmumuPlusJetsAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "AnalysisDataFormats/PFAnalyses/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/PFAnalyses/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <string>


class ZmumuPlusJetsAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit ZmumuPlusJetsAnalyzer(const edm::ParameterSet&);
  ~ZmumuPlusJetsAnalyzer();

  unsigned int jetID( const pat::Jet* jet, const double scale);

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;
  std::string fileName_;
  edm::InputTag diMuonTag_;
  edm::InputTag jetsTag_;
  const reco::GenParticle *tagQuark1_, *tagQuark2_;

  bool isMC_;

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* tagJetsIDP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDbyMjjP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDbyDEtaP4_; 

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* muonsP4_; 

  float chIsoLeg1_,nhIsoLeg1_,phIsoLeg1_;
  float chIsoLeg2_,nhIsoLeg2_,phIsoLeg2_;
  float run_,event_;
  float Zmass_,ZdeltaPhi_,MET_;
  float hltMu7_,hltMu9_,hltMu11_,hltMu15v1_,hltIsoMu13v3_,hltIsoMu13v4_;
  float numPV_;
  float dxy1_,dxy2_;

};


#endif
