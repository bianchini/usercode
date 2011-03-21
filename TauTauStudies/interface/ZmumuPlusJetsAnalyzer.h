#ifndef Bianchi_TauTauStudies_ZmumuPlusJetsAnalyzer_h
#define Bianchi_TauTauStudies_ZmumuPlusJetsAnalyzer_h

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

#include <string>


class ZmumuPlusJetsAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit ZmumuPlusJetsAnalyzer(const edm::ParameterSet&);
  ~ZmumuPlusJetsAnalyzer();

  unsigned int jetID( const pat::Jet* jet);

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;
  edm::InputTag diMuonTag_;
  edm::InputTag jetsTag_;

  bool isMC_;
  bool verbose_;
  float minCorrPt_;
  float minJetID_;

  std::vector< double >* jetsBtagHE_;
  std::vector< double >* jetsBtagHP_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDbyMjjP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* jetsIDbyDEtaP4_; 

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* muonsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diMuonP4_; 

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  >* METP4_;
  
  float run_,event_;
  float Zmass_,ZdeltaPhi_,sumEt_;
  float chIsoLeg1_,nhIsoLeg1_,phIsoLeg1_;
  float chIsoLeg2_,nhIsoLeg2_,phIsoLeg2_;
  float dxy1_,dxy2_;
  float MtLeg1_,MtLeg2_;
  float numPV_;

};


#endif
