#ifndef Bianchi_eToTaufakeRate_TauFakeRateAnalyzer_h
#define Bianchi_eToTaufakeRate_TauFakeRateAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"

#include <string>


class TauFakeRateAnalyzer : public edm::EDAnalyzer{


 public:

  explicit TauFakeRateAnalyzer(const edm::ParameterSet&);
  ~TauFakeRateAnalyzer();

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TTree* tree_;
  edm::InputTag tauTag_;
  edm::InputTag electronsTag_;
  bool isMC_;
  std::string matchTo_;
  double leadPFChargedHadrMva_;
  double leadPFChargedHadrHcalEnergy_;
  double leadPFChargedHadrEcalEnergy_;
  double leadPFChargedHadrTrackPt_;
  double leadPFChargedHadrTrackP_;
  double leadPFCandMva_;
  double leadPFCandHcalEnergy_;
  double leadPFCandEcalEnergy_;
  double leadPFCandPt_;
  double leadPFCandP_;
  double matchedID_;
  double fbrem_;
  double hasGsf_;
};


#endif
