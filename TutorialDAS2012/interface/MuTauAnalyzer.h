#ifndef Bianchi_TutorialDAS2012_MuTauAnalyzer_h
#define Bianchi_TutorialDAS2012_MuTauAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"



#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <TRandom3.h>


#include <string>
#include <utility>
#include <map>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

class MuTauAnalyzer : public edm::EDAnalyzer{


 public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit MuTauAnalyzer(const edm::ParameterSet&);
  ~MuTauAnalyzer();

  unsigned int jetID( const pat::Jet* jet );

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TFile* file_;
  TTree* tree_;

  edm::InputTag diTauTag_;
  edm::InputTag jetsTag_;
  edm::InputTag newJetsTag_;
  edm::InputTag rawMetTag_;
  edm::InputTag muonsTag_;
  edm::InputTag muonsRelTag_;
  edm::InputTag verticesTag_;
  edm::InputTag triggerResultsTag_;

  bool isMC_;
  bool verbose_;
  float minCorrPt_;
  float minJetID_;
  float deltaRLegJet_;

  std::vector< int >* tauXTriggers_;
  std::vector< int >* triggerBits_;

  std::vector< LV >* jetsIDP4_;
  std::vector< LV >* genJetsIDP4_; 
  std::vector< LV >* diTauVisP4_; 
  std::vector< LV >* diTauSVfitP4_; 
  std::vector< LV >* diTauLegsP4_; 
  std::vector< LV >* genDiTauLegsP4_; 
  std::vector< LV >* METP4_;
  std::vector< LV >* genMETP4_;
  std::vector< LV >* genVP4_;

  int genDecay_;
  unsigned long run_,event_,lumi_;
  float sumEt_;
  float chIsoLeg1v2_,nhIsoLeg1v2_,phIsoLeg1v2_,elecIsoLeg1v2_,muIsoLeg1v2_;
  float chIsoPULeg1v2_,nhIsoPULeg1v2_,phIsoPULeg1v2_;
  float chIsoLeg2_,nhIsoLeg2_,phIsoLeg2_;
  float dxy1_,dxy2_;
  float dz1_,dz2_;
  float MtLeg1_;
  float pZeta_;
  float pZetaVis_;
  float pZetaSig_;
  float mTauTauMin_;
  float numPV_;
  int numOfLooseIsoDiTaus_;
  int decayMode_;
  float diTauNSVfitMass_;
  float diTauNSVfitMassErrUp_;
  float diTauNSVfitMassErrDown_;
  float visibleTauMass_;
  int signalPFChargedHadrCands_;
  int signalPFGammaCands_;
  int tightestHPSDBWP_;
  int isTauLegMatched_;
  int isMuLegMatched_;
  int muFlag_;
  float diTauCharge_;
  float rhoFastJet_;
  int nPUVertices_;
  int nPUVerticesM1_;
  int nPUVerticesP1_;

};


#endif
