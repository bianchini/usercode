#include "Bianchi/TauTauStudies/interface/MuTauStreamAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <vector>
#include <utility>
#include <map>

using namespace std;
using namespace reco;

typedef std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more>::iterator CImap;

MuTauStreamAnalyzer::MuTauStreamAnalyzer(const edm::ParameterSet & iConfig){

  diTauTag_ = iConfig.getParameter<edm::InputTag>("diTaus");
  jetsTag_ = iConfig.getParameter<edm::InputTag>("jets");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults"); 
  isMC_ =  iConfig.getParameter<bool>("isMC");
  deltaRLegJet_  =  iConfig.getUntrackedParameter<double>("deltaRLegJet",0.3);
  minCorrPt_  =  iConfig.getUntrackedParameter<double>("minCorrPt",10.);
  minJetID_   =  iConfig.getUntrackedParameter<double>("minJetID",0.5);
  verbose_ =  iConfig.getUntrackedParameter<bool>("verbose",false);
}

void MuTauStreamAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","qqH tree");

  tRandom_ = new TRandom3();
 
  jetsBtagHE_  = new std::vector< double >();
  jetsBtagHP_  = new std::vector< double >();

  jetsChEfraction_  = new std::vector< float >();
  jetsChNfraction_  = new std::vector< float >();
  jetMoments_       = new std::vector< float >();

  gammadR_       = new std::vector< float >();
  gammaPt_       = new std::vector< float >();

  tauXTriggers_= new std::vector< int >();
  triggerBits_ = new std::vector< int >();

  jetsP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDP4_        = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDUpP4_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDDownP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDL1OffsetP4_= new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genJetsIDP4_     = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauVisP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauCAP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauICAP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfitP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genDiTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  METP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genMETP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genVP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  extraMuons_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  fpuweight_ = 0;//new PUWeight();

  weights2011_.push_back(0.0809561); 
  weights2011_.push_back(0.483674); 
  weights2011_.push_back(1.04601); 
  weights2011_.push_back(1.6949); 
  weights2011_.push_back(2.1009); 
  weights2011_.push_back(2.24013); 
  weights2011_.push_back(2.15938); 
  weights2011_.push_back(1.73219); 
  weights2011_.push_back(1.32508); 
  weights2011_.push_back(0.950115);
  weights2011_.push_back(0.6599); 
  weights2011_.push_back(0.438093); 
  weights2011_.push_back(0.284885); 
  weights2011_.push_back(0.180615); 
  weights2011_.push_back(0.112314); 
  weights2011_.push_back(0.0687845);
  weights2011_.push_back(0.0415087); 
  weights2011_.push_back(0.0247498); 
  weights2011_.push_back(0.0147688); 
  weights2011_.push_back(0.00865185); 
  weights2011_.push_back(0.00512891); 
  weights2011_.push_back(0.00291561); 
  weights2011_.push_back(0.00169447); 
  weights2011_.push_back(0.000998808);
  weights2011_.push_back(0.00087273); 
  weights2011_.push_back(0); 
  weights2011_.push_back(0); 
  

  tree_->Branch("jetsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsP4_);
  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("jetsIDUpP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDUpP4_);
  tree_->Branch("jetsIDDownP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDDownP4_);
  tree_->Branch("jetsIDL1OffsetP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDL1OffsetP4_);
  tree_->Branch("genJetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genJetsIDP4_);
  
  tree_->Branch("jetsBtagHE","std::vector<double>",&jetsBtagHE_);
  tree_->Branch("jetsBtagHP","std::vector<double>",&jetsBtagHP_);

  tree_->Branch("jetMoments","std::vector<float> ",&jetMoments_);

  tree_->Branch("gammadR","std::vector<float> ",&gammadR_);
  tree_->Branch("gammaPt","std::vector<float> ",&gammaPt_);

  tree_->Branch("jetsChEfraction","std::vector<float>",&jetsChEfraction_);
  tree_->Branch("jetsChNfraction","std::vector<float>",&jetsChNfraction_);

  tree_->Branch("tauXTriggers","std::vector<int>",&tauXTriggers_);
  tree_->Branch("triggerBits","std::vector<int>",&triggerBits_);

  tree_->Branch("diTauVisP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",  &diTauVisP4_);
  tree_->Branch("diTauCAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauCAP4_);
  tree_->Branch("diTauICAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauICAP4_);
  tree_->Branch("diTauSVfitP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfitP4_);

  tree_->Branch("diTauNSVfitMass",       &diTauNSVfitMass_,"diTauNSVfitMass/F");
  tree_->Branch("diTauNSVfitMassErrUp",  &diTauNSVfitMassErrUp_,"diTauNSVfitMassErrUp/F");
  tree_->Branch("diTauNSVfitMassErrDown",&diTauNSVfitMassErrDown_,"diTauNSVfitMassErrDown/F");

  tree_->Branch("diTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsP4_);
  tree_->Branch("genDiTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genDiTauLegsP4_);

  tree_->Branch("METP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&METP4_);
  tree_->Branch("genMETP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genMETP4_);
  tree_->Branch("genVP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genVP4_);
  tree_->Branch("genDecay",&genDecay_,"genDecay/I");

  tree_->Branch("extraMuons","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&extraMuons_);

  tree_->Branch("sumEt",&sumEt_,"sumEt/F");
  tree_->Branch("mTauTauMin",&mTauTauMin_,"mTauTauMin/F");
  tree_->Branch("MtLeg1",&MtLeg1_,"MtLeg1/F");
  tree_->Branch("pZeta",&pZeta_,"pZeta/F");
  tree_->Branch("pZetaVis",&pZetaVis_,"pZetaVis/F");
  tree_->Branch("pZetaSig",&pZetaSig_,"pZetaSig/F");

  tree_->Branch("chIsoLeg1v1",&chIsoLeg1v1_,"chIsoLeg1v1/F");
  tree_->Branch("nhIsoLeg1v1",&nhIsoLeg1v1_,"nhIsoLeg1v1/F");
  tree_->Branch("phIsoLeg1v1",&phIsoLeg1v1_,"phIsoLeg1v1/F");
  tree_->Branch("chIsoPULeg1v1",&chIsoPULeg1v1_,"chIsoPULeg1v1/F");
  tree_->Branch("nhIsoPULeg1v1",&nhIsoPULeg1v1_,"nhIsoPULeg1v1/F");
  tree_->Branch("phIsoPULeg1v1",&phIsoPULeg1v1_,"phIsoPULeg1v1/F");

  tree_->Branch("chIsoLeg1v2",&chIsoLeg1v2_,"chIsoLeg1v2/F");
  tree_->Branch("nhIsoLeg1v2",&nhIsoLeg1v2_,"nhIsoLeg1v2/F");
  tree_->Branch("phIsoLeg1v2",&phIsoLeg1v2_,"phIsoLeg1v2/F");
  tree_->Branch("chIsoPULeg1v2",&chIsoPULeg1v2_,"chIsoPULeg1v2/F");
  tree_->Branch("nhIsoPULeg1v2",&nhIsoPULeg1v2_,"nhIsoPULeg1v2/F");
  tree_->Branch("phIsoPULeg1v2",&phIsoPULeg1v2_,"phIsoPULeg1v2/F");

  tree_->Branch("chIsoLeg2",&chIsoLeg2_,"chIsoLeg2/F");
  tree_->Branch("nhIsoLeg2",&nhIsoLeg2_,"nhIsoLeg2/F");
  tree_->Branch("phIsoLeg2",&phIsoLeg2_,"phIsoLeg2/F");
  tree_->Branch("dxy1",&dxy1_,"dxy1/F");
  tree_->Branch("dxy2",&dxy2_,"dxy2/F");
  tree_->Branch("dz1",&dz1_,"dz1/F");
  tree_->Branch("dz2",&dz2_,"dz2/F");

  tree_->Branch("run",&run_,"run/l");
  tree_->Branch("event",&event_,"event/l");
  tree_->Branch("lumi",&lumi_,"lumi/l");
  tree_->Branch("numPV",&numPV_,"numPV/F");
  tree_->Branch("numOfDiTaus",&numOfDiTaus_,"numOfDiTaus/I");
  tree_->Branch("numOfLooseIsoDiTaus",&numOfLooseIsoDiTaus_,"numOfLooseIsoDiTaus/I");
  tree_->Branch("decayMode",&decayMode_,"decayMode/I");
  tree_->Branch("tightestHPSWP",&tightestHPSWP_,"tightestHPSWP/I");
  tree_->Branch("tightestHPSDBWP",&tightestHPSDBWP_,"tightestHPSDBWP/I");
  tree_->Branch("visibleTauMass",&visibleTauMass_,"visibleTauMass/F");
  tree_->Branch("leadPFChargedHadrCandTrackPt",&leadPFChargedHadrCandTrackPt_,"leadPFChargedHadrCandTrackPt/F");
  
  tree_->Branch("isTauLegMatched",&isTauLegMatched_,"isTauLegMatched/I");
  tree_->Branch("isMuLegMatched",&isMuLegMatched_,"isMuLegMatched/I");
  tree_->Branch("muFlag",&muFlag_,"muFlag/I");
  tree_->Branch("hasKft",&hasKft_,"hasKft/I");

  tree_->Branch("diTauCharge",&diTauCharge_,"diTauCharge/F");
  tree_->Branch("rhoFastJet",&rhoFastJet_,"rhoFastJet/F");
  tree_->Branch("rhoNeutralFastJet",&rhoNeutralFastJet_,"rhoNeutralFastJet/F");
  tree_->Branch("mcPUweight",&mcPUweight_,"mcPUweight/F");
  tree_->Branch("nPUVertices",&nPUVertices_,"nPUVertices/I");


}


MuTauStreamAnalyzer::~MuTauStreamAnalyzer(){
  delete jetsP4_; delete jetsIDP4_; delete jetsIDUpP4_; delete jetsIDDownP4_; 
  delete METP4_; delete diTauVisP4_; delete diTauCAP4_; delete diTauICAP4_; 
  delete diTauSVfitP4_; delete genVP4_;
  delete diTauLegsP4_; delete jetsBtagHE_; delete jetsBtagHP_; delete tauXTriggers_; delete triggerBits_;
  delete genJetsIDP4_; delete genDiTauLegsP4_; delete genMETP4_; delete extraMuons_; delete jetsIDL1OffsetP4_;
  delete tRandom_ ; /*delete fpuweight_*/ delete jetsChNfraction_; delete jetsChEfraction_;delete jetMoments_;
  delete gammadR_ ; delete gammaPt_;
}

void MuTauStreamAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){



  jetsP4_->clear();
  jetsIDP4_->clear();
  jetsIDUpP4_->clear();
  jetsIDDownP4_->clear();
  jetsIDL1OffsetP4_->clear();
  diTauVisP4_->clear();
  diTauCAP4_->clear();
  diTauICAP4_->clear();
  diTauSVfitP4_->clear();
  diTauLegsP4_->clear();
  METP4_->clear();
  genVP4_->clear();
  jetsChNfraction_->clear();
  jetsChEfraction_->clear();
  jetMoments_->clear();
  
  gammadR_->clear();
  gammaPt_->clear();

  genJetsIDP4_->clear();
  genDiTauLegsP4_->clear();
  genMETP4_->clear();
  extraMuons_->clear();
  jetsBtagHE_->clear();
  jetsBtagHP_->clear();
  tauXTriggers_->clear();
  triggerBits_->clear();
  
  edm::Handle<PATMuTauPairCollection> diTauHandle;
  iEvent.getByLabel(diTauTag_,diTauHandle);
  if( !diTauHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No diTau label available \n";
  const PATMuTauPairCollection* diTaus = diTauHandle.product();

  edm::Handle<pat::JetCollection> jetsHandle;
  iEvent.getByLabel(jetsTag_,jetsHandle);
  if( !jetsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No jets label available \n";
  const pat::JetCollection* jets = jetsHandle.product();

  edm::Handle<reco::VertexCollection> pvHandle;
  edm::InputTag pvTag("offlinePrimaryVertices");
  iEvent.getByLabel(pvTag,pvHandle);
  if( !pvHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No PV label available \n";
  const reco::VertexCollection* vertexes = pvHandle.product();

  std::vector<float> vtxZ;
  for(unsigned int k = 0; k<vertexes->size(); k++){
    vtxZ.push_back(((*vertexes)[k].position()).z());
  }


  if(verbose_){
    cout <<  "Run " << iEvent.run() << ", event " << (iEvent.eventAuxiliary()).event() 
	 << ", lumi " << iEvent.luminosityBlock() << endl;
    cout << "List of vertexes " << endl;
    for(unsigned int k = 0; k<vertexes->size(); k++){
      cout << "Vtx[" << k << "] (x,y,z) = (" << ((*vertexes)[k].position()).x()
	   << "," << ((*vertexes)[k].position()).y() << "," << ((*vertexes)[k].position()).z() << ")"
	   << endl;
    }
  }

  numPV_ = vertexes->size();

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel(edm::InputTag("patMETsPFlow"),metHandle);
  if( !metHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No MET label available \n";
  const pat::METCollection* met = metHandle.product();

  edm::Handle<pat::TriggerEvent> triggerHandle;
  iEvent.getByLabel(triggerResultsTag_, triggerHandle);
  if( !triggerHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No Trigger label available \n";
  const pat::TriggerEvent* trigger = triggerHandle.product();

  edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjsHandle;
  iEvent.getByLabel(edm::InputTag("patTrigger"),triggerObjsHandle);
  if( !triggerObjsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No Trigger Objects label available \n";
  const pat::TriggerObjectStandAloneCollection* triggerObjs = triggerObjsHandle.product();


  genDecay_ = -99;
  edm::Handle<reco::GenParticleCollection> genHandle;
  const reco::GenParticleCollection* genParticles = 0;
  if(isMC_){
    iEvent.getByLabel(edm::InputTag("genParticles"),genHandle);
    if( !genHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen particles label available \n";
    genParticles = genHandle.product();
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( !( (*genParticles)[k].pdgId() == 23 || abs((*genParticles)[k].pdgId()) == 24 || (*genParticles)[k].pdgId() == 25) || (*genParticles)[k].status()!=3)
	continue;
      if(verbose_) cout << "Boson status, pt,phi " << (*genParticles)[k].status() << "," << (*genParticles)[k].pt() << "," << (*genParticles)[k].phi() << endl;

      genVP4_->push_back( (*genParticles)[k].p4() );
      genDecay_ = (*genParticles)[k].pdgId();

      int breakLoop = 0;
      for(unsigned j = 0; j< ((*genParticles)[k].daughterRefVector()).size() && breakLoop!=1; j++){
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 11 ){
	  genDecay_ *= 11;
	  breakLoop = 1;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 13 ){
	  genDecay_ *= 13;
	  breakLoop = 1;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 15 ){
	  genDecay_ *= 15;
	  breakLoop = 1;
	}
      }
      if(verbose_) cout << "Decays to pdgId " << genDecay_/(*genParticles)[k].pdgId()  << endl;
      break;
    }
  }

  edm::Handle<reco::GenJetCollection> tauGenJetsHandle;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
  nPUVertices_= -99;
  nOOTPUVertices_=-99;
  float sum_nvtx = 0;

  const reco::GenJetCollection* tauGenJets = 0;
  if(isMC_){
    // tag gen jets
    iEvent.getByLabel(edm::InputTag("tauGenJetsSelectorAllHadrons"),tauGenJetsHandle);
    if( !tauGenJetsHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen jet label available \n";
    tauGenJets = tauGenJetsHandle.product();

    // PU infos
    iEvent.getByType(puInfoH);
    if(puInfoH.isValid()){
      for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
	//cout << "Bunc crossing " << it->getBunchCrossing() << endl;
	if(it->getBunchCrossing() ==0) nPUVertices_ = it->getPU_NumInteractions();
	else  nOOTPUVertices_ = it->getPU_NumInteractions();
	sum_nvtx += float(it->getPU_NumInteractions());
      }
    }
  }
  if(verbose_){
    cout << "Average num of int = " << sum_nvtx/3. << endl;
    cout << "Num of PU = " << nPUVertices_ << endl;
    cout << "Num of OOT PU = " << nOOTPUVertices_ << endl;
  }
  if( fpuweight_ ) mcPUweight_ = fpuweight_->GetWeight(nPUVertices_);
  //else mcPUweight_ =  int(sum_nvtx/3.+0.5) < int(weights2011_.size()) ? weights2011_[int(sum_nvtx/3.+0.5)] : weights2011_[(unsigned int)(weights2011_.size()-1)];
  else mcPUweight_ =  int(nPUVertices_+0.5) < int(weights2011_.size()) ? weights2011_[int(nPUVertices_+0.5)] : weights2011_[(unsigned int)(weights2011_.size()-1)];


  edm::Handle<double> rhoFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentral","rho", ""), rhoFastJetHandle);
  if( !rhoFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho label available \n";
  rhoFastJet_ = (*rhoFastJetHandle);

  edm::Handle<double> rhoNeutralFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsNeutral","rho", ""), rhoNeutralFastJetHandle);
  if( !rhoNeutralFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho neutral label available \n";
  rhoNeutralFastJet_ = (*rhoNeutralFastJetHandle);
  
  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel("muPtEtaID",muonsHandle);
  if( !muonsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No muons label available \n";
  const pat::MuonCollection* muons = muonsHandle.product();
  if(muons->size()<1){
    cout << " No muons !!! " << endl;
    return;
  } else if(muons->size()>1 && verbose_){
    cout << "WARNING: "<< muons->size() << "  muons found in the event !!! We will select only one" << endl;
  }
  edm::Handle<pat::MuonCollection> muonsRelHandle;
  iEvent.getByLabel("muPtEtaRelID",muonsRelHandle);
  if( !muonsRelHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No muonsRel label available \n";
  const pat::MuonCollection* muonsRel = muonsRelHandle.product();
  if(muonsRel->size()<1){
    cout << " No muonsRel !!! " << endl;
    return;
  } else if(muonsRel->size()>1 && verbose_){
    cout << "WARNING: "<< muonsRel->size() << "  muonsRel found in the event !!! We will select only one" << endl;
  }
  
  const PATMuTauPair *theDiTau = 0;
  if(diTaus->size()<1){
    cout << " No diTau !!! " << endl;
    return;
  } else if(diTaus->size()>1 && verbose_){
    cout << "WARNING: "<< diTaus->size() << "  diTaus found in the event !!! We will select only one" << endl;
  }
  numOfDiTaus_ = diTaus->size();

  unsigned int muIndex = 0;
  muFlag_=-99;

  if(  muons->size()==1 &&  muonsRel->size()==1 ){
    //muIndex = 0;
    muFlag_ = 0;
    if(verbose_) cout<< "Just one selected muon: flag= " << muFlag_ << endl;
  }
  else if(  muonsRel->size()>1 ){
    //muIndex = 0;
    muFlag_ = 3;
  
    bool found = false;
    for(unsigned int i=0; i<muonsRel->size(); i++){
      for(unsigned int j=0; j<muons->size(); j++){

	if(found) continue;

	if( Geom::deltaR( (*muonsRel)[i].p4(),(*muons)[j].p4())>0.3 
	    && (*muonsRel)[i].charge()*(*muons)[j].charge()<0
	    && (*muons)[j].userFloat("PFRelIsoDB04v2")<0.3 && (*muonsRel)[i].userFloat("PFRelIsoDB04v2")<0.3 ){
	  //muIndex = 0;
	  muFlag_ = 1;
	  if(verbose_) cout<< "Two muons failing diMu veto: flag= " << muFlag_ << endl;
	  found=true; 
	}
	else if( Geom::deltaR( (*muonsRel)[i].p4(),(*muons)[j].p4())>0.3 
		 && (*muonsRel)[i].charge()*(*muons)[j].charge()>0
		 && (*muons)[j].userFloat("PFRelIsoDB04v2")<0.3 && (*muonsRel)[i].userFloat("PFRelIsoDB04v2")<0.3 ){
	  //muIndex = 0;
	  muFlag_ = 2;
	  if(verbose_) cout<< "Two muons with SS: flag= " << muFlag_ << endl;
	  found=true;  
	}
      }
    }
  } 
  else{
    //muIndex = 0;
    muFlag_ = -1;
  }

  // select highest-pt electron that forms at least one di-tau 
  bool found2 = false;
  for(unsigned int j=0; j<muons->size(); j++){
    for(unsigned int i=0; i< diTaus->size(); i++){

      if(found2) continue;

      if( Geom::deltaR(((*diTaus)[i].leg1())->p4(),((*muons)[j]).p4())<0.01 ){
	muIndex=j;
	found2 = true;
      }
    }
  }
  if(verbose_) cout<< "Selected muon with index " << muIndex << endl;

  for(unsigned int i=0; i<muonsRel->size(); i++){
    if( Geom::deltaR( (*muonsRel)[i].p4(),(*muons)[muIndex].p4())>0.3) 
      extraMuons_->push_back( (*muonsRel)[i].p4() );
  }

  std::vector<unsigned int> selectedDiTausFromMu;
  for(unsigned int i=0; i< diTaus->size(); i++){
    if( Geom::deltaR(((*diTaus)[i].leg1())->p4(),((*muons)[muIndex]).p4())<0.01 ) 
      selectedDiTausFromMu.push_back(i);
  }
  if(verbose_) cout<< "Selection on Muons has selected " << selectedDiTausFromMu.size() << " diTaus...check tau leg now..." << endl;
  
  
  // choose the diTau with most isolated tau leg
  double sumIsoTau = 999.;
  //double highestPt = 0.;
  unsigned int index = 0;

  std::vector<unsigned int> identifiedTaus;
  std::vector<unsigned int> looseIsoTaus;
  //std::vector<unsigned int> not_identifiedTaus;

  for(unsigned int i=0; i<selectedDiTausFromMu.size(); i++){
    const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[selectedDiTausFromMu[i]].leg2()).get() );
    if(tau_i->tauID("decayModeFinding")<0.5) continue;
    identifiedTaus.push_back(selectedDiTausFromMu[i]);
    bool isIsolatedTau = false;
    //isIsolatedTau = ( (tau_i->isolationPFChargedHadrCandsPtSum() +
    //	       std::max( (tau_i->isolationPFGammaCandsEtSum() -
    //			  rhoFastJet_*TMath::Pi()*0.5*0.5), 0.0) )<2.); //new setting
    isIsolatedTau = tau_i->tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5;
    if(isIsolatedTau) looseIsoTaus.push_back(selectedDiTausFromMu[i]);
  }

  numOfLooseIsoDiTaus_ = looseIsoTaus.size();
  if(looseIsoTaus.size()>0 ){
    identifiedTaus.swap(looseIsoTaus);
    if(verbose_) cout << identifiedTaus.size() << "  isolated taus found..." << endl;
    for(unsigned int i=0; i<identifiedTaus.size(); i++){
      if(verbose_) cout << "Testing isolation of " << i << "th tau" << endl;
      const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[ identifiedTaus[i] ].leg2()).get() );
      double sumIsoTau_i = 0.;
      sumIsoTau_i += tau_i->isolationPFChargedHadrCandsPtSum();
      //sumIsoTau_i += tau_i->isolationPFGammaCandsEtSum();
      sumIsoTau_i += std::max( (tau_i->isolationPFGammaCandsEtSum() -
				rhoFastJet_*TMath::Pi()*0.5*0.5), 0.0);
      //sumIsoTau_i += tau_i->isolationPFNeutrHadrCandsEtSum();
      if(sumIsoTau_i<sumIsoTau){
	index = identifiedTaus[i];
	sumIsoTau = sumIsoTau_i;
      } 
    }
  } 
  else if(identifiedTaus.size()>0 ) {
    index = identifiedTaus[tRandom_->Integer( identifiedTaus.size() )];
    if(verbose_) cout << "Random selection has chosen index " << index << endl;   
  }

   
  if(verbose_) cout << "Chosen index " << index << endl;
  identifiedTaus.clear(); /*not_identifiedTaus.clear(); */ looseIsoTaus.clear();

  theDiTau = &(*diTaus)[index];

  
  diTauCharge_ = theDiTau->charge();
  METP4_->push_back((*met)[0].p4());
  if(isMC_) genMETP4_->push_back( (*met)[0].genMET()->p4() );
  sumEt_  = (*met)[0].sumEt();
  MtLeg1_ =  theDiTau->mt1MET();
  pZeta_     =  theDiTau->pZeta();
  pZetaVis_  =  theDiTau->pZetaVis();
  pZetaSig_  =  theDiTau->pZetaSig();
  mTauTauMin_=  theDiTau->mTauTauMin();

  isMuLegMatched_  = 0;
  isTauLegMatched_ = 0;

  const pat::Muon* leg1 = dynamic_cast<const pat::Muon*>( (theDiTau->leg1()).get() );
  const pat::Tau*  leg2 = dynamic_cast<const pat::Tau*>(  (theDiTau->leg2()).get() );

  vector<string> triggerPaths;
  vector<string> XtriggerPaths;
  vector<string> HLTfiltersMu;
  vector<string> HLTfiltersTau;

  if(isMC_){

    // X-triggers 
    XtriggerPaths.push_back("HLT_Mu15_LooseIsoPFTau20_v*");
    XtriggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v*");
    XtriggerPaths.push_back("HLT_Mu15_v*");
    XtriggerPaths.push_back("HLT_IsoMu12_v*");

    // Single Mu triggers
    triggerPaths.push_back("HLT_Mu15_LooseIsoPFTau20_v2");
    triggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v2");
    triggerPaths.push_back("HLT_Mu15_v2");
    triggerPaths.push_back("HLT_IsoMu12_v1");

    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered12");
    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered15");

    HLTfiltersTau.push_back("hltOverlapFilterIsoMu12IsoPFTau10");
  }
  else{

    XtriggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v*");
    XtriggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v*");
    XtriggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau20_v*");
    XtriggerPaths.push_back("HLT_IsoMu12_v*");


    triggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v1");
    triggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v2");
    triggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v4");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v2");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v4");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v5");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v6");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v8");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v9");
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau20_v2");
    triggerPaths.push_back("HLT_IsoMu12_v1");
    triggerPaths.push_back("HLT_IsoMu12_v2");

    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered12");
    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered15");

    HLTfiltersTau.push_back("hltOverlapFilterIsoMu12IsoPFTau10");
    HLTfiltersTau.push_back("hltOverlapFilterIsoMu15IsoPFTau15");
    HLTfiltersTau.push_back("hltOverlapFilterIsoMu15IsoPFTau20");

  }

  for(unsigned int i=0;i<triggerPaths.size();i++){
    if(!trigger){
      continue;
      cout << "Invalid triggerEvent !" << endl;
    }
    const pat::TriggerPath *triggerPath =  trigger->path(triggerPaths[i]);

    if(verbose_){
      cout<<  "Testing " << triggerPaths[i] << endl;
      if(triggerPath) cout << "Is there..." << endl;
      if(triggerPath && triggerPath->wasRun()) cout << "Was run..." << endl;
      if(triggerPath && triggerPath->wasRun() && triggerPath->wasAccept()) cout << "Was accepted..." << endl;
    }

    if(triggerPath && triggerPath->wasRun() && 
       triggerPath->wasAccept() && 
       triggerPath->prescale()==1) triggerBits_->push_back(1);
    else if (triggerPath && triggerPath->wasRun() && 
	     triggerPath->wasAccept() && 
	     triggerPath->prescale()!=1) triggerBits_->push_back(2);
    else triggerBits_->push_back(0);
  }

  /*
  for(unsigned int m = 0; m<XtriggerPaths.size(); m++){
    if((leg1->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()!=0 && 
       (leg2->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()!=0) tauXTriggers_->push_back(1);
    else if((leg1->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()!=0 && 
	    (leg2->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()==0)  tauXTriggers_->push_back(2);
    else if((leg1->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()==0 && 
	    (leg2->triggerObjectMatchesByPath(XtriggerPaths[m],false,false)).size()!=0)  tauXTriggers_->push_back(3);
    else tauXTriggers_->push_back(0);
  }
  */

  for(unsigned int i=0 ; i< HLTfiltersMu.size() ; i++){
    bool matched = false;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
      pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
      if(verbose_) {
	if( Geom::deltaR( aObj->triggerObject().p4(), leg1->p4() )<0.3 ){
	  for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
	    cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of muon" << endl;
	  }
	}
      }
      if( Geom::deltaR( aObj->triggerObject().p4(), leg1->p4() )<0.3  && aObj->hasFilterLabel(HLTfiltersMu[i]) ){
	matched = true;
      }
    }
    if(matched) 
      tauXTriggers_->push_back(1);
    else 
      tauXTriggers_->push_back(0);
    if(verbose_){
      if(matched) cout << "Muon matched within dR=0.3 with trigger object passing filter " << HLTfiltersMu[i] << endl;
      else cout << "!!! Muon is not trigger matched within dR=0.3 !!!" << endl;
    }
  }
  for(unsigned int i=0 ; i< HLTfiltersTau.size() ; i++){
    bool matched = false;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
      pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
      if(verbose_) {
        if( Geom::deltaR( aObj->triggerObject().p4(), leg2->p4() )<0.3 ){
          for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
            cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of tau" << endl;
          }
        }
      }
      if( Geom::deltaR( aObj->triggerObject().p4(), leg2->p4() )<0.3  && aObj->hasFilterLabel(HLTfiltersTau[i]) ){
	matched = true;
      }
    }
    if(matched) 
      tauXTriggers_->push_back(1);
    else 
      tauXTriggers_->push_back(0);
    if(verbose_){
      if(matched) cout << "Tau matched within dR=0.3 with trigger object passing filter " << HLTfiltersTau[i] << endl;
      else cout << "!!! Tau is not trigger matched within dR=0.3 !!!" << endl;
    }
  }


  // triggers Mu
  if(verbose_){
    const pat::TriggerObjectStandAloneCollection trColl = leg1->triggerObjectMatchesByType(trigger::TriggerMuon);
    cout << "Mu triggers" << endl;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false,false)).size(); k++){
	cout << (it->pathNames(false,false))[k] << " with filters: " << endl;
	for(unsigned int m = 0; m < (it->filterLabels()).size(); m++){
	  cout << "  - " << (it->filterLabels())[m] << endl;
	}
      }
    }
  }
  // triggers Tau
  if(verbose_){
    const pat::TriggerObjectStandAloneCollection trColl = leg2->triggerObjectMatchesByType(trigger::TriggerTau);
    cout << "Tau triggers" << endl;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false,false)).size(); k++){
	cout << (it->pathNames(false,false))[k]  << " with filters: " << endl;
	for(unsigned int m = 0; m < (it->filterLabels()).size(); m++){
	  cout << "  - " << (it->filterLabels())[m] << endl;
	}
      }
    }
  }

  diTauLegsP4_->push_back(leg1->p4());
  diTauLegsP4_->push_back(leg2->p4());
  
  if(isMC_){
    if( (leg1->genParticleById(13,0,true)).isNonnull() ||  (leg1->genParticleById(-13,0,true)).isNonnull() ){
      genDiTauLegsP4_->push_back( leg1->genParticleById(13,0,true)->p4() );
      isMuLegMatched_ = 1;
    }
    else{
      genDiTauLegsP4_->push_back( math::XYZTLorentzVectorD(0,0,0,0) );
      if(verbose_){
	for(unsigned int l = 0; l < leg1->genParticlesSize() ; l++){
	  if((leg1->genParticleRefs())[l]->pt() < 0.5 ) continue;
	  cout << "Mu leg matchged to particle " << (leg1->genParticleRefs())[l]->pdgId() 
	       << " with pt " << (leg1->genParticleRefs())[l]->pt()
	       << endl;
	}
      }
    }

    bool leg2IsFromMu = false;
    math::XYZTLorentzVectorD genMuP4(0,0,0,0);
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( abs((*genParticles)[k].pdgId()) != 11  || (*genParticles)[k].status()!=3 )
        continue;
      if(Geom::deltaR( (*genParticles)[k].p4(),leg2->p4())<0.15){
        leg2IsFromMu = true;
        genMuP4 = (*genParticles)[k].p4();
      }
    }

    if( leg2->genJet() !=0 )
      genDiTauLegsP4_->push_back(leg2->genJet()->p4());
    else if(leg2IsFromMu)
      genDiTauLegsP4_->push_back( genMuP4 );
    else{
      genDiTauLegsP4_->push_back( math::XYZTLorentzVectorD(0,0,0,0) );
      if(verbose_) cout << "WARNING: no genJet matched to the leg2 with eta,phi " << leg2->eta() << ", " << leg2->phi() << endl;
    }

    bool tauHadMatched = false;
    for(unsigned int k = 0; k < tauGenJets->size(); k++){
      if( Geom::deltaR( (*tauGenJets)[k].p4(),leg2->p4() ) < 0.15 ) tauHadMatched = true;
    }

    if( ( (leg2->genParticleById(15,0,true)).isNonnull() || (leg2->genParticleById(-15,0,true)).isNonnull() ) && tauHadMatched ) isTauLegMatched_ = 1;
    else if(verbose_){
      for(unsigned int l = 0; l < leg2->genParticlesSize() ; l++){
	if((leg2->genParticleRefs())[l]->pt() < 0.5 ) continue;
	cout << "Tau leg matchged to particle " << (leg2->genParticleRefs())[l]->pdgId() 
	     << " with pt " << (leg2->genParticleRefs())[l]->pt()
	     << endl;
      }
    }
  }

  if((leg2->signalPFChargedHadrCands()).size()==1 && (leg2->signalPFGammaCands()).size()==0) decayMode_ = 0; 
  else if((leg2->signalPFChargedHadrCands()).size()==1 && (leg2->signalPFGammaCands()).size()>0)  decayMode_ = 1; 
  else if((leg2->signalPFChargedHadrCands()).size()==3) decayMode_ = 2; 
  else  decayMode_ = -99;

  for(unsigned int k = 0 ; k < (leg2->signalPFGammaCands()).size() ; k++){
    reco::PFCandidateRef gamma = (leg2->signalPFGammaCands()).at(k);
    if( (leg2->leadPFChargedHadrCand()).isNonnull() )
      gammadR_->push_back( Geom::deltaR( gamma->p4(), leg2->leadPFChargedHadrCand()->p4() ) );
    else
      gammadR_->push_back( Geom::deltaR( gamma->p4(), leg2->p4() ) );
    gammaPt_->push_back(  gamma->pt() );
  }
  
  visibleTauMass_ = leg2->mass();
  if((leg2->leadPFChargedHadrCand()->trackRef()).isNonnull()){
    leadPFChargedHadrCandTrackPt_ = leg2->leadPFChargedHadrCand()->trackRef()->pt();
  } else if((leg2->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull()){
    leadPFChargedHadrCandTrackPt_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->pt();
  } else leadPFChargedHadrCandTrackPt_ = -99;

  tightestHPSWP_ = -1;
  if(leg2->tauID("byVLooseIsolation")>0.5) tightestHPSWP_++;
  if(leg2->tauID("byLooseIsolation")>0.5)  tightestHPSWP_++;
  if(leg2->tauID("byMediumIsolation")>0.5) tightestHPSWP_++;
  if(leg2->tauID("byTightIsolation")>0.5)  tightestHPSWP_++;
  tightestHPSDBWP_ = -1;
  if(leg2->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;
  if(leg2->tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;


  dxy1_ = vertexes->size()!=0 ? leg1->globalTrack()->dxy( (*vertexes)[0].position() ) : -99;
  dz1_ = vertexes->size()!=0 ? leg1->globalTrack()->dz( (*vertexes)[0].position() ) : -99;
  dxy2_ = -99;
  dz2_ = -99;
  hasKft_ = 0;
  if( vertexes->size()!=0 && (leg2->leadPFChargedHadrCand()).isNonnull() ){
    if( (leg2->leadPFChargedHadrCand()->trackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->trackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->trackRef()->dz( (*vertexes)[0].position() );
      hasKft_ = 1;
    }
    if( (leg2->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dz( (*vertexes)[0].position() );
    }
  }



  // isoDeposit definition: 2010
  isodeposit::AbsVetos vetos2010ChargedLeg1; 
  isodeposit::AbsVetos vetos2010NeutralLeg1; 
  isodeposit::AbsVetos vetos2010PhotonLeg1;
  // isoDeposit definition: 2011
  isodeposit::AbsVetos vetos2011ChargedLeg1; 
  isodeposit::AbsVetos vetos2011NeutralLeg1; 
  isodeposit::AbsVetos vetos2011PhotonLeg1;
 
  vetos2010ChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2010NeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.08));
  vetos2010NeutralLeg1.push_back(new isodeposit::ThresholdVeto(1.0));
  vetos2010PhotonLeg1.push_back( new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.05));
  vetos2010PhotonLeg1.push_back( new isodeposit::ThresholdVeto(1.0));

  vetos2011ChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.0));
  vetos2011NeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011NeutralLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2011PhotonLeg1.push_back( new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011PhotonLeg1.push_back( new isodeposit::ThresholdVeto(0.5));

  chIsoLeg1v1_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  nhIsoLeg1v1_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetos2010NeutralLeg1).first;
  phIsoLeg1v1_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetos2010PhotonLeg1).first;
  chIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010ChargedLeg1).first;
  nhIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010NeutralLeg1).first;
  phIsoPULeg1v1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2010PhotonLeg1).first;

  chIsoLeg1v2_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoLeg1v2_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoLeg1v2_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;
  chIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;


  chIsoLeg2_ = leg2->isolationPFChargedHadrCandsPtSum();
  //
  if(verbose_){
    cout << "Tau z position " << (leg2->vertex()).z() << endl;
    for(unsigned int i = 0; i < (leg2->isolationPFChargedHadrCands()).size(); i++){
      if( (((leg2->isolationPFChargedHadrCands()).at(i))->trackRef()).isNonnull() )
	cout << "Ch # " << i << " has z position " << (((leg2->isolationPFChargedHadrCands()).at(i))->trackRef()->referencePoint()).z() << " and pt " << ((leg2->isolationPFChargedHadrCands()).at(i))->pt() << endl;
    }
  }
  //
  nhIsoLeg2_ = 0.;
  for(unsigned int i = 0; i < (leg2->isolationPFNeutrHadrCands()).size(); i++){
    nhIsoLeg2_ += (leg2->isolationPFNeutrHadrCands()).at(i)->pt();
  }
  phIsoLeg2_ = leg2->isolationPFGammaCandsEtSum();
   

  // cleaning
  for(unsigned int i = 0; i <vetos2010ChargedLeg1.size(); i++){
    delete vetos2010ChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2010NeutralLeg1.size(); i++){
    delete vetos2010NeutralLeg1[i];
    delete vetos2010PhotonLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011ChargedLeg1.size(); i++){
    delete vetos2011ChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011NeutralLeg1.size(); i++){
    delete vetos2011NeutralLeg1[i];
    delete vetos2011PhotonLeg1[i];
  }
  //

  diTauVisP4_->push_back( theDiTau->p4Vis() );
  diTauCAP4_->push_back(  theDiTau->p4CollinearApprox() );
  diTauICAP4_->push_back( theDiTau->p4ImprovedCollinearApprox() );

  if(verbose_){
    if(theDiTau->collinearApproxIsValid() && theDiTau->ImprovedCollinearApproxIsValid())
      cout << "Found CA and ICA valid solutions" << endl;
    else if(!theDiTau->collinearApproxIsValid() && theDiTau->ImprovedCollinearApproxIsValid())
      cout << "ICA rescued the event!!!" << endl;
    else if(theDiTau->collinearApproxIsValid() && !theDiTau->ImprovedCollinearApproxIsValid())
      cout << "Strange!!! x1=" << theDiTau->x1CollinearApprox() << ", x2=" << theDiTau->x2CollinearApprox() << endl;
    if((theDiTau->p4CollinearApprox()).M()>0 && (theDiTau->p4ImprovedCollinearApprox()).M()<=0)
      cout << "Watch out! ICA gave mass=0 and CA mass !=0 ..." << endl;
  }

  math::XYZTLorentzVectorD nSVfitFitP4(0,0,0,(theDiTau->p4Vis()).M() );
  if( theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)!=0 /* &&theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->isValidSolution()*/ ){
    if(verbose_) cout << "Visible mass ==> " << nSVfitFitP4.E() << endl;
    nSVfitFitP4.SetPxPyPzE( 0,0,0, theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->mass() ) ;
    if(verbose_) cout << "SVFit fit solution ==> " << nSVfitFitP4.E() << endl;
  }
  diTauSVfitP4_->push_back( nSVfitFitP4  );

  diTauNSVfitMass_        = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->mass()        : -99; 
  diTauNSVfitMassErrUp_   = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->massErrUp()   : -99; 
  diTauNSVfitMassErrDown_ = (theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->isValidSolution() ) ? theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->massErrDown() : -99; 

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  lumi_  = iEvent.luminosityBlock();

  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJetsIDL1Offset;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJetsIDUp;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJetsIDDown;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedGenJetsID;
  std::map<double, std::pair<float,float> ,  MuTauStreamAnalyzer::more> bTaggers;
  std::map<double, std::pair<float,float> ,  MuTauStreamAnalyzer::more> jetPVassociation;
  std::map<double, std::pair<float,float> ,  MuTauStreamAnalyzer::more> jetMoments;
  


  for(unsigned int it = 0; it < jets->size() ; it++){

    math::XYZTLorentzVectorD leg2p4 = ( (leg2->pfJetRef()).isNonnull() ) ? leg2->pfJetRef()->p4() : leg2->p4();

    if( Geom::deltaR((*jets)[it].p4(),leg1->p4())<deltaRLegJet_ || 
	Geom::deltaR((*jets)[it].p4(), leg2p4 )<deltaRLegJet_ ){
      if(verbose_) cout << "The jet at (" <<(*jets)[it].pt()<<","<<(*jets)[it].eta()<<") is closer than "<<deltaRLegJet_ << " from one of the legs" << endl;  
      continue;
    }

    /////////////////////////////////////////////////////////////////////////
    //// try to use JES uncertainties
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameters;
    // get the jet corrector parameters collection from the global tag
    iSetup.get<JetCorrectionsRecord>().get("AK5PF", jetCorrParameters);
    // get the uncertainty parameters from the collection
    JetCorrectorParameters const & param = (*jetCorrParameters)["Uncertainty"];
    // instantiate the jec uncertainty object
    JetCorrectionUncertainty* deltaJEC = new JetCorrectionUncertainty(param);
    deltaJEC->setJetEta((*jets)[it].eta()); deltaJEC->setJetPt((*jets)[it].pt());
    float shift  = deltaJEC->getUncertainty( true );
    /////////////////////////////////////////////////////////////////////////

    if(verbose_){
      cout << "Analyzing jet with pt " << ((*jets)[it].p4()).Pt() 
	   << "and eta " << ((*jets)[it].p4()).Eta() << endl;
      cout << " ==> JEC uncertainty is " << shift*100 << " %" << endl;
      pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);
      for(unsigned int i = 0; i < (jet->availableJECSets()).size() ; i++ ){
	std::cout << (jet->availableJECSets())[i] << std::endl;
      }
      for(unsigned int i = 0; i < (jet->availableJECLevels()).size() ; i++ ){
	std::cout << (jet->availableJECLevels())[i] << std::endl;
      }
      std::cout << "L1FastJet ========> " << std::endl;
      std::cout << "Uncorrected " << jet->correctedJet("Uncorrected","none","patJetCorrFactors").pt() << std::endl;
      std::cout << "L1FastJet "   << jet->correctedJet("L1FastJet","none",  "patJetCorrFactors").pt() << std::endl;
      std::cout << "L2Relative "  << jet->correctedJet("L2Relative","none", "patJetCorrFactors").pt() << std::endl; 
      std::cout << "L3Absolute "  << jet->correctedJet("L3Absolute","none", "patJetCorrFactors").pt() << std::endl; 
      std::cout << "L1Offset ========> " << std::endl;
      std::cout << "Uncorrected " << jet->jecFactor("Uncorrected","none","patJetCorrFactorsL1Offset")*jet->pt() << std::endl;
      std::cout << "L1Offset"     << jet->jecFactor("L1Offset","none",   "patJetCorrFactorsL1Offset")*jet->pt() << std::endl;
      std::cout << "L2Relative "  << jet->jecFactor("L2Relative","none", "patJetCorrFactorsL1Offset")*jet->pt() << std::endl;
      std::cout << "L3Absolute"   << jet->jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*jet->pt() << std::endl;  
    }

    std::map<string,float> aMap;

    if( jetID( &(*jets)[it] , &((*vertexes)[0]), vtxZ, aMap ) < minJetID_ )  continue;

    sortedJets.insert( make_pair( (*jets)[it].correctedJet("Uncorrected").p4().Pt() ,(*jets)[it].correctedJet("Uncorrected").p4() ) );

    if((*jets)[it].p4().Pt() < minCorrPt_) continue;

    //add b-tag info
    bTaggers.insert( make_pair((*jets)[it].p4().Pt(), make_pair( (*jets)[it].bDiscriminator("trackCountingHighEffBJetTags"),(*jets)[it].bDiscriminator("trackCountingHighPurBJetTags")  ) ) );
    jetPVassociation.insert( make_pair( (*jets)[it].p4().Pt(), make_pair(aMap["chFracRawJetE"],aMap["chFracAllChargE"]) ) );
    jetMoments.insert( make_pair( (*jets)[it].p4().Pt(), make_pair( (*jets)[it].etaetaMoment(),(*jets)[it].phiphiMoment()) ) );
    
    //jetsBtagHE_->push_back((*jets)[it].bDiscriminator("trackCountingHighEffBJetTags"));
    //jetsBtagHP_->push_back((*jets)[it].bDiscriminator("trackCountingHighPurBJetTags"));

    if(isMC_) sortedJetsIDL1Offset.insert( make_pair( (*jets)[it].jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*(*jets)[it].pt() , (*jets)[it].jecFactor("L3Absolute","none", "patJetCorrFactorsL1Offset")*(*jets)[it].p4()) );   
    else sortedJetsIDL1Offset.insert( make_pair( (*jets)[it].jecFactor("L2L3Residual","none", "patJetCorrFactorsL1Offset")*(*jets)[it].pt() , (*jets)[it].jecFactor("L2L3Residual","none", "patJetCorrFactorsL1Offset")*(*jets)[it].p4()) ); 

    if(verbose_) cout << "Components: "
		      << "px=" << ((*jets)[it].p4()).Px() << " (" << (*jets)[it].px() << "), "
		      << "py=" << ((*jets)[it].p4()).Py() << " (" << (*jets)[it].py() << "), "
		      << "pz=" << ((*jets)[it].p4()).Pz() << " (" << (*jets)[it].pz() << "), "
		      << "E="  << ((*jets)[it].p4()).E()  << " (" << (*jets)[it].energy()  << ")"
		      << endl;
    sortedJetsID.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );
    sortedJetsIDUp.insert( make_pair( (*jets)[it].p4().Pt() ,  (*jets)[it].p4()*(1+shift) ) );
    sortedJetsIDDown.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4()*(1-shift) ) );

    if(isMC_){
      if((*jets)[it].genJet() != 0) sortedGenJetsID.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].genJet()->p4() ) );
      else sortedGenJetsID.insert( make_pair( (*jets)[it].p4().Pt() , math::XYZTLorentzVectorD(0,0,0,0) ) );
    }
     
  }
  
  for(CImap it = sortedJets.begin(); it != sortedJets.end() ; it++){
    jetsP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    jetsIDP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDUp.begin(); it != sortedJetsIDUp.end() ; it++){
    jetsIDUpP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDDown.begin(); it != sortedJetsIDDown.end() ; it++){
    jetsIDDownP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsIDL1Offset.begin(); it != sortedJetsIDL1Offset.end() ; it++){
    jetsIDL1OffsetP4_->push_back( it->second );
  }
  for(CImap it = sortedGenJetsID.begin(); it != sortedGenJetsID.end() ; it++){
    genJetsIDP4_->push_back( it->second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = bTaggers.begin(); it != bTaggers.end() ; it++){
    jetsBtagHE_->push_back( (it->second).first  );
    jetsBtagHP_->push_back( (it->second).second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = jetPVassociation.begin(); it != jetPVassociation.end() ; it++){
    jetsChEfraction_->push_back( (it->second).first  );
    jetsChNfraction_->push_back( (it->second).second );
  }
  for(std::map<double, std::pair<float,float> >::iterator it = jetMoments.begin(); it != jetMoments.end() ; it++){
    jetMoments_->push_back( (it->second).first  );
    jetMoments_->push_back( (it->second).second );
  }

  tree_->Fill();

}


unsigned int  MuTauStreamAnalyzer::jetID( const pat::Jet* jet, const reco::Vertex* vtx, std::vector<float> vtxZ, std::map<std::string,float>& map_){

  if( (jet->pt())<10 ) return 99; // always pass jet ID

  std::vector<reco::PFCandidatePtr> pfCandPtrs = jet->getPFConstituents();

  int nCharged = 0;
  int nPhotons = 0;
  int nNeutral = 0;
  int nConst = 0;

  float energyCharged = 0;
  float energyPhotons = 0;
  float energyNeutral = 0;
  float energyElectrons = 0;
 
  float totalEnergyFromConst = 0;

  float chEnFractionPV   = 0;
  float chEnFractionChPV = 0;
  float chFractionPV     = 0;

  for(unsigned i=0; i<pfCandPtrs.size(); ++i) {
    const reco::PFCandidate& cand = *(pfCandPtrs[i]);

    totalEnergyFromConst +=  cand.energy();
    nConst += 1;

    switch( cand.particleId() ) {
    case reco::PFCandidate::h: 
      nCharged++;
      energyCharged += cand.energy(); 

      if((cand.trackRef()).isNonnull()){
	bool isMatched = false;
	for(reco::Vertex::trackRef_iterator it = vtx->tracks_begin() ; it!=vtx->tracks_end() && !isMatched; it++){
	  if(  (*it).id() == (cand.trackRef()).id() && (*it).key() == (cand.trackRef()).key() ){
	    isMatched = true;
	    if(verbose_) cout << (*it).id() << ", " << (*it).key() << " is matched!" << endl;
	    chEnFractionPV += cand.energy();
	    chFractionPV+=1;
	  }
	}
	if(!isMatched){
	  float dist = vtxZ.size()>0 ? fabs(vtxZ[0]-((cand.trackRef())->referencePoint()).z()) : 999.;
	  float tmpDist = 999.;
	  for(unsigned k = 1; k < vtxZ.size(); k++){
	    if( fabs(((cand.trackRef())->referencePoint()).z() - vtxZ[k] ) < tmpDist )
	      tmpDist = fabs(((cand.trackRef())->referencePoint()).z() - vtxZ[k] );
	  }
	  if(tmpDist>dist){
	    isMatched = true;
	    if(verbose_) cout << "Matched by closest vtx in z!!" << endl;
	    chEnFractionPV += cand.energy();
	    chFractionPV+=1;
	  }
	}
	if(!isMatched && verbose_) {
	  cout << "Ch with pt " << cand.pt() << " and eta " << cand.eta() << " is not matched to the PV !!!" << endl;
	  cout << "z position of PV " << (vtx->position()).z() << ", z position of track " << ((cand.trackRef())->referencePoint()).z()  << ", x position of track " << ((cand.trackRef())->referencePoint()).x()   << ", y position of track " << ((cand.trackRef())->referencePoint()).y() << endl;
	}
      }

      break;
    case reco::PFCandidate::gamma:
      nPhotons++;
      energyPhotons += cand.energy();
      break;
    case reco::PFCandidate::h0:
      nNeutral++;
      energyNeutral += cand.energy();
      break;
    case reco::PFCandidate::e: 
      energyElectrons += cand.energy();

      if((cand.gsfTrackRef()).isNonnull()){
	bool isMatched = false;
	float dist = vtxZ.size()>0 ? fabs(vtxZ[0]-((cand.gsfTrackRef())->referencePoint()).z()) : 999.;
	float tmpDist = 999.;
	for(unsigned k = 1; k < vtxZ.size(); k++){
	  if( fabs(((cand.gsfTrackRef())->referencePoint()).z() - vtxZ[k] ) < tmpDist )
	    tmpDist = fabs(((cand.gsfTrackRef())->referencePoint()).z() - vtxZ[k] );
	}
	if(tmpDist>dist){
	  isMatched = true;
	  if(verbose_) cout << "Matched by closest vtx in z!!" << endl;
	  chEnFractionPV += cand.energy();
	  chFractionPV+=1;
	}
	if(!isMatched && verbose_) {
	  cout << "Ele with pt " << cand.pt() << " and eta " << cand.eta() << " is not matched to the PV !!!" << endl;
	  cout << "z position of PV " << (vtx->position()).z() << ", z position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).z()  << ", x position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).x()   << ", y position of gsfTrack " << ((cand.gsfTrackRef())->referencePoint()).y() << endl;
	}
      }
      
      break;
    case reco::PFCandidate::h_HF: // fill neutral
      nNeutral++;
      energyNeutral += cand.energy();
      break;
    case reco::PFCandidate::egamma_HF: // fill e/gamma
      nPhotons++;
      energyPhotons += cand.energy();
      break;
    default:
      break;
    }
  }

  chEnFractionChPV = chEnFractionPV;
  if(energyCharged>0 || energyElectrons>0)
    chEnFractionChPV /= (energyCharged+energyElectrons); 
  chEnFractionPV /= jet->correctedJet("Uncorrected").p4().E();
  chFractionPV   /= nCharged;

  //if(chFractionPV==0) cout << energyCharged << endl;

  map_["chMult"]          =  chFractionPV;
  map_["chFracRawJetE"]   =  chEnFractionPV;
  map_["chFracAllChargE"] =  chEnFractionChPV;
        
  bool loose=false;
  bool medium=false;
  bool tight=false;

  //loose id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.99 && 
       energyPhotons/totalEnergyFromConst<0.99 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<0.99 && 
       energyPhotons/totalEnergyFromConst<0.99 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<0.99
       )
      ) loose = true;
  // medium id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.95 && 
       energyPhotons/totalEnergyFromConst<0.95 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<1
       )
      ) medium = true;
  // tight id
  if( (TMath::Abs(jet->eta())>2.4 && 
       energyNeutral/totalEnergyFromConst<0.90 && 
       energyPhotons/totalEnergyFromConst<0.90 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<1
       )
      ) tight = true;
  
  if(loose && !medium && !tight) return 1;
  if(loose && medium && !tight)  return 2;
  if(loose && medium && tight)   return 3; 
  
  return 0;

}



void MuTauStreamAnalyzer::endJob(){}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(MuTauStreamAnalyzer);


