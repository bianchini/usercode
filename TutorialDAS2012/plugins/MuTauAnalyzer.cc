#include "Bianchi/TutorialDAS2012/interface/MuTauAnalyzer.h"
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

typedef std::map<double, math::XYZTLorentzVectorD ,MuTauAnalyzer::more>::iterator CImap;

MuTauAnalyzer::MuTauAnalyzer(const edm::ParameterSet & iConfig){

  diTauTag_          = iConfig.getParameter<edm::InputTag>("diTaus");
  jetsTag_           = iConfig.getParameter<edm::InputTag>("jets");
  rawMetTag_         = iConfig.getParameter<edm::InputTag>("rawMet");
  muonsTag_          = iConfig.getParameter<edm::InputTag>("muons");
  muonsRelTag_       = iConfig.getParameter<edm::InputTag>("muonsRel");
  verticesTag_       = iConfig.getParameter<edm::InputTag>("vertices");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults"); 
  isMC_              = iConfig.getParameter<bool>("isMC");
  deltaRLegJet_      = iConfig.getUntrackedParameter<double>("deltaRLegJet",0.3);
  minCorrPt_         = iConfig.getUntrackedParameter<double>("minCorrPt",10.);
  minJetID_          = iConfig.getUntrackedParameter<double>("minJetID",0.5);
  verbose_           = iConfig.getUntrackedParameter<bool>("verbose",false);
}

void MuTauAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","qqH tree");
 
  tauXTriggers_= new std::vector< int >();
  triggerBits_ = new std::vector< int >();

  jetsIDP4_       = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genJetsIDP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauVisP4_     = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfitP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauLegsP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genDiTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  METP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genMETP4_       = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genVP4_         = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("genJetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genJetsIDP4_);
  
  tree_->Branch("tauXTriggers","std::vector<int>",&tauXTriggers_);
  tree_->Branch("triggerBits","std::vector<int>",&triggerBits_);

  tree_->Branch("diTauVisP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",  &diTauVisP4_);
  tree_->Branch("diTauSVfitP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfitP4_);

  tree_->Branch("diTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsP4_);
  tree_->Branch("genDiTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genDiTauLegsP4_);

  tree_->Branch("METP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&METP4_);
  tree_->Branch("genMETP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genMETP4_);
  tree_->Branch("genVP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genVP4_);
  tree_->Branch("genDecay",&genDecay_,"genDecay/I");

  tree_->Branch("sumEt", &sumEt_,"sumEt/F");
  tree_->Branch("MtLeg1",&MtLeg1_,"MtLeg1/F");

  tree_->Branch("chIsoLeg1v2",  &chIsoLeg1v2_,"chIsoLeg1v2/F");
  tree_->Branch("nhIsoLeg1v2",  &nhIsoLeg1v2_,"nhIsoLeg1v2/F");
  tree_->Branch("phIsoLeg1v2",  &phIsoLeg1v2_,"phIsoLeg1v2/F");
  tree_->Branch("elecIsoLeg1v2",&elecIsoLeg1v2_,"elecIsoLeg1v2/F");
  tree_->Branch("muIsoLeg1v2",  &muIsoLeg1v2_ ,"muIsoLeg1v2/F");
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
  tree_->Branch("numOfLooseIsoDiTaus",&numOfLooseIsoDiTaus_,"numOfLooseIsoDiTaus/I");
  tree_->Branch("decayMode",&decayMode_,"decayMode/I");

  tree_->Branch("tightestHPSDBWP",&tightestHPSDBWP_,"tightestHPSDBWP/I");
  tree_->Branch("visibleTauMass",&visibleTauMass_,"visibleTauMass/F");

  tree_->Branch("signalPFChargedHadrCands",&signalPFChargedHadrCands_,"signalPFChargedHadrCands/I");
  tree_->Branch("signalPFGammaCands",&signalPFGammaCands_,"signalPFGammaCands/I");

  tree_->Branch("isTauLegMatched",&isTauLegMatched_,"isTauLegMatched/I");
  tree_->Branch("isMuLegMatched",&isMuLegMatched_,"isMuLegMatched/I");
  tree_->Branch("muFlag",&muFlag_,"muFlag/I");

  tree_->Branch("diTauCharge",&diTauCharge_,"diTauCharge/F");
  tree_->Branch("rhoFastJet",&rhoFastJet_,"rhoFastJet/F");

  tree_->Branch("nPUVertices",&nPUVertices_,"nPUVertices/I");
  tree_->Branch("nPUVerticesM1",&nPUVerticesM1_,"nPUVerticesM1/I");
  tree_->Branch("nPUVerticesP1",&nPUVerticesP1_,"nPUVerticesP1/I");



}


MuTauAnalyzer::~MuTauAnalyzer(){
  delete jetsIDP4_; 
  delete METP4_; delete diTauVisP4_; 
  delete diTauSVfitP4_; delete genVP4_;
  delete diTauLegsP4_; delete tauXTriggers_; delete triggerBits_;
  delete genJetsIDP4_; delete genDiTauLegsP4_; delete genMETP4_;
}

void MuTauAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  jetsIDP4_->clear();
  diTauVisP4_->clear();
  diTauSVfitP4_->clear();
  diTauLegsP4_->clear();
  METP4_->clear();
  genVP4_->clear();
  genJetsIDP4_->clear();
  genDiTauLegsP4_->clear();
  genMETP4_->clear();
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
  iEvent.getByLabel( verticesTag_ ,pvHandle);
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
  iEvent.getByLabel( rawMetTag_, metHandle);
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
  nPUVertices_       = -99;
  nPUVerticesP1_     = -99;
  nPUVerticesM1_     = -99;

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
	if(it->getBunchCrossing() == 0 ) 
	  nPUVertices_  = it->getPU_NumInteractions();
	else if(it->getBunchCrossing() == -1)  
	  nPUVerticesM1_= it->getPU_NumInteractions();
	else if(it->getBunchCrossing() == +1)  
	  nPUVerticesP1_= it->getPU_NumInteractions();
      }
    }
  }

  if(verbose_){
    cout << "Num of PU = "          << nPUVertices_ << endl;
    cout << "Num of OOT PU = "      << nPUVerticesM1_+nPUVerticesP1_<< endl;
  }

  edm::Handle<double> rhoFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentral","rho", ""), rhoFastJetHandle);
  if( !rhoFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho label available \n";
  rhoFastJet_ = (*rhoFastJetHandle);

  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel( muonsTag_ ,muonsHandle);
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
  iEvent.getByLabel(muonsRelTag_, muonsRelHandle);
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
  numOfLooseIsoDiTaus_ = diTaus->size();

  unsigned int muIndex = 0;
  muFlag_ = 0;

  bool found = false;
  for(unsigned int i=0; i<muonsRel->size(); i++){
    for(unsigned int j=0; j<muons->size(); j++){
      
      if(found) continue;
      
      if( Geom::deltaR( (*muonsRel)[i].p4(),(*muons)[j].p4())>0.3
	  && (*muonsRel)[i].charge()*(*muons)[j].charge()<0
	  && (*muons)[j].userFloat("PFRelIsoDB04v2")<0.3 && (*muonsRel)[i].userFloat("PFRelIsoDB04v2")<0.3 ){
	muFlag_ = 1;
	if(verbose_) cout<< "Two muons failing diMu veto: flag= " << muFlag_ << endl;
	found=true;
      }
      else if( Geom::deltaR( (*muonsRel)[i].p4(),(*muons)[j].p4())>0.3
	       && (*muonsRel)[i].charge()*(*muons)[j].charge()>0
	       && (*muons)[j].userFloat("PFRelIsoDB04v2")<0.3 && (*muonsRel)[i].userFloat("PFRelIsoDB04v2")<0.3 ){
	muFlag_ = 2;
	if(verbose_) cout<< "Two muons with SS: flag= " << muFlag_ << endl;
	found=true;
      }
    }
  }
  

  // select highest-pt muon that forms at least one di-tau
  // no ambiguity if "muons" are isolated ==> if extra
  // isolated muon exists the event is vetoed
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

  std::vector<unsigned int> selectedDiTausFromMu;
  for(unsigned int i=0; i< diTaus->size(); i++){
    if( Geom::deltaR(((*diTaus)[i].leg1())->p4(),((*muons)[muIndex]).p4())<0.01 ) 
      selectedDiTausFromMu.push_back(i);
  }
  if(verbose_) cout<< "Selection on Muons has selected " << selectedDiTausFromMu.size() << " diTaus...check tau leg now..." << endl;
  
  
  unsigned int index = 0;
  float maxTauPt     = 0;

  for(unsigned int i=0; i<selectedDiTausFromMu.size(); i++){
    if( ((*diTaus)[selectedDiTausFromMu[i]].leg2())->pt() > maxTauPt )
      index = selectedDiTausFromMu[i];
  }
  if(verbose_) cout << "Chosen index " << index << endl;

  theDiTau = &(*diTaus)[index];
  
  diTauCharge_ = theDiTau->charge();
  METP4_->push_back((*met)[0].p4());
  if(isMC_) genMETP4_->push_back( (*met)[0].genMET()->p4() );
  sumEt_  = (*met)[0].sumEt();
  MtLeg1_ =  theDiTau->mt1MET();

  isMuLegMatched_  = 0;
  isTauLegMatched_ = 0;

  const pat::Muon* leg1 = dynamic_cast<const pat::Muon*>( (theDiTau->leg1()).get() );
  const pat::Tau*  leg2 = dynamic_cast<const pat::Tau*>(  (theDiTau->leg2()).get() );

  diTauLegsP4_->push_back(leg1->p4());
  diTauLegsP4_->push_back(leg2->p4());

  vector<string> triggerPaths;
  vector<string> HLTfiltersMu;
  vector<string> HLTfiltersTau;

  if(isMC_){

    // for Fall11
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v9");
    triggerPaths.push_back("HLT_IsoMu15_v14");
    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered15");
    HLTfiltersTau.push_back("hltOverlapFilterIsoMu15IsoPFTau15");
  }
  else{
    triggerPaths.push_back("HLT_IsoMu15_LooseIsoPFTau15_v8");    
    HLTfiltersMu.push_back("hltSingleMuIsoL3IsoFiltered15");
    HLTfiltersTau.push_back("hltOverlapFilterIsoMu15IsoPFTau15");
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

  if(isMC_){
    if( (leg1->genParticleById(13,0,true)).isNonnull() ){
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
      if( abs((*genParticles)[k].pdgId()) != 13  || (*genParticles)[k].status()!=3 )
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

  
  visibleTauMass_ = leg2->mass();

   
  signalPFChargedHadrCands_ = (leg2->signalPFChargedHadrCands()).size();
  signalPFGammaCands_       = (leg2->signalPFGammaCands()).size();
   
  tightestHPSDBWP_ = -1;
  if(leg2->tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;
  if(leg2->tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP_++;
  if(leg2->tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP_++;


  dxy1_ = vertexes->size()!=0 ? leg1->globalTrack()->dxy( (*vertexes)[0].position() ) : -99;
  dz1_ = vertexes->size()!=0 ? leg1->globalTrack()->dz( (*vertexes)[0].position() ) : -99;
  dxy2_ = -99;
  dz2_ = -99;

  if( vertexes->size()!=0 && (leg2->leadPFChargedHadrCand()).isNonnull() ){
    if( (leg2->leadPFChargedHadrCand()->trackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->trackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->trackRef()->dz( (*vertexes)[0].position() );
    }
    if( (leg2->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull() ){
      dxy2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dxy( (*vertexes)[0].position() );
      dz2_ = leg2->leadPFChargedHadrCand()->gsfTrackRef()->dz( (*vertexes)[0].position() );
    }
  }

  // isoDeposit definition: 2011
  isodeposit::AbsVetos vetos2011ChargedLeg1; 
  isodeposit::AbsVetos vetos2011NeutralLeg1; 
  isodeposit::AbsVetos vetos2011PhotonLeg1;
 
  vetos2011ChargedLeg1.push_back(new isodeposit::ConeVeto(reco::isodeposit::Direction(leg1->eta(),leg1->phi()),0.0001));
  vetos2011ChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.0));
  vetos2011NeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011NeutralLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetos2011PhotonLeg1.push_back( new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.01));
  vetos2011PhotonLeg1.push_back( new isodeposit::ThresholdVeto(0.5));

  chIsoLeg1v2_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoLeg1v2_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoLeg1v2_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;
  elecIsoLeg1v2_ = 
   leg1->isoDeposit(pat::User3Iso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  muIsoLeg1v2_   = 
    leg1->isoDeposit(pat::User2Iso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  chIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011ChargedLeg1).first;
  nhIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011NeutralLeg1).first;
  phIsoPULeg1v2_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetos2011PhotonLeg1).first;

  chIsoLeg2_ = leg2->isolationPFChargedHadrCandsPtSum();
  nhIsoLeg2_ = 0.;
  for(unsigned int i = 0; i < (leg2->isolationPFNeutrHadrCands()).size(); i++){
    nhIsoLeg2_ += (leg2->isolationPFNeutrHadrCands()).at(i)->pt();
  }
  phIsoLeg2_ = leg2->isolationPFGammaCandsEtSum();
   
  for(unsigned int i = 0; i <vetos2011ChargedLeg1.size(); i++){
    delete vetos2011ChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetos2011NeutralLeg1.size(); i++){
    delete vetos2011NeutralLeg1[i];
    delete vetos2011PhotonLeg1[i];
  }
  //

  diTauVisP4_->push_back( theDiTau->p4Vis() );

  math::XYZTLorentzVectorD nSVfitFitP4(0,0,0,(theDiTau->p4Vis()).M() );
  if( theDiTau->hasNSVFitSolutions() && theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)!=0 /*&& theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->isValidSolution()*/ ){
    if(verbose_) cout << "Visible mass ==> " << nSVfitFitP4.E() << endl;
    nSVfitFitP4.SetPxPyPzE( 0,0,0, theDiTau->nSVfitSolution("psKine_MEt_logM_fit",0)->mass() ) ;
    if(verbose_) cout << "SVFit fit solution ==> " << nSVfitFitP4.E() << endl;
  }
  diTauSVfitP4_->push_back( nSVfitFitP4  );

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  lumi_  = iEvent.luminosityBlock();

  std::map<double, math::XYZTLorentzVectorD ,MuTauAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,MuTauAnalyzer::more> sortedGenJetsID;
 

  for(unsigned int it = 0; it < jets->size() ; it++){

    pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);

    math::XYZTLorentzVectorD leg2p4 = ( (leg2->pfJetRef()).isNonnull() ) ? leg2->pfJetRef()->p4() : leg2->p4();

    if( Geom::deltaR(jet->p4(), leg1->p4())<deltaRLegJet_ || 
	Geom::deltaR(jet->p4(), leg2p4 )<deltaRLegJet_ ){
      if(verbose_) cout << "The jet at (" <<jet->pt()<<","<<jet->eta()<<") is closer than "<<deltaRLegJet_ << " from one of the legs" << endl;  
      continue;
    }

    if(verbose_){
      cout << "Analyzing jet with pt " << (jet->p4()).Pt() 
	   << "and eta " << (jet->p4()).Eta() << endl;
      std::cout << "L1FastJet ========> " << std::endl;
      std::cout << "Uncorrected " << jet->correctedJet("Uncorrected","none","patJetCorrFactors").pt() << std::endl;
      std::cout << "L1FastJet "   << jet->correctedJet("L1FastJet","none",  "patJetCorrFactors").pt() << std::endl;
      std::cout << "L2Relative "  << jet->correctedJet("L2Relative","none", "patJetCorrFactors").pt() << std::endl; 
      std::cout << "L3Absolute "  << jet->correctedJet("L3Absolute","none", "patJetCorrFactors").pt() << std::endl; 
    }

    if( jetID( jet ) < minJetID_ ){
      if(verbose_) cout << "@@@@ Jet does not pass the Id" << endl;
      continue;
    }
    if(jet->p4().Pt() < minCorrPt_) continue;

    sortedJetsID.insert(     make_pair( jet->p4().Pt() ,  jet->p4() ) );

    if(isMC_){
      if(jet->genJet() != 0) sortedGenJetsID.insert( make_pair( jet->p4().Pt() ,jet->genJet()->p4() ) );
      else sortedGenJetsID.insert( make_pair( jet->p4().Pt() , math::XYZTLorentzVectorD(0,0,0,0) ) );
    }
     
  }
  
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    jetsIDP4_->push_back( it->second );
  }
  for(CImap it = sortedGenJetsID.begin(); it != sortedGenJetsID.end() ; it++){
    genJetsIDP4_->push_back( it->second );
  }
  
  tree_->Fill();

}


unsigned int  MuTauAnalyzer::jetID( const pat::Jet* jet ){

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

  for(unsigned i=0; i<pfCandPtrs.size(); ++i) {
    const reco::PFCandidate& cand = *(pfCandPtrs[i]);

    totalEnergyFromConst +=  cand.energy();
    nConst += 1;

    switch( cand.particleId() ) {
    case reco::PFCandidate::h: 
      nCharged++;
      energyCharged += cand.energy(); 
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
      break;
    case reco::PFCandidate::h_HF: // fill neutral
      nNeutral++;
      //energyNeutral += cand.energy();
      break;
    case reco::PFCandidate::egamma_HF: // fill e/gamma
      nPhotons++;
      energyPhotons += cand.energy();
      break;
    default:
      break;
    }
  }

  bool loose=false;
  bool medium=false;
  bool tight=false;

  if(verbose_){
    cout << "NeutrFRAC = " << energyNeutral/totalEnergyFromConst << endl;
    cout << "PhotoFRAC = " << energyPhotons/totalEnergyFromConst << endl;
    cout << "ChargFRAC = " << energyCharged/totalEnergyFromConst << endl;
    cout << "nConst = "    << nConst << endl;
    cout << "nCharged = "  << nCharged << endl;
    cout << "ElectFRAC ="  << energyElectrons/totalEnergyFromConst << endl;
  }


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


void MuTauAnalyzer::endJob(){}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(MuTauAnalyzer);


