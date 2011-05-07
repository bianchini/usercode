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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <vector>
#include <utility>
#include <map>

using namespace std;
using namespace reco;

typedef std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more>::iterator CImap;

MuTauStreamAnalyzer::MuTauStreamAnalyzer(const edm::ParameterSet & iConfig){

  diTauTag_ = iConfig.getParameter<edm::InputTag>("diTaus");
  jetsTag_ = iConfig.getParameter<edm::InputTag>("jets");
  isMC_ =  iConfig.getParameter<bool>("isMC");
  deltaRLegJet_  =  iConfig.getUntrackedParameter<double>("deltaRLegJet",0.3);
  minCorrPt_  =  iConfig.getUntrackedParameter<double>("minCorrPt",10.);
  minJetID_   =  iConfig.getUntrackedParameter<double>("minJetID",0.5);
  applyTauSignalSel_ = iConfig.getParameter<bool>("applyTauSignalSel");
  verbose_ =  iConfig.getUntrackedParameter<bool>("verbose",false);
}

void MuTauStreamAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","qqH tree");

  tRandom_ = new TRandom3();
 
  jetsBtagHE_  = new std::vector< double >();
  jetsBtagHP_  = new std::vector< double >();

  tauXTriggers_= new std::vector< int >();

  jetsP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDP4_        = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genJetsIDP4_       = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauVisP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauCAP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauICAP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfit1P4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfit2P4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauSVfit3P4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genDiTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  METP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genMETP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  fpuweight_ = new PUWeight();

  tree_->Branch("jetsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsP4_);
  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("genJetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genJetsIDP4_);
  
  tree_->Branch("jetsBtagHE","std::vector<double> ",&jetsBtagHE_);
  tree_->Branch("jetsBtagHP","std::vector<double> ",&jetsBtagHP_);

  tree_->Branch("tauXTriggers","std::vector<int> ",&tauXTriggers_);

  tree_->Branch("diTauVisP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",  &diTauVisP4_);
  tree_->Branch("diTauCAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauCAP4_);
  tree_->Branch("diTauICAP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",   &diTauICAP4_);
  tree_->Branch("diTauSVfit1P4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfit1P4_);
  tree_->Branch("diTauSVfit2P4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfit2P4_);
  tree_->Branch("diTauSVfit3P4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauSVfit3P4_);

  tree_->Branch("diTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsP4_);
  tree_->Branch("genDiTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genDiTauLegsP4_);

  tree_->Branch("METP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&METP4_);
  tree_->Branch("genMETP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genMETP4_);
  tree_->Branch("sumEt",&sumEt_,"sumEt/F");
  tree_->Branch("MtLeg1",&MtLeg1_,"MtLeg1/F");

  tree_->Branch("chIsoLeg1",&chIsoLeg1_,"chIsoLeg1/F");
  tree_->Branch("nhIsoLeg1",&nhIsoLeg1_,"nhIsoLeg1/F");
  tree_->Branch("phIsoLeg1",&phIsoLeg1_,"phIsoLeg1/F");
  tree_->Branch("chIsoPULeg1",&chIsoPULeg1_,"chIsoPULeg1/F");
  tree_->Branch("nhIsoPULeg1",&nhIsoPULeg1_,"nhIsoPULeg1/F");
  tree_->Branch("phIsoPULeg1",&phIsoPULeg1_,"phIsoPULeg1/F");

  tree_->Branch("chIsoLeg2",&chIsoLeg2_,"chIsoLeg2/F");
  tree_->Branch("nhIsoLeg2",&nhIsoLeg2_,"nhIsoLeg2/F");
  tree_->Branch("phIsoLeg2",&phIsoLeg2_,"phIsoLeg2/F");
  tree_->Branch("dxy1",&dxy1_,"dxy1/F");
  tree_->Branch("dxy2",&dxy2_,"dxy2/F");

  tree_->Branch("run",&run_,"run/F");
  tree_->Branch("event",&event_,"event/F");
  tree_->Branch("lumi",&lumi_,"lumi/F");
  tree_->Branch("numPV",&numPV_,"numPV/F");
  tree_->Branch("numOfDiTaus",&numOfDiTaus_,"numOfDiTaus/I");
  tree_->Branch("decayMode",&decayMode_,"decayMode/I");
  tree_->Branch("tightestHPSWP",&tightestHPSWP_,"tightestHPSWP/I");
  tree_->Branch("visibleTauMass",&visibleTauMass_,"visibleTauMass/F");
  tree_->Branch("leadPFChargedHadrCandTrackPt",&leadPFChargedHadrCandTrackPt_,"leadPFChargedHadrCandTrackPt/F");
  

  tree_->Branch("isTauLegMatched",&isTauLegMatched_,"isTauLegMatched/I");
  tree_->Branch("isMuLegMatched",&isMuLegMatched_,"isMuLegMatched/I");

  tree_->Branch("diTauCharge",&diTauCharge_,"diTauCharge/F");
  tree_->Branch("rhoFastJet",&rhoFastJet_,"rhoFastJet/F");
  tree_->Branch("mcPUweight",&mcPUweight_,"mcPUweight/F");


}


MuTauStreamAnalyzer::~MuTauStreamAnalyzer(){
  delete jetsP4_; delete jetsIDP4_; delete METP4_; delete diTauVisP4_; delete diTauCAP4_; delete diTauICAP4_; 
  delete diTauSVfit1P4_; delete diTauSVfit2P4_; delete diTauSVfit3P4_;
  delete diTauLegsP4_; delete jetsBtagHE_; delete jetsBtagHP_; delete tauXTriggers_;
  delete genJetsIDP4_; delete genDiTauLegsP4_; delete genMETP4_;
  delete tRandom_ ; delete fpuweight_;
}

void MuTauStreamAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){



  jetsP4_->clear();
  jetsIDP4_->clear();
  diTauVisP4_->clear();
  diTauCAP4_->clear();
  diTauICAP4_->clear();
  diTauSVfit1P4_->clear();
  diTauSVfit2P4_->clear();
  diTauSVfit3P4_->clear();
  diTauLegsP4_->clear();
  METP4_->clear();

  genJetsIDP4_->clear();
  genDiTauLegsP4_->clear();
  genMETP4_->clear();
  jetsBtagHE_->clear();
  jetsBtagHP_->clear();
  tauXTriggers_->clear();
  
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
  edm::InputTag pvTag("offlinePrimaryVerticesWithBS");
  iEvent.getByLabel(pvTag,pvHandle);
  if( !pvHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No PV label available \n";
  const reco::VertexCollection* vertexes = pvHandle.product();
  numPV_ = vertexes->size();

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel(edm::InputTag("patMETsPFlow"),metHandle);
  if( !metHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No MET label available \n";
  const pat::METCollection* met = metHandle.product();

  edm::Handle<reco::GenJetCollection> tauGenJetsHandle;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
  nPUVertices_= -99;
  nOOTPUVertices_=-99;

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
	if(it->getBunchCrossing() ==0) nPUVertices_ = it->getPU_NumInteractions();
	else  nOOTPUVertices_ = it->getPU_NumInteractions();
      }
    }
  }
  //cout << "Num of PU = " << nPUVertices_ << endl;
  //cout << "Num of OOT PU = " << nOOTPUVertices_ << endl;
  mcPUweight_ = fpuweight_->GetWeight(nPUVertices_);
  //cout << "Weight: " << weight << endl;


  edm::Handle<double> rhoFastJetHandle;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentral","rho", ""), rhoFastJetHandle);
  if( !rhoFastJetHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No rho label available \n";
  rhoFastJet_ = (*rhoFastJetHandle);
  
  
  const PATMuTauPair *theDiTau = 0;
  if(diTaus->size()<1){
    cout << " No diTau !!! " << endl;
    return;
  } else if(diTaus->size()>1 && verbose_){
    cout << "WARNING: "<< diTaus->size() << "  diTaus found in the event !!! We will select only one" << endl;
  }
  // choose the diTau with most isolated tau leg

  double sumIsoTau = 999.;
  //double highestPt = 0.;
  unsigned int index = 0;

  std::vector<unsigned int> identifiedTaus;
  std::vector<unsigned int> looseIsoTaus;
  //std::vector<unsigned int> not_identifiedTaus;

  for(unsigned int i=0; i<diTaus->size(); i++){
    const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[i].leg2()).get() );
    if(tau_i->tauID("leadingTrackFinding")<0.5) continue;
    identifiedTaus.push_back(i);
    if(tau_i->tauID("byLooseIsolation")>0.5) looseIsoTaus.push_back(i);
    //else not_identifiedTaus.push_back(i);
  }

  if(looseIsoTaus.size()>0 /*&& applyTauSignalSel_*/ ){
    identifiedTaus.swap(looseIsoTaus);
    if(verbose_) cout << identifiedTaus.size() << "  isolated taus found..." << endl;
    for(unsigned int i=0; i<identifiedTaus.size(); i++){
      if(verbose_) cout << "Testing isolation of " << i << "th tau" << endl;
      const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[ identifiedTaus[i] ].leg2()).get() );
      double sumIsoTau_i = 0.;
      sumIsoTau_i += tau_i->isolationPFChargedHadrCandsPtSum();
      sumIsoTau_i += tau_i->isolationPFGammaCandsEtSum();
      //sumIsoTau_i += tau_i->isolationPFNeutrHadrCandsEtSum();
      if(sumIsoTau_i<sumIsoTau){
	index = identifiedTaus[i];
	sumIsoTau = sumIsoTau_i;
      } 
    }
  } 
  else if(identifiedTaus.size()>0 /*&& !applyTauSignalSel_*/ ) {
    index = tRandom_->Integer( identifiedTaus.size() );
    if(verbose_) cout << "Random selection has chosen index " << index << endl;   
  }

    /*
      for(unsigned int i=0; identifiedTaus.size()>0 && !applyTauSignalSel_ && i<not_identifiedTaus.size(); i++){ 
      //if(identifiedTaus.size()>0) continue;
      if(verbose_) cout << "Testing max pt of charged hadron in the " << i << "th tau" << endl;
      const pat::Tau*  tau_i = dynamic_cast<const pat::Tau*>(  ((*diTaus)[ not_identifiedTaus[i] ].leg2()).get() );
      double highestPt_i = 0.;
      if(!(tau_i->pfTauTagInfoRef()).isNonnull()) continue;
      const PFCandidateRefVector pfChTau_i = tau_i->pfTauTagInfoRef()->PFChargedHadrCands();
      for(unsigned k = 0; k<pfChTau_i.size(); k++){
      PFCandidateRef hadron = pfChTau_i.at(k);
      if(hadron->pt()>highestPt_i){
      highestPt_i = hadron->pt();
      }
      }
      if(highestPt_i>highestPt){
      index =  not_identifiedTaus[i] ;
      highestPt = highestPt_i;
      }
      }
    */


  if(verbose_) cout << "Chosen index " << index << endl;
  identifiedTaus.clear(); /*not_identifiedTaus.clear(); */ looseIsoTaus.clear();

  theDiTau = &(*diTaus)[index];

  numOfDiTaus_ = diTaus->size();
  diTauCharge_ = theDiTau->charge();
  METP4_->push_back((*met)[0].p4());
  if(isMC_) genMETP4_->push_back( (*met)[0].genMET()->p4() );
  sumEt_  = (*met)[0].sumEt();
  MtLeg1_ =  theDiTau->mt1MET();
  isMuLegMatched_  = 0;
  isTauLegMatched_ = 0;

  const pat::Muon* leg1 = dynamic_cast<const pat::Muon*>( (theDiTau->leg1()).get() );
  const pat::Tau*  leg2 = dynamic_cast<const pat::Tau*>(  (theDiTau->leg2()).get() );

  vector<string> triggerPaths;
  if(isMC_){
    triggerPaths.push_back("HLT_Mu11_PFTau15_v*");
    triggerPaths.push_back("HLT_IsoMu9_PFTau15_v*");
  }
  else{
    triggerPaths.push_back("HLT_IsoMu12_LooseIsoPFTau10_v*");
    triggerPaths.push_back("HLT_Mu15_LooseIsoPFTau20_v*");
  }

  for(unsigned int m = 0; m<triggerPaths.size(); m++){
    if((leg1->triggerObjectMatchesByPath(triggerPaths[m],false)).size()!=0 && 
       (leg2->triggerObjectMatchesByPath(triggerPaths[m],false)).size()!=0) tauXTriggers_->push_back(1);
    else if((leg1->triggerObjectMatchesByPath(triggerPaths[m],false)).size()!=0 && 
	    (leg2->triggerObjectMatchesByPath(triggerPaths[m],false)).size()==0)  tauXTriggers_->push_back(2);
    else if((leg1->triggerObjectMatchesByPath(triggerPaths[m],false)).size()==0 && 
	    (leg2->triggerObjectMatchesByPath(triggerPaths[m],false)).size()!=0)  tauXTriggers_->push_back(3);
    else tauXTriggers_->push_back(0);
  }

  // triggers Mu
  if(verbose_){
    const pat::TriggerObjectStandAloneCollection trColl = leg1->triggerObjectMatchesByType(83);
    cout << "Mu triggers" << endl;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false)).size(); k++){
	cout << (it->pathNames(false))[k] << endl;
      }
    }
  }
  // triggers Tau
  if(verbose_){
    const pat::TriggerObjectStandAloneCollection trColl = leg2->triggerObjectMatchesByType(84);
    cout << "Tau triggers" << endl;
    for(pat::TriggerObjectStandAloneCollection::const_iterator it = trColl.begin(); it != trColl.end(); it++){
      for(unsigned int k = 0; k < (it->pathNames(false)).size(); k++){
	cout << (it->pathNames(false))[k] << endl;
      }
    }
  }

  diTauLegsP4_->push_back(leg1->p4());
  diTauLegsP4_->push_back(leg2->p4());
  
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

    if( leg2->genJet() !=0 ) genDiTauLegsP4_->push_back(leg2->genJet()->p4());
    else{
      genDiTauLegsP4_->push_back( math::XYZTLorentzVectorD(0,0,0,0) );
      if(verbose_) cout << "WARNING: no genJet matched to the leg2 with eta,phi " << leg2->eta() << ", " << leg2->phi() << endl;
    }

    bool tauHadMatched = false;
    for(unsigned int k = 0; k < tauGenJets->size(); k++){
      if( Geom::deltaR( (*tauGenJets)[k].p4(),leg2->p4() ) < 0.15 ) tauHadMatched = true;
    }

    if( (leg2->genParticleById(15,0,true)).isNonnull() && tauHadMatched ) isTauLegMatched_ = 1;
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
  leadPFChargedHadrCandTrackPt_ = (leg2->leadPFChargedHadrCand()->trackRef()).isNonnull() ?
    leg2->leadPFChargedHadrCand()->trackRef()->pt() : -99;

  tightestHPSWP_ = 0;
  if(leg2->tauID("byLooseIsolation")>0.5)  tightestHPSWP_++;
  if(leg2->tauID("byMediumIsolation")>0.5) tightestHPSWP_++;
  if(leg2->tauID("byTightIsolation")>0.5)  tightestHPSWP_++;


  dxy1_ = vertexes->size()!=0 ? leg1->globalTrack()->dxy( (*vertexes)[0].position() ) : -999;
  dxy2_ = (vertexes->size()!=0 && 
	   (leg2->leadPFChargedHadrCand()).isNonnull() 
	   && (leg2->leadPFChargedHadrCand()->trackRef()).isNonnull() ) ? 
    leg2->leadPFChargedHadrCand()->trackRef()->dxy( (*vertexes)[0].position() ) : -999;

  isodeposit::AbsVetos vetosChargedLeg1; 
  isodeposit::AbsVetos vetosNeutralLeg1; 
  isodeposit::AbsVetos vetosPhotonLeg1;
 
  vetosChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetosNeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.08));
  vetosNeutralLeg1.push_back(new isodeposit::ThresholdVeto(1.0));
  vetosPhotonLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.05));
  vetosPhotonLeg1.push_back(new isodeposit::ThresholdVeto(1.0));
  

  chIsoLeg1_   = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetosChargedLeg1).first;
  nhIsoLeg1_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetosNeutralLeg1).first;
  phIsoLeg1_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetosPhotonLeg1).first;
  chIsoPULeg1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetosChargedLeg1).first;
  nhIsoPULeg1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetosNeutralLeg1).first;
  phIsoPULeg1_ = 
    leg1->isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4,vetosPhotonLeg1).first;

  chIsoLeg2_ = -99;
  nhIsoLeg2_ = -99;
  phIsoLeg2_ = -99;
   

  // cleaning
  for(unsigned int i = 0; i <vetosChargedLeg1.size(); i++){
    delete vetosChargedLeg1[i];
  }
  for(unsigned int i = 0; i <vetosNeutralLeg1.size(); i++){
    delete vetosNeutralLeg1[i];
    delete vetosPhotonLeg1[i];
  }
  //

  diTauVisP4_->push_back( theDiTau->p4Vis() );
  diTauCAP4_->push_back( theDiTau->p4CollinearApprox() );
  diTauICAP4_->push_back( theDiTau->p4ImprovedCollinearApprox() );
  diTauSVfit1P4_->push_back( theDiTau->svFitSolution("psKine","",0)->p4()  );
  diTauSVfit2P4_->push_back( theDiTau->svFitSolution("psKine_MEt","",0)->p4()  );
  diTauSVfit3P4_->push_back( theDiTau->svFitSolution("psKine_MEt_ptBalance","",0)->p4()  );

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  lumi_ = iEvent.luminosityBlock();

  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,MuTauStreamAnalyzer::more> sortedGenJetsID;

  for(unsigned int it = 0; it < jets->size() ; it++){

    math::XYZTLorentzVectorD leg2p4 = ( (leg2->pfJetRef()).isNonnull() ) ? leg2->pfJetRef()->p4() : leg2->p4();

    if( Geom::deltaR((*jets)[it].p4(),leg1->p4())<deltaRLegJet_ || 
	Geom::deltaR((*jets)[it].p4(), leg2p4 )<deltaRLegJet_ ){
      if(verbose_) cout << "The jet at (" <<(*jets)[it].pt()<<","<<(*jets)[it].eta()<<") is closer than "<<deltaRLegJet_ << " from one of the legs" << endl;  
      continue;
    }

    if(verbose_){
      pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);
      //for(unsigned int i = 0; i < (jet->availableJECLevels()).size() ; i++ ){
      //std::cout << (jet->availableJECLevels())[i] << std::endl;
      //}
      std::cout << "Uncorrected " << jet->correctedJet("Uncorrected").pt() << std::endl;
      std::cout << "L1FastJet "   << jet->correctedJet("L1FastJet").pt() << std::endl;
      std::cout << "L2Relative "  << jet->correctedJet("L2Relative").pt() << std::endl; 
      std::cout << "L3Absolute "  << jet->correctedJet("L3Absolute").pt() << std::endl; 
    }

    if( jetID( &(*jets)[it] ) < minJetID_ )  continue;

    sortedJets.insert( make_pair( (*jets)[it].correctedJet("Uncorrected").p4().Pt() ,(*jets)[it].correctedJet("Uncorrected").p4() ) );

    if((*jets)[it].p4().Pt() < minCorrPt_) continue;

    //add b-tag info
    jetsBtagHE_->push_back((*jets)[it].bDiscriminator("trackCountingHighEffBJetTags"));
    jetsBtagHP_->push_back((*jets)[it].bDiscriminator("trackCountingHighPurBJetTags"));
                                
    sortedJetsID.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );
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
  for(CImap it = sortedGenJetsID.begin(); it != sortedGenJetsID.end() ; it++){
    genJetsIDP4_->push_back( it->second );
  }


  tree_->Fill();

}


unsigned int MuTauStreamAnalyzer::jetID( const pat::Jet* jet){

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


