#include "Bianchi/TauTauStudies/interface/ZmumuPlusJetsAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <vector>
#include <utility>
#include <map>

using namespace std;
using namespace reco;

typedef std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more>::iterator CImap;

ZmumuPlusJetsAnalyzer::ZmumuPlusJetsAnalyzer(const edm::ParameterSet & iConfig){

  diMuonTag_ = iConfig.getParameter<edm::InputTag>("diMuons");
  jetsTag_ = iConfig.getParameter<edm::InputTag>("jets");
  isMC_ =  iConfig.getParameter<bool>("isMC");
  minCorrPt_  =  iConfig.getUntrackedParameter<double>("minCorrPt",10.);
  minJetID_   =  iConfig.getUntrackedParameter<double>("minJetID",0.5);
  verbose_ =  iConfig.getUntrackedParameter<bool>("verbose",false);
}

void ZmumuPlusJetsAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","Z mumu plus jets tree");

  jetsBtagHE_  = new std::vector< double >();
  jetsBtagHP_  = new std::vector< double >();

  jetsP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDP4_        = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDbyMjjP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDbyDEtaP4_  = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
    

  muonsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diMuonP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();


  tree_->Branch("jetsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsP4_);
  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("jetsIDbyMjjP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDbyMjjP4_);
  tree_->Branch("jetsIDbyDEtaP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDbyDEtaP4_);
  
  tree_->Branch("jetsBtagHE","std::vector<double> ",&jetsBtagHE_);
  tree_->Branch("jetsBtagHP","std::vector<double> ",&jetsBtagHP_);

  tree_->Branch("muonsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&muonsP4_);
  tree_->Branch("diMuonP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diMuonP4_);

  tree_->Branch("METP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&METP4_);
  tree_->Branch("sumEt",&sumEt_,"sumEt/F");
  tree_->Branch("MtLeg1",&MtLeg1_,"MtLeg1/F");
  tree_->Branch("MtLeg2",&MtLeg2_,"MtLeg2/F");

  tree_->Branch("chIsoLeg1",&chIsoLeg1_,"chIsoLeg1/F");
  tree_->Branch("nhIsoLeg1",&nhIsoLeg1_,"nhIsoLeg1/F");
  tree_->Branch("phIsoLeg1",&phIsoLeg1_,"phIsoLeg1/F");
  tree_->Branch("chIsoLeg2",&chIsoLeg2_,"chIsoLeg2/F");
  tree_->Branch("nhIsoLeg2",&nhIsoLeg2_,"nhIsoLeg2/F");
  tree_->Branch("phIsoLeg2",&phIsoLeg2_,"phIsoLeg2/F");
  tree_->Branch("dxy1",&dxy1_,"dxy1/F");
  tree_->Branch("dxy2",&dxy2_,"dxy2/F");

  tree_->Branch("run",&run_,"run/F");
  tree_->Branch("event",&event_,"event/F");
  tree_->Branch("Zmass",&Zmass_,"Zmass/F");
  tree_->Branch("ZdeltaPhi",&ZdeltaPhi_,"ZdeltaPhi/F");
  tree_->Branch("numPV",&numPV_,"numPV/F");

}


ZmumuPlusJetsAnalyzer::~ZmumuPlusJetsAnalyzer(){
  delete jetsP4_; delete muonsP4_; delete jetsIDP4_; 
  delete jetsIDbyMjjP4_; delete jetsIDbyDEtaP4_; 
  delete METP4_; delete diMuonP4_;
}

void ZmumuPlusJetsAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){


  jetsP4_->clear();
  jetsIDP4_->clear();
  jetsIDbyMjjP4_->clear();
  jetsIDbyDEtaP4_->clear();
  jetsBtagHE_->clear();
  jetsBtagHP_->clear();
  muonsP4_->clear();
  diMuonP4_->clear();
  METP4_->clear();

  
  edm::Handle<CompositeCandidateCollection> diMuonHandle;
  iEvent.getByLabel(diMuonTag_,diMuonHandle);
  if( !diMuonHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No diMuon label available \n";
  const CompositeCandidateCollection* diMuons = diMuonHandle.product();

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

  std::map< double, const CompositeCandidate*, less<double> > visMassMap;
  const CompositeCandidate *theDiMuon = 0;
  for(unsigned int it = 0; it < diMuons->size() ; it++){
    visMassMap.insert( make_pair(TMath::Abs(((*diMuons)[it]).mass()-91.2),
				 &(*diMuons)[it]) );
  }

  theDiMuon = visMassMap.size()>0 ?  &(*(visMassMap.begin()->second)) : 0;
  if(theDiMuon==0 || theDiMuon->numberOfDaughters()<2){
    cout << " No valid diMuon !!! " << endl;
    return;
  }
  

  Zmass_ = theDiMuon->mass();
  ZdeltaPhi_ = abs(Geom::deltaPhi(theDiMuon->daughter(0)->phi(),theDiMuon->daughter(1)->phi()));
  METP4_->push_back((*met)[0].p4());
  sumEt_  = (*met)[0].sumEt();
  MtLeg1_ = TMath::Sqrt( (theDiMuon->daughter(0)->pt() + (*met)[0].pt() )*
			 (theDiMuon->daughter(0)->pt() + (*met)[0].pt() )- 
			 (theDiMuon->daughter(0)->p4() + (*met)[0].p4()).pt()*
			 (theDiMuon->daughter(0)->p4() + (*met)[0].p4()).pt()  );
  MtLeg2_ = TMath::Sqrt( (theDiMuon->daughter(1)->pt() + (*met)[0].pt() )*
			 (theDiMuon->daughter(1)->pt() + (*met)[0].pt() )- 
			 (theDiMuon->daughter(1)->p4() + (*met)[0].p4()).pt()*
			 (theDiMuon->daughter(1)->p4() + (*met)[0].p4()).pt()  );

  const pat::Muon* leg1 = dynamic_cast<const pat::Muon*>( (theDiMuon->daughter(0)->masterClone()).get() );
  const pat::Muon* leg2 = dynamic_cast<const pat::Muon*>( (theDiMuon->daughter(1)->masterClone()).get() );

  dxy1_ = vertexes->size()!=0 ? leg1->globalTrack()->dxy( (*vertexes)[0].position() ) : -999;
  dxy2_ = vertexes->size()!=0 ? leg2->globalTrack()->dxy( (*vertexes)[0].position() ) : -999;

  isodeposit::AbsVetos vetosChargedLeg1; 
  isodeposit::AbsVetos vetosNeutralLeg1; 
  isodeposit::AbsVetos vetosPhotonLeg1;
  isodeposit::AbsVetos vetosChargedLeg2; 
  isodeposit::AbsVetos vetosNeutralLeg2;
  isodeposit::AbsVetos vetosPhotonLeg2; 
  vetosChargedLeg1.push_back(new isodeposit::ThresholdVeto(0.5));
  vetosNeutralLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.08));
  vetosNeutralLeg1.push_back(new isodeposit::ThresholdVeto(1.0));
  vetosPhotonLeg1.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg1->eta(),leg1->phi()),0.05));
  vetosPhotonLeg1.push_back(new isodeposit::ThresholdVeto(1.0));

  vetosChargedLeg2.push_back(new isodeposit::ThresholdVeto(0.5)); 
  vetosNeutralLeg2.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg2->eta(),leg2->phi()),0.08));
  vetosNeutralLeg2.push_back(new isodeposit::ThresholdVeto(1.0));
  vetosPhotonLeg2.push_back(new isodeposit::ConeVeto(isodeposit::Direction(leg2->eta(),leg2->phi()),0.05));
  vetosPhotonLeg2.push_back(new isodeposit::ThresholdVeto(1.0));

  chIsoLeg1_ = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetosChargedLeg1).first;
  nhIsoLeg1_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetosNeutralLeg1).first;
  phIsoLeg1_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetosPhotonLeg1).first;
  chIsoLeg2_ = 
    leg2->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,vetosChargedLeg2).first;
  nhIsoLeg2_ = 
    leg2->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,vetosNeutralLeg2).first;
  phIsoLeg2_ = 
    leg2->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,vetosPhotonLeg2 ).first;

  // cleaning
  for(unsigned int i = 0; i <vetosChargedLeg1.size(); i++){
    delete vetosChargedLeg1[i];
    delete vetosChargedLeg2[i];
  }
  for(unsigned int i = 0; i <vetosNeutralLeg1.size(); i++){
    delete vetosNeutralLeg1[i];
    delete vetosNeutralLeg2[i];
    delete vetosPhotonLeg1[i];
    delete vetosPhotonLeg2[i];
  }
  //

  muonsP4_->push_back(theDiMuon->daughter(0)->p4());
  muonsP4_->push_back(theDiMuon->daughter(1)->p4());
  diMuonP4_->push_back(theDiMuon->daughter(0)->p4()+theDiMuon->daughter(1)->p4());

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJetsID;

  for(unsigned int it = 0; it < jets->size() ; it++){

    if(verbose_){
      pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);
      std::cout << "Raw " << jet->correctedJet("raw").pt() << std::endl;
      std::cout << "Rel " << jet->correctedJet("rel").pt() << std::endl;
      std::cout << "Abs " << jet->correctedJet("abs").pt() << std::endl; 
    }

    if((*jets)[it].p4().Pt() < minCorrPt_) continue;

    //add b-tag info
    jetsBtagHE_->push_back((*jets)[it].bDiscriminator("trackCountingHighEffBJetTags"));
    jetsBtagHP_->push_back((*jets)[it].bDiscriminator("trackCountingHighPurBJetTags"));
    
    sortedJets.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );
                                        
    if( jetID( &(*jets)[it] ) < minJetID_ )  continue;

    sortedJetsID.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );
     
  }
  
  for(CImap it = sortedJets.begin(); it != sortedJets.end() ; it++){
    jetsP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    jetsIDP4_->push_back( it->second );
  }

  double Mass = 0; double DEta = 0; CImap tag1,tag2;
  //by Mjj
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    for(CImap jt = it; jt != sortedJetsID.end() ; jt++){
      if( it!=jt && (it->second + jt->second).M() > Mass){
	Mass = (it->second + jt->second).M();
	tag1 = it; tag2 = jt;
      }
    }
  }
  if(sortedJetsID.size()>1){
    jetsIDbyMjjP4_->push_back( tag1->second );
    jetsIDbyMjjP4_->push_back( tag2->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    if(it==tag1 || it==tag2) continue;
    jetsIDbyMjjP4_->push_back( it->second );
  }
  //by delta eta
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    for(CImap jt = it; jt != sortedJetsID.end() ; jt++){
      if( it!=jt && TMath::Abs((it->second).eta()-(jt->second).eta()) > DEta){
	DEta =  TMath::Abs((it->second).eta()-(jt->second).eta());
	tag1 = it; tag2 = jt;
      }
    }
  }
  if(sortedJetsID.size()>1){
    jetsIDbyDEtaP4_->push_back( tag1->second );
    jetsIDbyDEtaP4_->push_back( tag2->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    if(it==tag1 || it==tag2) continue;
    jetsIDbyDEtaP4_->push_back( it->second );
  }


  tree_->Fill();

}


unsigned int ZmumuPlusJetsAnalyzer::jetID( const pat::Jet* jet){

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



void ZmumuPlusJetsAnalyzer::endJob(){}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(ZmumuPlusJetsAnalyzer);


