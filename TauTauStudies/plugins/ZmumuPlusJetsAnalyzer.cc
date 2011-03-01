#include "Bianchi/TauTauStudies/interface/ZmumuPlusJetsAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <utility>
#include <map>

using namespace std;
using namespace reco;

typedef std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more>::iterator CImap;

ZmumuPlusJetsAnalyzer::ZmumuPlusJetsAnalyzer(const edm::ParameterSet & iConfig){

  diMuonTag_ = iConfig.getParameter<edm::InputTag>("diMuons");
  jetsTag_ = iConfig.getParameter<edm::InputTag>("jets");
  isMC_ =  iConfig.getParameter<bool>("isMC");

}

void ZmumuPlusJetsAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","Z mumu plus jets tree");


  jetsP4_          = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDP4_        = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  tagJetsIDP4_     = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDbyMjjP4_   = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  jetsIDbyDEtaP4_  = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
    

  muonsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diMuonP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();


  tree_->Branch("jetsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsP4_);
  tree_->Branch("jetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDP4_);
  tree_->Branch("jetsIDbyMjjP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDbyMjjP4_);
  tree_->Branch("jetsIDbyDEtaP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&jetsIDbyDEtaP4_);
  tree_->Branch("tagJetsIDP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&tagJetsIDP4_);

  tree_->Branch("muonsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&muonsP4_);
  tree_->Branch("diMuonP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diMuonP4_);

  tree_->Branch("run",&run_,"run/F");
  tree_->Branch("event",&event_,"event/F");
  tree_->Branch("Zmass",&Zmass_,"Zmass/F");
  tree_->Branch("ZdeltaPhi",&ZdeltaPhi_,"ZdeltaPhi/F");
  tree_->Branch("MET",&MET_,"MET/F");
  tree_->Branch("numPV",&numPV_,"numPV/F");

}


ZmumuPlusJetsAnalyzer::~ZmumuPlusJetsAnalyzer(){
  delete jetsP4_; delete muonsP4_; delete jetsIDP4_; delete tagJetsIDP4_;  delete jetsIDbyMjjP4_; delete jetsIDbyDEtaP4_;
}

void ZmumuPlusJetsAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  //if true select as tag jets those matched to the tag quarks, 
  //else just check genJet matching
  bool isVBFH = false;

  jetsP4_->clear();
  jetsIDP4_->clear();
  tagJetsIDP4_->clear();
  jetsIDbyMjjP4_->clear();
  jetsIDbyDEtaP4_->clear();
  muonsP4_->clear();
  diMuonP4_->clear();
  tagQuark1_ = 0;
  tagQuark2_ = 0;
  
  const JetCorrector* corrector = JetCorrector::getJetCorrector("ak5PFResidual", iSetup);

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

  const reco::GenParticleCollection* genParticles = 0;
  edm::Handle<reco::GenParticleCollection> genHandle;
  edm::InputTag genTag("genParticles");
  iEvent.getByLabel(genTag,genHandle);
  if( isVBFH && isMC_ && genHandle.isValid() ){
    genParticles = genHandle.product();
    const reco::GenParticle *g1 = 0; 
    const reco::GenParticle *g2 = 0; 
    for(reco::GenParticleCollection::const_iterator ci = genParticles->begin(); ci!=genParticles->end(); ci++){
      //if( (TMath::Abs(ci->pdgId())==1 || TMath::Abs(ci->pdgId())==2 || TMath::Abs(ci->pdgId())==3 || TMath::Abs(ci->pdgId())==4 || TMath::Abs(ci->pdgId())==5)  ) 
      if(ci->pdgId() == 25){
	g1 = &(*(ci++));
	g2 = &(*((ci++)++));
	tagQuark1_ = g1->pt()>g2->pt() ? g1 : g2;
	tagQuark2_ = g1->pt()>g2->pt() ? g2 : g1;
	cout << "tag1 " << tagQuark1_->pdgId() << " tag2 " << tagQuark2_->pdgId() << endl;
	break;
      }
    }
  }

  numPV_ = vertexes->size();

  std::map< double, const CompositeCandidate*, less<double> > visMassMap;
  const CompositeCandidate *theDiMuon = 0;
  for(unsigned int it = 0; it < diMuons->size() ; it++){
    visMassMap.insert( make_pair(TMath::Abs(((*diMuons)[it]).mass()-91.2),
				 &(*diMuons)[it]) );
  }

  theDiMuon = visMassMap.size()>0 ?  &(*(visMassMap.begin()->second)) : 0;
  if(theDiMuon==0){
    cout << " No valid diMuon !!! " << endl;
    return;
  }
  
  if(theDiMuon->numberOfDaughters()<2){
    //std::cout << "less then 2 muons !!! " << std::endl;
    return;
  } //else {  std::cout << "passed !!! " << std::endl; }

  Zmass_ = theDiMuon->mass();
  ZdeltaPhi_ = Geom::deltaPhi(theDiMuon->daughter(0)->phi(),theDiMuon->daughter(1)->phi());
  MET_ = -999. ; //put MET


  const pat::Muon* mu = dynamic_cast<const pat::Muon*>( (theDiMuon->daughter(0)->masterClone()).get() );
  std::cout << mu->pt() << std::endl;

  muonsP4_->push_back(theDiMuon->daughter(0)->p4());
  muonsP4_->push_back(theDiMuon->daughter(1)->p4());
  diMuonP4_->push_back(theDiMuon->daughter(0)->p4()+theDiMuon->daughter(1)->p4());

  run_ = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedTagJetsID;

  for(unsigned int it = 0; it < jets->size() ; it++){
    double scale = isMC_ ? 1.0 : corrector->correction( (*jets)[it].p4() );
    cout << "scale factor for jet " << scale << endl;
   
    if((*jets)[it].p4().Pt()*scale < 10) continue;
    
    sortedJets.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );
    
    if( jetID( &(*jets)[it], scale ) > 0.5 ){
      
      sortedJetsID.insert( make_pair( (*jets)[it].p4().Pt()*scale ,scale*(*jets)[it].p4() ) );
     
      if( isVBFH && isMC_ && (*jets)[it].genJet()!=0 && tagQuark1_!=0 && tagQuark2_!=0 &&
	  (Geom::deltaR(tagQuark1_->p4(), ((*jets)[it].genJet())->p4()) < 0.3 ||
	   Geom::deltaR(tagQuark2_->p4(), ((*jets)[it].genJet())->p4()) < 0.3)
	  ) sortedTagJetsID.insert( make_pair( (*jets)[it].p4().Pt()*scale ,scale*(*jets)[it].p4() ) );
      if( !isVBFH && isMC_ && (*jets)[it].genJet()!=0 ) sortedTagJetsID.insert( make_pair( (*jets)[it].p4().Pt()*scale ,scale*(*jets)[it].p4() ) );
    }
  }
  
  for(CImap it = sortedJets.begin(); it != sortedJets.end() ; it++){
    jetsP4_->push_back( it->second );
  }
  for(CImap it = sortedJetsID.begin(); it != sortedJetsID.end() ; it++){
    jetsIDP4_->push_back( it->second );
  }
  for(CImap it = sortedTagJetsID.begin(); it != sortedTagJetsID.end() ; it++){
    tagJetsIDP4_->push_back( it->second );
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


unsigned int ZmumuPlusJetsAnalyzer::jetID( const pat::Jet* jet, const double scale){

  if( (jet->pt()*scale)<10 ) return 99; // always pass jet ID

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
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1) || 
      (TMath::Abs(jet->eta())<2.4 && 
       energyNeutral/totalEnergyFromConst<1 && 
       energyPhotons/totalEnergyFromConst<1 &&
       nConst > 1 &&
       energyCharged/totalEnergyFromConst>0 &&
       nCharged>0 &&
       energyElectrons/totalEnergyFromConst<1
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



void ZmumuPlusJetsAnalyzer::endJob(){

  //file_->cd();
  //tree_->Write();
  //file_->Close();

}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(ZmumuPlusJetsAnalyzer);


