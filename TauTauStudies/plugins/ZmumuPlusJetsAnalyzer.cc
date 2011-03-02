#include "Bianchi/TauTauStudies/interface/ZmumuPlusJetsAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

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
  applyResidualJEC_ =  iConfig.getParameter<bool>("applyResidualJEC");
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
  delete tagJetsIDP4_;  delete jetsIDbyMjjP4_; delete jetsIDbyDEtaP4_; 
  delete METP4_; delete diMuonP4_;
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
  jetsBtagHE_->clear();
  jetsBtagHP_->clear();
  muonsP4_->clear();
  diMuonP4_->clear();
  METP4_->clear();

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
  numPV_ = vertexes->size();

  const reco::GenParticleCollection* genParticles = 0;
  edm::Handle<reco::GenParticleCollection> genHandle;
  edm::InputTag genTag("genParticles");
  iEvent.getByLabel(genTag,genHandle);
  if( isVBFH && isMC_ && genHandle.isValid() ){
    genParticles = genHandle.product();
    const reco::GenParticle *g1 = 0; 
    const reco::GenParticle *g2 = 0; 
    for(reco::GenParticleCollection::const_iterator ci = genParticles->begin(); ci!=genParticles->end(); ci++){
      if(ci->pdgId() == 25){
	g1 = &(*(ci++));
	g2 = &(*((ci++)++));
	tagQuark1_ = g1->pt()>g2->pt() ? g1 : g2;
	tagQuark2_ = g1->pt()>g2->pt() ? g2 : g1;
	if(verbose_) cout << "tag1 " << tagQuark1_->pdgId() << " tag2 " << tagQuark2_->pdgId() << endl;
	break;
      }
    }
  }

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

  chIsoLeg1_ = 
    leg1->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  nhIsoLeg1_ = 
    leg1->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  phIsoLeg1_ = 
    leg1->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  chIsoLeg2_ = 
    leg2->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  nhIsoLeg2_ = 
    leg2->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  phIsoLeg2_ = 
    leg2->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
  
  muonsP4_->push_back(theDiMuon->daughter(0)->p4());
  muonsP4_->push_back(theDiMuon->daughter(1)->p4());
  diMuonP4_->push_back(theDiMuon->daughter(0)->p4()+theDiMuon->daughter(1)->p4());

  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJets;
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedJetsID;
  std::map<double, math::XYZTLorentzVectorD ,ZmumuPlusJetsAnalyzer::more> sortedTagJetsID;

  for(unsigned int it = 0; it < jets->size() ; it++){
    double scale = (!isMC_ && applyResidualJEC_) ?  corrector->correction( (*jets)[it].p4() ) : 1.0 ;

    if(verbose_){
      std::cout << "L2L3Residual corrections " << scale << std::endl;
      pat::Jet* jet = const_cast<pat::Jet*>(&(*jets)[it]);
      std::cout << "Raw " << jet->correctedJet("raw").pt() << std::endl;
      std::cout << "Rel " << jet->correctedJet("rel").pt() << std::endl;
      std::cout << "Abs " << jet->correctedJet("abs").pt() << std::endl; 
    }

    if((*jets)[it].p4().Pt()*scale < minCorrPt_) continue;
    
    sortedJets.insert( make_pair( (*jets)[it].p4().Pt() ,(*jets)[it].p4() ) );

    //add b-tag info
    jetsBtagHE_->push_back((*jets)[it].bDiscriminator("trackCountingHighEffBJetTags"));
    jetsBtagHP_->push_back((*jets)[it].bDiscriminator("trackCountingHighPurBJetTags"));
                                            
    if( jetID( &(*jets)[it], scale ) > minJetID_ ){
      
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


