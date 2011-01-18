#include "Bianchi/eToTaufakeRate/interface/TauFakeRateAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "Math/VectorUtil.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <utility>
#include <map>

using namespace std;
using namespace reco;


TauFakeRateAnalyzer::TauFakeRateAnalyzer(const edm::ParameterSet & iConfig){

  isMC_ =  iConfig.getParameter<bool>("isMC");
  matchTo_ =  iConfig.existsAs<string>("matchTo") ? iConfig.getParameter<string>("matchTo") : "default";
  tauTag_ = iConfig.getParameter<edm::InputTag>("tauTag");
  electronsTag_ = iConfig.getParameter<edm::InputTag>("electronsTag");
}

void TauFakeRateAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tau fake rate");
  tree_->Branch("leadPFChargedHadrMva",&leadPFChargedHadrMva_,"leadPFChargedHadrMva/D");
  tree_->Branch("leadPFChargedHadrHcalEnergy",&leadPFChargedHadrHcalEnergy_,"leadPFChargedHadrHcalEnergy/D");
  tree_->Branch("leadPFChargedHadrEcalEnergy",&leadPFChargedHadrEcalEnergy_,"leadPFChargedHadrEcalEnergy/D");
  tree_->Branch("leadPFChargedHadrTrackPt",&leadPFChargedHadrTrackPt_,"leadPFChargedHadrTrackPt/D");
  tree_->Branch("leadPFChargedHadrTrackP",&leadPFChargedHadrTrackP_,"leadPFChargedHadrTrackP/D");
  tree_->Branch("leadPFCandMva",&leadPFCandMva_,"leadPFCandMva/D");
  tree_->Branch("leadPFCandHcalEnergy",&leadPFCandHcalEnergy_,"leadPFCandHcalEnergy/D");
  tree_->Branch("leadPFCandEcalEnergy",&leadPFCandEcalEnergy_,"leadPFCandEcalEnergy/D");
  tree_->Branch("leadPFCandPt",&leadPFCandPt_,"leadPFCandPt/D");
  tree_->Branch("leadPFCandP",&leadPFCandP_,"leadPFCandP/D");
  tree_->Branch("matchedID",&matchedID_,"matchedID/D");
  tree_->Branch("fbrem",&fbrem_,"fbrem/D");
  tree_->Branch("hasGsf",&hasGsf_,"hasGsf/D");
  
}


TauFakeRateAnalyzer::~TauFakeRateAnalyzer(){}

void TauFakeRateAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  vector<math::XYZTLorentzVector> genMatchP4s;
  int ids[] = {95,90,85,80,70,60};

  edm::Handle<pat::TauCollection> tausHandle;
  iEvent.getByLabel(tauTag_,tausHandle);
  if( !tausHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No taus label available \n";
  const pat::TauCollection* taus = tausHandle.product();

  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByLabel(electronsTag_,electronsHandle);
  if( !electronsHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No electrons label available \n";
  const pat::ElectronCollection* electrons = electronsHandle.product();

  if(matchTo_.find("default")!=string::npos || matchTo_.find("tau")!=string::npos){
    edm::Handle<reco::GenJetCollection> tauGenJetsHandle;
    iEvent.getByLabel(edm::InputTag("genTauDecaysToHadrons"),tauGenJetsHandle);
    if( !tauGenJetsHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen jet label available \n";
    const reco::GenJetCollection* tauGenJets = tauGenJetsHandle.product();
    for(unsigned int j = 0; j<tauGenJets->size() ; j++){
      genMatchP4s.push_back( (*tauGenJets)[j].p4() );
    }
  } else if( matchTo_.find("electron")!=string::npos ){
    edm::Handle<reco::GenParticleCollection> genParticlesHandle;
    iEvent.getByLabel(edm::InputTag("genParticles"),genParticlesHandle);
    if( !genParticlesHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No genparticles label available \n";
    const reco::GenParticleCollection* genParticles = genParticlesHandle.product();
    for(unsigned int j = 0; j<genParticles->size() ; j++){
      if( TMath::Abs((*genParticles)[j].pdgId())!=11 ) continue;
      genMatchP4s.push_back( (*genParticles)[j].p4() );
    }
  } else edm::LogError("DataNotAvailable")
      << "not a valid matcher";

  for(unsigned int i = 0; i< taus->size(); i++){

    //cout << iEvent.id() << "  analyzing the " << i+1 << "  tau..."  << endl;
    if( ((*taus)[i].leadPFChargedHadrCand()).isNull() || ((*taus)[i].leadPFCand()).isNull() ){
      cout << "Null ref to charged or lead PF" << endl;
      if( ((*taus)[i].leadPFCand()).isNull() ) cout << "... no leadPF " << endl;
      continue;
    }

    // consider only 1-prong reco taus
    //if( ((*taus)[i].signalPFChargedHadrCands()).size()!=1) continue;

    bool matched = false;
    for(unsigned int j = 0; j<genMatchP4s.size() ; j++){
      if( !matched && ROOT::Math::VectorUtil::DeltaR( (*taus)[i].p4(), genMatchP4s[i] )<0.15){

	matched = true;
	matchedID_ = 1.0;

	// check if any gsfElectron is matched to the leadPFCahrgHadr by gsfTrack
	// if so, save the tightest cut-based WP, else 1.0
	for(unsigned int k = 0; k<electrons->size() ; k++){

	  if(((*electrons)[k].gsfTrack()).isNull()){
	    cout << "Null gsf trak for the electron!!!!" << endl;
	    continue;
	  }

	  // match by gsfTrack
	  bool matchedGsf        = ((((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()).isNull() ) ? false : (TMath::Abs(((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()->eta()-(*electrons)[k].gsfTrack()->eta()) < 1e-04 && TMath::Abs(((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()->phi()-(*electrons)[k].gsfTrack()->phi()) < 1e-04);
	  // match by track
	  bool matchedTrack      = (TMath::Abs(((*taus)[i].leadPFChargedHadrCand())->trackRef()->eta()-(*electrons)[k].gsfTrack()->eta()) < 1e-04 && TMath::Abs(((*taus)[i].leadPFChargedHadrCand())->trackRef()->phi()-(*electrons)[k].gsfTrack()->phi()) < 1e-04);
	  // match by ambiguous gsfTrack
	  bool matchedAmbGsf     = false;
	  for( reco::GsfTrackRefVector::const_iterator it = (*electrons)[k].ambiguousGsfTracksBegin() ; 
	       it!=(*electrons)[k].ambiguousGsfTracksEnd(); it++ ){
	    bool checkThis = ((((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()).isNull() ) ? false : ((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()==(*it);
	    matchedAmbGsf |= checkThis;
	  }
	  if( matchedGsf || matchedAmbGsf || matchedTrack ){
	    for(int id = 0; id<6; id++){
	      float match = (*electrons)[k].electronID(Form("simpleEleId%drelIso",ids[id]));
	      if( (match>0.5&&match<1.5)||
		  (match>2.5&&match<3.5)||
		  (match>4.5&&match<5.5)||
		  match>6.5 ) matchedID_ = 0.01*ids[id];
	    }
	  }
	}
	
	
	
	leadPFChargedHadrMva_        =   (*taus)[i].electronPreIDOutput() ;	
	leadPFChargedHadrHcalEnergy_ =  ((*taus)[i].leadPFChargedHadrCand())->hcalEnergy() ;
	leadPFChargedHadrEcalEnergy_ =  ((*taus)[i].leadPFChargedHadrCand())->ecalEnergy() ;
	leadPFChargedHadrTrackPt_    =  ((*taus)[i].leadPFChargedHadrCand())->trackRef()->pt();
	leadPFChargedHadrTrackP_     =  ((*taus)[i].leadPFChargedHadrCand())->trackRef()->p();
	leadPFCandMva_               =  ((*taus)[i].leadPFCand())->mva_e_pi() ;	
	leadPFCandHcalEnergy_        =  ((*taus)[i].leadPFCand())->hcalEnergy() ;
	leadPFCandEcalEnergy_        =  ((*taus)[i].leadPFCand())->ecalEnergy() ;
	leadPFCandPt_                =  ((*taus)[i].leadPFCand())->pt();
	leadPFCandP_                 =  ((*taus)[i].leadPFCand())->p();
	hasGsf_ = (((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()).isNonnull() ?
	  1 : 0;
	fbrem_ = (((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()).isNonnull() ?
	  ( ((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()->p() - 
	    ((*taus)[i].leadPFChargedHadrCand())->gsfTrackRef()->outerP() ) :
	  ( ((*taus)[i].leadPFChargedHadrCand())->trackRef()->p() - 
	    ((*taus)[i].leadPFChargedHadrCand())->trackRef()->outerP() ) ;

	tree_->Fill();
      }

    }
  }


}





void TauFakeRateAnalyzer::endJob(){
}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(TauFakeRateAnalyzer);


