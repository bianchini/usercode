#include "Bianchi/TauTauStudies/interface/MEtRecoilCorrectorProducer.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TLorentzVector.h"

#include <cstdio>
#include <map>

using namespace std;

MEtRecoilCorrectorProducer::MEtRecoilCorrectorProducer(const edm::ParameterSet & iConfig){

  metTag_          = iConfig.getParameter<edm::InputTag>("metTag");
  genParticlesTag_ = iConfig.getParameter<edm::InputTag>("genParticleTag");
  muonTag_         = iConfig.getParameter<edm::InputTag>("muonTag");
  electronTag_     = iConfig.getParameter<edm::InputTag>("electronTag");
  tauTag_          = iConfig.getParameter<edm::InputTag>("tauTag");
  jetTag_          = iConfig.getParameter<edm::InputTag>("jetTag");
  minJetPt_        = iConfig.getParameter<double>("minJetPt");

  verbose_  = iConfig.getParameter<bool>("verbose");
  isMC_     = iConfig.getParameter<bool>("isMC");
  
  inputFileNamezmm42X_ = iConfig.getParameter<edm::FileInPath>("inputFileNamezmm42X");
  inputFileNamedatamm_ = iConfig.getParameter<edm::FileInPath>("inputFileNamedatamm");
  inputFileNamewjets_  = iConfig.getParameter<edm::FileInPath>("inputFileNamewjets");
  inputFileNamezjets_  = iConfig.getParameter<edm::FileInPath>("inputFileNamezjets");
  inputFileNamehiggs_  = iConfig.getParameter<edm::FileInPath>("inputFileNamehiggs");
  
  eventCounter_ = -1;
  recoilCorr_   =  0;

  produces<pat::METCollection>("U");     
  produces<pat::METCollection>("N");     
  produces<pat::METCollection>("D");  

}

MEtRecoilCorrectorProducer::~MEtRecoilCorrectorProducer(){
  if(recoilCorr_!=0) delete recoilCorr_;
}

void MEtRecoilCorrectorProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  eventCounter_ ++;

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel(metTag_,metHandle);
  const pat::METCollection* mets = metHandle.product();

  double scaledMETPxN  = (*mets)[0].px();
  double scaledMETPyN  = (*mets)[0].py();
  double scaledMETPtN  = (*mets)[0].pt();
  double scaledMETPhiN = (*mets)[0].phi();

  double scaledMETPxU  = (*mets)[0].px();
  double scaledMETPyU  = (*mets)[0].py();
  double scaledMETPtU  = (*mets)[0].pt();
  double scaledMETPhiU = (*mets)[0].phi();

  double scaledMETPxD  = (*mets)[0].px();
  double scaledMETPyD  = (*mets)[0].py();
  double scaledMETPtD  = (*mets)[0].pt();
  double scaledMETPhiD = (*mets)[0].phi();

  std::auto_ptr< pat::METCollection > MEtRescaledCollU( new pat::METCollection() ) ;
  std::auto_ptr< pat::METCollection > MEtRescaledCollN( new pat::METCollection() ) ;
  std::auto_ptr< pat::METCollection > MEtRescaledCollD( new pat::METCollection() ) ;

  if(isMC_ && eventCounter_%1000==0){

    if(verbose_) cout << "Is MC, find out which process it is..." << endl;

    edm::Handle<reco::GenParticleCollection> genHandle;
    const reco::GenParticleCollection* genParticles = 0;
    iEvent.getByLabel(edm::InputTag("genParticles"),genHandle);
    if( !genHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen particles label available \n";
    genParticles = genHandle.product();
    
    int counterW = 0;
    int counterZ = 0;
    int counterH = 0;

    for(unsigned int k = 0; k < genParticles->size(); k ++){
      reco::GenParticle gen = (*genParticles)[k];
      if( gen.pdgId()== 23       && gen.status()==3 ) counterZ++;
      if( fabs(gen.pdgId())== 24 && gen.status()==3 ) counterW++;
      if( gen.pdgId()== 25       && gen.status()==3 ) counterH++;
    }
    if(verbose_) cout << "Zs = " << counterZ << ", Ws = " << counterW << ", Hs = " << counterH << endl;

    if(counterZ==1 && counterW==0 && counterH==0){
      genDecay_ = 23;
      recoilCorr_ = new  RecoilCorrector(inputFileNamezjets_.fullPath().data());
      recoilCorr_->addMCFile(inputFileNamezmm42X_.fullPath().data());
      recoilCorr_->addDataFile(inputFileNamedatamm_.fullPath().data());
      if(verbose_) cout << "Loading corrections for DY" << endl;
    }
    else if(counterZ==0 && counterW==1 && counterH==0){
      genDecay_ = 24;
      recoilCorr_ = new  RecoilCorrector(inputFileNamewjets_.fullPath().data());
      recoilCorr_->addMCFile(inputFileNamezmm42X_.fullPath().data());
      recoilCorr_->addDataFile(inputFileNamedatamm_.fullPath().data());
      if(verbose_) cout << "Loading corrections for W+jets" << endl;
    }
    else if(counterZ==0 && counterW==0 && counterH==1){
      genDecay_ = 25;
      recoilCorr_ = new  RecoilCorrector(inputFileNamehiggs_.fullPath().data());
      recoilCorr_->addMCFile(inputFileNamezmm42X_.fullPath().data());
      recoilCorr_->addDataFile(inputFileNamedatamm_.fullPath().data());
      if(verbose_) cout << "Loading corrections for Higgs" << endl;
    }
    else{
      genDecay_ = -99;
      if(verbose_) cout << "The event is not correctable!" << endl;
      recoilCorr_ = 0;
    }
  }


  if(recoilCorr_!=0){

    if(verbose_) cout << "Evaluate recoil..." << endl;

    edm::Handle<pat::MuonCollection> muonHandle;
    const pat::MuonCollection* muons = 0;
    if(muonTag_.label()!=""){
      iEvent.getByLabel(muonTag_,muonHandle);
      muons = muonHandle.product();
    }
    
    edm::Handle<pat::ElectronCollection> electronHandle;
    const pat::ElectronCollection* electrons = 0;
    if(electronTag_.label()!=""){
      iEvent.getByLabel(electronTag_,electronHandle);
      electrons = electronHandle.product();
    }
    
    edm::Handle<pat::TauCollection> tauHandle;
    const pat::TauCollection* taus = 0;
    if(tauTag_.label()!=""){
      iEvent.getByLabel(tauTag_,tauHandle);
      taus = tauHandle.product();
    }

    edm::Handle<reco::GenParticleCollection> genHandle;
    const reco::GenParticleCollection* genParticles = 0;
    iEvent.getByLabel(edm::InputTag("genParticles"),genHandle);
    if( !genHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No gen particles label available \n";
    genParticles = genHandle.product();
    
    edm::Handle<pat::JetCollection> jetsHandle;
    iEvent.getByLabel(jetTag_,jetsHandle);
    if( !jetsHandle.isValid() )  
      edm::LogError("DataNotAvailable")
	<< "No jets label available \n";
    const pat::JetCollection* jets = jetsHandle.product();
    
    int nJets30 = 0;
    for(unsigned int jet = 0 ; jet<jets->size(); jet++){
      if((*jets)[jet].pt()>minJetPt_ && TMath::Abs((*jets)[jet].eta())<4.5) 
	nJets30++;
    }
    if(verbose_) cout << "nJets30 = " << nJets30 << endl;

    math::XYZTLorentzVectorD genVP4;
    math::XYZTLorentzVectorD genLeg1P4(0,0,0,0);
    math::XYZTLorentzVectorD genLeg2P4(0,0,0,0);
    math::XYZTLorentzVectorD genVisLeg1P4(0,0,0,0);
    math::XYZTLorentzVectorD genVisLeg2P4(0,0,0,0);
  
    double u1 = 0; double u2 = 0;
    int eventDecay=-99;
    unsigned int bosonIndex=0;
    
   
    for(unsigned int k = 0; k < genParticles->size(); k ++){

      if( ! ((*genParticles)[k].pdgId() == 23 || abs((*genParticles)[k].pdgId()) == 24 ||(*genParticles)[k].pdgId() == 25) || 
	  (*genParticles)[k].status()!=3 )
	continue;

      if(verbose_) cout << "Boson status, pt,phi " << (*genParticles)[k].status() << "," << (*genParticles)[k].pt() << "," << (*genParticles)[k].phi() << endl;
      
      eventDecay = (*genParticles)[k].pdgId();
      
      int breakLoop = 0;
      for(unsigned j = 0; j< ((*genParticles)[k].daughterRefVector()).size() && breakLoop!=1; j++){
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 11 ){
	  eventDecay *= 11;
	  breakLoop = 1;
	  bosonIndex = k;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 13 ){
	  eventDecay *= 13;
	  breakLoop = 1;
	  bosonIndex = k;
	}
	if( abs(((*genParticles)[k].daughterRef(j))->pdgId()) == 15 ){
	  eventDecay *= 15;
	  breakLoop = 1;
	  bosonIndex = k;
	}
      }
      if(verbose_) cout << "Decays to pdgId " << eventDecay/(*genParticles)[k].pdgId()  << endl;
      break;
    }
    
    
    
    // Z/H->tautau, tau->tau_h, tau->l
    if(fabs(eventDecay)==23*15 || fabs(eventDecay)==25*15){

      if(verbose_) cout << " ==> Z->tautau event" << endl;
      
      genVP4 = ((*genParticles)[bosonIndex]).p4();

      for(unsigned j = 0; j< ((*genParticles)[bosonIndex].daughterRefVector()).size(); j++){

	reco::GenParticleRef daughter = (*genParticles)[bosonIndex].daughterRef(j);

	int daughterPdgId  = daughter->pdgId() ;
	int daughterStatus = daughter->status() ;

	if( daughterPdgId==-15 && daughterStatus==3){ // tau-

	  genLeg1P4    = daughter->p4();

	  if(verbose_) 
	    cout << " ==> tau- has pt,eta,phi= " << genLeg1P4.Pt() << ", " << genLeg1P4.Eta() << ", " << genLeg1P4.Phi() << endl;

	  genVisLeg1P4 = genLeg1P4;

	  for(unsigned k = 0; k<(daughter->daughterRefVector()).size(); k++){
	    if(daughter->daughterRef(k)->pdgId()==-15){
	      for(unsigned m = 0; m<(daughter->daughterRef(k)->daughterRefVector()).size(); m++){
		int nephewPdgId  = (daughter->daughterRef(k)->daughterRefVector())[m]->pdgId();
		int nephewStatus = (daughter->daughterRef(k)->daughterRefVector())[m]->status();
		if(verbose_) cout << "Nephew " << k << " => pdgId= "<< nephewPdgId 
				  << ", pt= " << (daughter->daughterRef(k)->daughterRefVector())[m]->pt()
				  << endl;
		if( (fabs(nephewPdgId)==12 || fabs(nephewPdgId)==14 || fabs(nephewPdgId)==16) && nephewStatus==1 ) 
		  genVisLeg1P4 -= daughter->daughterRef(k)->daughterRefVector()[m]->p4();
	      }
	    }
	  }

	  if(verbose_) 
	    cout << " ==> tau- has visible pt,eta,phi= " << genVisLeg1P4.Pt() << ", " << genVisLeg1P4.Eta() << ", " << genVisLeg1P4.Phi() << endl;
	}

	else if( daughterPdgId==15 && daughterStatus==3){
	  
	  genLeg2P4 = daughter->p4();
	  
	  if(verbose_) 
	    cout << " ==> tau+ has pt,eta,phi= " << genLeg2P4.Pt() << ", " << genLeg2P4.Eta() << ", " << genLeg2P4.Phi() << endl;
	  
	  genVisLeg2P4 = genLeg2P4;
	  
	  for(unsigned k = 0; k<(daughter->daughterRefVector()).size(); k++){
	    if(daughter->daughterRef(k)->pdgId()==15){
	      for(unsigned m = 0; m<(daughter->daughterRef(k)->daughterRefVector()).size(); m++){
		int nephewPdgId  = (daughter->daughterRef(k)->daughterRefVector())[m]->pdgId();
		int nephewStatus = (daughter->daughterRef(k)->daughterRefVector())[m]->status();
		if(verbose_) cout << "Nephew " << k << " => pdgId= "<< nephewPdgId 
				  << ", pt= " << (daughter->daughterRef(k)->daughterRefVector())[m]->pt()
				  << endl;
		if( (fabs(nephewPdgId)==12 || fabs(nephewPdgId)==14 || fabs(nephewPdgId)==16) && nephewStatus==1 ) 
		  genVisLeg1P4 -= daughter->daughterRef(k)->daughterRefVector()[m]->p4();
	      }
	    }
	  }
	  
	  if(verbose_) 
	    cout << " ==> tau+ has visible pt,eta,phi= " << genVisLeg1P4.Pt() << ", " << genVisLeg1P4.Eta() << ", " << genVisLeg1P4.Phi() << endl;
	}
      }
      
      if(verbose_) cout << "Now let's match to reco objects" << endl;

      bool leg1IsMatched       = false;
      bool leg1IsMatchedByLept = false;
      bool leg2IsMatched       = false;
      bool leg2IsMatchedByLept = false;
      math::XYZTLorentzVectorD recoLeg1P4(0,0,0,0);
      math::XYZTLorentzVectorD recoLeg2P4(0,0,0,0);
      
      for(unsigned int i = 0 ; muons!=0 && !leg1IsMatched && i <muons->size(); i++){
	if( Geom::deltaR((*muons)[i].p4(), genLeg1P4)<0.3 ){
	  leg1IsMatched       = true;
	  leg1IsMatchedByLept = true;
	  recoLeg1P4 = (*muons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; electrons!=0 && !leg1IsMatched && i <electrons->size(); i++){
	if( Geom::deltaR((*electrons)[i].p4(), genLeg1P4)<0.3){
	  leg1IsMatched       = true;
	  leg1IsMatchedByLept = true;
	  recoLeg1P4 = (*electrons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; taus!=0 && !leg1IsMatched && i <taus->size(); i++){
	if( Geom::deltaR((*taus)[i].p4(), genLeg1P4)<0.3){
	  leg1IsMatched       = true;
	  leg1IsMatchedByLept = false;
	  recoLeg1P4 = (*taus)[i].p4();
	}
      }
      
      for(unsigned int i = 0 ; muons!=0 && !leg2IsMatched && i <muons->size(); i++){
	if( Geom::deltaR((*muons)[i].p4(), genLeg2P4)<0.3 ){
	  leg2IsMatched       = true;
	  leg2IsMatchedByLept = true;
	  recoLeg2P4 = (*muons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; electrons!=0 && !leg2IsMatched && i <electrons->size(); i++){
	if( Geom::deltaR((*electrons)[i].p4(), genLeg2P4)<0.3){
	  leg2IsMatched       = true;
	  leg2IsMatchedByLept = true;
	  recoLeg2P4 = (*electrons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; taus!=0 && !leg2IsMatched && i <taus->size(); i++){
	if( Geom::deltaR((*taus)[i].p4(), genLeg2P4)<0.3){
	  leg2IsMatched       = true;
	  leg2IsMatchedByLept = false;
	  recoLeg2P4 = (*taus)[i].p4();
	}
      }
   
      double leptPt  = -99;
      double leptPhi = -99;

      if(leg1IsMatched && leg2IsMatched){
	if(verbose_) cout << "Both gen taus can be matched to reco objects" << endl;
	leptPt  = (recoLeg1P4+recoLeg2P4).Pt();
	leptPhi = (recoLeg1P4+recoLeg2P4).Phi();
      }
      else if(leg1IsMatched && !leg2IsMatched){
	if(verbose_) cout << "Only one gen tau can be matched to reco objects" << endl;
	leptPt  = (recoLeg1P4+genVisLeg2P4).Pt();
	leptPhi = (recoLeg1P4+genVisLeg2P4).Phi();
      }
      else if(!leg1IsMatched && leg2IsMatched){
	if(verbose_) cout << "Only one gen tau can be matched to reco objects" << endl;
	leptPt  = (genVisLeg1P4+recoLeg2P4).Pt();
	leptPhi = (genVisLeg1P4+recoLeg2P4).Phi();
      }
      else{
	if(verbose_) cout << "None of the gen taus can be matched to reco objects" << endl;
	leptPt  = (genVisLeg1P4+genVisLeg2P4).Pt();
	leptPhi = (genVisLeg1P4+genVisLeg2P4).Phi();
      }

      if(verbose_) cout << "Raw MEt " << scaledMETPtN << endl;

      if(verbose_) cout << "Evaluating Nominal recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtN,scaledMETPhiN,genVP4.Pt() ,genVP4.Phi(), 
			       leptPt,leptPhi, 
				u1, u2 , 0 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Nominal   = "    << scaledMETPtN << endl; 
      if(verbose_) cout << "Evaluating +1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtU,scaledMETPhiU,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 
				u1, u2 , +1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Up = " << scaledMETPtU << endl;
      if(verbose_) cout << "Evaluating -1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtD,scaledMETPhiD,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 			     
				u1, u2 , -1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Down = "    << scaledMETPtD << endl;

    }
   
    /////////////////////////////////////////////////////////////////////////
    else if(fabs(eventDecay)==24*11 || fabs(eventDecay)==24*13 || fabs(eventDecay)==24*15){

      if(verbose_) cout << " ==> W event" << endl;

      bool leg1IsMatched       = false;
      bool leg1IsMatchedByLept = false;
      math::XYZTLorentzVectorD recoLeg1P4(0,0,0,0);
      
      genVP4 = ((*genParticles)[bosonIndex]).p4();
      
      for(unsigned j = 0; j< ((*genParticles)[bosonIndex].daughterRefVector()).size(); j++){
	
	reco::GenParticleRef daughter = (*genParticles)[bosonIndex].daughterRef(j);
	
	int daughterPdgId  = daughter->pdgId() ;
	int daughterStatus = daughter->status() ;

	if(fabs(daughterPdgId)==11 && daughterStatus==3){

	  if(verbose_) cout << "W decays to electron + neutrino" << endl;

	  genLeg1P4    = daughter->p4();
	  genVisLeg1P4 = genLeg1P4;

	  for(unsigned int i = 0 ; electrons!=0 && !leg1IsMatched && i <electrons->size(); i++){
	    if( Geom::deltaR((*electrons)[i].p4(), genLeg1P4)<0.3){
	      leg1IsMatched       = true;
	      leg1IsMatchedByLept = true;
	      recoLeg1P4 = (*electrons)[i].p4();
	    }
	  }
	}
	else if(fabs(daughterPdgId)==13 && daughterStatus==3){

	  if(verbose_) cout << "W decays to muon + neutrino" << endl;

	  genLeg1P4    = daughter->p4();
	  genVisLeg1P4 = genLeg1P4;

	  for(unsigned int i = 0 ; muons!=0 && !leg1IsMatched && i <muons->size(); i++){
	    if( Geom::deltaR((*muons)[i].p4(), genLeg1P4)<0.3 ){
	      leg1IsMatched       = true;
	      leg1IsMatchedByLept = true;
	      recoLeg1P4 = (*muons)[i].p4();
	    }
	  }
	}
	else if(fabs(daughterPdgId)==15 && daughterStatus==3){

	  if(verbose_) cout << "W decays to tau + neutrino" << endl;

	  genLeg1P4    = daughter->p4();
	  genVisLeg1P4 = genLeg1P4;

	  for(unsigned k = 0; k<(daughter->daughterRefVector()).size(); k++){
	    if(fabs(daughter->daughterRef(k)->pdgId())==15){
	      for(unsigned m = 0; m<(daughter->daughterRef(k)->daughterRefVector()).size(); m++){
		int nephewPdgId  = (daughter->daughterRef(k)->daughterRefVector())[m]->pdgId();
		int nephewStatus = (daughter->daughterRef(k)->daughterRefVector())[m]->status();
		if(verbose_) cout << "Nephew " << k << " => pdgId= "<< nephewPdgId 
				  << ", pt= " << (daughter->daughterRef(k)->daughterRefVector())[m]->pt()
				  << endl;
		if( (fabs(nephewPdgId)==12 || fabs(nephewPdgId)==14 || fabs(nephewPdgId)==16) && nephewStatus==1 ) 
		  genVisLeg1P4 -= daughter->daughterRef(k)->daughterRefVector()[m]->p4();
	      }
	    }
	  }

	  if(verbose_) 
	    cout << " ==> tau has visible pt,eta,phi= " << genVisLeg1P4.Pt() << ", " << genVisLeg1P4.Eta() << ", " << genVisLeg1P4.Phi() << endl;

	  for(unsigned int i = 0 ; electrons!=0 && !leg1IsMatched && i <electrons->size(); i++){
	    if( Geom::deltaR((*electrons)[i].p4(), genLeg1P4)<0.3){
	      leg1IsMatched       = true;
	      leg1IsMatchedByLept = true;
	      recoLeg1P4 = (*electrons)[i].p4();
	    }
	  }
	  
	  for(unsigned int i = 0 ; muons!=0 && !leg1IsMatched && i <muons->size(); i++){
	    if( Geom::deltaR((*muons)[i].p4(), genLeg1P4)<0.3 ){
	      leg1IsMatched       = true;
	      leg1IsMatchedByLept = true;
	      recoLeg1P4 = (*muons)[i].p4();
	    }
	  }
	}
      }

      double leptPt  = -99;
      double leptPhi = -99;
      
      if(leg1IsMatched){
	if(verbose_) cout << "The lepton from W decay is matched to reco objects" << endl;
	leptPt  = recoLeg1P4.Pt();
	leptPhi = recoLeg1P4.Phi();
      }
      else{
	if(verbose_) cout << "The lepton from W decay is not matched to reco objects" << endl;
	leptPt  = genVisLeg1P4.Pt();
	leptPhi = genVisLeg1P4.Phi();
      }

      if(verbose_) cout << "Raw MEt " << scaledMETPtN << endl;

      if(verbose_) cout << "Evaluating Nominal recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtN,scaledMETPhiN,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 
				u1, u2 , 0 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Nominal   = "    << scaledMETPtN << endl; 
      if(verbose_) cout << "Evaluating +1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtU,scaledMETPhiU,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 
				u1, u2 , +1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Up = " << scaledMETPtU << endl;
      if(verbose_) cout << "Evaluating -1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtD,scaledMETPhiD,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 			     
				u1, u2 , -1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Down = "    << scaledMETPtD << endl;
      
    }

    // Z->ll + fakes
    else if(fabs(eventDecay)==23*11 || fabs(eventDecay)==23*13){

      if(verbose_) cout << " ==> Z->ll + fakes event" << endl;
      
      genVP4 = ((*genParticles)[bosonIndex]).p4();

      for(unsigned j = 0; j< ((*genParticles)[bosonIndex].daughterRefVector()).size(); j++){

	reco::GenParticleRef daughter = (*genParticles)[bosonIndex].daughterRef(j);

	int daughterPdgId  = daughter->pdgId() ;
	int daughterStatus = daughter->status() ;

	if( (daughterPdgId==-11 || daughterPdgId==-13 ) && daughterStatus==3){
	  genLeg1P4    = daughter->p4();
	  if(verbose_) 
	    cout << " ==> lept- has pt,eta,phi= " << genLeg1P4.Pt() << ", " << genLeg1P4.Eta() << ", " << genLeg1P4.Phi() << endl;
	  genVisLeg1P4 = genLeg1P4;
	}
	if( (daughterPdgId==+11 || daughterPdgId==+13 ) && daughterStatus==3){
	  genLeg2P4    = daughter->p4();
	  if(verbose_) 
	    cout << " ==> lept+ has pt,eta,phi= " << genLeg2P4.Pt() << ", " << genLeg2P4.Eta() << ", " << genLeg2P4.Phi() << endl;
	  genVisLeg1P4 = genLeg1P4;
	}
      }


      if(verbose_) cout << "Now let's match to reco objects" << endl;

      bool leg1IsMatched       = false;
      bool leg2IsMatched       = false;
      math::XYZTLorentzVectorD recoLeg1P4(0,0,0,0);
      math::XYZTLorentzVectorD recoLeg2P4(0,0,0,0);
      
      for(unsigned int i = 0 ; muons!=0 && !leg1IsMatched && i <muons->size(); i++){
	if( Geom::deltaR((*muons)[i].p4(), genLeg1P4)<0.3 ){
	  leg1IsMatched       = true;
	  recoLeg1P4 = (*muons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; electrons!=0 && !leg1IsMatched && i <electrons->size(); i++){
	if( Geom::deltaR((*electrons)[i].p4(), genLeg1P4)<0.3){
	  leg1IsMatched       = true;
	  recoLeg1P4 = (*electrons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; taus!=0 && !leg1IsMatched && i <taus->size(); i++){
	if( Geom::deltaR((*taus)[i].p4(), genLeg1P4)<0.3){
	  leg1IsMatched       = true;
	  recoLeg1P4 = (*taus)[i].p4();
	}
      }
      for(unsigned int i = 0 ; muons!=0 && !leg2IsMatched && i <muons->size(); i++){
	if( Geom::deltaR((*muons)[i].p4(), genLeg2P4)<0.3 ){
	  leg2IsMatched       = true;
	  recoLeg2P4 = (*muons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; electrons!=0 && !leg2IsMatched && i <electrons->size(); i++){
	if( Geom::deltaR((*electrons)[i].p4(), genLeg2P4)<0.3){
	  leg2IsMatched       = true;
	  recoLeg2P4 = (*electrons)[i].p4();
	}
      }
      for(unsigned int i = 0 ; taus!=0 && !leg2IsMatched && i <taus->size(); i++){
	if( Geom::deltaR((*taus)[i].p4(), genLeg2P4)<0.3){
	  leg2IsMatched       = true;
	  recoLeg2P4 = (*taus)[i].p4();
	}
      }

      double leptPt  = -99;
      double leptPhi = -99;
      
      if(leg1IsMatched && leg2IsMatched){
	if(verbose_) cout << "Both gen leptons can be matched to reco objects" << endl;
	leptPt  = (recoLeg1P4+recoLeg2P4).Pt();
	leptPhi = (recoLeg1P4+recoLeg2P4).Phi();
      }
      else if(leg1IsMatched && !leg2IsMatched){
	if(verbose_) cout << "Only one gen lepton can be matched to reco objects" << endl;
	leptPt  = (recoLeg1P4+genLeg2P4).Pt();
	leptPhi = (recoLeg1P4+genLeg2P4).Phi();
      }
      else if(!leg1IsMatched && leg2IsMatched){
	if(verbose_) cout << "Only one gen lepton can be matched to reco objects" << endl;
	leptPt  = (genLeg1P4+recoLeg2P4).Pt();
	leptPhi = (genLeg1P4+recoLeg2P4).Phi();
      }
      else{
	if(verbose_) cout << "None of the gen leptons can be matched to reco objects" << endl;
	leptPt  = (genLeg1P4+genLeg2P4).Pt();
	leptPhi = (genLeg1P4+genLeg2P4).Phi();
      }

      if(verbose_) cout << "Raw MEt " << scaledMETPtN << endl;

      if(verbose_) cout << "Evaluating Nominal recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtN,scaledMETPhiN,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 
				u1, u2 , 0 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Nominal   = "    << scaledMETPtN << endl; 
      if(verbose_) cout << "Evaluating +1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtU,scaledMETPhiU,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 
				u1, u2 , +1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Up = " << scaledMETPtU << endl;
      if(verbose_) cout << "Evaluating -1 recoil-corrected MEt" << endl;
      recoilCorr_->CorrectType1(scaledMETPtD,scaledMETPhiD,genVP4.Pt() ,genVP4.Phi(), 
				leptPt,leptPhi, 			     
				u1, u2 , -1 , TMath::Min(nJets30,2) );
      if(verbose_) cout << "corrected MET Down = "    << scaledMETPtD << endl;

    }

  }
  


  TLorentzVector LVscaledMETN;
  LVscaledMETN.SetPtEtaPhiM(scaledMETPtN,0,scaledMETPhiN,0);
  scaledMETPxN = LVscaledMETN.Px();
  scaledMETPyN = LVscaledMETN.Py();

  TLorentzVector LVscaledMETU;
  LVscaledMETU.SetPtEtaPhiM(scaledMETPtU,0,scaledMETPhiU,0);
  scaledMETPxU = LVscaledMETU.Px();
  scaledMETPyU = LVscaledMETU.Py();

  TLorentzVector LVscaledMETD;
  LVscaledMETD.SetPtEtaPhiM(scaledMETPtD,0,scaledMETPhiD,0);
  scaledMETPxD = LVscaledMETD.Px();
  scaledMETPyD = LVscaledMETD.Py();


  pat::MET scaledMETU( reco::MET(((*mets)[0]).sumEt(), 
				 reco::MET::LorentzVector(scaledMETPxU, scaledMETPyU, 0, sqrt(scaledMETPxU*scaledMETPxU+scaledMETPyU*scaledMETPyU)), 
				 reco::MET::Point(0,0,0)) );
  pat::MET scaledMETN( reco::MET(((*mets)[0]).sumEt(), 
				 reco::MET::LorentzVector(scaledMETPxN, scaledMETPyN, 0, sqrt(scaledMETPxN*scaledMETPxN+scaledMETPyN*scaledMETPyN)), 
				 reco::MET::Point(0,0,0)) );
  pat::MET scaledMETD( reco::MET(((*mets)[0]).sumEt(), 
				 reco::MET::LorentzVector(scaledMETPxD, scaledMETPyD, 0, sqrt(scaledMETPxD*scaledMETPxD+scaledMETPyD*scaledMETPyD)), 
				 reco::MET::Point(0,0,0)) );

  MEtRescaledCollU->push_back(scaledMETU);  
  MEtRescaledCollN->push_back(scaledMETN);  
  MEtRescaledCollD->push_back(scaledMETD);  

  iEvent.put(MEtRescaledCollU, "U");
  iEvent.put(MEtRescaledCollN, "N");
  iEvent.put(MEtRescaledCollD, "D");
  
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MEtRecoilCorrectorProducer);


