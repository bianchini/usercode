#include "Bianchi/TauTauStudies/interface/MuonsUserEmbedded.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

MuonsUserEmbedded::MuonsUserEmbedded(const edm::ParameterSet & iConfig){

  muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
 
  produces<pat::MuonCollection>("");

}

MuonsUserEmbedded::~MuonsUserEmbedded(){
}

void MuonsUserEmbedded::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel(muonTag_,muonsHandle);
  const pat::MuonCollection* muons = muonsHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  const reco::VertexCollection* vertexes = vertexHandle.product();

  std::auto_ptr< pat::MuonCollection > muonsUserEmbeddedColl( new pat::MuonCollection() ) ;

  for(unsigned int i = 0; i < muons->size(); i++){
    pat::Muon aMuon( (*muons)[i] );
   
    double dxyWrtPV =  -99.;
    double dzWrtPV  =  -99.;

    if(vertexes->size()!=0 && aMuon.isGlobalMuon()){
      dxyWrtPV = (aMuon.globalTrack())->dxy( (*vertexes)[0].position() ) ;
      dzWrtPV  = (aMuon.globalTrack())->dxy( (*vertexes)[0].position() ) ;
    }
    else if (vertexes->size()!=0 && aMuon.isTrackerMuon()){
      dxyWrtPV = (aMuon.innerTrack())->dxy( (*vertexes)[0].position() ) ;
      dzWrtPV  = (aMuon.innerTrack())->dxy( (*vertexes)[0].position() ) ;
    }

    aMuon.addUserFloat("dxyWrtPV",dxyWrtPV);
    aMuon.addUserFloat("dzWrtPV",dzWrtPV);

    reco::isodeposit::AbsVetos vetosCharged; 
    reco::isodeposit::AbsVetos vetosNeutral;  
    reco::isodeposit::AbsVetos vetosPhotons;
    vetosCharged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
    vetosNeutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(aMuon.eta(),aMuon.phi()),0.08));
    vetosNeutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));
    vetosPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(aMuon.eta(),aMuon.phi()),0.05));
    vetosPhotons.push_back(new reco::isodeposit::ThresholdVeto(1.0));

    float chIso03 = 
      aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosCharged).first;
    float nhIso03 = 
      aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosNeutral).first;
    float phIso03 = 
      aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosPhotons).first;
    float chIso04 = 
      aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosCharged).first;
    float nhIso04 = 
      aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosNeutral).first;
    float phIso04 = 
      aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosPhotons).first;

    aMuon.addUserFloat("PFRelIso04",(chIso04+nhIso04+phIso04)/aMuon.pt());
    aMuon.addUserFloat("PFRelIso03",(chIso03+nhIso03+phIso03)/aMuon.pt());

    aMuon.addUserFloat("isInRun",iEvent.run());

    muonsUserEmbeddedColl->push_back(aMuon);
    
    // cleaning
    for(unsigned int i = 0; i <vetosCharged.size(); i++){
      delete vetosCharged[i];
    }
    for(unsigned int i = 0; i <vetosNeutral.size(); i++){
      delete vetosNeutral[i];
      delete vetosPhotons[i];
    }

  }


  iEvent.put( muonsUserEmbeddedColl );
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
 

DEFINE_FWK_MODULE(MuonsUserEmbedded);


