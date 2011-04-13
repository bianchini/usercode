#include "Bianchi/TauTauStudies/interface/ElectronsUserEmbedded.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

ElectronsUserEmbedded::ElectronsUserEmbedded(const edm::ParameterSet & iConfig){

  electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
  vertexTag_   = iConfig.getParameter<edm::InputTag>("vertexTag");
 
  produces<pat::ElectronCollection>("");

}

ElectronsUserEmbedded::~ElectronsUserEmbedded(){
}

void ElectronsUserEmbedded::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByLabel(electronTag_,electronsHandle);
  const pat::ElectronCollection* electrons = electronsHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  const reco::VertexCollection* vertexes = vertexHandle.product();

  std::auto_ptr< pat::ElectronCollection > electronsUserEmbeddedColl( new pat::ElectronCollection() ) ;

  for(unsigned int i = 0; i < electrons->size(); i++){
    pat::Electron aElectron( (*electrons)[i] );
   
    double dxyWrtPV =  -99.;

    if(vertexes->size()!=0 && (aElectron.gsfTrack()).isNonnull() )  
      dxyWrtPV = (aElectron.gsfTrack())->dxy( (*vertexes)[0].position() ) ;
    else if (vertexes->size()!=0 && (aElectron.track()).isNonnull() )
      dxyWrtPV = (aElectron.track())->dxy( (*vertexes)[0].position() ) ;

    aElectron.addUserFloat("dxyWrtPV",dxyWrtPV);

    reco::isodeposit::AbsVetos vetosCharged; 
    reco::isodeposit::AbsVetos vetosNeutral;  
    reco::isodeposit::AbsVetos vetosPhotons;
    vetosCharged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
    vetosNeutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(aElectron.eta(),aElectron.phi()),0.08));
    vetosNeutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));
    vetosPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(aElectron.eta(),aElectron.phi()),0.05));
    vetosPhotons.push_back(new reco::isodeposit::ThresholdVeto(1.0));

    float chIso03 = 
      aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosCharged).first;
    float nhIso03 = 
      aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosNeutral).first;
    float phIso03 = 
      aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosPhotons).first;
    float chIso04 = 
      aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosCharged).first;
    float nhIso04 = 
      aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosNeutral).first;
    float phIso04 = 
      aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosPhotons).first;

    aElectron.addUserFloat("PFRelIso04",(chIso04+nhIso04+phIso04)/aElectron.pt());
    aElectron.addUserFloat("PFRelIso03",(chIso03+nhIso03+phIso03)/aElectron.pt());

    aElectron.addUserFloat("isInRun",iEvent.run());

    electronsUserEmbeddedColl->push_back(aElectron);
    
    // cleaning
    for(unsigned int i = 0; i <vetosCharged.size(); i++){
      delete vetosCharged[i];
    }
    for(unsigned int i = 0; i <vetosNeutral.size(); i++){
      delete vetosNeutral[i];
      delete vetosPhotons[i];
    }

  }


  iEvent.put( electronsUserEmbeddedColl );
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronsUserEmbedded);


