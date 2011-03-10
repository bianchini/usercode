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

    if(vertexes->size()!=0 && aMuon.isGlobalMuon())  
      dxyWrtPV = (aMuon.globalTrack())->dxy( (*vertexes)[0].position() ) ;
    else if (vertexes->size()!=0 && aMuon.isTrackerMuon())
      dxyWrtPV = (aMuon.innerTrack())->dxy( (*vertexes)[0].position() ) ;

    aMuon.addUserFloat("dxyWrtPV",dxyWrtPV);

    float chIso = 
      aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
    float nhIso = 
      aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;
    float phIso = 
      aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4,reco::IsoDeposit::Vetos(),false ).first;

    aMuon.addUserFloat("PFRelIso03",(chIso+nhIso+phIso)/aMuon.pt());

    aMuon.addUserFloat("isInRun",iEvent.run());

    muonsUserEmbeddedColl->push_back(aMuon);

  }


  iEvent.put( muonsUserEmbeddedColl );
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
 

DEFINE_FWK_MODULE(MuonsUserEmbedded);


