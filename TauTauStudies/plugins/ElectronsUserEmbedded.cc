#include "Bianchi/TauTauStudies/interface/ElectronsUserEmbedded.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

using namespace edm;
using namespace std;
using namespace reco;


ElectronsUserEmbedded::ElectronsUserEmbedded(const edm::ParameterSet & iConfig){

  electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
  vertexTag_   = iConfig.getParameter<edm::InputTag>("vertexTag");
  isMC_        = iConfig.getParameter<bool>("isMC");
 
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

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel("scalersRawToDigi", dcsHandle);
  float evt_bField;

  if (!isMC_) {
    // scale factor = 3.801/18166.0 which are
    // average values taken over a stable two
    // week period
    float currentToBFieldScaleFactor = 2.09237036221512717e-04;
    float current = -9999/currentToBFieldScaleFactor;
    if( dcsHandle.isValid() && (*dcsHandle).size() > 0 ) {
      current = (*dcsHandle)[0].magnetCurrent();
    }
      
    evt_bField = current*currentToBFieldScaleFactor;
  }
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
  }

  //Get the CTF tracks
  Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
  
  //get GSF Tracks
  Handle<reco::GsfTrackCollection> gsftracks_h;
  iEvent.getByLabel("electronGsfTracks", gsftracks_h);
  
  std::auto_ptr< pat::ElectronCollection > electronsUserEmbeddedColl( new pat::ElectronCollection() ) ;

  for(unsigned int i = 0; i < electrons->size(); i++){

    pat::Electron aElectron( (*electrons)[i] );
    const reco::GsfElectron* aGsf = static_cast<reco::GsfElectron*>(&aElectron); 

    const reco::Track *el_track = (const reco::Track*)((aElectron).gsfTrack().get());  
    const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
    float nHits = p_inner.numberOfHits();

    ConversionFinder convFinder;
    vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(aElectron.core()), tracks_h, gsftracks_h, evt_bField);
    ConversionInfo convInfo  = convFinder.getConversionInfo(*aGsf, tracks_h, gsftracks_h, evt_bField);
    double els_conv_dist     = convInfo.dist();
    double els_conv_dcot     = convInfo.dcot();
    double els_conv_radius   = convInfo.radiusOfConversion();
    math::XYZPoint els_conv_Point = convInfo.pointOfConversion(); 
    TrackRef els_conv_ctfRef = convInfo.conversionPartnerCtfTk(); 
    GsfTrackRef els_conv_gsfRef = convInfo.conversionPartnerGsfTk();
    double els_conv_delMissHits =  convInfo.deltaMissingHits();

    float dPhi  = aElectron.deltaPhiSuperClusterTrackAtVtx();
    float dEta  = aElectron.deltaEtaSuperClusterTrackAtVtx();
    float sihih = aElectron.sigmaIetaIeta();
    float HoE   = aElectron.hadronicOverEm();
    //cout << "dEta " << dEta << " dPhi " << dPhi << " -- dcot " << els_conv_dcot << " -- nHits " << nHits << endl;

    aElectron.addUserFloat("nHits",nHits);
    aElectron.addUserFloat("dist",abs(els_conv_dist));
    aElectron.addUserFloat("dcot",abs(els_conv_dcot));
    aElectron.addUserFloat("dPhi",abs(dPhi));
    aElectron.addUserFloat("dEta",abs(dEta));
    aElectron.addUserFloat("sihih",sihih);
    aElectron.addUserFloat("HoE",HoE);

    double dxyWrtPV =  -99.;
    double dzWrtPV =  -99.;

    if(vertexes->size()!=0 && (aElectron.gsfTrack()).isNonnull() ){
      dxyWrtPV = (aElectron.gsfTrack())->dxy( (*vertexes)[0].position() ) ;
      dzWrtPV  = (aElectron.gsfTrack())->dz( (*vertexes)[0].position() ) ;
    }
    else if (vertexes->size()!=0 && (aElectron.track()).isNonnull() ){
      dxyWrtPV = (aElectron.track())->dxy( (*vertexes)[0].position() ) ;
      dzWrtPV  = (aElectron.track())->dz( (*vertexes)[0].position() ) ;
    }

    aElectron.addUserFloat("dxyWrtPV",dxyWrtPV);
    aElectron.addUserFloat("dzWrtPV",dzWrtPV);

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


