#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TLorentzVector.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include <vector>
#include <utility>
#include <map>

using namespace edm;
using namespace std;

class UserDefinedVariables : public edm::EDProducer {
public:
  explicit UserDefinedVariables(const edm::ParameterSet & iConfig);
  virtual ~UserDefinedVariables() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  virtual void beginJob();

private:

  typedef std::vector<double> vdouble;

  edm::InputTag objects_;
  edm::InputTag objects2_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag met_;            
  bool isMC_;
  edm::LumiReWeighting LumiWeights2011A_;
  edm::LumiReWeighting LumiWeights2011B_;
  vdouble TrueDist2011A_f_;
  vdouble TrueDist2011B_f_;
  vdouble MCDist_f_;
};

UserDefinedVariables::UserDefinedVariables(const edm::ParameterSet & iConfig) :
  objects_(iConfig.getParameter<edm::InputTag>("objects")),
  objects2_(iConfig.getParameter<edm::InputTag>("objects2")),
  triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")), 
  met_(iConfig.getParameter<edm::InputTag>("met")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  TrueDist2011A_f_(iConfig.getParameter<vdouble>("TrueDist2011A")),
  TrueDist2011B_f_(iConfig.getParameter<vdouble>("TrueDist2011B")),
  MCDist_f_(iConfig.getParameter<vdouble>("MCDist"))
{
  produces<edm::ValueMap<float> >("Mt");
  produces<edm::ValueMap<float> >("eleIdWP70");
  produces<edm::ValueMap<float> >("numEleWP80");
  produces<edm::ValueMap<float> >("numEleWP95");
  produces<edm::ValueMap<float> >("puMCWeightRun2011A");
  produces<edm::ValueMap<float> >("puMCWeightRun2011B");
  produces<edm::ValueMap<float> >("triggerBit");
  produces<edm::ValueMap<float> >("tauXTriggersLeg1");
  produces<edm::ValueMap<float> >("tauXTriggersLeg2");
  produces<edm::ValueMap<float> >("triggerBitSingleEle");
  produces<edm::ValueMap<float> >("eleXTriggersLeg1");
}


UserDefinedVariables::~UserDefinedVariables()
{

}

void UserDefinedVariables::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

    // read input
    Handle<View<pat::Electron> > objects;
    iEvent.getByLabel(objects_, objects);

    Handle<View<reco::Candidate> > objects2;
    iEvent.getByLabel(objects2_, objects2);

    edm::Handle<pat::METCollection> metHandle;
    iEvent.getByLabel(met_,metHandle);
    const pat::METCollection* met = metHandle.product();

    int nPUVertices = -99;
    //int nOOTPUVertices = -99;
    float mcPUweight2011A = 1;
    float mcPUweight2011B = 1;
    if(isMC_){

	edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

	std::vector<PileupSummaryInfo>::const_iterator PVI;

	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

	  	nPUVertices = PVI->getTrueNumInteractions();
	}

    	mcPUweight2011A = LumiWeights2011A_.weight(nPUVertices);
    	mcPUweight2011B = LumiWeights2011B_.weight(nPUVertices);

    }

    //cout << mcPUweight << " -- " << nPUVertices << endl;

    // prepare vector for output   
    std::vector<float> values;
    std::vector<float> values2;
    std::vector<float> values3;
    std::vector<float> values4;

    float numWP80 = 0;
    float numWP95 = 0;
    std::vector<float> valuesWP80;
    std::vector<float> valuesWP95;

    View<pat::Electron>::const_iterator object; 
    for (object = objects->begin(); object != objects->end(); ++object) {

      float scalarSumPt = (object->p4()).Pt() + ((*met)[0].p4()).Pt();
      float vectorSumPt = (object->p4() + (*met)[0].p4()).Pt() ;
      float Mt = TMath::Sqrt( scalarSumPt*scalarSumPt - vectorSumPt*vectorSumPt ) ;

      values.push_back(Mt);
      values2.push_back(mcPUweight2011A);
      values3.push_back(mcPUweight2011B);

      float simpleCutsWP70 = (object->userFloat("nHits")==0 && object->userInt("antiConv")>0.5 &&

                               ((object->isEB() && object->userFloat("sihih")<0.010 && object->userFloat("dPhi")<0.03 &&
                                                   object->userFloat("dEta")< 0.004 && object->userFloat("HoE") <0.025) ||
                                (object->isEE() && object->userFloat("sihih")<0.030 && object->userFloat("dPhi")<0.020 &&
                                                   object->userFloat("dEta") <0.005 && object->userFloat("HoE") <0.025))
			     );

      float simpleCutsWP80 = (object->userFloat("nHits")==0 && object->userInt("antiConv")>0.5 &&

                               ((object->isEB() && object->userFloat("sihih")<0.010 && object->userFloat("dPhi")<0.06 &&
                                                   object->userFloat("dEta")< 0.004 && object->userFloat("HoE") <0.04) ||
                                (object->isEE() && object->userFloat("sihih")<0.030 && object->userFloat("dPhi")<0.030 &&
                                                   object->userFloat("dEta") <0.007 && object->userFloat("HoE") <0.025))
			     );

      float simpleCutsWP95 = (object->userFloat("nHits")<=1 &&

                               ((object->isEB() && object->userFloat("sihih")<0.010 && object->userFloat("dPhi")<0.80 &&
                                                   object->userFloat("dEta")< 0.007 && object->userFloat("HoE") <0.15) ||
                                (object->isEE() && object->userFloat("sihih")<0.030 && object->userFloat("dPhi")<0.70 &&
                                                   object->userFloat("dEta") <0.010 && object->userFloat("HoE") <0.07))
			     );

      if(simpleCutsWP80) ++numWP80;
      if(simpleCutsWP95) ++numWP95;

      values4.push_back(simpleCutsWP70);

    }

    //std::cout<<numWP80<<" "<<numWP95<<std::endl;

    int sizeEle = objects->size();
    for(int i=0; i<sizeEle; i++){

	valuesWP80.push_back(numWP80);
	valuesWP95.push_back(numWP95);

    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(objects, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, "Mt");

    std::auto_ptr<ValueMap<float> > valMap9(new ValueMap<float>());
    ValueMap<float>::Filler filler9(*valMap9);
    filler9.insert(objects, values4.begin(), values4.end());
    filler9.fill();
    iEvent.put(valMap9, "eleIdWP70");

    std::auto_ptr<ValueMap<float> > valMap2(new ValueMap<float>());
    ValueMap<float>::Filler filler2(*valMap2);
    filler2.insert(objects, values2.begin(), values2.end());
    filler2.fill();
    iEvent.put(valMap2, "puMCWeightRun2011A");

    std::auto_ptr<ValueMap<float> > valMap3(new ValueMap<float>());
    ValueMap<float>::Filler filler3(*valMap3);
    filler3.insert(objects, values3.begin(), values3.end());
    filler3.fill();
    iEvent.put(valMap3, "puMCWeightRun2011B");

    std::auto_ptr<ValueMap<float> > valMapWP80(new ValueMap<float>());
    ValueMap<float>::Filler fillerWP80(*valMapWP80);
    fillerWP80.insert(objects, valuesWP80.begin(), valuesWP80.end());
    fillerWP80.fill();
    iEvent.put(valMapWP80, "numEleWP80");

    std::auto_ptr<ValueMap<float> > valMapWP95(new ValueMap<float>());
    ValueMap<float>::Filler fillerWP95(*valMapWP95);
    fillerWP95.insert(objects, valuesWP95.begin(), valuesWP95.end());
    fillerWP95.fill();
    iEvent.put(valMapWP95, "numEleWP95");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //vector<string> XtriggerPaths;
    string triggerPaths;
    string HLTfiltersElec;
    string HLTfiltersTau;

    std::vector<float> tauXTriggersLeg1_;
    std::vector<float> tauXTriggersLeg2_;
    std::vector<float> triggerBits_;

    string triggerPathsSingleEle;
    string HLTfiltersElecSingleEle;

    std::vector<float> eleXTriggersLeg1_;
    std::vector<float> triggerBitsSingleEle_;

    int runNumber = iEvent.id().run();

    if(isMC_){

    	// X-triggers Non esistono nel fall11
    	//XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*");
    	//XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");
    	// for Fall11
    	//triggerPaths.push_back("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");
    	//triggerPaths.push_back("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1");
	triggerPaths = "HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1";

    	//HLTfiltersElec.push_back("hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    	//HLTfiltersElec.push_back("hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle18MediumIsoPFTau20");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle20MediumIsoPFTau20");
	HLTfiltersElec = "hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
	HLTfiltersTau = "hltOverlapFilterIsoEle18MediumIsoPFTau20";

        triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v8";
    	HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";

    }

    else{
    
    	// X-triggers
    	//XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*");
    	//XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");
    	//XtriggerPaths.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v*");
    	//XtriggerPaths.push_back("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*");

    	// Single Electron triggers + X-triggers

	if(runNumber >= 160404 && runNumber <= 161176){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau15";
	}
	else if(runNumber >= 161216 && runNumber <= 163261){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau15";
	}
	else if(runNumber >= 163269 && runNumber <= 163869){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau15";
	}
	else if(runNumber >= 165088 && runNumber <= 165633){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau20";
	}
	else if(runNumber >= 165970 && runNumber <= 166967){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau20";
	}
	else if(runNumber >= 167039 && runNumber <= 167913){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v9";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15IsoPFTau20";
	}
	else if(runNumber >= 170249 && runNumber <= 173198){
    		triggerPaths = "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2";
    		HLTfiltersElec = "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle15TightIsoPFTau20";
	}
	else if(runNumber >= 173236 && runNumber <= 178380){
    		triggerPaths = "HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1";
    		HLTfiltersElec = "hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter";
    		HLTfiltersTau = "hltOverlapFilterIsoEle18MediumIsoPFTau20";
	}
	else if(runNumber >= 178420  && runNumber <= 179889){
    		triggerPaths = "HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v5";
    		HLTfiltersElec = "hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20";
    		HLTfiltersTau = "hltOverlapFilterIsoEle20MediumIsoPFTau20";
	}
	else if(runNumber >= 179959 && runNumber <= 180252){
    		triggerPaths = "HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6";
    		HLTfiltersElec = "hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20";
    		HLTfiltersTau = "hltOverlapFilterIsoEle20MediumIsoPFTau20";
	}
                              
    	// watch out! name convention changed with time
    	//HLTfiltersElec.push_back("hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter"); //Non esiste
    	//HLTfiltersElec.push_back("hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
    	// watch out! name convention changed with time
    	//HLTfiltersElec.push_back("hltEle18CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter"); //Non esiste
    	//HLTfiltersElec.push_back("hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
 
    	//HLTfiltersElec.push_back("hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20");
 
   
    	/////HLTfiltersTau.push_back("hltPFTau15TrackLooseIso"); //Boh saranno esistiti in qualche momento della storia umana
    	/////HLTfiltersTau.push_back("hltPFTau20TrackLooseIso"); //Boh saranno esistiti in qualche momento della storia umana
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau15");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle15IsoPFTau20");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle15TightIsoPFTau20");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle18MediumIsoPFTau20");
    	//HLTfiltersTau.push_back("hltOverlapFilterIsoEle20MediumIsoPFTau20");

	if(runNumber >= 160404 && runNumber <= 160877){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 160888 && runNumber <= 163261){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 163269 && runNumber <= 163593){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 165364 && runNumber <= 165208){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 165970 && runNumber <= 166346){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v5";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 167039 && runNumber <= 167913){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v6";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 170249 && runNumber <= 173198){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v7";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 173236 && runNumber <= 175875){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v8";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 178420  && runNumber <= 179889){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v9";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}
	else if(runNumber >= 179959 && runNumber <= 180252){
        	triggerPathsSingleEle = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v10";
    		HLTfiltersElecSingleEle = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter";
	}

    }

    //std::cout<<"triggerPaths: "<<triggerPaths.c_str()<<" HLTfiltersElec: "<<HLTfiltersElec.c_str()<<" HLTfiltersTau: "<<HLTfiltersTau.c_str()<<std::endl;
    //std::cout<<"triggerPathsSingleEle: "<<triggerPathsSingleEle.c_str()<<" HLTfiltersElecSingleEle: "<<HLTfiltersElec.c_str()<<std::endl;

    edm::Handle<pat::TriggerEvent> triggerHandle;
    iEvent.getByLabel(triggerResultsTag_, triggerHandle);
    const pat::TriggerEvent* trigger = triggerHandle.product();

    edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjsHandle;
    iEvent.getByLabel(edm::InputTag("patTrigger"),triggerObjsHandle);
    const pat::TriggerObjectStandAloneCollection* triggerObjs = triggerObjsHandle.product();
    

    for (object = objects->begin(); object != objects->end(); ++object) {

	if(trigger){

		const pat::TriggerPath *triggerPath =  trigger->path(triggerPaths);
		const pat::TriggerPath *triggerPathSingleEle =  trigger->path(triggerPathsSingleEle);

	    	/*if(verbose_){
	      		cout<<  "Testing " << triggerPaths[i] << endl;
	      		if(triggerPath) cout << "Is there..." << endl;
	      		if(triggerPath && triggerPath->wasRun()) cout << "Was run..." << endl;
	      		if(triggerPath && triggerPath->wasRun() && triggerPath->wasAccept()) cout << "Was accepted..." << endl;
	    	}*/
	    
	   	if(triggerPath && triggerPath->wasRun() && triggerPath->wasAccept() && triggerPath->prescale() ==1 ) triggerBits_.push_back(1);
	    	else if (triggerPath && triggerPath->wasRun() && triggerPath->wasAccept() && triggerPath->prescale()!=1) triggerBits_.push_back(2);
	   	else triggerBits_.push_back(0);

	   	if(triggerPathSingleEle && triggerPathSingleEle->wasRun() && triggerPathSingleEle->wasAccept() && triggerPathSingleEle->prescale() ==1 ) 
			triggerBitsSingleEle_.push_back(1);
	    	else if (triggerPathSingleEle && triggerPathSingleEle->wasRun() && triggerPathSingleEle->wasAccept() && triggerPathSingleEle->prescale()!=1) 
			triggerBitsSingleEle_.push_back(2);
	   	else triggerBitsSingleEle_.push_back(0);

	}

	bool matched = false;
	bool matchedSingleEle = false;
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));

	      	/*if(verbose_) {
			if( Geom::deltaR( aObj->triggerObject().p4(), leg1->p4() )<0.3 ){
		  		for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
		    			cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of electron" << endl;
		  		}
			}
	      	}*/

	      	if( Geom::deltaR( aObj->triggerObject().p4(), object->p4() ) < 0.3 && aObj->hasFilterLabel(HLTfiltersElec) ){
			matched = true;
	      	}
	      	if( Geom::deltaR( aObj->triggerObject().p4(), object->p4() ) < 0.3 && aObj->hasFilterLabel(HLTfiltersElecSingleEle) ){
			matchedSingleEle = true;
	      	}
	}
	if(matched) tauXTriggersLeg1_.push_back(1);
	else tauXTriggersLeg1_.push_back(0);
	if(matchedSingleEle) eleXTriggersLeg1_.push_back(1);
	else eleXTriggersLeg1_.push_back(0);
	/*if(verbose_){
	      	if(matched) cout << "Electron matched within dR=0.3 with trigger object passing filter " << HLTfiltersElec[i] << endl;
	      	else cout << "!!! Electron is not trigger matched within dR=0.3 !!!" << endl;
	}*/

    }

    View<reco::Candidate>::const_iterator object2; 
    for (object2 = objects2->begin(); object2 != objects2->end(); ++object2) {

	bool matched = false;
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));

	      	/*if(verbose_) {
			if( Geom::deltaR( aObj->triggerObject().p4(), leg2->p4() )<0.3 ){
		  		for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
		    			cout << "Object passing " << (aObj->filterLabels())[k] << " within 0.3 of tau" << endl;
		  		}
			}
	      	}*/

	      	if( Geom::deltaR( aObj->triggerObject().p4(), object2->p4() ) < 0.3  && aObj->hasFilterLabel(HLTfiltersTau) ){
			matched = true;
	      	}
	}

	if(matched) tauXTriggersLeg2_.push_back(1);
	else tauXTriggersLeg2_.push_back(0);
	/*if(verbose_){
	      	if(matched) cout << "Tau matched within dR=0.3 with trigger object passing filter " << HLTfiltersTau[i] << endl;
	      	else cout << "!!! Tau is not trigger matched within dR=0.3 !!!" << endl;
	}*/

    }

    //std::cout<<"O1 Size: "<<objects->size()<<" triggerBits: "<<triggerBits_.size()<<" tauXTriggersLeg1: "<<tauXTriggersLeg1_.size()<<std::endl;
    //std::cout<<"O2 Size: "<<objects2->size()<<" tauXTriggersLeg2: "<<tauXTriggersLeg2_.size()<<std::endl;

    std::auto_ptr<ValueMap<float> > valMap6(new ValueMap<float>());
    ValueMap<float>::Filler filler6(*valMap6);
    filler6.insert(objects, triggerBits_.begin(), triggerBits_.end());
    filler6.fill();
    iEvent.put(valMap6, "triggerBit");

    std::auto_ptr<ValueMap<float> > valMap4(new ValueMap<float>());
    ValueMap<float>::Filler filler4(*valMap4);
    filler4.insert(objects, tauXTriggersLeg1_.begin(), tauXTriggersLeg1_.end());
    filler4.fill();
    iEvent.put(valMap4, "tauXTriggersLeg1");

    std::auto_ptr<ValueMap<float> > valMap5(new ValueMap<float>());
    ValueMap<float>::Filler filler5(*valMap5);
    filler5.insert(objects2, tauXTriggersLeg2_.begin(), tauXTriggersLeg2_.end());
    filler5.fill();
    iEvent.put(valMap5, "tauXTriggersLeg2");

    std::auto_ptr<ValueMap<float> > valMap7(new ValueMap<float>());
    ValueMap<float>::Filler filler7(*valMap7);
    filler7.insert(objects, triggerBitsSingleEle_.begin(), triggerBitsSingleEle_.end());
    filler7.fill();
    iEvent.put(valMap7, "triggerBitSingleEle");

    std::auto_ptr<ValueMap<float> > valMap8(new ValueMap<float>());
    ValueMap<float>::Filler filler8(*valMap8);
    filler8.insert(objects, eleXTriggersLeg1_.begin(), eleXTriggersLeg1_.end());
    filler8.fill();
    iEvent.put(valMap8, "eleXTriggersLeg1");

}

void
UserDefinedVariables::beginJob()
{

  if(isMC_){

	  std::vector< float > MCDist ;
	  std::vector< float > TrueDist2011A;
	  std::vector< float > TrueDist2011B;

	  int sizeMCDist_f_ = MCDist_f_.size();

	  for( int i=0; i<sizeMCDist_f_; ++i) {
	      TrueDist2011A.push_back(TrueDist2011A_f_[i]);
	      MCDist.push_back(MCDist_f_[i]);
	  }

	  for( int i=0; i<sizeMCDist_f_; ++i) {
	      TrueDist2011B.push_back(TrueDist2011B_f_[i]);
	  }

	  //std::cout<<MCDist_f_.size()<<" "<<TrueDist2011A_f_.size()<<" "<<TrueDist2011B_f_.size()<<std::endl;
	  //std::cout<<MCDist.size()<<" "<<TrueDist2011A.size()<<" "<<TrueDist2011B.size()<<std::endl;

	  LumiWeights2011A_ = edm::LumiReWeighting(MCDist, TrueDist2011A);
	  LumiWeights2011B_ = edm::LumiReWeighting(MCDist, TrueDist2011B);

  }

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(UserDefinedVariables);

