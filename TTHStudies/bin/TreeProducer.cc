#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "Bianchi/TTHStudies/interface/Samples.h"

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


typedef struct
{
  float bmass;
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass;
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass;
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;
} genTopInfo;

typedef struct
{
  int index;
  float pt;
  float eta;
  float phi;
  float mass;
  float csv;
  float topB;
  float topW;
  float atopB;
  float atopW;
  float higgsB;
  float flavor;
  float unc;
  void reset(){
    index = -99; pt = -99; eta = -99; phi = -99; mass = -99; csv = -99; topB = -99; topW = -99;  atopW = -99;  atopB = -99;
    higgsB = -99; flavor = -99; unc = -99;
  }
  void set(int index_,  float pt_, float eta_, float phi_, float mass_,float csv_, float topB_, float topW_, float atopB_,
	   float atopW_, float higgsB_, float flavor_, float unc_){
    index = index_; 
    pt = pt_; 
    eta = eta_; 
    phi = phi_; 
    mass = mass_; 
    csv = csv_; 
    topB = topB_; 
    topW = topW_;  
    atopW = atopW_;  
    atopB = atopB_;
    higgsB = higgsB_; 
    flavor = flavor_; 
    unc = unc_;
  }
} JetByPt;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;



Float_t deltaR(TLorentzVector j1, TLorentzVector j2){
  float res = TMath::Sqrt(TMath::Power(j1.Eta()-j2.Eta(),2) + TMath::Power( TMath::ACos(TMath::Cos(j1.Phi()-j2.Phi())) ,2));
  return res;
}

Int_t match(TLorentzVector j1, TLorentzVector j2){
  return deltaR(j1,j2) < GENJETDR;
}



int main(int argc, const char* argv[])
{

  std::cout << "TreeProducer" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  fwlite::TFileService fs = fwlite::TFileService("./root/treeProducer.root");
  TTree* genTree          = fs.make<TTree>("genTree","event tree");
  TTree* genJetLightTree  = fs.make<TTree>("genJetLightTree","event tree");
  TTree* genEventTree     = fs.make<TTree>("genEventTree","event tree");
  TTree* genJetHeavyTree  = fs.make<TTree>("genJetHeavyTree","event tree");

  float X1;
  float X2;
  float PtTTH  ;
  float PzTTH  ;
  float PhiTTH  ;
  float BetaTTH  ;
  float GammaTTH    ;
  float DeltaTTH    ;
  float PtHStar ;
  float DphiHStar ;
  float GammaTT;
  float MassTT;
  float BetaW  ;
  float GammaW  ;
  float DeltaW  ;
  float BetaWLep  ;
  float GammaWLep  ;
  float DeltaWLep  ;
  float PxT, PyT, PzT;
  float PxH, PyH, PzH;
  //float M12;

  float ptRecoLight;
  float phiRecoLight;
  float etaRecoLight;
  float csvRecoLight;
  float ptLight;
  float genPtLight;
  float etaLight;
  float phiLight;
  float massRecoLight;

  float sumEt;
  float et;
  float etReco;
  float phi;
  float phiReco;
  float impEtReco;
  float impPhiReco;
  float impSumEt;

  float ptRecoHeavy;
  float phiRecoHeavy;
  float etaRecoHeavy;
  float csvRecoHeavy;
  float ptHeavy;
  float genPtHeavy;
  float etaHeavy;
  float phiHeavy;
  float massRecoHeavy;


  genTree->Branch("X1",       &X1,     "X1/F");
  genTree->Branch("X2",       &X2,     "X2/F");
  genTree->Branch("PtTTH",    &PtTTH,  "PtTTH/F");
  genTree->Branch("PzTTH",    &PzTTH,  "PzTTH/F");
  genTree->Branch("PhiTTH",   &PhiTTH, "PhiTTH/F");
  genTree->Branch("BetaTTH",  &BetaTTH,     "BetaTTH/F");
  genTree->Branch("GammaTTH", &GammaTTH,     "GammaTTH/F");
  genTree->Branch("GammaTT",  &GammaTT,      "GammaTT/F");
  genTree->Branch("DeltaTTH", &DeltaTTH,     "DeltaTTH/F");
  genTree->Branch("PtHStar",  &PtHStar,       "PtHStar/F");
  genTree->Branch("DphiHStar",&DphiHStar,    "DphiHStar/F");
  genTree->Branch("MassTT",   &MassTT,       "MassTT/F");
  genTree->Branch("BetaW",    &BetaW,        "BetaW/F");
  genTree->Branch("GammaW",   &GammaW,       "GammaW/F");
  genTree->Branch("DeltaW",   &DeltaW,       "DeltaW/F");
  genTree->Branch("BetaWLep", &BetaWLep,     "BetaWLep/F");
  genTree->Branch("GammaWLep",&GammaWLep,    "GammaWLep/F");
  genTree->Branch("DeltaWLep",&DeltaWLep,    "DeltaWLep/F");
  genTree->Branch("PxT",   &PxT,       "PxT/F");
  genTree->Branch("PyT",   &PyT,       "PyT/F");
  genTree->Branch("PzT",   &PzT,       "PzT/F");
  genTree->Branch("PxH",   &PxH,       "PxH/F");
  genTree->Branch("PyH",   &PyH,       "PyH/F");
  genTree->Branch("PzH",   &PzH,       "PzH/F");

  genJetLightTree->Branch("ptReco",   &ptRecoLight, "ptReco/F");
  genJetLightTree->Branch("etaReco",  &etaRecoLight, "etaReco/F");
  genJetLightTree->Branch("phiReco",  &phiRecoLight, "phiReco/F");
  genJetLightTree->Branch("csvReco",   &csvRecoLight, "csvReco/F");
  genJetLightTree->Branch("pt",       &ptLight,     "pt/F");
  genJetLightTree->Branch("ptGen",    &genPtLight,     "ptGen/F");
  genJetLightTree->Branch("eta",      &etaLight,    "eta/F");
  genJetLightTree->Branch("phi",      &phiLight,    "phi/F");
  genJetLightTree->Branch("massReco", &massRecoLight,   "massReco/F");

  genEventTree->Branch("sumEt",    &sumEt,   "sumEt/F");
  genEventTree->Branch("et",       &et,    "et/F");
  genEventTree->Branch("etReco",   &etReco,"etReco/F");
  genEventTree->Branch("phi",      &phi,   "phi/F");
  genEventTree->Branch("phiReco",  &phiReco,"phiReco/F");
  genEventTree->Branch("impEtReco",  &impEtReco, "impEtReco/F");
  genEventTree->Branch("impPhiReco", &impPhiReco,"impPhiReco/F");
  genEventTree->Branch("impSumEt",    &impSumEt,   "impSumEt/F");


  genJetHeavyTree->Branch("ptReco",   &ptRecoHeavy, "ptReco/F");
  genJetHeavyTree->Branch("etaReco",  &etaRecoHeavy, "etaReco/F");
  genJetHeavyTree->Branch("phiReco",  &phiRecoHeavy, "phiReco/F");
  genJetHeavyTree->Branch("csvReco",   &csvRecoHeavy, "csvReco/F");
  genJetHeavyTree->Branch("pt",       &ptHeavy,     "pt/F");
  genJetHeavyTree->Branch("ptGen",    &genPtHeavy,     "ptGen/F");
  genJetHeavyTree->Branch("eta",      &etaHeavy,    "eta/F");
  genJetHeavyTree->Branch("phi",      &phiHeavy,    "phi/F");
  genJetHeavyTree->Branch("massReco",     &massRecoHeavy,   "massReco/F");


  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi              (in.getParameter<double>("lumi") );
  bool verbose             (in.getParameter<bool>("verbose") );
  bool computeCMSVariables (in.getParameter<bool>("computeCMSVariables") );
  bool printP4             (in.getParameter<bool>("printP4") );
  int toPrint              (in.getParameter<int>("toPrint") );

  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }

  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    float hJetsgenpt [999];
    float aJetsgenpt [999];
    float hJetsgeneta[999];
    float aJetsgeneta[999];
    float hJetsgenphi[999];
    float aJetsgenphi[999];
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;

    currentTree->SetBranchAddress("jet1",   &jet1);
    currentTree->SetBranchAddress("jet2",   &jet2);
    currentTree->SetBranchAddress("jet3",   &jet3);
    currentTree->SetBranchAddress("jet4",   &jet4);
    currentTree->SetBranchAddress("jet5",   &jet5);
    currentTree->SetBranchAddress("jet6",   &jet6);
    currentTree->SetBranchAddress("jet7",   &jet7);
    currentTree->SetBranchAddress("jet8",   &jet8);
    currentTree->SetBranchAddress("jet9",   &jet9);
    currentTree->SetBranchAddress("jet10",  &jet10);
    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("hJet_genPt",  hJetsgenpt);
    currentTree->SetBranchAddress("aJet_genPt",  aJetsgenpt);
    currentTree->SetBranchAddress("hJet_genEta",  hJetsgeneta);
    currentTree->SetBranchAddress("aJet_genEta",  aJetsgeneta);
    currentTree->SetBranchAddress("hJet_genPhi",  hJetsgenphi);
    currentTree->SetBranchAddress("aJet_genPhi",  aJetsgenphi);


    int printed = 0;
    Long64_t nentries = currentTree->GetEntries();
    
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);

      //continue;
      
      TLorentzVector topBLV   ;
      TLorentzVector topW1LV  ; 
      TLorentzVector topW2LV  ;
      TLorentzVector atopBLV  ; 
      TLorentzVector atopW1LV ; 
      TLorentzVector atopW2LV ;
      TLorentzVector genBLV   ;
      TLorentzVector genBbarLV;
      if(genTop.bmass>0){
	topBLV.SetPtEtaPhiM(  genTop.bpt,genTop.beta,genTop.bphi,genTop.bmass );
	topW1LV.SetPtEtaPhiM( genTop.wdau1pt,genTop.wdau1eta, genTop.wdau1phi,genTop.wdau1mass);
	topW2LV.SetPtEtaPhiM( genTop.wdau2pt,genTop.wdau2eta, genTop.wdau2phi,genTop.wdau2mass);
      }
      if(genTbar.bmass>0){
	atopBLV.SetPtEtaPhiM(  genTbar.bpt,genTbar.beta,genTbar.bphi,genTbar.bmass );
	atopW1LV.SetPtEtaPhiM( genTbar.wdau1pt,genTbar.wdau1eta, genTbar.wdau1phi,genTbar.wdau1mass);
	atopW2LV.SetPtEtaPhiM( genTbar.wdau2pt,genTbar.wdau2eta,genTbar.wdau2phi,genTbar.wdau2mass);
     }
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
      

      //cout << topBLV.Pt() << ", " << topW1LV.Pt() << ", " << topW2LV.Pt() << ", " << atopBLV.Pt() << ", " << atopW1LV.Pt() << ", " << atopW2LV.Pt() << ", " << genBLV.Pt() << ", " << genBbarLV.Pt() << endl;

      vector<TLorentzVector> myJets;
      std::map<unsigned int, JetByPt> map;
      vector<TLorentzVector> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
      map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
      map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
      map[8] = jet9;  map[9] = jet10; 
    
      if(jet1.pt>0){
	TLorentzVector jet1LV;
	jet1LV.SetPtEtaPhiM(jet1.pt, jet1.eta, jet1.phi, jet1.mass);
	myJets.push_back( jet1LV );
      }
      if(jet2.pt>0){
	TLorentzVector jet2LV;
	jet2LV.SetPtEtaPhiM(jet2.pt, jet2.eta, jet2.phi, jet2.mass);
	myJets.push_back( jet2LV );
      }
      if(jet3.pt>0){
	TLorentzVector jet3LV;
	jet3LV.SetPtEtaPhiM(jet3.pt, jet3.eta, jet3.phi, jet3.mass);
	myJets.push_back( jet3LV );
      }   
      if(jet4.pt>0){
	TLorentzVector jet4LV;
	jet4LV.SetPtEtaPhiM(jet4.pt, jet4.eta, jet4.phi, jet4.mass);
	myJets.push_back( jet4LV );
      }
      if(jet5.pt>0){
	TLorentzVector jet5LV;
	jet5LV.SetPtEtaPhiM(jet5.pt, jet5.eta, jet5.phi, jet5.mass);
	myJets.push_back( jet5LV );
      }
      if(jet6.pt>0){
	TLorentzVector jet6LV;
	jet6LV.SetPtEtaPhiM(jet6.pt, jet6.eta, jet6.phi, jet6.mass);
	myJets.push_back( jet6LV );
      }
      if(jet7.pt>0){
	TLorentzVector jet7LV;
	jet7LV.SetPtEtaPhiM(jet7.pt, jet7.eta, jet7.phi, jet7.mass);
	myJets.push_back( jet7LV );
      }
      if(jet8.pt>0){
	TLorentzVector jet8LV;
	jet8LV.SetPtEtaPhiM(jet8.pt, jet8.eta, jet8.phi, jet8.mass);
	myJets.push_back( jet8LV );
      }
      if(jet9.pt>0){
	TLorentzVector jet9LV;
	jet9LV.SetPtEtaPhiM(jet9.pt, jet9.eta, jet9.phi, jet9.mass);
	myJets.push_back( jet9LV );
      }
      if(jet10.pt>0){
	TLorentzVector jet10LV;
	jet10LV.SetPtEtaPhiM(jet10.pt, jet10.eta, jet10.phi, jet10.mass);
	myJets.push_back( jet10LV );
      }
      
      unsigned int iter = 0;    
      for(unsigned int k = 0; k<myJets.size() ; k++){
	
	if(myJets[k].Pt()>20 && TMath::Abs(myJets[k].Eta()) < 2.5 ){
	  
	  myJetsFilt.push_back( myJets[k] );
	  
	  switch(k){
	  case 0:
	    mapFilt[iter] = jet1;
	    break;
	  case 1:
	    mapFilt[iter] = jet2;
	    break;
	  case 2:
	    mapFilt[iter] = jet3;
	    break;
	  case 3:
	    mapFilt[iter] = jet4;
	    break;
	  case 4:
	    mapFilt[iter] = jet5;
	    break;
	  case 5:
	    mapFilt[iter] = jet6;
	    break;
	  case 6:
	    mapFilt[iter] = jet7;
	    break;
	  case 7:
	    mapFilt[iter] = jet8;
	    break;
	  case 8:
	    mapFilt[iter] = jet9;
	    break;
	  case 9:
	    mapFilt[iter] = jet10;
	    break;
	  default:
	    break;
	  }
	  
	  iter++;
	}
      }
      
      int numJets20     =0;     
      int numJets30     =0; 
      int numJets60     =0; 
      int numJets20BtagL=0; 
      int numJets20BtagM=0; 
      int numJets20BtagT=0;
      int numJets30BtagL=0; 
      int numJets30BtagM=0; 
      int numJets30BtagT=0;
      int numJets60BtagL=0; 
      int numJets60BtagM=0; 
      int numJets60BtagT=0;

      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
	
	if( pt_k>20 && csv_k>0.244 ) numJets20BtagL++;
	if( pt_k>20 && csv_k>0.679 ) numJets20BtagM++;
	if( pt_k>20 && csv_k>0.898 ) numJets20BtagT++;
	if( pt_k>30 && csv_k>0.244 ) numJets30BtagL++;
	if( pt_k>30 && csv_k>0.679 ) numJets30BtagM++;
	if( pt_k>30 && csv_k>0.898 ) numJets30BtagT++;
	if( pt_k>60 && csv_k>0.244 ) numJets60BtagL++;
	if( pt_k>60 && csv_k>0.679 ) numJets60BtagM++;
	if( pt_k>60 && csv_k>0.898 ) numJets60BtagT++;
	if( pt_k>20     ) numJets20++;
	if( pt_k>30     ) numJets30++;
	if( pt_k>60     ) numJets60++;
	
	//bool isB    = genBLV.Pt()>0 ?    deltaR( myJetsFilt[k], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt()<0.30        : false;
	//bool isBbar = genBbarLV.Pt()>0 ? deltaR( myJetsFilt[k], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.30  : false;      
      }

      if( myJetsFilt.size()<4 /*|| numJets30<6 || numJets30BtagM<4*/) continue;


      if(  !((abs(genTop.wdau1id)<6 &&  abs(genTbar.wdau1id)>6) ||  
	     (abs(genTop.wdau1id)>6 &&  abs(genTbar.wdau1id)<6))  )
	continue; // semileptonic decays 

      ///////////////////////////////////////////////////////////////////////////

      ptRecoLight=-99;
      etaRecoLight=-99;
      phiRecoLight=-99;
      csvRecoLight=-99;
      ptLight=-99;
      etaLight=-99;
      phiLight=-99;
      massRecoLight=-99;
      ptRecoHeavy=-99;
      etaRecoHeavy=-99;
      phiRecoHeavy=-99;
      csvRecoHeavy=-99;
      ptHeavy=-99;
      etaHeavy=-99;
      phiHeavy=-99;
      massRecoHeavy=-99;

      //////////////////////////////////////////////////
      // met
      //////////////////////////////////////////////////

      sumEt   = METtype1p2corr.sumet;
      etReco  = METtype1p2corr.et;
      phiReco = METtype1p2corr.phi;


      TVector3 recoMEt( etReco*TMath::Cos(phiReco), etReco*TMath::Sin(phiReco), 0.0) ;
      //recoMEt.SetPtEtaPhiM( etReco, 0.0, phiReco, 0.0);
      impSumEt   = sumEt;
      //impEtReco  = etReco;
      //impPhiReco = phiReco;

      if(abs(genTop.wdau1id)==12 || abs(genTop.wdau1id)==14 || abs(genTop.wdau1id)==16){
	et  = genTop.wdau1pt;
	phi = genTop.wdau1phi;
      }
      else if(abs(genTop.wdau2id)==12 || abs(genTop.wdau2id)==14 || abs(genTop.wdau2id)==16){
	et  = genTop.wdau2pt;
	phi = genTop.wdau2phi;
      }
      else if(abs(genTbar.wdau1id)==12 || abs(genTbar.wdau1id)==14 || abs(genTbar.wdau1id)==16){
	et  = genTbar.wdau1pt;
	phi = genTbar.wdau1phi;
      }
      else if(abs(genTbar.wdau2id)==12 || abs(genTbar.wdau2id)==14 || abs(genTbar.wdau2id)==16){
	et  = genTbar.wdau2pt;
	phi = genTbar.wdau2phi;
      }
      else{
	et  = -99;
	phi = -99;
      }
 
      //genEventTree->Fill();

      ///////////////////////////////////////////////////
      // Heavy
      ///////////////////////////////////////////////////
      
      unsigned int indexB    = 999;
      unsigned int indexBbar = 999;
      unsigned int indexH1   = 999;
      unsigned int indexH2   = 999;
      unsigned int indexW1   = 999;
      unsigned int indexW2   = 999;

      int whichTopHad = abs(genTop.wdau1id)<6 ? 0 : 1;
      TLorentzVector wCand1 = whichTopHad==0 ? topW1LV : atopW1LV;
      TLorentzVector wCand2 = whichTopHad==0 ? topW2LV : atopW2LV;
      
      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	if(deltaR(myJetsFilt[k],topBLV)<0.3)
	  indexB = k;
	if(deltaR(myJetsFilt[k],atopBLV)<0.3)
	  indexBbar = k;
	if(genBLV.E()>0 && deltaR(myJetsFilt[k],genBLV)<0.3)
	  indexH1 = k;
	if(genBbarLV.E()>0 && deltaR(myJetsFilt[k],genBbarLV)<0.3)
	  indexH2 = k;
	if(deltaR(myJetsFilt[k],wCand1)<0.3)
	  indexW1 = k;
	if(deltaR(myJetsFilt[k],wCand2)<0.3)
	  indexW2 = k;
      }


      genPtHeavy = -99; 
      if(indexB!=999 && indexB!=indexBbar && indexB!=indexH1 && indexB!=indexH2 && indexB!=indexW1 && indexB!=indexW2){
	// CHECK!!!!!!
	//ptRecoHeavy    = mapFilt[indexB].index>=0 ? hJetsgenpt[mapFilt[indexB].index] :  aJetsgenpt[-mapFilt[indexB].index-1] ;	 
	ptRecoHeavy    = myJetsFilt[indexB].E();
	etaRecoHeavy   = myJetsFilt[indexB].Eta();
	phiRecoHeavy   = myJetsFilt[indexB].Phi();
	massRecoHeavy  = myJetsFilt[indexB].M();
	csvRecoHeavy   = mapFilt[indexB].csv>0 ? mapFilt[indexB].csv : 0.0 ;
	
	float dPx = mapFilt[indexB].index>=0 ? 
	  hJetsgenpt[mapFilt[indexB].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexB].index] ) - myJetsFilt[indexB].Px():
	  aJetsgenpt[-mapFilt[indexB].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexB].index-1] ) - myJetsFilt[indexB].Px();
	float dPy = mapFilt[indexB].index>=0 ? 
	  hJetsgenpt[mapFilt[indexB].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexB].index] ) - myJetsFilt[indexB].Py():
	  aJetsgenpt[-mapFilt[indexB].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexB].index-1] ) - myJetsFilt[indexB].Py();

	if(mapFilt[indexB].index>=0)
	  genPtHeavy =  hJetsgenpt[mapFilt[indexB].index];
	else if(mapFilt[indexB].index!=-99)
	  genPtHeavy  = aJetsgenpt[-mapFilt[indexB].index-1];
	else{}

	TVector3 shift( dPx, dPy , 0.0);
	recoMEt  -= shift;
	impSumEt += (topBLV.E()-myJetsFilt[indexB].E());
	
      }else{ 
	ptRecoHeavy=-99;
	etaRecoHeavy=-99;
	phiRecoHeavy=-99;
	csvRecoHeavy=-99;
	massRecoHeavy=-99;
      }
      ptHeavy       = topBLV.E();
      etaHeavy      = topBLV.Eta();
      phiHeavy      = topBLV.Phi();
      if(topBLV.E()>0) genJetHeavyTree->Fill();
      

      genPtHeavy = -99; 
      if(indexBbar!=999 && indexBbar!=indexB && indexBbar!=indexH1 && indexBbar!=indexH2 && indexBbar!=indexW1 && indexBbar!=indexW2){
	//ptRecoHeavy    = mapFilt[indexBbar].index>=0 ? hJetsgenpt[mapFilt[indexBbar].index] :  aJetsgenpt[-mapFilt[indexBbar].index-1] ;	 
	ptRecoHeavy    = myJetsFilt[indexBbar].E();
	massRecoHeavy  = myJetsFilt[indexBbar].M();
	csvRecoHeavy   = mapFilt[indexBbar].csv>0 ? mapFilt[indexBbar].csv : 0.0 ;
	etaRecoHeavy   = myJetsFilt[indexBbar].Eta();
	phiRecoHeavy   = myJetsFilt[indexBbar].Phi();

	float dPx = mapFilt[indexBbar].index>=0 ? 
	  hJetsgenpt[mapFilt[indexBbar].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexBbar].index] ) - myJetsFilt[indexBbar].Px():
	  aJetsgenpt[-mapFilt[indexBbar].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexBbar].index-1] ) - myJetsFilt[indexBbar].Px();
	float dPy = mapFilt[indexBbar].index>=0 ? 
	  hJetsgenpt[mapFilt[indexBbar].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexBbar].index] ) - myJetsFilt[indexBbar].Py():
	  aJetsgenpt[-mapFilt[indexBbar].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexBbar].index-1] ) - myJetsFilt[indexBbar].Py();

	if(mapFilt[indexBbar].index>=0)
	  genPtHeavy =  hJetsgenpt[mapFilt[indexBbar].index];
	else if(mapFilt[indexBbar].index!=-99)
	  genPtHeavy  = aJetsgenpt[-mapFilt[indexBbar].index-1];
	else{}

	TVector3 shift( dPx , dPy , 0.0);

	if(indexB!=indexBbar) {
	  recoMEt  -= shift;
	  impSumEt += (atopBLV.E()-myJetsFilt[indexBbar].E());
	}

      }else{
	ptRecoHeavy=-99;
	csvRecoHeavy=-99;
	massRecoHeavy=-99;
	etaRecoHeavy=-99;
	phiRecoHeavy=-99;
      }
      ptHeavy       = atopBLV.E();
      etaHeavy      = atopBLV.Eta();
      phiHeavy      = atopBLV.Phi();
      if(atopBLV.E()>0) genJetHeavyTree->Fill();


      genPtHeavy = -99; 
      if(indexH1!=999 && indexH1!=indexB && indexH1!=indexBbar && indexH1!=indexH2 && indexH1!=indexW1 && indexH1!=indexW2){
	//ptRecoHeavy    = mapFilt[indexH1].index>=0 ? hJetsgenpt[mapFilt[indexH1].index] :  aJetsgenpt[-mapFilt[indexH1].index-1] ;	 
	ptRecoHeavy    = myJetsFilt[indexH1].E();
	massRecoHeavy  = myJetsFilt[indexH1].M();
	csvRecoHeavy   = mapFilt[indexH1].csv>0 ? mapFilt[indexH1].csv : 0.0 ;
	etaRecoHeavy   = myJetsFilt[indexH1].Eta();
	phiRecoHeavy   = myJetsFilt[indexH1].Phi();
	
	float dPx = mapFilt[indexH1].index>=0 ? 
	  hJetsgenpt[mapFilt[indexH1].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexH1].index] ) - myJetsFilt[indexH1].Px():
	  aJetsgenpt[-mapFilt[indexH1].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexH1].index-1] ) - myJetsFilt[indexH1].Px();
	float dPy = mapFilt[indexH1].index>=0 ? 
	  hJetsgenpt[mapFilt[indexH1].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexH1].index] ) - myJetsFilt[indexH1].Py():
	  aJetsgenpt[-mapFilt[indexH1].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexH1].index-1] ) - myJetsFilt[indexH1].Py();

	if(mapFilt[indexH1].index>=0)
	  genPtHeavy =  hJetsgenpt[mapFilt[indexH1].index];
	else if(mapFilt[indexH1].index!=-99)
	  genPtHeavy  = aJetsgenpt[-mapFilt[indexH1].index-1];
	else{}


	TVector3 shift( dPx , dPy , 0.0);

	if(indexH1!=indexBbar && indexH1!=indexB) {
	  recoMEt  -= shift;
	  impSumEt += (genBLV.E()-myJetsFilt[indexH1].E());
	}
	
      }
      else{
	  ptRecoHeavy=-99;
	  csvRecoHeavy=-99;
	  massRecoHeavy=-99;
	  etaRecoHeavy=-99;
	  phiRecoHeavy=-99;
	}
      ptHeavy       = genBLV.E()>0 ? genBLV.E() : -999;
      etaHeavy      = genBLV.E()>0 ? genBLV.Eta(): -999;
      phiHeavy      = genBLV.E()>0 ? genBLV.Phi(): -999;
      if(genBLV.E()>0) genJetHeavyTree->Fill();
	
      genPtHeavy = -99; 
      if(indexH2!=999 && indexH2!=indexB && indexH2!=indexBbar && indexH2!=indexH1 && indexH2!=indexW1 && indexH2!=indexW2){
	//ptRecoHeavy    = mapFilt[indexH2].index>=0 ? hJetsgenpt[mapFilt[indexH2].index] :  aJetsgenpt[-mapFilt[indexH2].index-1] ;	 
	ptRecoHeavy    = myJetsFilt[indexH2].E();
	massRecoHeavy  = myJetsFilt[indexH2].M();
	csvRecoHeavy   = mapFilt[indexH2].csv>0 ? mapFilt[indexH2].csv : 0.0 ;
	etaRecoHeavy   = myJetsFilt[indexH2].Eta();
	phiRecoHeavy   = myJetsFilt[indexH2].Phi();

	float dPx = mapFilt[indexH2].index>=0 ? 
	  hJetsgenpt[mapFilt[indexH2].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexH2].index] ) - myJetsFilt[indexH2].Px():
	  aJetsgenpt[-mapFilt[indexH2].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexH2].index-1] ) - myJetsFilt[indexH2].Px();
	float dPy = mapFilt[indexH2].index>=0 ? 
	  hJetsgenpt[mapFilt[indexH2].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexH2].index] ) - myJetsFilt[indexH2].Py():
	  aJetsgenpt[-mapFilt[indexH2].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexH2].index-1] ) - myJetsFilt[indexH2].Py();
	
	if(mapFilt[indexH2].index>=0)
	  genPtHeavy =  hJetsgenpt[mapFilt[indexH2].index];
	else if(mapFilt[indexH2].index!=-99)
	  genPtHeavy  = aJetsgenpt[-mapFilt[indexH2].index-1];
	else{}

	TVector3 shift( dPx, dPy, 0.0);

	if(indexH2!=indexBbar && indexH2!=indexB && indexH2!=indexH1) {
	  recoMEt  -= shift;
	  impSumEt += (genBbarLV.E()-myJetsFilt[indexH2].E());
	}
	
      }
      else{
	ptRecoHeavy=-99;
	csvRecoHeavy=-99;
	massRecoHeavy=-99;
	etaRecoHeavy=-99;
	phiRecoHeavy=-99;
      }
      ptHeavy       = genBbarLV.E()>0 ? genBbarLV.E() : -999;
      etaHeavy      = genBbarLV.E()>0 ? genBbarLV.Eta(): -999;
      phiHeavy      = genBbarLV.E()>0 ? genBbarLV.Phi(): -999;
      if(genBbarLV.E()>0) genJetHeavyTree->Fill();
      
      ///////////////////////////////////////////////////
      // Light
      ///////////////////////////////////////////////////
      
      genPtLight = -99;
      if(indexW1!=999 && indexW1!=indexB && indexW1!=indexBbar && indexW1!=indexH1 && indexW1!=indexH2 && indexW1!=indexW2){  // CHECK!!!!!!
	//ptRecoLight    = mapFilt[indexW1].index>=0 ? hJetsgenpt[mapFilt[indexW1].index] :  aJetsgenpt[-mapFilt[indexW1].index-1] ;
	ptRecoLight    = myJetsFilt[indexW1].E();
	etaRecoLight   = myJetsFilt[indexW1].Eta();
	phiRecoLight   = myJetsFilt[indexW1].Phi();
	massRecoLight  = myJetsFilt[indexW1].M();
	csvRecoLight   = mapFilt[indexW1].csv>0 ? mapFilt[indexW1].csv : 0.0 ;
	
	
	float dPx = mapFilt[indexW1].index>=0 ? 
	  hJetsgenpt[mapFilt[indexW1].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexW1].index] ) - myJetsFilt[indexW1].Px():
	  aJetsgenpt[-mapFilt[indexW1].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexW1].index-1] ) - myJetsFilt[indexW1].Px();
	float dPy = mapFilt[indexW1].index>=0 ? 
	  hJetsgenpt[mapFilt[indexW1].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexW1].index] ) - myJetsFilt[indexW1].Py():
	  aJetsgenpt[-mapFilt[indexW1].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexW1].index-1] ) - myJetsFilt[indexW1].Py();
	
	if(mapFilt[indexW1].index>=0)
	  genPtLight =  hJetsgenpt[mapFilt[indexW1].index];
	else if(mapFilt[indexW1].index!=-99)
	  genPtLight  = aJetsgenpt[-mapFilt[indexW1].index-1];
	else{}

	  TVector3 shift(dPx ,dPy , 0.0);
	  
	  if(indexW1!=indexBbar && indexW1!=indexB && indexW1!=indexH1 &&  indexW1!=indexH2) {
	    recoMEt  -= shift;
	    impSumEt += (wCand1.E()-myJetsFilt[indexW1].E());
	  }


	}
	else{
	  ptRecoLight=-99;
	  csvRecoLight=-99;
	  massRecoLight=-99;
	  etaRecoLight=-99;
	  phiRecoLight=-99;
	}
	ptLight       = wCand1.E();
	etaLight      = wCand1.Eta();
	phiLight      = wCand1.Phi();
	if(wCand1.E()>0) genJetLightTree->Fill();
	
	genPtLight = -99;
	if(indexW2!=999 && indexW2!=indexB && indexW2!=indexBbar && indexW2!=indexH1 && indexW2!=indexH2 && indexW2!=indexW1){
	  //ptRecoLight    = mapFilt[indexW2].index>=0 ? hJetsgenpt[mapFilt[indexW2].index] :  aJetsgenpt[-mapFilt[indexW2].index-1] ;
	  ptRecoLight    = myJetsFilt[indexW2].E();
	  massRecoLight  = myJetsFilt[indexW2].M();
	  etaRecoLight   = myJetsFilt[indexW2].Eta();
	  phiRecoLight   = myJetsFilt[indexW2].Phi();
	  csvRecoLight   =  mapFilt[indexW2].csv>0 ? mapFilt[indexW2].csv : 0.0 ;

	  float dPx = mapFilt[indexW2].index>=0 ? 
	    hJetsgenpt[mapFilt[indexW2].index]*TMath::Cos(  hJetsgenphi[mapFilt[indexW2].index] ) - myJetsFilt[indexW2].Px():
	    aJetsgenpt[-mapFilt[indexW2].index-1]*TMath::Cos(  aJetsgenphi[-mapFilt[indexW2].index-1] ) - myJetsFilt[indexW2].Px();
	  float dPy = mapFilt[indexW2].index>=0 ? 
	    hJetsgenpt[mapFilt[indexW2].index]*TMath::Sin(  hJetsgenphi[mapFilt[indexW2].index] ) - myJetsFilt[indexW2].Py():
	    aJetsgenpt[-mapFilt[indexW2].index-1]*TMath::Sin(  aJetsgenphi[-mapFilt[indexW2].index-1] ) - myJetsFilt[indexW2].Py();

	  if(mapFilt[indexW2].index>=0)
	    genPtLight =  hJetsgenpt[mapFilt[indexW2].index];
	  else if(mapFilt[indexW2].index!=-99)
	    genPtLight  = aJetsgenpt[-mapFilt[indexW2].index-1];
	  else{}

	  TVector3 shift( dPx, dPy, 0.0);

	  if(indexW2!=indexBbar && indexW2!=indexB && indexW2!=indexH2 &&  indexW2!=indexH1 &&  indexW2!=indexW1) {
	    recoMEt  -= shift;
	    impSumEt += (wCand2.E()-myJetsFilt[indexW2].E());
	  }

	}else{
	  ptRecoLight=-99;
	  csvRecoLight=-99;
	  massRecoLight=-99;
	  etaRecoLight=-99;
	  phiRecoLight=-99;
	}	
	ptLight       = wCand2.E();
	etaLight      = wCand2.Eta();
	phiLight      = wCand2.Phi();
	if(wCand2.E()>0) genJetLightTree->Fill();
	
      
	impEtReco  = recoMEt.Pt();
	impPhiReco = recoMEt.Phi();
	genEventTree->Fill();




      X1 = -99;
      X2 = -99;
      PtTTH   = -99;
      PhiTTH   = -99;
      BetaTTH   = -99;
      GammaTTH     = -99;
      GammaTT      = -99;
      DeltaTTH     = -99;
      PtHStar  = -99;
      MassTT   = -99;
      DphiHStar  = -99;
      BetaW   = -99;
      GammaW   = -99;
      DeltaW   = -99;
      BetaWLep = -99;
      GammaWLep = -99;
      DeltaWLep = -99;

      if(computeCMSVariables && genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0){

	TLorentzVector TOPLEP;
	TLorentzVector TOPHAD;
	TLorentzVector TOPHADW1;
	TLorentzVector TOPHADW2;
	TLorentzVector TOPHADB;
	TLorentzVector TOPLEPW1;
	TLorentzVector TOPLEPW2;
	TLorentzVector TOPLEPB;
	TLorentzVector HIGGS;

	int properEvent = 1;

	if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	  TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	  TOPLEPW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
	}
	else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	  TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	  TOPLEPW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	}
	
	else{properEvent=0;}
	HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());

	if(properEvent){


	  if(printP4 && printed<toPrint){
	    cout << "jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta() << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << endl;
	    cout << "jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta() << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << endl;
	    cout << "jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << "," <<  TOPLEPB.Eta() << "," << TOPLEPB.Phi() << "," << TOPLEPB.M() << endl;
	    cout << "jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta() << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << endl;
	    cout << "jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta() << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << endl;
	    cout << "jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << "," <<  TOPHADB.Eta() << "," << TOPHADB.Phi() << "," << TOPHADB.M() << endl;
	    cout << "jet7.SetPtEtaPhiM(" << genBLV.Pt() << "," <<  genBLV.Eta() << "," << genBLV.Phi() << "," << genBLV.M() << endl;
	    cout << "jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << "," <<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << endl;
	    printed++;
	  }

	  TLorentzVector TOT = TOPLEP+TOPHAD+HIGGS;
	  TVector3 boostToCMS       = TOT.BoostVector();
	  TVector3 boostToTopHadCMS = TOPHAD.BoostVector();
	  TVector3 boostToTopLepCMS = TOPLEP.BoostVector();

	  PtTTH  = TOT.Pt();
	  PhiTTH = TMath::Cos(TOT.Phi());
	  PzTTH  = TOT.Pz();
	  X1 = (  TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);
	  X2 = ( -TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);
	  
	  //trasform
	  float tmp1 = X1;
	  float tmp2 = X2;
	  X1 = TMath::Sqrt(tmp1*tmp2);
	  X2 = tmp1-tmp2;
	  
	  TOPLEP.Boost(-boostToCMS);
	  TOPHAD.Boost(-boostToCMS);
	  HIGGS .Boost(-boostToCMS);
	  
	  PxT = TOPLEP.Px();
	  PyT = TOPLEP.Py();
	  PzT = TOPLEP.Pz();
	  PxH = HIGGS.Px();
	  PyH = HIGGS.Py();
	  PzH = HIGGS.Pz();
	  
	  if(TOPLEP.Mag()<=0 || TOPHAD.Mag()<=0) cout << "Error: null momentum" << endl;
	  
	  TVector3 decayPlane    = ((TOPHAD.Vect()).Cross(TOPLEP.Vect())).Unit();
	  TVector3 aux1 = (boostToCMS.Cross( TVector3(0,0,1) )).Unit();
	  TVector3 aux2 = (boostToCMS.Cross( aux1 )).Unit();
	  TVector3 aux3 = boostToCMS.Unit();
	  
	  TVector3 aux4( aux1.Px()*TMath::Cos(decayPlane.Angle(aux1)) , aux1.Py()*TMath::Cos(decayPlane.Angle(aux1)), aux1.Pz()*TMath::Cos(decayPlane.Angle(aux1)));
	  TVector3 aux5( aux2.Px()*TMath::Cos(decayPlane.Angle(aux2)) , aux2.Py()*TMath::Cos(decayPlane.Angle(aux2)), aux2.Pz()*TMath::Cos(decayPlane.Angle(aux2)));
	  
	  //TVector3 decayPlanePhi = decayPlane.Cross( TVector3(0,0,1) );
	  
	  BetaTTH  = TMath::Cos(decayPlane.Angle(aux3));
	  DeltaTTH = TMath::Cos((aux4+aux5).Angle(aux1));
	  GammaTTH = TMath::Cos((HIGGS.Vect()).Angle( boostToCMS ));
	  
	  PtHStar   = HIGGS.P()/TOT.M();
	  DphiHStar = TMath::Cos((HIGGS.Vect()).Angle(TOPLEP.Vect()));
	  
	  TOPLEP.Boost( -(TOPLEP+TOPHAD).BoostVector() ) ;
	  TOPHAD.Boost( -(TOPLEP+TOPHAD).BoostVector() ) ;
	  
	  GammaTT  = TMath::Cos( (TOPLEP.Vect()).Angle( (TOPLEP+TOPHAD).BoostVector() ) );
	  MassTT   = (TOPLEP+TOPHAD).M();
	  
	  //////////////////////////
	  int chLep = 0;
	  if(abs(genTop.wdau1id)==11  || abs(genTop.wdau1id)==13  || abs(genTop.wdau1id)==15 ||
	     abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15  )
	    chLep = 1;
	  else if(abs(genTop.wdau2id) ==11  || abs(genTop.wdau2id)==13  || abs(genTop.wdau2id)==15 ||
		  abs(genTbar.wdau2id)==11  || abs(genTbar.wdau2id)==13 || abs(genTbar.wdau2id)==15  )
	    chLep = 2;
	  else{}
	  
	  if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	    TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	    TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	    TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	    TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	    TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	    TOPLEPW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	    TOPLEPW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	    TOPLEPB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	  }
	  else{
	    TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	    TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	    TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	    TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	    TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	    TOPLEPW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	    TOPLEPW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	    TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	  }
	  
	  int whichW = TOPHADW1.Pt() > TOPHADW2.Pt() ? 1 : 2;

	  //TOPHAD.Boost(-boostToTopHadCMS);
	  TOPHADW1.Boost(-boostToTopHadCMS);
	  TOPHADW2.Boost(-boostToTopHadCMS);
	  TOPHADB.Boost( -boostToTopHadCMS);
	  
	  TVector3 decayTPlane = (TOPHADB.Vect()).Cross( (TOPHADW1).Vect()  ).Unit();
	  
	  TVector3 boostToWHadCMS = (TOPHADW1+TOPHADW2).BoostVector();
	  TOPHADW1.Boost(-boostToWHadCMS);
	  TOPHADW2.Boost(-boostToWHadCMS);
	  
	  //cout << TOPHADB.P() << endl;
	  //cout << "W1: " << TOPHADW1.E() << " W1: " << TOPHADW2.E() << endl;
	  
	  BetaW     = TMath::Cos((TOPHADB.Vect()).Angle( boostToTopHadCMS  ));
	  DeltaW    = TMath::Cos( decayTPlane.Angle( boostToTopHadCMS )  );
	  GammaW    = i%2==0 ? TMath::Cos( (TOPHADW1.Vect()).Angle(boostToWHadCMS) ) : TMath::Cos( (TOPHADW2.Vect()).Angle(boostToWHadCMS) );  //average over W flavor (unobserved)
	  

	  ///////////////////////////

	  TOPLEPW1.Boost(-boostToTopLepCMS);
	  TOPLEPW2.Boost(-boostToTopLepCMS);
	  TOPLEPB.Boost( -boostToTopLepCMS);
	  
	  TVector3 boostToWLepCMS = (TOPLEPW1+TOPLEPW2).BoostVector();
	  TOPLEPW1.Boost(-boostToWLepCMS);
	  TOPLEPW2.Boost(-boostToWLepCMS);
	  
	  //cout << TOPLEPB.P() << endl;
	  //cout << "W1: " << TOPLEPW1.E() << " W1: " << TOPLEPW2.E() << endl;

	  BetaWLep     = TMath::Cos( (TOPLEPB.Vect()).Angle( boostToTopLepCMS  ));
	  DeltaWLep    = 0.0;
	  if( chLep==1)
	    GammaWLep    = TMath::Cos( (TOPLEPW1.Vect()).Angle(boostToWLepCMS) );
	  else if( chLep==2)
	    GammaWLep    = TMath::Cos( (TOPLEPW2.Vect()).Angle(boostToWLepCMS) );
	  else
	    GammaWLep = -99;
	  
	  ///////
	  
	  genTree->Fill();
	}
      }
      
    } // event loop

    mySamples->GetFile( currentName )->Close();
    
  } // samples loop


  return 0;

}
