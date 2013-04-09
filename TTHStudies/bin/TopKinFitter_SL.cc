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

#include "TopQuarkAnalysis/TopKinFitter/interface/TtSemiLepKinFitter.h"
#include "TopQuarkAnalysis/TopKinFitter/plugins/TtSemiLepKinFitProducer.h"

#include "Bianchi/TTHStudies/interface/Samples.h"



#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;

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


struct sorter {
  bool operator() (Float_t i, Float_t j) const { return (i>j);}
};


Float_t deltaR(LV j1, LV j2){
  return TMath::Sqrt(TMath::Power(j1.eta()-j2.eta(),2) + TMath::Power(j1.phi()-j2.phi(),2));
}

Int_t match(LV j1, LV j2){
  return deltaR(j1,j2) < GENJETDR;
}

bool bookKeeper(unsigned int index1, unsigned int index2, 
		unsigned int index3, unsigned int index4, 
		vector<unsigned int>& visited){
  
  unsigned int combination1 = 1000*index1 + 100*index2 + 10*index3 + 1*index4;
  unsigned int combination2 = 1000*index2 + 100*index1 + 10*index3 + 1*index4;
  bool isThere = false;
  
  for(unsigned int k = 0; k < visited.size() ; k++ )
    if( combination1 == visited[k] || combination2 == visited[k] ) isThere = true;
  
  if(!isThere){
    visited.push_back( combination1 );
    visited.push_back( combination2 );
    return true;
  }
  else{
    return false;
  }

}


int main(int argc, const char* argv[])
{

  std::cout << "TopKinFitter_SL" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  std::vector< edm::ParameterSet > udscResolutions(      in.getParameter<edm::VParameterSet>("udscResolutions") );
  std::vector< edm::ParameterSet > bResolutions(         in.getParameter<edm::VParameterSet>("bResolutions") );
  std::vector< edm::ParameterSet > lepResolutions(       in.getParameter<edm::VParameterSet>("lepResolutions") );
  std::vector<double> metResolutionsCoeff(               in.getParameter<std::vector<double> >("metResolutionsCoeff") );
  std::vector<double> helicityCoeff(                     in.getParameter<std::vector<double> >("helicityCoeff") );
  if(helicityCoeff.size()!=9){
    cout << "Nine parameters are required for the helicity angle param... return" << endl;
    return 0;
  }

  double lumi     (in.getParameter<double>("lumi") );
  bool verbose    (in.getParameter<bool>("verbose") );
  double udscPtCut(in.getParameter<double>("udscPtCut") );
  double bPtCut   (in.getParameter<double>("bPtCut") );
  std::string likelihoodFile(in.getParameter<std::string>("likelihoodFile") );
  vector<std::string> likelihoods(in.getParameter<vector<std::string> >("likelihoods") );

  int useKin  = 0;
  int useBtag = 0;
  int useHel  = 0;

  for(unsigned int l = 0; l < likelihoods.size() ; l++){
    if( likelihoods[l].find("KIN") !=string::npos ) useKin =1;
    if( likelihoods[l].find("BTAG")!=string::npos ) useBtag=1;
    if( likelihoods[l].find("HEL") !=string::npos ) useHel =1;
  }

  TFile* fLikelihoods          = TFile::Open(likelihoodFile.c_str(),"READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");

  TF1* fHelTopLept = new TF1("fHelTopLept",Form("%f + %f*x + %f*x*x + %f*x*x*x",              helicityCoeff[0], helicityCoeff[1], helicityCoeff[2], helicityCoeff[3]),-1,1);
  TF1* fHelTopHadr = new TF1("fHelTopHadr",Form("%f + %f*x + %f*x*x + %f*x*x*x + %f*x*x*x*x", helicityCoeff[4], helicityCoeff[5], helicityCoeff[6], helicityCoeff[7], helicityCoeff[8]),-1,1);


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
    
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    genParticleInfo genB, genBbar;
    Float_t vLepton_mass[99];
    Float_t vLepton_pt  [99];
    Float_t vLepton_eta [99];
    Float_t vLepton_phi [99];
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
    currentTree->SetBranchAddress("vLepton_mass",vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt",vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta",vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi",vLepton_phi);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);


    Long64_t nentries = currentTree->GetEntries();
    
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%100==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      LV genBLV   (0.,0.,0.,0.);
      LV genBbarLV(0.,0.,0.,0.);
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPt(  genB.pt );
	genBLV.SetEta( genB.eta );
	genBLV.SetPhi( genB.phi );
	genBLV.SetM(   genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPt(  genBbar.pt );
	genBbarLV.SetEta( genBbar.eta );
	genBbarLV.SetPhi( genBbar.phi );
	genBbarLV.SetM(   genBbar.mass );
      }
      
      vector<LV> myJets;
      std::map<unsigned int, JetByPt> map;
      vector<LV> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
      map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
      map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
      map[8] = jet9;  map[9] = jet10; 
    
      myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
      myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
      myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
      myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
      myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
      myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
      myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
      myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
      myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
      myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));
      
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
      
      int numJets20 =0;     
      int numJets30 =0; 
      int numJets60 =0; 
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


      ////////////////////////////  TOP KIN FITTER    ////////////////////////////

      std::vector< TtSemiLepKinFitter::Constraint > constraints; 
      constraints.push_back(TtSemiLepKinFitter::kWHadMass);
      constraints.push_back(TtSemiLepKinFitter::kWLepMass);
      constraints.push_back(TtSemiLepKinFitter::kTopHadMass);
      constraints.push_back(TtSemiLepKinFitter::kTopLepMass);
      //constraints.push_back(TtSemiLepKinFitter::kEqualTopMasses);
      //constraints.push_back(TtSemiLepKinFitter::kNeutrinoMass); // causes crash

      std::vector< double > jetEnergyResolutionScaleFactors;
      std::vector< double > jetEnergyResolutionEtaBinning;
      jetEnergyResolutionScaleFactors.push_back(1.);
      jetEnergyResolutionEtaBinning.push_back(0.);
      jetEnergyResolutionEtaBinning.push_back(-1.);


      std::vector< edm::ParameterSet > metResolutions;
      edm::ParameterSet paramSetMet;  
      paramSetMet.addParameter("et",string(Form("sqrt(%f + %f*%.5f)",      metResolutionsCoeff[0],  metResolutionsCoeff[1],  METtype1p2corr.sumet)));
      paramSetMet.addParameter("eta",string("999"));
      paramSetMet.addParameter("phi",string(Form("sqrt(%f + %f*%.5f)/%.5f",metResolutionsCoeff[2],  metResolutionsCoeff[3],  METtype1p2corr.sumet,METtype1p2corr.et)));
      metResolutions .push_back(paramSetMet);
	
      int maxNbIter = 200;
      TtSemiLepKinFitter *fitter = 
	new TtSemiLepKinFitter(TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,
			       maxNbIter,5e-5,1e-4,constraints,80.4,173.,
			       &udscResolutions,&bResolutions,&lepResolutions,&metResolutions,
			       &jetEnergyResolutionScaleFactors,&jetEnergyResolutionEtaBinning);


      LV leptonLV(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);

      TLorentzVector p4Lepton(leptonLV.Px(),leptonLV.Py(),leptonLV.Pz(),leptonLV.E());
      TLorentzVector p4Neutrino(nuPx,nuPy,0,nuE);
      TLorentzVector p4HadB(0.,0.,0.,0.);
      TLorentzVector p4HadP(0.,0.,0.,0.);
      TLorentzVector p4HadQ(0.,0.,0.,0.);
      TLorentzVector p4LepB(0.,0.,0.,0.);

      unsigned int indices[myJetsFilt.size()];
      for(unsigned int k = 0; k<myJetsFilt.size() ; k++)
	indices[k]=k; 
      vector<unsigned int> visited;

      if( myJetsFilt.size()<4 ) continue;

      int counterPer = 0;
      float minLogP  = 999;
      vector<unsigned int> bestRank; 
      do {

	float logp = 999;
	bestRank.push_back(0);	bestRank.push_back(0);
	bestRank.push_back(0);	bestRank.push_back(0);

	unsigned int indexWdau1Cand = indices[0];
	unsigned int indexWdau2Cand = indices[1];
	unsigned int indexbCand     = indices[2];
	unsigned int indexbbarCand  = indices[3];

	if( !((mapFilt[indexbCand ]).csv     > 0. && (mapFilt[indexbbarCand ]).csv  > 0. &&
	    (mapFilt[indexWdau1Cand ]).csv <  0.679 && (mapFilt[indexWdau2Cand ]).csv <  0.679) ) 
	  continue;
	if( !( (mapFilt[indexWdau1Cand ]).pt >  udscPtCut && (mapFilt[indexWdau2Cand ]).pt > udscPtCut  &&
	       (mapFilt[indexbCand ]).pt     >  bPtCut    && (mapFilt[indexbbarCand ]).pt  > bPtCut ) )
	  continue;
	

	if( !bookKeeper( indexWdau1Cand, indexWdau2Cand, indexbCand, indexbbarCand, visited) ) 
	  continue;

	float WMass   = (myJetsFilt[indexWdau1Cand]+myJetsFilt[indexWdau2Cand]).M();
	float TopMass = (myJetsFilt[indexWdau1Cand]+myJetsFilt[indexWdau2Cand]+myJetsFilt[indexbCand]).M();

	if( WMass  <  40 || WMass  > 120 ) continue;
	if( TopMass< 120 || TopMass> 230 ) continue;
	
	counterPer++;
	
	p4HadB.SetPx( (myJetsFilt[indexbCand]).Px() );
	p4HadB.SetPy( (myJetsFilt[indexbCand]).Py() );
	p4HadB.SetPz( (myJetsFilt[indexbCand]).Pz() );
	p4HadB.SetE(  (myJetsFilt[indexbCand]).E()  );

	p4HadP.SetPx( (myJetsFilt[indexWdau1Cand]).Px() );
	p4HadP.SetPy( (myJetsFilt[indexWdau1Cand]).Py() );
	p4HadP.SetPz( (myJetsFilt[indexWdau1Cand]).Pz() );
	p4HadP.SetE(  (myJetsFilt[indexWdau1Cand]).E()  );

	p4HadQ.SetPx( (myJetsFilt[indexWdau2Cand]).Px() );
	p4HadQ.SetPy( (myJetsFilt[indexWdau2Cand]).Py() );
	p4HadQ.SetPz( (myJetsFilt[indexWdau2Cand]).Pz() );
	p4HadQ.SetE(  (myJetsFilt[indexWdau2Cand]).E()  );

	p4LepB.SetPx( (myJetsFilt[indexbbarCand]).Px() );
	p4LepB.SetPy( (myJetsFilt[indexbbarCand]).Py() );
	p4LepB.SetPz( (myJetsFilt[indexbbarCand]).Pz() );
	p4LepB.SetE(  (myJetsFilt[indexbbarCand]).E()  );

	int status = fitter->fit(p4HadP, p4HadQ, p4HadB, p4LepB, p4Lepton, p4Neutrino, -1, CovarianceMatrix::kMuon);
	if(status!=0) continue;

	float chi2 = fitter->fitS();
	float p    = fitter->fitProb();
	if(useKin && p>0.)
	  logp =  -TMath::Log(p);
	
	if(useBtag && fLikelihoodCsvB && fLikelihoodCsvNonB){
	  float btag1   = fLikelihoodCsvB->Eval( (mapFilt[indexbCand    ]).csv);
	  float btag2   = fLikelihoodCsvB->Eval( (mapFilt[indexbbarCand ]).csv);
	  float buntag1 = fLikelihoodCsvNonB->Eval( (mapFilt[indexWdau1Cand ]).csv);
	  float buntag2 = fLikelihoodCsvNonB->Eval( (mapFilt[indexWdau2Cand ]).csv);
	  if( useKin && btag1>0. && btag2>0. && buntag1>0. && buntag2>0.){
	    logp += (-TMath::Log(btag1));
	    logp += (-TMath::Log(btag2));
	    logp += (-TMath::Log(buntag1));
	    logp += (-TMath::Log(buntag2));
	  }
	  else if( !useKin && btag1>0. && btag2>0. && buntag1>0. && buntag2>0.){
	    logp =  (-TMath::Log(btag1));
	    logp += (-TMath::Log(btag2));
	    logp += (-TMath::Log(buntag1));
	    logp += (-TMath::Log(buntag2));
	  }
	  else
	    logp = 999;	  
	}

	

	pat::Particle resHadP     = fitter->fittedHadP();
	pat::Particle resHadQ     = fitter->fittedHadQ();
	pat::Particle resHadB     = fitter->fittedHadB();
	pat::Particle resLepB     = fitter->fittedLepB();
	pat::Particle resLepton   = fitter->fittedLepton();
	pat::Particle resNeutrino = fitter->fittedNeutrino();
	
	TLorentzVector lvHadP(resHadP.px(),resHadP.py(),resHadP.pz(),resHadP.energy());
	TLorentzVector lvHadQ(resHadQ.px(),resHadQ.py(),resHadQ.pz(),resHadQ.energy());
	TLorentzVector lvHadB(resHadB.px(),resHadB.py(),resHadB.pz(),resHadB.energy());
	TLorentzVector lvLepB(resLepB.px(),resLepB.py(),resLepB.pz(),resLepB.energy());
	TLorentzVector lvLepton(resLepton.px(),resLepton.py(),resLepton.pz(),resLepton.energy());
	TLorentzVector lvNeutrino(resNeutrino.px(),resNeutrino.py(),resNeutrino.pz(),resNeutrino.energy());
	TLorentzVector lvHadW   = lvHadP   + lvHadQ ;
	TLorentzVector lvHadTop = lvHadW   + lvHadB ;
	TLorentzVector lvLepW   = lvLepton + lvNeutrino ;
	TLorentzVector lvLepTop = lvLepW   + lvLepB ;

	if(useHel && fHelTopLept && fHelTopHadr){
	  
	  TLorentzVector bHad;
	  TLorentzVector w1Had;
	  TLorentzVector w2Had;
	  w1Had.SetPxPyPzE( lvHadP.Px(), lvHadP.Py(), lvHadP.Pz(), lvHadP.E());
	  w2Had.SetPxPyPzE( lvHadQ.Px(), lvHadQ.Py(), lvHadQ.Pz(), lvHadQ.E());
	  bHad.SetPxPyPzE(  lvHadB.Px(), lvHadB.Py(), lvHadB.Pz(), lvHadB.E());
	  TVector3 boostWHad(lvHadW.Px()/lvHadW.E(), lvHadW.Py()/lvHadW.E(), lvHadW.Pz()/lvHadW.E());
	  w1Had.Boost(-boostWHad);
	  w2Had.Boost(-boostWHad);
	  bHad.Boost( -boostWHad);
	  double cosChiStarTopHad = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
	  float helHad = fHelTopHadr->Eval( cosChiStarTopHad );

	  TLorentzVector bLep;
	  TLorentzVector w1Lep;
	  TLorentzVector w2Lep;
	  w1Lep.SetPxPyPzE( lvLepton.Px(), lvLepton.Py(), lvLepton.Pz(), lvLepton.E());
	  w2Lep.SetPxPyPzE( lvNeutrino.Px(), lvNeutrino.Py(), lvNeutrino.Pz(), lvNeutrino.E());
	  bLep.SetPxPyPzE(  lvLepB.Px(), lvLepB.Py(), lvLepB.Pz(), lvLepB.E());
	  TVector3 boostWLep(lvLepW.Px()/lvLepW.E(), lvLepW.Py()/lvLepW.E(), lvLepW.Pz()/lvLepW.E());
	  w1Lep.Boost(-boostWLep);
	  w2Lep.Boost(-boostWLep);
	  bLep.Boost( -boostWLep);
	  double cosChiStarTopLep = TMath::Cos(w1Lep.Vect().Angle( -bLep.Vect() ));
	  float helLep = fHelTopLept->Eval( cosChiStarTopLep );

	  if( (useKin || useBtag) && helHad>0. && helLep>0. ){
	    logp += (-TMath::Log(helHad));
	    logp += (-TMath::Log(helLep));
	  }
	  else if( !(useKin || useBtag) && helHad>0. && helLep>0.){
	    logp = (-TMath::Log(helHad));
	    logp += (-TMath::Log(helLep));
	  }
	  else
	    logp = 999;
	}


	if(logp<minLogP){
	  minLogP = logp;
	  bestRank.clear();
	  bestRank.push_back( indexWdau1Cand);
	  bestRank.push_back( indexWdau2Cand);
	  bestRank.push_back( indexbCand);
	  bestRank.push_back( indexbbarCand);
	}

	if(verbose) cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << "] => log(p) = " << logp << endl;

      }
      while ( std::next_permutation( indices, indices+myJetsFilt.size()  ) );

      if(verbose || true){
	cout << "Number of permutation:  " << counterPer << endl;
	cout << "Smallest logLikelihood: " << minLogP
	     << " [" << bestRank[0] << "," <<  bestRank[1] << "," << bestRank[2] << "," <<  bestRank[3] << "]"
	     << endl;
      }


      //int matchB    = 0;
      //int matchBbar = 0;
      //for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
      //if( minLogP>10 ) continue;
      //if( k==counterPer[0] || k==counterPer[1] ||  k==counterPer[2]  |  k==counterPer[3]   ) continue;
      //bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.3;
      //bool isBbar = deltaR( myJetsFilt[k], genBbarLV) < 0.3;
      //if( isB)     matchB++;
      //if( isBbar)  matchBbar++;
      //}
      

      delete fitter;
    } // event loop

    mySamples->GetFile( currentName )->Close();
    
  } // samples loop


  return 0;

}
