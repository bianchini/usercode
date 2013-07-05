#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMatrixT.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/AllIntegrationTypes.h"

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


#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooBernstein.h"
#include "RooNDKeysPdf.h"
#include "RooChiSquarePdf.h"

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/CalcME.h"

#define GENJETDR 0.3
#define VERBOSE  false
#define PI TMath::Pi()

using namespace std;
using namespace RooFit;

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
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;




pair<double,double> getMaxValue( TH1F* hMassProb){

  float est  = -99;
  float prob = -99;
  float a = -99;
  float b = -99;
  float c = -99;
  
  int maxBin = hMassProb->GetMaximumBin();
  
  if( maxBin<  hMassProb->GetNbinsX() && maxBin>1 ){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
    
    TMatrixD A(3,3);
    const double elements[9] = 
      { 1, xD,  xD*xD,
	1, xC,  xC*xC,
	1, xU,  xU*xU 
      };
    A.SetMatrixArray(elements, "C");
    TMatrixD AInv(3,3);
    double det;
    AInv = A.Invert(&det);
    
    TMatrixD Y(3,1);
    const double yPos[3] = 
      { yD, yC, yU
      };
    Y.SetMatrixArray(yPos, "C");
    
    TMatrixD C(3,1);
    const double dummy[3] = 
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est  = -b/2/a ;
    prob = a*est*est + b*est + c;
  }
  else if(maxBin== hMassProb->GetNbinsX()){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-2);
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-2);
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    
    TMatrixD A(3,3);
    const double elements[9] = 
	    { 1, xD,  xD*xD,
	      1, xC,  xC*xC,
	      1, xU,  xU*xU 
	    };
    A.SetMatrixArray(elements, "C");
    TMatrixD AInv(3,3);
    double det;
    AInv = A.Invert(&det);
    
    TMatrixD Y(3,1);
    const double yPos[3] = 
      { yD, yC, yU
      };
    Y.SetMatrixArray(yPos, "C");
    
    TMatrixD C(3,1);
    const double dummy[3] = 
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est = -b/2/a ;
    prob = a*est*est + b*est + c;	  
  }
  else if(maxBin==1){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+2);
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+2);
    
    TMatrixD A(3,3);
    const double elements[9] = 
      { 1, xD,  xD*xD,
	1, xC,  xC*xC,
	1, xU,  xU*xU 
	    };
    A.SetMatrixArray(elements, "C");
    TMatrixD AInv(3,3);
    double det;
    AInv = A.Invert(&det);
    
    TMatrixD Y(3,1);
    const double yPos[3] = 
      { yD, yC, yU
      };
    Y.SetMatrixArray(yPos, "C");
    
    TMatrixD C(3,1);
    const double dummy[3] = 
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est = -b/2/a ;
    prob = a*est*est + b*est + c;
  }
  else{
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
  }

  return make_pair(est, prob);

}


int main(int argc, const char* argv[])
{

  std::cout << "MEValidator" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName  ( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile   ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering     ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF     ( in.getParameter<std::string>  ("pathToTF"   ) );
  int vegasPoints          ( in.getParameter<int>   ("vegasPoints"   ) );
  bool verbose             ( in.getParameter<bool>  ("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met")        );

  float pertW1    ( in.getParameter<double> ("pertW1")     );
  float pertW2    ( in.getParameter<double> ("pertW2")     );
  float pertBHad  ( in.getParameter<double> ("pertBHad")   );
  float pertBLep  ( in.getParameter<double> ("pertBLep")   );
  float pertB1    ( in.getParameter<double> ("pertW1")     );
  float pertB2    ( in.getParameter<double> ("pertW2")     );

  float scaleH  ( in.getParameter<double> ("scaleH")   );
  float scaleL  ( in.getParameter<double> ("scaleL")   );
  float scaleMET( in.getParameter<double> ("scaleMET") );

  int   useME     ( in.getParameter<int>    ("useME")      );
  int   useJac    ( in.getParameter<int>    ("useJac")     );
  int   useMET    ( in.getParameter<int>    ("useMET")     );
  int   useTF     ( in.getParameter<int>    ("useTF")      );
  int   usePDF    ( in.getParameter<int>    ("usePDF")     );

  int   doParton      ( in.getParameter<int>    ("doParton")     );
  int   doSmear       ( in.getParameter<int>    ("doSmear")     );
  int   doMassScan    ( in.getParameter<int>    ("doMassScan")     );
  int   doPermutations( in.getParameter<int>    ("doPermutations") );

  int   printP4         ( in.getParameter<int>    ("printP4")          );
  int   mode            ( in.getUntrackedParameter<int>    ("mode", 0) );
  int   norm            ( in.getUntrackedParameter<int>    ("norm", 0) );

  vector<double> masses       ( in.getParameter<vector<double> >  ("masses")        );
  vector<string> functions    ( in.getParameter<vector<string> >  ("functions")     );

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];

  gSystem->Exec(("rm "+outFileName).c_str());

  TStopwatch* clock = new TStopwatch();
  TRandom3*   ran   = new TRandom3();
  ran->SetSeed(4321);


  const int nMassPoints  = masses.size();
  double mH[nMassPoints];
  for( unsigned int m = 0; m < masses.size() ; m++)
    mH[m] = masses[m];

  //{ /*50,  55,  */60,  65,  70,  75,  
  // 80 , 85,  90, 95, 100, 105, 110, 
  // 115, 120, 125, 130, 135, 140, 145, 150, 155, 
  // 160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250 
  //}; 


  ///////////////////////////////////////////////////////


  TH1F*  hMass         = new TH1F("hMass","",        nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMassProb     = new TH1F("hMassProb","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);

  TH1F*  hBestMass_gen     = new TH1F("hBestMass_gen","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMass_rec     = new TH1F("hBestMass_rec","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb_gen = new TH1F("hBestMassProb_gen","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb_rec = new TH1F("hBestMassProb_rec","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);

  TTree* tree  = new TTree("tree","");

  //TTree* tGenGoodComb  = new TTree("tree","");
  //TTree* tGenAllComb   = new TTree("tree","");
  //TTree* tRecGoodComb  = new TTree("tree","");
  //TTree* tRecAllComb   = new TTree("tree","");
 
  int counter_;

  float probGG_;
  float probAtSgnGG_;
  float massGG_;
  float evalCpuGG_,evalReaGG_;

  float probRG_;
  float probAtSgnRG_;
  float massRG_;
  float evalCpuRG_,evalReaRG_;

  float probGA_;
  float probAtSgnGA_;
  float probAtGoodGA_;
  float massGA_;
  float evalCpuGA_,evalReaGA_;

  float probRA_;
  float probAtSgnRA_;
  float probAtGoodRA_;
  float massRA_;
  float evalCpuRA_,evalReaRA_;

  tree->Branch("counter",                  &counter_,      "counter/I");

  tree->Branch("prob_gen_good",            &probGG_,      "prob_gen_good/F");
  tree->Branch("probAtSgn_gen_good",       &probAtSgnGG_, "probAtSgn_gen_good/F");
  tree->Branch("mass_gen_good",            &massGG_,      "mass_gen_good/F");
  tree->Branch("evalCpu_gen_good",         &evalCpuGG_,   "evalCpu_gen_good/F");
  tree->Branch("evalRea_gen_good",         &evalReaGG_,   "evalRea_gen_good/F");

  tree->Branch("prob_rec_good",            &probRG_,      "prob_rec_good/F");
  tree->Branch("probAtSgn_rec_good",       &probAtSgnRG_, "probAtSgn_rec_good/F");
  tree->Branch("mass_rec_good",            &massRG_,      "mass_rec_good/F");
  tree->Branch("evalCpu_rec_good",         &evalCpuRG_,   "evalCpu_rec_good/F");
  tree->Branch("evalRea_rec_good",         &evalReaRG_,   "evalRea_rec_good/F");

  tree->Branch("prob_gen_all",            &probGA_,      "prob_gen_all/F");
  tree->Branch("probAtSgn_gen_all",       &probAtSgnGA_, "probAtSgn_gen_all/F");
  tree->Branch("probAtGood_gen_all",      &probAtGoodGA_, "probAtGood_gen_all/F");
  tree->Branch("mass_gen_all",            &massGA_,      "mass_gen_all/F");
  tree->Branch("evalCpu_gen_all",         &evalCpuGA_,   "evalCpu_gen_all/F");
  tree->Branch("evalRea_gen_all",         &evalReaGA_,   "evalRea_gen_all/F");

  tree->Branch("prob_rec_all",            &probRA_,      "prob_rec_all/F");
  tree->Branch("probAtSgn_rec_all",       &probAtSgnRA_, "probAtSgn_rec_all/F");
  tree->Branch("probAtGood_rec_all",      &probAtGoodRA_, "probAtGood_rec_all/F");
  tree->Branch("mass_rec_all",            &massRA_,      "mass_rec_all/F");
  tree->Branch("evalCpu_rec_all",         &evalCpuRA_,   "evalCpu_rec_all/F");
  tree->Branch("evalRea_rec_all",         &evalReaRA_,   "evalRea_rec_all/F");

  ///////////////////////////////////////////////////////



  int par;
  switch(mode){
  case 0:
    par = 4;
    break;
  case 1:
    par = 6;
    break;
  case 2:
    par = 19;
    break;
  case 3:
    par = 19;
    break;
  default:
    par = 4;
    break;
  }


  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , par, int(verbose));
  if(mode == 0)
    meIntegrator->setIntType( MEIntegratorNew::SL2wj );
  else if(mode == 1)
    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
  else if(mode == 2)
    meIntegrator->setIntType( MEIntegratorNew::SLXSec );
  else if(mode == 3)
    meIntegrator->setIntType( MEIntegratorNew::SLAcc );
  else{
    cout << "Unsupported mode... exit" << endl;
    delete meIntegrator;
    return 0;
  }
  if( norm == 0)
    meIntegrator->setWeightNorm( MEIntegratorNew::None );
  else if( norm==1 )
    meIntegrator->setWeightNorm( MEIntegratorNew::xSec );
  else if( norm==2 )
    meIntegrator->setWeightNorm( MEIntegratorNew::Acc );
  else{
    cout << "Unsupported normalization... exit" << endl;
    delete meIntegrator;
    return 0;
  }
  meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
				 TString(functions[1].c_str()),  
				 TString(functions[2].c_str()));

  meIntegrator->setTopMass( 174.3 , 80.19);
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);


 ////////////////////////////////////////////////////////


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


  ////////////////////////////////////////////////////////



  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;

    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;
    Float_t vLepton_charge [99];

    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    

    int counter = 0;
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
  
      TLorentzVector TOPLEP;
      TLorentzVector TOPHAD;
      TLorentzVector TOPHADW1;
      TLorentzVector TOPHADW2;
      TLorentzVector TOPHADB;
      TLorentzVector TOPLEPW1;
      TLorentzVector TOPLEPW2;
      TLorentzVector TOPLEPB;
      TLorentzVector HIGGS;

      float neutCosTheta;

      bool properEvent = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( topW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  neutCosTheta = TMath::Cos( topW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{
	properEvent=false;
      }


      if(mode==0)
	properEvent = ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
			TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
			TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			);
      else if(mode==1)
	properEvent = (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			 (TOPHADW1.Pt() < 30 || TMath::Abs(TOPHADW1.Eta()) >2.5) &&
			 TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
			 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			 ) ||
		       ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			 TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
			 (TOPHADW2.Pt() < 30 || TMath::Abs(TOPHADW2.Eta()) >2.5) &&
			 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			 ) );
      else{ }	
      

      if(!properEvent) continue;
      counter++;

      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      counter_ = counter;

      vector<TLorentzVector> jets;
      jets.push_back( TOPLEPW1  );  
      jets.push_back( TOPLEPW2  );  
      jets.push_back( TOPLEPB   );
      jets.push_back( TOPHADW1  );  
      jets.push_back( TOPHADW2  );  
      jets.push_back( TOPHADB   );
      jets.push_back( genBLV    );
      jets.push_back( genBbarLV );

      if(mode==1){
	TLorentzVector w1 = jets[3];
	TLorentzVector w2 = jets[4];
	if( (w1.Pt()<30 || TMath::Abs(w1.Eta())>2.5) && (w2.Pt()>30 && TMath::Abs(w2.Eta())<2.5)){
	  jets[3] = w2;
	  jets[4] = w1;
	}
	else if( (w2.Pt()<30 || TMath::Abs(w2.Eta())>2.5) && (w1.Pt()>30 && TMath::Abs(w1.Eta())<2.5)){
	  jets[3] = w1;
	  jets[4] = w2;
	}
	else{
	  cout << "Inconsistentcy of mode and jet selections" << endl;
	  delete meIntegrator;
	  return 1;
	}
      }

      
      if(printP4){
	cout << "******* INPUT ******" << endl;      
	cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	cout << "lep:  jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta() << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << ")" << endl;
	cout << "met:  jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta() << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << ")" << endl;
	cout << "blep: jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << "," <<  TOPLEPB.Eta() << "," << TOPLEPB.Phi() << "," << TOPLEPB.M() << ")" << endl;
	cout << "w1:   jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta() << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	cout << "w2:   jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta() << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	cout << "bhad: jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << "," <<  TOPHADB.Eta() << "," << TOPHADB.Phi() << "," << TOPHADB.M() << ")" << endl;
	cout << "h1:   jet7.SetPtEtaPhiM(" << genBLV.Pt() << "," <<  genBLV.Eta() << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	cout << "h2:   jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << "," <<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	cout << "Top Lep mass = " << (TOPLEPW1+TOPLEPW2+TOPLEPB).M() << " <=> neutrino eta=0!!!" << endl;
	cout << "Top Had mass = " << (TOPHADW1+TOPHADW2+TOPHADB).M() << endl;
	cout << "Higgs mass = "   << (genBLV+genBbarLV).M() << endl;
	cout << "Neut cosTheta = " << neutCosTheta << endl;
	}


      meIntegrator->setJets(&jets);
      meIntegrator->setSumEt( METtype1p2corr.sumet );
      meIntegrator->setMEtCov(-99,-99,0);
      meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );

      probGG_      = 0;
      probAtSgnGG_ = 0;
      massGG_      = 0;
      evalCpuGG_   = 0;
      evalReaGG_   = 0;
      probGA_      = 0;
      massGA_      = 0;
      evalCpuGA_   = 0;
      evalReaGA_   = 0;
      probAtSgnGA_ = 0;
      probAtGoodGA_= 0;
      probRG_      = 0;
      probAtSgnRG_ = 0;
      massRG_      = 0;
      evalCpuRG_   = 0;
      evalReaRG_   = 0;	  	
      probRA_      = 0;
      massRA_      = 0;
      evalCpuRA_   = 0;
      evalReaRA_   = 0;
      probAtSgnRA_ = 0;
      probAtGoodRA_= 0;



      if( doParton){

	if( doMassScan ){
	
	  meIntegrator->initVersors(1);
		
	  pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	  pair<double, double> range_x1 =  make_pair(-1,1);
	  pair<double, double> range_x2 =  make_pair(-PI,PI);
	  pair<double, double> range_x3 =  make_pair(-1,1);
	  pair<double, double> range_x4 = (meIntegrator->getNuPhiCI(0.95));
	  pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	  
	  double x0L = range_x0.first;
	  double x0U = range_x0.second;
	  double x1L = range_x1.first;
	  double x1U = range_x1.second;
	  double x2L = range_x2.first;
	  double x2U = range_x2.second;
	  double x3L = range_x3.first;
	  double x3U = range_x3.second;
	  double x4L = range_x4.first;
	  double x4U = range_x4.second;
	  double x5L = range_x5.first;
	  double x5U = range_x5.second;
	
	  double xLmode0[4] = {x0L, x3L, x4L, x5L};
	  double xUmode0[4] = {x0U, x3U, x4U, x5U};
	  
	  double xLmode1[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	  double xUmode1[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	  
	  if( printP4 ){
	    cout << "Integration range: " << endl;
	    for(unsigned int k = 0; k < 4; k++){
	      cout << "Var " << k << ": [" << xLmode0[k] << "," <<  xUmode0[k] << "]" << endl;
	    }	   
	  }
	
	  clock->Start();
	  
	  hMass->Reset();
	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
	    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
	    ig2.SetFunction(toIntegrate);
	    meIntegrator->SetPar(par);
	  
	    double p = 0.;
	    if(mode==0)
	      p = ig2.Integral(xLmode0, xUmode0);
	    else if(mode==1)
	      p = ig2.Integral(xLmode1, xUmode1);
	    else{ }
	    
	    clock->Stop();
	    evalCpuGG_ += clock->CpuTime();
	    evalReaGG_ += clock->RealTime();
	    if(printP4) cout << "Mass " << mH[m] << " => prob  = " << p << endl;
	    hMass->Fill(mH[m], p);
	    if( mH[m]<met+0.5 && mH[m]>met-0.5) probAtSgnGG_ += p;
	    clock->Reset();
	  }
	  evalCpuGG_ /= nMassPoints;
	  evalReaGG_ /= nMassPoints;
	  
	  pair<double,double> bestMass = getMaxValue(hMass);
	  hBestMass_gen->Fill( bestMass.first );
	  massGG_ = bestMass.first;
	  probGG_ = bestMass.second;	 
	}
	
	if(doPermutations){
	
	  hMassProb->Reset();
	
	  double pTotAtSgn = 0.;
	  int combScanned = 0;
	
	  int permutations[12] = 
	    {234567, 534267,     // CORRECT 
	     634725, 734625,     // ALL WRONG
	     534762, 734562,     // ONE WRONG
	     234765, 734265,     // ONE WRONG
	     634572, 534672,     // ONE WRONG
	     234675, 634275      // ONE WRONG
	    };
	  
	  
	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	    
	    for(unsigned int pos = 0; pos < 12; pos++){
	      meIntegrator->initVersors( permutations[pos] );
	      
	      if( !(meIntegrator->compatibilityCheck(0.95)) ){
		if(printP4) cout << "M=" << mH[m] << " incompatible with tested combination " << permutations[pos] << endl;		  
		continue;
	      }
	      
	      combScanned++;
	      
	      pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	      pair<double, double> range_x1 =  make_pair(-1,1);
	      pair<double, double> range_x2 =  make_pair(-PI,PI);
	      pair<double, double> range_x3 =  make_pair(-1,1);
	      pair<double, double> range_x4 = (meIntegrator->getNuPhiCI(0.95));
	      pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	      
	      double x0L = range_x0.first;
	      double x0U = range_x0.second;
	      double x1L = range_x1.first;
	      double x1U = range_x1.second;
	      double x2L = range_x2.first;
	      double x2U = range_x2.second;
	      double x3L = range_x3.first;
	      double x3U = range_x3.second;
	      double x4L = range_x4.first;
	      double x4U = range_x4.second;
	      double x5L = range_x5.first;
	      double x5U = range_x5.second;
	    
	      double xLmode0[4] = {x0L, x3L, x4L, x5L};
	      double xUmode0[4] = {x0U, x3U, x4U, x5U};
	      
	      double xLmode1[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	      double xUmode1[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	      
	      if( printP4 ){
		cout << "Integration range: " << endl;
		for(unsigned int k = 0; k < 4; k++){
		  cout << "Var " << k << ": [" << xLmode0[k] << "," <<  xUmode0[k] << "]" << endl;
		}	   
	      }
	    
	      clock->Start();

	      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
	      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
	      ig2.SetFunction(toIntegrate);
	      meIntegrator->SetPar(par);
	      
	      double p = 0.;
	      if(mode==0) 
		p = ig2.Integral(xLmode0, xUmode0);
	      else if(mode==1)
		p = ig2.Integral(xLmode1, xUmode1);
	      else{ }
	      
	      if( pos==0 ) pTotAtSgn += p;
	      
	      clock->Stop();
	      evalCpuGA_ += clock->CpuTime();
	      evalReaGA_ += clock->RealTime();
	      
	      clock->Reset();

	      if( mH[m]<met+0.5 && mH[m]>met-0.5) probAtSgnGA_ += p;
	      
	      float old = hMassProb->GetBinContent( hMassProb->FindBin(mH[m]) );
	      hMassProb->SetBinContent( hMassProb->FindBin(mH[m]), old+p );
	    } // permutations
	  }
	  
	  evalCpuGA_ /= combScanned;
	  evalReaGA_ /= combScanned;
	  probAtGoodGA_ = pTotAtSgn;
	
	  pair<double,double> bestMass = getMaxValue(hMassProb);
	  hBestMassProb_gen->Fill( bestMass.first );
	  massGA_ = bestMass.first;
	  probGA_ = bestMass.second;

	}

      }


      if( doSmear ){

	meIntegrator->initVersors(1);

	if( meIntegrator->smearByTF(30.) ){
	  
	  if( doMassScan ){

	    pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	    pair<double, double> range_x1 =  make_pair(-1,1);
	    pair<double, double> range_x2 =  make_pair(-PI,PI);
	    pair<double, double> range_x3 =  make_pair(-1,1);
	    pair<double, double> range_x4 = (meIntegrator->getNuPhiCI(0.95));
	    pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	    
	    double x0L = range_x0.first;
	    double x0U = range_x0.second;
	    double x1L = range_x1.first;
	    double x1U = range_x1.second;
	    double x2L = range_x2.first;
	    double x2U = range_x2.second;
	    double x3L = range_x3.first;
	    double x3U = range_x3.second;
	    double x4L = range_x4.first;
	    double x4U = range_x4.second;
	    double x5L = range_x5.first;
	    double x5U = range_x5.second;
	    
	    double xLmode0[4] = {x0L, x3L, x4L, x5L};
	    double xUmode0[4] = {x0U, x3U, x4U, x5U};
	    
	    double xLmode1[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	    double xUmode1[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	    
	    if( printP4 ){
	      cout << "Integration range: " << endl;
	      for(unsigned int k = 0; k < 4; k++){
		cout << "Var " << k << ": [" << xLmode0[k] << "," <<  xUmode0[k] << endl;
	      }	   
	    }
	    
	    clock->Start();
	    
	    hMass->Reset();
	    for(int m = 0; m < nMassPoints ; m++){
	      meIntegrator->setMass( mH[m] );
	      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
	      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
	      ig2.SetFunction(toIntegrate);
	      meIntegrator->SetPar(par);
	      
	      double p = 0.;
	      if(mode==0)
		p = ig2.Integral(xLmode0, xUmode0);
	      else if(mode==1)
		p = ig2.Integral(xLmode1, xUmode1);
	      else{ }
	      
	      clock->Stop();
	      evalCpuRG_ += clock->CpuTime();
	      evalReaRG_ += clock->RealTime();
	      if(printP4) cout << "Mass " << mH[m] << " => prob  = " << p << endl;
	      hMass->Fill(mH[m], p);
	      if( mH[m]<met+0.5 && mH[m]>met-0.5) probAtSgnRG_ += p;
	      clock->Reset();
	    }
	    evalCpuRG_ /= nMassPoints;
	    evalReaRG_ /= nMassPoints;
	    
	    pair<double,double> bestMass = getMaxValue(hMass);
	    hBestMass_rec->Fill( bestMass.first );
	    massRG_ = bestMass.first;
	    probRG_ = bestMass.second;	 
	  }
	
	  
	  if(doPermutations){
	    
	    hMassProb->Reset();
	    
	    double pTotAtSgn = 0.;
	    int combScanned = 0;
	    
	    int permutations[12] = 
	      {234567, 534267,      // CORRECT 
	       634725, 734625,     // ALL WRONG
	       534762, 734562,     // ONE WRONG
	       234765, 734265,     // ONE WRONG
	       634572, 534672,     // ONE WRONG
	       234675, 634275      // ONE WRONG
	      };
	    
	    
	    for(int m = 0; m < nMassPoints ; m++){
	      meIntegrator->setMass( mH[m] );
	      
	      for(unsigned int pos = 0; pos < 12; pos++){
		meIntegrator->initVersors( permutations[pos] );
		
		if( !(meIntegrator->compatibilityCheck(0.95)) ){
		  if(printP4) cout << "M=" << mH[m] << " incompatible with tested combination " << permutations[pos] << endl;		  
		  continue;
		}
		
		combScanned++;
		
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 = (meIntegrator->getNuPhiCI(0.95));
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		
		double x0L = range_x0.first;
		double x0U = range_x0.second;
		double x1L = range_x1.first;
		double x1U = range_x1.second;
		double x2L = range_x2.first;
		double x2U = range_x2.second;
		double x3L = range_x3.first;
		double x3U = range_x3.second;
		double x4L = range_x4.first;
		double x4U = range_x4.second;
		double x5L = range_x5.first;
		double x5U = range_x5.second;
		
		double xLmode0[4] = {x0L, x3L, x4L, x5L};
		double xUmode0[4] = {x0U, x3U, x4U, x5U};
		
		double xLmode1[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode1[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
		
		if( printP4 ){
		  cout << "Integration range: " << endl;
		  for(unsigned int k = 0; k < 4; k++){
		    cout << "Var " << k << ": [" << xLmode0[k] << "," <<  xUmode0[k] << endl;
		  }	   
		}
		
		clock->Start();
		
		ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
		ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
		ig2.SetFunction(toIntegrate);
		meIntegrator->SetPar(par);
		
		double p = 0.;
		if(mode==0) 
		  p = ig2.Integral(xLmode0, xUmode0);
		else if(mode==1)
		p = ig2.Integral(xLmode1, xUmode1);
		else{ }
		
		if( pos==0 ) pTotAtSgn += p;
		
		clock->Stop();
		evalCpuRA_ += clock->CpuTime();
		evalReaRA_ += clock->RealTime();
		
		clock->Reset();
		
		if( mH[m]<met+0.5 && mH[m]>met-0.5) probAtSgnGA_ += p;
		
		float old = hMassProb->GetBinContent( hMassProb->FindBin(mH[m]) );
		hMassProb->SetBinContent( hMassProb->FindBin(mH[m]), old+p );
	      } // permutations
	    }
	    
	    evalCpuRA_ /= combScanned;
	    evalReaRA_ /= combScanned;
	    probAtGoodRA_ = pTotAtSgn;
	    
	    pair<double,double> bestMass = getMaxValue(hMassProb);
	    hBestMassProb_rec->Fill( bestMass.first );
	    massRA_ = bestMass.first;
	    probRA_ = bestMass.second;
	    
	  }
	}

      }

      tree->Fill();      
    } // nentries
    
  } // samples

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
  hBestMass_gen->Write("hBestMass_gen",        TObject::kOverwrite);
  hBestMassProb_gen->Write("hBestMassProb_gen",TObject::kOverwrite);
  hBestMass_rec->Write("hBestMass_rec",        TObject::kOverwrite);
  hBestMassProb_rec->Write("hBestMassProb_rec",TObject::kOverwrite);

  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete ran; delete clock;
  cout << "Finished!!!" << endl;
  
  return 0;
}

