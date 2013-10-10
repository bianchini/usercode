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
#include "THStack.h"
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
#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"

 
#define GENJETDR  0.3
#define MAX_REEVAL_TRIES 3
#define VERBOSE  false
#define PI TMath::Pi()

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
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;
 

struct sorterByPt {
  bool operator() (float i,float j) const { return (i>j);}
};


typedef struct
{
  TLorentzVector p4;
  float csv;
} JetObservable;


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

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  std::cout << "MEAnalysis" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ CONFIGURATION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  // SAMPLES
  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering   ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF   ( in.getParameter<std::string>  ("pathToTF"   ) );
  std::string pathToCP   ( in.getParameter<std::string>  ("pathToCP"   ) );
  bool   verbose         ( in.getParameter<bool>  ("verbose") );
  
  // PARAMETERS
  double lumi            ( in.getParameter<double>("lumi") );
  float  MH              ( in.getUntrackedParameter<double> ("MH",     125.));
  float  MT              ( in.getUntrackedParameter<double> ("MT",    174.3));
  float  MW              ( in.getUntrackedParameter<double> ("MW",    80.19));
  float  MwL             ( in.getUntrackedParameter<double> ("MwL",      60));
  float  MwH             ( in.getUntrackedParameter<double> ("MwH",     100));
  //float  MhL             ( in.getUntrackedParameter<double> ("MhL",     100));
  //float  MhH             ( in.getUntrackedParameter<double> ("MhH",     140));
  double maxChi2_        ( in.getUntrackedParameter<double> ("maxChi2", 2.5));
  vector<double> massesH ( in.getParameter<vector<double> > ("massesH"));
  vector<double> massesT ( in.getParameter<vector<double> > ("massesT"));

  // FLAGS
  int   switchoffOL      ( in.getUntrackedParameter<int>    ("switchoffOL",     0));
  int   speedup          ( in.getUntrackedParameter<int>    ("speedup",         0));
  int   doubleGaussianB  ( in.getUntrackedParameter<int>    ("doubleGaussianB", 1));
  int   useBtag          ( in.getUntrackedParameter<int>    ("useBtag",         0));
  int   doType0          ( in.getUntrackedParameter<int>    ("doType0", 0));
  int   doType1          ( in.getUntrackedParameter<int>    ("doType1", 0));
  int   doType2          ( in.getUntrackedParameter<int>    ("doType2", 0));
  int   doType3          ( in.getUntrackedParameter<int>    ("doType3", 0));
  //int   doType4          ( in.getUntrackedParameter<int>    ("doType4", 0));
  int   doType6          ( in.getUntrackedParameter<int>    ("doType6", 0));
  int   doType7          ( in.getUntrackedParameter<int>    ("doType7", 0));
  int   useME            ( in.getParameter<int>             ("useME")     );
  int   useJac           ( in.getParameter<int>             ("useJac")    );
  int   useMET           ( in.getParameter<int>             ("useMET")    );
  int   useTF            ( in.getParameter<int>             ("useTF")     );
  int   usePDF           ( in.getParameter<int>             ("usePDF")    );
  //int   printP4          ( in.getParameter<int>             ("printP4")   );
  int   norm             ( in.getUntrackedParameter<int>    ("norm",      0));
  int   hypo             ( in.getUntrackedParameter<int>    ("hypo",      0));
  int   SoB              ( in.getUntrackedParameter<int>    ("SoB",       1));
  int   doCSVup          ( in.getUntrackedParameter<int>    ("doCSVup",   0));
  int   doCSVdown        ( in.getUntrackedParameter<int>    ("doCSVdown", 0));
  int   doJECup          ( in.getUntrackedParameter<int>    ("doJECup",   0));
  int   doJECdown        ( in.getUntrackedParameter<int>    ("doJECdown", 0));
  int   fixNumEvJob      ( in.getUntrackedParameter<int>    ("fixNumEvJob",1));
  vector<string> functions( in.getParameter<vector<string> >("functions"));
  vector<int>    evLimits( in.getParameter<vector<int> >   ("evLimits"));



  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ INITIALIZE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

  // upper and lower bounds to be processed
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];

  TStopwatch* clock = new TStopwatch();

  // file with b-tag pdf
  TFile* fCP = TFile::Open(pathToCP.c_str(),"READ");

  // b-tag pdf for b-quark ('b'), c-quark ('c'), and light jets ('l')
  map<string,TH1F*> btagger; 
  if( useBtag && fCP!=0 ){
    btagger["b_Bin0"] = fCP->Get("csv_b_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin0__csvReco") : 0;
    btagger["b_Bin1"] = fCP->Get("csv_b_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin1__csvReco") : 0;
    btagger["c_Bin0"] = fCP->Get("csv_c_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin0__csvReco") : 0;
    btagger["c_Bin1"] = fCP->Get("csv_c_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin1__csvReco") : 0;
    btagger["l_Bin0"] = fCP->Get("csv_l_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin0__csvReco") : 0;
    btagger["l_Bin1"] = fCP->Get("csv_l_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin1__csvReco") : 0;
  }
  else if( useBtag && fCP!=0 ){
    cout << "Cound not find " << pathToCP << ": exit" << endl;
    return 0;
  }

  // Higgs mass values for scan
  const int nMassPoints  = massesH.size();
  double mH[nMassPoints];
  for( unsigned int m = 0; m < massesH.size() ; m++)
    mH[m] = massesH[m];

  // Top mass values for scan
  const int nTopMassPoints  = massesT.size();
  double mT[nTopMassPoints]; 
  for( unsigned int m = 0; m < massesT.size() ; m++)
    mT[m] = massesT[m];


  int permutations_SL2wj[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };

  int permutations_SL1wj[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };

  int permutations_DL[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };
  

  // configure MEIntegrator
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , 4 , int(verbose));
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

  // set normalization formulas ( not used if norm==0 )
  meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
				 TString(functions[1].c_str()),  
				 TString(functions[2].c_str()),
				 TString(functions[3].c_str()),  
				 TString(functions[4].c_str()),  
				 TString(functions[5].c_str())
				 );
  
  // initialize top and W mass
  meIntegrator->setTopMass( MT , MW );

  // configure ME calculation
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);

  // use double-gaussian for b quark energy TF
  meIntegrator->setUseRefinedTF(doubleGaussianB);

  // use nominal TF from ROOT file
  meIntegrator->initTFparameters(1.0,1.0,1.0,1.0, 1.0);
  if(switchoffOL){
    meIntegrator->switchOffOL(); 
    cout << "*** Switching off OpenLoops to speed-up the calculation ***" << endl;
  }



  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  // clean output file (if any)
  gSystem->Exec(("rm "+outFileName).c_str());

  // output file
  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");

  // total event counter for normalization
  TH1F*  hcounter = new TH1F("hcounter","",1,0,1);
  int events_     = 0;
  
  // per-event probability vs mass
  TH1F*  h_prob_vs_higgs_mass     = new TH1F("h_prob_vs_higgs_mass",   "", nMassPoints,    mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  h_prob_vs_top_mass       = new TH1F("h_prob_vs_top_mass",     "", nTopMassPoints, mT[0]-2.5, mT[nTopMassPoints-1]+2.5);
  TH1F*  h_prob_vs_top_mass_ttbb  = new TH1F("h_prob_vs_top_mass_ttbb","", nTopMassPoints, mT[0]-2.5, mT[nTopMassPoints-1]+2.5);
  TH1F*  h_prob_vs_top_mass_ttbj  = new TH1F("h_prob_vs_top_mass_ttbj","", nTopMassPoints, mT[0]-2.5, mT[nTopMassPoints-1]+2.5);
  TH1F*  h_prob_vs_top_mass_ttcc  = new TH1F("h_prob_vs_top_mass_ttcc","", nTopMassPoints, mT[0]-2.5, mT[nTopMassPoints-1]+2.5);
  TH1F*  h_prob_vs_top_mass_ttjj  = new TH1F("h_prob_vs_top_mass_ttjj","", nTopMassPoints, mT[0]-2.5, mT[nTopMassPoints-1]+2.5);

  // output tree
  TTree* tree  = new TTree("tree","");

  // counts how many events have been analyzed (to have fixed-size jobs)
  int counter_;
  // number of jet-quark permutations
  int nPermut_;
  // number of event interpretations
  int nInter_;
  // total number of jet-quark permutations (nPermut*nInter)
  int nTotPermut_;
  // number of matches to higgs quarks among tagged jets
  int matchesH_;
  // number of matches to W quarks among un-tagged jets
  int matchesW_;
  // number of matches to top quarks among tagged jets
  int matchesT_;
  // number of matches to higgs quarks among all jets
  int matchesHAll_;
  // number of matches to W quarks among all jets
  int matchesWAll_;
  // number of matches to top quarks among all jets
  int matchesTAll_;
  // count how many quarks from W decay overlap by dR<0.5
  int overlapLight_;
  // count how many b-quarks overlap by dR<0.5
  int overlapHeavy_;
  // integration type
  int type_;
  // num. of b-hadrons and c-quarks
  int nSimBs_, nC_, nCTop_;
  // num of b-hadrons inside jets
  int nMatchSimBs_;
  // type-dependent flags
  int flag_type0_;
  int flag_type1_;
  int flag_type2_;
  int flag_type3_;
  int flag_type4_;
  int flag_type6_;
  // event-wise probability (summed over permutations)
  float probAtSgn_;
  float probAtSgn_alt_;
  float probAtSgn_ttbb_;
  float probAtSgn_alt_ttbb_;
  float probAtSgn_alt_ttjj_;
  float probAtSgn_alt_ttbj_;
  float probAtSgn_alt_ttcc_;
  // best mass (from likelihood scan)
  float bestMH_;
  float bestMT_;
  float bestMTbb_;
  float bestMTbj_;
  float bestMTcc_;
  float bestMTjj_;
  // per-permutation probability
  float probAtSgn_permut_[999];
  float probAtSgn_alt_permut_[999];
  float probAtSgn_bb_permut_[999];
  float probAtSgn_bj_permut_[999];
  float probAtSgn_cc_permut_[999];
  float probAtSgn_jj_permut_[999];
  // event-dependent weight (for normalization)
  float weight_;
  // cpu time
  float time_;
  // event information
  EventInfo EVENT_;

  tree->Branch("counter",      &counter_,       "counter/I");
  tree->Branch("nPermut",      &nPermut_,       "nPermut/I");
  tree->Branch("nTotPermut",   &nTotPermut_,    "nTotPermut/I");
  tree->Branch("nInter",       &nInter_,        "nInter/I");  
  tree->Branch("matchesH",     &matchesH_,      "matchesH/I");
  tree->Branch("matchesW",     &matchesW_,      "matchesW/I");
  tree->Branch("matchesT",     &matchesT_,      "matchesT/I");
  tree->Branch("matchesHAll",  &matchesHAll_,   "matchesHAll/I");
  tree->Branch("matchesWAll",  &matchesWAll_,   "matchesWAll/I");
  tree->Branch("matchesTAll",  &matchesTAll_,   "matchesTAll/I");
  tree->Branch("overlapLight", &overlapLight_,  "overlapLight/I");
  tree->Branch("overlapHeavy", &overlapHeavy_,  "overlapHeavy/I");
  tree->Branch("type",         &type_,          "type/I");
  tree->Branch("nSimBs",       &nSimBs_,        "nSimBs/I");
  tree->Branch("nMatchSimBs",  &nMatchSimBs_,   "nMatchSimBs/I");
  tree->Branch("nC",           &nC_,            "nC/I");
  tree->Branch("nCTop",        &nCTop_,         "nCTop/I");
  tree->Branch("weight",       &weight_,        "weight/F");
  tree->Branch("time",         &time_,          "time/F");
  tree->Branch("flag_type0",   &flag_type0_,    "flag_type0/I");
  tree->Branch("flag_type1",   &flag_type1_,    "flag_type1/I");
  tree->Branch("flag_type2",   &flag_type2_,    "flag_type2/I");
  tree->Branch("flag_type3",   &flag_type3_,    "flag_type3/I");
  tree->Branch("flag_type4",   &flag_type4_,    "flag_type4/I");
  tree->Branch("flag_type6",   &flag_type6_,    "flag_type6/I");
  tree->Branch("bestMH",       &bestMH_,        "bestMH/F");
  tree->Branch("bestMT",       &bestMT_,        "bestMT/F");
  tree->Branch("bestMTbb",     &bestMTbb_,      "bestMTbb/F");
  tree->Branch("bestMTbj",     &bestMTbj_,      "bestMTbj/F");
  tree->Branch("bestMTcc",     &bestMTcc_,      "bestMTcc/F");
  tree->Branch("bestMTjj",     &bestMTjj_,      "bestMTjj/F");
  tree->Branch("EVENT",        &EVENT_,         "run/I:lumi/I:event/I:json/I");

  tree->Branch(Form("p_%d_all_s",     int(MH)),   &probAtSgn_,           Form("p_%d_all_s/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_b",     int(MH)),   &probAtSgn_alt_,       Form("p_%d_all_b/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_s_ttbb",int(MH)),   &probAtSgn_ttbb_,      Form("p_%d_all_s_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbb",int(MH)),   &probAtSgn_alt_ttbb_,  Form("p_%d_all_b_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttjj",int(MH)),   &probAtSgn_alt_ttjj_,  Form("p_%d_all_b_ttjj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbj",int(MH)),   &probAtSgn_alt_ttbj_,  Form("p_%d_all_b_ttbj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttcc",int(MH)),   &probAtSgn_alt_ttcc_,  Form("p_%d_all_b_ttcc/F",         int(MH)) );
  tree->Branch(Form("p_%d_s",         int(MH)),   probAtSgn_permut_,     Form("p_%d_s[nTotPermut]/F",  int(MH)) );
  tree->Branch(Form("p_%d_b",         int(MH)),   probAtSgn_alt_permut_, Form("p_%d_b[nTotPermut]/F",  int(MH)) );
  tree->Branch(Form("p_%d_bb",        int(MH)),   probAtSgn_bb_permut_,  Form("p_%d_bb[nTotPermut]/F", int(MH)) );
  tree->Branch(Form("p_%d_bj",        int(MH)),   probAtSgn_bj_permut_,  Form("p_%d_bj[nTotPermut]/F", int(MH)) );
  tree->Branch(Form("p_%d_cc",        int(MH)),   probAtSgn_cc_permut_,  Form("p_%d_cc[nTotPermut]/F", int(MH)) );
  tree->Branch(Form("p_%d_jj",        int(MH)),   probAtSgn_jj_permut_,  Form("p_%d_jj[nTotPermut]/F", int(MH)) );

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ OPEN FILES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
 

  // read input files  
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

  // loop over input files
  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    float scaleFactor        = mySamples->GetWeight(currentName);

    // variables to be used from the input files
    genParticleInfo genB, genBbar;
    genTopInfo      genTop, genTbar;
    metInfo         METtype1p2corr;
    EventInfo       EVENT;
    int             nvlep, nSimBs, nC, nCTop, nhJets, naJets;
    Int_t   vLepton_type      [999];
    Float_t vLepton_mass      [999];
    Float_t vLepton_pt        [999];
    Float_t vLepton_eta       [999];
    Float_t vLepton_phi       [999];
    Float_t vLepton_charge    [999];
    Float_t vLepton_pfCorrIso [999];
    Float_t hJet_pt           [999];
    Float_t hJet_eta          [999];
    Float_t hJet_phi          [999];
    Float_t hJet_e            [999];
    Float_t hJet_puJetIdL     [999];
    Float_t hJet_csv_nominal  [999];
    Float_t hJet_csv_upBC     [999];
    Float_t hJet_csv_downBC   [999];
    Float_t hJet_csv_upL      [999];
    Float_t hJet_csv_downL    [999];
    Float_t hJet_JECUnc       [999];
    Float_t hJet_genPt        [999];
    Float_t hJet_genEta       [999];
    Float_t hJet_genPhi       [999];
    Float_t aJet_pt           [999];
    Float_t aJet_eta          [999];
    Float_t aJet_phi          [999];
    Float_t aJet_e            [999];
    Float_t aJet_puJetIdL     [999];
    Float_t aJet_csv_nominal  [999];
    Float_t aJet_csv_upBC     [999];
    Float_t aJet_csv_downBC   [999];
    Float_t aJet_csv_upL      [999];
    Float_t aJet_csv_downL    [999];
    Float_t aJet_JECUnc       [999];
    Float_t aJet_genPt        [999];
    Float_t aJet_genEta       [999];
    Float_t aJet_genPhi       [999];
    float SimBsmass           [999];
    float SimBspt             [999];
    float SimBseta            [999];
    float SimBsphi            [999];
    //int nSvs;
    //float SvmassSv [999];
    //float Svpt     [999];
    //float Sveta    [999];
    //float Svphi    [999];
  
    currentTree->SetBranchAddress("EVENT",       &EVENT);
    currentTree->SetBranchAddress("nhJets",      &nhJets);
    currentTree->SetBranchAddress("naJets",      &naJets);
    currentTree->SetBranchAddress("nSimBs",      &nSimBs);
    currentTree->SetBranchAddress("nvlep",       &nvlep);
    currentTree->SetBranchAddress("nC",          &nC);
    currentTree->SetBranchAddress("nCTop",       &nCTop);
    currentTree->SetBranchAddress("genB",            &genB);
    currentTree->SetBranchAddress("genBbar",         &genBbar);
    currentTree->SetBranchAddress("genTop",          &genTop);
    currentTree->SetBranchAddress("genTbar",         &genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",  &METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_mass"  ,   vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,   vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,   vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,   vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_pfCorrIso",vLepton_pfCorrIso);
    currentTree->SetBranchAddress("vLepton_type",     vLepton_type);
    currentTree->SetBranchAddress("hJet_pt",          hJet_pt);    
    currentTree->SetBranchAddress("hJet_eta",         hJet_eta);    
    currentTree->SetBranchAddress("hJet_phi",         hJet_phi);    
    currentTree->SetBranchAddress("hJet_e",           hJet_e);    
    currentTree->SetBranchAddress("hJet_puJetIdL",    hJet_puJetIdL);
    currentTree->SetBranchAddress("hJet_csv_nominal", hJet_csv_nominal);
    currentTree->SetBranchAddress("hJet_csv_upBC",    hJet_csv_upBC);
    currentTree->SetBranchAddress("hJet_csv_downBC",  hJet_csv_downBC);
    currentTree->SetBranchAddress("hJet_csv_upL",     hJet_csv_upL);
    currentTree->SetBranchAddress("hJet_csv_downL",   hJet_csv_downL);
    currentTree->SetBranchAddress("hJet_JECUnc",      hJet_JECUnc);
    currentTree->SetBranchAddress("hJet_genPt",       hJet_genPt);
    currentTree->SetBranchAddress("hJet_genEta",      hJet_genEta);
    currentTree->SetBranchAddress("hJet_genPhi",      hJet_genPhi);
    currentTree->SetBranchAddress("aJet_pt",          aJet_pt);    
    currentTree->SetBranchAddress("aJet_eta",         aJet_eta);    
    currentTree->SetBranchAddress("aJet_phi",         aJet_phi);    
    currentTree->SetBranchAddress("aJet_e",           aJet_e);    
    currentTree->SetBranchAddress("aJet_puJetIdL",    hJet_puJetIdL);
    currentTree->SetBranchAddress("aJet_csv_nominal", aJet_csv_nominal);
    currentTree->SetBranchAddress("aJet_csv_upBC",    aJet_csv_upBC);
    currentTree->SetBranchAddress("aJet_csv_downBC",  aJet_csv_downBC);
    currentTree->SetBranchAddress("aJet_csv_upL",     aJet_csv_upL);
    currentTree->SetBranchAddress("aJet_csv_downL",   aJet_csv_downL);
    currentTree->SetBranchAddress("aJet_JECUnc",      aJet_JECUnc);
    currentTree->SetBranchAddress("aJet_genPt",       aJet_genPt);
    currentTree->SetBranchAddress("aJet_genEta",      aJet_genEta);
    currentTree->SetBranchAddress("aJet_genPhi",      aJet_genPhi);
    currentTree->SetBranchAddress("SimBs_mass",   SimBsmass);
    currentTree->SetBranchAddress("SimBs_pt",     SimBspt);
    currentTree->SetBranchAddress("SimBs_eta",    SimBseta);
    currentTree->SetBranchAddress("SimBs_phi",    SimBsphi);
    //currentTree->SetBranchAddress("nSvs",        &nSvs);
    //currentTree->SetBranchAddress("Sv_massSv",   SvmassSv);
    //currentTree->SetBranchAddress("Sv_pt",       Svpt);
    //currentTree->SetBranchAddress("Sv_eta",      Sveta);
    //currentTree->SetBranchAddress("Sv_phi",      Svphi);

 
    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
    /* @@@@@@@@@@@@@@@@@@@@@@@@@ EVENT LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
    
    // loop over entries
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();

    if( evHigh<0 ) evHigh = nentries;
    for (Long64_t i = 0; i < nentries ; i++){
      
      // if fixed-size job and above upper bound, continue...
      if(counter>evHigh && fixNumEvJob) continue;

      // otherwise, if outside the event window, continue...
      if(!fixNumEvJob && !(i>=evLow && i<evHigh) ) continue;
      events_++;

      if(i%5000==0) cout << i << endl;

      // read event...
      currentTree->GetEntry(i);
      
      // reset variables
      probAtSgn_          =  0.;
      probAtSgn_alt_      =  0.;
      probAtSgn_ttbb_     =  0.;
      probAtSgn_alt_ttbb_ =  0.;
      probAtSgn_alt_ttbj_ =  0.;
      probAtSgn_alt_ttcc_ =  0.;
      probAtSgn_alt_ttjj_ =  0.;
      bestMH_             = -99;
      bestMT_             = -99;
      bestMTbb_           = -99;
      bestMTbj_           = -99;
      bestMTcc_           = -99;
      bestMTjj_           = -99;
      matchesH_           = -99;
      matchesW_           = -99;
      matchesT_           = -99;
      matchesHAll_        = -99;
      matchesWAll_        = -99;
      matchesTAll_        = -99;
      overlapHeavy_       = -99;
      overlapLight_       = -99;
      type_               = -99;
      nPermut_            = +99;
      nTotPermut_         = +99;
      nInter_             =   1;
      flag_type0_         = -99;
      flag_type1_         = -99;
      flag_type2_         = -99;
      flag_type3_         = -99;
      flag_type4_         = -99;
      flag_type6_         = -99;
      nSimBs_             = -99;
      nMatchSimBs_        = -99;
      nC_                 = -99;
      nCTop_              = -99;
      time_               = -99;


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN PARTICLES @@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // define top decay products;      
      TLorentzVector topBLV   (1,0,0,1);
      TLorentzVector topW1LV  (1,0,0,1); 
      TLorentzVector topW2LV  (1,0,0,1);
      TLorentzVector atopBLV  (1,0,0,1); 
      TLorentzVector atopW1LV (1,0,0,1); 
      TLorentzVector atopW2LV (1,0,0,1);
      TLorentzVector genBLV   (1,0,0,1);
      TLorentzVector genBbarLV(1,0,0,1);

      // set-up top decay products (if available from input file...)
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
      if(genB.mass>0 && (genB.momid==25 || genB.momid==23)){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && (genBbar.momid==25 || genBbar.momid==23)){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
  
      // define LV for the 8 particles in TTH events
      TLorentzVector TOPLEP  (1,0,0,1);
      TLorentzVector TOPHAD  (1,0,0,1);
      TLorentzVector TOPHADW1(1,0,0,1);
      TLorentzVector TOPHADW2(1,0,0,1);
      TLorentzVector TOPHADB (1,0,0,1);
      TLorentzVector TOPLEPW1(1,0,0,1);
      TLorentzVector TOPLEPW2(1,0,0,1);
      TLorentzVector TOPLEPB (1,0,0,1);

      // dummy cut (for the moment)
      bool properEventSL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);
      bool properEventDL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      if(!(properEventSL || properEventDL)){
	cout << "A dummy cut has failed..." << endl;
	continue;
      }

      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
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
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else if(abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPHADW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	}
	else{
	  TOPHADW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPHADW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	}
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{}


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN BHADRONS @@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      // find out number of b-hadrons in the event...
      nMatchSimBs_ = 0;
      for(int l = 0; l<nSimBs; l++){
	TLorentzVector Bs(1,0,0,1);
	Bs.SetPtEtaPhiM( SimBspt[l], SimBseta[l], SimBsphi[l], SimBsmass[l]);	    
	
	if( topBLV.Pt()>10 && deltaR(topBLV,  Bs)<0.5 ) continue;
	if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs)<0.5 ) continue;
	
	for(int hj = 0; hj<nhJets; hj++){
	  TLorentzVector hJLV(1,0,0,1);
	  if(hJet_genPt[hj]>10) 
	    hJLV.SetPtEtaPhiM( hJet_genPt[hj], hJet_genEta[hj], hJet_genPhi[hj], 0.0);
	  if( hJLV.Pt()>20 && TMath::Abs(hJLV.Eta())<5 && deltaR(Bs, hJLV)<0.5 ) nMatchSimBs_++;
	}
	for(int aj = 0; aj<naJets; aj++){
	  TLorentzVector aJLV(1,0,0,1);
	  if(aJet_genPt[aj]>10) 
	    aJLV.SetPtEtaPhiM( aJet_genPt[aj], aJet_genEta[aj], aJet_genPhi[aj], 0.0);
	  if( aJLV.Pt()>20 && TMath::Abs(aJLV.Eta())<5 && deltaR(Bs, aJLV)<0.5 ) nMatchSimBs_++;
	}	
      }
      if( nSimBs>=2){
	for(int l = 0; l<nSimBs-1; l++){
	  TLorentzVector Bs1(1,0,0,1);
	  Bs1.SetPtEtaPhiM( SimBspt[l], SimBseta[l], SimBsphi[l], SimBsmass[l]);	    
	  if( topBLV.Pt()>10 && deltaR(topBLV,  Bs1)<0.5 ) continue;
	  if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs1)<0.5 ) continue;
	  for(int m = l+1; m<nSimBs; m++){
	    TLorentzVector Bs2(1,0,0,1);
		Bs2.SetPtEtaPhiM( SimBspt[m], SimBseta[m], SimBsphi[m], SimBsmass[m]);	    	    
		if( topBLV.Pt()>10 && deltaR(topBLV,  Bs2)<0.5 ) continue;
		if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs2)<0.5 ) continue;
		if( deltaR(Bs1,Bs2)<0.50 ) nMatchSimBs_--;
	  }
	}
      }

      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@ LEPTON SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */
      

      // charged leptons and MET
      TLorentzVector leptonLV, leptonLV2;
      TLorentzVector neutrinoLV;


      ///////////////////////////////////
      //         SL events             //
      ///////////////////////////////////

      properEventSL = false;      
      if( nvlep==1 ){

	// first lepton...
	leptonLV.SetPtEtaPhiM(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);      
	// cut on lepton (SL)
	properEventSL = 
	  (vLepton_type[0]==13 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.1 && vLepton_pfCorrIso[0]<0.10) ||
	  (vLepton_type[0]==11 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.5 && !(TMath::Abs(leptonLV.Eta())>1.442 &&  TMath::Abs(leptonLV.Eta())<1.566) && vLepton_pfCorrIso[0]<0.10) ;	

      }


      ///////////////////////////////////
      //         DL events             //
      ///////////////////////////////////

      properEventDL = false;
      if( nvlep>=2 ){
	
	// first lepton...
	leptonLV.SetPtEtaPhiM (vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
	// second lepton...
	leptonLV2.SetPtEtaPhiM(vLepton_pt[1],vLepton_eta[1],vLepton_phi[1],vLepton_mass[1]);

	// cut on leptons (DL)
	properEventDL = (
			 ( vLepton_type[0]==13 && vLepton_type[1]==13 && 
			   ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.1 && vLepton_pfCorrIso[0]<0.10 &&
			      leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.4 && vLepton_pfCorrIso[1]<0.20) ||
			     (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.4 && vLepton_pfCorrIso[0]<0.20 &&
			      leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.1 && vLepton_pfCorrIso[1]<0.10) )
			   ) ||
			 ( vLepton_type[0]==11 && vLepton_type[0]==11 && 
			   ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.10 &&
			      leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.20) ||
			     (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.20 &&
			      leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.10) )
			   )
			 )
	  &&  vLepton_charge[0]*vLepton_charge[1]<0 ;
      }

      ////////////////////////////////////////////////////////////////////////

      // MET
      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
      neutrinoLV.SetPxPyPzE(nuPx,nuPy,0. ,nuE);

      // continue if leptons do not satisfy cuts
      if( !(properEventSL || properEventDL) ) continue;


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // this container orders the jets by decreasing pt
      std::map<double, JetObservable , sorterByPt> jet_map;
      std::map<double, JetObservable >::iterator   jet_map_it;

      // loop over jet collections
      for(int coll = 0 ; coll < 2 ; coll++){

	// loop over jets
	for(int hj = 0; hj < (coll==0 ? nhJets : naJets); hj++){

	  float pt     = (coll==0) ? hJet_pt [hj]  : aJet_pt [hj];
	  float eta    = (coll==0) ? hJet_eta[hj]  : aJet_eta[hj];
	  float phi    = (coll==0) ? hJet_phi[hj]  : aJet_phi[hj];
	  float e      = (coll==0) ? hJet_e  [hj]  : aJet_e  [hj];
	  float m2     = e*e - pt*pt*TMath::CosH(eta)*TMath::CosH(eta);
	  if(m2<0) m2 = 0.; 
	  float m      = TMath::Sqrt( m2 ); 

	  int id       = (coll==0) ? hJet_puJetIdL[hj] : aJet_puJetIdL[hj];
	  float JECUnc = (coll==0) ? hJet_JECUnc  [hj] : aJet_JECUnc  [hj];

	  // for JEC systematics
	  float shift     = 1.0;
	  if     ( doJECup   )  shift *= (1+JECUnc);
	  else if( doJECdown )  shift *= (1-JECUnc);
	  else{}
	  pt *= shift;
	  m  *= shift;

	  //cout << "Jet #" << coll << "-" << hj << " => (" << pt << "," << eta << "," << phi << "," << m << "), ID=" <<
	  //id << endl;

	  // only jets in acceptance...
	  if( TMath::Abs(eta)> 2.5 ) continue;
	  // only jets passing ID...
	  //if( id < 0.5 ) continue;	
	  // only jets above pt cut...
	  if( pt < 30  ) continue;	  

	  TLorentzVector p4;
	  p4.SetPtEtaPhiM( pt, eta, phi, m );

	  // for csv systematics
	  float csv_nominal =  (coll==0) ? hJet_csv_nominal[hj] : aJet_csv_nominal[hj];
	  float csv_upBC    =  (coll==0) ? hJet_csv_upBC   [hj] : aJet_csv_upBC   [hj];
	  float csv_downBC  =  (coll==0) ? hJet_csv_downBC [hj] : aJet_csv_downBC [hj];
	  float csv_upL     =  (coll==0) ? hJet_csv_upL    [hj] : aJet_csv_upL    [hj];
	  float csv_downL   =  (coll==0) ? hJet_csv_downL  [hj] : aJet_csv_downL  [hj];
	  float csv = csv_nominal;
	  if     ( doCSVup  ) csv =  TMath::Max(csv_upBC,   csv_upL);
	  else if( doCSVdown) csv =  TMath::Min(csv_downBC, csv_downL);
	  else{}

	  // the jet observables (p4 and csv)
	  JetObservable myJet;
	  myJet.p4  = p4;
	  myJet.csv = csv; 
	  
	  // use pt to order jet collection
	  jet_map[ p4.Pt() ] = myJet;
	  
	}
      }

   
      // fill arrays of jets
      std::vector<TLorentzVector>  jets_p4;
      std::vector<double>          jets_csv;
      std::vector<double>          jets_csv_prob_b;
      std::vector<double>          jets_csv_prob_c;
      std::vector<double>          jets_csv_prob_j;
      int jetsAboveCut = 0;

      for( jet_map_it = jet_map.begin() ; jet_map_it != jet_map.end(); jet_map_it++){

	//cout << "Map: " << ((jet_map_it->second).p4).Pt() << "," << (jet_map_it->second).csv << endl;

	TLorentzVector p4 = (jet_map_it->second).p4;
	float csv         = (jet_map_it->second).csv;
	// (Min needed because in csvUp, csv can exceed 1...)
	csv               =  TMath::Min( csv, float(0.999999) );

	// count jets above 40 GeV
	if( p4.Pt()>40 ) jetsAboveCut++;

	// store jet p4...
	jets_p4.push_back ( p4  );
	// store csv
	jets_csv.push_back( csv );
	
	// needed to find appropriate csv PDF
	string bin = "";
	if( TMath::Abs( p4.Eta() ) <= 1.0 ) 
	  bin = "Bin0";
	if( TMath::Abs( p4.Eta() ) >  1.0 ) 
	  bin = "Bin1";

	// store PDF(csv)
	jets_csv_prob_b.push_back( btagger["b_"+bin]!=0 ? btagger["b_"+bin]->GetBinContent( btagger["b_"+bin]->FindBin( csv ) ) : 1.);
	jets_csv_prob_c.push_back( btagger["c_"+bin]!=0 ? btagger["c_"+bin]->GetBinContent( btagger["c_"+bin]->FindBin( csv ) ) : 1.);
	jets_csv_prob_j.push_back( btagger["l_"+bin]!=0 ? btagger["l_"+bin]->GetBinContent( btagger["l_"+bin]->FindBin( csv ) ) : 1.);

      }
      // continue if not enough jets
      if( jetsAboveCut<4 ){
	//cout << "Less then 4 jets.. continue" << endl;
	continue;
      }

      
      // jet multiplicity
      int numJets30UntagM = 0; 
      int numJets30UntagL = 0; 
      int numJets30BtagL  = 0; 
      int numJets30BtagM  = 0; 
      int numJets30BtagT  = 0;

      // indices of jets_p4 for the first five un-tagged jets (CSVM)
      unsigned int w1=999;
      unsigned int w2=999;
      unsigned int w3=999;
      unsigned int w4=999;
      unsigned int w5=999;
      // indices of jets_p4 for the first four tagged jets (CSVM)
      unsigned int b1=999;
      unsigned int b2=999;
      unsigned int b3=999;
      unsigned int b4=999;
      // indices of jets_p4 for the first five un-tagged jets (CSVL)
      unsigned int wL1=999;
      unsigned int wL2=999;
      unsigned int wL3=999;
      unsigned int wL4=999;
      unsigned int wL5=999;
      // indices of jets_p4 for the first four tagged jets (CSVL)
      unsigned int bL1=999;
      unsigned int bL2=999;
      unsigned int bL3=999;
      unsigned int bL4=999;

      for(unsigned int k = 0; k < jets_p4.size(); k++){  

	float csv_k = jets_csv[k];
	float pt_k  = (jets_p4[k]).Pt();

	// passes CSVL...
	int btag_L = csv_k>0.244 ;	
	// passes CSVM...
	int btag_M = csv_k>0.679 ;
	// passes CSVT...
	int btag_T = csv_k>0.898 ;	

	if( pt_k>30 &&  btag_L )  numJets30BtagL ++;
	if( pt_k>30 &&  btag_M )  numJets30BtagM ++;
	if( pt_k>30 &&  btag_T )  numJets30BtagT ++;
	if( pt_k>30 && !btag_M )  numJets30UntagM++;
	if( pt_k>30 && !btag_L )  numJets30UntagL++;

	if     (  btag_M && b1==999 && b2==999 && b3==999 && b4==999)            b1 = k;	   
	else if(  btag_M && b1!=999 && b2==999 && b3==999 && b4==999)            b2 = k;	    
	else if(  btag_M && b1!=999 && b2!=999 && b3==999 && b4==999)            b3 = k;	    
	else if(  btag_M && b1!=999 && b2!=999 && b3!=999 && b4==999)            b4 = k;	    
	else if( !btag_M && w1==999 && w2==999 && w3==999 && w4==999 && w5==999) w1 = k;
	else if( !btag_M && w1!=999 && w2==999 && w3==999 && w4==999 && w5==999) w2 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3==999 && w4==999 && w5==999) w3 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3!=999 && w4==999 && w5==999) w4 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3!=999 && w4!=999 && w5==999) w5 = k; 
	else{}

	if     (  btag_L && bL1==999 && bL2==999 && bL3==999 && bL4==999)             bL1 = k;	    
	else if(  btag_L && bL1!=999 && bL2==999 && bL3==999 && bL4==999)             bL2 = k;	    
	else if(  btag_L && bL1!=999 && bL2!=999 && bL3==999 && bL4==999)             bL3 = k;	    
	else if(  btag_L && bL1!=999 && bL2!=999 && bL3!=999 && bL4==999)             bL4 = k;	    
	else if( !btag_L && wL1==999 && wL2==999 && wL3==999 && wL4==999 && wL5==999) wL1 = k;
	else if( !btag_L && wL1!=999 && wL2==999 && wL3==999 && wL4==999 && wL5==999) wL2 = k; 
	else if( !btag_L && wL1!=999 && wL2!=999 && wL3==999 && wL4==999 && wL5==999) wL3 = k; 
	else if( !btag_L && wL1!=999 && wL2!=999 && wL3!=999 && wL4==999 && wL5==999) wL4 = k; 
	else if( !btag_L && wL1!=999 && wL2!=999 && wL3!=999 && wL4!=999 && wL5==999) wL5 = k; 
	else{}

      }


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ EVENT SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */


      ////////////////////////////////////////////////////
      // SEMILEPTONIC EVENTS WITH ==4 BTAG && >=1 UNTAG //
      // FULLLEPTONIC EVENTS WITH ==4 BTAG              //
      ////////////////////////////////////////////////////

      //cout << "Recap: numJets30BtagM=" << numJets30BtagM << ", numJets30UntagM=" << numJets30UntagM << endl;

      bool analyze_type0 = properEventSL && numJets30BtagM==4 && numJets30UntagM==2 && doType0;
      bool analyze_type1 = properEventSL && numJets30BtagM==4 && numJets30UntagM==2 && doType1;
      bool analyze_type2 = properEventSL && numJets30BtagM==4 && numJets30UntagM==1 && doType2;
      bool analyze_type3 = properEventSL && numJets30BtagM==4 && numJets30UntagM >2 && doType3;
      bool analyze_type6 = properEventDL && numJets30BtagM==4 && doType6;
      bool analyze_type7 = properEventDL && numJets30BtagM==3 && numJets30BtagL==4  && doType7;


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN MATCH @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      vector<unsigned int> btag_indices;
      if( analyze_type0 || analyze_type1 || analyze_type2 || analyze_type3 || analyze_type6){
	if(b1!=999) btag_indices.push_back( b1 );
	if(b2!=999) btag_indices.push_back( b2 );
	if(b3!=999) btag_indices.push_back( b3 );
	if(b4!=999) btag_indices.push_back( b4 );
      }
      else if(analyze_type7){
	if(bL1!=999) btag_indices.push_back( bL1 );
	if(bL2!=999) btag_indices.push_back( bL2 );
	if(bL3!=999) btag_indices.push_back( bL3 );
	if(bL4!=999) btag_indices.push_back( bL4 );
      }
      vector<unsigned int> buntag_indices;
      if( analyze_type0 || analyze_type1 || analyze_type2 || analyze_type3){
	if(w1!=999) buntag_indices.push_back( w1 );
	if(w2!=999) buntag_indices.push_back( w2 );
	if(w3!=999) buntag_indices.push_back( w3 );
	if(w4!=999) buntag_indices.push_back( w4 );
	if(w5!=999) buntag_indices.push_back( w5 );
      }
      
      // temporary variables
      int hMatches = 0;
      int tMatches = 0;
      int wMatches = 0;

      for( unsigned int w = 0; w<btag_indices.size(); w++){

	if     (  genBLV.Py()>0    && deltaR(jets_p4[ btag_indices[w] ],genBLV)   <GENJETDR ) hMatches++;
	else if(  genBbarLV.Py()>0 && deltaR(jets_p4[ btag_indices[w] ],genBbarLV)<GENJETDR)  hMatches++;    
	if     (  TOPHADB.Py()>0   && deltaR(jets_p4[ btag_indices[w] ],TOPHADB)  <GENJETDR ) tMatches++;
	else if(  TOPLEPB.Py()>0   && deltaR(jets_p4[ btag_indices[w] ],TOPLEPB)  <GENJETDR)  tMatches++;
	if     (  TOPHADW1.Py()>0  && deltaR(jets_p4[ btag_indices[w] ],TOPHADW1) <GENJETDR ) wMatches++;
	else if(  TOPHADW2.Py()>0  && deltaR(jets_p4[ btag_indices[w] ],TOPHADW2) <GENJETDR)  wMatches++;
    
      }
      matchesH_    = hMatches;
      matchesHAll_ = hMatches;
      matchesT_    = tMatches;
      matchesTAll_ = tMatches;
      matchesWAll_ = wMatches;  
      wMatches = 0;

      for( unsigned int w = 0; w<buntag_indices.size(); w++){

	if     (  genBLV.Py()>0    && deltaR(jets_p4[ buntag_indices[w] ],genBLV)   <GENJETDR ) matchesHAll_++;
	else if(  genBbarLV.Py()>0 && deltaR(jets_p4[ buntag_indices[w] ],genBbarLV)<GENJETDR ) matchesHAll_++;
	if     (  TOPHADB.Py()>0   && deltaR(jets_p4[ buntag_indices[w] ],TOPHADB)  <GENJETDR ) matchesTAll_++;
	else if(  TOPLEPB.Py()>0   && deltaR(jets_p4[ buntag_indices[w] ],TOPLEPB)  <GENJETDR)  matchesTAll_++;
	if     (  TOPHADW1.Py()>0  && deltaR(jets_p4[ buntag_indices[w] ],TOPHADW1) <GENJETDR ) wMatches++;
	else if(  TOPHADW2.Py()>0  && deltaR(jets_p4[ buntag_indices[w] ],TOPHADW2) <GENJETDR)  wMatches++;
    
      }
      matchesW_    =  wMatches; 
      matchesWAll_ += wMatches;

      // gen level overlap
      vector<TLorentzVector> genHeavy;
      if(genBLV.Py()>0)    genHeavy.push_back( genBLV);
      if(genBbarLV.Py()>0) genHeavy.push_back( genBbarLV);
      if(TOPHADB.Py()>0)   genHeavy.push_back( TOPHADB);
      if(TOPLEPB.Py()>0)   genHeavy.push_back( TOPLEPB);
      int overlapH = 0;
      if(genHeavy.size()>1){
	for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	  for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	    if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	  }	    
	}
      }
      overlapHeavy_ = overlapH;
      vector<TLorentzVector> genLight;
      if( TOPHADW1.Py()>0) genLight.push_back( TOPHADW1);
      if( TOPHADW2.Py()>0) genLight.push_back( TOPHADW2);
      int overlapL = 0;
      for(unsigned int k = 0; k < genLight.size(); k++){  
	for(unsigned int l = 0; l < genHeavy.size(); l++){  
	  if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	}	    
      }
      if( genLight.size()>1 && deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
      overlapLight_  = overlapL;      
      


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ANALYSIS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      //  input 4-vectors
      vector<TLorentzVector> jets;

      // internal map: [ position in "jets" ] -> [ position in "jets_p4" ]
      std::map< unsigned int, unsigned int> pos_to_index;


      if( analyze_type0 || analyze_type1 || analyze_type2 || analyze_type3 || analyze_type6 || analyze_type7){	
	
	vector<unsigned int> untaggedjets;
	if( w1!=999 ) untaggedjets.push_back(w1);
	if( w2!=999 ) untaggedjets.push_back(w2);
	if( w3!=999 ) untaggedjets.push_back(w3);
	if( w4!=999 ) untaggedjets.push_back(w4);
	if( w5!=999 ) untaggedjets.push_back(w5);
	
	// find out which two untagged jets come from W->qq'
	unsigned int ind1 = 999;
	unsigned int ind2 = 999;
	
	if( (analyze_type0 || analyze_type1) ){
	  
	  // sanity check:
	  if(untaggedjets.size() != 2){
	    cout << "Inconsistency found for (analyze_type0 || analyze_type1)... continue" << endl;
	    continue;
	  }

	  // use untagged mass to assign to type 0 OR type 1
	  float WMass = (jets_p4[w1]+jets_p4[w2]).M();
	  
	  // set index for untagged jets
	  ind1 = w1;
	  ind2 = w2;
	  
	  if( (WMass>MwL && WMass<MwH)  && analyze_type0 ){
	    
	    counter++;      
	    if( fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	    cout << "Processing event # " << counter << " (TYPE 0), mW=" << WMass << " GeV." << " Event ID " << EVENT.event << endl;
	    
	    /////////////////////////////////////////////////////	      
	    type_       =  0;
	    nPermut_    = 12;
	    nInter_     =  1;
	    nTotPermut_ = 12;
	    meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	    /////////////////////////////////////////////////////
	  }
	  else if( !( WMass>MwL && WMass<MwH) && analyze_type1){
	    
	    counter++;      
	    if( fixNumEvJob && !(counter>=evLow && counter<=evHigh)  ) continue;
	    cout << "Processing event # " << counter << " (TYPE 1), mW=" << WMass << " GeV." << " Event ID " << EVENT.event << endl;	   
	    
	    /////////////////////////////////////////////////////	  
	    type_       =  1;
	    nPermut_    = 12;
	    nInter_     =  2;
	    nTotPermut_ = 24;
	    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	    /////////////////////////////////////////////////////
	  }
	  else{ continue; }

	}
	else if( analyze_type2 ){	   
	  
	  // sanity check:
	  if(untaggedjets.size() != 1){
	    cout << "Inconsistency found for analyze_type2... continue" << endl;
	    continue;
	  }

	  counter++;
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (TYPE 1)." << " Event ID " << EVENT.event << endl;
	  
	  // set index for untagged jet
	  ind1 = w1;
	  ind2 = w1;
	  
	  /////////////////////////////////////////////////////
	  type_       =  2;
	  nPermut_    = 12;
	  nInter_     =  1;
	  nTotPermut_ = 12;
	  meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	  /////////////////////////////////////////////////////
	}
	else if( analyze_type3 ){
	  
	  // sanity check:
	  if(untaggedjets.size() <= 2 ){
	    cout << "Inconsistency found for analyze_type3... continue" << endl;
	    continue;
	  }

	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (TYPE 3)." << " Event ID " << EVENT.event << endl;	   
	  
	  // find out which are ind1 and ind2...
	  float minDiff     = 99999.;
	  for(unsigned int uj1 = 0; uj1<untaggedjets.size()-1; uj1++){
	    for(unsigned int uj2 = uj1+1; uj2<untaggedjets.size(); uj2++){
	      
	      float WMass12 = (jets_p4[ untaggedjets[uj1] ]+jets_p4[ untaggedjets[uj2] ]).M();
	      if( TMath::Abs(WMass12-80.1)<minDiff ){
		minDiff = TMath::Abs(WMass12-MW);
		ind1 = untaggedjets[uj1];
		ind2 = untaggedjets[uj2];
	      }
	      
	    }
	  }	  
	  /////////////////////////////////////////////////////
	  type_       =  3;
	  nPermut_    = 12;
	  nInter_     =  1;
	  nTotPermut_ = 12;
	  meIntegrator->setIntType( MEIntegratorNew::SL2wj );
	  /////////////////////////////////////////////////////	  
	}
	else if( analyze_type6 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (TYPE 6)." << " Event ID " << EVENT.event << endl;	   

	  /////////////////////////////////////////////////////
	  type_       =  6;
	  nPermut_    = 12;
	  nInter_     =  1;
	  nTotPermut_ = 12;
	  meIntegrator->setIntType( MEIntegratorNew::DL );
	  /////////////////////////////////////////////////////	  

	}
	else if( analyze_type7 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (TYPE 7)." << " Event ID " << EVENT.event << endl;	   
	  
	  /////////////////////////////////////////////////////
	  type_       =  7;
	  nPermut_    = 12;
	  nInter_     =  1;
	  nTotPermut_ = 12;
	  meIntegrator->setIntType( MEIntegratorNew::DL );
	  /////////////////////////////////////////////////////	  

	}
	else{ 
	  cout << "Inconsistency in the analysis... continue." << endl; 
	  continue; 
	}


	// sanity-check
	if(type_<=3 && (ind1==999 || ind2==999)){
	  cout << "Inconsistency found: ind1 or ind2 are not set...continue." << endl;
	  continue;
	}
	
	
	// setup jet collection
	jets.clear();

	// b1,...,w1,w2 are indices for jets_p4 collection;
	// This is a map between the internal ordering bLep=2, W1Had=3, ..., higgs2 = 7, and jets_p4
	pos_to_index.clear();

	if( type_<=3 ){
	  jets.push_back( leptonLV     );  
	  jets.push_back( neutrinoLV   );  
	  jets.push_back( jets_p4[b1]  );
	  jets.push_back( jets_p4[ind1]);  
	  jets.push_back( jets_p4[ind2]);  
	  jets.push_back( jets_p4[b2]  );
	  jets.push_back( jets_p4[b3]  );
	  jets.push_back( jets_p4[b4]  );
	  
	  pos_to_index[2] = b1;
	  pos_to_index[3] = ind1;
	  pos_to_index[4] = ind2;
	  pos_to_index[5] = b2;
	  pos_to_index[6] = b3;
	  pos_to_index[7] = b4;
	}
	else if( type_==6 ){
	  jets.push_back( leptonLV    );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( jets_p4[b1] );
	  jets.push_back( leptonLV2   );  
	  jets.push_back( neutrinoLV  ); // dummy  
	  jets.push_back( jets_p4[b2] );
	  jets.push_back( jets_p4[b3] );
	  jets.push_back( jets_p4[b4] );
	}
	else if( type_==7 ){
	  jets.push_back( leptonLV    );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( jets_p4[bL1] );
	  jets.push_back( leptonLV2   );  
	  jets.push_back( neutrinoLV  ); // dummy  
	  jets.push_back( jets_p4[bL2] );
	  jets.push_back( jets_p4[bL3] );
	  jets.push_back( jets_p4[bL4] );
	}
	
	// init particle momenta...
	meIntegrator->setJets(&jets);
	
	
	// set all prob. to 0.0;
	for(int v = 0 ; v < nTotPermut_; v++){
	  probAtSgn_permut_    [v] = 0.;
	  probAtSgn_alt_permut_[v] = 0.;
	  probAtSgn_bb_permut_ [v] = 0.;
	  probAtSgn_bj_permut_ [v] = 0.;
	  probAtSgn_cc_permut_ [v] = 0.;
	  probAtSgn_jj_permut_ [v] = 0.;
	}
	
	/////////////////////////////////////////////////////////////
	
	// check if there is a tag-untag pair that satisfies the "cs-tag" 
	vector<unsigned int> btag_indices;
	if( type_<6 ){
	  btag_indices.push_back( b1 );
	  btag_indices.push_back( b2 );
	  btag_indices.push_back( b3 );
	  btag_indices.push_back( b4 );
	}
	for( unsigned int w = 0; w<btag_indices.size(); w++){
	  
	  float m1 = ( jets_p4[btag_indices[w]] + jets_p4[ind1] ).M();
	  float m2 = ( jets_p4[btag_indices[w]] + jets_p4[ind2] ).M();
	  
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	    flag_type0_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type0_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type0_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type0_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type0_ = 4;  
	  }
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	    flag_type1_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type1_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type1_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type1_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type1_ = 4;  
	  }
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 2 ){
	    flag_type2_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type2_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type2_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type2_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type2_ = 4;  
	  }
	}
	// for type 3, the W-tag is different...
	float WMass = type_==3 ? (jets_p4[ ind1 ]+jets_p4[ ind2 ]).M() : -999.;
	if( WMass>MwL && WMass<MwH )  flag_type3_ = 1;
	
	
	/////////////////////////////////////////////////////////////
	
	// init MET stuff
	meIntegrator->setSumEt( METtype1p2corr.sumet );
	meIntegrator->setMEtCov(-99,-99,0);
	
	// specify if topLep has pdgid +6 or -6
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	
	// type1 has two iterations: 
	//  > iter0 := treat w1 as W jet and w2 as a gluon
	//  > iter1 := treat w2 as W jet and w1 as a gluon	  
	
	for(int iter=0; iter<nInter_; iter++){
	  	  
	  // if type1, swap the light jets
	  if( type_ == 1 && iter==1){
	    jets.clear();
	    jets.push_back( leptonLV      );  
	    jets.push_back( neutrinoLV    );  
	    jets.push_back( jets_p4[b1]   );
	    jets.push_back( jets_p4[ind2] ); // <-- here
	    jets.push_back( jets_p4[ind1] ); // <-- here 
	    jets.push_back( jets_p4[b2]   );
	    jets.push_back( jets_p4[b3]   );
	    jets.push_back( jets_p4[b4]   );
	    
	    // update the map
	    pos_to_index.clear();
	    pos_to_index[2] = b1;
	    pos_to_index[3] = ind2;  // <-- here
	    pos_to_index[4] = ind1;  // <-- here
	    pos_to_index[5] = b2;
	    pos_to_index[6] = b3;
	    pos_to_index[7] = b4;
	    
	    // reset the jets
	    meIntegrator->setJets(&jets);	
	  }
	  
	  // for each event, clean the histograms
	  h_prob_vs_higgs_mass   ->Reset();
	  h_prob_vs_top_mass     ->Reset();
	  h_prob_vs_top_mass_ttbb->Reset();
	  h_prob_vs_top_mass_ttbj->Reset();
	  h_prob_vs_top_mass_ttcc->Reset();
	  h_prob_vs_top_mass_ttjj->Reset();
	    
	  // start the clock...	  
	  clock->Start();
	  
	  // loop over Higgs mass values...
	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	    
	    // loop over Top mass values...
	    for(int t = 0; t < nTopMassPoints ; t++){
	      meIntegrator->setTopMass( mT[t] , MW );
	      
	      // these are used for bookkeeping
	      double maxP_s = 0.;
	      double maxP_b = 0.;
		
	      // loop over permutations
	      for(unsigned int pos = 0; pos < (unsigned int)nPermut_ ; pos++){
		
		// consider permutation #pos
		if( type_==0 || type_==3)
		  meIntegrator->initVersors( permutations_SL2wj[pos] );
		if( type_==1 || type_==2)
		  meIntegrator->initVersors( permutations_SL1wj[pos] );
		if( type_>=6 )
		  meIntegrator->initVersors( permutations_DL[pos] );

		// check invariant mass of jet system:
		double mass, massLow, massHigh;
		bool skip        = !( meIntegrator->compatibilityCheck    (0.95, /*printP4*/ 0, mass, massLow, massHigh ) );
		bool skip_WHad   = false;
		bool skip_TopHad = false;
		if((type_==0 || type_==3)){
		  skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		  skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		}

		// if type 0/3 and incompatible with MW or MT (and we are not scanning vs MT) continue
		// ( this applies to both hypothesis )
		if( nTopMassPoints==1 && (skip_WHad || skip_TopHad) ){
		  cout << "Skip                                 Perm. #" << pos << endl;
		  continue;
		}
		
		// retrieve intergration boundaries from meIntegrator
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);	    
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 =  useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
		
		double x0L = range_x0.first; double x0U = range_x0.second;
		double x1L = range_x1.first; double x1U = range_x1.second;
		double x2L = range_x2.first; double x2U = range_x2.second;
		double x3L = range_x3.first; double x3U = range_x3.second;
		double x4L = range_x4.first; double x4U = range_x4.second;
		double x5L = range_x5.first; double x5U = range_x5.second;
		double x6L = range_x6.first; double x6U = range_x6.second;
		
		// these hold for the sgn integration and type0...
		double xLmode0_s[4] = {x0L, x3L, x4L, x5L};
		double xUmode0_s[4] = {x0U, x3U, x4U, x5U};
		// these hold for the bkg integration and type0...
		double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};		 
		
		// these hold for the sgn integration and type1...
		double xLmode1_s[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode1_s[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
		// these hold for the bkg integration and type1...
		double xLmode1_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode1_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
		
		// these hold for the sgn integration and type2...
		double xLmode2_s[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode2_s[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
		// these hold for the bkg integration and type2...	      
		double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};	     	     	   	       
		
		// these hold for the sgn integration and type3...
		double xLmode3_s[4] = {x0L, x3L, x4L, x5L};
		double xUmode3_s[4] = {x0U, x3U, x4U, x5U};
		// these hold for the bkg integration and type3...
		double xLmode3_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode3_b[5] = {x0U, x3U, x4U, x5U, x6U};
		
		// these hold for the sgn integration and type6...
		double xLmode6_s[5] = {x1L, x2L, x1L, x2L, x5L};
		double xUmode6_s[5] = {x1U, x2U, x1U, x2U, x5U};
		// these hold for the bkg integration and type6...
		double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
		double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};

		// these hold for the sgn integration and type7...
		double xLmode7_s[5] = {x1L, x2L, x1L, x2L, x5L};
		double xUmode7_s[5] = {x1U, x2U, x1U, x2U, x5U};
		// these hold for the bkg integration and type7...
		double xLmode7_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
		double xUmode7_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};

		// number of integration variables (TTH hypothesis)
		int nParam;
		if     ( type_==0 )  nParam = 4;
		else if( type_==1 )  nParam = 6;
		else if( type_==2 )  nParam = 6;
		else if( type_==3 )  nParam = 4;
		else if( type_==6 )  nParam = 5;
		else if( type_==7 )  nParam = 5;
		else{
		  cout << "No type match...contine." << endl;
		  continue;
		}
		
		// the per-permutation probability...
		double p = 0.;	       	
		  
		// loop over hypothesis [ TTH, TTbb ]
		for(int hyp = 0 ; hyp<2;  hyp++){
		  
		  // if doing higgs mass scan, don't consider bkg hypo
		  if( nMassPoints>1 && hyp==1 ) continue;

		  // if doing top mass scan, don't consider sgn hypo
		  if( nTopMassPoints>1 && hyp==0 ) continue;

		  // if consider only one hypothesis (SoB=0) 
		  // and the current hypo is not the desired one, continue...
		  if( SoB==0 && hyp!=hypo) continue;

		  // if current hypo is TTH, but M(b1b2) incompatible with 125
		  // (and we are not scanning vs MH) continue...
		  if( hyp==0 && nMassPoints==1 && skip){
		    cout << "Skip    hypo " << (hyp==0 ? "ttH " : "ttbb") 
			 << ", MH=" << mH[m] << ", MT=" << mT[t]
			 << ". Perm. #" << pos;
		    if(type_==1) 
		      cout << " (inter. " << iter << ")";
		    cout << " => p=" << p << endl;
		    continue;
		  }
		  
		  // setup hypothesis
		  cout << "Testing hypo " << (hyp==0 ? "ttH " : "ttbb") 
		       << ", MH=" << mH[m] << ", MT=" << mT[t]
		       << ". Perm. #" << pos;
		  if(type_==1) 
		    cout << " (inter. " << iter << ")";
		    
		  meIntegrator->setHypo(hyp);
		  
		  // initial number of function calles
		  int intPoints = 4000;
		  if( type_==0 )  intPoints =  2000;
		  if( type_==1 )  intPoints =  4000;
		  if( type_==2 )  intPoints =  4000;
		  if( type_==3 )  intPoints =  2000;
		  if( type_==6 )  intPoints = 10000;
		  if( type_==7 )  intPoints = 10000;
		  
		  // count how many time the integration is rerun per permutation
		  int ntries = 0;
		  
		  // skip ME calculation... for debugging
		  if(speedup==0){
		    
		    // refinement: redo integration if it returned a bad chi2
		    while( ntries < MAX_REEVAL_TRIES){
		      
		      // integrand
		      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, nParam+hyp);
		      // VEGAS integrator
		      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		      ig2.SetFunction(toIntegrate);
		      // setup # of parameters
		      meIntegrator->SetPar( nParam+hyp );	 
		      
		      // the integration ranges depend on hyp and type
		      if     ( type_==0 ){
			p = (hyp==0 ? ig2.Integral(xLmode0_s, xUmode0_s) : ig2.Integral(xLmode0_b, xUmode0_b));
		      }
		      else if( type_==1 ){
			p = (hyp==0 ? ig2.Integral(xLmode1_s, xUmode1_s) : ig2.Integral(xLmode1_b, xUmode1_b));
		      }
		      else if( type_==2 ){
			p = (hyp==0 ? ig2.Integral(xLmode2_s, xUmode2_s) : ig2.Integral(xLmode2_b, xUmode2_b));
		      }
		      else if( type_==3 ){
			p = (hyp==0 ? ig2.Integral(xLmode3_s, xUmode3_s) : ig2.Integral(xLmode3_b, xUmode3_b));
		      }
		      else if( type_==6 ){
			p = (hyp==0 ? ig2.Integral(xLmode6_s, xUmode6_s) : ig2.Integral(xLmode6_b, xUmode6_b));
		      }
		      else if( type_==7 ){
			p = (hyp==0 ? ig2.Integral(xLmode7_s, xUmode7_s) : ig2.Integral(xLmode7_b, xUmode7_b));
		      }
		      else{ /* ... */ }		    
		      
		      // chi2 of the integration
		      double chi2 =  ig2.ChiSqr();
		      
		      // check if the actual permutation returned a small or large number...
		      // if the value is less than 10% of the largest found that far, skip 
		      // the improvement
		      if( hyp==0 ){
			if( p>maxP_s ) maxP_s = p;		  
			else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		      }
		      else{
			if( p>maxP_b ) maxP_b = p;		  
			else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		      }
		      
		      // if the chi2 is bad, increse # of points and repeat the integration
		      if( chi2 > maxChi2_ ){
			ntries++;
			intPoints *= 1.5;
		      }
		      // otherwise, just go to the next permutation...
		      else 
			ntries = MAX_REEVAL_TRIES+1;			
		    }	
		  }
		  else{
		    // can still be interested in b-tagging, so set p=1...
		    p = 1.;
		  }
		  
		  // to avoid problems
		  if( TMath::IsNaN(p) ) p = 0.;
		  cout << " => p=" << p << endl;
		  
		  /////////////////////////////////////////////////////////////
		  
		  if( useBtag ){
		    
		    int bLep_pos      = (type_==0 || type_==3 || type_>=6) ? (permutations_SL2wj[pos])/100000   : (permutations_SL1wj[pos])/100000;
		    int bHad_pos      = (type_==0 || type_==3 || type_>=6) ? (permutations_SL2wj[pos])%1000/100 : (permutations_SL1wj[pos])%1000/100;
		    int b1_pos        = (type_==0 || type_==3 || type_>=6) ? (permutations_SL2wj[pos])%100/10   : (permutations_SL1wj[pos])%100/10;
		    int b2_pos        = (type_==0 || type_==3 || type_>=6) ? (permutations_SL2wj[pos])%10       : (permutations_SL1wj[pos])%10;
		    
		    double p_b_bLep =  jets_csv_prob_b[ pos_to_index[bLep_pos] ];
		    double p_b_bHad =  jets_csv_prob_b[ pos_to_index[bHad_pos] ];
		    double p_b_b1   =  jets_csv_prob_b[ pos_to_index[b1_pos] ];
		    double p_c_b1   =  jets_csv_prob_c[ pos_to_index[b1_pos] ];
		    double p_j_b1   =  jets_csv_prob_j[ pos_to_index[b1_pos] ];
		    double p_b_b2   =  jets_csv_prob_b[ pos_to_index[b2_pos] ];
		    double p_c_b2   =  jets_csv_prob_c[ pos_to_index[b2_pos] ];
		    double p_j_b2   =  jets_csv_prob_j[ pos_to_index[b2_pos] ];
		      
		    // total and per-permutation ME*btag probability for nominal MH and MT
		    if( mH[m]<MH+0.5 && mH[m]>MH-0.5 && mT[t]<MT+0.5 && mT[t]>MT-0.5){
		      if(hyp==0){
			probAtSgn_ttbb_         += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 );
		      }
		      else{
			probAtSgn_alt_ttbb_     += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 );
			probAtSgn_alt_ttbj_     += ( p * p_b_bLep * p_b_bHad * (p_b_b1 * p_j_b2 + p_j_b1 * p_b_b2 )*0.5 );
			probAtSgn_alt_ttcc_     += ( p * p_b_bLep * p_b_bHad * p_c_b1 * p_c_b2 );
			probAtSgn_alt_ttjj_     += ( p * p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 );
		      }
		      probAtSgn_bb_permut_[(unsigned int)(pos+iter*12)] =  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2;
		      probAtSgn_bj_permut_[(unsigned int)(pos+iter*12)] =  p_b_bLep * p_b_bHad * (p_b_b1 * p_j_b2 + p_j_b1 * p_b_b2 )*0.5;
		      probAtSgn_cc_permut_[(unsigned int)(pos+iter*12)] =  p_b_bLep * p_b_bHad * p_c_b1 * p_c_b2;
		      probAtSgn_jj_permut_[(unsigned int)(pos+iter*12)] =  p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2;
		    }
		      
		    // total ME probability vs MT*btag for nominal MH under B hypo
		    if( hyp>0 && mH[m]<MH+0.5 && mH[m]>MH-0.5){
		      h_prob_vs_top_mass_ttbb->Fill( mT[t], ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 ) );
		      h_prob_vs_top_mass_ttbj->Fill( mT[t], ( p * p_b_bLep * p_b_bHad * (p_b_b1 * p_j_b2 + p_j_b1 * p_b_b2 )*0.5 ));
		      h_prob_vs_top_mass_ttcc->Fill( mT[t], ( p * p_b_bLep * p_b_bHad * p_c_b1 * p_c_b2 ) );
		      h_prob_vs_top_mass_ttjj->Fill( mT[t], ( p * p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 ) );
		    }			
		    
		  }

		  /////////////////////////////////////////////////////////////
		    
		  
		  // total and per-permutation ME probability for nominal MH and MT
		  if( mH[m]<MH+0.5 && mH[m]>MH-0.5&& mT[t]<MT+0.5 &&  mT[t]>MT-0.5){
		    if(hyp==0){
		      probAtSgn_     += p;
			probAtSgn_permut_[(unsigned int)(pos+iter*nPermut_)]      = p;
		    }
		    else{
		      probAtSgn_alt_ += p;
		      probAtSgn_alt_permut_[(unsigned int)(pos+iter*nPermut_)]  = p;
		    }
		  }
		  
		  // total ME probability vs MH for nominal MT under S hypo
		  if( hyp==0 && mT[t]<MT+0.5 && mT[t]>MT-0.5){
		    h_prob_vs_higgs_mass->Fill( mH[m], p );
		  }
		  // total ME probability vs MT for nominal MH under B hypo
		  if( hyp>0 && mH[m]<MH+0.5 && mH[m]>MH-0.5){
		    h_prob_vs_top_mass->Fill( mT[t],   p );
		  }
		  
		  /////////////////////////////////////////////////////////////
		  
		}  // hypothesis		  
	      }  // permutations
	    }  // nTopMAssPoints
	  }  // nMassPoints	    
	}  // iterations
	
	// get best mass value by quadratic interpolation of the mass histograms
	bestMH_   = (getMaxValue(h_prob_vs_higgs_mass))   .first;
	bestMT_   = (getMaxValue(h_prob_vs_top_mass))     .first;
	bestMTbb_ = (getMaxValue(h_prob_vs_top_mass_ttbb)).first;
	bestMTbj_ = (getMaxValue(h_prob_vs_top_mass_ttbj)).first;
	bestMTcc_ = (getMaxValue(h_prob_vs_top_mass_ttcc)).first;
	bestMTjj_ = (getMaxValue(h_prob_vs_top_mass_ttjj)).first;
	
	// stop clock and reset
	clock->Stop();
	time_ = clock->RealTime();
	clock->Reset();	    
	
      }

      ///////////////////////////////////////////////////
      // ALL THE REST...                               //
      ///////////////////////////////////////////////////

      else{
	continue;
      }
      
      counter_   = counter;
      nSimBs_    = nSimBs;
      nC_        = nC;
      nCTop_     = nCTop;
      weight_    = scaleFactor;

      EVENT_.run   = EVENT.run;
      EVENT_.lumi  = EVENT.lumi;
      EVENT_.event = EVENT.event;
      EVENT_.json  = EVENT.json;
      
      tree->Fill();      

    } // nentries

    hcounter->SetBinContent(1,float(events_)/nentries);    

  } // samples



  fout_tmp->cd();

  hcounter->Write               ("",TObject::kOverwrite );
  h_prob_vs_higgs_mass->Write   ("",TObject::kOverwrite );
  h_prob_vs_top_mass->Write     ("",TObject::kOverwrite );
  h_prob_vs_top_mass_ttbb->Write("",TObject::kOverwrite );
  h_prob_vs_top_mass_ttbj->Write("",TObject::kOverwrite );
  h_prob_vs_top_mass_ttcc->Write("",TObject::kOverwrite );
  h_prob_vs_top_mass_ttjj->Write("",TObject::kOverwrite );

  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock;
  cout << "Finished!!!" << endl;
  
  return 0;
}
