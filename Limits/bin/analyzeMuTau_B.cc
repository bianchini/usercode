#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"

#include "HiggsAnalysis/CombinedLimit/interface/TH1Keys.h"

#define VERBOSE          true
#define SAVE             true
#define addVH            true
#define EMBEDDEDSAMPLES  true
#define W3JETS           false
#define LOOSEISO         true
#define kFactorSM         1.0


///////////////////////////////////////////////////////////////////////////////////////////////

void makeHistoFromDensity(TH1* hDensity, TH1* hHistogram){

  if(hDensity->GetNbinsX() != hHistogram->GetNbinsX()){
    cout << "makeHistoFromDensity: different binning" << endl;
    return;
  }

  for(int k = 1 ; k <= hDensity->GetNbinsX(); k++){
    float bink   = hDensity->GetBinContent(k);
    float widthk = hHistogram->GetBinWidth(k);
    hDensity->SetBinContent(k, bink*widthk );
  }
  hDensity->Scale(hHistogram->Integral()/hDensity->Integral());
}

///////////////////////////////////////////////////////////////////////////////////////////////

void plotMuTau( Int_t mH_           = 120,
		Int_t useEmbedding_ = 0,
		string selection_   = "inclusive",
		string analysis_    = "",		  
		TString variable_   = "diTauVisMass",
		TString XTitle_     = "full mass",
		TString Unities_    = "GeV",
		TString outputDir   = "./",
		Int_t nBins_ = 80, Float_t xMin_=0, Float_t xMax_=400,
		Float_t magnifySgn_ = 1.0,
		Float_t hltEff_     = 1.0,
		Int_t logy_         = 0,
		Float_t maxY_       = 1.2
		) 
{   

  string postfix_ = "";

  ofstream out(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/yieldsMuTau_mH%d_%s_%s.txt",outputDir.Data(),mH_,selection_.c_str(), analysis_.c_str() ),ios_base::out); 
  out.precision(5);
  out<< " => " << selection_ << endl;

  
  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/bins/bins_muTau_%s_%s.txt",variable_.Data(), selection_.c_str())); 
  if(nBins_<0 &&  !is.good()){
    cout << "Bins file not found" << endl;
    return;
  }

  int nBinsFromFile = 0;
  while (is.good())     
    {
      is.getline(c,999,',');     
      if (is.good()){
	nBinsFromFile++;
	//cout << c << endl;
      }
    }

  // choose the number of bins
  int nBins =  nBins_>0 ? nBins_ : nBinsFromFile-1 ;
  TArrayF bins(nBins+1);
  cout << "Making histograms with " << nBins << " bins:" << endl;

  is.close();
  is.open(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/bins/bins_muTau_%s_%s.txt",variable_.Data(), selection_.c_str())); 
  
  nBinsFromFile = 0;

  if(nBins_>0){
    for( ; nBinsFromFile <= nBins ; nBinsFromFile++){
      bins[nBinsFromFile] =  xMin_ + nBinsFromFile*(xMax_-xMin_)/nBins_;
    }
  }
  else{
    while (is.good())  
      {
	is.getline(c,999,',');     
	if (is.good() && nBinsFromFile<=nBins) {
	  bins[nBinsFromFile] = atof(c);
	  cout << bins[nBinsFromFile] << ", " ;
	}
	nBinsFromFile++;
      }
    cout << endl;
  }

  // Luminosity analyzed & parameters
  //float Lumi = (-47.4 + 215.3 + 930.7 + 410.6 + (450.6+212.7) + (735.2+254.8+778.2+682.0) )*1.00;
  //float Lumi = (-47.4 + 221.7 + 970 + 390 + 706 + 2741)*1.00;
  // from lumiPixel
  float Lumi   = (-47.4 + 215.6 + 955.3 + 389.9 + 706.719 + 2714)*(1-0.056);

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  float OStoSSRatioQCD            = 1.11;
  float SSIsoToSSAIsoRatioQCD     = 1.00;

  float MutoTauCorrectionFactor   = 1.00;
  float JtoTauCorrectionFactor    = 1.00;

  float embeddedMEtCutEff         = 1.00;
  float madgraphMEtCutEff         = 1.00;

  // Fall11_06Dec2011
  float WcorrectionFactorOS        = 0.92;  
  float WcorrectionFactorSS        = 1.08; 
  float ExtrapolationFactorZ       = 1.0;
  float ExtrapolationFactorZDataMC = 1.0;
  float ErrorExtrapolationFactorZ  = 1.0;

  //float NoVbfExtrapolationFactorZ = 0.997;
  //float VbfExtrapolationFactorZ   = 1.37;
  //float BoostExtrapolationFactorZ = 0.98;

  float VbfExtrapolationFactorW   = 1.00;
  float BoostExtrapolationFactorW = 1.00;


  /////////////////  change SVfit mass here ///////////////////

  //string variableStr = "";
  //TString variable(variableStr.c_str());
  TString variable = variable_;

  //////////////////////////////////////////////////////////////

  bool useMt      = true;
  string antiWcut = useMt ? "MtLeg1MVA" : "-(pZetaMVA-1.5*pZetaVisMVA)" ; //<------------------
  //string antiWcut = useMt ? "MtLeg1" : "-(pZetaMVA-1.5*pZetaVisMVA)" ;
  float antiWsgn  = useMt ? 40. :  20. ; 
  float antiWsdb  = useMt ? 60. :  40. ; 

  bool use2Dcut   = false;
  if( use2Dcut ){
    antiWcut = "!(MtLeg1MVA<40 && (pZetaMVA-1.5*pZetaVisMVA)>-20)";
    antiWsgn = 0.5;
    antiWsdb = 0.5;
  }

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);

  TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);
 
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(logy_);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TLegend* leg = new TLegend(0.63,0.48,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=7 TeV}{%.1f fb^{-1} #tau_{#mu}#tau_{had}}", Lumi/1000. ));

  THStack* aStack = new THStack("aStack","");

  TH1F* hSiml    = new TH1F( "hSiml"   ,"all"               , nBins , bins.GetArray());
  TH1F* hSgn     = new TH1F( "hSgn "   ,"vbf+ggf"           , nBins , bins.GetArray());
  TH1F* hSgn1    = new TH1F( "hSgn1"   ,"vbf"               , nBins , bins.GetArray());
  TH1F* hSgn2    = new TH1F( "hSgn2"   ,"ggf"               , nBins , bins.GetArray());
  TH1F* hSgn3    = new TH1F( "hSgn3"   ,"vh"                , nBins , bins.GetArray());
  TH1F* hData    = new TH1F( "hData"   ,"        "          , nBins , bins.GetArray());
  TH1F* hDataEmb = new TH1F( "hDataEmb","Embedded"          , nBins , bins.GetArray());
  TH1F* hW       = new TH1F( "hW"      ,"W+jets"            , nBins , bins.GetArray());
  TH1F* hW3Jets  = new TH1F( "hW3Jets" ,"W+3jets"           , nBins , bins.GetArray());
  TH1F* hEWK     = new TH1F( "hEWK"    ,"EWK"               , nBins , bins.GetArray());
  TH1F* hZtt     = new TH1F( "hZtt"    ,"Ztautau"           , nBins , bins.GetArray());
  TH1F* hZmm     = new TH1F( "hZmm"    ,"Z+jets, mu to tau" , nBins , bins.GetArray());
  TH1F* hZmj     = new TH1F( "hZmj"    ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hZfakes  = new TH1F( "hZfakes" ,"Z+jets, jet to tau", nBins , bins.GetArray());
  TH1F* hTTb     = new TH1F( "hTTb"    ,"ttbar"             , nBins , bins.GetArray());
  TH1F* hQCD     = new TH1F( "hQCD"    ,"QCD"               , nBins , bins.GetArray());
  TH1F* hLooseIso= new TH1F( "hLooseIso","Loose Iso"        , nBins , bins.GetArray());
  TH1F* hAntiIso = new TH1F( "hAntiIso","Anti Iso"          , nBins , bins.GetArray());
  TH1F* hVV      = new TH1F( "hVV"     ,"Diboson"           , nBins , bins.GetArray());

  TH1*  hW3JetsKeys   = 0;
  TH1*  hWKeys        = 0;
  TH1*  hLooseIsoKeys = 0;
  TH1*  hAntiIsoKeys  = 0;
  TH1*  hZmmKeys      = 0;
  TH1*  hZmjKeys      = 0;
  TH1*  hVVKeys       = 0;

 
  // pZeta OS, N pZ sideband OS, pZeta SS, N sideband SS, N QCD SS, OS/SS
  TH1F* hParameters   = new TH1F( "hParameters", "" ,12, 0, 12);

  // Open the files
  TFile *fData              
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleRun2011-MuTau-All_run_Open_MuTauStream.root", "READ");  
  TFile *fDataLooseIso  ///////////////////            
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleRun2011-MuTau-All_run_Open_MuTauStream.root", "READ");  
  TFile *fDataEmbedded              
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleRun2011-MuTau-Embedded-All_run_Open_MuTauStream.root", "READ");  
  TFile *fSignalVBF         
    = new TFile(Form("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleVBFH%d-MuTau-powheg-PUS6_run_Open_MuTauStream.root",mH_) ,"READ");  
  TFile *fSignalGGH         
    = new TFile(Form("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleGGFH%d-MuTau-powheg-PUS6_run_Open_MuTauStream.root",mH_),"READ"); 
  TFile *fSignalVH         
    = new TFile(Form("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleVH%d-MuTau-pythia-PUS6_run_Open_MuTauStream.root",mH_),"READ");  
  TFile *fBackgroundDY
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root","READ"); 
  TFile *fBackgroundWJets   
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleWJets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root","READ"); 
  TFile *fBackgroundW3Jets   
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleW3Jets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root","READ"); 
  TFile *fBackgroundTTbar  
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleTTJets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root","READ"); 
  TFile *fBackgroundOthers  
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval//nTupleOthers-MuTau-PUS6_run_Open_MuTauStream.root","READ"); 

  // choose the analysis: Nominal "", jet up/Down "JetUp/Down" , elec up/down "MuUp/Down" , tau up/down "TauUp/Down"
  TString tree         = "outTreePtOrd"+postfix_+analysis_;
  TString treeEmbedded = "outTreePtOrd"+postfix_;
  if(analysis_.find("TauUp")  !=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("TauDown")!=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("MuUp")  !=string::npos) 
    treeEmbedded = tree;
  if(analysis_.find("MuDown")!=string::npos) 
    treeEmbedded = tree;

  TTree *data                = (TTree*)fData->Get(("outTreePtOrd"+postfix_).c_str());
  TTree *dataLooseIso        = LOOSEISO  ? (TTree*)fDataLooseIso->Get(("outTreePtOrd"+postfix_).c_str()) : 0;
  TTree *dataEmbedded        = EMBEDDEDSAMPLES ? (TTree*)fDataEmbedded->Get(treeEmbedded) : 0;
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);
  TTree *signalVH            = addVH ? (TTree*)fSignalVH->Get(tree) : 0;

  // split the DY->ll into l=e/mu and l=tau (MC level) ===> temporary, need fix !!!
  TFile* dummy1 = new TFile("dummy2.root","RECREATE");
  cout << "Now copying g/Z -> tau+ tau- " << endl;
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)==(23*15)");                 // g/Z -> tau+ tau-
  cout << "Now copying g/Z -> mu+mu- mu->tau" << endl;
  TTree *backgroundDYMutoTau = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) &&  leptFakeTau"); // g/Z -> mu+mu- mu->tau
  cout << "Now copying g/Z -> mu+mu- jet->tau" << endl;
  TTree *backgroundDYJtoTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) && !leptFakeTau"); // g/Z -> mu+mu- jet->tau


  cout << backgroundDYTauTau->GetEntries()  << " come from DY->tautau"         << endl;
  cout << backgroundDYMutoTau->GetEntries() << " come from DY->mumu, mu->tau"  << endl;
  cout << backgroundDYJtoTau->GetEntries()  << " come from DY->mumu, jet->tau" << endl;

  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundW3Jets    = W3JETS ? (TTree*)fBackgroundW3Jets->Get(tree) : 0;
  TTree *backgroundOthers    = (TTree*)fBackgroundOthers->Get(tree);
 

  ///// LEPT PT ///////
  TCut lpt("ptL1>17 && isPFMuon && isTightMuon");
  TCut tpt("ptL2>20");

  if(selection_.find("High")!=string::npos)
    tpt = tpt&&TCut("ptL2>40");
  else if(selection_.find("Low")!=string::npos)
    tpt = tpt&&TCut("ptL2<40");

  ////// TAU ISO //////
  TCut tiso("tightestHPSMVAWP>=0");
  //TCut tiso("tightestHPSDBWP>0 && tightestHPSMVAWP>=0"); // <-----------------------

  ////// MU ISO ///////
  TCut liso("combRelIsoLeg1DBetav2<0.10");
  TCut laiso("combRelIsoLeg1DBetav2>0.30 && combRelIsoLeg1DBetav2<0.50");
  TCut lliso("combRelIsoLeg1DBetav2<0.30");

 
  ////// EVENT WISE //////
  TCut lveto("muFlag==0");
  TCut SS("diTauCharge!=0");
  TCut OS("diTauCharge==0");
  TCut pZ( Form("((%s)<%f)",antiWcut.c_str(),antiWsgn));
  TCut apZ(Form("((%s)>%f)",antiWcut.c_str(),antiWsdb));
  TCut hltevent("pairIndex<1 && HLTx==1 && ( run>=163269 || run==1)");
  TCut hltmatch("HLTmatch==1");


  ////// CATEGORIES ///
  TCut oneJet("nJets30>=1");
  TCut twoJets("nJets30>=2");

  //TCut vbf("pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1");     // <--- BASELINE
  //TCut vbf("pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1 && jet1PUWP>0.5 && jet2PUWP>0.5 && (jetVetoPUWP>0.5 && jetVetoPUWP<0)");
  TCut vbf("pt1>30 && pt2>30 && isVetoInJets!=1 && MVAvbf>0.80");

  TCut vh("pt1>30 && pt2>30 && Mjj>70 && Mjj<120 && diJetPt>150 && MVAvbf<0.80 && nJets20BTagged<1");

  TCut boost("pt1>30 && nJets20BTagged<1"); // <--- NEW
  boost = boost && !vbf && !vh;

  TCut boost2("pt1>100 && pt1<150 && !(pt2>30 && eta1*eta2<0 && Mjj>400 && Deta>4.0 && isVetoInJets!=1)");  

  TCut bTag("nJets30<2 && nJets20BTagged>0");
  TCut nobTag("nJets30<2 && nJets20BTagged==0");

  TCut novbf = !vbf && !vh && !boost && !bTag;
  //TCut novbf("pt1<9999");   


  TCut sbin; TCut sbinEmbedding; TCut sbinEmbeddingPZetaRel; TCut sbinPZetaRel; TCut sbinSS; TCut sbinPZetaRelSS; TCut sbinPZetaRev; TCut sbinPZetaRevSS; TCut sbinSSaIso; TCut sbinSSlIso;

  TCut sbinInclusive;
  sbinInclusive            = lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch;
  TCut sbinEmbeddingInclusive;
  sbinEmbeddingInclusive   = lpt && tpt && tiso && liso && lveto && OS && pZ                         ;
  TCut sbinPZetaRelSSInclusive;
  sbinPZetaRelSSInclusive  = lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch;
  TCut sbinSSInclusive;
  sbinSSInclusive          = lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch;
  TCut sbinSSaIsoInclusive;
  sbinSSaIsoInclusive      = lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch;

  if(selection_.find("inclusive")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                         ;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                               ;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch;
  }
  else if(selection_.find("oneJet")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && oneJet;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && oneJet;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && oneJet;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && oneJet;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && oneJet;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && oneJet;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && oneJet;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && oneJet;
  }
  else if(selection_.find("twoJets")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && twoJets;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && twoJets;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && twoJets;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && twoJets;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && twoJets;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && twoJets;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && twoJets;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && twoJets;
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && vbf;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && vbf;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && vbf;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && vbf;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && vbf;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && vbf;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && vbf;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && vbf;
    
  }
  else if(selection_.find("vh")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && vh;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && vh;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && vh;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && vh;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && vh;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && vh;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && vh;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && vh;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && vh;

  }
  else if(selection_.find("novbf")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && novbf;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && novbf;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && novbf;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && novbf;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && novbf;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && novbf;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && novbf;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && novbf;
  }
  else if(selection_.find("boost")!=string::npos && selection_.find("boost2")==string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && boost;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && boost;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && boost;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && boost;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && boost;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && boost;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && boost;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost;
  }
  else if(selection_.find("boost2")!=string::npos ){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && boost2;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && boost2;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && boost2;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && boost2;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && boost2;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && boost2;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && boost2;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && boost2;
  }
  else if(selection_.find("bTag")!=string::npos && selection_.find("nobTag")==string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && bTag;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && bTag;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && bTag;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && bTag;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && bTag;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && bTag;    
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && bTag;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && bTag;
  }
  else if(selection_.find("nobTag")!=string::npos){
    sbin                   =  lpt && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch && nobTag;
    sbinEmbedding          =  lpt && tpt && tiso && liso && lveto && OS && pZ                          && nobTag;
    sbinEmbeddingPZetaRel  =  lpt && tpt && tiso && liso && lveto && OS                                && nobTag;
    sbinPZetaRel           =  lpt && tpt && tiso && liso && lveto && OS        && hltevent && hltmatch && nobTag;
    sbinPZetaRev           =  lpt && tpt && tiso && liso && lveto && OS && apZ && hltevent && hltmatch && nobTag;
    sbinPZetaRevSS         =  lpt && tpt && tiso && liso && lveto && SS && apZ && hltevent && hltmatch && nobTag;
    sbinSS                 =  lpt && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinPZetaRelSS         =  lpt && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch && nobTag;
    sbinSSaIso             =  lpt && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
    sbinSSlIso             =  lpt && tpt && tiso && lliso&& lveto && SS && pZ  && hltevent && hltmatch && nobTag;
  }


  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  cout << "******** Extrapolation factors for Z->tautau normalization ********" << endl;
  // inclusive DY->tautau:
  TH1F* hExtrap = new TH1F("hExtrap","",nBins , bins.GetArray());
  backgroundDYTauTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinInclusive);
  float ExtrapDYInclusive = hExtrap->Integral()*Lumi*hltEff_/1000.;
  hExtrap->Reset();
  cout << "All Z->tautau = " << ExtrapDYInclusive << endl; 

  TCut sbinInclusiveEmbeddedCut = sbinEmbeddingInclusive;
  TCut sbinEmbeddedCut          = sbinEmbedding;

  // if VBF, minimize ttbar contamination asking for 0 btag jet:
  if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    sbinInclusiveEmbeddedCut = sbinInclusiveEmbeddedCut && TCut("nJets20BTagged<1");
    sbinEmbeddedCut          = sbinEmbeddedCut          && TCut("nJets20BTagged<1");
  }

  dataEmbedded->Draw(variable+">>hExtrap", "(HLTTau*HLTMu*embeddingWeight)"*sbinInclusiveEmbeddedCut);
  float ExtrapEmbedDen =  hExtrap->Integral();
  hExtrap->Reset();
  dataEmbedded->Draw(variable+">>hExtrap", "(HLTTau*HLTMu*embeddingWeight)"*sbinEmbeddedCut);
  float ExtrapEmbedNum =  hExtrap->Integral();
  hExtrap->Reset();

  ExtrapolationFactorZ = ExtrapEmbedNum/ExtrapEmbedDen;

  ErrorExtrapolationFactorZ = TMath::Sqrt(ExtrapolationFactorZ*(1-ExtrapolationFactorZ)/ExtrapEmbedDen);
  cout << "Extrap. factor using embedded sample: " << ExtrapolationFactorZ << " +/- " << ErrorExtrapolationFactorZ << endl;
  backgroundDYTauTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbin);
  float ExtrapolationFactorMadGraph = hExtrap->Integral()*Lumi*hltEff_/1000./ExtrapDYInclusive;
  cout << "MadGraph prediction = " << ExtrapolationFactorMadGraph << endl;
  ExtrapolationFactorZDataMC  = ExtrapolationFactorZ/ExtrapolationFactorMadGraph;
  cout << " ==> data/MC = " << ExtrapolationFactorZDataMC << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  cout << "******** Extrapolation factors for QCD normalization ********" << endl;
  hExtrap->Reset();
  backgroundWJets->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&pZ));
  float ExtrapSSWinSignalRegionMC   = hExtrap->Integral();
  hExtrap->Reset();
  backgroundWJets->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapSSWinSidebandRegionMC = hExtrap->Integral();
  float ExtrapscaleFactorSS         = ExtrapSSWinSignalRegionMC>0 ? ExtrapSSWinSidebandRegionMC/ExtrapSSWinSignalRegionMC : 1.0;
  cout << " Extrapolation factor W SS (inclusive) " << ExtrapscaleFactorSS << endl;

  hExtrap->Reset();
  backgroundTTbar->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapttbarExtrSS    = hExtrap->Integral()*Lumi/1000*hltEff_;
  hExtrap->Reset();
  backgroundOthers->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapothersExtrSS   = hExtrap->Integral()*Lumi/1000*hltEff_;
  hExtrap->Reset();
  backgroundDYJtoTau->Draw(variable+">>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSSInclusive&&apZ));
  float ExtrapdyjtotauExtrSS = hExtrap->Integral()*Lumi/1000*hltEff_;

  hExtrap->Reset();
  data->Draw(variable+">>hExtrap", sbinPZetaRelSSInclusive&&apZ);
  float ExtrapSSWinSignalRegionDATA = hExtrap->Integral();
  cout << "Extrapolation for QCD (inclusive): total data events in sideband " << ExtrapSSWinSignalRegionDATA << endl;
  ExtrapSSWinSignalRegionDATA -= ExtrapttbarExtrSS;
  ExtrapSSWinSignalRegionDATA -= ExtrapothersExtrSS;
  ExtrapSSWinSignalRegionDATA -= ExtrapdyjtotauExtrSS;
  ExtrapSSWinSignalRegionDATA /= ExtrapscaleFactorSS;
  cout << "Extrapolation for QCD (inclusive): W+jets in SS signal region (inclusive) is estimated to be " << ExtrapSSWinSignalRegionDATA << endl;

  hExtrap->Reset();
  data->Draw(variable+">>hExtrap", sbinSSInclusive);
  float SSeventsExtrap = hExtrap->Integral();
  cout << "Extrapolation for SS events in data (inclusive) " << hExtrap->GetEntries() << endl;
  cout << "Subtracting W+jets (SS)..." << endl;
  SSeventsExtrap  -= ExtrapSSWinSignalRegionDATA;

  hExtrap->Reset();
  backgroundTTbar->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000*hltEff_;
  hExtrap->Reset();
  backgroundDYMutoTau->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000*hltEff_*MutoTauCorrectionFactor;

  hExtrap->Reset();
  backgroundDYJtoTau->Draw(variable+">>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSSInclusive);
  SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000*hltEff_*JtoTauCorrectionFactor;
  hExtrap->Reset();

  SSeventsExtrap *= OStoSSRatioQCD;

  dataLooseIso->Draw(variable+">>hExtrap", sbinSSaIsoInclusive);
  float SSeventsExtrapAiso = hExtrap->GetEntries();
  SSIsoToSSAIsoRatioQCD = SSeventsExtrap/SSeventsExtrapAiso ;
  cout << "The extrapolation factor Iso>0.3 / Iso<0.1 is " << SSIsoToSSAIsoRatioQCD << endl;

  cout << "************** END extrapolation *******************" << endl;
  delete hExtrap;
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // estimate the W+jets in the selection bin using pZeta extrapolation

  //TH1F* hWMt = new TH1F("hWMt","",1,-10,10);
  TH1F* hWMt = new TH1F("hWMt","",nBins , bins.GetArray());

  ///////////////////////////////////////// Doing OS first...
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the W+jets OS sample" << endl;
  float OSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_/1000.;
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float OSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_/1000.;
  float scaleFactorOS = OSWinSignalRegionMC>0 ? OSWinSidebandRegionMC/OSWinSignalRegionMC : 1.0 ;

  if(useMt)
    cout << "Extrapolation factor for W OS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorOS << endl;
  else
    cout << "Extrapolation factor for W OS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorOS << endl;    
 
  hWMt->Reset();
  cout << "Estimating cobtribution from Ztt, ttbar and others in OS low pZeta tail..." << endl;
  backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float ttbarExtrOS  = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from ttbar in OS is " << ttbarExtrOS << endl;
  hWMt->Reset();
  backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float othersExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from single-t and di-boson in OS is " << othersExtrOS << endl;
  hWMt->Reset();
  backgroundDYTauTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dytautauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->tautau in OS is " << dytautauExtrOS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dyjtotauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->mumu, jet->tau in OS is " << dyjtotauExtrOS << endl;
  hWMt->Reset();
  backgroundDYMutoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRel&&apZ));
  float dymutotauExtrOS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->mumu, mu->tau in OS is " << dymutotauExtrOS << endl;
  hWMt->Reset();

  data->Draw(variable+">>hWMt", sbinPZetaRev);
  float OSWinSignalRegionDATA = hWMt->Integral();
  cout << "Selected events in data in low pZeta/low Mt tail " << OSWinSignalRegionDATA << endl;
  OSWinSignalRegionDATA -= ttbarExtrOS;
  OSWinSignalRegionDATA -= othersExtrOS;
  OSWinSignalRegionDATA -= dytautauExtrOS;
  OSWinSignalRegionDATA -= dyjtotauExtrOS;
  OSWinSignalRegionDATA -= dymutotauExtrOS;
  OSWinSignalRegionDATA /= scaleFactorOS;
  cout << "W+jets in signal region is estimated to be "  
       << OSWinSignalRegionDATA*scaleFactorOS << "/" << scaleFactorOS << " = " 
       << OSWinSignalRegionDATA <<  " +/- " << sqrt(OSWinSignalRegionDATA/scaleFactorOS)/scaleFactorOS << endl;
  cout << "  ===> the MC prediction was " << OSWinSignalRegionMC << endl;

  hParameters->SetBinContent(1, 1./scaleFactorOS );
  hParameters->SetBinContent(2, OSWinSignalRegionDATA*scaleFactorOS );

  ///////////////////////////////////////// Doing SS last...
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&pZ));
  cout << "Using  " << hWMt->GetEntries() << " entries from the SS W+jets sample" << endl;
  float SSWinSignalRegionMC   = hWMt->Integral()*Lumi*hltEff_/1000.;
  hWMt->Reset();
  backgroundWJets->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float SSWinSidebandRegionMC = hWMt->Integral()*Lumi*hltEff_/1000.;
  float scaleFactorSS = SSWinSignalRegionMC>0 ? SSWinSidebandRegionMC/SSWinSignalRegionMC : 1.0;
 
  if(useMt)
    cout << "Extrapolation factor for W SS : P(MtCorr>" << antiWsdb << ")/P(MtCorr<" << antiWsgn << ") ==> " << scaleFactorSS << endl;
  else
    cout << "Extrapolation factor for W SS : P(pZetaCorr<- "<< antiWsdb << ")/P(pZetaCorr>"<< antiWsgn << ") ==> " << scaleFactorSS << endl;    

  hWMt->Reset();
  cout << "Estimating cobtribution Ztt,from ttbar and others in SS low pZeta tail..." << endl;
  backgroundTTbar->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float ttbarExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from ttbar in SS is " << ttbarExtrSS << endl;
  hWMt->Reset();
  backgroundOthers->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float othersExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from single-t and di-boson in SS is " << othersExtrSS << endl;
  hWMt->Reset();
  backgroundDYJtoTau->Draw(variable+">>hWMt","(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelSS&&apZ));
  float dyjtotauExtrSS = hWMt->Integral()*Lumi*hltEff_/1000.;
  cout << "Contribution from DY->mumu, jet->tau in SS is " << dyjtotauExtrSS << endl;
  hWMt->Reset();

  data->Draw(variable+">>hWMt",sbinPZetaRevSS);
  float SSWinSignalRegionDATA = hWMt->Integral();
  cout << "Selected events in data in low pZeta/low Mt tail " << SSWinSignalRegionDATA << endl;
  SSWinSignalRegionDATA -= ttbarExtrSS;
  SSWinSignalRegionDATA -= othersExtrSS;
  SSWinSignalRegionDATA -= dyjtotauExtrSS;
  SSWinSignalRegionDATA /= scaleFactorSS;
  cout << "W+jets in SS signal region is estimated to be "  
       << SSWinSignalRegionDATA*scaleFactorSS << "/" << scaleFactorSS << " = " 
       << SSWinSignalRegionDATA <<  " +/- " << sqrt(SSWinSignalRegionDATA/scaleFactorSS)/scaleFactorSS << endl;
  cout << "  ===> the MC prediction was " << SSWinSignalRegionMC << endl;

  hParameters->SetBinContent(3, 1./scaleFactorSS );
  hParameters->SetBinContent(4, SSWinSignalRegionDATA*scaleFactorSS );

  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("Data");
  if(dataLooseIso){
    samples.push_back("LooseIso");
    samples.push_back("AntiIso");
  }
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  if(signalVH)
    samples.push_back("VH115");
  samples.push_back("Others");
  samples.push_back("TTbar");
  samples.push_back("SS");
  samples.push_back("WJets");
  if(backgroundW3Jets)
    samples.push_back("W3Jets");
  samples.push_back("DYMutoTau");
  samples.push_back("DYJtoTau");
  samples.push_back("DYToTauTau");
  if(dataEmbedded)
    samples.push_back("Embedded");
  

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]         = data;
  tMap["LooseIso"]     = dataLooseIso;
  tMap["AntiIso"]      = dataLooseIso;
  tMap["Embedded"]     = dataEmbedded;
  tMap["ggH115"]       = signalGGH;
  tMap["qqH115"]       = signalVBF;
  tMap["VH115"]        = signalVH;
  tMap["DYToTauTau"]   = backgroundDYTauTau;
  tMap["DYMutoTau"]    = backgroundDYMutoTau;
  tMap["DYJtoTau"]     = backgroundDYJtoTau;
  tMap["WJets"]        = backgroundWJets;
  tMap["W3Jets"]       = backgroundW3Jets;
  tMap["Others"]       = backgroundOthers;
  tMap["TTbar"]        = backgroundTTbar;
  tMap["SS"]           = data;


  
  std::map<TString,Float_t> vMap;


  for( unsigned iter=0; iter<samples.size(); iter++){

    cout << "Dealing with sample " << samples[iter] << endl;
    
    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,"" , nBins , bins.GetArray());

    TTree* currentTree = 0;
    
    if((it->first).find("SS")!=string::npos){
      
      cout << "Remove W contamination from SS data sample ... " << endl;
      currentTree = (it->second);

      float error2OnQCD = 0.0;
      
      TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
      hHelp->Reset();
      currentTree->Draw(variable+">>hHelp", sbinSS);
      int SSevents = hHelp->GetEntries();
      cout << "Selected SS events in data " << hHelp->GetEntries() << endl;
      h1->Add(hHelp,1);

      hHelp->Reset();
      backgroundWJets->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from W+jets (from " << hHelp->GetEntries() << " entries)" << endl;
      float sFWSS = ( selection_.find("novbf") !=string::npos || 
		      selection_.find("bTag")  !=string::npos || 
		      selection_.find("boost") !=string::npos || 
		      selection_.find("inclusive")!=string::npos) ? 
	SSWinSignalRegionDATA/SSWinSignalRegionMC : WcorrectionFactorSS; // from the extrapolation factor DATA/MC

      if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos) 
	sFWSS *= VbfExtrapolationFactorW;
      else if(selection_.find("boost")!=string::npos)
	sFWSS *= BoostExtrapolationFactorW;

      hHelp->Scale(sFWSS*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from W+jets by extrapolating" << endl;
      cout << " ==> removing W+jets from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow( hHelp->Integral()/hHelp->GetEntries(), 2)*hHelp->GetEntries(); // error on MC W+jets SS events
      error2OnQCD +=  pow(WcorrectionFactorSS*0.06,2)*pow(hHelp->GetEntries(),2);        // error on W+jets extrapolation factor ==> 6% according to Artur
      cout << sqrt(error2OnQCD) << " <==  W" << endl;      

      hHelp->Reset();
      backgroundTTbar->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from TTbar (from " << hHelp->GetEntries() << " entries)" << endl;
      hHelp->Scale(1.0*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from TTbar" << endl;
      cout << " ==> removing TTbar from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC TTbar SS events
      cout << sqrt(error2OnQCD) << " <== W + TTb" << endl;      

      hHelp->Reset();
      backgroundDYMutoTau->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from DY->mumu, mu->jet" << endl;
      hHelp->Scale(MutoTauCorrectionFactor*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->mumu, mu->tau" << endl;
      cout << " ==> removing DY->mumu, mu->tau from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC DY->mumu, mu->tau events
      cout << sqrt(error2OnQCD) << " <== W + TTb + DY(1)" << endl;      

      hHelp->Reset();
      backgroundDYJtoTau->Draw(variable+">>hHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFTau*SFMu)"*sbinSS);
      cout << "We expect " << hHelp->Integral()*Lumi/1000*hltEff_ << " SS events from DY->mumu, jet->tau" << endl;
      hHelp->Scale(JtoTauCorrectionFactor*Lumi/1000*hltEff_);
      cout << "We estimate " << hHelp->Integral() << " SS events from DY->mumu, jet->tau" << endl;
      cout << " ==> removing DY->mumu, mu->jet from SS...." << endl;
      h1->Add(hHelp, -1 );
      if(hHelp->GetEntries()>0) error2OnQCD +=  pow(hHelp->Integral()/hHelp->GetEntries(),2)*hHelp->GetEntries();   // error on MC DY->mumu, jet->tau events
      cout << sqrt(error2OnQCD) << " <== W + TTb + DY(1,2)" << endl;      

      //  OS/SS ratio
      h1->Scale(OStoSSRatioQCD);
      cout << "After removing the expected contribution from W+jets and rescaling by " << OStoSSRatioQCD << " we expect " 
	   << h1->Integral() << " events from QCD processes" << endl;

      hParameters->SetBinContent(5, SSevents);
      hParameters->SetBinContent(6, h1->Integral()/SSevents);

      cout << "Total unceratinty from bkg subtraction in SS region is " << sqrt(error2OnQCD) << endl;
      float totalRelErrorOnQCD = 0.02 + sqrt(error2OnQCD)/h1->Integral(); //0.02 ==> uncertainty on OS/SS ratio
      hParameters->SetBinContent(7,totalRelErrorOnQCD);

      ////////////////////////////////////////////////

      hParameters->SetBinContent(8,ExtrapolationFactorZ);
      hParameters->SetBinContent(9,ErrorExtrapolationFactorZ);
      hParameters->SetBinContent(10,ExtrapolationFactorZDataMC);
      hParameters->SetBinContent(11,SSIsoToSSAIsoRatioQCD);
      hParameters->SetBinContent(12,embeddedMEtCutEff/madgraphMEtCutEff);


    }
    else{

      currentTree = (it->second);

      if((it->first).find("Embed")==string::npos){

	if((it->first).find("DYToTauTau")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbinPZetaRel);  
	  float madgraphNoMEtCut = h1->Integral();
	  h1->Reset();
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau)"*sbin);
	  madgraphMEtCutEff = h1->Integral()/madgraphNoMEtCut;
	  cout << "Efficiency of antiW cut on madgraph " << madgraphMEtCutEff << endl;
	}
	else if((it->first).find("W3Jets")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hW3JetsKeys = new TH1Keys("hW3JetsKeys","W+3jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hW3JetsKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for W3Jets filled with integral " << hW3JetsKeys->Integral() << " and entries " << hW3JetsKeys->GetEntries() << endl;
	}
	else if((it->first).find("WJets")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hWKeys = new TH1Keys("hWKeys","W+jets smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hWKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for WJets filled with integral " << hWKeys->Integral() << " and entries " << hWKeys->GetEntries() << endl;
	}
	else if((it->first).find("DYMutoTau")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hZmmKeys = new TH1Keys("hZmmKeys","Z+jets, mu to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hZmmKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmm filled with integral " << hZmmKeys->Integral() << " and entries " << hZmmKeys->GetEntries() << endl;
	}
	else if((it->first).find("DYJtoTau")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hZmjKeys = new TH1Keys("hZmjKeys","Z+jets, jet to tau smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hZmjKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for Zmj filled with integral " << hZmjKeys->Integral() << " and entries " << hZmjKeys->GetEntries() << endl;
	}
	else if((it->first).find("Others")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,     "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  hVVKeys = new TH1Keys("hVVKeys","Others smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hVVKeys", "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
	  cout << "Keys for VV filled with integral " << hVVKeys->Integral() << " and entries " << hVVKeys->GetEntries() << endl;
	}
	else if((it->first).find("LooseIso")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,    sbinSSlIso);
	  hLooseIsoKeys = new TH1Keys("hLooseIsoKeys","Loose Iso smoothed", nBins , bins.GetArray());
	  if(  ((selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos) || 
		selection_.find("boost")!=string::npos ||
		selection_.find("vh")!=string::npos) )
	    currentTree->Draw(variable+">>hLooseIsoKeys", sbinSSlIso);
	  cout << "Keys for LooseIso filled with integral " << hLooseIsoKeys->Integral() << " and entries " << hLooseIsoKeys->GetEntries() << endl;
	}
	else if((it->first).find("AntiIso")!=string::npos){
	  currentTree->Draw(variable+">>"+h1Name,    sbinSSaIso);
	  hAntiIsoKeys = new TH1Keys("hAntiIsoKeys","Anti Iso smoothed", nBins , bins.GetArray());
	  currentTree->Draw(variable+">>hAntiIsoKeys", sbinSSaIso);
	  cout << "Keys for AntiIso filled with integral " << hAntiIsoKeys->Integral() << " and entries " << hAntiIsoKeys->GetEntries() << endl;
	}
	else
	  currentTree->Draw(variable+">>"+h1Name, "(sampleWeight*puWeight*HLTweightTau*HLTweightMu*SFMu*SFTau*HqTWeight)"*sbin);
      }
      else{
	currentTree->Draw(variable+">>"+h1Name, "(HLTTau*HLTMu*embeddingWeight)"*sbinEmbeddingPZetaRel);
	float embeddedNoMEtCut = h1->Integral();
	h1->Reset();
	currentTree->Draw(variable+">>"+h1Name, "(HLTTau*HLTMu*embeddingWeight)"*sbinEmbedding);

	embeddedMEtCutEff =  h1->Integral()/embeddedNoMEtCut;
	cout << "Efficiency of antiW cut on embedded " << embeddedMEtCutEff << endl;
      }


      // scale by correction factors
      if( ! ((it->first).find("Data")!=string::npos || 
	     (it->first).find("LooseIso")!=string::npos ||
	     (it->first).find("AntiIso")!=string::npos) ) 
	h1->Scale(Lumi/1000*hltEff_);

      // if W+jets, scale by extrapolation
      float sFWOS = ( selection_.find("novbf")!=string::npos  || 
		      selection_.find("boost")!=string::npos  || 
		      selection_.find("bTag")!=string::npos   || 
		      selection_.find("inclusive")!=string::npos) ? 
	OSWinSignalRegionDATA/OSWinSignalRegionMC : WcorrectionFactorOS;
      if((it->first).find("WJets")!=string::npos){

	if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	  sFWOS *= VbfExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << VbfExtrapolationFactorW << " according to the Z->mumu+j+vbf/Z->mumu+j ratio" << endl;
	}
	else if(selection_.find("boost")!=string::npos){
	  sFWOS *= BoostExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << BoostExtrapolationFactorW << " according to the Z->mumu+j+vbf/Z->mumu+j ratio" << endl;
	}	else if(selection_.find("boost")!=string::npos){
	  sFWOS *= BoostExtrapolationFactorW;
	  cout << "Wjets will be rescaled by " << BoostExtrapolationFactorW << " according to the Z->mumu+j+vbf/Z->mumu+j ratio" << endl;
	}

	h1->Scale( sFWOS );
	hW->Add(h1,1.0);
      }

      // if DY->tautau, and vbf scale by ratio data/MC
      if((it->first).find("DYToTauTau")!=string::npos){


	h1->Scale( ExtrapolationFactorZDataMC );

	//if(selection_.find("novbf")!=string::npos){
	//cout << "DY->tautau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//h1->Scale( NoVbfExtrapolationFactorZ );
	//}
	//else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	//cout << "DY->tautau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//h1->Scale( VbfExtrapolationFactorZ );
	//}
	//else if(selection_.find("boost")!=string::npos){
	//cout << "DY->tautau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//h1->Scale( BoostExtrapolationFactorZ );
	//}
      }

      // if DY->mumu, mu->tau, scale by fake-rate
      if((it->first).find("DYMutoTau")!=string::npos){

	float sF = MutoTauCorrectionFactor;

	sF *= ExtrapolationFactorZDataMC;

	//if(selection_.find("novbf")!=string::npos){
	//sF *= NoVbfExtrapolationFactorZ;
	//cout << "DY->tautau, mu->tau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//}
	//else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	//sF *= VbfExtrapolationFactorZ;
	//cout << "DY->tautau, mu->tau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//}
	//else if(selection_.find("boost")!=string::npos){
	//cout << "DY->tautau, mu->tau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//sF *= BoostExtrapolationFactorZ;
	//}

	h1->Scale(sF);
	hZmm->Add(h1,1.0);
	hZfakes->Add(h1,1.0);
      }

      // if DY->mumu, jet->tau, scale by fake-rate
      if((it->first).find("DYJtoTau")!=string::npos){

	float sF = JtoTauCorrectionFactor;

	sF *= ExtrapolationFactorZDataMC;

	//if(selection_.find("novbf")!=string::npos){
	//sF *= NoVbfExtrapolationFactorZ;
	//cout << "DY->tautau, jet->tau will be rescaled by " << NoVbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//}
	//else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
	//sF *= VbfExtrapolationFactorZ;
	//cout << "DY->tautau, jet->tau will be rescaled by " << VbfExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//}
	//else if(selection_.find("boost")!=string::npos){
	//cout << "DY->tautau, jet->tau will be rescaled by " << BoostExtrapolationFactorZ << " according to the Z->mumu+vbf/Z->mumu ratio" << endl;
	//sF *=BoostExtrapolationFactorZ;
	//}

	h1->Scale(sF);
	hZmj->Add(h1,1.0);
	hZfakes->Add(h1,1.0);
      }

      if((it->first).find("ggH115")!=string::npos && selection_.find("boost")!=string::npos){
	h1->Scale(kFactorSM);
      }


    }
  
    /////////////////////////////////////////////////////////////////////////////////////

    if( (it->first).find("DYMutoTau")!=string::npos ||  
	(it->first).find("DYJtoTau")!=string::npos || 
	(it->first).find("WJets")!=string::npos || 
	(it->first).find("Others")!=string::npos ) {
      hEWK->SetFillColor(kRed-3);
      hEWK->Add(h1,1.0);

      if( (it->first).find("Others")!=string::npos ){
	hVV->Add(h1,1.0);
	if(hVVKeys->Integral()>0) hVVKeys->Scale(hVV->Integral()/hVVKeys->Integral());
	hVVKeys->SetFillColor(kRed-3);
      }
    }

    if( (it->first).find("DYToTauTau")!=string::npos ) {
      hZtt->Add(h1,1.0);
      hZtt->SetFillColor(kYellow-9);
    }
    if( (it->first).find("Embedded")!=string::npos ) {
      //if(hZtt->Integral()>0) h1->Scale(hZtt->Integral()/h1->Integral());
      h1->Scale( (ExtrapolationFactorZ*ExtrapDYInclusive)/h1->Integral());
      h1->Scale(embeddedMEtCutEff/madgraphMEtCutEff);

      hDataEmb->Add(h1,1.0);
      hDataEmb->SetFillColor(kYellow-9);
    }
    if( (it->first).find("DYMutoTau")!=string::npos ) {
      if(hZmmKeys->Integral()>0) hZmmKeys->Scale(hZmm->Integral()/hZmmKeys->Integral());
      hZmmKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("DYJtoTau")!=string::npos ) {
      if(hZmjKeys->Integral()>0) hZmjKeys->Scale(hZmj->Integral()/hZmjKeys->Integral());
      hZmjKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("WJets")!=string::npos ) {
      if(hWKeys->Integral()>0) hWKeys->Scale(hW->Integral()/hWKeys->Integral());
      hWKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("W3Jets")!=string::npos ) {
      hW3Jets->Add(h1,1.0);
      if(hW3Jets->Integral()>0) hW3Jets->Scale(hW->Integral()/hW3Jets->Integral());
      hW3Jets->SetFillColor(kRed-3);
      if(hW3JetsKeys->Integral()>0) hW3JetsKeys->Scale(hW->Integral()/hW3JetsKeys->Integral());
      hW3JetsKeys->SetFillColor(kRed-3);
    }
    if( (it->first).find("LooseIso")!=string::npos ) {
      hLooseIso->Add(h1,1.0);
      if(hLooseIsoKeys->Integral()>0) hLooseIsoKeys->Scale(hLooseIso->Integral()/hLooseIsoKeys->Integral());
      hLooseIsoKeys->SetFillColor(0);
    }
    if( (it->first).find("AntiIso")!=string::npos ) {
      float sF = SSIsoToSSAIsoRatioQCD*OStoSSRatioQCD;
      cout << "Anti-isolated SS data scaled by " <<  SSIsoToSSAIsoRatioQCD << " ratio measured in inclusive sample ==>" << endl;
      cout << "   SS anti-isolated events = " << h1->GetEntries() << " ==> " << h1->Integral()*SSIsoToSSAIsoRatioQCD << " predicted in signal region" << endl;
      h1->Scale(sF);
      hAntiIso->Add(h1,1.0);
      if(hAntiIsoKeys->Integral()>0) hAntiIsoKeys->Scale(hAntiIso->Integral()/hAntiIsoKeys->Integral());
      hAntiIsoKeys->SetFillColor(0);
    }
    if( (it->first).find("TTbar")!=string::npos ) {
      hTTb->Add(h1,1.0);
      hTTb->SetFillColor(kBlue-2);     
    }
    if( (it->first).find("SS")!=string::npos ) {
      hQCD->Add(h1,1.0);
      hQCD->SetFillColor(kMagenta-9);
    }
    if((it->first).find("qqH115")!=string::npos){
      hSgn1->Add(h1,1.0);
      hSgn1->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("ggH115")!=string::npos){
      hSgn2->Add(h1,1.0);
      hSgn2->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("VH115")!=string::npos){
      hSgn3->Add(h1,1.0);
      hSgn3->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn3->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
    }
    if((it->first).find("Data")!=string::npos){
      hData->Add(h1,1.0);
      hData->Sumw2();
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.2);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetXTitle(XTitle_+" ("+Unities_+")");
      hData->SetYTitle(Form(" Events/(%.1f %s)", hData->GetBinWidth(1), Unities_.Data() ) );
      hData->SetTitleSize(0.04,"X");
      hData->SetTitleSize(0.05,"Y");
      hData->SetTitleOffset(0.95,"Y");
    }

    if(VERBOSE) cout<<(it->first) << " ==> " 
		    << h1->Integral() << " +/- " 
		    << TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries())
		    << endl;

    //out << (it->first) << "  " << h1->Integral() << " $\\pm$ " <<  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()) << endl;
    char* c = new char[50];
    if(h1->Integral()>=10) 
      sprintf(c,"$%.0f\\pm%.0f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=1)
      sprintf(c,"$%.1f\\pm%.1f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.1)
      sprintf(c,"$%.2f\\pm%.2f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else if(h1->Integral()>=0.01)
      sprintf(c,"$%.3f\\pm%.3f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    else
      sprintf(c,"$%.5f\\pm%.5f$",h1->Integral(),  TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries()));
    out << string(c) << "  //" << (it->first) << endl;
    delete c;

  }

  out.close();

  // all signal summed together:
  hSgn->Add(hSgn1,hSgn2,1,1);
  hSgn->Add(hSgn3);
  hSgn->SetFillColor(0);
  hSgn->SetLineColor(kBlue);
  hSgn->SetLineWidth(2);
  hSgn->SetLineStyle(kDashed);

  // adding alltogether
  hSiml->Add(hQCD,1.0);
  hSiml->Add(hEWK,1.0);
  hSiml->Add(hTTb,1.0);
  if(useEmbedding_)
    hSiml->Add(hDataEmb,1.0);
  else
    hSiml->Add(hZtt,1.0);


  //Adding to the stack
  aStack->Add(hQCD);
  aStack->Add(hEWK);
  aStack->Add(hTTb);
  if(useEmbedding_)
    aStack->Add(hDataEmb);
  else
    aStack->Add(hZtt);
  if(!logy_)
    aStack->Add(hSgn);

  // legend
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hSgn,Form("(%.0fx) H#rightarrow#tau#tau m_{H}=%d",magnifySgn_,mH_),"F");
  if(useEmbedding_)
    leg->AddEntry(hDataEmb,"Z#rightarrow#tau#tau (embedded)","F");
  else
    leg->AddEntry(hZtt,"Z#rightarrow#tau#tau","F"); 
  leg->AddEntry(hTTb,"t#bar{t}","F");
  leg->AddEntry(hEWK,"Electroweak","F");
  leg->AddEntry(hQCD,"QCD","F");
  
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  if(!logy_)
    hData->SetAxisRange(0.0, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  else
    hData->SetAxisRange(0.1, TMath::Max( hData->GetMaximum(), hSiml->GetMaximum() )*maxY_ ,"Y");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  if(logy_)
    hSgn->Draw("HISTSAME");
  

  leg->Draw();

  pad2->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  //TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA-MC)}{#sqrt{DATA}}",
  //		  hStack->GetNbinsX(), 
  //		  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());

  TH1F* hRatio = new TH1F( "hRatio" ," ; ; #frac{(DATA-MC)}{#sqrt{DATA}}" , nBins , bins.GetArray());
  hRatio->Reset();
  hRatio->SetXTitle("");
  hRatio->SetYTitle("#frac{(DATA-MC)}{#sqrt{DATA}}");

  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerSize(0.8);
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");
  TH1F* hError = (TH1F*)hRatio->Clone("hError");
  hError->Reset();
  hError->SetFillColor(kRed);
  hError->SetFillStyle(3004);
  hError->SetMarkerStyle(kDot);

  float uncertZtt = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertZtt += (0.06  * 0.06) ; // Tau-ID 
    uncertZtt += (0.035 * 0.035); // Lumi 
    uncertZtt += (0.10  * 0.10);  // Extrap. factor 
  }
  uncertZtt = TMath::Sqrt(uncertZtt);
  float uncertTTb = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertTTb += (0.075 * 0.075) ; // xsection 
  }
  uncertTTb = TMath::Sqrt(uncertTTb);
  float uncertEWK = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertEWK += (0.07 * 0.07) ; // extrapolation 
    uncertEWK += (0.10 * 0.10) ; // extrapolation  
  }
  uncertEWK = TMath::Sqrt(uncertEWK);
  float uncertQCD = 0;
  if(selection_.find("novbf")!=string::npos || selection_.find("inclusive")!=string::npos){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
  }
  else if(selection_.find("vbf")!=string::npos && selection_.find("novbf")==string::npos){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
    uncertQCD += (0.10 * 0.10) ; // extrapolation 
  }
  else if(selection_.find("boost")!=string::npos){
    uncertQCD += (0.02 * 0.02) ; // extrapolation 
    uncertQCD += (0.10 * 0.10) ; // extrapolation  
  }
  uncertQCD = TMath::Sqrt(uncertQCD);

  float maxPull = 0.;
  for(int k = 0 ; k < hRatio->GetNbinsX(); k++){
    float pull = hData->GetBinContent(k) - hSiml->GetBinContent(k);
    if(hData->GetBinContent(k)>0)
      pull /= TMath::Sqrt(hData->GetBinContent(k));
    hRatio->SetBinContent(k, pull);
    hError->SetBinContent(k, 0.0);
    float error2 = 0.0;
    error2 += pow(hZtt->GetBinContent(k) * uncertZtt, 2);
    error2 += pow(hTTb->GetBinContent(k) * uncertTTb, 2);
    error2 += pow(hEWK->GetBinContent(k) * uncertEWK, 2);
    error2 += pow(hQCD->GetBinContent(k) * uncertQCD, 2);
    if(hData->GetBinContent(k)>0)
      hError->SetBinError(k, TMath::Sqrt(error2/hData->GetBinContent(k)) );
    else
      hError->SetBinError(k,0);
    if(TMath::Abs(pull) > maxPull)
      maxPull = TMath::Abs(pull);
  }
  hRatio->SetAxisRange(-1.2*maxPull,1.2*maxPull,"Y");
  hRatio->Draw("P");
  hError->Draw("E3SAME");

  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hStack->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");
  
  //return;

  c1->SaveAs(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/plots/%s/plot_muTau_mH%d_%s_%s_%s.png",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()));
  c1->SaveAs(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/plots/%s/plot_muTau_mH%d_%s_%s_%s.pdf",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()));

  // templates for fitting
  TFile* fout = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/muTau_mH%d_%s_%s_%s.root",outputDir.Data(), mH_,selection_.c_str(),analysis_.c_str(),variable_.Data()),"RECREATE");
  fout->cd();

  hSiml->Write();
  hQCD->Write();
  hZmm->Write();
  hZmj->Write();
  hZfakes->Write();
  hTTb->Write();
  hZtt->Write();
  hDataEmb->Write();
  hLooseIso->Write();
  hAntiIso->Write();
  hW->Write();
  if(hWKeys) hW3Jets->Write();

  TH1* hWKeysHisto = hWKeys ? ((TH1Keys*)hWKeys)->GetHisto() : 0;
  if(hWKeysHisto){
    makeHistoFromDensity(hWKeysHisto,hW);
    hWKeysHisto->SetName("hWKeys");
    hWKeysHisto->Write();
  }
  TH1* hW3JetsKeysHisto = hW3JetsKeys ? ((TH1Keys*)hW3JetsKeys)->GetHisto() : 0;
  if(hW3JetsKeysHisto){
    makeHistoFromDensity(hW3JetsKeysHisto,hW3Jets);
    hW3JetsKeysHisto->SetName("hW3JetsKeys");
    hW3JetsKeysHisto->Write();
  }
  TH1* hZmmKeysHisto = hZmmKeys ? ((TH1Keys*)hZmmKeys)->GetHisto() : 0;
  if(hZmmKeysHisto){
    makeHistoFromDensity(hZmmKeysHisto,hZmm);
    hZmmKeysHisto->SetName("hZmmKeys");
    hZmmKeysHisto->Write();
  }

  TH1* hZmjKeysHisto = hZmjKeys ? ((TH1Keys*)hZmjKeys)->GetHisto() : 0;
  if(hZmjKeysHisto){
    makeHistoFromDensity(hZmjKeysHisto,hZmj);
    hZmjKeysHisto->SetName("hZmjKeys");
    hZmjKeysHisto->Write();
  }
  TH1* hLooseIsoKeysHisto = hLooseIsoKeys ? ((TH1Keys*)hLooseIsoKeys)->GetHisto() : 0;
  if(hLooseIsoKeysHisto){
    makeHistoFromDensity(hLooseIsoKeysHisto,hLooseIso);
    hLooseIsoKeysHisto->SetName("hLooseIsoKeys");
    hLooseIsoKeysHisto->Write();
  }
  TH1* hAntiIsoKeysHisto = hAntiIsoKeys ? ((TH1Keys*)hAntiIsoKeys)->GetHisto() : 0;
  if(hAntiIsoKeysHisto){
    makeHistoFromDensity(hAntiIsoKeysHisto,hAntiIso);
    hAntiIsoKeysHisto->SetName("hAntiIsoKeys");
    hAntiIsoKeysHisto->Write();
  }
  TH1* hVVKeysHisto = hVVKeys ? ((TH1Keys*)hVVKeys)->GetHisto() : 0;
  if(hVVKeysHisto){
    makeHistoFromDensity(hVVKeysHisto, hVV);
    hVVKeysHisto->SetName("hVVKeys");
    hVVKeysHisto->Write();
  }

  hVV->Write();
  hEWK->Write();
  hSgn1->Write();
  hSgn2->Write();
  hSgn3->Write();
  hData->Write();
  hParameters->Write();
  fout->Write();
  fout->Close();

  delete hQCD; delete hZmm; delete hZmj; delete hZfakes; delete hTTb; delete hZtt; delete hW; delete hW3Jets; delete hLooseIso; delete hAntiIso;
  if(hW3JetsKeys)   delete hW3JetsKeys;
  if(hWKeys)        delete hWKeys;
  if(hZmmKeys)      delete hZmmKeys;
  if(hZmjKeys)      delete hZmjKeys;
  if(hVVKeys)       delete hVVKeys;
  if(hLooseIsoKeys) delete hLooseIsoKeys;
  if(hAntiIsoKeys)  delete hAntiIsoKeys;
  delete hVV; delete hSgn1; delete hSgn2; delete hSgn3; delete hData; delete hParameters;
  delete hWMt; delete aStack;  delete hEWK; delete hSiml; delete hDataEmb; delete hSgn; delete hRatio; delete line;
  delete fout;


  fSignalGGH->Close();       delete fSignalGGH; 
  fSignalVBF->Close();       delete fSignalVBF;
  fSignalVH->Close();        delete fSignalVH; 
  fBackgroundOthers->Close();delete fBackgroundOthers;
  fBackgroundTTbar->Close(); delete fBackgroundTTbar;
  fBackgroundWJets->Close(); delete fBackgroundWJets;
  fData->Close();            delete fData; 
  dummy1->Close();           delete dummy1;
  fBackgroundDY->Close();    delete fBackgroundDY;

}


///////////////////////////////////////////////////////////////////////////////////////////////



void plotMuTauAll( Int_t useEmbedded = 1, TString outputDir = "May2012/Reload_PreApproval"){

  vector<string> variables;
  vector<int> mH;

  //variables.push_back("diTauVisMass");
  variables.push_back("diTauNSVfitMass");
 
  //mH.push_back(105);
  //mH.push_back(110);
  //mH.push_back(115);
  mH.push_back(120);
  //mH.push_back(125);
  //mH.push_back(130);
  //mH.push_back(135);
  //mH.push_back(140);
  //mH.push_back(145);
  //mH.push_back(160);

  //plotMuTau(120,1,"inclusive",""   ,"diTauVisMass","visible mass","GeV" ,outputDir,50,0,200,5.0,1.0,0,1.2);
  //plotMuTau(120,0,"inclusive",""   ,"MtLeg1MVA   ","M_{T}","GeV" ,       outputDir,40,0,160,5.0,1.0,0,1.2);
  //plotMuTau(120,1,"inclusive",""   ,"ptL2","#tau p_{T}","GeV"           ,outputDir,30,0, 90,5.0,1.0,0,1.2);
  //plotMuTau(120,1,"inclusive",""   ,"ptL1","#mu p_{T}", "GeV"           ,outputDir,30,0, 90,5.0,1.0,0,1.2);

  //plotMuTau(120,0,"inclusive",""   ,"numPV","reconstructed vertexes","units" ,outputDir,30,0,30,5.0,1.0,0,1.5);

  //plotMuTau(120,0,"inclusive",""   ,"etaL1","#mu #eta", "units"                          ,outputDir,25,-2.5, 2.5,5.0,1.0,0,2.);
  //plotMuTau(120,0,"inclusive",""   ,"etaL2","#tau #eta","units"                          ,outputDir,25,-2.5, 2.5,5.0,1.0,0,2.);
  //plotMuTau(120,0,"inclusive",""   ,"nJets30","jet multiplicity","units"                 ,outputDir,10,0, 10,5.0,1.0,1,10);
  //plotMuTau(120,0,"inclusive",""   ,"nJets20BTagged","b-tagged jet multiplicity","units" ,outputDir,5,0, 5,5.0,1.0,1,10);
  //plotMuTau(120,1,"oneJet",""      ,"pt1","leading jet p_{T}","GeV"       ,outputDir,36,20, 200,5.0,1.0,1,100);
  //plotMuTau(120,1,"oneJet",""      ,"eta1","leading jet #eta","units"     ,outputDir,21,-5, 5,5.0,1.0,0,2.);
  //plotMuTau(120,1,"twoJets",""     ,"pt1","leading jet p_{T}","GeV"       ,outputDir,36,20, 200,5.0,1.0,1,100);
  //plotMuTau(120,1,"twoJets",""     ,"pt2","trailing jet p_{T}","GeV"      ,outputDir,36,20, 200,5.0,1.0,1,100);
  //plotMuTau(120,1,"twoJets",""     ,"eta1","leading jet #eta","units"     ,outputDir,21,-5, 5,5.0,1.0,0,2.);
  //plotMuTau(120,1,"twoJets",""     ,"eta2","trailing jet #eta","units"    ,outputDir,21,-5, 5,5.0,1.0,0,2.);
  //plotMuTau(120,1,"twoJets",""     ,"Deta","|#Delta#eta|_{jj}","units"    ,outputDir,20,0, 8,   5.0,1.0,0,1.5);
  //plotMuTau(120,1,"twoJets",""     ,"Mjj","M_{jj}","GeV"                  ,outputDir,20,0, 1000,5.0,1.0,1,100);
  //plotMuTau(120,1,"twoJets",""     ,"MVAvbf","BDT output","units"         ,outputDir,20,-1, 1,5.0,1.0,1,100);
  
  //return;

  for(unsigned int i = 0 ; i < variables.size(); i++){
    for(unsigned j = 0; j < mH.size(); j++){

      //plotMuTau(mH[j],useEmbedded,"inclusive",""   ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      
      plotMuTau(mH[j],useEmbedded,"novbfHigh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
 
      return;

      plotMuTau(mH[j],useEmbedded,"novbfLow",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfLow","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfLow","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfLow","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfLow","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"novbfHigh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfHigh","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfHigh","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfHigh","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbfHigh","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"boostLow",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostLow","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostLow","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostLow","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostLow","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"boostHigh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostHigh","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostHigh","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostHigh","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boostHigh","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"bTagLow",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagLow","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagLow","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagLow","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagLow","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"bTagHigh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagHigh","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagHigh","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagHigh","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"bTagHigh","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"vbf",""         ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","TauUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","TauDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","JetUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","JetDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"vh",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vh","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vh","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vh","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vh","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);



      /*
      plotMuTau(mH[j],useEmbedded,"novbf","Raw"       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbf","RawTauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbf","RawTauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbf","RawJetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"novbf","RawJetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      
      plotMuTau(mH[j],useEmbedded,"vbf","Raw"         ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","RawTauUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","RawTauDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","RawJetUp"    ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"vbf","RawJetDown"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      plotMuTau(mH[j],useEmbedded,"boost","Raw"       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boost","RawTauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boost","RawTauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boost","RawJetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      plotMuTau(mH[j],useEmbedded,"boost","RawJetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      */




      //plotMuTau(mH[j],useEmbedded,"boost2",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"boost2","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"boost2","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"boost2","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"boost2","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,1.0,1.0,0,1.2);

      //plotMuTau(mH[j],useEmbedded,"twoJets",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"twoJets","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"twoJets","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"twoJets","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"twoJets","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      
      //plotMuTau(mH[j],useEmbedded,"oneJet",""       ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"oneJet","TauUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"oneJet","TauDown",variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"oneJet","JetUp"  ,variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      //plotMuTau(mH[j],useEmbedded,"oneJet","JetDown",variables[i],"mass","GeV",outputDir,-1,0,100,5.0,1.0,0,1.2);
      
    }
  }
  


}




int main(int argc, const char* argv[])
{

  std::cout << "plotMuTau()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  plotMuTauAll();
  //plotMuTauA();
  //plotMuTau(120,1,"",""   ,"MtLeg1Corr","M_{T}","GeV" ,"./",40,0,160,5.0,1.0,0,1.2);

  //plotMuTau(120,1,"inclusive","","diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotMuTau(120,1,"novbf","","diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotMuTau(120,1,"boost","","diTauVisMass","","","Dec2011/iter3",10,0,7000);
  //plotMuTau(120,1,"vbf","","diTauVisMass","","","Dec2011/iter3",10,0,7000);

}
