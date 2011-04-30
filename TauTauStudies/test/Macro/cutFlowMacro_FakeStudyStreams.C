#include <cstdlib>
#include <iostream> 
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
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

#define VERBOSE true

void cutFlowStudyElec( TString weightFile = "TMVAClassificationPtOrd_qqH115vsWZttQCD_Cuts.weights.xml",
		       Double_t effS_ = 0.3) 
{
  
  ofstream out("tmp.txt");
  out.precision(4);

  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   

  Float_t pt1, pt2;
  Float_t Deta, Mjj;
  Float_t eta1,eta2;

  reader->AddVariable( "pt1", &pt1);
  reader->AddVariable( "pt2", &pt2);
  reader->AddVariable( "Deta",&Deta);
  reader->AddVariable( "Mjj", &Mjj);
  reader->AddSpectator("eta1",&eta1);
  reader->AddSpectator("eta2",&eta2);
  reader->BookMVA( "Cuts", TString("/home/llr/cms/lbianchini/CMSSW_3_9_9/src/Bianchi/TauTauStudies/test/Macro/weights/")+weightFile ); 
 
  TFile *fFullSignalVBF           = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_VBFH115-powheg-PUS1.root","READ"); 
  TFile *fFullSignalGGH           = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_GGFH115-powheg-PUS1.root","READ");  
  TFile *fFullBackgroundDYTauTau  = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_DYToTauTau-20-PUS1.root","READ"); 
  TFile *fFullBackgroundDYEleEle  = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_DYToEE-20-PUS1.root","READ"); 
  TFile *fFullBackgroundWJets     = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_WJets-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundQCD       = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_QCD.root","READ"); 
  TFile *fFullBackgroundG1J       = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_G1Jet.root","READ"); 
  TFile *fFullBackgroundTTbar     = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_TTJets-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundDiBoson   = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011/treeElecTauStream_DiBoson.root","READ"); 

  // OpenNTuples
  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleVBFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleGGFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDYToTauTau-20-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYEleEle  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDYToEE-20-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleWJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameQCD       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleQCD_Open_ElecTauStream.root";
  TString fBackgroundNameG1J       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleG1Jet_Open_ElecTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleTTJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDiBoson   = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDiBoson_Open_ElecTauStream.root";

  TFile *fSignalVBF(0); 
  TFile *fSignalGGH(0); 
  TFile *fBackgroundDYTauTau(0);
  TFile *fBackgroundDYEleEle(0);
  TFile *fBackgroundWJets(0);
  TFile *fBackgroundQCD(0);
  TFile *fBackgroundG1J(0);
  TFile *fBackgroundTTbar(0);
  TFile *fBackgroundDiBoson(0);
 
  fSignalVBF          = TFile::Open( fSignalNameVBF ); 
  fSignalGGH          = TFile::Open( fSignalNameGGH ); 
  fBackgroundDYTauTau = TFile::Open( fBackgroundNameDYTauTau ); 
  fBackgroundDYEleEle = TFile::Open( fBackgroundNameDYEleEle ); 
  fBackgroundWJets    = TFile::Open( fBackgroundNameWJets ); 
  fBackgroundQCD      = TFile::Open( fBackgroundNameQCD ); 
  fBackgroundG1J      = TFile::Open( fBackgroundNameG1J ); 
  fBackgroundTTbar    = TFile::Open( fBackgroundNameTTbar ); 
  fBackgroundDiBoson  = TFile::Open( fBackgroundNameDiBoson ); 

  if(!fSignalVBF || !fBackgroundDYTauTau || !fBackgroundWJets || !fBackgroundQCD || !fBackgroundTTbar ||
     !fSignalGGH || !fBackgroundDYEleEle || !fBackgroundG1J || !fBackgroundDiBoson ){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";

  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);
  TTree *backgroundDYTauTau  = (TTree*)fBackgroundDYTauTau->Get(tree);
  TTree *backgroundDYEleEle  = (TTree*)fBackgroundDYEleEle->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundQCD       = (TTree*)fBackgroundQCD->Get(tree);
  TTree *backgroundG1J       = (TTree*)fBackgroundG1J->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundDiBoson   = (TTree*)fBackgroundDiBoson->Get(tree);

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["Zeleele"]=backgroundDYEleEle;
  tMap["Wjets"]=backgroundWJets;
  tMap["G1J"]=backgroundG1J;
  tMap["QCD"]=backgroundQCD;
  tMap["TTbar"]=backgroundTTbar;
  tMap["DiBoson"]=backgroundDiBoson;

  std::map<std::string,TTree*>::iterator jt;

  Float_t pt1_, pt2_;
  Float_t Deta_, Mjj_;
  Float_t Dphi,diTauSVFitPt,diTauSVFitEta,diTauVisMass,diTauSVFitMass,ptL1,ptL2,etaL1,etaL2,diTauCharge,MtLeg1,numPV,combRelIsoLeg1,sampleWeight,ptVeto;
  Int_t tightestHPSWP;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  samples.push_back("DiBoson");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("Zeleele");
  samples.push_back("Ztautau");
  samples.push_back("G1J");
  samples.push_back("QCD");

  std::map<std::string,float> crossSec;
  crossSec["ggH115"]=( 7.65e-02 * 18.13 );
  crossSec["qqH115"]=( 0.1012);
  crossSec["DiBoson"]=( -1  );
  crossSec["TTbar"]=( 157.5 );
  crossSec["Wjets"]=( 31314.0);
  crossSec["Zeleele"]=( 1666  );
  crossSec["Ztautau"]=( 1666  );
  crossSec["G1J"]=( -1 );
  crossSec["QCD"]=( -1 );

  float Lumi = 1000;


  // here I choose the order in the stack
  std::vector<string> filters;
  filters.push_back("Total");
  filters.push_back("Good vertex");
  filters.push_back("Exactly 1 WP95 electron, p_{T}>15 GeV");
  filters.push_back("No VBTF muon, p_{T}>15 GeV");
  filters.push_back("One WP80 electron");
  filters.push_back("One HPS-id tau, p_{T}>20 GeV");
  filters.push_back("#Delta R ele-tau>0.3");
  filters.push_back("electron isolation");
  filters.push_back("tau isolation");
  filters.push_back("M_{T}(electron,MET)<40 GeV");
  filters.push_back("OS");
  filters.push_back("2 jets, E_{T}>20,15 GeV, #eta_{1}*#eta_{2}<0");
  filters.push_back("VBF cuts");

  // here I define the map between a sample name and its file ptr
  std::map<std::string,TFile*> fullMap;
  fullMap["ggH115"]   = fFullSignalGGH;
  fullMap["qqH115"]   = fFullSignalVBF;
  fullMap["Ztautau"]  = fFullBackgroundDYTauTau;
  fullMap["Zeleele"]  = fFullBackgroundDYEleEle;
  fullMap["Wjets"]    = fFullBackgroundWJets;
  fullMap["G1J"]      = fFullBackgroundG1J;
  fullMap["QCD"]      = fFullBackgroundQCD;
  fullMap["TTbar"]    = fFullBackgroundTTbar;
  fullMap["DiBoson"]  = fFullBackgroundDiBoson;

  std::map<std::string,TFile*>::iterator it;

  std::map<std::string,float> cutMap_allEventsFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    float tot = totalEvents;
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_allEventsFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_vertexScrapingFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("vertexScrapingFilter/totalEvents");
    float tot =  allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_vertexScrapingFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_oneElectronFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("oneElectronFilter/totalEvents");
    float tot =  allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_oneElectronFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_noMuonFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("noMuonFilter/totalEvents");
    float tot =  allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_noMuonFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_electronLegFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("electronLegFilter/totalEvents");
    float tot =  allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_electronLegFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_tauLegFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("tauLegFilter/totalEvents");
    float tot =  allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_tauLegFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_atLeastOneDiTauFilter;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    allEvents = (TH1F*)(it->second)->Get("atLeastOneDiTauFilter/totalEvents");
    float tot = allEvents->GetBinContent(1);
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_atLeastOneDiTauFilter[it->first] = tot;
  }
  std::map<std::string,float> cutMap_ElecIso;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    cout<<it->first<<endl;
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (crossSec[it->first]>0) ?  "(chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1"
      : "weight*((chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1)";
    ((TTree*) (it->second->Get("elecTauStreamAnalyzer/tree")) )->Draw("diTauLegsP4[0].Eta()>>h1",cut);
    float tot = h1->Integral();
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_ElecIso[it->first] = tot;
    delete h1;
  }
  std::map<std::string,float> cutMap_TauIso;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    cout<<it->first<<endl;
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (crossSec[it->first]>0) ?  "(chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0"
      : "weight*((chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0)";
    ((TTree*) (it->second->Get("elecTauStreamAnalyzer/tree")) )->Draw("diTauLegsP4[0].Eta()>>h1",cut);
    float tot = h1->Integral();
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_TauIso[it->first] = tot;
    delete h1;
  }
  std::map<std::string,float> cutMap_OS;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    cout<<it->first<<endl;
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (crossSec[it->first]>0) ?  "(chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0 && diTauCharge==0"
      : "weight*((chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0 && diTauCharge==0)";
    ((TTree*) (it->second->Get("elecTauStreamAnalyzer/tree")) )->Draw("diTauLegsP4[0].Eta()>>h1",cut);
    float tot = h1->Integral();
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_OS[it->first] = tot;
    delete h1;
  }
  std::map<std::string,float> cutMap_Mt;
  for(it = fullMap.begin(); it != fullMap.end(); it++){
    cout<<it->first<<endl;
    TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (crossSec[it->first]>0) ?  "(chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0 && diTauCharge==0 && MtLeg1<40"
      : "weight*((chIsoLeg1+nhIsoLeg1+phIsoLeg1)/diTauLegsP4[0].Pt()<0.1 && tightestHPSWP>0 && diTauCharge==0 && MtLeg1<40)";
    ((TTree*) (it->second->Get("elecTauStreamAnalyzer/tree")) )->Draw("diTauLegsP4[0].Eta()>>h1",cut);
    float tot = h1->Integral();
    if(crossSec[it->first]>0) tot *= (Lumi/ (totalEvents/crossSec[it->first])  ) ;
    cutMap_Mt[it->first] = tot;
    delete h1;
  }
  std::map<std::string,float> cutMap_VBFPre;
  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*(pt1>0 && combRelIsoLeg1<0.1 && tightestHPSWP>0 && diTauCharge==0 && MtLeg1<40)"; 
    jt->second->Draw("etaL1>>h1",cut);
    float tot = h1->Integral();
    cutMap_VBFPre[jt->first] = tot;
    delete h1;
  }
  std::map<std::string,float> cutMap_VBF;
  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TCut cut =  "(pt1>0 && combRelIsoLeg1<0.1 && tightestHPSWP>0 && diTauCharge==0 && MtLeg1<40)"; 

    TFile* dummy = new TFile("dummy.root","RECREATE");  
    TTree* currentTree = (TTree*)(jt->second)->CopyTree(cut);
    float tot = 0;

    currentTree->SetBranchAddress( "pt1", &pt1_ );
    currentTree->SetBranchAddress( "pt2", &pt2_ );
    currentTree->SetBranchAddress( "Deta",&Deta_ );
    currentTree->SetBranchAddress( "Mjj", &Mjj_ );
    currentTree->SetBranchAddress( "diTauSVFitPt",&diTauSVFitPt);
    //currentTree->SetBranchAddress( "diTauSVFitEta",&diTauSVFitEta);
    currentTree->SetBranchAddress( "diTauSVFitMass",&diTauSVFitMass);
    currentTree->SetBranchAddress( "diTauVisMass",&diTauVisMass);
    currentTree->SetBranchAddress( "ptL1", &ptL1 );
    currentTree->SetBranchAddress( "ptL2",  &ptL2 );
    currentTree->SetBranchAddress( "etaL1", &etaL1 );
    currentTree->SetBranchAddress( "etaL2", &etaL2 );
    currentTree->SetBranchAddress( "combRelIsoLeg1",&combRelIsoLeg1);
    currentTree->SetBranchAddress( "tightestHPSWP",&tightestHPSWP);
    currentTree->SetBranchAddress( "diTauCharge",&diTauCharge);
    currentTree->SetBranchAddress( "MtLeg1",&MtLeg1);
    currentTree->SetBranchAddress( "numPV",&numPV);
    currentTree->SetBranchAddress( "sampleWeight",&sampleWeight);
    currentTree->SetBranchAddress( "ptVeto",&ptVeto);

    for (Long64_t ievt=0; ievt<currentTree->GetEntries();ievt++) {

      currentTree->GetEntry(ievt);

      if (ievt%10000 == 0){
	std::cout << (jt->first) << " ---> processing event: " << ievt << " ..." <<std::endl;
      }

      pt1  = pt1_;
      pt2  = pt2_;
      Deta = Deta_;
      Mjj  = Mjj_;

      bool pass = effS_>0 ? reader->EvaluateMVA( "Cuts", effS_ ) : (pt1>0);

      if(pass)	tot+=sampleWeight;
    
    }// end loop   

    cutMap_VBF[jt->first] = tot;
  }


  
 

  std::vector< std::map<std::string,float> > allFilters;
  allFilters.push_back(cutMap_allEventsFilter);
  allFilters.push_back(cutMap_vertexScrapingFilter);
  allFilters.push_back(cutMap_oneElectronFilter);
  allFilters.push_back(cutMap_noMuonFilter);
  allFilters.push_back(cutMap_electronLegFilter);
  allFilters.push_back(cutMap_tauLegFilter);
  allFilters.push_back(cutMap_atLeastOneDiTauFilter);
  allFilters.push_back(cutMap_ElecIso);
  allFilters.push_back(cutMap_TauIso);
  allFilters.push_back(cutMap_Mt);
  allFilters.push_back(cutMap_OS);
  allFilters.push_back(cutMap_VBFPre);
  allFilters.push_back(cutMap_VBF);

  out<<"\\begin{center}"<<endl;
  out<<"\\begin{tabular}[!htbp]{|c";
  for(int k = 0 ; k < samples.size(); k++) out<<"|c";
  out<<"|} \\hline"<<endl;
  out<< "Cut & ";
  for(int k = 0 ; k < samples.size(); k++){
    out << (fullMap.find(samples[k]))->first;
    if(k!=samples.size()-1) out <<" & " ;
    else out << " \\hline" << endl;
  }

  
  for(int i = 0; i < allFilters.size(); i++){
    out << filters[i] << " & ";
    for(int k = 0 ; k < samples.size(); k++){
      out << (allFilters[i].find(samples[k]))->second ;
      if(k!=samples.size()-1) out <<" & " ;
      else out << " \\hline" << endl;
    }

  }
  
  out<<"\\end{tabular}"<<endl;
  out<<"\\end{center}"<<endl;
 

  return;

  /*
  TMVA::Tools::Instance();

  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   

  Float_t pt1, pt2;
  Float_t Deta, Mjj;
  Float_t eta1,eta2;

  reader->AddVariable( "pt1", &pt1);
  reader->AddVariable( "pt2", &pt2);
  reader->AddVariable( "Deta",&Deta);
  reader->AddVariable( "Mjj", &Mjj);

  reader->AddSpectator("eta1",&eta1);
  reader->AddSpectator("eta2",&eta2);

  reader->BookMVA( "Cuts", TString("weights/")+weightFile ); 

  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleVBFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleGGFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDYToTauTau-20-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYEleEle  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDYToEE-20-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleWJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameQCD       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleQCD_Open_ElecTauStream.root";
  TString fBackgroundNameG1J       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleG1Jet_Open_ElecTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleTTJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDiBoson   = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011/nTupleDiBoson_Open_ElecTauStream.root";

  TFile *fSignalVBF(0); 
  TFile *fSignalGGH(0); 
  TFile *fBackgroundDYTauTau(0);
  TFile *fBackgroundDYEleEle(0);
  TFile *fBackgroundWJets(0);
  TFile *fBackgroundQCD(0);
  TFile *fBackgroundG1J(0);
  TFile *fBackgroundTTbar(0);
  TFile *fBackgroundDiBoson(0);
 

  fSignalVBF          = TFile::Open( fSignalNameVBF ); 
  fSignalGGH          = TFile::Open( fSignalNameGGH ); 
  fBackgroundDYTauTau = TFile::Open( fBackgroundNameDYTauTau ); 
  fBackgroundDYEleEle = TFile::Open( fBackgroundNameDYEleEle ); 
  fBackgroundWJets    = TFile::Open( fBackgroundNameWJets ); 
  fBackgroundQCD      = TFile::Open( fBackgroundNameQCD ); 
  fBackgroundG1J      = TFile::Open( fBackgroundNameG1J ); 
  fBackgroundTTbar    = TFile::Open( fBackgroundNameTTbar ); 
  fBackgroundDiBoson  = TFile::Open( fBackgroundNameDiBoson ); 

  if(!fSignalVBF || !fBackgroundDYTauTau || !fBackgroundWJets || !fBackgroundQCD || !fBackgroundTTbar ||
     !fSignalGGH || !fBackgroundDYEleEle || !fBackgroundG1J || !fBackgroundDiBoson ){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";

  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);
  TTree *backgroundDYTauTau  = (TTree*)fBackgroundDYTauTau->Get(tree);
  TTree *backgroundDYEleEle  = (TTree*)fBackgroundDYEleEle->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundQCD       = (TTree*)fBackgroundQCD->Get(tree);
  TTree *backgroundG1J       = (TTree*)fBackgroundG1J->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundDiBoson   = (TTree*)fBackgroundDiBoson->Get(tree);
 
  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  samples.push_back("DiBoson");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("Zeleele");
  samples.push_back("Ztautau");
  samples.push_back("G1J");
  samples.push_back("QCD");

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["Zeleele"]=backgroundDYEleEle;
  tMap["Wjets"]=backgroundWJets;
  tMap["G1J"]=backgroundG1J;
  tMap["QCD"]=backgroundQCD;
  tMap["TTbar"]=backgroundTTbar;
  tMap["DiBoson"]=backgroundDiBoson;

  Float_t pt1_, pt2_;
  Float_t Deta_, Mjj_;
  Float_t Dphi,diTauSVFitPt,diTauSVFitEta,diTauVisMass,diTauSVFitMass,ptL1,ptL2,etaL1,etaL2,diTauCharge,MtLeg1,numPV,combRelIsoLeg1,sampleWeight,ptVeto;
  Int_t tightestHPSWP;

  std::map<TString,Float_t> vMap;

  for( unsigned iter=0; iter<samples.size(); iter++){

    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,"" , nBins_ ,xMin_ , xMax_);

    TFile* dummy = new TFile("dummy.root","RECREATE");  
    TTree* currentTree = (TTree*)(it->second)->CopyTree(Cuts_);
    Float_t counter = 0;

    currentTree->SetBranchAddress( "pt1", &pt1_ );
    currentTree->SetBranchAddress( "pt2", &pt2_ );
    currentTree->SetBranchAddress( "Deta",&Deta_ );
    currentTree->SetBranchAddress( "Mjj", &Mjj_ );
    currentTree->SetBranchAddress( "diTauSVFitPt",&diTauSVFitPt);
    //currentTree->SetBranchAddress( "diTauSVFitEta",&diTauSVFitEta);
    currentTree->SetBranchAddress( "diTauSVFitMass",&diTauSVFitMass);
    currentTree->SetBranchAddress( "diTauVisMass",&diTauVisMass);
    currentTree->SetBranchAddress( "ptL1", &ptL1 );
    currentTree->SetBranchAddress( "ptL2",  &ptL2 );
    currentTree->SetBranchAddress( "etaL1", &etaL1 );
    currentTree->SetBranchAddress( "etaL2", &etaL2 );
    currentTree->SetBranchAddress( "combRelIsoLeg1",&combRelIsoLeg1);
    currentTree->SetBranchAddress( "tightestHPSWP",&tightestHPSWP);
    currentTree->SetBranchAddress( "diTauCharge",&diTauCharge);
    currentTree->SetBranchAddress( "MtLeg1",&MtLeg1);
    currentTree->SetBranchAddress( "numPV",&numPV);
    currentTree->SetBranchAddress( "sampleWeight",&sampleWeight);
    currentTree->SetBranchAddress( "ptVeto",&ptVeto);

    for (Long64_t ievt=0; ievt<currentTree->GetEntries();ievt++) {

      currentTree->GetEntry(ievt);

      if (ievt%10000 == 0){
	std::cout << (it->first) << " ---> processing event: " << ievt << " ..." <<std::endl;
      }

      pt1  = pt1_;
      pt2  = pt2_;
      Deta = Deta_;
      Mjj  = Mjj_;

      vMap["Dphi"]= Dphi;
      vMap["diTauSVFitPt"]= diTauSVFitPt;
      //vMap["diTauSVFitEta"]= diTauSVFitEta;
      vMap["diTauSVFitMass"]= diTauSVFitMass;
      vMap["diTauVisMass"]= diTauVisMass;
      vMap["ptL1"]= ptL1;
      vMap["ptL2"]= ptL2;
      vMap["etaL1"]= etaL1;
      vMap["etaL2"]= etaL2;
      vMap["diTauCharge"]= Float_t(diTauCharge);
      vMap["MtLeg1"]= MtLeg1;
      vMap["numPV"]= numPV;
      vMap["combRelIsoLeg1"]= combRelIsoLeg1;
      vMap["sampleWeight"]= sampleWeight;
      vMap["pt1_"]= pt1;
      vMap["pt2_"]= pt2;
      vMap["Deta_"]= Deta;
      vMap["Mjj_"]= Mjj;
      vMap["ptVeto"]= ptVeto;

      // if efficiency is more than 0 use the optimized cuts
      // else use VBF  preselection
      bool pass = effS_>0 ? reader->EvaluateMVA( "Cuts", effS_ ) : (pt1>0);

      if(pass){
	counter+=sampleWeight;
      }
    
    }// end loop   

    cout << it->first << " ==> " << counter  << endl;

  }
  

  delete reader;
  */
}


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


void fakeStudyMu( TString weightFile = "TMVAClassificationPtOrd_qqH115vsWZttQCD_Cuts.weights.xml",
		  Double_t effS_ = 0.5,
		  TCut Cuts_ = "pt1>0",
		  TString variable_ = "MtLeg1",
		  TString XTitle_ = "M_{T}(#mu,MET)",
		  TString Unities_ = "GeV",
		  Int_t nBins_ = 12, Float_t xMin_=0, Float_t xMax_=120,
		  Float_t magnifySgn_ = 10,
		  Int_t logy_ = 0 ) 
{   
  
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);

  TLegend* leg = new TLegend(0.60,0.47,0.90,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader("#mu+#tau_{had}");

  THStack* aStack = new THStack("aStack",Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}", 1000.));
  TH1F* hSiml = new TH1F();
  TH1F* hSgn   = new TH1F();
  TH1F* hSgn1  = new TH1F();
  TH1F* hSgn2  = new TH1F();


  TMVA::Tools::Instance();

  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   

  Float_t pt1, pt2;
  Float_t Deta, Mjj;
  Float_t eta1,eta2;

  reader->AddVariable( "pt1", &pt1);
  reader->AddVariable( "pt2", &pt2);
  reader->AddVariable( "Deta",&Deta);
  reader->AddVariable( "Mjj", &Mjj);

  reader->AddSpectator("eta1",&eta1);
  reader->AddSpectator("eta2",&eta2);

  reader->BookMVA( "Cuts", TString("weights/")+weightFile ); 

  // define histos

  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleVBFH115-Mu-powheg-PUS1_Open_MuTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleGGFH115-Mu-powheg-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleDYToTauTau-Mu-20-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameDYMuMu    = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleDYToMuMu-20-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleWJets-Mu-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameQCD       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleQCDmu_Open_MuTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleTTJets-Mu-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameDiBoson   = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011/nTupleDiBoson-Mu_Open_MuTauStream.root";


  TFile *fSignalVBF(0); 
  TFile *fSignalGGH(0); 
  TFile *fBackgroundDYTauTau(0);
  TFile *fBackgroundDYMuMu(0);
  TFile *fBackgroundWJets(0);
  TFile *fBackgroundQCD(0);
  TFile *fBackgroundTTbar(0);
  TFile *fBackgroundDiBoson(0);
 

  fSignalVBF          = TFile::Open( fSignalNameVBF ); 
  fSignalGGH          = TFile::Open( fSignalNameGGH ); 
  fBackgroundDYTauTau = TFile::Open( fBackgroundNameDYTauTau ); 
  fBackgroundDYMuMu   = TFile::Open( fBackgroundNameDYMuMu ); 
  fBackgroundWJets    = TFile::Open( fBackgroundNameWJets ); 
  fBackgroundQCD      = TFile::Open( fBackgroundNameQCD ); 
  fBackgroundTTbar    = TFile::Open( fBackgroundNameTTbar ); 
  fBackgroundDiBoson  = TFile::Open( fBackgroundNameDiBoson ); 

  if(!fSignalVBF || !fBackgroundDYTauTau || !fBackgroundWJets || !fBackgroundQCD || !fBackgroundTTbar ||
     !fSignalGGH || !fBackgroundDYMuMu || !fBackgroundDiBoson ){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";

  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);
  TTree *backgroundDYTauTau  = (TTree*)fBackgroundDYTauTau->Get(tree);
  TTree *backgroundDYMuMu    = (TTree*)fBackgroundDYMuMu->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundQCD       = (TTree*)fBackgroundQCD->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundDiBoson   = (TTree*)fBackgroundDiBoson->Get(tree);
 
  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  samples.push_back("DiBoson");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("Zmumu");
  samples.push_back("Ztautau");
  samples.push_back("QCD");

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["Zmumu"]=backgroundDYMuMu;
  tMap["Wjets"]=backgroundWJets;
  tMap["QCD"]=backgroundQCD;
  tMap["TTbar"]=backgroundTTbar;
  tMap["DiBoson"]=backgroundDiBoson;

  Float_t pt1_, pt2_;
  Float_t Deta_, Mjj_;
  Float_t Dphi,diTauSVFitPt,diTauSVFitEta,diTauVisMass,diTauSVFitMass,ptL1,ptL2,etaL1,etaL2,diTauCharge,MtLeg1,numPV,combRelIsoLeg1,sampleWeight,ptVeto;
  Int_t tightestHPSWP;

  std::map<TString,Float_t> vMap;

  for( unsigned iter=0; iter<samples.size(); iter++){

    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,"" , nBins_ ,xMin_ , xMax_);

    TFile* dummy = new TFile("dummy.root","RECREATE");  
    TTree* currentTree = (TTree*)(it->second)->CopyTree(Cuts_);
    Int_t counter = 0;

    currentTree->SetBranchAddress( "pt1", &pt1_ );
    currentTree->SetBranchAddress( "pt2", &pt2_ );
    currentTree->SetBranchAddress( "Deta",&Deta_ );
    currentTree->SetBranchAddress( "Mjj", &Mjj_ );
    currentTree->SetBranchAddress( "diTauSVFitPt",&diTauSVFitPt);
    //currentTree->SetBranchAddress( "diTauSVFitEta",&diTauSVFitEta);
    currentTree->SetBranchAddress( "diTauSVFitMass",&diTauSVFitMass);
    currentTree->SetBranchAddress( "diTauVisMass",&diTauVisMass);
    currentTree->SetBranchAddress( "ptL1", &ptL1 );
    currentTree->SetBranchAddress( "ptL2",  &ptL2 );
    currentTree->SetBranchAddress( "etaL1", &etaL1 );
    currentTree->SetBranchAddress( "etaL2", &etaL2 );
    currentTree->SetBranchAddress( "combRelIsoLeg1",&combRelIsoLeg1);
    currentTree->SetBranchAddress( "tightestHPSWP",&tightestHPSWP);
    currentTree->SetBranchAddress( "diTauCharge",&diTauCharge);
    currentTree->SetBranchAddress( "MtLeg1",&MtLeg1);
    currentTree->SetBranchAddress( "numPV",&numPV);
    currentTree->SetBranchAddress( "sampleWeight",&sampleWeight);
    currentTree->SetBranchAddress( "ptVeto",&ptVeto);

    for (Long64_t ievt=0; ievt<currentTree->GetEntries();ievt++) {

      currentTree->GetEntry(ievt);

      if (ievt%10000 == 0){
	std::cout << (it->first) << " ---> processing event: " << ievt << " ..." <<std::endl;
      }

      pt1  = pt1_;
      pt2  = pt2_;
      Deta = Deta_;
      Mjj  = Mjj_;

      vMap["Dphi"]= Dphi;
      vMap["diTauSVFitPt"]= diTauSVFitPt;
      //vMap["diTauSVFitEta"]= diTauSVFitEta;
      vMap["diTauSVFitMass"]= diTauSVFitMass;
      vMap["diTauVisMass"]= diTauVisMass;
      vMap["ptL1"]= ptL1;
      vMap["ptL2"]= ptL2;
      vMap["etaL1"]= etaL1;
      vMap["etaL2"]= etaL2;
      vMap["diTauCharge"]= Float_t(diTauCharge);
      vMap["MtLeg1"]= MtLeg1;
      vMap["numPV"]= numPV;
      vMap["combRelIsoLeg1"]= combRelIsoLeg1;
      vMap["sampleWeight"]= sampleWeight;
      vMap["pt1_"]= pt1;
      vMap["pt2_"]= pt2;
      vMap["Deta_"]= Deta;
      vMap["Mjj_"]= Mjj;
      vMap["ptVeto"]= ptVeto;

      // if efficiency is more than 0 use the optimized cuts
      // else use VBF  preselection
      bool pass = effS_>0 ? reader->EvaluateMVA( "Cuts", effS_ ) : (pt1>0);

      if(pass){
	counter++;
	h1->Fill( vMap[variable_], sampleWeight);
      }
    
    }// end loop


    if( (it->first).find("Zmumu")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z#rightarrow #mu#mu","F");
    }
    if( (it->first).find("Ztautau")!=string::npos ) {
      h1->SetFillColor(kRed);
      leg->AddEntry(h1,"Z#rightarrow #tau#tau","F");
    }
    if( (it->first).find("TTbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"t#bar{t}+jets","F");
    }
    if( (it->first).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"W+jets","F");
    }
    if( (it->first).find("DiBoson")!=string::npos ) {
      h1->SetFillColor(38);
      leg->AddEntry(h1,"Di-Boson","F");
    }
    if( (it->first).find("QCD")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"QCD-multijets","F");
    }
    if((it->first).find("qqH115")!=string::npos){
      hSgn1 = (TH1F*)h1->Clone("hSgn1");
      hSgn1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("VBF H(115) X %.0f",magnifySgn_),"F");
    }
    if((it->first).find("ggH115")!=string::npos){
      hSgn2 = (TH1F*)h1->Clone("hSgn2");
      hSgn2->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("GGF H(115) X %.0f",magnifySgn_),"F");
    }

    if(iter==0) hSiml=(TH1F*)h1->Clone("hSiml");
    else hSiml->Add(h1);

    aStack->Add(h1);

    if(VERBOSE) cout<<(it->first) << " ==> " 
		   << h1->Integral() << " +/- " 
		   << TMath::Sqrt(h1->GetEntries())*(h1->Integral()/h1->GetEntries())
		   << endl;
  }

  // all signal summed together:
  hSgn = (TH1F*)hSgn1->Clone("hSgn");
  hSgn->Add(hSgn1,hSgn2,1,1);

  if(VERBOSE) cout<< "S/sqrt(B) ==> " 
		  << hSgn->Integral()/ TMath::Sqrt(hSiml->Integral()) << " +/- " 
		  << (1./2)*TMath::Sqrt(hSiml->GetEntries())*(hSiml->GetSumOfWeights())/hSiml->Integral()*( hSgn->Integral()/ TMath::Sqrt(hSiml->Integral())  )
		  << endl;

  aStack->Draw("HIST");
  //hSgn1->Draw("HISTSAME");
  //hSgn2->Draw("HISTSAME");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  leg->Draw();


  delete reader;

}



