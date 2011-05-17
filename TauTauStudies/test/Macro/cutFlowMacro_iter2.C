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
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"


#define VERBOSE true

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void cutFlowStudyElec( Double_t hltEff_=0.8 ) 
{
  
  ofstream out("cutFlow_ElecTauStream_iter2.txt");
  out.precision(4);


  TFile *fFullData                = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_Run2011-Elec.root","READ"); 
  TFile *fFullSignalVBF           = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_VBFH115-powheg-PUS1.root","READ"); 
  TFile *fFullSignalGGH           = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_GGFH115-powheg-PUS1.root","READ");  
  TFile *fFullBackgroundDYTauTau  = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_DYJets-50-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundDYEleEle  = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_DYJets-50-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundWJets     = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_WJets-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundQCD       = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_QCD.root","READ"); 
  TFile *fFullBackgroundG1J       = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_G1Jet.root","READ"); 
  TFile *fFullBackgroundTTbar     = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_TTJets-madgraph-PUS1.root","READ");
  TFile *fFullBackgroundSingleTop = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_SingleTop.root","READ");
  TFile *fFullBackgroundDiBoson   = new TFile("/data_CMS/cms/lbianchini//ElecTauStream2011_iter2/treeElecTauStream_DiBoson.root","READ"); 



  // OpenNTuples
  TString fDataName                = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleRun2011-Elec_Open_ElecTauStream.root";
  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleVBFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleGGFH115-powheg-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleDYJets-50-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameDYEleEle  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleDYJets-50-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleWJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameQCD       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleQCD_Open_ElecTauStream.root";
  TString fBackgroundNameG1J       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleG1Jet_Open_ElecTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleTTJets-madgraph-PUS1_Open_ElecTauStream.root";
  TString fBackgroundNameSingleTop = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleSingleTop_Open_ElecTauStream.root";
  TString fBackgroundNameDiBoson   = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStream2011_iter2/Inclusive/nTupleDiBoson_Open_ElecTauStream.root";

  TFile *fData(0); 
  TFile *fSignalVBF(0); 
  TFile *fSignalGGH(0); 
  TFile *fBackgroundDYTauTau(0);
  TFile *fBackgroundDYEleEle(0);
  TFile *fBackgroundWJets(0);
  TFile *fBackgroundQCD(0);
  TFile *fBackgroundG1J(0);
  TFile *fBackgroundTTbar(0);
  TFile *fBackgroundSingleTop(0);
  TFile *fBackgroundDiBoson(0);
 
  fData               = TFile::Open( fDataName ); 
  fSignalVBF          = TFile::Open( fSignalNameVBF ); 
  fSignalGGH          = TFile::Open( fSignalNameGGH ); 
  fBackgroundDYTauTau = TFile::Open( fBackgroundNameDYTauTau ); 
  fBackgroundDYEleEle = TFile::Open( fBackgroundNameDYEleEle ); 
  fBackgroundWJets    = TFile::Open( fBackgroundNameWJets ); 
  fBackgroundQCD      = TFile::Open( fBackgroundNameQCD );
  fBackgroundG1J      = TFile::Open( fBackgroundNameG1J ); 
  fBackgroundTTbar    = TFile::Open( fBackgroundNameTTbar ); 
  fBackgroundSingleTop= TFile::Open( fBackgroundNameSingleTop ); 
  fBackgroundDiBoson  = TFile::Open( fBackgroundNameDiBoson ); 

  if(!fSignalVBF || !fBackgroundDYTauTau || !fBackgroundWJets || !fBackgroundQCD || !fBackgroundG1J || !fBackgroundTTbar || !fBackgroundSingleTop ||
     !fSignalGGH || !fBackgroundDYEleEle || !fBackgroundDiBoson || !fData){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";

  TTree *data                = (TTree*)fData->Get(tree);
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);

  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDYTauTau->Get(tree))->CopyTree("isTauLegMatched>0.5");
  TTree *backgroundDYEleEle    = ((TTree*)fBackgroundDYEleEle->Get(tree))->CopyTree("isTauLegMatched<0.5");
  cout <<backgroundDYTauTau->GetEntries() << " -- " << backgroundDYEleEle->GetEntries() << endl;

  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundQCD       = (TTree*)fBackgroundQCD->Get(tree);
  TTree *backgroundG1J       = (TTree*)fBackgroundG1J->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundSingleTop = (TTree*)fBackgroundSingleTop->Get(tree);
  TTree *backgroundDiBoson   = (TTree*)fBackgroundDiBoson->Get(tree);

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]=data;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["ZfakeTau"]=backgroundDYEleEle;
  tMap["Wjets"]=backgroundWJets;
  tMap["QCD"]=backgroundQCD;
  tMap["G1J"]=backgroundG1J;
  tMap["TTbar"]=backgroundTTbar;
  tMap["SingleTop"]=backgroundSingleTop;
  tMap["DiBoson"]=backgroundDiBoson;

  std::map<std::string,TTree*>::iterator jt;

  
  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  samples.push_back("DiBoson");
  samples.push_back("SingleTop");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("G1J");
  samples.push_back("QCD");
  samples.push_back("ZfakeTau");
  samples.push_back("Ztautau");
  samples.push_back("Data");

  std::map<std::string,float> crossSec;
  crossSec["ggH115"]=( 7.65e-02 * 18.13 );
  crossSec["qqH115"]=( 0.1012);
  crossSec["DiBoson"]=( -1  );
  crossSec["TTbar"]=( 157.5 );
  crossSec["SingleTop"]=( -1 );
  crossSec["Wjets"]=( 31314.0);
  crossSec["ZfakeTau"]=( 3048  );
  crossSec["Ztautau"]=( 3048  );
  crossSec["G1J"]=( -1 );
  crossSec["QCD"]=( 296600000*0.0002855 );
  crossSec["Data"]=( 0 );

  float Lumi = 185.37;


  // here I choose the order in the stack
  std::vector<string> filters;
  filters.push_back("total");
  filters.push_back("vertex");
  filters.push_back("at least 1 e-tau");
  filters.push_back("e pt-eta");
  filters.push_back("e-ID");
  filters.push_back("tau pt-eta");
  filters.push_back("tau-ID");
  filters.push_back("tau against-mu");
  filters.push_back("tau against-e");
  filters.push_back("tau-iso");
  filters.push_back("e-iso");
  filters.push_back("di-e veto");
  filters.push_back("Mt");
  filters.push_back("OS");
  filters.push_back("2-jets");
  filters.push_back("2-tag jets");
  filters.push_back("VBF cuts");
  filters.push_back("jet-veto");
  filters.push_back("anti-b tag");
  filters.push_back("HLT");

  // here I define the map between a sample name and its file ptr
  std::map<std::string,TFile*> fullMap;
  fullMap["Data"]     = fFullData;
  fullMap["ggH115"]   = fFullSignalGGH;
  fullMap["qqH115"]   = fFullSignalVBF;
  fullMap["Ztautau"]  = fFullBackgroundDYTauTau;
  fullMap["ZfakeTau"] = fFullBackgroundDYEleEle;
  fullMap["Wjets"]    = fFullBackgroundWJets;
  fullMap["G1J"]      = fFullBackgroundG1J;
  fullMap["QCD"]      = fFullBackgroundQCD;
  fullMap["TTbar"]    = fFullBackgroundTTbar;
  fullMap["SingleTop"]= fFullBackgroundSingleTop;
  fullMap["DiBoson"]  = fFullBackgroundDiBoson;

  std::map<std::string,TFile*>::iterator it;

  std::map<std::string,float> cutMap_allEventsFilter;
  std::map<std::string,float> cutMap_allEventsFilterE;

  std::map<std::string,float> cutMap_vertexScrapingFilter;
  std::map<std::string,float> cutMap_vertexScrapingFilterE;

  std::map<std::string,float> cutMap_atLeastOneElecTauFilter;
  std::map<std::string,float> cutMap_atLeastOneElecTauFilterE;

  std::map<std::string,float> cutMap_elecPtEtaFilter;
  std::map<std::string,float> cutMap_elecPtEtaFilterE;

  std::map<std::string,float> cutMap_elecPtEtaIDFilter;
  std::map<std::string,float> cutMap_elecPtEtaIDFilterE;

  std::map<std::string,float> cutMap_tauPtEtaFilter;
  std::map<std::string,float> cutMap_tauPtEtaFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDAgMuFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDAgMuFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDAgMuAgElecFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDAgMuAgElecFilterE;

  std::map<std::string,float> cutMap_TauIso;
  std::map<std::string,float> cutMap_TauIsoE;

  std::map<std::string,float> cutMap_ElecIso;
  std::map<std::string,float> cutMap_ElecIsoE;

  std::map<std::string,float> cutMap_DiElecVeto;
  std::map<std::string,float> cutMap_DiElecVetoE;

  std::map<std::string,float> cutMap_Mt;
  std::map<std::string,float> cutMap_MtE;

  std::map<std::string,float> cutMap_OS;
  std::map<std::string,float> cutMap_OSE;

  std::map<std::string,float> cutMap_2jets;
  std::map<std::string,float> cutMap_2jetsE;

  std::map<std::string,float> cutMap_VBFPre;
  std::map<std::string,float> cutMap_VBFPreE;

  std::map<std::string,float> cutMap_VBF;
  std::map<std::string,float> cutMap_VBFE;

  std::map<std::string,float> cutMap_JetVeto;
  std::map<std::string,float> cutMap_JetVetoE;

  std::map<std::string,float> cutMap_antiBtag;
  std::map<std::string,float> cutMap_antiBtagE;

  std::map<std::string,float> cutMap_HLT;
  std::map<std::string,float> cutMap_HLTE;

  std::vector< std::map<std::string,float> > allFilters;
  allFilters.push_back(cutMap_allEventsFilter);
  allFilters.push_back(cutMap_vertexScrapingFilter);
  allFilters.push_back(cutMap_atLeastOneElecTauFilter);
  allFilters.push_back(cutMap_elecPtEtaFilter);
  allFilters.push_back(cutMap_elecPtEtaIDFilter);
  allFilters.push_back(cutMap_tauPtEtaFilter);
  allFilters.push_back(cutMap_tauPtEtaIDFilter);
  allFilters.push_back(cutMap_tauPtEtaIDAgMuFilter);
  allFilters.push_back(cutMap_tauPtEtaIDAgMuAgElecFilter);
 

  std::vector< std::map<std::string,float> > allFiltersE;
  allFiltersE.push_back(cutMap_allEventsFilterE);
  allFiltersE.push_back(cutMap_vertexScrapingFilterE);
  allFiltersE.push_back(cutMap_atLeastOneElecTauFilterE);
  allFiltersE.push_back(cutMap_elecPtEtaFilterE);
  allFiltersE.push_back(cutMap_elecPtEtaIDFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDAgMuFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDAgMuAgElecFilterE);
 

  std::vector< string > filtersByName;
  filtersByName.push_back("allEventsFilter");
  filtersByName.push_back("vertexScrapingFilter");
  filtersByName.push_back("atLeastOneElecTauFilter");
  filtersByName.push_back("elecPtEtaFilter");
  filtersByName.push_back("elecPtEtaIDFilter");
  filtersByName.push_back("tauPtEtaFilter");
  filtersByName.push_back("tauPtEtaIDFilter");
  filtersByName.push_back("tauPtEtaIDAgMuFilter");
  filtersByName.push_back("tauPtEtaIDAgMuAgElecFilter");
  
  for(unsigned int k = 0 ; k <filtersByName.size() ; k++){

    for(it = fullMap.begin(); it != fullMap.end(); it++){
      TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
      float totalEvents = allEvents->GetBinContent(1);
      allEvents = (TH1F*)(it->second)->Get( (filtersByName[k]+"/totalEvents").c_str() );
      float tot =  allEvents->GetBinContent(1);
      float totalEquivalentEvents = allEvents->GetEffectiveEntries();
      if(crossSec[it->first]>0){
	tot *= (Lumi/ (totalEvents/crossSec[it->first]) * hltEff_  ) ;
      } else if (crossSec[it->first]<0) tot *= Lumi/1000*hltEff_;
      (allFilters[k])[it->first]  = tot;
      (allFiltersE[k])[it->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    }
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_TauIso[jt->first] = tot;
    cutMap_TauIsoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_ElecIso[jt->first] = tot;
    cutMap_ElecIsoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_DiElecVeto[jt->first] = tot;
    cutMap_DiElecVetoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_Mt[jt->first] = tot;
    cutMap_MtE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_OS[jt->first] = tot;
    cutMap_OSE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>20 && pt2>15)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_2jets[jt->first] = tot;
    cutMap_2jetsE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1  && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>20 && pt2>15 && eta1*eta2<0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_VBFPre[jt->first] = tot;
    cutMap_VBFPreE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_VBF[jt->first] = tot;
    cutMap_VBFE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_JetVeto[jt->first] = tot;
    cutMap_JetVetoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20 && jetsBtagHE1<2.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_antiBtag[jt->first] = tot;
    cutMap_antiBtagE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (jt->first).find("Data")!=string::npos ? "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20  && jetsBtagHE1<2.1 && HLTmatch==1 && HLTx==1)" : "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1  && elecFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20 && jetsBtagHE1<2.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_HLT[jt->first] = tot;
    cutMap_HLTE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }
  
  allFilters.push_back(cutMap_TauIso);
  allFilters.push_back(cutMap_ElecIso);
  allFilters.push_back(cutMap_DiElecVeto);
  allFilters.push_back(cutMap_Mt);
  allFilters.push_back(cutMap_OS);
  allFilters.push_back(cutMap_2jets);
  allFilters.push_back(cutMap_VBFPre);
  allFilters.push_back(cutMap_VBF);
  allFilters.push_back(cutMap_JetVeto);
  allFilters.push_back(cutMap_antiBtag);
  allFilters.push_back(cutMap_HLT);

  allFiltersE.push_back(cutMap_TauIsoE);
  allFiltersE.push_back(cutMap_ElecIsoE);
  allFiltersE.push_back(cutMap_DiElecVetoE);
  allFiltersE.push_back(cutMap_MtE);
  allFiltersE.push_back(cutMap_OSE);
  allFiltersE.push_back(cutMap_2jetsE);
  allFiltersE.push_back(cutMap_VBFPreE);
  allFiltersE.push_back(cutMap_VBFE);
  allFiltersE.push_back(cutMap_JetVetoE);
  allFiltersE.push_back(cutMap_antiBtagE);
  allFiltersE.push_back(cutMap_HLTE);

  std::vector<std::pair<float,float> > SMsums;

  for(unsigned int i = 0; i < allFilters.size(); i++){
    float SM=0,SME2=0;
    int helper = 0;
    for(int k = 0 ; k < samples.size()-1; k++){
      if( !(samples[k].find("Z")!=string::npos && i<9 && helper!=0) ) 
	SM+=(allFilters[i].find(samples[k]))->second;
      if( !(samples[k].find("Z")!=string::npos && i<9 && helper!=0) ) 
	SME2+=(allFiltersE[i].find(samples[k]))->second*(allFiltersE[i].find(samples[k]))->second;
      if( samples[k].find("Z")!=string::npos ) helper++;
    }
    SMsums.push_back( make_pair(SM,SME2 ));
  }


  //out<<"\\begin{center}"<<endl;
  out<<"\\begin{tabular}[!htbp]{|c";
  for(int k = 0 ; k < samples.size(); k++) out<<"|c";
  out<<"|c|} \\hline"<<endl;
  out<< "Cut & ";
  for(int k = 0 ; k < samples.size(); k++){
    out << (fullMap.find(samples[k]))->first;
    if(k!=samples.size()-1) out <<" & " ;
    else out << " & $\\Sigma$ SM \\\\ " << endl;
  }
  out <<  " \\hline" << endl;

  
  for(int i = 0; i < allFilters.size(); i++){
    out << filters[i] << " & ";
    for(int k = 0 ; k < samples.size(); k++){
      if(samples[k].find("Data")==string::npos) 
	out << (allFilters[i].find(samples[k]))->second << " $\\pm$ " << (allFiltersE[i].find(samples[k]))->second;
      else
	out << (allFilters[i].find(samples[k]))->second;
      if(k!=samples.size()-1) out <<" & " ;
      else out << " & " << SMsums[i].first << "$\\pm$" <<  sqrt(SMsums[i].second) << " \\\\ " << endl;
    }
    out <<  " \\hline" << endl;

  }
  
  out<<"\\end{tabular}"<<endl;
  //out<<"\\end{center}"<<endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *c1 = new TCanvas("c1CutFlowMass","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);
 
  TLegend* leg = new TLegend(0.55,0.55,0.80,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader( "e+#tau" );

  THStack* aStack = new THStack("aStack",Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}",Lumi));

  std::vector<TH1F*> histos;
  for(unsigned int k = 0 ; k < samples.size(); k++){
    TH1F* h1 = new TH1F( ("h1_"+samples[k]).c_str(), Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1} ; ; Events",Lumi) , filters.size(), 0, filters.size());

    if( (samples[k]).find("ZfakeTau")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z+jets, fake-#tau","F");
    }
    if( (samples[k]).find("Ztautau")!=string::npos ) {
      h1->SetFillColor(kRed);
      leg->AddEntry(h1,"Z+jets, genuine #tau","F");
    }
    if( (samples[k]).find("TTbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"t#bar{t}+jets","F");
    }
    if( (samples[k]).find("SingleTop")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"single-t","F");
    }
    if( (samples[k]).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"W+jets","F");
    }
    if( (samples[k]).find("DiBoson")!=string::npos ) {
      h1->SetFillColor(38);
      leg->AddEntry(h1,"Di-Boson","F");
    }
    if( (samples[k]).find("QCD")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"QCD-multijets","F");
    }
    if( (samples[k]).find("G1J")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"#gamma+jet","F");
    }
    if((samples[k]).find("qqH115")!=string::npos){
      h1->SetLineWidth(2);
      h1->SetLineColor(kBlack);
      h1->SetFillColor(kBlack);     
      h1->SetFillStyle(3004);
      leg->AddEntry(h1,"qqH(115)","F");
    }
    if((samples[k]).find("ggH115")!=string::npos){
      h1->SetLineWidth(2);
      h1->SetLineColor(kBlack);  
      h1->SetFillColor(kBlack);     
      h1->SetFillStyle(3005);
      leg->AddEntry(h1,"ggH(115)","F");
    }
    

    //h1->SetLineWidth(1.4);
    for(int v = 0; v < filters.size(); v++){
      h1->GetXaxis()->SetBinLabel(v+1,filters[v].c_str());
      h1->SetBinContent(v+1, (allFilters[v].find(samples[k]))->second );
    }
    if(samples[k].find("Data")!=string::npos){
      h1->Sumw2();
      h1->SetMarkerStyle(20);
      h1->SetMarkerSize(1.2);
      h1->SetMarkerColor(kBlack);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,"DATA","P");
    }
    //else leg->AddEntry(h1,samples[k].c_str(),"L");
    histos.push_back(h1);
  }

  for(unsigned int k = 0 ; k < histos.size(); k++){
    string histoName =  std::string( histos[k]->GetName() );
    if( histoName.find("Data")== string::npos) aStack->Add(histos[k]);
  }
  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  aStack->GetYaxis()->SetRangeUser(0.01,10000000);

  for(unsigned int k = 0 ; k < histos.size(); k++){
    string histoName =  std::string( histos[k]->GetName() );
    if( histoName.find("Data")!= string::npos){
      histos[k]->Sumw2();
      histos[k]->Draw("PSAME");
    }
  }

  leg->Draw();

  return;

}

///////////////////////////////////////////////////////////////////////////////////////////////////

void cutFlowStudyMu( Double_t hltEff_=0.9218 ) 
{
  
  ofstream out("cutFlow_MuTauStream_iter2.txt");
  out.precision(4);

 
  TFile *fFullData                = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_Run2011-Mu.root","READ"); 
  TFile *fFullSignalVBF           = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_VBFH115-Mu-powheg-PUS1.root","READ"); 
  TFile *fFullSignalGGH           = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_GGFH115-Mu-powheg-PUS1.root","READ");  
  TFile *fFullBackgroundDYTauTau  = new TFile("/data_CMS/cms/lbianchini/MuTauStream2011_iter2/treeMuTauStream_DYJets-Mu-50-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundDYMuMu    = new TFile("/data_CMS/cms/lbianchini/MuTauStream2011_iter2/treeMuTauStream_DYJets-Mu-50-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundWJets     = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_WJets-Mu-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundQCD       = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_QCDmu.root","READ"); 
  TFile *fFullBackgroundTTbar     = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_TTJets-Mu-madgraph-PUS1.root","READ"); 
  TFile *fFullBackgroundSingleTop = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_SingleTop-Mu.root","READ"); 
  TFile *fFullBackgroundDiBoson   = new TFile("/data_CMS/cms/lbianchini//MuTauStream2011_iter2/treeMuTauStream_DiBoson-Mu.root","READ"); 

  // OpenNTuples
  TString fDataName                = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleRun2011-Mu_Open_MuTauStream.root";
  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleVBFH115-Mu-powheg-PUS1_Open_MuTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleGGFH115-Mu-powheg-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleDYJets-Mu-50-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameDYMuMu    = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleDYJets-Mu-50-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleWJets-Mu-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameQCD       = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleQCDmu_Open_MuTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleTTJets-Mu-madgraph-PUS1_Open_MuTauStream.root";
  TString fBackgroundNameSingleTop = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleSingleTop-Mu_Open_MuTauStream.root";
  TString fBackgroundNameDiBoson   = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleDiBoson-Mu_Open_MuTauStream.root";

  TFile *fData(0); 
  TFile *fSignalVBF(0); 
  TFile *fSignalGGH(0); 
  TFile *fBackgroundDYTauTau(0);
  TFile *fBackgroundDYMuMu(0);
  TFile *fBackgroundWJets(0);
  TFile *fBackgroundQCD(0);
  TFile *fBackgroundTTbar(0);
  TFile *fBackgroundSingleTop(0);
  TFile *fBackgroundDiBoson(0);
 
  fData               = TFile::Open( fDataName ); 
  fSignalVBF          = TFile::Open( fSignalNameVBF ); 
  fSignalGGH          = TFile::Open( fSignalNameGGH ); 
  fBackgroundDYTauTau = TFile::Open( fBackgroundNameDYTauTau ); 
  fBackgroundDYMuMu   = TFile::Open( fBackgroundNameDYMuMu ); 
  fBackgroundWJets    = TFile::Open( fBackgroundNameWJets ); 
  fBackgroundQCD      = TFile::Open( fBackgroundNameQCD ); 
  fBackgroundTTbar    = TFile::Open( fBackgroundNameTTbar ); 
  fBackgroundSingleTop= TFile::Open( fBackgroundNameSingleTop ); 
  fBackgroundDiBoson  = TFile::Open( fBackgroundNameDiBoson ); 

  if(!fSignalVBF || !fBackgroundDYTauTau || !fBackgroundWJets || !fBackgroundQCD || !fBackgroundTTbar || !fBackgroundSingleTop ||
     !fSignalGGH || !fBackgroundDYMuMu || !fBackgroundDiBoson || !fData){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";

  TTree *data                = (TTree*)fData->Get(tree);
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);

  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDYTauTau->Get(tree))->CopyTree("isTauLegMatched>0.5");
  TTree *backgroundDYMuMu    = ((TTree*)fBackgroundDYMuMu->Get(tree))->CopyTree("isTauLegMatched<0.5");
  cout <<backgroundDYTauTau->GetEntries() << " -- " << backgroundDYMuMu->GetEntries() << endl;

  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundQCD       = (TTree*)fBackgroundQCD->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundSingleTop = (TTree*)fBackgroundSingleTop->Get(tree);
  TTree *backgroundDiBoson   = (TTree*)fBackgroundDiBoson->Get(tree);

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]=data;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["ZfakeTau"]=backgroundDYMuMu;
  tMap["Wjets"]=backgroundWJets;
  tMap["QCD"]=backgroundQCD;
  tMap["TTbar"]=backgroundTTbar;
  tMap["SingleTop"]=backgroundSingleTop;
  tMap["DiBoson"]=backgroundDiBoson;

  std::map<std::string,TTree*>::iterator jt;

  
  // here I choose the order in the stack
  std::vector<string> samples;
  samples.push_back("ggH115");
  samples.push_back("qqH115");
  samples.push_back("DiBoson");
  samples.push_back("SingleTop");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("QCD");
  samples.push_back("ZfakeTau");
  samples.push_back("Ztautau");
  samples.push_back("Data");

  std::map<std::string,float> crossSec;
  crossSec["ggH115"]=( 7.65e-02 * 18.13 );
  crossSec["qqH115"]=( 0.1012);
  crossSec["DiBoson"]=( -1  );
  crossSec["TTbar"]=( 157.5 );
  crossSec["SingleTop"]=( -1 );
  crossSec["Wjets"]=( 31314.0);
  crossSec["ZfakeTau"]=( 3048  );
  crossSec["Ztautau"]=( 3048  );
  crossSec["QCD"]=( 296600000*0.0002855 );
  crossSec["Data"]=( 0 );

  float Lumi = 24.86+159.15;


  // here I choose the order in the stack
  std::vector<string> filters;
  filters.push_back("total");
  filters.push_back("vertex");
  filters.push_back("at least 1 mu-tau");
  filters.push_back("mu pt-eta");
  filters.push_back("mu-ID");
  filters.push_back("tau pt-eta");
  filters.push_back("tau-ID");
  filters.push_back("tau against-mu");
  filters.push_back("tau against-e");
  filters.push_back("tau-iso");
  filters.push_back("mu-iso");
  filters.push_back("di-mu veto");
  filters.push_back("Mt");
  filters.push_back("OS");
  filters.push_back("2-jets");
  filters.push_back("2-tag jets");
  filters.push_back("VBF cuts");
  filters.push_back("jet-veto");
  filters.push_back("anti-b tag");
  filters.push_back("HLT");

  // here I define the map between a sample name and its file ptr
  std::map<std::string,TFile*> fullMap;
  fullMap["Data"]     = fFullData;
  fullMap["ggH115"]   = fFullSignalGGH;
  fullMap["qqH115"]   = fFullSignalVBF;
  fullMap["Ztautau"]  = fFullBackgroundDYTauTau;
  fullMap["ZfakeTau"] = fFullBackgroundDYMuMu;
  fullMap["Wjets"]    = fFullBackgroundWJets;
  fullMap["QCD"]      = fFullBackgroundQCD;
  fullMap["TTbar"]    = fFullBackgroundTTbar;
  fullMap["SingleTop"]= fFullBackgroundSingleTop;
  fullMap["DiBoson"]  = fFullBackgroundDiBoson;

  std::map<std::string,TFile*>::iterator it;

  std::map<std::string,float> cutMap_allEventsFilter;
  std::map<std::string,float> cutMap_allEventsFilterE;

  std::map<std::string,float> cutMap_vertexScrapingFilter;
  std::map<std::string,float> cutMap_vertexScrapingFilterE;

  std::map<std::string,float> cutMap_atLeastOneMuTauFilter;
  std::map<std::string,float> cutMap_atLeastOneMuTauFilterE;

  std::map<std::string,float> cutMap_muPtEtaFilter;
  std::map<std::string,float> cutMap_muPtEtaFilterE;

  std::map<std::string,float> cutMap_muPtEtaIDFilter;
  std::map<std::string,float> cutMap_muPtEtaIDFilterE;

  std::map<std::string,float> cutMap_tauPtEtaFilter;
  std::map<std::string,float> cutMap_tauPtEtaFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDAgMuFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDAgMuFilterE;

  std::map<std::string,float> cutMap_tauPtEtaIDAgMuAgElecFilter;
  std::map<std::string,float> cutMap_tauPtEtaIDAgMuAgElecFilterE;

  std::map<std::string,float> cutMap_TauIso;
  std::map<std::string,float> cutMap_TauIsoE;

  std::map<std::string,float> cutMap_MuIso;
  std::map<std::string,float> cutMap_MuIsoE;

  std::map<std::string,float> cutMap_DiMuVeto;
  std::map<std::string,float> cutMap_DiMuVetoE;

  std::map<std::string,float> cutMap_Mt;
  std::map<std::string,float> cutMap_MtE;

  std::map<std::string,float> cutMap_OS;
  std::map<std::string,float> cutMap_OSE;

  std::map<std::string,float> cutMap_2jets;
  std::map<std::string,float> cutMap_2jetsE;

  std::map<std::string,float> cutMap_VBFPre;
  std::map<std::string,float> cutMap_VBFPreE;

  std::map<std::string,float> cutMap_VBF;
  std::map<std::string,float> cutMap_VBFE;

  std::map<std::string,float> cutMap_JetVeto;
  std::map<std::string,float> cutMap_JetVetoE;

  std::map<std::string,float> cutMap_antiBtag;
  std::map<std::string,float> cutMap_antiBtagE;

  std::map<std::string,float> cutMap_HLT;
  std::map<std::string,float> cutMap_HLTE;

  std::vector< std::map<std::string,float> > allFilters;
  allFilters.push_back(cutMap_allEventsFilter);
  allFilters.push_back(cutMap_vertexScrapingFilter);
  allFilters.push_back(cutMap_atLeastOneMuTauFilter);
  allFilters.push_back(cutMap_muPtEtaFilter);
  allFilters.push_back(cutMap_muPtEtaIDFilter);
  allFilters.push_back(cutMap_tauPtEtaFilter);
  allFilters.push_back(cutMap_tauPtEtaIDFilter);
  allFilters.push_back(cutMap_tauPtEtaIDAgMuFilter);
  allFilters.push_back(cutMap_tauPtEtaIDAgMuAgElecFilter);
 

  std::vector< std::map<std::string,float> > allFiltersE;
  allFiltersE.push_back(cutMap_allEventsFilterE);
  allFiltersE.push_back(cutMap_vertexScrapingFilterE);
  allFiltersE.push_back(cutMap_atLeastOneMuTauFilterE);
  allFiltersE.push_back(cutMap_muPtEtaFilterE);
  allFiltersE.push_back(cutMap_muPtEtaIDFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDAgMuFilterE);
  allFiltersE.push_back(cutMap_tauPtEtaIDAgMuAgElecFilterE);
 

  std::vector< string > filtersByName;
  filtersByName.push_back("allEventsFilter");
  filtersByName.push_back("vertexScrapingFilter");
  filtersByName.push_back("atLeastOneMuTauFilter");
  filtersByName.push_back("muPtEtaFilter");
  filtersByName.push_back("muPtEtaIDFilter");
  filtersByName.push_back("tauPtEtaFilter");
  filtersByName.push_back("tauPtEtaIDFilter");
  filtersByName.push_back("tauPtEtaIDAgMuFilter");
  filtersByName.push_back("tauPtEtaIDAgMuAgElecFilter");
  
  for(unsigned int k = 0 ; k <filtersByName.size() ; k++){

    for(it = fullMap.begin(); it != fullMap.end(); it++){
      TH1F* allEvents = (TH1F*)(it->second)->Get("allEventsFilter/totalEvents");
      float totalEvents = allEvents->GetBinContent(1);
      allEvents = (TH1F*)(it->second)->Get( (filtersByName[k]+"/totalEvents").c_str() );
      float tot =  allEvents->GetBinContent(1);
      float totalEquivalentEvents = allEvents->GetEffectiveEntries();
      if(crossSec[it->first]>0){
	tot *= (Lumi/ (totalEvents/crossSec[it->first]) * hltEff_  ) ;
      } else if (crossSec[it->first]<0) tot *= Lumi/1000*hltEff_;
      (allFilters[k])[it->first]  = tot;
      (allFiltersE[k])[it->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    }
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_TauIso[jt->first] = tot;
    cutMap_TauIsoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_MuIso[jt->first] = tot;
    cutMap_MuIsoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_DiMuVeto[jt->first] = tot;
    cutMap_DiMuVetoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_Mt[jt->first] = tot;
    cutMap_MtE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_OS[jt->first] = tot;
    cutMap_OSE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>20 && pt2>15)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_2jets[jt->first] = tot;
    cutMap_2jetsE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1  && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>20 && pt2>15 && eta1*eta2<0)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_VBFPre[jt->first] = tot;
    cutMap_VBFPreE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_VBF[jt->first] = tot;
    cutMap_VBFE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_JetVeto[jt->first] = tot;
    cutMap_JetVetoE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20 && jetsBtagHE1<2.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_antiBtag[jt->first] = tot;
    cutMap_antiBtagE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }

  for(jt = tMap.begin(); jt != tMap.end(); jt++){
    cout<<jt->first<<endl;
    TH1F* h1 = new TH1F("h1","",1,-10,10); 
    TCut cut =  (jt->first).find("Data")!=string::npos ? "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20  && jetsBtagHE1<2.1 && HLTmatch==1 && ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)))" : "sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1  && muFlag==0 && MtLeg1<40 && diTauCharge==0 && pt1>38 && pt2>25 && eta1*eta2<0 && Deta>2.7 && Mjj>520 && ptVeto<20 && jetsBtagHE1<2.1)"; 
    jt->second->Draw("etaL1>>h1",cut);
    if((jt->first).find("Data")==string::npos) h1->Scale(Lumi/1000*hltEff_); 
    float tot = h1->Integral();
    float totalEquivalentEvents = h1->GetEffectiveEntries();
    cutMap_HLT[jt->first] = tot;
    cutMap_HLTE[jt->first] = totalEquivalentEvents>0 ? sqrt(totalEquivalentEvents)*(tot/totalEquivalentEvents) : 0;
    delete h1;
  }
  
  allFilters.push_back(cutMap_TauIso);
  allFilters.push_back(cutMap_MuIso);
  allFilters.push_back(cutMap_DiMuVeto);
  allFilters.push_back(cutMap_Mt);
  allFilters.push_back(cutMap_OS);
  allFilters.push_back(cutMap_2jets);
  allFilters.push_back(cutMap_VBFPre);
  allFilters.push_back(cutMap_VBF);
  allFilters.push_back(cutMap_JetVeto);
  allFilters.push_back(cutMap_antiBtag);
  allFilters.push_back(cutMap_HLT);

  allFiltersE.push_back(cutMap_TauIsoE);
  allFiltersE.push_back(cutMap_MuIsoE);
  allFiltersE.push_back(cutMap_DiMuVetoE);
  allFiltersE.push_back(cutMap_MtE);
  allFiltersE.push_back(cutMap_OSE);
  allFiltersE.push_back(cutMap_2jetsE);
  allFiltersE.push_back(cutMap_VBFPreE);
  allFiltersE.push_back(cutMap_VBFE);
  allFiltersE.push_back(cutMap_JetVetoE);
  allFiltersE.push_back(cutMap_antiBtagE);
  allFiltersE.push_back(cutMap_HLTE);

  std::vector<std::pair<float,float> > SMsums;

  for(unsigned int i = 0; i < allFilters.size(); i++){
    float SM=0,SME2=0;
    int helper = 0;
    for(int k = 0 ; k < samples.size()-1; k++){
      if( !(samples[k].find("Z")!=string::npos && i<9 && helper!=0) ) 
	SM+=(allFilters[i].find(samples[k]))->second;
      if( !(samples[k].find("Z")!=string::npos && i<9 && helper!=0) ) 
	SME2+=(allFiltersE[i].find(samples[k]))->second*(allFiltersE[i].find(samples[k]))->second;
      if( samples[k].find("Z")!=string::npos ) helper++;
    }
    SMsums.push_back( make_pair(SM,SME2 ));
  }


  //out<<"\\begin{center}"<<endl;
  out<<"\\begin{tabular}[!htbp]{|c";
  for(int k = 0 ; k < samples.size(); k++) out<<"|c";
  out<<"|c|} \\hline"<<endl;
  out<< "Cut & ";
  for(int k = 0 ; k < samples.size(); k++){
    out << (fullMap.find(samples[k]))->first;
    if(k!=samples.size()-1) out <<" & " ;
    else out << " & $\\Sigma$ SM \\\\ " << endl;
  }
  out <<  " \\hline" << endl;

  
  for(int i = 0; i < allFilters.size(); i++){
    out << filters[i] << " & ";
    for(int k = 0 ; k < samples.size(); k++){
      if(samples[k].find("Data")==string::npos) 
	out << (allFilters[i].find(samples[k]))->second << " $\\pm$ " << (allFiltersE[i].find(samples[k]))->second;
      else
	out << (allFilters[i].find(samples[k]))->second;
      if(k!=samples.size()-1) out <<" & " ;
      else out << " & " << SMsums[i].first << "$\\pm$" <<  sqrt(SMsums[i].second) << " \\\\ " << endl;
    }
    out <<  " \\hline" << endl;

  }
  
  out<<"\\end{tabular}"<<endl;
  //out<<"\\end{center}"<<endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *c1 = new TCanvas("c1CutFlowMass","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);
 
  TLegend* leg = new TLegend(0.55,0.55,0.80,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader( "#mu+#tau" );

  THStack* aStack = new THStack("aStack",Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}",Lumi));

  std::vector<TH1F*> histos;
  for(unsigned int k = 0 ; k < samples.size(); k++){
    TH1F* h1 = new TH1F( ("h1_"+samples[k]).c_str(), Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1} ; ; Events",Lumi) , filters.size(), 0, filters.size());

    if( (samples[k]).find("ZfakeTau")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z+jets, fake-#tau","F");
    }
    if( (samples[k]).find("Ztautau")!=string::npos ) {
      h1->SetFillColor(kRed);
      leg->AddEntry(h1,"Z+jets, genuine #tau","F");
    }
    if( (samples[k]).find("TTbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"t#bar{t}+jets","F");
    }
    if( (samples[k]).find("SingleTop")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"single-t","F");
    }
    if( (samples[k]).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"W+jets","F");
    }
    if( (samples[k]).find("DiBoson")!=string::npos ) {
      h1->SetFillColor(38);
      leg->AddEntry(h1,"Di-Boson","F");
    }
    if( (samples[k]).find("QCD")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"QCD-multijets","F");
    }
    if((samples[k]).find("qqH115")!=string::npos){
      h1->SetLineWidth(2);
      h1->SetLineColor(kBlack);
      h1->SetFillColor(kBlack);     
      h1->SetFillStyle(3004);
      leg->AddEntry(h1,"qqH(115)","F");
    }
    if((samples[k]).find("ggH115")!=string::npos){
      h1->SetLineWidth(2);
      h1->SetLineColor(kBlack);  
      h1->SetFillColor(kBlack);     
      h1->SetFillStyle(3005);
      leg->AddEntry(h1,"ggH(115)","F");
    }
    

    //h1->SetLineWidth(1.4);
    for(int v = 0; v < filters.size(); v++){
      h1->GetXaxis()->SetBinLabel(v+1,filters[v].c_str());
      h1->SetBinContent(v+1, (allFilters[v].find(samples[k]))->second );
    }
    if(samples[k].find("Data")!=string::npos){
      h1->Sumw2();
      h1->SetMarkerStyle(20);
      h1->SetMarkerSize(1.2);
      h1->SetMarkerColor(kBlack);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,"DATA","P");
    }
    //else leg->AddEntry(h1,samples[k].c_str(),"L");
    histos.push_back(h1);
  }

  for(unsigned int k = 0 ; k < histos.size(); k++){
    string histoName =  std::string( histos[k]->GetName() );
    if( histoName.find("Data")== string::npos) aStack->Add(histos[k]);
  }
  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  aStack->GetYaxis()->SetRangeUser(0.01,10000000);

  for(unsigned int k = 0 ; k < histos.size(); k++){
    string histoName =  std::string( histos[k]->GetName() );
    if( histoName.find("Data")!= string::npos){
      histos[k]->Sumw2();
      histos[k]->Draw("PSAME");
    }
  }

  leg->Draw();

  return;

}
