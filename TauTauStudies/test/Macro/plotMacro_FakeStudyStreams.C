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

void fakeStudyElec( TString weightFile = "TMVAClassificationPtOrd_qqH115vsWZttQCD_Cuts.weights.xml",
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
  leg->SetHeader("e+#tau_{had}");

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

      if(reader->EvaluateMVA( "Cuts", effS_ )){
	counter++;
	h1->Fill( vMap[variable_], sampleWeight);
      }
    
    }// end loop


    if( (it->first).find("Zeleele")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z#rightarrow ee","F");
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
    if( (it->first).find("G1J")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"#gamma+jet","F");
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

      if(reader->EvaluateMVA( "Cuts", effS_ )){
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



