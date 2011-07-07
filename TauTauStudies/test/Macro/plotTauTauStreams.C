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
#include "TCut.h"

#define VERBOSE true
#define SAVE true


////////////////////////////////////////////////////////

void plotElecTau( TCut Cuts_ = "combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && diTauCharge==0 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20",
		  Bool_t doVBFPreRegion_ = false,
		  TString variable_ = "diTauVisMass",
		  TString XTitle_ = "mass",
		  TString Unities_ = "GeV",
		  Int_t nBins_ = 40, Float_t xMin_=0, Float_t xMax_=200,
		  Float_t magnifySgn_ = 1,
		  Float_t hltEff_ = 1.1,
		  Int_t enableHLTmatching_ = 0,
		  Int_t logy_ = 0 ) 
{   
  
  
  
  float Lumi = 123.;
  
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);

  TLegend* leg = new TLegend(0.55,0.47,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader("e+#tau_{had}");

  THStack* aStack = new THStack("aStack",Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}", Lumi ));
  TH1F* hSiml  = new TH1F();
  TH1F* hSgn   = new TH1F();
  TH1F* hSgn1  = new TH1F();
  TH1F* hSgn2  = new TH1F();
  TH1F* hData  = new TH1F();

  // OpenNTuples
  TString fDataName                = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleRun2011-Elec_Open_ElecTauStream.root";
  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleVBFH125-powheg-PUS4_run_Open_ElecTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleGGFH125-powheg-PUS4_run_Open_ElecTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleDYJets-50-madgraph-PUS4_run_Open_ElecTauStream.root";
  TString fBackgroundNameDYElecElec    = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleDYJets-50-madgraph-PUS4_run_Open_ElecTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleWJets-madgraph-PUS4_run_Open_ElecTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamSummer11_iter1//nTupleTT-pythia-PUS3_run_Open_ElecTauStream.root";

  TFile *fData              = new TFile( fDataName, "READ");  
  TFile *fSignalVBF         = new TFile(fSignalNameVBF ,"READ");  
  TFile *fSignalGGH         = new TFile(fSignalNameGGH,"READ");  
  TFile *fBackgroundDYTauTau= new TFile(fBackgroundNameDYTauTau,"READ"); 
  TFile *fBackgroundDYElecElec  = new TFile(fBackgroundNameDYElecElec,"READ"); 
  TFile *fBackgroundWJets   = new TFile(fBackgroundNameWJets,"READ"); 
  TFile *fBackgroundTTbar   = new TFile(fBackgroundNameTTbar,"READ"); 

  if(!fSignalVBF || 
     !fBackgroundDYTauTau || 
     !fBackgroundWJets || 
     !fBackgroundTTbar ||
     !fSignalGGH || 
     !fBackgroundDYElecElec ||
     !fData ){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";
  TCut inclusiveRegionSS("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && HLTx==1)");
  TCut vbfPreRegionSS("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1)");
  TCut inclusiveRegionSSMtSide("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)<-40 && diTauCharge!=0 && HLTx==1)");
  TCut vbfPreRegionSSMtSide("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)<-40 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1)");
  TCut inclusiveRegionSSmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && HLTx==1 )");
  TCut vbfPreRegionSSmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1 )");
  TCut inclusiveRegionSSMtmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && diTauCharge!=0 && HLTx==1 )");
  TCut vbfPreRegionSSMtmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && elecFlag==0 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1 )");

  TTree *data                = (TTree*)fData->Get(tree);
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);

  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDYTauTau->Get(tree))->CopyTree("abs(genDecay)==(23*15)");
  TTree *backgroundDYElecElec    = ((TTree*)fBackgroundDYElecElec->Get(tree))->CopyTree("abs(genDecay)==(23*11)");
  cout <<backgroundDYTauTau->GetEntries() << " -- " << backgroundDYElecElec->GetEntries() << endl;

  TTree* dataSS              = 0;
  TTree* WJetsSS             = 0;
  TTree* WJetsSSFullMt       = 0;
  dataSS  = doVBFPreRegion_ ? 
    ((TTree*)fData->Get(tree))->CopyTree(vbfPreRegionSS) : 
    ((TTree*)fData->Get(tree))->CopyTree(inclusiveRegionSS);
  WJetsSS = doVBFPreRegion_ ? 
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(vbfPreRegionSSmc) :
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(inclusiveRegionSSmc); 
  WJetsSSFullMt = doVBFPreRegion_ ? 
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(vbfPreRegionSSMtmc) :
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(inclusiveRegionSSMtmc);
  
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);

  // here I choose the order in the stack
  std::vector<string> samples;
  //samples.push_back("ggH115");
  //samples.push_back("qqH115");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("SS");
  samples.push_back("Zelecelec");
  samples.push_back("Ztautau");
  samples.push_back("Data");
  

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]=data;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["Zelecelec"]=backgroundDYElecElec;
  tMap["Wjets"]=backgroundWJets;
  tMap["TTbar"]=backgroundTTbar;
  tMap["SS"]=dataSS;


  Float_t pt1, pt2, Deta, Mjj, eta1, eta2;
  Float_t diTauNSVfitMass,diTauNSVfitMassErrUp,diTauNSVfitMassErrDown,diTauVisPt,diTauVisEta,diTauVisMass,ptL1,ptL2,etaL1,etaL2,diTauCharge,MtLeg1,numPV,combRelIsoLeg1DBeta,combRelIsoLeg1,pZeta,pZetaVis, MEt, pZetaCorr,pZetaVisCorr, MEtCorr, MtLeg1Corr, ptVeto, leadTrackPt, jetsBtagHE1, jetsBtagHE2;
  Float_t sampleWeight,puWeight, HLTweight;

  std::map<TString,Float_t> vMap;


  for( unsigned iter=0; iter<samples.size(); iter++){

    cout << samples[iter] << endl;
    
    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}", Lumi) , nBins_ ,xMin_ , xMax_);
    
    TFile* dummy = new TFile("dummy.root","RECREATE");  
    TCut Cuts = Cuts_;
    TCut hltData   = "(HLTx==1)";   
    TCut hltMC     = "(HLTx==1)";
    if((it->first).find("Data")!=string::npos) Cuts = Cuts_ && hltData;
    else Cuts = Cuts_ && hltMC;
    
    TTree* currentTree = 0;
    
    if((it->first).find("SS")!=string::npos){
      
      cout << "Dealing with SS" << endl;
      
      currentTree = (it->second);
      
      //cout << currentTree->GetEntries() << endl;
      TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
      hHelp->Reset();
      TString variable = ( (std::string(variable_.Data())).find("pZeta")!=string::npos ) ? "pZetaCorr-1.5*pZetaVisCorr" : variable_;
      currentTree->Draw(variable+">>hHelp");
      //cout << hHelp->GetEntries() << endl;
      h1->Add(hHelp,1);
      //cout << h1->Integral() << endl;
      TH1F* hWMt = new TH1F("hWMt","",200,-200,200);
      WJetsSSFullMt->Draw("(pZetaCorr-1.5*pZetaVisCorr)>>hWMt");
      float scaleFactor = (hWMt->Integral(0,80))/(hWMt->Integral(90,200));
      cout << "Scale factor = " << scaleFactor << endl;
      TH1F* hW = (TH1F*)h1->Clone("hW");
      hW->Reset();
      float SSWinSignalRegion = doVBFPreRegion_ ?
	data->GetEntries(inclusiveRegionSSMtSide) :  
	data->GetEntries(vbfPreRegionSSMtSide);
      SSWinSignalRegion *= (1./scaleFactor);
      WJetsSS->Draw(variable+">>hW");
      hW->Scale(SSWinSignalRegion/hW->Integral());
      h1->Add(hW, -1 );
      h1->Scale(1.06);
      //cout << h1->Integral() << endl;
    }
    else{

      currentTree =  (TTree*)(it->second)->CopyTree(Cuts);
      
      Int_t counter = 0;
      
      currentTree->SetBranchAddress( "pt1", &pt1 );
      currentTree->SetBranchAddress( "pt2", &pt2 );
      currentTree->SetBranchAddress( "Deta",&Deta );
      currentTree->SetBranchAddress( "Mjj", &Mjj );
      currentTree->SetBranchAddress( "diTauNSVfitMass",&diTauNSVfitMass);
      currentTree->SetBranchAddress( "diTauNSVfitMassErrUp",&diTauNSVfitMassErrUp);
      currentTree->SetBranchAddress( "diTauNSVfitMassErrDown",&diTauNSVfitMassErrDown);

      currentTree->SetBranchAddress( "diTauVisPt",&diTauVisPt);
      currentTree->SetBranchAddress( "diTauVisEta",&diTauVisEta);
      currentTree->SetBranchAddress( "diTauVisMass",&diTauVisMass);
      currentTree->SetBranchAddress( "ptL1",  &ptL1 );
      currentTree->SetBranchAddress( "ptL2",  &ptL2 );
      currentTree->SetBranchAddress( "etaL1", &etaL1 );
      currentTree->SetBranchAddress( "etaL2", &etaL2 );
      currentTree->SetBranchAddress( "combRelIsoLeg1DBeta",&combRelIsoLeg1DBeta);
      currentTree->SetBranchAddress( "diTauCharge",&diTauCharge);
      currentTree->SetBranchAddress( "MtLeg1",&MtLeg1);
      currentTree->SetBranchAddress( "pZeta",&pZeta);
      currentTree->SetBranchAddress( "pZetaVis",&pZetaVis);
      currentTree->SetBranchAddress( "MtLeg1Corr",&MtLeg1Corr);
      currentTree->SetBranchAddress( "pZetaCorr",&pZetaCorr);
      currentTree->SetBranchAddress( "pZetaVisCorr",&pZetaVisCorr);
      currentTree->SetBranchAddress( "numPV",&numPV);
      currentTree->SetBranchAddress( "sampleWeight",&sampleWeight);
      currentTree->SetBranchAddress( "puWeight",&puWeight);
      currentTree->SetBranchAddress( "ptVeto",&ptVeto);
      currentTree->SetBranchAddress( "MEt",&MEt);
      currentTree->SetBranchAddress( "MEtCorr",&MEtCorr);
      currentTree->SetBranchAddress( "jetsBtagHE1",&jetsBtagHE1);
      currentTree->SetBranchAddress( "jetsBtagHE2",&jetsBtagHE2);
      
      cout << "starting loop..." << endl;

      for (Long64_t ievt=0; ievt<currentTree->GetEntries();ievt++) {
	
	currentTree->GetEntry(ievt);
	
	if (ievt%10000 == 0){
	  std::cout << (it->first) << " ---> processing event: " << ievt << " ..." <<std::endl;
	}
	
	
	vMap["diTauNSVfitMass"]  = diTauNSVfitMass;
	vMap["diTauNSVfitMassErrUp"]  = diTauNSVfitMassErrUp;
	vMap["diTauNSVfitMassErrDown"]  = diTauNSVfitMassErrDown;
	vMap["diTauVisPt"]    = diTauVisPt;
	vMap["diTauVisEta"]   = diTauVisEta;
	vMap["diTauVisMass"]  = diTauVisMass;
	vMap["ptL1"]  = ptL1;
	vMap["ptL2"]  = ptL2;
	vMap["etaL1"] = etaL1;
	vMap["etaL2"] = etaL2;
	vMap["diTauCharge"]= Float_t(diTauCharge);
	vMap["MEt"]= MEt;
	vMap["MtLeg1"]= MtLeg1;
	vMap["pZeta"] = pZeta;
	vMap["pZetaVis"] = pZetaVis;
	vMap["pZetaCutVar"] = pZeta-1.5*pZetaVis;
	vMap["MEtCorr"]= MEtCorr;
	vMap["MtLeg1Corr"]= MtLeg1Corr;
	vMap["pZetaCorr"] = pZetaCorr;
	vMap["pZetaVisCorr"] = pZetaVisCorr;
	vMap["pZetaCutVarCorr"] = pZetaCorr-1.5*pZetaVisCorr;

	vMap["numPV"]= numPV;
	vMap["combRelIsoLeg1"]= combRelIsoLeg1;
	vMap["combRelIsoLeg1DBeta"]= combRelIsoLeg1DBeta;
	vMap["sampleWeight"]= sampleWeight;
	vMap["puWeight"]= puWeight;
	vMap["puWeight"]= puWeight;
	vMap["jetsBtagHE1"]= jetsBtagHE1;
	vMap["jetsBtagHE2"]= jetsBtagHE2;
	vMap["leadTrackPt"]= leadTrackPt;
	
	vMap["pt1"]= pt1;
	vMap["pt2"]= pt2;
	vMap["Deta"]= Deta;
	vMap["Mjj"]= Mjj;
	vMap["ptVeto"]= ptVeto;
	
	if((it->first).find("Data")==string::npos) sampleWeight*=(Lumi/1000*hltEff_*puWeight);
	
	counter++;
	h1->Fill( vMap[variable_], sampleWeight);
	
      }// end loop
    }

    if( (it->first).find("Zelecelec")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z#rightarrowll","F");
    }
    if( (it->first).find("Ztautau")!=string::npos ) {
      h1->SetFillColor(kWhite);
      leg->AddEntry(h1,"Z#rightarrow#tau#tau","F");
    }
    if( (it->first).find("TTbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"t#bar{t}+jets","F");
    }
    if( (it->first).find("Top")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"Top","F");
    }
    if( (it->first).find("SingleTop")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"single-t","F");
    }
    if( (it->first).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"W+jets","F");
    }
    if( (it->first).find("EWK")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"EWK","F");
    }
    if( (it->first).find("DiBoson")!=string::npos ) {
      h1->SetFillColor(38);
      leg->AddEntry(h1,"Di-Boson","F");
    }
    if( (it->first).find("QCD")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"QCD-multijets","F");
    }
    if( (it->first).find("SS")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"SS \\ W","F");
    }
    if((it->first).find("qqH115")!=string::npos){
      hSgn1 = (TH1F*)h1->Clone("hSgn1");
      hSgn1->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("VBF H(115) X %.0f",magnifySgn_),"F");
    }
    if((it->first).find("ggH115")!=string::npos){
      hSgn2 = (TH1F*)h1->Clone("hSgn2");
      hSgn2->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("GGF H(115) X %.0f",magnifySgn_),"F");
    }
    if((it->first).find("Data")!=string::npos){
      hData = (TH1F*)h1->Clone("hData");
      hData->Sumw2();
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.2);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetXTitle(XTitle_+" ("+Unities_+")");
      hData->SetYTitle(Form(" Events/(%.0f %s)", hData->GetBinWidth(1), Unities_.Data() ) );
      hData->SetTitleSize(0.04,"X");
      hData->SetTitleSize(0.05,"Y");
      hData->SetTitleOffset(0.95,"Y");
      leg->AddEntry(hData,"DATA","P");
    }

    if(iter==0) hSiml=(TH1F*)h1->Clone("hSiml");
    else if((it->first).find("Data")==string::npos) hSiml->Add(h1);

    if((it->first).find("Data")==string::npos) aStack->Add(h1);

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

  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  leg->Draw();

  //c1->SaveAs("tmp.png");

}


///////////////////////////////////////////////////////////////////////////////////////////////

void plotMuTau( TCut Cuts_ = "combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && diTauCharge==0 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20",
		  Bool_t doVBFPreRegion_ = false,
		  TString variable_ = "diTauVisMass",
		  TString XTitle_ = "mass",
		  TString Unities_ = "GeV",
		  Int_t nBins_ = 40, Float_t xMin_=0, Float_t xMax_=200,
		  Float_t magnifySgn_ = 1,
		  Float_t hltEff_ = 1.1,
		  Int_t enableHLTmatching_ = 0,
		  Int_t logy_ = 0 ) 
{   
  
  
  
  float Lumi = 145.;
  
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy_);

  TLegend* leg = new TLegend(0.55,0.47,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader("#mu+#tau_{had}");

  THStack* aStack = new THStack("aStack",Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}", Lumi ));
  TH1F* hSiml  = new TH1F();
  TH1F* hSgn   = new TH1F();
  TH1F* hSgn1  = new TH1F();
  TH1F* hSgn2  = new TH1F();
  TH1F* hData  = new TH1F();

  // OpenNTuples
  TString fDataName                = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleRun2011-Mu_Open_MuTauStream.root";
  TString fSignalNameVBF           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleVBFH125-Mu-powheg-PUS4_run_Open_MuTauStream.root";
  TString fSignalNameGGH           = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleGGFH125-Mu-powheg-PUS4_run_Open_MuTauStream.root";
  TString fBackgroundNameDYTauTau  = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleDYJets-Mu-50-madgraph-PUS4_run_Open_MuTauStream.root";
  TString fBackgroundNameDYMuMu    = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleDYJets-Mu-50-madgraph-PUS4_run_Open_MuTauStream.root";
  TString fBackgroundNameWJets     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleWJets-Mu-madgraph-PUS4_run_Open_MuTauStream.root";
  TString fBackgroundNameTTbar     = "/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamSummer11_iter1//nTupleTT-Mu-pythia-PUS3_run_Open_MuTauStream.root";

  TFile *fData              = new TFile( fDataName, "READ");  
  TFile *fSignalVBF         = new TFile(fSignalNameVBF ,"READ");  
  TFile *fSignalGGH         = new TFile(fSignalNameGGH,"READ");  
  TFile *fBackgroundDYTauTau= new TFile(fBackgroundNameDYTauTau,"READ"); 
  TFile *fBackgroundDYMuMu  = new TFile(fBackgroundNameDYMuMu,"READ"); 
  TFile *fBackgroundWJets   = new TFile(fBackgroundNameWJets,"READ"); 
  TFile *fBackgroundTTbar   = new TFile(fBackgroundNameTTbar,"READ"); 

  if(!fSignalVBF || 
     !fBackgroundDYTauTau || 
     !fBackgroundWJets || 
     !fBackgroundTTbar ||
     !fSignalGGH || 
     !fBackgroundDYMuMu ||
     !fData ){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTreePtOrd";
  TCut inclusiveRegionSS("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && ( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) ))");
  TCut vbfPreRegionSS("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && ( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) ))");
  TCut inclusiveRegionSSMtSide("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)<-40 && diTauCharge!=0 && ( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) ))");
  TCut vbfPreRegionSSMtSide("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)<-40 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && ( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) ))");
  TCut inclusiveRegionSSmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && HLTx==1 )");
  TCut vbfPreRegionSSmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && (pZetaCorr-1.5*pZetaVisCorr)>-20 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1 )");
  TCut inclusiveRegionSSMtmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && diTauCharge!=0 && HLTx==1 )");
  TCut vbfPreRegionSSMtmc("(combIsoLeg2<2 && combRelIsoLeg1DBeta<0.10 && muFlag==0 && diTauCharge!=0 && pt1>30 && pt2>30 && eta1*eta2<0 && Mjj>350 && Deta>3.5 && HLTx==1 )");

  TTree *data                = (TTree*)fData->Get(tree);
  TTree *signalVBF           = (TTree*)fSignalVBF->Get(tree);
  TTree *signalGGH           = (TTree*)fSignalGGH->Get(tree);

  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDYTauTau->Get(tree))->CopyTree("abs(genDecay)==(23*15)");
  TTree *backgroundDYMuMu    = ((TTree*)fBackgroundDYMuMu->Get(tree))->CopyTree("abs(genDecay)==(23*13)");
  cout <<backgroundDYTauTau->GetEntries() << " -- " << backgroundDYMuMu->GetEntries() << endl;

  TTree* dataSS              = 0;
  TTree* WJetsSS             = 0;
  TTree* WJetsSSFullMt       = 0;
  dataSS  = doVBFPreRegion_ ? 
    ((TTree*)fData->Get(tree))->CopyTree(vbfPreRegionSS) : 
    ((TTree*)fData->Get(tree))->CopyTree(inclusiveRegionSS);
  WJetsSS = doVBFPreRegion_ ? 
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(vbfPreRegionSSmc) :
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(inclusiveRegionSSmc); 
  WJetsSSFullMt = doVBFPreRegion_ ? 
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(vbfPreRegionSSMtmc) :
    ((TTree*)fBackgroundWJets->Get(tree))->CopyTree(inclusiveRegionSSMtmc);
  
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);

  // here I choose the order in the stack
  std::vector<string> samples;
  //samples.push_back("ggH115");
  //samples.push_back("qqH115");
  samples.push_back("TTbar");
  samples.push_back("Wjets");
  samples.push_back("SS");
  //samples.push_back("Zmumu");
  samples.push_back("Ztautau");
  samples.push_back("Data");
  

  // here I define the map between a sample name and its tree
  std::map<std::string,TTree*> tMap;
  tMap["Data"]=data;
  tMap["ggH115"]=signalGGH;
  tMap["qqH115"]=signalVBF;
  tMap["Ztautau"]=backgroundDYTauTau;
  tMap["Zmumu"]=backgroundDYMuMu;
  tMap["Wjets"]=backgroundWJets;
  tMap["TTbar"]=backgroundTTbar;
  tMap["SS"]=dataSS;


  Float_t pt1, pt2, Deta, Mjj, eta1, eta2;
  Float_t diTauNSVfitMass,diTauNSVfitMassErrUp,diTauNSVfitMassErrDown,diTauVisPt,diTauVisEta,diTauVisMass,ptL1,ptL2,etaL1,etaL2,diTauCharge,MtLeg1,numPV,combRelIsoLeg1DBeta,combRelIsoLeg1,pZeta,pZetaVis, MEt, pZetaCorr,pZetaVisCorr, MEtCorr, MtLeg1Corr, ptVeto, leadTrackPt, jetsBtagHE1, jetsBtagHE2;
  Float_t sampleWeight,puWeight, HLTweight;

  std::map<TString,Float_t> vMap;


  for( unsigned iter=0; iter<samples.size(); iter++){

    cout << samples[iter] << endl;
    
    std::map<std::string,TTree*>::iterator it = tMap.find(samples[iter]);

    TString h1Name = "h1_"+it->first;
    TH1F* h1 = new TH1F( h1Name ,Form("CMS Preliminary 2011    #sqrt{s}=7 TeV L=%.0f pb^{-1}", Lumi) , nBins_ ,xMin_ , xMax_);
    
    TFile* dummy = new TFile("dummy.root","RECREATE");  
    TCut Cuts = Cuts_;
    TCut hltData   = "( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) )";   
    TCut hltMC     = "(HLTx==1)";
    if((it->first).find("Data")!=string::npos) Cuts = Cuts_ && hltData;
    else Cuts = Cuts_ && hltMC;
    
    TTree* currentTree = 0;
    
    if((it->first).find("SS")!=string::npos){
      
      cout << "Dealing with SS" << endl;
      
      currentTree = (it->second);
      
      //cout << currentTree->GetEntries() << endl;
      TH1F* hHelp = (TH1F*)h1->Clone("hHelp");
      hHelp->Reset();
      TString variable = ( (std::string(variable_.Data())).find("pZeta")!=string::npos ) ? "pZetaCorr-1.5*pZetaVisCorr" : variable_;
      currentTree->Draw(variable+">>hHelp");
      //cout << hHelp->GetEntries() << endl;
      h1->Add(hHelp,1);
      //cout << h1->Integral() << endl;
      TH1F* hWMt = new TH1F("hWMt","",200,-200,200);
      WJetsSSFullMt->Draw("(pZetaCorr-1.5*pZetaVisCorr)>>hWMt");
      float scaleFactor = (hWMt->Integral(0,80))/(hWMt->Integral(90,200));
      cout << "Scale factor = " << scaleFactor << endl;
      TH1F* hW = (TH1F*)h1->Clone("hW");
      hW->Reset();
      float SSWinSignalRegion = doVBFPreRegion_ ?
	data->GetEntries(inclusiveRegionSSMtSide) :  
	data->GetEntries(vbfPreRegionSSMtSide);
      SSWinSignalRegion *= (1./scaleFactor);
      WJetsSS->Draw(variable+">>hW");
      hW->Scale(SSWinSignalRegion/hW->Integral());
      h1->Add(hW, -1 );
      h1->Scale(1.06);
      //cout << h1->Integral() << endl;
    }
    else{

      currentTree =  (TTree*)(it->second)->CopyTree(Cuts);
      
      Int_t counter = 0;
      
      currentTree->SetBranchAddress( "pt1", &pt1 );
      currentTree->SetBranchAddress( "pt2", &pt2 );
      currentTree->SetBranchAddress( "Deta",&Deta );
      currentTree->SetBranchAddress( "Mjj", &Mjj );
      currentTree->SetBranchAddress( "diTauNSVfitMass",&diTauNSVfitMass);
      currentTree->SetBranchAddress( "diTauNSVfitMassErrUp",&diTauNSVfitMassErrUp);
      currentTree->SetBranchAddress( "diTauNSVfitMassErrDown",&diTauNSVfitMassErrDown);

      currentTree->SetBranchAddress( "diTauVisPt",&diTauVisPt);
      currentTree->SetBranchAddress( "diTauVisEta",&diTauVisEta);
      currentTree->SetBranchAddress( "diTauVisMass",&diTauVisMass);
      currentTree->SetBranchAddress( "ptL1",  &ptL1 );
      currentTree->SetBranchAddress( "ptL2",  &ptL2 );
      currentTree->SetBranchAddress( "etaL1", &etaL1 );
      currentTree->SetBranchAddress( "etaL2", &etaL2 );
      currentTree->SetBranchAddress( "combRelIsoLeg1DBeta",&combRelIsoLeg1DBeta);
      currentTree->SetBranchAddress( "diTauCharge",&diTauCharge);
      currentTree->SetBranchAddress( "MtLeg1",&MtLeg1);
      currentTree->SetBranchAddress( "pZeta",&pZeta);
      currentTree->SetBranchAddress( "pZetaVis",&pZetaVis);
      currentTree->SetBranchAddress( "MtLeg1Corr",&MtLeg1Corr);
      currentTree->SetBranchAddress( "pZetaCorr",&pZetaCorr);
      currentTree->SetBranchAddress( "pZetaVisCorr",&pZetaVisCorr);
      currentTree->SetBranchAddress( "numPV",&numPV);
      currentTree->SetBranchAddress( "sampleWeight",&sampleWeight);
      currentTree->SetBranchAddress( "puWeight",&puWeight);
      currentTree->SetBranchAddress( "ptVeto",&ptVeto);
      currentTree->SetBranchAddress( "MEt",&MEt);
      currentTree->SetBranchAddress( "MEtCorr",&MEtCorr);
      currentTree->SetBranchAddress( "jetsBtagHE1",&jetsBtagHE1);
      currentTree->SetBranchAddress( "jetsBtagHE2",&jetsBtagHE2);
      
      cout << "starting loop..." << endl;

      for (Long64_t ievt=0; ievt<currentTree->GetEntries();ievt++) {
	
	currentTree->GetEntry(ievt);
	
	if (ievt%10000 == 0){
	  std::cout << (it->first) << " ---> processing event: " << ievt << " ..." <<std::endl;
	}
	
	
	vMap["diTauNSVfitMass"]  = diTauNSVfitMass;
	vMap["diTauNSVfitMassErrUp"]  = diTauNSVfitMassErrUp;
	vMap["diTauNSVfitMassErrDown"]  = diTauNSVfitMassErrDown;
	vMap["diTauVisPt"]    = diTauVisPt;
	vMap["diTauVisEta"]   = diTauVisEta;
	vMap["diTauVisMass"]  = diTauVisMass;
	vMap["ptL1"]  = ptL1;
	vMap["ptL2"]  = ptL2;
	vMap["etaL1"] = etaL1;
	vMap["etaL2"] = etaL2;
	vMap["diTauCharge"]= Float_t(diTauCharge);
	vMap["MEt"]= MEt;
	vMap["MtLeg1"]= MtLeg1;
	vMap["pZeta"] = pZeta;
	vMap["pZetaVis"] = pZetaVis;
	vMap["pZetaCutVar"] = pZeta-1.5*pZetaVis;
	vMap["MEtCorr"]= MEtCorr;
	vMap["MtLeg1Corr"]= MtLeg1Corr;
	vMap["pZetaCorr"] = pZetaCorr;
	vMap["pZetaVisCorr"] = pZetaVisCorr;
	vMap["pZetaCutVarCorr"] = pZetaCorr-1.5*pZetaVisCorr;

	vMap["numPV"]= numPV;
	vMap["combRelIsoLeg1"]= combRelIsoLeg1;
	vMap["combRelIsoLeg1DBeta"]= combRelIsoLeg1DBeta;
	vMap["sampleWeight"]= sampleWeight;
	vMap["puWeight"]= puWeight;
	vMap["puWeight"]= puWeight;
	vMap["jetsBtagHE1"]= jetsBtagHE1;
	vMap["jetsBtagHE2"]= jetsBtagHE2;
	vMap["leadTrackPt"]= leadTrackPt;
	
	vMap["pt1"]= pt1;
	vMap["pt2"]= pt2;
	vMap["Deta"]= Deta;
	vMap["Mjj"]= Mjj;
	vMap["ptVeto"]= ptVeto;
	
	if((it->first).find("Data")==string::npos) sampleWeight*=(Lumi/1000*hltEff_*puWeight);
	
	counter++;
	h1->Fill( vMap[variable_], sampleWeight);
	
      }// end loop
    }

    if( (it->first).find("Zmumu")!=string::npos ) {
      h1->SetFillColor(7);
      leg->AddEntry(h1,"Z#rightarrowll","F");
    }
    if( (it->first).find("Ztautau")!=string::npos ) {
      h1->SetFillColor(kWhite);
      leg->AddEntry(h1,"Z#rightarrow#tau#tau","F");
    }
    if( (it->first).find("TTbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"t#bar{t}+jets","F");
    }
    if( (it->first).find("Top")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"Top","F");
    }
    if( (it->first).find("SingleTop")!=string::npos ) {
      h1->SetFillColor(29);
      leg->AddEntry(h1,"single-t","F");
    }
    if( (it->first).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"W+jets","F");
    }
    if( (it->first).find("EWK")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"EWK","F");
    }
    if( (it->first).find("DiBoson")!=string::npos ) {
      h1->SetFillColor(38);
      leg->AddEntry(h1,"Di-Boson","F");
    }
    if( (it->first).find("QCD")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"QCD-multijets","F");
    }
    if( (it->first).find("SS")!=string::npos ) {
      h1->SetFillColor(42);
      leg->AddEntry(h1,"SS \\ W","F");
    }
    if((it->first).find("qqH115")!=string::npos){
      hSgn1 = (TH1F*)h1->Clone("hSgn1");
      hSgn1->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn1->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3004);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("VBF H(115) X %.0f",magnifySgn_),"F");
    }
    if((it->first).find("ggH115")!=string::npos){
      hSgn2 = (TH1F*)h1->Clone("hSgn2");
      hSgn2->Scale(magnifySgn_);
      h1->Scale(magnifySgn_);
      hSgn2->SetLineWidth(2);
      h1->SetFillColor(kBlack);
      h1->SetFillStyle(3005);
      h1->SetLineColor(kBlack);
      leg->AddEntry(h1,Form("GGF H(115) X %.0f",magnifySgn_),"F");
    }
    if((it->first).find("Data")!=string::npos){
      hData = (TH1F*)h1->Clone("hData");
      hData->Sumw2();
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.2);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      hData->SetXTitle(XTitle_+" ("+Unities_+")");
      hData->SetYTitle(Form(" Events/(%.0f %s)", hData->GetBinWidth(1), Unities_.Data() ) );
      hData->SetTitleSize(0.04,"X");
      hData->SetTitleSize(0.05,"Y");
      hData->SetTitleOffset(0.95,"Y");
      leg->AddEntry(hData,"DATA","P");
    }

    if(iter==0) hSiml=(TH1F*)h1->Clone("hSiml");
    else if((it->first).find("Data")==string::npos) hSiml->Add(h1);

    if((it->first).find("Data")==string::npos) aStack->Add(h1);

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

  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle(XTitle_+" ("+Unities_+")");
  hStack->SetYTitle(Form(" Events/(%.0f %s)", hStack->GetBinWidth(1), Unities_.Data() ) );
  hStack->SetTitleSize(0.04,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.95,"Y");
  leg->Draw();

  //c1->SaveAs("tmp.png");

}












