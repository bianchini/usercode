#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TBenchmark.h"

#include <vector>

using namespace std;

typedef std::vector< pair<TFile*, pair<std::string,float> > > FileList; 

void plotDEta(  FileList fileList_  , float Lumi_ = 30 ){

  TCanvas *c1 = new TCanvas("c1DEta","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad("pad1","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2","",0.05,0.02,0.96,0.26);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(1);

  TLegend* leg = new TLegend(0.60,0.50,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader(">1 jet, p^{jet}_{T}>30 GeV");
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,Form(" ; |#Delta #eta| ; Events/(%.0f units)", (10./25.) ), 25 ,0,10);

    if( ((fileList_[i].second).first).find("Zjets")!=string::npos ) {
      h1->SetFillColor(kRed);
      leg->AddEntry(h1,"MadGraph Z+jets","F");
    }
    if( ((fileList_[i].second).first).find("ttbar")!=string::npos ) {
      h1->SetFillColor(kBlue);
      leg->AddEntry(h1,"MadGraph t#bar{t}+jets","F");
    }
    if( ((fileList_[i].second).first).find("Wjets")!=string::npos ) {
      h1->SetFillColor(kGreen);
      leg->AddEntry(h1,"MadGraph W+jets","F");
    }
    if( ((fileList_[i].second).first).find("tW")!=string::npos ){
      h1->SetFillColor(kYellow);
      leg->AddEntry(h1,"MadGraph single-t","F");
    }
    if( ((fileList_[i].second).first).find("QCD")!=string::npos ) {
      h1->SetFillColor(kBlack);
      leg->AddEntry(h1,"PYTHIA QCD","F");
    }

    currentTree->Draw( Form("abs(jetsIDP4[0].eta()-jetsIDP4[1].eta())>>%s",h1Name.c_str()), 
		       "jetsIDP4@.size()>1 && jetsIDP4[jetsIDP4@.size()-1].pt()>30" );

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("|#Delta #eta|");
  hStack->SetYTitle(Form(" Events/(%.2f units)", (10./25.) ) );
  hStack->SetTitleSize(0.05,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.75,"Y");
  leg->Draw();

  pad2->cd();
  TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA-MC)}{MC}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");

  hRatio->Draw();
}


void mainPlot(){

  TFile* fDYJets = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_DYJets-madgraph-50-PU-v2.root","READ");
  TFile* fTT     = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_TT-madgraph-PU-v2.root","READ");
  TFile* fWJets  = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_WJets-madgraph-PU-v2.root","READ");
  TFile* fT      = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_TToBLNu-tW-madhraph-PU-v2.root","READ");
  TFile* fQCD    = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_QCD-pythia-PU-v2.root","READ");

  FileList fileList;
  fileList.push_back( make_pair(fQCD,    make_pair("QCD",   349988.0 )  ));
  fileList.push_back( make_pair(fWJets,  make_pair("Wjets",  31314.0 )  ));
  fileList.push_back( make_pair(fT,      make_pair("tW",        10.6 )  ));
  fileList.push_back( make_pair(fTT,     make_pair("ttbar",    157.5 )  ));
  fileList.push_back( make_pair(fDYJets, make_pair("Zjets",   3048.0 )  ));

  plotDEta( fileList , 30);


}


