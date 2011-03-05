#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TBenchmark.h"
#include "TGraph.h"
#include "TVectorT.h"


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


  TPad* pad1 = new TPad("pad1DEta","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.26);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(1);

  TLegend* leg = new TLegend(0.60,0.55,0.85,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader("N^{jet}>1, p^{jet}_{T}>30 GeV");
  
  THStack* aStack = new THStack("aStack",Form("CMS Preliminary    #sqrt{s}=7 TeV L=%.0f pb^{-1}",Lumi_));
  TH1F* hData = new TH1F();
  TH1F* hSiml = new TH1F();

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 16 ,0 , 8);

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
    if( ((fileList_[i].second).first).find("data")!=string::npos ) {
      h1->SetMarkerColor(kBlack);
      h1->SetMarkerStyle(kFullCircle);
      h1->SetMarkerSize(0.8);
      leg->AddEntry(h1,"DATA","lep");
    }

    currentTree->Draw( Form("abs(jetsIDP4[0].eta()-jetsIDP4[1].eta())>>%s",h1Name.c_str()), 
		       "jetsIDP4@.size()>1 && jetsIDP4[jetsIDP4@.size()-1].pt()>30" );

    if(((fileList_[i].second).first).find("data")!=string::npos){
      hData = h1;
      hData->Sumw2();
      continue;
    }

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );
   
    if(i==0) hSiml=(TH1F*)h1->Clone("hSiml");
    else hSiml->Add(h1);

    aStack->Add(h1);

  }

  aStack->Draw("HIST");
  hData->Draw("PSAME");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("|#Delta #eta_{jj}|");
  hStack->SetYTitle(Form(" Events/(%.1f units)", hStack->GetBinWidth(1) ) );
  hStack->SetTitleSize(0.05,"X");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleOffset(0.75,"Y");
  leg->Draw();

  pad2->cd();
  TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA-MC)}{MC}",
			  hStack->GetNbinsX(), 
			  hStack->GetXaxis()->GetXmin(), hStack->GetXaxis()->GetXmax());
  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerSize(0.8);
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");
  TH1F* hDataClone = (TH1F*)hData->Clone("hDataClone");
  hDataClone->Add(hSiml,-1);
  hRatio->Divide( hDataClone ,hSiml,1.0,1.0);
  hRatio->SetAxisRange(-3,3,"Y");

  hRatio->Draw();
  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hStack->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

}


void plotJetMultiplicity(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30.){

  TCanvas *c1 = new TCanvas("c1JetMultiplicity","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad("pad1JetMultiplicity","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2JetMultiplicity","",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("p^{jet}_{T}>%.0f GeV/c",cutOff_) );
  
  THStack* aStack = new THStack("aStack",Form("CMS Preliminary    #sqrt{s}=7 TeV L=%.0f pb^{-1}",Lumi_));
  TH1F* hData = new TH1F();
  TH1F* hSiml = new TH1F();

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 7 ,-0.5, 6.5);

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
    if( ((fileList_[i].second).first).find("data")!=string::npos ) {
      h1->SetMarkerColor(kBlack);
      h1->SetMarkerStyle(kFullCircle);
      h1->SetMarkerSize(0.8);
      leg->AddEntry(h1,"DATA","lep");
    }

    currentTree->Draw( Form("jetsIDP4@.size()>>%s",h1Name.c_str()), 
		       Form("jetsIDP4@.size()<1 || (jetsIDP4@.size()>0 && jetsIDP4[jetsIDP4@.size()-1].pt()>%f)",cutOff_) );

    if(((fileList_[i].second).first).find("data")!=string::npos){
      hData = h1;
      hData->Sumw2();
      continue;
    }

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("Exclusive jet number");
  hStack->SetYTitle("Number of events");
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


void plotBtag(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30.){

  TCanvas *c1 = new TCanvas("c1Btag","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad("pad1Btag","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2Btag","",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("p^{jet}_{T}>%.0f GeV/c",cutOff_) );
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 40 ,-20, 20);

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

    currentTree->Draw( Form("jetsBtagHE>>%s",h1Name.c_str()), 
		       Form("jetsIDP4@.size()<1 || (jetsIDP4@.size()>0 && jetsIDP4[jetsIDP4@.size()-1].pt()>%f)",cutOff_) );

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    h1->SetBinContent(1,h1->GetBinContent(0));
    h1->SetBinContent( h1->GetNbinsX() , h1->GetBinContent(h1->GetNbinsX()+1));

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("TCHE output");
  hStack->SetYTitle("Number of events");
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


void plotBtag2(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30.){

  TCanvas *c1 = new TCanvas("c1Btag2","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad("pad1Btag2","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2Btag2","",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("p^{jet}_{T}>%.0f GeV/c",cutOff_) );
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 40 ,-20, 20);

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

    currentTree->Draw( Form("jetsBtagHP>>%s",h1Name.c_str()), 
		       Form("jetsIDP4@.size()<1 || (jetsIDP4@.size()>0 && jetsIDP4[jetsIDP4@.size()-1].pt()>%f)",cutOff_) );

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    h1->SetBinContent(1,h1->GetBinContent(1)+h1->GetBinContent(0));
    h1->SetBinContent( h1->GetNbinsX() , h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("TCHP output");
  hStack->SetYTitle("Number of events");
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


void plotLeadingJetEt(  FileList fileList_  , float Lumi_ = 30 , float etaCutLow_=0.0, float etaCutHigh_=5.0){

  TCanvas *c1 = new TCanvas(Form("c1LeadingJetEt%f",etaCutLow_),"",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad(Form("pad1LeadingJetEt%f",etaCutLow_),"",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad(Form("pad2LeadingJetEt%f",etaCutLow_),"",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("%.1f<|#eta^{jet}|<%.1f",etaCutLow_,etaCutHigh_) );
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 36 ,20, 200);

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

    currentTree->Draw( Form("jetsIDP4[0].Et()>>%s",h1Name.c_str()), 
		       Form("jetsIDP4@.size()>0 && abs(jetsIDP4[0].eta())<%f && abs(jetsIDP4[0].eta())>%f",etaCutHigh_,etaCutLow_ ) );
    
    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("E_{T} leading jet");
  hStack->SetYTitle("Events/(5 GeV/c)");
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


void plotZpt(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30., int jetMult_ = 99, int nBins_ = 40){

  TCanvas *c1 = new TCanvas(Form("c1Zpt%d",jetMult_),"",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad(Form("pad1Zpt%d",jetMult_),"",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad(Form("pad2Zpt%d",jetMult_),"",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("p^{jet}_{T}>%.0f GeV/c, N^{jets}=%d ",cutOff_,jetMult_) );
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", nBins_ ,0, 200);

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

    currentTree->Draw( Form("diMuonP4[0].Pt()>>%s",h1Name.c_str()), 
		       Form("(jetsIDP4@.size()>%f && jetsIDP4@.size()<%f && ( (jetsIDP4@.size()<1) || (jetsIDP4@.size()>0 && jetsIDP4[jetsIDP4@.size()-1].pt()>%f)) )",(float)(jetMult_-0.5),(float)(jetMult_+0.5),cutOff_) );

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    h1->SetBinContent( h1->GetNbinsX() , h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("p_{T}^{#mu#mu}");
  hStack->SetYTitle( Form("Events/(%.0f GeV/c)",hStack->GetBinWidth(1)) );
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


void plotZeppZ(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30.){

  TCanvas *c1 = new TCanvas("c1ZeppZ","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);


  TPad* pad1 = new TPad("pad1ZeppZ","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2ZeppZ","",0.05,0.02,0.96,0.26);
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
  leg->SetHeader( Form("p^{jet}_{T}>%.0f GeV/c",cutOff_) );
  
  THStack* aStack = new THStack("aStack",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_));

  for(unsigned int i = 0 ; i < fileList_.size() ; i++){

    TFile* currentFile = (TFile*)fileList_[i].first ;
    if( currentFile->IsZombie() ) continue;
    TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
    float totalEvents = allEvents->GetBinContent(1);

    TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
    string h1Name = "h1_"+(fileList_[i].second).first;
    TH1F* h1 = new TH1F( h1Name.c_str() ,"", 28 ,-7, 7);

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

    currentTree->Draw( Form("diMuonP4[0].Eta()-(jetsIDP4[0].Eta()+jetsIDP4[1].Eta())/2.>>%s",h1Name.c_str()), 
		       Form("jetsIDP4@.size()>1 && jetsIDP4[jetsIDP4@.size()-1].pt()>%f",cutOff_) );

    h1->Scale( Lumi_ / (totalEvents/((fileList_[i].second).second)) );

    aStack->Add(h1);

  }


  aStack->Draw("HIST");
  TH1F* hStack = (TH1F*)aStack->GetHistogram();
  hStack->SetXTitle("#eta^{#mu#mu}-#frac{#eta_{j1}+#eta_{j2}}{2}");
  hStack->SetYTitle( Form("Events/(%.1f unities)",hStack->GetBinWidth(1)) );
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





void plotJetVeto(  FileList fileList_  , float Lumi_ = 30 , float cutOff_ = 30. , float minDEta_ = 1.0){

  TCanvas *c1 = new TCanvas("c1JetVeto","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(0);

  TPad* pad1 = new TPad("pad1JetVeto","",0.05,0.27,0.96,0.97);
  TPad* pad2 = new TPad("pad2JetVeto","",0.05,0.02,0.96,0.26);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(0);

  TLegend* leg = new TLegend(0.15,0.60,0.35,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader( Form("p^{j1,j2}_{T}>%.0f GeV/c ; |#eta_{j1}-#eta_{j2}|>%.1f ; (#eta_{j1,j2}+ 0.5) < #eta_{j3} < (#eta_{j2,j1}- 0.5)",cutOff_,minDEta_) );

  char* ptCut[3] = {"25","30","35"};
 
  char* numPV[7] = {"1","2","3","4","5","6","7"};
  float vxPV[7]={1,2,3,4,5,6,7};
  float vexPV[7]={0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float vyPV15[7];
  float vyPV25[7];
  float vyPV35[7];
  float veyPV15[7];
  float veyPV25[7];
  float veyPV35[7];

  for(unsigned int i = 0 ; i< 3; i++){

    for(unsigned int j = 0 ; j< 7; j++){
      
      float pass_i_j = 0;
      float norm_i_j = 0; 
      float pass_i_j_num = 0; 
      float norm_i_j_num = 0;

      for(unsigned int f = 0 ; f < fileList_.size() ; f++){
	
	TFile* currentFile = (TFile*)fileList_[f].first ;
	if( currentFile->IsZombie() ) continue;
	TH1F* allEvents = (TH1F*)currentFile->Get("allEventsFilter/totalEvents");
	float totalEvents = allEvents->GetBinContent(1);
	
	TTree* currentTree = (TTree*)currentFile->Get("zPlusJetsAnalyzer/tree");
	string h1Name = "h1_"+(fileList_[f].second).first;
	TH1F* h1Pass = new TH1F( "h1Pass" ,"", 1 ,55, 125);
	TH1F* h1Norm = new TH1F( "h1Norm" ,"", 1 ,55, 125);
	
	currentTree->Draw("Zmass>>h1Pass", 
			  Form("((jetsIDP4@.size()>1 && jetsIDP4@.size()<3 && abs(jetsIDP4[0].eta()-jetsIDP4[1].eta())>%f && jetsIDP4[1].Et()>%f && numPV>(%s-0.5) && numPV<(%s+0.5)) || (jetsIDP4@.size()>2 && abs(jetsIDP4[0].eta()-jetsIDP4[1].eta())>%f && jetsIDP4[1].Et()>%f && jetsIDP4[2].Et()<%s && ( (jetsIDP4[2].eta()>jetsIDP4[1].eta()+0.5 &&  jetsIDP4[2].eta()<jetsIDP4[0].eta()-0.5) || (jetsIDP4[2].eta()>jetsIDP4[0].eta()+0.5 &&  jetsIDP4[2].eta()<jetsIDP4[1].eta()-0.5) ) && numPV>(%s-0.5) && numPV<(%s+0.5)))",minDEta_,cutOff_,numPV[j],numPV[j],minDEta_,cutOff_,ptCut[i],numPV[j],numPV[j] ) );
	
	currentTree->Draw("Zmass>>h1Norm", 
			  Form("(jetsIDP4@.size()>1 && abs(jetsIDP4[0].eta()-jetsIDP4[1].eta())>%f && jetsIDP4[1].Et()>%f && numPV>(%s-0.5) && numPV<(%s+0.5))",minDEta_,cutOff_,numPV[j],numPV[j]) );
	
	h1Pass->Scale(  Lumi_ / (totalEvents/((fileList_[f].second).second)) );
	h1Norm->Scale(  Lumi_ / (totalEvents/((fileList_[f].second).second)) );

	pass_i_j+=h1Pass->Integral();
	pass_i_j_num+=h1Pass->GetEntries();
	norm_i_j+=h1Norm->Integral();
	norm_i_j_num+=h1Norm->GetEntries();

	delete h1Norm;
	delete h1Pass;
	
      }//end f

      cout << pass_i_j << "   /   " << norm_i_j << endl;

      if(ptCut[i]=="25"){
	vyPV15[j]  =  pass_i_j/norm_i_j;
	veyPV15[j] = sqrt( vyPV15[j]*(1-vyPV15[j]) / norm_i_j_num);
      }
      if(ptCut[i]=="30"){
	vyPV25[j] =  pass_i_j/norm_i_j;
	veyPV25[j] = sqrt( vyPV25[j]*(1-vyPV25[j]) / norm_i_j_num);
      }
      if(ptCut[i]=="35") {
	vyPV35[j] =  pass_i_j/norm_i_j;
	veyPV35[j] = sqrt( vyPV35[j]*(1-vyPV35[j]) / norm_i_j_num);
      }

    }//end j

  }//end i

 
  TH2F* hMaster = new TH2F("hMaster",Form("#sqrt{s}=7 TeV L=%.0f pb^{-1}   CMS Preliminary",Lumi_),8,0.5,7.5,100,0,1);
  hMaster->SetXTitle("num PV");
  hMaster->SetYTitle( "jet veto efficiency" );
  hMaster->SetTitleSize(0.05,"X");
  hMaster->SetTitleSize(0.05,"Y");
  hMaster->SetTitleOffset(0.85,"Y");
  hMaster->Draw();

  TVectorF TvyPV15(7,vyPV15);
  TVectorF TvxPV15(7,vxPV);
  TVectorF TvyPV25(7,vyPV25);
  TVectorF TvxPV25(7,vxPV);
  TVectorF TvyPV35(7,vyPV35);
  TVectorF TvxPV35(7,vxPV);
  TVectorF TveyPV15(7,veyPV15);
  TVectorF TvexPV15(7,vexPV);
  TVectorF TveyPV25(7,veyPV25);
  TVectorF TvexPV25(7,vexPV);
  TVectorF TveyPV35(7,veyPV35);
  TVectorF TvexPV35(7,vexPV);

  TGraphErrors* graph15 = new TGraphErrors(TvxPV15,TvyPV15,TvexPV15,TveyPV15);
  graph15->SetMarkerStyle(kOpenCircle);
  graph15->SetMarkerSize(1.2);
  graph15->SetMarkerColor(kRed);

  TGraphErrors* graph25 = new TGraphErrors(TvxPV25,TvyPV25,TvexPV25,TveyPV25);
  graph25->SetMarkerStyle(kOpenTriangleUp);
  graph25->SetMarkerSize(1.2);
  graph25->SetMarkerColor(kBlue);

  TGraphErrors* graph35 = new TGraphErrors(TvxPV35,TvyPV35,TvexPV35,TveyPV35);
  graph35->SetMarkerStyle(kOpenSquare);
  graph35->SetMarkerSize(1.2);
  graph35->SetMarkerColor(kMagenta);

  graph15->Draw("P");
  graph25->Draw("P");
  graph35->Draw("P");

  leg->AddEntry(graph15,"E_{T}^{veto}=25 GeV","P");
  leg->AddEntry(graph25,"E_{T}^{veto}=30 GeV","P");
  leg->AddEntry(graph35,"E_{T}^{veto}=35 GeV","P");

  leg->Draw();


  pad2->cd();
  TH1F* hRatio = new TH1F("hRatio", " ; ; #frac{(DATA-MC)}{MC}",
			  hMaster->GetNbinsX(), 
			  hMaster->GetXaxis()->GetXmin(), hMaster->GetXaxis()->GetXmax());
  hRatio->SetLabelSize(0.12,"X");
  hRatio->SetLabelSize(0.10,"Y");
  hRatio->SetTitleSize(0.12,"Y");
  hRatio->SetTitleOffset(0.36,"Y");

  hRatio->Draw();
}



void mainPlot(){

  TFile* fData   = new TFile("/data_CMS/cms/lbianchini/ZmumuPlusJetsStudy/treeZmumuPlusJets_Mu-Run2010AB.root","READ");
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
  fileList.push_back( make_pair(fData,    make_pair("data",       -99 )  ));
  
  plotDEta( fileList , 36.);
  /*
  
  plotBtag( fileList, 36. , 30.);
  plotBtag2( fileList, 36. , 30.);
  
  plotLeadingJetEt( fileList, 36. , 0.0, 2.4);
  plotLeadingJetEt( fileList, 36. , 2.4, 4.5);
 
  plotZpt( fileList, 36. , 30., 0, 40);
  plotZpt( fileList, 36. , 30., 1, 20);
  plotZpt( fileList, 36. , 30., 2, 20);
  plotZpt( fileList, 36. , 30., 3, 15);
  plotZpt( fileList, 36. , 30., 4, 10);
  plotZpt( fileList, 36. , 30., 5,  5);
  plotZpt( fileList, 36. , 30., 6,  5);
  

  plotZeppZ(  fileList, 36. , 30.);

  

  plotJetMultiplicity( fileList, 36. , 30.);
  plotJetVeto(  fileList, 36. , 30., 2.0 ); 
  */

}


