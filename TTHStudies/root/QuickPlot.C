#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"

#include "TPaveText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TLine.h"

#include "TVector3.h"
#include "TLorentzVector.h"






void plot_analysis(TString cat   ="SL", 
		   string header="type==3",
		   TString title = "",
		   TString fname= "type0",
		   float fact = 1, float fact2 = 1
		   ){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  //TTJetsSemiLept
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_newType3_TTH125.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_newType3_TTJetsSemiLept.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTH125.root");

  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTJetsSemiLept.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTJetsSemiLept.root");

  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTH125.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTH125.root");

  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTJetsSemiLept.root");
  
  TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTJetsSemiLept.root");
  TFile* fB = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTJetsSemiLept.root");



  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",2, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",2, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();

  //tS->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hS", fact));
  //tB->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hB", fact));

  //tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), TCut((header+" && nSimBs>2").c_str()));
  //tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), TCut((header+" && nSimBs==2").c_str()));

  tS->Draw(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hS", fact), TCut((header+" && nSimBs>2  && p_125_all_s_ttbb/(p_125_all_s_ttbb+0.02*((1-0.5)*p_125_all_b_ttbb+0.5*p_125_all_b_ttjj))<0.1").c_str()));
  tB->Draw(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hB", fact), TCut((header+" && nSimBs==2 && p_125_all_s_ttbb/(p_125_all_s_ttbb+0.02*((1-0.5)*p_125_all_b_ttbb+0.5*p_125_all_b_ttjj))<0.1").c_str()));

  //tS->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hS", fact, fact2, fact2), TCut((header+" && nSimBs>2" ).c_str()));
  //tB->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hB", fact, fact2, fact2), TCut((header+" && nSimBs==2").c_str()));

  cout << "hS=" << hS->Integral() << endl;
  cout << "hB=" << hB->Integral() << endl;


  if(hS->GetMaximum()>=hB->GetMaximum()){
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hS->SetMaximum(1.0);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hB->SetMaximum(1.0);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hS,"ttH");
  leg->AddEntry(hB,"tt+jets");
  leg->Draw();

  TH1F* hRatio = (TH1F*) hS->Clone("hRatio");
  hRatio->Reset();
  hRatio->SetFillColor(0);
  hRatio->Sumw2();
  hRatio->Divide(hS,hB,1,1,"B");
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);
  //hRatio->Draw("HISTE");
  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");

  //c1->SaveAs("plot_SL_CompMadWeight.png");
  //c1->SaveAs("ttbb_vs_ttjj.png");


}

void plot_syst(TString cat   = "SL",
	       TString syst  = "csv", 
	       string header = "type==0",
	       TString fname = "type0",
	       TString title = "",
	       int doRatio = 0,
	       TString sample = "v6_TTJetsSemiLept",
	       float fact = 0.02, float fact2 = 0.5
	       ){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 


  TFile* fN  = TFile::Open("MEAnalysis_"+cat+"_nominal_"+sample+".root");
  TFile* fU  = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Up_"+sample+".root");
  TFile* fD  = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Down_"+sample+".root");


  TTree* tN = (TTree*)fN->Get("tree");
  TTree* tU = (TTree*)fU->Get("tree");
  TTree* tD = (TTree*)fD->Get("tree");

  TH1F* hN = new TH1F("hN","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hN->SetLineColor(kBlack);
  hN->SetLineWidth(3);
  hN->Sumw2();
  TH1F* hU = new TH1F("hU","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hU->SetLineColor(kBlue);
  hU->SetLineWidth(3);
  hU->Sumw2();
  TH1F* hD = new TH1F("hD","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hD->SetLineColor(kRed);
  hD->SetLineWidth(3);
  hD->Sumw2();

  tN->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hN", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));
  tU->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hU", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));
  tD->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hD", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));

  hN->SetMaximum(  hN->GetMaximum()*1.8);
  hN->SetMinimum( 0 );

  hN->Draw("HISTE");
  hU->Draw("HISTESAME");
  hD->Draw("HISTESAME");
 

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hN,"nominal");
  leg->AddEntry(hU,syst+" up");
  leg->AddEntry(hD,syst+" down");
  leg->Draw();

  TH1F* hRatioN = (TH1F*) hN->Clone("hRatioN");
  hRatioN->Reset();
  hRatioN->SetFillColor(0);
  hRatioN->Sumw2();
  hRatioN->Divide(hN,hN,1,1,"B");
  TH1F* hRatioU = (TH1F*) hN->Clone("hRatioU");
  hRatioU->Reset();
  hRatioU->SetFillColor(0);
  hRatioU->Sumw2();
  hRatioU->Divide(hU,hN,1,1,"B");
  TH1F* hRatioD = (TH1F*) hN->Clone("hRatioD");
  hRatioD->Reset();
  hRatioD->SetFillColor(0);
  hRatioD->Sumw2();
  hRatioD->Divide(hD,hN,1,1,"B");

  if(doRatio){
    hRatioN->SetMinimum(0.4);
    hRatioN->SetMaximum(1.6);
    hRatioN->Draw("HISTE");
    hRatioU->Draw("HISTESAME");
    hRatioD->Draw("HISTESAME");
  }

  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");

}



void draw(vector<float> param, TTree* t = 0, TString var = "", TH1F* h = 0, TCut cut = ""){

  h->Reset();
  if(param.size()==0){
    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    return;
  }
  else{
    if( param.size()!=2 ){
      cout << "Error in draw()" << endl;
      return;
    }

    float cutval = h->GetBinLowEdge( 3 );
    TCut    cut1( Form("%s<%f", var.Data(), cutval) );    
    TString var2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))",    param[0]));

    TCut    cut2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))<%f", param[0],param[1]));    
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&cut2) );
    float binL    = h->Integral();
    float binLErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&(!cut2)) );
    float binH    = h->Integral();
    float binHErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();

    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    h->SetBinContent(1, binL);  h->SetBinError(1, binLErr);
    h->SetBinContent(2, binH);  h->SetBinError(2, binHErr);
  }

  return;
}


void plot_category(string cat = "SL_VType2_nominal", 
		   string  cut = "type==0",
		   float fact1 = 0.02,
		   float fact2 = 0.,
		   float lumiScale = 3.2,
		   string header = "",
		   string fname  = "",
		   int newvar = 1,
		   int nBins  = 6,
		   int splitFirstBin = 0
		   ){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.52,0.50,0.77,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  THStack* aStack = new THStack("aStack","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1};  L_{S}/(L_{S}+L_{B}) ; events ");
  //TH1F* hStack = (TH1F*)aStack->GetHistogram();
  //hStack->SetTitle(cutName.c_str());
  //hStack->SetXTitle(xTitle.c_str());
  //hStack->SetYTitle(yTitle.c_str());

 //  vector<string> samples;
//   samples.push_back("SingleT");
//   samples.push_back("TTV");
//   samples.push_back("EWK");
//   samples.push_back("DiBoson");
//   if( cat.find("VType2")!=string::npos || cat.find("VType3")!=string::npos ){
//     samples.push_back("TTJetsSemiLept");
//     samples.push_back("TTJetsSemiLept");
//   }
//   if( cat.find("VType0")!=string::npos || cat.find("VType1")!=string::npos ){
//     samples.push_back("TTJetsFullLept");
//     samples.push_back("TTJetsFullLept");
//   }
//   samples.push_back("TTH125");

  string version = cat.find("SL")!=string::npos ? "_v6" : "_v5";


  TString var("");
  if(newvar)
    var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, fact2, fact2));
  else
    var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));
  vector<float> param;
  param.clear();
  if(splitFirstBin!=0){
    param.push_back(1.0);
    param.push_back(0.5);
  }

  TArrayF bins (nBins+1);
  TArrayF bins2(nBins+2);
  cout << "Making histograms with " << nBins << " bins:" << endl;
  if( param.size()==0){
    for(int b = 0; b < nBins+1; b++)
      bins[b] = b*1.0/nBins;
  }
  else{
    for(int b = 0; b < nBins+2; b++){
      if(b<=2) bins2[b] = b*0.5/(nBins);
      else     bins2[b] = (b-1)*1.0/(nBins);
      cout <<  bins2[b] << ", ";
    }
    cout << endl;
  }



  vector<string> samples;
  if( cat.find("SL")!=string::npos ){
    samples.push_back("nominal"+version+"_SingleT");
    samples.push_back("nominal"+version+"_TTV");

    //samples.push_back("nominal"+version+"_EWK");
    //samples.push_back("nominal"+version+"_DiBoson");
    samples.push_back("nominal"+version+"_TTJetsSemiLept");
    samples.push_back("nominal"+version+"_TTJetsSemiLept");
    samples.push_back("nominal"+version+"_TTH125");

    //samples.push_back("VType2_nominal_wtag_TTJetsSemiLept");
    //samples.push_back("VType2_nominal_wtag_TTJetsSemiLept");
    //samples.push_back("VType2_nominal_wtag_TTH125");
  }
  if( cat.find("DL")!=string::npos ){
    samples.push_back("nominal"+version+"_SingleT");
    samples.push_back("nominal"+version+"_TTV");
    //samples.push_back("nominal"+version+"_EWK");
    //samples.push_back("nominal"+version+"_DiBoson");
    samples.push_back("nominal"+version+"_TTJetsFullLept");
    samples.push_back("nominal"+version+"_TTJetsFullLept");
    samples.push_back("nominal"+version+"_TTH125");
  }

  int countTTJets  = 0;
  int countSingleT = 0;
  int countTTV     = 0;
  int countEWK     = 0;
  int countDiBoson = 0;
  int countTTH     = 0;

  TH1F* hS = new TH1F("hS","Simulation 8 TeV",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
  hS->SetLineWidth(3);
  hS->SetLineStyle(kDashed);
  hS->SetLineColor(kRed);

  TH1F* hErr = new TH1F("hErr","Simulation 8 TeV L=19.5 fb^{-1}" , param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
  hErr->Sumw2();
  leg->AddEntry(hErr, "MC unc. (stat.)", "L");

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;
    TFile* f = TFile::Open(("MEAnalysis_"+cat+"_"+samples[sample]+".root").c_str());
    if(f==0 || f->IsZombie() ) continue;
    TTree* tree     = (TTree*)f->Get("tree");
    //TH1F*  hcounter = (TH1F*) f->Get("hcounter");

    TCut tcut(cut.c_str());

    //tcut = tcut && TCut("p_125_all_b_ttbb/((p_125_all_b_ttbb+2*p_125_all_b_ttjj))>0.18");

    if( samples[sample].find("TTJets") != string::npos && countTTJets==0){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs>2");
    }
    else if( samples[sample].find("TTJets") != string::npos && countTTJets==1){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs<=2");
    }

    TH1F* h = new TH1F(("h_"+samples[sample]).c_str(),"Simulation 8 TeV L=19.5 fb^{-1}; events; L_{S}/{L_{S}+L_{B}}",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
    h->Sumw2();
    draw( param, tree , var, h , tcut );

    //if(newvar)
    //tree->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>h_%s", fact1, fact2, fact2, samples[sample].c_str()), "weight"*tcut);
    //else
    //tree->Draw(Form("p_125_all_s/(p_125_all_s+%f*p_125_all_b)>>h_%s", fact1, samples[sample].c_str()), "weight"*tcut);
    //h->Scale(1./hcounter->GetBinContent(1) * lumiScale);

    h->Scale(lumiScale);
    cout << "Events = " << h->Integral() << endl;

    //for(int b = 1; b<= h->GetNbinsX(); b++ ){
      //h->SetBinContent( b, h->GetBinContent(b)/h->GetBinWidth(b) );
      //h->SetBinError  ( b, h->GetBinError  (b)/h->GetBinWidth(b) );
      //if( param.size()>0 ){
      //if(b==1)
      //h->GetXaxis()->SetBinLabel(b, "L_{LH}<0.5" );
      //else if(b==2) 
      //  h->GetXaxis()->SetBinLabel(b, "L_{LH}>0.5" );
      //else
      //  h->GetXaxis()->SetBinLabel(b, Form("%.1f",  h->GetBinCenter(b)) );
      //h->GetXaxis()->SetLabelSize(0.10);
      //h->GetXaxis()->SetNdivisions(nBins*100 + 2*nBins,kTRUE);	  
      //}
    //}


    if(samples[sample].find("TTH125") == string::npos){
      hErr->Add( h, 1.0);
    }
   
    if( samples[sample].find("TTH125") != string::npos ){
      h->SetLineColor( kRed );
      h->SetFillColor( kRed );
      if(countTTH==0) leg->AddEntry(h, "t#bar{t}H", "F");
      countTTH++;
      hS->Add( h, 5.);
    }
    if( samples[sample].find("TTJets") != string::npos ){
      if(countTTJets==1){
	leg->AddEntry(h, "t#bar{t}+jets HF", "F");
	h->SetLineColor( 17 );
	h->SetFillColor( 17 );
      }
      else if(countTTJets==2){
	h->SetLineColor( 18 );
	//h->SetLineWidth( 0 );
	leg->AddEntry(h, "t#bar{t}+jets LF", "F");
	h->SetFillColor( 18 );
      }
    }
    if( samples[sample].find("SingleT") != string::npos ){
      h->SetLineColor( kMagenta );
      h->SetFillColor( kMagenta );
      if(countSingleT==0) leg->AddEntry(h, "Single top", "F");
      countSingleT++;
    }
    if( samples[sample].find("EWK") != string::npos ){
      h->SetLineColor( kGreen );
      h->SetFillColor( kGreen );
      if(countEWK==0) leg->AddEntry(h, "V+jets", "F");
      countEWK++;
    }
    if( samples[sample].find("DiBoson") != string::npos ){
      h->SetLineColor( kYellow );
      h->SetFillColor( kYellow );
      if(countDiBoson==0) leg->AddEntry(h, "VV", "F");
      countDiBoson++;
    }
    if( samples[sample].find("TTV") != string::npos ){
      h->SetLineColor( 30 );
      h->SetFillColor( 30 );
      if(countTTV==0) leg->AddEntry(h, "t#bar{t}V", "F");
      countTTV++;
    }

    aStack->Add( h );
  }

  hErr->GetYaxis()->SetTitle("Events");
  hErr->GetXaxis()->SetTitle("L_{S} = P_{tth} / ( P_{tth} + P_{ttbb} + P_{ttjj} )");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  hErr->GetMaximum()*1.35;
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);
  hErr->Draw("HISTE1");
  aStack->Draw("HISTSAME");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 5", "L");
  hS->Draw("HISTSAME");
  hErr->Draw("HISTE1SAME");
  leg->Draw();

  TLine* line = new TLine(hErr->GetBinLowEdge(3), max , hErr->GetBinLowEdge(3), 0.);
  line->SetLineWidth(4);
  line->SetLineStyle(kSolid);
  line->SetLineColor(kBlack);
  if( param.size()>0 ) line->Draw("SAME");

  TPaveText *pt1 = new TPaveText(0.101, 0.839161, 0.198142, 0.895105,"brNDC");
  pt1->SetFillStyle(1001);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(kWhite);
  pt1->SetTextSize(0.03); 
  pt1->AddText("L_{bb}<0.5");

  TPaveText *pt2 = new TPaveText(0.191, 0.839161, 0.294118, 0.895105,"brNDC");
  pt2->SetFillStyle(1001);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(kWhite);
  pt2->SetTextSize(0.03); 
  pt2->AddText("L_{bb}>0.5");
 
  if( param.size()>0 ){
    pt1->Draw();
    pt2->Draw();
  }

  cout << "Signal = " << hS->Integral()/5. << endl;

  c1->SaveAs(  (fname+version+".png").c_str() );

}



void plotAll(){


  
  plot_category( "SL", 
		 "type==0",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat1 (4b2j W-tag)",
		 "Plot_SL_Cat1",
		 1,
		 4,
		 1
		 );
  
  
  plot_category( "SL", 
		 "type==1",
		 0.018, 0.70, 
		 19.5/12.1,
		 "Cat2 (4b2j !W-tag)",
		 "Plot_SL_Cat2",
		 1,
		 5,
		 1
		 );


  plot_category( "SL", 
		 "type==2 && flag_type2>0",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat3 (4b1j W-tag)",
		 "Plot_SL_Cat3",
		 1,
		 5,
		 1
		 );
  

  plot_category( "SL", 
		 "type==2 && flag_type2<=0",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat4 (4b1j !W-tag)",
		 "Plot_SL_Cat4",
		 1,
		 5,
		 1
		 );


  plot_category( "SL", 
		 "type==3",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat5 (4b3j)",
		 "Plot_SL_Cat5",
		 1,
		 5,
		 1
		 );

  return;

  
  plot_category( "DL", 
		 "type==6",
		 0.02, 0,
		 19.5/12.1 * 2,
		 "Cat6 (4b), tight",
		 "Plot_DL_Cat6",
		 0,
		 4,
		 0
		 );

  plot_category( "DL", 
		 "type==7",
		 0.02, 0,
		 19.5/12.1 * 2,
		 "Cat7 (4b), loose",
		 "Plot_DL_Cat7",
		 0,
		 6,
		 0
		 );
  
  
  

}



