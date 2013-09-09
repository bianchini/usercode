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

#include "TVector3.h"
#include "TLorentzVector.h"






void plot_analysis(TString cat   ="SL_VType2_nominal", 
		   string header="type==0",
		   TString title = "",
		   TString fname= "type0",
		   float fact = 1.5
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
  TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_TTJetsSemiLept.root");
  TFile* fB = TFile::Open("MEAnalysis_"+cat+"_JECUp_TTJetsSemiLept.root");

  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",6, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",6, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();

  //tS->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hS", fact));
  //tB->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hB", fact));

  tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>hS", fact), TCut(header.c_str()));
  tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>hB", fact), TCut(header.c_str()));

  if(hS->GetMaximum()>=hB->GetMaximum()){
    //hS->Scale( 1./hS->Integral());
    //hB->Scale( 1./hB->Integral());
    //hS->SetMaximum(0.5);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    //hS->Scale( 1./hS->Integral());
    //hB->Scale( 1./hB->Integral());
    //hB->SetMaximum(0.5);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hS,"Signal (ttH)");
  leg->AddEntry(hB,"Background (tt+jets)");
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

}

void plot_syst(TString cat   = "SL_VType2",
	       TString syst  = "csv", 
	       string header = "type==0",
	       int doRatio = 0,
	       TString sample = "TTJetsSemiLept",
	       float fact = 1.5
	       ){

  TString title = "";
  TString fname = "type0";

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
  TFile* fU = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Up_"+sample+".root");
  TFile* fD = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Down_"+sample+".root");


  TTree* tN = (TTree*)fN->Get("tree");
  TTree* tU = (TTree*)fU->Get("tree");
  TTree* tD = (TTree*)fD->Get("tree");

  TH1F* hN = new TH1F("hN","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",6, 0,1);
  hN->SetLineColor(kBlack);
  hN->SetLineWidth(3);
  hN->Sumw2();
  TH1F* hU = new TH1F("hU","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",6, 0,1);
  hU->SetLineColor(kBlue);
  hU->SetLineWidth(3);
  hU->Sumw2();
  TH1F* hD = new TH1F("hD","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",6, 0,1);
  hD->SetLineColor(kRed);
  hD->SetLineWidth(3);
  hD->Sumw2();

  //tS->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hN", fact));
  //tB->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hD", fact));

  tN->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>hN", fact), TCut(header.c_str()));
  tU->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>hU", fact), TCut(header.c_str()));
  tD->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>hD", fact), TCut(header.c_str()));

  hN->SetMaximum(  hN->GetMaximum()*1.8);
  hN->SetMinimum( 0 );

  hN->Draw("HISTE");
  hU->Draw("HISTESAME");
  hD->Draw("HISTESAME");
 

  leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
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




void plot_category(string cat = "SL_VType2_nominal", 
		   string  cut = "type==0",
		   float fact = 1.5,
		   float lumiScale = 3.2,
		   string header = "",
		   string fname  = ""
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

  TH1F* hS = new TH1F("hS","Simulation 8 TeV", 6,0,1);
  hS->SetLineWidth(3);
  hS->SetLineStyle(kDashed);
  hS->SetLineColor(kRed);

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


  vector<string> samples;
  if( cat.find("SL")!=string::npos ){
    samples.push_back("nominal_SingleT");
    samples.push_back("nominal_TTV");
    samples.push_back("nominal_EWK");
    samples.push_back("nominal_DiBoson");
    samples.push_back("nominal_TTJetsSemiLept");
    samples.push_back("nominal_TTJetsSemiLept");
    samples.push_back("nominal_TTH125");
  }
  if( cat.find("DL")!=string::npos ){
    samples.push_back("nominal_SingleT");
    samples.push_back("nominal_TTV");
    samples.push_back("nominal_EWK");
    samples.push_back("nominal_DiBoson");
    samples.push_back("nominal_TTJetsFullLept");
    samples.push_back("nominal_TTJetsFullLept");
    samples.push_back("nominal_TTH125");
  }

  int countTTJets  = 0;
  int countSingleT = 0;
  int countTTV     = 0;
  int countEWK     = 0;
  int countDiBoson = 0;
  int countTTH     = 0;

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;
    TFile* f = TFile::Open(("MEAnalysis_"+cat+"_"+samples[sample]+".root").c_str());
    if(f==0 || f->IsZombie() ) continue;
    TTree* tree     = (TTree*)f->Get("tree");
    //TH1F*  hcounter = (TH1F*) f->Get("hcounter");

    TCut tcut(cut.c_str());

    if( samples[sample].find("TTJets") != string::npos && countTTJets==0){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs>2");
    }
    else if( samples[sample].find("TTJets") != string::npos && countTTJets==1){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs<=2");
    }

    TH1F* h = new TH1F(("h_"+samples[sample]).c_str(),"Simulation 8 TeV L=19.5 fb^{-1}; events; L_{S}/{L_{S}+L_{B}}", 6,0,1);
    tree->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)>>h_%s", fact, samples[sample].c_str()), "weight"*tcut);

    //h->Scale(1./hcounter->GetBinContent(1) * lumiScale);

    h->Scale(lumiScale);
    cout << "Events = " << h->Integral() << endl;

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
	h->SetLineColor( 19 );
	//h->SetLineWidth( 0 );
	leg->AddEntry(h, "t#bar{t}+jets LF", "F");
	h->SetFillColor( 19 );
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


  aStack->Draw("HIST");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 5", "L");
  hS->Draw("HISTSAME");
  leg->Draw();

  cout << "Signal = " << hS->Integral()/5. << endl;

  c1->SaveAs(  (fname+".png").c_str() );

}



void plotAll(){

  plot_category( "SL", 
		 "type==0",
		 0.5,
		 19.5/12.1,
		 "1l + 4b + 2j, W-tag",
		 "Plot_SL_type0"
		 );

  plot_category( "SL", 
		 "type==1",
		 0.5,
		 19.5/12.1,
		 "1l + 4b + 2j, !W-tag",
		 "Plot_SL_type1"
		 );

  plot_category( "SL", 
		 "type==2 && flag_type2>0",
		 0.5,
		 19.5/12.1,
		 "1l + 4b + 1j, W-tag",
		 "Plot_SL_type2_1"
		 );

  plot_category( "SL", 
		 "type==2 && flag_type2<0",
		 0.5,
		 19.5/12.1,
		 "1l + 4b + 1j, !W-tag",
		 "Plot_SL_type2_2"
		 );

  plot_category( "SL", 
		 "type==3",
		 0.5,
		 19.5/12.1,
		 "1l + 4b + 3j",
		 "Plot_SL_type3"
		 );

  plot_category( "DL", 
		 "type==6",
		 0.5,
		 19.5/12.1 * 2,
		 "2l + 4b, tight",
		 "Plot_SL_type6"
		 );

  plot_category( "DL", 
		 "type==7",
		 0.5,
		 19.5/12.1 * 2,
		 "2l + 4b, loose",
		 "Plot_SL_type7"
		 );


}
