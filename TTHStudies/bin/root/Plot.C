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

#define SAVEPLOTS 1


void plot_BestMass(TString fname  = "gen_default",
		   TString header = "Signal (parton-level)",
		   TString hname  = "hBestMass",
		   int color = 2
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

  TFile *f = TFile::Open("TestMENew_"+fname+".root","READ");


  TH1F* hBestMass = (TH1F*)f->Get(hname);
  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }

  //hBestMass->Rebin(3);

  hBestMass->SetFillStyle(3004);
  hBestMass->SetFillColor(color);
  hBestMass->SetLineWidth(2);
  hBestMass->SetLineColor(color);
  hBestMass->SetTitle("");
  hBestMass->SetXTitle("M (GeV)");
  hBestMass->SetYTitle("P(M)");
  
  TLegend* leg = new TLegend(0.48,0.65,0.78,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  string str_header(hname.Data());
  if( str_header.find("Prob")!=string::npos)
    leg->AddEntry(hBestMass, Form("#splitline{All permutations}{#mu=%.1f, RMS=%.1f}",  
				  hBestMass->GetMean(), hBestMass->GetRMS() ),"F");
  else
    leg->AddEntry(hBestMass, Form("#splitline{Good permutation}{#mu=%.1f, RMS=%.1f}",  
				  hBestMass->GetMean(), hBestMass->GetRMS() ),"F");


  hBestMass->DrawNormalized("HIST");
  leg->Draw();

  if(SAVEPLOTS){
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+".png");
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+".pdf");
  }

  return;
}




void plot_TF1d(TString fname  = "gen_default",
	       TString header = "Light jet TF",
	       TString unitsX = "jet E (GeV)",
	       TString unitsY = "TF(E|#hat{E})",
	       TString hname  = "tfWjet1",
	       int event = 1,
	       TString extraLabel = ""){

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

  TFile *f = TFile::Open("TestMENew_"+fname+".root","READ");
  TH1F* hBestMass = (TH1F*)f->Get(TString(Form("Event_%d/",event))+hname);


  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }

  hBestMass->SetLineWidth(2);
  hBestMass->SetLineColor(kRed);
  hBestMass->SetTitle("");
  hBestMass->SetXTitle( unitsX );
  hBestMass->SetYTitle( "units" );
  hBestMass->SetTitleSize(0.04,"X");
  hBestMass->SetTitleSize(0.04,"Y");

  TLegend* leg = new TLegend(0.65,0.65,0.80,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 
  leg->AddEntry(hBestMass, unitsY);

  hBestMass->Draw("HIST");
  leg->Draw();

  if(SAVEPLOTS){
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+extraLabel+".png");
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+extraLabel+".pdf");
  }


  return;
}



void plot_TF2d(TString fname  = "gen_default",
	       TString header = "MEt TF",
	       TString unitsX = "p_{T} (GeV)",
	       TString unitsY = "#phi",
	       TString unitsZ = "TF(#phi|#hat{#phi},p_{T})",
	       TString hname  = "tfMetPhi",
	       int event = 1
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

  TFile *f = TFile::Open("TestMENew_"+fname+".root","READ");
  TH2F* hBestMass = (TH2F*)f->Get(TString(Form("Event_%d/",event))+hname);

  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }

  hBestMass->SetLineWidth(2);
  hBestMass->SetLineColor(kRed);
  hBestMass->SetTitle("");
  hBestMass->SetXTitle( unitsX );
  hBestMass->SetYTitle( unitsY );

  TLegend* leg = new TLegend(0.65,0.65,0.80,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  leg->AddEntry(hBestMass, unitsZ);

  hBestMass->Draw("LEGO");
  leg->Draw();

  if(SAVEPLOTS){
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+".png");
    c1->SaveAs("plots/Plot_"+hname+"_"+fname+".pdf");
  }


  return;
}



void plot_param(TString header = "Light jet TF",	      
		TString unitsX = "jet E (GeV)",
		TString unitsY = "Acceptance",
		TString hname  = "accLightBin0",
		float xLow  = 20,
		float xHigh = 100,
		float yLow  = -999,
		float yHigh = -999,
		TString extraLabel = ""
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

  TFile *f = TFile::Open("ControlPlots.root","READ");
  TH1F* hBestMass = (TH1F*)f->Get(hname);


  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }


  TString fname;
  TIter iter( hBestMass->GetListOfFunctions() );
  TObject* obj = iter();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    //cout << func->GetName() << endl;
    fname = func->GetName();
    obj = iter();
  }
  TF1 *func = hBestMass->GetFunction(fname);


  hBestMass->SetMarkerStyle(kFullCircle);
  hBestMass->SetMarkerColor(kBlue);
  hBestMass->SetTitle("");
  hBestMass->SetXTitle( unitsX );
  hBestMass->SetYTitle( unitsY );
  hBestMass->SetTitleSize(0.05,"X");
  hBestMass->SetTitleSize(0.05,"Y");

  hBestMass->SetTitleOffset(0.90,"X");
  hBestMass->SetTitleOffset(0.83,"Y");

  hBestMass->GetXaxis()->SetRangeUser(xLow,xHigh);
  if(yLow!=-999 && yHigh!=-999)
    hBestMass->GetYaxis()->SetRangeUser(yLow,yHigh);

  TLegend* leg = new TLegend(0.15,0.65,0.55,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 
  leg->AddEntry(func, "Fit", "L");

  hBestMass->Draw("PE");
  leg->Draw();

  if(SAVEPLOTS){
    c1->SaveAs("plots/Plot_"+hname+"_tf_"+extraLabel+".png");
    c1->SaveAs("plots/Plot_"+hname+"_tf_"+extraLabel+".pdf");
  }


  return;
}


void plot_genreco(string mode = "SL2wj", string norm = "acc",
		  string mass = "bestInt",
		  string cat  = "SL(4,2)",
		  int  nBins  = 20, 
		  float xLow  = 50,
		  float xHigh = 250){

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


  TLegend* leg = new TLegend(0.52,0.60,0.82,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  TPaveText *pt = new TPaveText(.52,.22,.82,.35,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04); 


  TFile *fgen = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+".root").c_str(),"READ");
  TFile *frec = TFile::Open(("MEValidator_"+mode+"_rec_"+norm+".root").c_str(),"READ");

  TTree* tgen = (TTree*)fgen->Get("tree");
  TTree* trec = (TTree*)frec->Get("tree");

  TH1F* hgen_good = new TH1F("hgen_good","; M (GeV); units", nBins, xLow, xHigh);
  hgen_good->SetLineColor(kRed);
  hgen_good->SetLineStyle(kSolid);
  hgen_good->SetLineWidth(3);
  hgen_good->SetFillStyle(3002);
  hgen_good->SetFillColor(kRed);
  TH1F* hgen_bad  = new TH1F("hgen_bad", "; M (GeV); units", nBins, xLow, xHigh);
  hgen_bad->SetLineColor(kRed+2);
  hgen_bad->SetLineStyle(kDashed);
  hgen_bad->SetLineWidth(2);
  hgen_bad->SetFillStyle(3002);
  hgen_bad->SetFillColor(kRed+2);

  tgen->Draw( ("m_"+mass+"_all_gen>>hgen_good").c_str(),  ("match_"+mass+"_all_gen==1").c_str());
  tgen->Draw( ("m_"+mass+"_all_gen>>hgen_bad" ).c_str(),  ("match_"+mass+"_all_gen==0").c_str());
  float normgen = hgen_good->Integral() + hgen_bad->Integral();
  hgen_good->Scale(1./normgen);
  hgen_bad ->Scale(1./normgen);

  THStack*  sgen =  new THStack("gen","");
  sgen->Add(hgen_bad);
  sgen->Add(hgen_good);

  TH1F* hrec_good = new TH1F("hrec_good","; M (GeV); units", nBins, xLow, xHigh);
  hrec_good->SetLineColor(kBlue);
  hrec_good->SetLineStyle(kSolid);
  hrec_good->SetLineWidth(3);
  hrec_good->SetFillStyle(3354);
  hrec_good->SetFillColor(kBlue);
  TH1F* hrec_bad  = new TH1F("hrec_bad", "; M (GeV); units", nBins, xLow, xHigh);
  hrec_bad->SetLineColor(kBlue+2);
  hrec_bad->SetLineStyle(kDashed);
  hrec_bad->SetLineWidth(2);
  hrec_bad->SetFillStyle(3354);
  hrec_bad->SetFillColor(kBlue+2);

  trec->Draw( ("m_"+mass+"_all_rec>>hrec_good").c_str(),  ("match_"+mass+"_all_rec==1").c_str());
  trec->Draw( ("m_"+mass+"_all_rec>>hrec_bad" ).c_str(),  ("match_"+mass+"_all_rec==0").c_str());
  float normrec = hrec_good->Integral() + hrec_bad->Integral();
  hrec_good->Scale(1./normrec);
  hrec_bad ->Scale(1./normrec);


  THStack*  srec =  new THStack("rec","");
  srec->Add(hrec_bad);
  srec->Add(hrec_good);

  sgen->Draw("HIST");
  srec->Draw("HISTSAME");
  hgen_bad->Draw("HISTSAME");
  
  TH1F* hStack = (TH1F*)sgen->GetHistogram();
  hStack->SetXTitle("M (GeV)");
  hStack->SetYTitle("units");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleSize(0.05,"X");
  hStack->SetTitleOffset(0.92,"X");
  hStack->SetTitleOffset(0.98,"Y");
  hStack->SetTitle(Form("Match efficiency: %.0f%% (partons), %.0f%% (smeared)", hgen_good->Integral()*100, hrec_good->Integral()*100  ));

  leg->SetHeader(("t#bar{t}H(125), "+cat).c_str());
  leg->AddEntry(hgen_good, "partons, good" ,  "F");
  leg->AddEntry(hgen_bad,  "partons, wrong",  "F");
  leg->AddEntry(hrec_good, "smeared, good",   "F");
  leg->AddEntry(hrec_bad,  "smeared, wrong",  "F");

  leg->Draw();
  //sgen->Draw("HISTSAME");

  pt->AddText("Good matching:");
  pt->AddText(Form("partons: %.0f%%",hgen_good->Integral()*100 ))->SetTextColor(kRed);
  pt->AddText(Form("smeared: %.0f%%",hrec_good->Integral()*100 ))->SetTextColor(kBlue);
  pt->SetTextAlign(31);

  pt->Draw();

  if(SAVEPLOTS){
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass+".png").c_str());
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass+".pdf").c_str());
  }


  return;
}



void plot_masses(string mode = "SL2wj", 
		 string norm = "acc",
		 string level = "gen",
		 string mass1 = "bestInt",  string mass2 = "bestMax",
		 string cat  = "SL(4,2)",
		 int  nBins  = 20, 
		 float xLow  = 50,
		 float xHigh = 250){

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


  TLegend* leg = new TLegend(0.50,0.60,0.80,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  TPaveText *pt = new TPaveText(.52,.22,.82,.35,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04); 


  TFile *fgen = TFile::Open(("MEValidator_"+mode+"_"+level+"_"+norm+".root").c_str(),"READ");
  TFile *frec = TFile::Open(("MEValidator_"+mode+"_"+level+"_"+norm+".root").c_str(),"READ");

  TTree* tgen = (TTree*)fgen->Get("tree");
  TTree* trec = (TTree*)frec->Get("tree");

  TH1F* hgen_good = new TH1F("hgen_good","; M (GeV); units", nBins, xLow, xHigh);
  hgen_good->SetLineColor(kRed);
  hgen_good->SetLineStyle(kSolid);
  hgen_good->SetLineWidth(3);
  hgen_good->SetFillStyle(3002);
  hgen_good->SetFillColor(kRed);
 
  tgen->Draw( ("m_"+mass1+"_all_"+level+">>hgen_good").c_str(),  ("m_"+mass1+"_all_"+level+">0").c_str());

  float normgen = hgen_good->Integral() ;
  hgen_good->Scale(1./normgen);

  THStack*  sgen =  new THStack("gen","");
  sgen->Add(hgen_good);

  TH1F* hrec_good = new TH1F("hrec_good","; M (GeV); units", nBins, xLow, xHigh);
  hrec_good->SetLineColor(kBlue);
  hrec_good->SetLineStyle(kSolid);
  hrec_good->SetLineWidth(3);
  hrec_good->SetFillStyle(3354);
  hrec_good->SetFillColor(kBlue);

  trec->Draw( ("m_"+mass2+"_all_"+level+">>hrec_good").c_str(),  ("m_"+mass2+"_all_"+level+">0").c_str());

  float normrec = hrec_good->Integral() ;
  hrec_good->Scale(1./normrec);


  THStack*  srec =  new THStack("rec","");
  srec->Add(hrec_good);

  sgen->Draw("HIST");
  srec->Draw("HISTSAME");
   
  TH1F* hStack = (TH1F*)sgen->GetHistogram();
  hStack->SetXTitle("M (GeV)");
  hStack->SetYTitle("units");
  hStack->SetTitleSize(0.05,"Y");
  hStack->SetTitleSize(0.05,"X");
  hStack->SetTitleOffset(0.92,"X");
  hStack->SetTitleOffset(0.98,"Y");
  hStack->SetTitle(Form("Match efficiency: %.0f%% (partons), %.0f%% (smeared)", hgen_good->Integral()*100, hrec_good->Integral()*100  ));

  string what = "dummy";
  if(level.find("rec")!=string::npos) what = "smeared";
  if(level.find("gen")!=string::npos) what = "partons";

  string mass_1 = "dummy";
  if(mass1.find("Int")!=string::npos) 
    mass_1 = "max_{p}#left{#intP(M|Y_{p})dM#right}";
  else if(mass1.find("Max")!=string::npos) 
    mass_1 = "max_{p,M}#left{P(M|Y_{p})#right}";
  else
    mass_1 = "max_{M}#left{#sumP(M|Y_{p})#right}";

 
  string mass_2 = "dummy";
  if(mass2.find("Int")!=string::npos) 
    mass_2 = "max_{p}#left{#intP(M|Y_{p})dM#right}";
  else if(mass2.find("Max")!=string::npos) 
    mass_2 = "max_{p,M}#left{P(M|Y_{p})#right}";
  else
    mass_2 = "max_{M}#left{#sumP(M|Y_{p})#right}";



  leg->SetHeader(("t#bar{t}H(125), "+cat).c_str());
  leg->AddEntry(hgen_good, (mass_1).c_str() ,  "F");
  leg->AddEntry(hrec_good, (mass_2).c_str(),   "F");

  leg->Draw();
  //sgen->Draw("HISTSAME");

  pt->AddText("Good matching:");
  pt->AddText(Form("partons: %.0f%%",hgen_good->Integral()*100 ))->SetTextColor(kRed);
  pt->AddText(Form("smeared: %.0f%%",hrec_good->Integral()*100 ))->SetTextColor(kBlue);
  pt->SetTextAlign(31);

  //pt->Draw();

  if(SAVEPLOTS){
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass1+"_vs_"+mass2+".png").c_str());
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass1+"_vs_"+mass2+".pdf").c_str());
  }


  return;
}




void plot_systematics(int syst = 0, string mode = "SL2wj", string norm = "acc",
		      string mass = "bestInt",
		      string cat  = "SL(4,2)",
		      int  nBins  = 20, 
		      float xLow  = 50,
		      float xHigh = 250){

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


  TLegend* leg = new TLegend(0.50,0.60,0.80,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->SetHeader(("t#bar{t}H(125), "+cat).c_str());

  TFile *fgen_nominal = syst<6 ? TFile::Open(("MEValidator_"+mode+"_gen_"+norm+".root").c_str(),"READ") : TFile::Open(("MEValidator_"+mode+"_rec_"+norm+".root").c_str(),"READ");
  TFile *fgen_noME    = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+"_noME.root").c_str(),"READ");
  TFile *fgen_noMET   = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+"_noMET.root").c_str(),"READ");
  TFile *fgen_noJac   = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+"_noJac.root").c_str(),"READ");
  TFile *fgen_noTF    = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+"_noTF.root").c_str(),"READ");
  TFile *fgen_noPDF   = TFile::Open(("MEValidator_"+mode+"_gen_"+norm+"_noPDF.root").c_str(),"READ");
  TFile *fgen_jUp     = TFile::Open(("MEValidator_"+mode+"_rec_"+norm+"_jUp.root").c_str(),"READ");
  TFile *fgen_bUp     = TFile::Open(("MEValidator_"+mode+"_rec_"+norm+"_bUp.root").c_str(),"READ");
  TFile *fgen_METUp   = TFile::Open(("MEValidator_"+mode+"_rec_"+norm+"_METUp.root").c_str(),"READ");

  vector<TFile*> files;
  files.push_back(fgen_nominal);
  files.push_back(fgen_noME);
  files.push_back(fgen_noMET);
  files.push_back(fgen_noJac);
  files.push_back(fgen_noTF);
  files.push_back(fgen_noPDF);
  files.push_back(fgen_jUp);
  files.push_back(fgen_bUp);
  files.push_back(fgen_METUp);

  vector<TH1F*> histos;

  TH1F* h = new TH1F("h","; M (GeV); units", nBins, xLow, xHigh);
  h->SetTitleSize(0.05,"Y");
  h->SetTitleSize(0.05,"X");
  h->SetTitleOffset(0.92,"X");
  h->SetTitleOffset(0.98,"Y");

  for(unsigned int k = 0; k < files.size(); k++){
    TTree* t = (TTree*)files[k]->Get("tree");
    TH1F*  h2 = (TH1F*)h->Clone(Form("h_%d",k));
    if(k<6 && syst<6)
      t->Draw( Form("m_%s_all_gen>>h_%d", mass.c_str(), k) ,  ("m_"+mass+"_all_gen>0").c_str());
    else if( k<6 && syst>5)
      t->Draw( Form("m_%s_all_rec>>h_%d", mass.c_str(), k) ,  ("m_"+mass+"_all_rec>0").c_str());
    else if(k>5 && syst>5)
      t->Draw( Form("m_%s_all_rec>>h_%d", mass.c_str(), k) ,  ("m_"+mass+"_all_rec>0").c_str());
    else{}

    h2->Scale(1./h2->Integral());
    if(k==0){
      h2->SetLineColor(kRed);
      h2->SetLineStyle(kSolid);
      h2->SetLineWidth(3);
      h2->SetFillStyle(3354);
      h2->SetFillColor(kRed);
    }
    else{
      h2->SetLineColor(kBlue);
      h2->SetLineStyle(kDashed);
      h2->SetLineWidth(3);
      h2->SetFillStyle(0);
      h2->SetFillColor(0);
    }
    histos.push_back( h2 );
  }

  float max = 0;
  for( unsigned int k = 0; k < histos.size(); k++ ){
    if( histos[k]->GetMaximum()>max) max =  histos[k]->GetMaximum();
  }

  string name = "dummy";

  for(unsigned int k = 0; k < histos.size(); k++){
    histos[k]->GetYaxis()->SetRangeUser(0, max*1.2);

    if(k==0) 
      histos[k]->Draw("HIST");
    else if(k>0 && k==syst)
      histos[k]->Draw("HISTSAME");
    else{};

    if( (string(files[k]->GetName())).find("noME.")!=string::npos )
      name = "w/o ME";
    else if( (string(files[k]->GetName())).find("noMET.")!=string::npos )
      name = "w/o E_{T}^{miss} TF";
    else if( (string(files[k]->GetName())).find("noTF.")!=string::npos )
      name = "w/o jet TF";
    else if( (string(files[k]->GetName())).find("noJac.")!=string::npos )
      name = "w/o Jacobians";
    else if( (string(files[k]->GetName())).find("noPDF.")!=string::npos )
      name = "w/o PDF";
    else if( (string(files[k]->GetName())).find("jUp.")!=string::npos )
      name = "1.25 #times #sigma_{j}";
    else if( (string(files[k]->GetName())).find("bUp.")!=string::npos )
      name = "1.25 #times #sigma_{b}";
    else if( (string(files[k]->GetName())).find("METUp.")!=string::npos )
      name = "1.25 #times #sigma_{MET}";
    else
      name = "Nominal";

    if(k==0 || (k>0 && k==syst)) leg->AddEntry( histos[k], name.c_str(), "F" );
  }
 
  leg->Draw();

  if(SAVEPLOTS){
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass+"_"+string(Form("%d",syst))+".png").c_str());
    c1->SaveAs(("plots/Plot_"+mode+"_"+norm+"_"+mass+"_"+string(Form("%d",syst))+".pdf").c_str());
  }


  return;

}


void makeAll(){

  plot_genreco("SL2wj",   "acc","bestInt","SL(4,2)",           25,50,250);
  plot_genreco("SL1wj",   "acc","bestInt","SL(4,2)",           25,50,250);
  plot_genreco("SLNoBHad","acc","bestInt","SL(3,2) no b_{q}",  25,50,250);
  plot_genreco("SLNoBLep","acc","bestInt","SL(3,2) no b_{l}",  25,50,250);
  plot_genreco("DL",      "acc","bestInt","DL(4,N)",           25,50,250);
  
  plot_genreco("SL2wj",   "unnorm","bestInt","SL(4,2)",           25,50,250);
  plot_genreco("SL1wj",   "unnorm","bestInt","SL(4,2)",           25,50,250);
  plot_genreco("SLNoBHad","unnorm","bestInt","SL(3,2) no b_{q}",  25,50,250);
  plot_genreco("SLNoBLep","unnorm","bestInt","SL(3,2) no b_{l}",  25,50,250);
  plot_genreco("DL",      "unnorm","bestInt","DL(4,N)",           25,50,250);

  plot_systematics(1, "SL2wj",  "acc", "bestInt", "SL(4,2) partons", 50, 50, 250);
  plot_systematics(2, "SL2wj",  "acc", "bestInt", "SL(4,2) partons", 50, 50, 250);
  plot_systematics(3, "SL2wj",  "acc", "bestInt", "SL(4,2) partons", 50, 50, 250);
  plot_systematics(4, "SL2wj",  "acc", "bestInt", "SL(4,2) partons", 50, 50, 250);
  plot_systematics(5, "SL2wj",  "acc", "bestInt", "SL(4,2) partons", 50, 50, 250);
  plot_systematics(6, "SL2wj",  "acc", "bestInt", "SL(4,2) smeared", 15, 50, 250);
  plot_systematics(7, "SL2wj",  "acc", "bestInt", "SL(4,2) smeared", 15, 50, 250);
  plot_systematics(8, "SL2wj",  "acc", "bestInt", "SL(4,2) smeared", 15, 50, 250);

  plot_masses("SL2wj","acc","gen","bestInt","best","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","acc","gen","bestInt","best","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","acc","gen","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","acc","gen","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("DL","acc","gen","bestInt","best","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","unnorm","gen","bestInt","best","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","unnorm","gen","bestInt","best","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","unnorm","gen","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","unnorm","gen","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("DL","unnorm","gen","bestInt","best","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","acc","rec","bestInt","best","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","acc","rec","bestInt","best","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","acc","rec","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","acc","rec","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("DL","acc","rec","bestInt","best","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","unnorm","rec","bestInt","best","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","unnorm","rec","bestInt","best","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","unnorm","rec","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","unnorm","rec","bestInt","best","SL(3,2) partons",20,50,250);
  plot_masses("DL","unnorm","rec","bestInt","best","DL(4,N) partons",20,50,250);


  plot_masses("SL2wj","acc","gen","bestInt","bestMax","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","acc","gen","bestInt","bestMax","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","acc","gen","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","acc","gen","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("DL","acc","gen","bestInt","bestMax","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","unnorm","gen","bestInt","bestMax","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","unnorm","gen","bestInt","bestMax","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","unnorm","gen","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","unnorm","gen","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("DL","unnorm","gen","bestInt","bestMax","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","acc","rec","bestInt","bestMax","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","acc","rec","bestInt","bestMax","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","acc","rec","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","acc","rec","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("DL","acc","rec","bestInt","bestMax","DL(4,N) partons",20,50,250);

  plot_masses("SL2wj","unnorm","rec","bestInt","bestMax","SL(4,2) partons",20,50,250);
  plot_masses("SL1wj","unnorm","rec","bestInt","bestMax","SL(4,1) partons",20,50,250);
  plot_masses("SLNoBHad","unnorm","rec","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("SLNoBLep","unnorm","rec","bestInt","bestMax","SL(3,2) partons",20,50,250);
  plot_masses("DL","unnorm","rec","bestInt","bestMax","DL(4,N) partons",20,50,250);


  //plot_param("udcsg 0<|#eta|<1.5", "parton E (GeV)",   "#sigma^{j}(E)", "resolLightBin0", 30,  180, 0, 30, "");
  //plot_param("udcsg 1.5<|#eta|<2.5", "parton E (GeV)", "#sigma^{j}(E)", "resolLightBin1", 30,  180, 0, 30, "");
  //plot_param("b-quarks 0<|#eta|<1.5", "parton E (GeV)",   "#sigma^{j}(E)", "resolHeavyBin0", 30,  180, 0, 30, "");
  //plot_param("b-quarks 1.5<|#eta|<2.5", "parton E (GeV)", "#sigma^{j}(E)", "resolHeavyBin1", 30,  180, 0, 30, "");

  //plot_param("udcsg 0<|#eta|<1.5",   "parton E (GeV)", "#mu^{j}(E)", "respLightBin0", 30,  200, 30, 200, "");
  //plot_param("udcsg 1.5<|#eta|<2.5", "parton E (GeV)", "#mu^{j}(E)", "respLightBin1", 30,  200, 30, 180, "");
  //plot_param("b-quarks 0<|#eta|<1.5",   "parton E (GeV)", "#mu^{j}(E)", "respHeavyBin0", 30,  200, 30, 200, "");
  //plot_param("b-quarks 1.5<|#eta|<2.5", "parton E (GeV)", "#mu^{j}(E)", "respHeavyBin1", 30,  200, 30, 200, "");

  //plot_param("E^{miss}_{x,y}", "#sumE_{T} (GeV)",   "#sigma^{MET}_{x,y}", "hWidthsPxResol", 0,  8000, 0, 50, "");

  //plot_param("MET, #sum E_T<1200 (GeV)", "#nu p_{T} (GeV)", "#mu^{#nu_{T}}", "hMeanEt_0", 10,  299, -20, 30, "");
  //plot_param("MET, 1200<#sum E_T<1800 (GeV)", "#nu p_{T} (GeV)", "#mu^{#nu_{T}}", "hMeanEt_1", 10,  299, -20, 30, "");
  //plot_param("MET, 1800<#sum E_T (GeV)", "#nu p_{T} (GeV)", "#mu^{#nu_{T}}", "hMeanEt_2", 10,  299, -20, 30, "");

  //plot_param("MET, #sum E_T<1200 (GeV)", "#nu p_{T} (GeV)", "#sigma^{#nu_{T}}", "hWidthEt_0", 10,  299, 0, 60, "");
  //plot_param("MET, 1200<#sum E_T<1800 (GeV)", "#nu p_{T} (GeV)", "#sigma^{#nu_{T}}", "hWidthEt_1", 10,  299, 0, 60, "");
  //plot_param("MET, 1800<#sum E_T (GeV)", "#nu p_{T} (GeV)", "#sigma^{#nu_{T}}", "hWidthEt_2", 10,  299, 0, 60, "");

  return;

  /*

  plot_BestMass("rec_default_ttjets_largerLR_range", "t#bar{t}+jets", "hBestMassProb", 4);
  */

  //plot_TF1d("gen_default","Light jet TF", "parton E (GeV)","TF(E|#hat{E})","tfWjet1");
  //plot_TF1d("gen_default","Heavy jet TF", "parton E (GeV)","TF(E|#hat{E})","tfbHad");
  //plot_TF1d("gen_default","MET TF", "#nu p_{T} (GeV)","TF(E|#hat{E})","tfMetPt");
  //plot_TF2d("gen_default","MET TF", "#nu p_{T} (GeV)","#phi", "TF(#phi|#hat{#phi},p_{T})","tfMetPhi");

  /*
  plot_TF1d("rec_default","P(M|Y_{good})", "M (GeV)","signal","hMass", 1);
  plot_TF1d("rec_default","P(M|Y)", "M (GeV)","signal","hMassProb", 1, "_2peak_BAD");
  plot_TF1d("rec_default","P(M|Y)", "M (GeV)","signal","hMassProb", 3, "_2peak_GOOD");
  plot_TF1d("rec_default","P(M|Y)", "M (GeV)","signal","hMassProb", 4, "_1peak_GOOD");
  plot_TF1d("rec_default","P(M|Y)", "M (GeV)","signal","hMassProb", 34,"_3peak_BAD");

  return;

  plot_BestMass("gen_default", "Signal (parton)", "hBestMass");
  plot_BestMass("gen_default", "Signal (parton)", "hBestMassProb");

  plot_BestMass("rec_default", "Signal (nominal res.)", "hBestMass");
  plot_BestMass("rec_default", "Signal (nominal res.)", "hBestMassProb");
  plot_BestMass("rec_METup", "Signal (+50% MET)",      "hBestMassProb");
  plot_BestMass("rec_ScaleLup", "Signal (+50% udcsg)", "hBestMassProb");
  plot_BestMass("rec_ScaleBup", "Signal (+50% b)",     "hBestMassProb");

  plot_BestMass("gen_noJac", "Signal (parton) no Jac.", "hBestMassProb");
  plot_BestMass("gen_noPDF", "Signal (parton) no PDF", "hBestMassProb");
  plot_BestMass("gen_noMET", "Signal (parton) no E_{T}^{miss} TF", "hBestMassProb");
  plot_BestMass("gen_noME",  "Signal (parton) no PS", "hBestMassProb");
  plot_BestMass("gen_noTF",  "Signal (parton) no jet TF", "hBestMassProb");

  plot_BestMass("gen_default_largerLR_range", "Signal (parton) [60,210]", "hBestMassProb");
  plot_BestMass("rec_default_largerLR_range", "Signal (res.) [60,210]", "hBestMassProb");
  */
}
