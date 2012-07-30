#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMatrixT.h"
#include "TMatrixTBase.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TCut.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TBenchmark.h"
#include "TGraph.h"
#include "TVectorT.h"
#include "TMultiGraph.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TStyle.h"


#include <vector>
#include <utility>
#include <map>
#include <algorithm>

float pileupWeight2( int intimepileup_ ){

  float weights[53];
  weights[0]=1.0;
  weights[1]=0.241725;
  weights[2]=0.497953;
  weights[3]=0.750906;
  weights[4]=0.958778;
  weights[5]=1.15292;
  weights[6]=1.26729;
  weights[7]=1.33763;
  weights[8]=1.39117;
  weights[9]=1.38692;
  weights[10]=1.41177;
  weights[11]=1.37077;
  weights[12]=1.34191;
  weights[13]=1.27619;
  weights[14]=1.20034;
  weights[15]=1.1264;
  weights[16]=1.00513;
  weights[17]=0.898563;
  weights[18]=0.783283;
  weights[19]=0.660026;
  weights[20]=0.545681;
  weights[21]=0.444979;
  weights[22]=0.355539;
  weights[23]=0.278989;
  weights[24]=0.214793;
  weights[25]=0.161305;
  weights[26]=0.12141;
  weights[27]=0.089384;
  weights[28]=0.0655027;
  weights[29]=0.0470954;
  weights[30]=0.033824;
  weights[31]=0.0241277;
  weights[32]=0.0168523;
  weights[33]=0.0118342;
  weights[34]=0.00831188;
  weights[35]=0.00574736;
  weights[36]=0.00395389;
  weights[37]=0.00270099;
  weights[38]=0.00184071;
  weights[39]=0.00126892;
  weights[40]=0.000799038;
  weights[41]=0.000568358;
  weights[42]=0.000366065;
  weights[43]=0.000241041;
  weights[44]=0.000152796;
  weights[45]=5.53181e-05;
  weights[46]=5.53181e-05;
  weights[47]=5.53181e-05;
  weights[48]=5.53181e-05;
  weights[49]=5.53181e-05;
  weights[50]=5.53181e-05;
  weights[51]=5.53181e-05;
  weights[52]=5.53181e-05;

  if(intimepileup_<52)
    return weights[intimepileup_+1];
  else
    return 5.53181e-05;
}

double mc( int intimepileup_  ){

  Double_t Fall2011[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

  if(intimepileup_<50)
    return Fall2011[intimepileup_];
  else
    return Fall2011[49];


}

void plot(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.68,0.55,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  TH1F* h = new TH1F("h","; number of pile-up interactions; units",35,0,35);
  //h->SetMarkerStyle(kFullSquare);
  //h->SetMarkerColor(kRed);
  //h->SetMarkerSize(1.5);
  h->SetLineWidth(3);
  h->SetLineStyle(kSolid);
  h->SetLineColor(kRed);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  //for(int i = 1; i < 54; i++) h->SetBinContent(i, pileupWeight2(i) );
  for(int i = 0; i < 50; i++) h->SetBinContent(i+1, pileupWeight2(i) );

  h->Draw("HIST");
  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  //leg->AddEntry(h,"MC pu-weight","P");
  leg->AddEntry(h,"MC pu distribution","L");
  leg->Draw();


}


void plot2(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.68,0.55,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);
  
  TH1F* hLeg = new TH1F("hLeg","",1,0,1);
  hLeg->SetLineWidth(3);
  hLeg->SetLineStyle(kSolid);
  hLeg->SetLineColor(kBlue);
  
  //TFile f("../../../TauTauStudies/test/Macro/pileUp/Run2011PileUp.root","READ");
  TFile f("PileUp2011.root","READ");
  
  TH1F* h =  (TH1F*)(((TH1F*)f.Get("pileup"))->Clone("h"));
  h->SetLineWidth(3);
  h->SetLineStyle(kSolid);
  h->SetLineColor(kBlue);
  h->SetTitle("");
  h->SetXTitle("number of pile-up interactions");
  h->SetYTitle("units");
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  
  h->DrawNormalized("HIST");
  
  
  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(hLeg,"#splitline{Data-driven pu}{distribution}","L");
  leg->Draw();


}

void plotMuTrig(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  TH1F* h = new TH1F("h","; muon p_{T} (GeV); efficiency",70,10,80);
  //h->SetLineWidth(3);
  //h->SetLineColor(kRed);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->SetAxisRange(0.30,1.1,"Y");

  TFile corrections("../../../Utilities/data/corrections/llrCorrections.root");

  TF1 *turnOnMuAllBL         = (TF1*)corrections.Get("turnOnMuAllBL");
  TF1 *turnOnMuAllEC         = (TF1*)corrections.Get("turnOnMuAllEC");
  TF1 *turnOnTauMuTauAllBL   = (TF1*)corrections.Get("turnOnTauMuTauAllBL");  
  TF1 *turnOnTauMuTauAllEC   = (TF1*)corrections.Get("turnOnTauMuTauAllEC");  
  TF1 *turnOnTauElecTauAllBL = (TF1*)corrections.Get("turnOnTauElecTauAllBL");
  TF1 *turnOnTauElecTauAllEC = (TF1*)corrections.Get("turnOnTauElecTauAllEC");
  TF1 *turnOnElecAllBL       = (TF1*)corrections.Get("turnOnElecAllBL");
  TF1 *turnOnElecAllEC       = (TF1*)corrections.Get("turnOnElecAllEC");

  h->Draw("");

  turnOnMuAllBL->SetLineWidth(4);
  turnOnMuAllBL->SetLineStyle(kSolid);
  turnOnMuAllBL->SetLineColor(kRed);
  turnOnMuAllBL->Draw("SAME");

  turnOnMuAllEC->SetLineWidth(4);
  turnOnMuAllEC->SetLineStyle(kDashed);
  turnOnMuAllEC->SetLineColor(kBlue);
  turnOnMuAllEC->Draw("SAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(turnOnMuAllBL,"|#eta(#mu)|<1.5","L");
  leg->AddEntry(turnOnMuAllEC,"|#eta(#mu)|>1.5","L");
  leg->Draw();

  gPad->SaveAs("muonHLT.pdf");



}


void plotElecTrig(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  TH1F* h = new TH1F("h","; electron p_{T} (GeV); efficiency",62,18,80);
  //h->SetLineWidth(3);
  //h->SetLineColor(kRed);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->SetAxisRange(0.30,1.1,"Y");

  TFile corrections("../../../Utilities/data/corrections/llrCorrections.root");

  TF1 *turnOnMuAllBL         = (TF1*)corrections.Get("turnOnMuAllBL");
  TF1 *turnOnMuAllEC         = (TF1*)corrections.Get("turnOnMuAllEC");
  TF1 *turnOnTauMuTauAllBL   = (TF1*)corrections.Get("turnOnTauMuTauAllBL");  
  TF1 *turnOnTauMuTauAllEC   = (TF1*)corrections.Get("turnOnTauMuTauAllEC");  
  TF1 *turnOnTauElecTauAllBL = (TF1*)corrections.Get("turnOnTauElecTauAllBL");
  TF1 *turnOnTauElecTauAllEC = (TF1*)corrections.Get("turnOnTauElecTauAllEC");
  TF1 *turnOnElecAllBL       = (TF1*)corrections.Get("turnOnElecAllBL");
  TF1 *turnOnElecAllEC       = (TF1*)corrections.Get("turnOnElecAllEC");

  h->Draw("");

  turnOnElecAllBL->SetLineWidth(4);
  turnOnElecAllBL->SetLineStyle(kSolid);
  turnOnElecAllBL->SetLineColor(kRed);
  turnOnElecAllBL->Draw("SAME");

  turnOnElecAllEC->SetLineWidth(4);
  turnOnElecAllEC->SetLineStyle(kDashed);
  turnOnElecAllEC->SetLineColor(kBlue);
  turnOnElecAllEC->Draw("SAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(turnOnElecAllBL,"|#eta(e)|<1.5","L");
  leg->AddEntry(turnOnElecAllEC,"|#eta(e)|>1.5","L");
  leg->Draw();

  gPad->SaveAs("eleHLT.pdf");



}


void plotTauMuTrig(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  TH1F* h = new TH1F("h","; tau p_{T} (GeV); efficiency",65,15,80);
  //h->SetLineWidth(3);
  //h->SetLineColor(kRed);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->SetAxisRange(0.30,1.1,"Y");

  TFile corrections("../../../Utilities/data/corrections/llrCorrections.root");

  TF1 *turnOnMuAllBL         = (TF1*)corrections.Get("turnOnMuAllBL");
  TF1 *turnOnMuAllEC         = (TF1*)corrections.Get("turnOnMuAllEC");
  TF1 *turnOnTauMuTauAllBL   = (TF1*)corrections.Get("turnOnTauMuTauAllBL");  
  TF1 *turnOnTauMuTauAllEC   = (TF1*)corrections.Get("turnOnTauMuTauAllEC");  
  TF1 *turnOnTauElecTauAllBL = (TF1*)corrections.Get("turnOnTauElecTauAllBL");
  TF1 *turnOnTauElecTauAllEC = (TF1*)corrections.Get("turnOnTauElecTauAllEC");
  TF1 *turnOnElecAllBL       = (TF1*)corrections.Get("turnOnElecAllBL");
  TF1 *turnOnElecAllEC       = (TF1*)corrections.Get("turnOnElecAllEC");

  h->Draw("");

  turnOnTauMuTauAllBL->SetLineWidth(4);
  turnOnTauMuTauAllBL->SetLineStyle(kSolid);
  turnOnTauMuTauAllBL->SetLineColor(kRed);
  turnOnTauMuTauAllBL->Draw("SAME");

  turnOnTauMuTauAllEC->SetLineWidth(4);
  turnOnTauMuTauAllEC->SetLineStyle(kDashed);
  turnOnTauMuTauAllEC->SetLineColor(kBlue);
  turnOnTauMuTauAllEC->Draw("SAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(turnOnTauMuTauAllBL,"|#eta(#tau)|<1.5 (#mu+#tau)","L");
  leg->AddEntry(turnOnTauMuTauAllEC,"|#eta(#tau)|>1.5 (#mu+#tau)","L");
  leg->Draw();

  gPad->SaveAs("tauHLT_muTau.pdf");



}


void plotTauElecTrig(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  TH1F* h = new TH1F("h","; tau p_{T} (GeV); efficiency",65,15,80);
  //h->SetLineWidth(3);
  //h->SetLineColor(kRed);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->SetAxisRange(0.30,1.1,"Y");

  TFile corrections("../../../Utilities/data/corrections/llrCorrections.root");

  TF1 *turnOnMuAllBL         = (TF1*)corrections.Get("turnOnMuAllBL");
  TF1 *turnOnMuAllEC         = (TF1*)corrections.Get("turnOnMuAllEC");
  TF1 *turnOnTauMuTauAllBL   = (TF1*)corrections.Get("turnOnTauMuTauAllBL");  
  TF1 *turnOnTauMuTauAllEC   = (TF1*)corrections.Get("turnOnTauMuTauAllEC");  
  TF1 *turnOnTauElecTauAllBL = (TF1*)corrections.Get("turnOnTauElecTauAllBL");
  TF1 *turnOnTauElecTauAllEC = (TF1*)corrections.Get("turnOnTauElecTauAllEC");
  TF1 *turnOnElecAllBL       = (TF1*)corrections.Get("turnOnElecAllBL");
  TF1 *turnOnElecAllEC       = (TF1*)corrections.Get("turnOnElecAllEC");

  h->Draw("");

  turnOnTauElecTauAllBL->SetLineWidth(4);
  turnOnTauElecTauAllBL->SetLineStyle(kSolid);
  turnOnTauElecTauAllBL->SetLineColor(kRed);
  turnOnTauElecTauAllBL->Draw("SAME");

  turnOnTauElecTauAllEC->SetLineWidth(4);
  turnOnTauElecTauAllEC->SetLineStyle(kDashed);
  turnOnTauElecTauAllEC->SetLineColor(kBlue);
  turnOnTauElecTauAllEC->Draw("SAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(turnOnTauElecTauAllBL,"|#eta(#tau)|<1.5 (e+#tau)","L");
  leg->AddEntry(turnOnTauElecTauAllEC,"|#eta(#tau)|>1.5 (e+#tau)","L");
  leg->Draw();

  gPad->SaveAs("tauHLT_eleTau.pdf");



}


void plotMuIso(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.055);

  //TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleGGFH120-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
  TTree* tree = (TTree*)f->Get("outTreePtOrdRaw");

  TH1F* hDen = new TH1F("hDen","; number of vertices; isolation efficiency",30, 0, 30 );
  TH1F* hNum = new TH1F("hNum","; number of vertices; isolation efficiency",30, 0, 30 );
 
  hNum->SetMarkerSize(1.4);
  hNum->SetMarkerStyle(kFullCircle);
  hNum->SetMarkerColor(kRed);
  hNum->GetXaxis()->SetTitleSize(0.055);
  hNum->GetXaxis()->SetTitleOffset(0.80);
  hNum->GetYaxis()->SetTitleSize(0.055);
  hNum->GetYaxis()->SetTitleOffset(0.80);
  hNum->SetAxisRange(0.0,1.2,"Y");

  TH1F* hNum2 = new TH1F("hNum2","; number of vertices; isolation efficiency",30, 0, 30 );
 
  hNum2->SetMarkerSize(1.4);
  hNum2->SetMarkerStyle(kOpenCircle);
  hNum2->SetMarkerColor(kBlue);
  hNum2->GetXaxis()->SetTitleSize(0.055);
  hNum2->GetXaxis()->SetTitleOffset(0.80);
  hNum2->GetYaxis()->SetTitleSize(0.055);
  hNum2->GetYaxis()->SetTitleOffset(0.80);
  hNum2->SetAxisRange(0.0,1.2,"Y");

  TCut cutDen("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && MtLeg1<40 && tightestHPSMVAWP>=0");
  TCut cutNum("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && MtLeg1<40 && tightestHPSMVAWP>=0 && combRelIsoLeg1DBetav2<0.10");
  TCut cutNum2("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && MtLeg1<40 && tightestHPSMVAWP>=0 && combRelIsoLeg1<0.10");

  hDen->Sumw2();
  hNum->Sumw2();
  hNum2->Sumw2();

  tree->Draw("numPV>>hDen","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen);
  tree->Draw("numPV>>hNum","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutNum);
  tree->Draw("numPV>>hNum2","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutNum2);


  hNum->Divide(hDen);
  hNum2->Divide(hDen);

  hNum->Draw("PE");
  hNum2->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Simulation #sqrt{s}=7 TeV}");
  leg->AddEntry(hNum2,"w/o #Delta#beta","P");
  leg->AddEntry(hNum,"w/ #Delta#beta","P");
  leg->Draw();

}


void plotTauIso(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.050);

  TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
  TTree* tree = (TTree*)f->Get("outTreePtOrdRaw");

  TH1F* hDen = new TH1F("hDen","; gen #tau p_{T} (GeV); isolation efficiency", 15, 15, 90 );
  TH1F* hNum = new TH1F("hNum","; gen #tau p_{T} (GeV) ;isolation efficiency",15, 15, 90 );
 
  hNum->SetMarkerSize(1.4);
  hNum->SetMarkerStyle(kFullCircle);
  hNum->SetMarkerColor(kRed);
  hNum->GetXaxis()->SetTitleSize(0.055);
  hNum->GetXaxis()->SetTitleOffset(0.80);
  hNum->GetYaxis()->SetTitleSize(0.055);
  hNum->GetYaxis()->SetTitleOffset(0.80);
  hNum->SetAxisRange(0.0,1.1,"Y");

  TH1F* hNum2 = new TH1F("hNum2","; gen #tau p_{T} (GeV); isolation efficiency",15, 15, 90 );
 
  hNum2->SetMarkerSize(1.4);
  hNum2->SetMarkerStyle(kOpenCircle);
  hNum2->SetMarkerColor(kBlue);
  hNum2->GetXaxis()->SetTitleSize(0.055);
  hNum2->GetXaxis()->SetTitleOffset(0.80);
  hNum2->GetYaxis()->SetTitleSize(0.055);
  hNum2->GetYaxis()->SetTitleOffset(0.80);
  hNum2->SetAxisRange(0.0,1.1,"Y");

  TH1F* hNum3 = new TH1F("hNum3","; number of vertices; isolation efficiency",15, 15, 90 );
 
  hNum3->SetMarkerSize(1.4);
  hNum3->SetMarkerStyle(kOpenSquare);
  hNum3->SetMarkerColor(kMagenta);
  hNum3->GetXaxis()->SetTitleSize(0.055);
  hNum3->GetXaxis()->SetTitleOffset(0.80);
  hNum3->GetYaxis()->SetTitleSize(0.055);
  hNum3->GetYaxis()->SetTitleOffset(0.80);
  hNum3->SetAxisRange(0.0,1.1,"Y");

  TCut cutDen( "isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && isTauLegMatched && TMath::Abs(genTauEta)<2.3");
  TCut cutNumL("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && isTauLegMatched && TMath::Abs(genTauEta)<2.3 && tightestHPSMVAWP>=0");
  TCut cutNumM("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && isTauLegMatched && TMath::Abs(genTauEta)<2.3 && tightestHPSDBWP>0");//ightestHPSMVAWP>=1
  TCut cutNumT("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && isTauLegMatched && TMath::Abs(genTauEta)<2.3 && tightestHPSMVAWP>=2");

  hDen->Sumw2();
  hNum->Sumw2();
  hNum2->Sumw2();
  hNum3->Sumw2();

  tree->Draw("genTauPt>>hDen", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen);
  tree->Draw("genTauPt>>hNum", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutNumL);
  tree->Draw("genTauPt>>hNum2","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutNumM);
  tree->Draw("genTauPt>>hNum3","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutNumT);


  hNum->Divide(hDen);
  hNum2->Divide(hDen);
  hNum3->Divide(hDen);

  hNum->Draw("PE");
  hNum2->Draw("PESAME");
  //hNum3->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Z#rightarrow#tau#tau MC #sqrt{s}=7 TeV}");
  leg->AddEntry(hNum, "MVA Loose","P");
  leg->AddEntry(hNum2,"#Delta#beta Loose","P");
  //leg->AddEntry(hNum3,"Tight","P");
  leg->Draw();

}


void plotTauFakeIso(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.050);

  TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleRun2011-MuTau-All_run_Open_MuTauStream.root");
  TTree* tree = (TTree*)f->Get("outTreePtOrdRaw");

  TH1F* hDen = new TH1F("hDen","; #tau-jet p_{T} (GeV); fake rate", 15, 15, 90 );
  TH1F* hNum = new TH1F("hNum","; #tau-jet p_{T} (GeV); fake rate", 15, 15, 90 );
 
  hNum->SetMarkerSize(1.4);
  hNum->SetMarkerStyle(kFullCircle);
  hNum->SetMarkerColor(kRed);
  hNum->GetXaxis()->SetTitleSize(0.055);
  hNum->GetXaxis()->SetTitleOffset(0.80);
  hNum->GetYaxis()->SetTitleSize(0.055);
  hNum->GetYaxis()->SetTitleOffset(0.80);
  hNum->SetAxisRange(0.001,1.1,"Y");

  TH1F* hNum2 = new TH1F("hNum2","; #tau-jet p_{T} (GeV); fake rate",15, 15, 90 );
 
  hNum2->SetMarkerSize(1.4);
  hNum2->SetMarkerStyle(kOpenCircle);
  hNum2->SetMarkerColor(kBlue);
  hNum2->GetXaxis()->SetTitleSize(0.055);
  hNum2->GetXaxis()->SetTitleOffset(0.80);
  hNum2->GetYaxis()->SetTitleSize(0.055);
  hNum2->GetYaxis()->SetTitleOffset(0.80);
  hNum2->SetAxisRange(0.001,1.1,"Y");

  TH1F* hNum3 = new TH1F("hNum3","; #tau-jet p_{T} (GeV); fake rate",15, 15, 90 );
 
  hNum3->SetMarkerSize(1.4);
  hNum3->SetMarkerStyle(kOpenSquare);
  hNum3->SetMarkerColor(kMagenta);
  hNum3->GetXaxis()->SetTitleSize(0.055);
  hNum3->GetXaxis()->SetTitleOffset(0.80);
  hNum3->GetYaxis()->SetTitleSize(0.055);
  hNum3->GetYaxis()->SetTitleOffset(0.80);
  hNum3->SetAxisRange(0.001,1.1,"Y");

  TCut cutDen( "pairIndex<1 && isTightMuon && isPFMuon && muFlag==0 && diTauCharge!=0 && HLTmatch && HLTx && MtLeg1Corr<40 && nJets20BTagged<1 && combRelIsoLeg1DBetav2>0.30");
  TCut cutNumL("pairIndex<1 && isTightMuon && isPFMuon && muFlag==0 && diTauCharge!=0 && HLTmatch && HLTx && MtLeg1Corr<40 && nJets20BTagged<1 && combRelIsoLeg1DBetav2>0.30 && tightestHPSMVAWP>=0");
  TCut cutNumM("pairIndex<1 && isTightMuon && isPFMuon && muFlag==0 && diTauCharge!=0 && HLTmatch && HLTx && MtLeg1Corr<40 && nJets20BTagged<1 && combRelIsoLeg1DBetav2>0.30 && tightestHPSDBWP>0");//tightestHPSMVAWP>=1
  TCut cutNumT("pairIndex<1 && isTightMuon && isPFMuon && muFlag==0 && diTauCharge!=0 && HLTmatch && HLTx && MtLeg1Corr<40 && nJets20BTagged<1 && combRelIsoLeg1DBetav2>0.30 && tightestHPSMVAWP>=2");

  hDen->Sumw2();
  hNum->Sumw2();
  hNum2->Sumw2();
  hNum3->Sumw2();

  tree->Draw("pfJetPt>>hDen", cutDen);
  tree->Draw("pfJetPt>>hNum", cutNumL);
  tree->Draw("pfJetPt>>hNum2",cutNumM);
  tree->Draw("pfJetPt>>hNum3",cutNumT);


  hNum->Divide(hDen);
  hNum2->Divide(hDen);
  hNum3->Divide(hDen);

  hNum->Draw("PE");
  hNum2->Draw("PESAME");
  //hNum3->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Data #sqrt{s}=7 TeV}");
  leg->AddEntry(hNum, "MVA Loose","P");
  leg->AddEntry(hNum2,"#Delta#beta Loose","P");
  //leg->AddEntry(hNum3,"Tight","P");
  leg->Draw();

}

void plotJetID(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.050);

  TFile* f = new TFile("/data_CMS/cms/lbianchini/MuTauStreamFall11_04May2012_Approval/treeMuTauStream_VBFH120-MuTau-powheg-PUS6_run.root");
  TTree* tree = (TTree*)f->Get("muTauStreamAnalyzerRaw/tree");
  TH1F* hAll = new TH1F("hAll","",15,0,5.0);
  TH1F* hAll2 = new TH1F("hAll2","",15,0,5.0);
  

  TH1F* hPassL = new TH1F("hPassL","; jet |#eta|; pile-up jet-ID efficiency",15,0.0,5.0);
  hPassL->SetMarkerSize(1.4);
  hPassL->SetMarkerStyle(kFullCircle);
  hPassL->SetMarkerColor(kRed);
  hPassL->GetXaxis()->SetTitleSize(0.055);
  hPassL->GetXaxis()->SetTitleOffset(0.80);
  hPassL->GetYaxis()->SetTitleSize(0.055);
  hPassL->GetYaxis()->SetTitleOffset(0.80);
  hPassL->SetAxisRange(0.0,1.2,"Y");

  TH1F* hPassL2 = new TH1F("hPassL2","; jet #eta; pile-up jet-ID efficiency",15,0,5.0);
  hPassL2->SetMarkerSize(1.4);
  hPassL2->SetMarkerStyle(kOpenCircle);
  hPassL2->SetMarkerColor(kBlue);
  hPassL2->GetXaxis()->SetTitleSize(0.055);
  hPassL2->GetXaxis()->SetTitleOffset(0.80);
  hPassL2->GetYaxis()->SetTitleSize(0.055);
  hPassL2->GetYaxis()->SetTitleOffset(0.80);
  hPassL2->SetAxisRange(0.0,1.2,"Y");


  hAll->Sumw2(); 
  hPassL->Sumw2();
  hPassL2->Sumw2();

  tree->Draw("TMath::Abs(jetsIDP4[0].Eta())>>hAll",   "isTauLegMatched && jetsIDP4@.size()>0 && jetsIDP4[0].Pt()>30 && genJetsIDP4[0].Pt()>30 && TMath::Abs(genJetsIDP4[0].Eta()-jetsIDP4[0].Eta())<0.2 && TMath::Abs(genJetsIDP4[0].Phi()-jetsIDP4[0].Phi())<0.2 && jetPUWP[0]>-99 && numPV<=10");
  tree->Draw("TMath::Abs(jetsIDP4[0].Eta())>>hAll2",  "isTauLegMatched && jetsIDP4@.size()>0 && jetsIDP4[0].Pt()>30 && genJetsIDP4[0].Pt()>30 && TMath::Abs(genJetsIDP4[0].Eta()-jetsIDP4[0].Eta())<0.2 && TMath::Abs(genJetsIDP4[0].Phi()-jetsIDP4[0].Phi())<0.2 && jetPUWP[0]>-99 && numPV>10");
  tree->Draw("TMath::Abs(jetsIDP4[0].Eta())>>hPassL", "isTauLegMatched && jetsIDP4@.size()>0 && jetsIDP4[0].Pt()>30 && genJetsIDP4[0].Pt()>30 && TMath::Abs(genJetsIDP4[0].Eta()-jetsIDP4[0].Eta())<0.2 && TMath::Abs(genJetsIDP4[0].Phi()-jetsIDP4[0].Phi())<0.2 && jetPUWP[0]>0.5 && numPV<=10");
  tree->Draw("TMath::Abs(jetsIDP4[0].Eta())>>hPassL2","isTauLegMatched && jetsIDP4@.size()>0 && jetsIDP4[0].Pt()>30 && genJetsIDP4[0].Pt()>30 && TMath::Abs(genJetsIDP4[0].Eta()-jetsIDP4[0].Eta())<0.2 && TMath::Abs(genJetsIDP4[0].Phi()-jetsIDP4[0].Phi())<0.2 && jetPUWP[0]>0.5 && numPV>10");


  hPassL->Divide(hAll);
  hPassL2->Divide(hAll2);
  hPassL->Draw("PE");
  hPassL2->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Simulation #sqrt{s}=7 TeV}");
  leg->AddEntry(hPassL,  "Loose #PV#leq10 ", "P");
  leg->AddEntry(hPassL2, "Loose #PV>10", "P");
  leg->Draw();


}


void plotRecoil(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.38,0.25,0.68,0.55,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.050);

  //TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
  //TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleVBFH120-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleWJets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root");

  TTree* tree = (TTree*)f->Get("outTreePtOrdRaw");

  TH1F* hRaw = new TH1F("hRaw","; E_{T}^{miss} (GeV); units",18, 0, 120 );
  TH1F* hCor = new TH1F("hCor","; number of vertices; units",18, 0, 120 );

  TH1F* hRaw2 = new TH1F("hRaw2","; number of vertices; units",18, 0, 120 );
  TH1F* hCor2 = new TH1F("hCor2","; number of vertices; units",18, 0, 120 );

  TH1F* hRaw3 = new TH1F("hRaw3","; number of vertices; units",18, 0, 120 );
  TH1F* hCor3 = new TH1F("hCor3","; number of vertices; units",18, 0, 120 );

  hRaw->SetLineWidth(3);
  hRaw->SetLineColor(kRed);
  hRaw->GetXaxis()->SetTitleSize(0.055);
  hRaw->GetXaxis()->SetTitleOffset(0.80);
  hRaw->GetYaxis()->SetTitleSize(0.055);
  hRaw->GetYaxis()->SetTitleOffset(0.85);
  //hRaw->SetAxisRange(0.0,0.16,"Y");

  hCor->SetLineWidth(3);
  hCor->SetLineStyle(kDashed);
  hCor->SetLineColor(kRed);
  hCor->GetXaxis()->SetTitleSize(0.055);
  hCor->GetXaxis()->SetTitleOffset(0.80);
  hCor->GetYaxis()->SetTitleSize(0.055);
  hCor->GetYaxis()->SetTitleOffset(0.80);
  //hCor->SetAxisRange(0.0,1.2,"Y");

  hRaw2->SetLineWidth(3);
  hRaw2->SetLineColor(kBlue);
  hRaw2->GetXaxis()->SetTitleSize(0.055);
  hRaw2->GetXaxis()->SetTitleOffset(0.80);
  hRaw2->GetYaxis()->SetTitleSize(0.055);
  hRaw2->GetYaxis()->SetTitleOffset(0.80);
  //hRaw->SetAxisRange(0.0,1.2,"Y");

  hCor2->SetLineWidth(3);
  hCor2->SetLineStyle(kDashed);
  hCor2->SetLineColor(kBlue);
  hCor2->GetXaxis()->SetTitleSize(0.055);
  hCor2->GetXaxis()->SetTitleOffset(0.80);
  hCor2->GetYaxis()->SetTitleSize(0.055);
  hCor2->GetYaxis()->SetTitleOffset(0.80);
  //hCor->SetAxisRange(0.0,1.2,"Y");

  hRaw3->SetLineWidth(3);
  hRaw3->SetLineColor(kMagenta);
  hRaw3->GetXaxis()->SetTitleSize(0.055);
  hRaw3->GetXaxis()->SetTitleOffset(0.80);
  hRaw3->GetYaxis()->SetTitleSize(0.055);
  hRaw3->GetYaxis()->SetTitleOffset(0.80);
  //hRaw->SetAxisRange(0.0,1.2,"Y");

  hCor3->SetLineWidth(3);
  hCor3->SetLineStyle(kDashed);
  hCor3->SetLineColor(kMagenta);
  hCor3->GetXaxis()->SetTitleSize(0.055);
  hCor3->GetXaxis()->SetTitleOffset(0.80);
  hCor3->GetYaxis()->SetTitleSize(0.055);
  hCor3->GetYaxis()->SetTitleOffset(0.80);
  //hCor->SetAxisRange(0.0,1.2,"Y");


  TCut cutDen("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && tightestHPSMVAWP>=0 && combRelIsoLeg1DBetav2<0.10 && !isTauLegMatched && nJets30==0");
  TCut cutDen2("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && tightestHPSMVAWP>=0 && combRelIsoLeg1DBetav2<0.10 && nJets30==1 && !isTauLegMatched");
  TCut cutDen3("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && tightestHPSMVAWP>=0 && combRelIsoLeg1DBetav2<0.10 && nJets30>=2 && !isTauLegMatched");

  tree->Draw("MEt>>hRaw","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen);
  tree->Draw("MEtCorr>>hCor","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen);

  tree->Draw("MEt>>hRaw2","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen2);
  tree->Draw("MEtCorr>>hCor2","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen2);

  tree->Draw("MEt>>hRaw3","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen3);
  tree->Draw("MEtCorr>>hCor3","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cutDen3);

  hRaw->DrawNormalized();
  hCor->DrawNormalized("SAME");

  hRaw2->DrawNormalized("SAME");
  hCor2->DrawNormalized("SAME");

  hRaw3->DrawNormalized("SAME");
  hCor3->DrawNormalized("SAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{W#rightarrow l#nu  #sqrt{s}=7 TeV}");
  leg->AddEntry(hRaw,"Raw N_{jets}=0","L");
  leg->AddEntry(hCor,"Corr N_{jets}=0","L");
  leg->AddEntry(hRaw2,"Raw N_{jets}=1","L");
  leg->AddEntry(hCor2,"Corr N_{jets}=1","L");
  leg->AddEntry(hRaw3,"Raw N_{jets}>1","L");
  leg->AddEntry(hCor3,"Corr N_{jets}>1","L");
  leg->Draw();
}


void plotVBF(TString var = "Deta", int nBins = 10, float xLow = 0, float xHigh = 10, TString title = "#Delta#eta", int sgn=0){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.45,0.61,0.65,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.050);


  TFile* fSgn = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleVBFH120-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fBkg = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");

  TTree* treeSgn = (TTree*)fSgn->Get("outTreePtOrdRaw");
  TTree* treeBkg = (TTree*)fBkg->Get("outTreePtOrdRaw");

  TH1F* hSgn = new TH1F("hSgn","; "+title+"; units",nBins, xLow, xHigh );
  TH1F* hBkg = new TH1F("hBkg","; "+title+"; units",nBins, xLow, xHigh );

  hSgn->SetLineWidth(1);
  hSgn->SetLineColor(kBlue);
  hSgn->SetFillStyle(3004);
  hSgn->SetFillColor(kBlue);
  hSgn->GetXaxis()->SetTitleSize(0.055);
  hSgn->GetXaxis()->SetTitleOffset(0.80);
  hSgn->GetYaxis()->SetTitleSize(0.055);
  hSgn->GetYaxis()->SetTitleOffset(0.85);
  //hRaw->SetAxisRange(0.0,0.16,"Y");

  hBkg->SetLineWidth(3);
  hBkg->SetLineStyle(kDashed);
  hBkg->SetLineColor(kRed);
  hBkg->GetXaxis()->SetTitleSize(0.055);
  hBkg->GetXaxis()->SetTitleOffset(0.80);
  hBkg->GetYaxis()->SetTitleSize(0.055);
  hBkg->GetYaxis()->SetTitleOffset(0.80);



  TCut cut("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && tightestHPSMVAWP>=0 && combRelIsoLeg1DBetav2<0.10 && nJets30>1 && (ptVeto<30 || isVetoInJets!=1) && isTauLegMatched");

  treeSgn->Draw(var+">>hSgn","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cut);
  treeBkg->Draw(var+">>hBkg","puWeight2*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*cut);

  if(sgn){
    hSgn->DrawNormalized();
    hBkg->DrawNormalized("SAME");
  }
  else{
    hBkg->DrawNormalized();
    hSgn->DrawNormalized("SAME");
  }

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Simulation #sqrt{s}=7 TeV}");
  leg->AddEntry(hSgn,"qq #rightarrow H #rightarrow #tau#tau","F");
  leg->AddEntry(hBkg,"Z#rightarrow #tau#tau","L");
  leg->Draw();

}

void plotVBFAll(){


  plotVBF("diJetPt",20,0,500,"p_{T}(jj) (GeV)",0);
  gPad->SaveAs("diJetPt.pdf");
  plotVBF("Dphi",10,0,3.2,"#Delta#phi(jj)",0);
  gPad->SaveAs("Dphi.pdf");
  plotVBF("Mjj",30,0,1500,"M(jj) (GeV)",0);
  gPad->SaveAs("Mjj.pdf");
  plotVBF("Deta",15,0,8,"#Delta#eta(jj)",0);
  gPad->SaveAs("Deta.pdf");
  plotVBF("diTauRecoPt",20,0,500,"p_{T}(H) (GeV)",0);
  gPad->SaveAs("diTauRecoPt.pdf");
  plotVBF("dPhiHjet",20,0,3.2,"#Delta#phi(Hjj)",1);
  gPad->SaveAs("dPhiHjet.pdf");
  plotVBF("c2",20,0,300,"p_{T}(#tau#tau) (GeV)",0);
  gPad->SaveAs("c2.pdf");
  plotVBF("c1",20,0,5,"#Delta#eta(min)",0);
  gPad->SaveAs("c1.pdf");
  plotVBF("Dphi",8,0,3.2,"#Delta#phi(jj)",0);
  gPad->SaveAs("Dphi.pdf");




}



void plotWshape(string ptL2 = ">20", string nJets20BTagged = "<1", float nJets30 = 0, string label = " W-Sdb",
		TString title = "Wshape0jetPt2040Wsdb.pdf"){


 TFile *fData = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleRun2011-MuTau-All_run_Open_MuTauStream.root", "READ");
 TTree *data = (TTree*)fData->Get("outTreePtOrdRaw");

 TFile *fMadGraph = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleW3Jets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root", "READ");
 TTree *madGraph = (TTree*)fMadGraph->Get("outTreePtOrdRaw");

 TFile *fMadGraphDY = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root", "READ");
 TTree *madGraphDYIncl = (TTree*)fMadGraphDY->Get("outTreePtOrdRaw");

 TFile *fMadGraphIncl = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleWJets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root", "READ");
 TTree *madGraphIncl = (TTree*)fMadGraphIncl->Get("outTreePtOrdRaw");

 TFile *fTTIncl = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleTTJets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root", "READ");
 TTree *ttIncl = (TTree*)fTTIncl->Get("outTreePtOrdRaw");




 TCut lpt("ptL1>17 && isPFMuon && isTightMuon");
 TCut tpt(Form("ptL2%s",ptL2.c_str()));

 ////// TAU ISO //////
 TCut tiso("tightestHPSMVAWP>=0");
 //TCut tiso("hpsMVA<0.80");     // <-----------------------
 TCut ltiso("tightestHPSMVAWP>-99");
 TCut mtiso("hpsMVA>0.50");

 ////// MU ISO ///////
 TCut liso("combRelIsoLeg1DBetav2<0.10");
 TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
 TCut lliso("combRelIsoLeg1DBetav2<0.30");


 ////// EVENT WISE //////
 TCut lveto("muFlag==0"); //<----------------------------------
 TCut SS("diTauCharge!=0");
 TCut OS("diTauCharge==0");
 TCut apZ(Form("((MtLeg1Corr)>%f)",70.));
 TCut hltevent("pairIndex<1 && HLTx==1 && ( run>=163269 || run==1)");
 TCut hltmatch("HLTmatch==1");

 TCut vbf(Form("nJets30>=%f && nJets20BTagged%s",nJets30,nJets20BTagged.c_str()));
 //TCut vbf("nJets30>1  && MVAvbf>0.5 &&  (ptVeto<30 || isVetoInJets!=1)");
 //TCut vbf("nJets30>1  && Mjj>400 && Deta>3 &&  (ptVeto<30 || isVetoInJets!=1) && nJets30<999");


 TCut sbinPZetaRelInclusive;
 sbinPZetaRelInclusive     = vbf  && lpt && tpt && tiso && liso  && lveto && OS && hltevent && hltmatch;
 TCut sbinPZetaRelAIsoInclusive;
 sbinPZetaRelAIsoInclusive = vbf  && lpt && tpt && tiso && laiso && lveto && OS && hltevent && hltmatch; 

 TCanvas *c11 = new TCanvas("c11","",5,30,650,600);
 c11->SetGrid(0,0);
 c11->SetFillStyle(4000);
 c11->SetFillColor(10);
 c11->SetTicky();
 c11->SetObjectStat(0);

 TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
 TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);
 
 pad1->SetFillColor(0);
 pad2->SetFillColor(0);
 pad1->Draw();
 pad2->Draw();
 
 pad1->cd();
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


 TH1F* hMassData = new TH1F("hMassData","",20,0,400);
 hMassData->Sumw2();
 data->Draw("diTauNSVfitMass>>hMassData", sbinPZetaRelInclusive&&apZ);
 hMassData->SetMarkerStyle(kFullCircle);
 hMassData->SetMarkerSize(1.3);

 TH1F* hMassMadGraph = new TH1F("hMassMadGraph"," ; SVfit mass (GeV) ; units",20,0,400);
 hMassMadGraph->Sumw2();
 madGraphIncl->Draw("diTauNSVfitMass>>hMassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 hMassMadGraph->SetLineColor(kRed);




 TH1F* h2MassMadGraph = new TH1F("h2MassMadGraph"," ; SVfit mass (GeV) ; units",20,0,400);
 h2MassMadGraph->SetLineWidth(1);
 h2MassMadGraph->SetFillStyle(3003);
 h2MassMadGraph->SetFillColor(kBlue);
 h2MassMadGraph->Sumw2();
 h2MassMadGraph->GetXaxis()->SetTitleSize(0.055);
 h2MassMadGraph->GetXaxis()->SetTitleOffset(0.80);
 h2MassMadGraph->GetYaxis()->SetTitleSize(0.055);
 h2MassMadGraph->GetYaxis()->SetTitleOffset(0.85);

 TH1F* h2MassMadGraphUp = (TH1F*)h2MassMadGraph->Clone("h2MassMadGraphUp");
 h2MassMadGraphUp->SetFillStyle(0);
 h2MassMadGraphUp->SetLineColor(kRed);
 h2MassMadGraphUp->SetLineWidth(2);
 h2MassMadGraphUp->SetLineStyle(kDashed);
 TH1F* h2MassMadGraphDown = (TH1F*)h2MassMadGraph->Clone("h2MassMadGraphDown");
 h2MassMadGraphDown->SetFillStyle(0);
 h2MassMadGraphDown->SetLineColor(kRed);
 h2MassMadGraphDown->SetLineWidth(2);
 h2MassMadGraphDown->SetLineStyle(kDotted);

 madGraphIncl->Draw("diTauNSVfitMass>>h2MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 madGraphIncl->Draw("diTauNSVfitMass*(1+0.03)>>h2MassMadGraphUp", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 madGraphIncl->Draw("diTauNSVfitMass*(1-0.03)>>h2MassMadGraphDown", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));

 float Wincl = h2MassMadGraph->Integral();
 cout << "W = " << Wincl << endl;

 TH1F* h3MassMadGraph = new TH1F("h3MassMadGraph","",20,0,400);
 h3MassMadGraph->Sumw2();

 TH1F* h3MassMadGraphUp = (TH1F*)h3MassMadGraph->Clone("h3MassMadGraphUp");
 h3MassMadGraphUp->SetFillStyle(0);
 h3MassMadGraphUp->SetLineColor(kRed);
 h3MassMadGraphUp->SetLineWidth(2);
 h3MassMadGraphUp->SetLineStyle(kDashed);
 TH1F* h3MassMadGraphDown = (TH1F*)h3MassMadGraph->Clone("h3MassMadGraphDown");
 h3MassMadGraphDown->SetFillStyle(0);
 h3MassMadGraphDown->SetLineColor(kRed);
 h3MassMadGraphDown->SetLineWidth(2);
 h3MassMadGraphDown->SetLineStyle(kDotted);

 ttIncl->Draw("diTauNSVfitMass>>h3MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 ttIncl->Draw("diTauNSVfitMass*(1+0.03)>>h3MassMadGraphUp", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 ttIncl->Draw("diTauNSVfitMass*(1-0.03)>>h3MassMadGraphDown", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));

 float TTincl = h3MassMadGraph->Integral();
 cout << "TT = " << TTincl << endl;

 TH1F* h4MassMadGraph = new TH1F("h4MassMadGraph","",20,0,400);
 h4MassMadGraph->Sumw2();
 madGraphDYIncl->Draw("diTauNSVfitMass>>h4MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightMu*SFTau*SFMu)"*(sbinPZetaRelInclusive&&apZ));
 float DYincl = h4MassMadGraph->Integral();
 cout << "DY = " << DYincl << endl;

 TH1F* h5MassMadGraph = new TH1F("h5MassMadGraph","",20,0,400);
 data->Draw("diTauNSVfitMass>>h5MassMadGraph", (sbinPZetaRelAIsoInclusive&&apZ));
 float AIsoincl = h5MassMadGraph->Integral();
 cout << "AIso = " << AIsoincl << endl;


 h2MassMadGraph->Add(h3MassMadGraph,1.0);
 h2MassMadGraph->Add(h4MassMadGraph,1.0);

 h2MassMadGraphUp->Add(h3MassMadGraphUp,1.0);
 h2MassMadGraphUp->Add(h4MassMadGraph,1.0);

 h2MassMadGraphDown->Add(h3MassMadGraphDown,1.0);
 h2MassMadGraphDown->Add(h4MassMadGraph,1.0);


 h2MassMadGraph->SetLineColor(kBlue);
 h2MassMadGraph->Scale(1./h2MassMadGraph->Integral());
 h2MassMadGraphUp->Scale(1./h2MassMadGraphUp->Integral());
 h2MassMadGraphDown->Scale(1./h2MassMadGraphDown->Integral());

 hMassData->Scale(1./hMassData->Integral());

 h2MassMadGraph->Draw("EHIST");
 h2MassMadGraphUp->Draw("HISTSAME");
 h2MassMadGraphDown->Draw("HISTSAME");

 hMassData->Draw("PESAME");
 //hMassMadGraph->DrawNormalized("ESAME");

 TLegend* leg = new TLegend(0.50,0.51,0.70,0.85,NULL,"brNDC");
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);
 leg->SetFillColor(10);
 leg->SetTextSize(0.050);
 leg->SetHeader(Form("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV %s}",label.c_str()));

 leg->AddEntry(hMassData,"Data","P");
 leg->AddEntry(h2MassMadGraph,"Background","L");
 leg->AddEntry(h2MassMadGraphUp,"Fake-Tau 3% Up","L");
 leg->AddEntry(h2MassMadGraphDown,"Fake-Tau 3% Down","L");
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

 TH1F* hRatio = new TH1F( "hRatio" ," ; ; #frac{DATA}{MC}" , 20,0,400);
 hRatio->Reset();
 hRatio->SetXTitle("");
 hRatio->SetYTitle("#frac{(DATA)}{MC}");
 hRatio->GetXaxis()->SetTitleSize(0.1);
 hRatio->GetXaxis()->SetTitleOffset(0.80);
 hRatio->GetYaxis()->SetTitleSize(0.1);
 hRatio->GetYaxis()->SetLabelSize(0.09);
 hRatio->GetXaxis()->SetLabelSize(0.09);
 hRatio->GetYaxis()->SetTitleOffset(0.85);
 hRatio->Sumw2();
 hRatio->SetMarkerStyle(kFullCircle);
 hRatio->SetMarkerSize(1.3);
 hRatio->Divide(hMassData,h2MassMadGraph,1.0,1.0);
 hRatio->SetAxisRange(0.5,1.5,"Y");
 hRatio->Draw("PE");

 TH1F* hRatioUp = (TH1F*)hRatio->Clone("hRatioUp");
 hRatioUp->Reset();
 hRatioUp->SetLineColor(kRed);
 hRatioUp->SetLineWidth(3);
 //hRatioUp->SetFillColor(kRed);
 //hRatioUp->SetFillStyle(3003);
 hRatioUp->SetLineStyle(kDashed);
 hRatioUp->Divide(h2MassMadGraphUp,h2MassMadGraph,1.0,1.0);
 hRatioUp->Draw("HISTSAME");

 TH1F* hRatioDown = (TH1F*)hRatio->Clone("hRatioDown");
 hRatioDown->Reset();
 hRatioDown->SetLineColor(kRed);
 //hRatioDown->SetFillColor(kRed);
 //hRatioDown->SetFillStyle(3003);
 hRatioDown->SetLineWidth(3);
 hRatioDown->SetLineStyle(kDotted);
 hRatioDown->Divide(h2MassMadGraphDown,h2MassMadGraph,1.0,1.0);
 hRatioDown->Draw("HISTSAME");

 TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
 line->SetLineStyle(3);
 line->SetLineWidth(1.5);
 line->SetLineColor(kBlack);
 line->Draw("SAME");

 c11->cd();
 c11->SaveAs(title);


}

void plotWshapeAll(){

  plotWshape("<40", "<1", 0, " W-Sdb","Wshape0jetPt2040WsdbMuTau.pdf");
  plotWshape(">40", "<1", 0, " W-Sdb","Wshape0jetPt40WsdbMuTau.pdf");
  plotWshape("<40", "<1", 1, " W-Sdb","Wshape1jetPt2040WsdbMuTau.pdf");
  plotWshape(">40", "<1", 1, " W-Sdb","Wshape1jetPt40WsdbMuTau.pdf");
  plotWshape(">20", "<1", 2, " W-Sdb","Wshape2jetPt20WsdbMuTau.pdf");

  plotWshape("<40", ">1", 0, " t#bar{t}-Sdb","Wshape0jetPt2040TTsdbMuTau.pdf");
  plotWshape(">40", ">1", 0, " t#bar{t}-Sdb","Wshape0jetPt40TTsdbMuTau.pdf");
  plotWshape("<40", ">1", 1, " t#bar{t}-Sdb","Wshape1jetPt2040TTsdbMuTau.pdf");
  plotWshape(">40", ">1", 1, " t#bar{t}-Sdb","Wshape1jetPt40TTsdbMuTau.pdf");
  plotWshape(">20", ">1", 2, " t#bar{t}-Sdb","Wshape2jetPt20TTsdbMuTau.pdf");



}



void plotWshapeElecTau(string ptL2 = ">20", string nJets20BTagged = "<1", float nJets30 = 0, string label = " W-Sdb",
		       TString title = "Wshape0jetPt2040WsdbElecTau.pdf"){


 TFile *fData = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_Approval_thesis/nTupleRun2011-ElecTau-All_run_Open_ElecTauStream.root", "READ");
 TTree *data = (TTree*)fData->Get("outTreePtOrdRaw");

 TFile *fMadGraph = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_Approval_thesis/nTupleW3Jets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root", "READ");
 TTree *madGraph = (TTree*)fMadGraph->Get("outTreePtOrdRaw");

 TFile *fMadGraphDY = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_Approval_thesis/nTupleDYJets-ElecTau-50-madgraph-PUS6_run_Open_ElecTauStream.root", "READ");
 TTree *madGraphDYIncl = (TTree*)fMadGraphDY->Get("outTreePtOrdRaw");

 TFile *fMadGraphIncl = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_Approval_thesis/nTupleWJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root", "READ");
 TTree *madGraphIncl = (TTree*)fMadGraphIncl->Get("outTreePtOrdRaw");

 TFile *fTTIncl = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_Approval_thesis/nTupleTTJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root", "READ");
 TTree *ttIncl = (TTree*)fTTIncl->Get("outTreePtOrdRaw");

 TCut lpt("ptL1>20 && TMath::Abs(etaL1)<2.1"); 
 TCut lID("((TMath::Abs(etaL1)<0.80 && mvaPOGNonTrig>0.925) || (TMath::Abs(etaL1)<1.479 && TMath::Abs(etaL1)>0.80 && mvaPOGNonTrig>0.975) || (TMath::Abs(etaL1)>1.479 && mvaPOGNonTrig>0.985))");
 lpt = lpt && lID;
 TCut tpt(Form("ptL2%s  && TMath::Abs(etaL2)<2.3",ptL2.c_str()));


 TCut tiso("tightestHPSMVAWP>=0 && tightestAntiECutWP>1"); 
 TCut tAiso("tightestHPSMVAWP>=0 && tightestAntiECutWP<1"); 
 TCut ltiso("tightestHPSMVAWP>-99 && tightestAntiECutWP>1"); 
 TCut mtiso("hpsMVA>0.0");
 
 ////// ELEC ISO ///////
 TCut liso("combRelIsoLeg1DBetav2<0.05");
 TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
 TCut lliso("combRelIsoLeg1DBetav2<0.30");
 
 ////// EVENT WISE //////
 TCut lveto("elecFlag==0");
 TCut lAveto("elecFlag==1");
 TCut SS("diTauCharge!=0");
 TCut OS("diTauCharge==0");
 TCut apZ(Form("((MtLeg1Corr)>%f)",70.));
 TCut hltevent("pairIndex<1 && HLTx==1 && (run>=163269 || run==1)");
 TCut hltmatch("HLTmatch==1");

 TCut vbf(Form("nJets30>=%f && nJets20BTagged%s",nJets30,nJets20BTagged.c_str()));
 //TCut vbf("nJets30>1  && MVAvbf>0.5 &&  (ptVeto<30 || isVetoInJets!=1)");
 //TCut vbf("nJets30>1  && Mjj>400 && Deta>3 &&  (ptVeto<30 || isVetoInJets!=1) && nJets30<999");


 TCut sbinPZetaRelInclusive;
 sbinPZetaRelInclusive     = vbf  && lpt && tpt && tiso && liso  && lveto && OS && hltevent && hltmatch;
 TCut sbinPZetaRelAIsoInclusive;
 sbinPZetaRelAIsoInclusive = vbf  && lpt && tpt && tiso && laiso && lveto && OS && hltevent && hltmatch; 

 TCanvas *c11 = new TCanvas("c11","",5,30,650,600);
 c11->SetGrid(0,0);
 c11->SetFillStyle(4000);
 c11->SetFillColor(10);
 c11->SetTicky();
 c11->SetObjectStat(0);

 TPad* pad1 = new TPad("pad1DEta","",0.05,0.22,0.96,0.97);
 TPad* pad2 = new TPad("pad2DEta","",0.05,0.02,0.96,0.20);
 
 pad1->SetFillColor(0);
 pad2->SetFillColor(0);
 pad1->Draw();
 pad2->Draw();
 
 pad1->cd();
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


 TH1F* hMassData = new TH1F("hMassData","",20,0,400);
 hMassData->Sumw2();
 data->Draw("diTauNSVfitMass>>hMassData", sbinPZetaRelInclusive&&apZ);
 hMassData->SetMarkerStyle(kFullCircle);
 hMassData->SetMarkerSize(1.3);

 TH1F* hMassMadGraph = new TH1F("hMassMadGraph"," ; SVfit mass (GeV) ; units",20,0,400);
 hMassMadGraph->Sumw2();
 madGraphIncl->Draw("diTauNSVfitMass>>hMassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 hMassMadGraph->SetLineColor(kRed);




 TH1F* h2MassMadGraph = new TH1F("h2MassMadGraph"," ; SVfit mass (GeV) ; units",20,0,400);
 h2MassMadGraph->SetLineWidth(1);
 h2MassMadGraph->SetFillStyle(3003);
 h2MassMadGraph->SetFillColor(kBlue);
 h2MassMadGraph->Sumw2();
 h2MassMadGraph->GetXaxis()->SetTitleSize(0.055);
 h2MassMadGraph->GetXaxis()->SetTitleOffset(0.80);
 h2MassMadGraph->GetYaxis()->SetTitleSize(0.055);
 h2MassMadGraph->GetYaxis()->SetTitleOffset(0.85);

 TH1F* h2MassMadGraphUp = (TH1F*)h2MassMadGraph->Clone("h2MassMadGraphUp");
 h2MassMadGraphUp->SetFillStyle(0);
 h2MassMadGraphUp->SetLineColor(kRed);
 h2MassMadGraphUp->SetLineWidth(2);
 h2MassMadGraphUp->SetLineStyle(kDashed);
 TH1F* h2MassMadGraphDown = (TH1F*)h2MassMadGraph->Clone("h2MassMadGraphDown");
 h2MassMadGraphDown->SetFillStyle(0);
 h2MassMadGraphDown->SetLineColor(kRed);
 h2MassMadGraphDown->SetLineWidth(2);
 h2MassMadGraphDown->SetLineStyle(kDotted);

 madGraphIncl->Draw("diTauNSVfitMass>>h2MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 madGraphIncl->Draw("diTauNSVfitMass*(1+0.03)>>h2MassMadGraphUp", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 madGraphIncl->Draw("diTauNSVfitMass*(1-0.03)>>h2MassMadGraphDown", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));

 float Wincl = h2MassMadGraph->Integral();
 cout << "W = " << Wincl << endl;

 TH1F* h3MassMadGraph = new TH1F("h3MassMadGraph","",20,0,400);
 h3MassMadGraph->Sumw2();

 TH1F* h3MassMadGraphUp = (TH1F*)h3MassMadGraph->Clone("h3MassMadGraphUp");
 h3MassMadGraphUp->SetFillStyle(0);
 h3MassMadGraphUp->SetLineColor(kRed);
 h3MassMadGraphUp->SetLineWidth(2);
 h3MassMadGraphUp->SetLineStyle(kDashed);
 TH1F* h3MassMadGraphDown = (TH1F*)h3MassMadGraph->Clone("h3MassMadGraphDown");
 h3MassMadGraphDown->SetFillStyle(0);
 h3MassMadGraphDown->SetLineColor(kRed);
 h3MassMadGraphDown->SetLineWidth(2);
 h3MassMadGraphDown->SetLineStyle(kDotted);

 ttIncl->Draw("diTauNSVfitMass>>h3MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 ttIncl->Draw("diTauNSVfitMass*(1+0.03)>>h3MassMadGraphUp", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 ttIncl->Draw("diTauNSVfitMass*(1-0.03)>>h3MassMadGraphDown", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));

 float TTincl = h3MassMadGraph->Integral();
 cout << "TT = " << TTincl << endl;

 TH1F* h4MassMadGraph = new TH1F("h4MassMadGraph","",20,0,400);
 h4MassMadGraph->Sumw2();
 madGraphDYIncl->Draw("diTauNSVfitMass>>h4MassMadGraph", "(sampleWeight*puWeight2*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelInclusive&&apZ));
 float DYincl = h4MassMadGraph->Integral();
 cout << "DY = " << DYincl << endl;

 TH1F* h5MassMadGraph = new TH1F("h5MassMadGraph","",20,0,400);
 data->Draw("diTauNSVfitMass>>h5MassMadGraph", (sbinPZetaRelAIsoInclusive&&apZ));
 float AIsoincl = h5MassMadGraph->Integral();
 cout << "AIso = " << AIsoincl << endl;


 h2MassMadGraph->Add(h3MassMadGraph,1.0);
 h2MassMadGraph->Add(h4MassMadGraph,1.0);

 h2MassMadGraphUp->Add(h3MassMadGraphUp,1.0);
 h2MassMadGraphUp->Add(h4MassMadGraph,1.0);

 h2MassMadGraphDown->Add(h3MassMadGraphDown,1.0);
 h2MassMadGraphDown->Add(h4MassMadGraph,1.0);


 h2MassMadGraph->SetLineColor(kBlue);
 h2MassMadGraph->Scale(1./h2MassMadGraph->Integral());
 h2MassMadGraphUp->Scale(1./h2MassMadGraphUp->Integral());
 h2MassMadGraphDown->Scale(1./h2MassMadGraphDown->Integral());

 hMassData->Scale(1./hMassData->Integral());

 h2MassMadGraph->Draw("EHIST");
 h2MassMadGraphUp->Draw("HISTSAME");
 h2MassMadGraphDown->Draw("HISTSAME");

 hMassData->Draw("PESAME");
 //hMassMadGraph->DrawNormalized("ESAME");

 TLegend* leg = new TLegend(0.50,0.51,0.70,0.85,NULL,"brNDC");
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);
 leg->SetFillColor(10);
 leg->SetTextSize(0.050);
 leg->SetHeader(Form("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV %s}",label.c_str()));

 leg->AddEntry(hMassData,"Data","P");
 leg->AddEntry(h2MassMadGraph,"Background","L");
 leg->AddEntry(h2MassMadGraphUp,"Fake-Tau 3% Up","L");
 leg->AddEntry(h2MassMadGraphDown,"Fake-Tau 3% Down","L");
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

 TH1F* hRatio = new TH1F( "hRatio" ," ; ; #frac{DATA}{MC}" , 20,0,400);
 hRatio->Reset();
 hRatio->SetXTitle("");
 hRatio->SetYTitle("#frac{(DATA)}{MC}");
 hRatio->GetXaxis()->SetTitleSize(0.1);
 hRatio->GetXaxis()->SetTitleOffset(0.80);
 hRatio->GetYaxis()->SetTitleSize(0.1);
 hRatio->GetYaxis()->SetLabelSize(0.09);
 hRatio->GetXaxis()->SetLabelSize(0.09);
 hRatio->GetYaxis()->SetTitleOffset(0.85);
 hRatio->Sumw2();
 hRatio->SetMarkerStyle(kFullCircle);
 hRatio->SetMarkerSize(1.3);
 hRatio->Divide(hMassData,h2MassMadGraph,1.0,1.0);
 hRatio->SetAxisRange(0.5,1.5,"Y");
 hRatio->Draw("PE");

 TH1F* hRatioUp = (TH1F*)hRatio->Clone("hRatioUp");
 hRatioUp->Reset();
 hRatioUp->SetLineColor(kRed);
 hRatioUp->SetLineWidth(3);
 //hRatioUp->SetFillColor(kRed);
 //hRatioUp->SetFillStyle(3003);
 hRatioUp->SetLineStyle(kDashed);
 hRatioUp->Divide(h2MassMadGraphUp,h2MassMadGraph,1.0,1.0);
 hRatioUp->Draw("HISTSAME");

 TH1F* hRatioDown = (TH1F*)hRatio->Clone("hRatioDown");
 hRatioDown->Reset();
 hRatioDown->SetLineColor(kRed);
 //hRatioDown->SetFillColor(kRed);
 //hRatioDown->SetFillStyle(3003);
 hRatioDown->SetLineWidth(3);
 hRatioDown->SetLineStyle(kDotted);
 hRatioDown->Divide(h2MassMadGraphDown,h2MassMadGraph,1.0,1.0);
 hRatioDown->Draw("HISTSAME");

 TF1* line = new TF1("line","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
 line->SetLineStyle(3);
 line->SetLineWidth(1.5);
 line->SetLineColor(kBlack);
 line->Draw("SAME");

 c11->cd();
 c11->SaveAs(title);


}

void plotWshapeElecTauAll(){

  plotWshapeElecTau("<40", "<1", 0, " W-Sdb","Wshape0jetPt2040WsdbElecTau.pdf");
  plotWshapeElecTau(">40", "<1", 0, " W-Sdb","Wshape0jetPt40WsdbElecTau.pdf");
  plotWshapeElecTau("<40", "<1", 1, " W-Sdb","Wshape1jetPt2040WsdbElecTau.pdf");
  plotWshapeElecTau(">40", "<1", 1, " W-Sdb","Wshape1jetPt40WsdbElecTau.pdf");
  plotWshapeElecTau(">20", "<1", 2, " W-Sdb","Wshape2jetPt20WsdbElecTau.pdf");

  plotWshapeElecTau("<40", ">1", 0, " t#bar{t}-Sdb","Wshape0jetPt2040TTsdbElecTau.pdf");
  plotWshapeElecTau(">40", ">1", 0, " t#bar{t}-Sdb","Wshape0jetPt40TTsdbElecTau.pdf");
  plotWshapeElecTau("<40", ">1", 1, " t#bar{t}-Sdb","Wshape1jetPt2040TTsdbElecTau.pdf");
  plotWshapeElecTau(">40", ">1", 1, " t#bar{t}-Sdb","Wshape1jetPt40TTsdbElecTau.pdf");
  plotWshapeElecTau(">20", ">1", 2, " t#bar{t}-Sdb","Wshape2jetPt20TTsdbElecTau.pdf");

}

void plotWtemplateVBF(){



 TFile *fMadGraph = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleW3Jets-MuTau-madgraph-PUS6_run_Open_MuTauStream.root", "READ");
 TTree *madGraph = (TTree*)fMadGraph->Get("outTreePtOrdRaw");

 TCut lpt("ptL1>17 && isPFMuon && isTightMuon");
 TCut tpt("ptL2>20");

 ////// TAU ISO //////
 TCut tiso("tightestHPSMVAWP>=0");
 //TCut tiso("hpsMVA<0.80");     // <-----------------------
 TCut ltiso("tightestHPSMVAWP>-99");
 TCut mtiso("hpsMVA>0.00");

 ////// MU ISO ///////
 TCut liso("combRelIsoLeg1DBetav2<0.10");
 TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
 TCut lliso("combRelIsoLeg1DBetav2<0.30");


 ////// EVENT WISE //////
 TCut lveto("muFlag==0"); //<----------------------------------
 TCut SS("diTauCharge!=0");
 TCut OS("diTauCharge==0");
 TCut apZ(Form("((MtLeg1Corr)>%f)",70.));
 TCut pZ(Form("((MtLeg1Corr)<%f)",40.));
 TCut hltevent("pairIndex<1 && HLTx==1 && ( run>=163269 || run==1)");
 TCut hltmatch("HLTmatch==1");

 TCut vbf("nJets30>1  && MVAvbf>0.5 &&  (ptVeto<30 || isVetoInJets!=1)");
 TCut vbfLoose("nJets30>1  && MVAvbf>-0.30 &&  (ptVeto<30 || isVetoInJets!=1)");
 //TCut vbf("nJets30>1  && Mjj>400 && Deta>3 &&  (ptVeto<30 || isVetoInJets!=1) && nJets30<999");


 TCut sbinPZetaRelInclusive;
 sbinPZetaRelInclusive      = vbf       && lpt && tpt && tiso && liso  && lveto && OS && hltevent && hltmatch;
 TCut sbinPZetaRelInclusive2;
 sbinPZetaRelInclusive2     = vbfLoose  && lpt && tpt && tiso && liso  && lveto && OS && hltevent && hltmatch;
 TCut sbinPZetaRelInclusive3;
 sbinPZetaRelInclusive3     = vbf       && lpt && tpt && mtiso && liso  && lveto && OS && hltevent && hltmatch;
 TCut sbinPZetaRelInclusive4;
 sbinPZetaRelInclusive4     = vbfLoose  && lpt && tpt && mtiso && liso  && lveto && OS && hltevent && hltmatch;

 TCanvas *c11 = new TCanvas("c11","",5,30,650,600);
 c11->SetGrid(0,0);
 c11->SetFillStyle(4000);
 c11->SetFillColor(10);
 c11->SetTicky();
 c11->SetObjectStat(0);

 TH1F* hMass = new TH1F("hMass"," ; SVfit mass (GeV) ; units",10,0,400);
 hMass->Sumw2();
 madGraph->Draw("diTauNSVfitMass>>hMass", sbinPZetaRelInclusive&&pZ);
 //hMass->SetMarkerStyle(kFullCircle);
 //hMass->SetMarkerSize(1.3);
 hMass->GetXaxis()->SetTitleSize(0.055);
 hMass->GetXaxis()->SetTitleOffset(0.80);
 hMass->GetYaxis()->SetTitleSize(0.055);
 hMass->GetYaxis()->SetTitleOffset(0.85);
 hMass->SetLineWidth(2);
 hMass->SetLineColor(kRed);

 TH1F* hMass2 = new TH1F("hMass2","; SVfit mass (GeV) ; units",10,0,400);
 hMass2->Sumw2();
 madGraph->Draw("diTauNSVfitMass>>hMass2", sbinPZetaRelInclusive2&&pZ);
 //hMass2->SetMarkerStyle(kOpenCircle);
 //hMass2->SetMarkerSize(1.3);
 hMass2->GetXaxis()->SetTitleSize(0.055);
 hMass2->GetXaxis()->SetTitleOffset(0.80);
 hMass2->GetYaxis()->SetTitleSize(0.055);
 hMass2->GetYaxis()->SetTitleOffset(0.85);
 hMass2->SetLineWidth(3);
 hMass2->SetLineColor(kBlue);

 TH1F* hMass3 = new TH1F("hMass3","; SVfit mass (GeV) ; units",10,0,400);
 hMass3->Sumw2();
 madGraph->Draw("diTauNSVfitMass>>hMass3", sbinPZetaRelInclusive3&&pZ);
 hMass3->GetXaxis()->SetTitleSize(0.055);
 hMass3->GetXaxis()->SetTitleOffset(0.80);
 hMass3->GetYaxis()->SetTitleSize(0.055);
 hMass3->GetYaxis()->SetTitleOffset(0.85);
 hMass3->SetMarkerStyle(kOpenSquare);
 hMass3->SetMarkerSize(1.3);

 TH1F* hMass4 = new TH1F("hMass4","; SVfit mass (GeV) ; units",10,0,400);
 hMass4->Sumw2();
 madGraph->Draw("diTauNSVfitMass>>hMass4", sbinPZetaRelInclusive4&&pZ);
 hMass4->GetXaxis()->SetTitleSize(0.055);
 hMass4->GetXaxis()->SetTitleOffset(0.80);
 hMass4->GetYaxis()->SetTitleSize(0.055);
 hMass4->GetYaxis()->SetTitleOffset(0.85);
 hMass4->SetMarkerStyle(kFullSquare);
 hMass4->SetMarkerSize(1.3);


 //hMass->DrawNormalized("E");
 //hMass2->DrawNormalized("E");
 hMass3->DrawNormalized("PE");
 hMass4->DrawNormalized("PESAME");

 TLegend* leg = new TLegend(0.50,0.51,0.70,0.85,NULL,"brNDC");
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);
 leg->SetFillColor(10);
 leg->SetTextSize(0.050);
 leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV W+3 jets}");

 leg->AddEntry(hMass4,"Relaxed VBF, relaxed #tau-iso",   "L");
 leg->AddEntry(hMass3,"Full VBF, relaxed #tau-iso","P");
 leg->Draw();



}


void plotBestFit(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  
  TLegend* leg = new TLegend(0.14,0.45,0.41,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetFillColor(0);

  //combine -M MaxLikelihoodFit combined_SM_mH110.txt --stepSize 0.0001 --rMin -20 --rMax 20

  float X[]        = {110,115,120,125,130,135,140,145};
  
  //float obsY[]     = {0.514915 , 0.354472 , 0.217412 ,0.0     ,0.       ,0.     ,0.,     0.};
  float obsY[]     = {0.610128 , 0.380753 , 0.244691 , -0.1626 , -0.49405  ,-0.56743     , -0.8491 ,  -1.3524};

  float expX1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  //float expY1sL[]  = {0.514915 , 0.354472 , 0.217412 ,0.      ,0.       ,0.     ,0.,     0.};
  //float expY1sH[]  = {1.10714  , 1.00143  , 0.9799   ,0.8444  ,0.66112  ,0.7357 , 0.7792, 0.9388};
  float expY1sL[]  = {0.98753 , 0.96267 , 0.9441 , 0.99270,  1.0513,  1.1820,  1.3806,  1.83771};
  float expY1sH[]  = {1.008   , 0.9754  , 0.9523 , 0.9917 ,  1.036,   1.164 ,  1.371, 1.8209};



  //for(int i = 0 ; i < 8 ; i++){
  //expY1sH[i] = TMath::Abs(expY1sH[i]-obsY[i]);
  //expY1sL[i] = TMath::Abs(expY1sL[i]-obsY[i]);
  //}


  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");
  
  TGraphAsymmErrors* observed = new TGraphAsymmErrors(8, X, obsY, expX1sL ,expX1sL , expX1sL, expX1sL);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(8, X, obsY, expX1sL, expX1sL,  expY1sL, expY1sH);

  oneSigma->SetMarkerColor(kBlack);
  oneSigma->SetMarkerStyle(kFullCircle);
  oneSigma->SetFillColor(7);
  oneSigma->SetFillStyle(1001);

  observed->SetMarkerColor(kBlack);
  observed->SetMarkerSize(1.5);
  observed->SetMarkerStyle(kFullCircle);
  observed->SetLineColor(kBlack);
  observed->SetLineWidth(4);

  mg->Add(oneSigma);
  mg->Add(observed);

  mg->Draw("ALP3");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->SetLimits(110,145);
  mg->GetYaxis()->SetTitleOffset(0.85);
  mg->GetXaxis()->SetTitleOffset(0.82);
  mg->SetMinimum(-4.);
  mg->SetMaximum(3);
  mg->GetXaxis()->SetTitle("m_{H} (GeV)");
  mg->GetYaxis()->SetTitle("best fit value of #mu");
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.055);

  TString legend("#tau_{l}#tau_{h}");

  leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=7 TeV}{4.9 fb^{-1}, %s}",legend.Data()));
  leg->SetTextSize(0.05);

  leg->AddEntry(observed,"Observed","P");
  leg->AddEntry(oneSigma,"#pm 1#sigma","F");

  leg->Draw();

  TF1 *line = new TF1("line","1",100,150);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  line->SetLineWidth(3);

  line->Draw("SAME");

}




void plotSVfit(){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TFile* fNew90  = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
  TFile* fNew110 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH110-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew115 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH115-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew120 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH120-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew125 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH125-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew130 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH130-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew135 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH135-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew140 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH140-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew145 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH145-MuTau-powheg-PUS6_run_Open_MuTauStream.root");

  TTree* tNew90  = (TTree*)fNew90->Get("outTreePtOrdRaw");
  TTree* tNew110 = (TTree*)fNew110->Get("outTreePtOrdRaw");
  TTree* tNew115 = (TTree*)fNew115->Get("outTreePtOrdRaw");
  TTree* tNew120 = (TTree*)fNew120->Get("outTreePtOrdRaw");
  TTree* tNew125 = (TTree*)fNew125->Get("outTreePtOrdRaw");
  TTree* tNew130 = (TTree*)fNew130->Get("outTreePtOrdRaw");
  TTree* tNew135 = (TTree*)fNew135->Get("outTreePtOrdRaw");
  TTree* tNew140 = (TTree*)fNew140->Get("outTreePtOrdRaw");
  TTree* tNew145 = (TTree*)fNew145->Get("outTreePtOrdRaw");

  vector<int> mass;
  mass.push_back(91);
  mass.push_back(110);
  mass.push_back(115);
  mass.push_back(120);
  mass.push_back(125);
  mass.push_back(130);
  mass.push_back(135);
  mass.push_back(140);
  mass.push_back(145);

  float X[]        = {91,110,115,120,125,130,135,140,145};
  float Y[]        = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float Y2[]        = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expX1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expX2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.};


  TH1F* mean = new TH1F("meanNew","; generated mass (GeV); mean", 80,  80, 160);
 

  for(unsigned int i = 0 ; i < mass.size() ; i++){
  
    TH1F* h = new TH1F("h","",400,0,400);
    TH1F* h2 = new TH1F("h2","",400,0,400);
 
    cout << i << endl;

    if( i == 0 ){
      tNew90->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew90->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0 && pt1>150)");
      //h->Reset();
    }
    else if( i == 1 ){
      tNew110->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew110->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0 &&  pt1>150)");
      //h->Reset();
    }
    else if( i == 2 ){
      tNew115->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew115->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0 && pt1>150)");
      //h->Reset();
    }
    else if( i == 3 ){
      tNew120->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew120->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }
    else if( i == 4 ){
      tNew125->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew125->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }
    else if( i == 5 ){
      tNew130->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew130->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }
    else if( i == 6 ){
      tNew135->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew135->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }
    else if( i == 7 ){
      tNew140->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew140->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }
    else if( i == 8 ){
      tNew145->Draw("diTauNSVfitMass>>h", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0)");
      tNew145->Draw("diTauNSVfitMass>>h2", "puWeight2*HLTweightTau*HLTweightMu*(isTauLegMatched && ptL1>17 && HLTx && MtLeg1Corr<40 && diTauCharge==0 && combRelIsoLeg1DBetav2<0.1 && muFlag==0  && pt1>150)");
      //h->Reset();
    }

    Double_t xq[1];  // position where to compute the quantiles in [0,1]
    Double_t yq[1]; 
    xq[0] = 0.5; 
    h->GetQuantiles(1,yq,xq);


    X[i]         = mass[i];
    Y[i]         = h->GetMean();
    Y2[i]        = h2->GetMean();
    //Y[i]         = yq[0];
    //Y[i]         = h->GetBinCenter(h->GetMaximumBin());
    expY1sL[i]   = h->GetRMS();
    expY2sL[i]   = h2->GetRMS();

    cout << "Mean ==> " << Y[i] << endl;
    cout << "Mean 2 ==> " << Y2[i] << endl;

    delete h;

  }
    

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");
  TGraphAsymmErrors* observed = new TGraphAsymmErrors(9, X,  Y , expX1sL ,expX1sL , expY1sL, expY1sL);
  observed->SetMarkerColor(kBlue);
  observed->SetMarkerSize(1.6);
  observed->SetMarkerStyle(kFullSquare);
  observed->SetLineColor(kBlue);
  observed->SetLineWidth(4);

  TGraphAsymmErrors* observed2 = new TGraphAsymmErrors(9, X,  Y2 , expX2sL ,expX2sL , expY2sL, expY2sL);
  observed2->SetMarkerColor(kRed);
  observed2->SetMarkerSize(1.6);
  observed2->SetMarkerStyle(kOpenCircle);
  observed2->SetLineColor(kRed);
  observed2->SetLineWidth(3);

  mg->Add(observed);
  mg->Add(observed2);

  //mean->Draw();
  mg->Draw("AP");
  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->SetLimits(60,160);
  mg->GetYaxis()->SetTitleOffset(0.85);
  mg->GetXaxis()->SetTitleOffset(0.82);
  mg->SetMinimum(60);
  mg->SetMaximum(160);
  mg->GetXaxis()->SetTitle("generated mass (GeV)");
  mg->GetYaxis()->SetTitle("reconstructed mass (GeV)");
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.055);

  TLegend* leg = new TLegend(0.55,0.67,0.78,0.82,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Preliminary #sqrt{s}=7 TeV}{Simulation #tau_{#mu}#tau_{h}}");
  leg->SetTextSize(0.05);

  leg->AddEntry(observed,"Mean #pm RMS (inclusive)","PL");
  leg->AddEntry(observed,"Mean #pm RMS (>0 jets)","PL");
  leg->Draw();



}



void plotSVfitComp(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.60,0.47,0.90,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05);

  leg->SetHeader( "#splitline{CMS Preliminary}{#sqrt{s}=7 TeV #tau_{#mu}#tau_{h}}" );

  TFile* fOld = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
  TTree* tOld = (TTree*)fOld->Get("outTreePtOrd");
  TH1F* hOld = new TH1F("hOld"," ; SVfit mass (GeV) ; units",70,0,420);
  hOld->SetLineColor(kBlack);
  hOld->SetLineWidth(3);
  tOld->Draw("diTauNSVfitMass>>hOld","(HLTx && MtLeg1<40 && ptL1>17 && isTauLegMatched)");

  cout << hOld->GetMean() << " " << hOld->GetRMS() << endl;

  //TFile* fNew = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval_thesis/nTupleGGFH125-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TFile* fNew = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval/nTupleSUSYGG120-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TTree* tNew = (TTree*)fNew->Get("outTreePtOrd");
  TH1F* hNew = new TH1F("hNew"," ; SVfit mass (GeV) ; units",70,0,420);
  hNew->SetLineColor(kRed);
  hNew->SetLineWidth(3);
  hNew->SetFillStyle(3003);
  hNew->SetFillColor(kRed);
  tNew->Draw("diTauNSVfitMass>>hNew","(HLTx && MtLeg1<40 && isTauLegMatched && ptL1>17)");

  cout << hNew->GetMean() << " " << hNew->GetRMS() << endl;

  TFile* fNew2 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval/nTupleSUSYGG180-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TTree* tNew2 = (TTree*)fNew2->Get("outTreePtOrd");
  TH1F* hNew2 = new TH1F("hNew2"," ; SVfit mass (GeV) ; units",70,0,420);
  hNew2->SetLineColor(kMagenta);
  hNew2->SetLineWidth(3);
  hNew2->SetFillStyle(3005);
  hNew2->SetFillColor(kMagenta);
  tNew2->Draw("diTauNSVfitMass>>hNew2","(HLTx && MtLeg1<40 && isTauLegMatched && ptL1>17)");

  TFile* fNew3 = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_PreApproval/nTupleSUSYGG250-MuTau-powheg-PUS6_run_Open_MuTauStream.root");
  TTree* tNew3 = (TTree*)fNew3->Get("outTreePtOrd");
  TH1F* hNew3 = new TH1F("hNew3"," ; SVfit mass (GeV) ; units",70,0,420);
  hNew3->SetLineColor(kBlue);
  hNew3->SetLineWidth(3);
  hNew3->SetFillStyle(3004);
  hNew3->SetFillColor(kBlue);
  tNew3->Draw("diTauNSVfitMass>>hNew3","(HLTx && MtLeg1<40 && isTauLegMatched && ptL1>17)");


  hNew->GetYaxis()->SetTitleOffset(0.85);
  hNew->GetXaxis()->SetTitleOffset(0.82);
  hNew->GetXaxis()->SetTitleSize(0.055);
  hNew->GetYaxis()->SetTitleSize(0.055);
  hOld->GetYaxis()->SetTitleOffset(0.85);
  hOld->GetXaxis()->SetTitleOffset(0.82);
  hOld->GetXaxis()->SetTitleSize(0.055);
  hOld->GetYaxis()->SetTitleSize(0.055);

  hOld->DrawNormalized();
  hNew->DrawNormalized("SAME");
  hNew2->DrawNormalized("SAME");
  hNew3->DrawNormalized("SAME");

  leg->AddEntry( hOld, "Z#rightarrow #tau#tau","L");
  leg->AddEntry( hNew, "H(120)#rightarrow #tau#tau","F");
  leg->AddEntry( hNew2, "H(180)#rightarrow #tau#tau","F");
  leg->AddEntry( hNew3, "H(250)#rightarrow #tau#tau","F");
  leg->Draw();

  //gPad->SaveAs("ZTT_SVfit.png");
  //gPad->SaveAs("ZTT_SVfit.pdf");


}
