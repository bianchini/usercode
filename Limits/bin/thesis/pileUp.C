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

  TH1F* h = new TH1F("h","; number of pile-up interactions; weight",35,0,35);
  h->SetMarkerStyle(kFullSquare);
  h->SetMarkerColor(kRed);
  h->SetMarkerSize(1.5);
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleOffset(0.80);
  for(int i = 1; i < 54; i++) h->SetBinContent(i, pileupWeight2(i) );

  h->Draw("P");
  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV}");
  leg->AddEntry(h,"MC pu-weight","P");
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
  
  TFile f("../../../TauTauStudies/test/Macro/pileUp/Run2011PileUp.root","READ");
  
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

  TFile* f = new TFile("../batch/nTupleDYJets-MuTau-50-madgraph-PUS6_run_Open_MuTauStream.root");
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
  TCut cutNumM("isTightMuon && isPFMuon && muFlag==0 && diTauCharge==0 && HLTmatch && HLTx && isTauLegMatched && TMath::Abs(genTauEta)<2.3 && tightestHPSMVAWP>=1");
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
  hNum3->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Z#rightarrow#tau#tau MC #sqrt{s}=7 TeV}");
  leg->AddEntry(hNum, "Loose","P");
  leg->AddEntry(hNum2,"Medium","P");
  leg->AddEntry(hNum3,"Tight","P");
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

  TFile* f = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStreamFall11_04May2012_Approval/nTupleRun2011-MuTau-All_run_Open_MuTauStream.root");
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
  TCut cutNumM("pairIndex<1 && isTightMuon && isPFMuon && muFlag==0 && diTauCharge!=0 && HLTmatch && HLTx && MtLeg1Corr<40 && nJets20BTagged<1 && combRelIsoLeg1DBetav2>0.30 && tightestHPSMVAWP>=1");
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
  hNum3->Draw("PESAME");

  leg->SetHeader("#splitline{CMS Preliminary 2011}{Data #sqrt{s}=7 TeV}");
  leg->AddEntry(hNum, "Loose","P");
  leg->AddEntry(hNum2,"Medium","P");
  leg->AddEntry(hNum3,"Tight","P");
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
