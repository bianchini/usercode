#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TText.h"
#include "TFrame.h"
#include "TObject.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TList.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TVectorT.h"

#include <sstream> 
#include <strstream> 
#include <vector>
#include <string>
#include <utility>
#include <iostream>

void plot(){

  TFile* f_Zee     = new TFile("../tagAndProbe/testNewWriteFromPAT_DYToEE-PYTHIA-TAUEFF.root");
  TFile* f_Ztautau = new TFile("../tagAndProbe/testNewWriteFromPAT_DYToTauTau-PYTHIA-TAUEFF.root");

  f_Zee->cd();
  TTree* tree_Zee     = (TTree*) gDirectory->Get("etoTauMargLoose90/fitter_tree"); 
  f_Ztautau->cd();
  TTree* tree_Ztautau = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF[45];
  float vyF[45];

  float ZeeAll     = (float)tree_Zee->GetEntries("mcTrue");
  float ZtautauAll = (float)tree_Ztautau->GetEntries("");

  std::cout << ZtautauAll << std::endl;

  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("electronPreIDOutput<=%f &&  mcTrue",-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    std::cout << ZtautauCut << std::endl;
    vxF[i]=ZeeCut/ZeeAll;
    vyF[i]=ZtautauCut/ZtautauAll;
      //h2->Fill(ZeeCut/ZeeAll,ZtautauCut/ZtautauAll);
    std::cout << Form("electronPreIDOutput<=%f &&  mcTrue: ",-0.1+cut) <<ZeeCut/ZeeAll << " --- " << ZtautauCut/ZtautauAll << std::endl;
  }

  h2->SetXTitle("efficiency of mva<X cut on electrons");
  h2->SetYTitle("efficiency of mva<X cut on taus");
  h2->SetAxisRange(0,1.0,"X");
  h2->SetAxisRange(0.95,1.02,"Y");
  h2->Draw();

  TVectorF vx(45,vxF);
  TVectorF vy(45,vyF);
  TGraph* graph = new TGraph(vx,vy);
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(kRed);
  graph->Draw("P");

}

void plotVBTFWP(const string tau_ = "HPS"){


  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogx(1);

  TLegend* leg = new TLegend(0.36,0.15,0.80,0.5,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile* f_Zee     = new TFile("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-TAUEFF.root");
  TFile* f_Ztautau = new TFile("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-TAUEFF.root");


  f_Zee->cd();
  TTree* tree_Zee     = (TTree*) gDirectory->Get( ("tauFakeRateAnalyzer"+tau_+"/tree").c_str() ); 
  f_Ztautau->cd();
  TTree* tree_Ztautau = (TTree*) gDirectory->Get( ("tauFakeRateAnalyzer"+tau_+"/tree").c_str() ); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF[45];
  float vyF[45];

  float vxF_VBTF[6];
  float vyF_VBTF[6];

  float ZeeAll     = (float)tree_Zee->GetEntries("");
  float ZtautauAll = (float)tree_Ztautau->GetEntries("");

  std::cout << ZtautauAll << std::endl;

  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    std::cout << ZtautauCut << std::endl;
    vxF[i]=ZeeCut/ZeeAll;
    vyF[i]=ZtautauCut/ZtautauAll;
      //h2->Fill(ZeeCut/ZeeAll,ZtautauCut/ZtautauAll);
    std::cout << Form("leadPFChargedHadrMva<=%f: ",-0.1+cut) <<ZeeCut/ZeeAll << " --- " << ZtautauCut/ZtautauAll << std::endl;
  }

  h2->SetXTitle("e#rightarrow #tau_{had} fake-rate");
  h2->SetYTitle("#tau_{had} efficiency");
  h2->SetAxisRange(0.02,1.0,"X");
  h2->SetAxisRange(0.95,1.01,"Y");
  h2->Draw();

  float ids[] = {1.0,.95,.90,.85,.80,.70,.60};
  for(int i = 0; i<6; i++){
    float ZeeCut_i     = (float) tree_Zee->GetEntries( Form("matchedID>%f",ids[i]-0.025));
    float ZtautauCut_i = (float) tree_Ztautau->GetEntries(Form("matchedID>%f",ids[i]-0.025));
    vxF_VBTF[i]=ZeeCut_i/ZeeAll;
    vyF_VBTF[i]=ZtautauCut_i/ZtautauAll;
    std::cout << Form("matchedID>%f",ids[i]-0.025) << "    " << ZeeCut_i/ZeeAll << " <-- VBTF Zee ----  VBTF Ztautau ---> " << ZtautauCut_i/ZtautauAll << std::endl;
  }

  TVectorF vx_VBTF(6,vxF_VBTF);
  TVectorF vy_VBTF(6,vyF_VBTF);

  TVectorF vx(45,vxF);
  TVectorF vy(45,vyF);
  TGraph* graph = new TGraph(vx,vy);
  TGraph* graph_VBTF = new TGraph(vx_VBTF,vy_VBTF);
  if( tau_.find("HPS")!=string::npos) graph->SetMarkerStyle(kOpenCircle);
  else  graph->SetMarkerStyle(kOpenSquare);
  graph->SetMarkerSize(1.2);
  if( tau_.find("HPS")!=string::npos) graph->SetMarkerColor(kRed);
  else graph->SetMarkerColor(kBlue);
  graph_VBTF->SetMarkerStyle(kFullStar);
  graph_VBTF->SetMarkerSize(1.8);
  graph_VBTF->SetMarkerColor(kBlack);
  graph->Draw("P");
  graph_VBTF->Draw("P");

  string tau = tau_.find("HPS")!=string::npos ? "HPS" : "Shrinking Cone";
  leg->SetHeader( ("#splitline{Simulation: "+tau+" #tau_{had}-candidates}{passing tau-ID and loose isolation}").c_str() );
  leg->AddEntry(graph,"#splitline{discriminator by #xi^{lch}}{-0.1#leq #xi^{lch}_{cut} #leq1.0}","P");
  leg->AddEntry(graph_VBTF,"#splitline{cut-based discriminator}{WP95,90,85,80,70,60 (ID-only)}","P");
  leg->Draw();
}




void plot2D(string type = "HPS",  float WP = 0.3){

  
  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile* f_Zee     = new TFile("../testNewWriteFromPAT_Zee.root");
  TFile* f_Ztautau = new TFile("../testNewWriteFromPAT_Ztautau.root");

  f_Zee->cd();
  TTree* tree_Zee     = (TTree*) gDirectory->Get("etoTauMargLoose90/fitter_tree"); 
  f_Ztautau->cd();
  TTree* tree_Ztautau = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);
  h2->SetAxisRange(0.16,0.30,"X");
  h2->SetAxisRange(0.99,1.00,"Y");

  float ZeeAll     = (float)tree_Zee->GetEntries("mcTrue");
  float ZtautauAll = (float)tree_Ztautau->GetEntries("");


  float vxF1D[45];
  float vyF1D[45];

  float vxF2D_HoP[61];
  float vyF2D_HoP[61];
  float vxF2D_EoP[61];
  float vyF2D_EoP[61];
  float vxF2D_HoPEoP[961];
  float vyF2D_HoPEoP[961];
  float vxF2D_EMF[61];
  float vyF2D_EMF[61];
  
  std::cout << "MVA" << std::endl;
  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("electronPreIDOutput<=%f &&  mcTrue",-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    vxF1D[i]=ZeeCut/ZeeAll;
    vyF1D[i]=ZtautauCut/ZtautauAll;
  }

  std::cout << "HoP" << std::endl;
  // HoP WP = -0.1
  for(int i = 0; i<=60; i++){
    float cut = 0.02*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("(electronPreIDOutput<=%f || (electronPreIDOutput>%f && hcalEnergy/leadPFChargedHadrCandTrackP>=%f) )  && mcTrue",WP, WP, cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f || (leadPFChargedHadrMva>%f && leadPFChargedHadrHcalEnergy/leadPFChargedHadrTrackP>=%f)",WP, WP,cut));
    vxF2D_HoP[i]=ZeeCut/ZeeAll;
    vyF2D_HoP[i]=ZtautauCut/ZtautauAll;
  }
  std::cout << "EoP" << std::endl;
  // EoP WP = -0.1
  for(int i = 0; i<=60; i++){
    float cut = 0.02*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("(electronPreIDOutput<=%f || (electronPreIDOutput>%f && ecalEnergy/leadPFChargedHadrCandTrackP<=%f) )  && mcTrue",WP, WP, cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f || (leadPFChargedHadrMva>%f && leadPFChargedHadrEcalEnergy/leadPFChargedHadrTrackP<=%f)",WP, WP,cut));
    vxF2D_EoP[i]=ZeeCut/ZeeAll;
    vyF2D_EoP[i]=ZtautauCut/ZtautauAll;
  }
  std::cout << "HoPEoP" << std::endl;
  // HoPEoP WP = -0.1
  for(int i = 0; i<=30; i++){
    float cuti = 0.04*i;
    for(int j = 0; j<=30; j++){
      float cutj = 0.04*j;
      float ZeeCut = (float) tree_Zee->GetEntries(Form("(electronPreIDOutput<=%f || (electronPreIDOutput>%f && (hcalEnergy/leadPFChargedHadrCandTrackP>=%f || ecalEnergy/leadPFChargedHadrCandTrackP<=%f)) )  && mcTrue",WP, WP, cuti, cutj));
      float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f || (leadPFChargedHadrMva>%f && (leadPFChargedHadrHcalEnergy/leadPFChargedHadrTrackP>=%f || leadPFChargedHadrEcalEnergy/leadPFChargedHadrTrackP<=%f))",WP, WP,cuti,cutj));
      vxF2D_HoPEoP[i*31+j]=ZeeCut/ZeeAll;
      vyF2D_HoPEoP[i*31+j]=ZtautauCut/ZtautauAll;
    }
  }
  std::cout << "EMF" << std::endl;
  // EMF WP = -0.1
  for(int i = 0; i<=60; i++){
    float cut = 0.02*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("(electronPreIDOutput<=%f || (electronPreIDOutput>%f && hcalEnergy/(ecalEnergy+hcalEnergy)>%f) )  && mcTrue",WP, WP, cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f || (leadPFChargedHadrMva>%f && leadPFChargedHadrHcalEnergy/(leadPFChargedHadrEcalEnergy+leadPFChargedHadrHcalEnergy)>=%f)",WP, WP,cut));
    vxF2D_EMF[i]=ZeeCut/ZeeAll;
    vyF2D_EMF[i]=ZtautauCut/ZtautauAll;
  }

  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.0);
  h2->SetXTitle("electrons fake-rate");
  h2->SetYTitle("tau-jet efficiency");
  h2->Draw();

  TVectorF vx1D(45,vxF1D);
  TVectorF vy1D(45,vyF1D);
  TGraph* graph1D = new TGraph(vx1D,vy1D);

  TVectorF vx2D_EoP(61,vxF2D_EoP);
  TVectorF vy2D_EoP(61,vyF2D_EoP);
  TGraph* graph2D_EoP = new TGraph(vx2D_EoP,vy2D_EoP);
  TVectorF vx2D_HoP(61,vxF2D_HoP);
  TVectorF vy2D_HoP(61,vyF2D_HoP);
  TGraph* graph2D_HoP = new TGraph(vx2D_HoP,vy2D_HoP);
  TVectorF vx2D_HoPEoP(961,vxF2D_HoPEoP);
  TVectorF vy2D_HoPEoP(961,vyF2D_HoPEoP);
  TGraph* graph2D_HoPEoP = new TGraph(vx2D_HoPEoP,vy2D_HoPEoP);
  TVectorF vx2D_EMF(61,vxF2D_EMF);
  TVectorF vy2D_EMF(61,vyF2D_EMF);
  TGraph* graph2D_EMF = new TGraph(vx2D_EMF,vy2D_EMF);

  graph1D->SetMarkerStyle(21);
  graph1D->SetMarkerSize(1.0);
  graph1D->SetMarkerColor(kBlue);
  graph1D->Draw("P");

  graph2D_EoP->SetMarkerStyle(kFullCircle);
  graph2D_EoP->SetMarkerSize(1.0);
  graph2D_EoP->SetMarkerColor(kGreen);
  graph2D_EoP->Draw("P");
  graph2D_HoP->SetMarkerStyle(kFullSquare);
  graph2D_HoP->SetMarkerSize(1.0);
  graph2D_HoP->SetMarkerColor(kRed);
  graph2D_HoP->Draw("P");
  graph2D_HoPEoP->SetMarkerStyle(kFullTriangleUp);
  graph2D_HoPEoP->SetMarkerSize(1.0);
  graph2D_HoPEoP->SetMarkerColor(5);
  graph2D_HoPEoP->Draw("P");
  graph2D_EMF->SetMarkerStyle(kOpenCircle);
  graph2D_EMF->SetMarkerSize(1.0);
  graph2D_EMF->SetMarkerColor(kBlack);
  graph2D_EMF->Draw("P");

  leg->SetHeader((type+" passing tau-ID and loose iso").c_str());
  leg->AddEntry(graph1D,"(#epsilon_{e},#epsilon_{#tau}) for mva<X, -0.1<X<1.0","P");
  leg->AddEntry(graph2D_HoP,Form("(#epsilon_{e},#epsilon_{#tau}) for mva<%.2f || (mva>%.2f && H/p>X), 0<X<1.2",WP,WP),"P");
  leg->AddEntry(graph2D_EoP,Form("(#epsilon_{e},#epsilon_{#tau}) for mva<%.2f || (mva>%.2f && E/p<X), 0<X<1.2 ",WP,WP),"P");
  leg->AddEntry(graph2D_HoPEoP,Form("(#epsilon_{e},#epsilon_{#tau}) for mva<%.2f || (mva>-%.2f && (H/p>X||E/p<Y)), 0<X,Y<1.2",WP,WP),"P");
  leg->AddEntry(graph2D_EMF,Form("(#epsilon_{e},#epsilon_{#tau}) for mva<%.2f || (mva>%.2f && H/(E+H)>X), 0<X<1.2 ",WP,WP),"P");
  leg->Draw();

}


void plotMVALeadVsLeadCharg(){

  TFile* f_Ztautau = new TFile("../");
  f_Ztautau->cd();
  TTree* tree_Ztautau = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 
  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF[45];
  float vyF[45];

  float ZtautauAll = (float)tree_Ztautau->GetEntries("");

  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZtautauCutLeadPi    = (float) tree_Ztautau->GetEntries(Form("leadPFCandMva<=%f",-0.1+cut));
    float ZtautauCutLeadCh    = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    vxF[i]=ZtautauCutLeadCh/ZtautauAll;
    vyF[i]=ZtautauCutLeadPi/ZtautauAll;
  }

  h2->SetXTitle("efficiency of mva<X cut on lead charged");
  h2->SetYTitle("efficiency of mva<X cut on lead candidate");
  h2->SetAxisRange(0.95,1.02,"X");
  h2->SetAxisRange(0.95,1.02,"Y");
  
  h2->Draw();

  TVectorF vx(45,vxF);
  TVectorF vy(45,vyF);
  TGraph* graph = new TGraph(vx,vy);
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kRed);  
  graph->SetMarkerSize(1.0);
  graph->Draw("P");

  TF1* line = new TF1("line","x",0.9,1.2);
  line->Draw("SAME");

}


void plotMVALeadVsMVAChargROC(){

  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile* f_Zee     = new TFile("../testNewWriteFromPAT_Zee.root");
  TFile* f_Ztautau = new TFile("../testNewWriteFromPAT_Ztautau.root");

  f_Zee->cd();
  TTree* tree_Zee     = (TTree*) gDirectory->Get("etoTauSCMarg90/fitter_tree"); 
  f_Ztautau->cd();
  TTree* tree_Ztautau = (TTree*) gDirectory->Get("tauFakeRateAnalyzerSHC/tree"); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF_L[45];
  float vyF_L[45];
  float vxF_C[45];
  float vyF_C[45];

  float ZeeAll     = (float)tree_Zee->GetEntries("mcTrue");
  float ZtautauAll = (float)tree_Ztautau->GetEntries("");

  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("electronPreIDOutput<=%f &&  mcTrue",-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFChargedHadrMva<=%f",-0.1+cut));
    vxF_C[i]=ZeeCut/ZeeAll;
    vyF_C[i]=ZtautauCut/ZtautauAll;
  }
  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee->GetEntries(Form("electronPreIDOutputLeadPFCand<=%f &&  mcTrue",-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau->GetEntries(Form("leadPFCandMva<=%f",-0.1+cut));
    vxF_L[i]=ZeeCut/ZeeAll;
    vyF_L[i]=ZtautauCut/ZtautauAll;
  }

  h2->SetXTitle("electrons fake-rate");
  h2->SetYTitle("tau-jet efficiency");
  h2->SetAxisRange(0,1.0,"X");
  h2->SetAxisRange(0.95,1.02,"Y");
  h2->Draw();

  TVectorF vx_L(45,vxF_L);
  TVectorF vy_L(45,vyF_L);
  TGraph* graph_L = new TGraph(vx_L,vy_L);
  TVectorF vx_C(45,vxF_C);
  TVectorF vy_C(45,vyF_C);
  TGraph* graph_C = new TGraph(vx_C,vy_C);
  graph_C->SetMarkerStyle(kFullCircle);
  graph_C->SetMarkerSize(1.0);
  graph_C->SetMarkerColor(kRed);
  graph_C->Draw("P");
  graph_L->SetMarkerStyle(kFullSquare);
  graph_L->SetMarkerSize(1.0);
  graph_L->SetMarkerColor(kBlue);
  graph_L->Draw("P");

  leg->SetHeader("PF-Shrinking Cone passing tau-ID and loose iso");
  leg->AddEntry(graph_C,"discr by mva_e_pi of LeadChargHadr","P");
  leg->AddEntry(graph_L,"discr by mva_e_pi of LeadPFCand","P");
  leg->Draw();

}



void plotSHCvsHPS(const string pion_ = "ch"){

  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogx(1);

  TLegend* leg = new TLegend(0.36,0.15,0.80,0.5,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);

  TFile* f_Zee     = new TFile("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-TAUEFF.root");
  TFile* f_Ztautau = new TFile("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-TAUEFF.root");


  f_Zee->cd();
  TTree* tree_Zee_SHC     = (TTree*) gDirectory->Get("tauFakeRateAnalyzerSHC/tree"); 
  TTree* tree_Zee_HPS     = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 
  f_Ztautau->cd();
  TTree* tree_Ztautau_SHC = (TTree*) gDirectory->Get("tauFakeRateAnalyzerSHC/tree"); 
  TTree* tree_Ztautau_HPS = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF_L[45];
  float vyF_L[45];
  float vxF_C[45];
  float vyF_C[45];

  float ZeeAll_SHC     = (float)tree_Zee_SHC->GetEntries("leadPFChargedHadrTrackPt>3");
  float ZtautauAll_SHC = (float)tree_Ztautau_SHC->GetEntries("leadPFChargedHadrTrackPt>3");
  float ZeeAll_HPS     = (float)tree_Zee_HPS->GetEntries("leadPFChargedHadrTrackPt>3");
  float ZtautauAll_HPS = (float)tree_Ztautau_HPS->GetEntries("leadPFChargedHadrTrackPt>3");

  string pion = pion_.find("ch")!=string::npos ? "leadPFChargedHadrMva" : "leadPFCandMva"; 

  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float)     tree_Zee_SHC->GetEntries(Form("%s<=%f && leadPFChargedHadrTrackPt>3",pion.c_str(),-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau_SHC->GetEntries(Form("%s<=%f  && leadPFChargedHadrTrackPt>3",pion.c_str(),-0.1+cut));
    vxF_C[i]=ZeeCut/ZeeAll_SHC;
    vyF_C[i]=ZtautauCut/ZtautauAll_SHC;
  }
  for(int i = 0; i<=44; i++){
    float cut = 0.025*i;
    float ZeeCut = (float) tree_Zee_HPS->GetEntries(Form("%s<=%f && leadPFChargedHadrTrackPt>3",pion.c_str(),-0.1+cut));
    float ZtautauCut = (float) tree_Ztautau_HPS->GetEntries(Form("%s<=%f && leadPFChargedHadrTrackPt>3",pion.c_str(),-0.1+cut));
    vxF_L[i]=ZeeCut/ZeeAll_HPS;
    vyF_L[i]=ZtautauCut/ZtautauAll_HPS;
  }

  h2->SetXTitle("e#rightarrow#tau_{had} fake-rate");
  h2->SetYTitle("#tau_{had} efficiency");
  h2->SetAxisRange(0.02,1.0,"X");
  h2->SetAxisRange(0.95,1.01,"Y");
  h2->Draw();

  TVectorF vx_L(45,vxF_L);
  TVectorF vy_L(45,vyF_L);
  TGraph* graph_L = new TGraph(vx_L,vy_L);
  TVectorF vx_C(45,vxF_C);
  TVectorF vy_C(45,vyF_C);
  TGraph* graph_C = new TGraph(vx_C,vy_C);
  graph_C->SetMarkerStyle(kOpenSquare);
  graph_C->SetMarkerSize(1.2);
  graph_C->SetMarkerColor(kBlue);
  graph_C->Draw("P");
  graph_L->SetMarkerStyle(kOpenCircle);
  graph_L->SetMarkerSize(1.2);
  graph_L->SetMarkerColor(kRed);
  graph_L->Draw("P");

  string xi = pion_.find("ch")!=string::npos ? "#xi^{lch}" : "#xi^{lh}" ;
  leg->SetHeader(("Simulation: discriminator by "+xi).c_str());
  leg->AddEntry(graph_L,"#splitline{HPS #tau_{had}-candidates}{passing tau-ID and loose isolation}","P");
  leg->AddEntry(graph_C,"#splitline{Shrinking Cone #tau_{had}-candidates}{passing tau-ID and loose isolation}","P");
  leg->Draw();

}




void plotDiffCuts(){

  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile* f_Zee     = new TFile("../testNewWriteFromPAT_Zee_debug.root");
  TFile* f_Ztautau = new TFile("../testNewWriteFromPAT_Ztautau_debug.root");

  f_Zee->cd();
  TTree* tree_Zee_SHC     = (TTree*) gDirectory->Get("tauFakeRateAnalyzerSHC/tree"); 
  TTree* tree_Zee_HPS     = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 
  f_Ztautau->cd();
  TTree* tree_Ztautau_SHC = (TTree*) gDirectory->Get("tauFakeRateAnalyzerSHC/tree"); 
  TTree* tree_Ztautau_HPS = (TTree*) gDirectory->Get("tauFakeRateAnalyzerHPS/tree"); 

  TH2F* h2 = new TH2F("h2","",220,0,1.1,220,0,1.1);

  float vxF_L[45];
  float vyF_L[45];
  float vxF_C[45];
  float vyF_C[45];

  float ZeeAll_SHC     = (float)tree_Zee_SHC->GetEntries("");
  float ZtautauAll_SHC = (float)tree_Ztautau_SHC->GetEntries("");
  float ZeeAll_HPS     = (float)tree_Zee_HPS->GetEntries("");
  float ZtautauAll_HPS = (float)tree_Ztautau_HPS->GetEntries("");

  for(int i = 0; i<=44; i++){
    float cut = 0.005*i;
    float ZeeCut = (float)     tree_Zee_HPS->GetEntries(Form("leadPFChargedHadrMva<=-0.1 && fbrem/leadPFChargedHadrTrackP<%f",cut));
    float ZtautauCut = (float) tree_Ztautau_HPS->GetEntries(Form("leadPFChargedHadrMva<=-0.1 && fbrem/leadPFChargedHadrTrackP<%f",cut));
    cout << Form("leadPFChargedHadrMva<=-0.1 && fbrem/leadPFChargedHadrTrackP<%f",cut) << endl;
    cout << ZeeCut/ZeeAll_HPS <<" ---- "<<ZtautauCut/ZtautauAll_HPS  << endl;
    vxF_C[i]=ZeeCut/ZeeAll_HPS;
    vyF_C[i]=ZtautauCut/ZtautauAll_HPS;
  }
  for(int i = 0; i<=44; i++){
    float cut = 0.005*i;
    float ZeeCut = (float) tree_Zee_HPS->GetEntries(Form("leadPFChargedHadrMva<=%f",cut));
    float ZtautauCut = (float) tree_Ztautau_HPS->GetEntries(Form("leadPFChargedHadrMva<=%f",cut));
    vxF_L[i]=ZeeCut/ZeeAll_HPS;
    vyF_L[i]=ZtautauCut/ZtautauAll_HPS;
  }

  h2->SetXTitle("electrons fake-rate");
  h2->SetYTitle("tau-jet efficiency");
  h2->SetAxisRange(0,1.0,"X");
  h2->SetAxisRange(0.70,1.02,"Y");
  h2->Draw();

  TVectorF vx_L(45,vxF_L);
  TVectorF vy_L(45,vyF_L);
  TGraph* graph_L = new TGraph(vx_L,vy_L);
  TVectorF vx_C(45,vxF_C);
  TVectorF vy_C(45,vyF_C);
  TGraph* graph_C = new TGraph(vx_C,vy_C);
  graph_C->SetMarkerStyle(kFullCircle);
  graph_C->SetMarkerSize(1.0);
  graph_C->SetMarkerColor(kRed);
  graph_C->Draw("P");
  graph_L->SetMarkerStyle(kFullSquare);
  graph_L->SetMarkerSize(1.0);
  graph_L->SetMarkerColor(kBlue);
  graph_L->Draw("P");

  leg->SetHeader("HPS, passing tau-ID and loose iso");
  leg->AddEntry(graph_C,"HPC: discr by mva_e_pi and fbrem<10%","P");
  leg->AddEntry(graph_L,"HPS: discr by mva_e_pi","P");
  leg->Draw();

}
