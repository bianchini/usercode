#include <cstdlib>
#include <iostream> 
#include <fstream>
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
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"


void bbb( TH1F* hin, TH1F* hout_Down, TH1F* hout_Up, int bin){

  float bin_down = hin->GetBinContent( bin ) - hin->GetBinError( bin );
  float bin_up   = hin->GetBinContent( bin ) + hin->GetBinError( bin );

  hout_Down->Reset();
  hout_Down->Add(hin,1.0);
  hout_Down->SetBinContent(bin, TMath::Max(bin_down,0.01));

  hout_Up->Reset();
  hout_Up->Add(hin,1.0);
  hout_Up->SetBinContent  (bin, TMath::Max(bin_up,0.01));

  return;

}


void produce( TString fname = "SL_VType2", string cut = "type==0", TString category = "cat0", float fact1 = 0.02, float fact2 = 0., float lumiScale = 20./12.1,
	      int newvar = 1, int nBins=6){


  string version = "_v5";
  //string version = "";

  cout << "Doing version " << version << " and category " << category << endl;

  float scaleTTH      = 1.0;
  float scaleTTJetsLF = 1.0;
  float scaleTTJetsHF = 1.0;

  TFile* fout = new TFile(fname+".root","UPDATE");
  TDirectory* dir =  fout->GetDirectory( fname+"_"+category); 
  if( !dir) dir = fout->mkdir( fname+"_"+category ) ;

  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units", nBins, 0,1 );
 
  TString var("");
  if(newvar)
    var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, fact2, fact2));
  else
    var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));

  string basecut = cut;
  TString ttjets = (string(fname.Data())).find("SL")!=string::npos ? "TTJetsSemiLept" : "TTJetsFullLept";

  ////////////////////////////// TTH
  
  // NOMINAL
  TFile* f_TTH125 = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_TTH125.root");
  TTree* t_TTH125 = (TTree*)f_TTH125->Get("tree");
  TH1F* h_TTH125 = (TH1F*)h->Clone("h_TTH125");
  h_TTH125->Reset();
  h_TTH125->Sumw2();
  t_TTH125->Draw(var+">>h_TTH125",   "weight"*TCut(cut.c_str()) );
  h_TTH125->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125->Write("TTH125", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    TH1F* h_TTH125_b_up   = (TH1F*)h->Clone(Form("h_TTH125_%d_Up",bin));
    TH1F* h_TTH125_b_down = (TH1F*)h->Clone(Form("h_TTH125_%d_Down",bin));
    bbb( h_TTH125, h_TTH125_b_up, h_TTH125_b_down, bin);

    dir->cd();
    h_TTH125_b_up  ->Write(Form("TTH125_TTH125%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTH125_b_down->Write(Form("TTH125_TTH125%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTH125_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp"+version+"_TTH125.root");
  TTree* t_TTH125_JECUp = (TTree*)f_TTH125_JECUp->Get("tree");
  TH1F* h_TTH125_JECUp = (TH1F*)h->Clone("h_TTH125_JECUp");
  h_TTH125_JECUp->Reset();
  h_TTH125_JECUp->Sumw2();
  t_TTH125_JECUp->Draw(var+">>h_TTH125_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECUp->Write("TTH125_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTH125_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown"+version+"_TTH125.root");
  TTree* t_TTH125_JECDown = (TTree*)f_TTH125_JECDown->Get("tree");
  TH1F* h_TTH125_JECDown = (TH1F*)h->Clone("h_TTH125_JECDown");
  h_TTH125_JECDown->Reset();
  h_TTH125_JECDown->Sumw2();
  t_TTH125_JECDown->Draw(var+">>h_TTH125_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECDown->Write("TTH125_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTH125_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp"+version+"_TTH125.root");
  TTree* t_TTH125_csvUp = (TTree*)f_TTH125_csvUp->Get("tree");
  TH1F* h_TTH125_csvUp = (TH1F*)h->Clone("h_TTH125_csvUp");
  h_TTH125_csvUp->Reset();
  h_TTH125_csvUp->Sumw2();
  t_TTH125_csvUp->Draw(var+">>h_TTH125_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvUp->Write("TTH125_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTH125_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown"+version+"_TTH125.root");
  TTree* t_TTH125_csvDown = (TTree*)f_TTH125_csvDown->Get("tree");
  TH1F* h_TTH125_csvDown = (TH1F*)h->Clone("h_TTH125_csvDown");
  h_TTH125_csvDown->Reset();
  h_TTH125_csvDown->Sumw2();
  t_TTH125_csvDown->Draw(var+">>h_TTH125_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvDown->Write("TTH125_csvDown", TObject::kOverwrite);

  ////////////////////////////// TTJets: HF

  cut = basecut+"&&nSimBs>2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HF = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF = (TTree*)f_TTJetsSemiLept_HF->Get("tree");
  TH1F* h_TTJetsSemiLept_HF = (TH1F*)h->Clone("h_TTJetsSemiLept_HF");
  h_TTJetsSemiLept_HF->Reset();
  h_TTJetsSemiLept_HF->Sumw2();
  t_TTJetsSemiLept_HF->Draw(var+">>h_TTJetsSemiLept_HF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF->Write("TTJetsHF", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HF->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HF_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HF_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HF_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HF_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HF, h_TTJetsSemiLept_HF_b_up, h_TTJetsSemiLept_HF_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HF_b_up  ->Write(Form("TTJetsHF_TTJetsHF%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HF_b_down->Write(Form("TTJetsHF_TTJetsHF%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HF_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECUp = (TTree*)f_TTJetsSemiLept_HF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECUp");
  h_TTJetsSemiLept_HF_JECUp->Reset();
  h_TTJetsSemiLept_HF_JECUp->Sumw2();
  t_TTJetsSemiLept_HF_JECUp->Draw(var+">>h_TTJetsSemiLept_HF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_JECUp->Write("TTJetsHF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HF_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECDown = (TTree*)f_TTJetsSemiLept_HF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECDown");
  h_TTJetsSemiLept_HF_JECDown->Reset();
  h_TTJetsSemiLept_HF_JECDown->Sumw2();
  t_TTJetsSemiLept_HF_JECDown->Draw(var+">>h_TTJetsSemiLept_HF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_JECDown->Write("TTJetsHF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HF_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvUp = (TTree*)f_TTJetsSemiLept_HF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvUp");
  h_TTJetsSemiLept_HF_csvUp->Reset();
  h_TTJetsSemiLept_HF_csvUp->Sumw2();
  t_TTJetsSemiLept_HF_csvUp->Draw(var+">>h_TTJetsSemiLept_HF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_csvUp->Write("TTJetsHF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HF_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvDown = (TTree*)f_TTJetsSemiLept_HF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvDown");
  h_TTJetsSemiLept_HF_csvDown->Reset();
  h_TTJetsSemiLept_HF_csvDown->Sumw2();
  t_TTJetsSemiLept_HF_csvDown->Draw(var+">>h_TTJetsSemiLept_HF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_csvDown->Write("TTJetsHF_csvDown", TObject::kOverwrite);


  ////////////////////////////// TTJets: LF

  cut = basecut+"&&nSimBs==2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_LF = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF = (TTree*)f_TTJetsSemiLept_LF->Get("tree");
  TH1F* h_TTJetsSemiLept_LF = (TH1F*)h->Clone("h_TTJetsSemiLept_LF");
  h_TTJetsSemiLept_LF->Reset();
  h_TTJetsSemiLept_LF->Sumw2();
  t_TTJetsSemiLept_LF->Draw(var+">>h_TTJetsSemiLept_LF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF->Write("TTJetsLF", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_LF_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_LF_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Down",bin));
    bbb( h_TTJetsSemiLept_LF, h_TTJetsSemiLept_LF_b_up, h_TTJetsSemiLept_LF_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_LF_b_up  ->Write(Form("TTJetsLF_TTJetsLF%sbin%dUp",  category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_LF_b_down->Write(Form("TTJetsLF_TTJetsLF%sbin%dDown",category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_LF_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECUp = (TTree*)f_TTJetsSemiLept_LF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECUp");
  h_TTJetsSemiLept_LF_JECUp->Reset();
  h_TTJetsSemiLept_LF_JECUp->Sumw2();
  t_TTJetsSemiLept_LF_JECUp->Draw(var+">>h_TTJetsSemiLept_LF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECUp->Write("TTJetsLF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_LF_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECDown = (TTree*)f_TTJetsSemiLept_LF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECDown");
  h_TTJetsSemiLept_LF_JECDown->Reset();
  h_TTJetsSemiLept_LF_JECDown->Sumw2();
  t_TTJetsSemiLept_LF_JECDown->Draw(var+">>h_TTJetsSemiLept_LF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECDown->Write("TTJetsLF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_LF_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvUp = (TTree*)f_TTJetsSemiLept_LF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvUp");
  h_TTJetsSemiLept_LF_csvUp->Reset();
  h_TTJetsSemiLept_LF_csvUp->Sumw2();
  t_TTJetsSemiLept_LF_csvUp->Draw(var+">>h_TTJetsSemiLept_LF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvUp->Write("TTJetsLF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_LF_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvDown = (TTree*)f_TTJetsSemiLept_LF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvDown");
  h_TTJetsSemiLept_LF_csvDown->Reset();
  h_TTJetsSemiLept_LF_csvDown->Sumw2();
  t_TTJetsSemiLept_LF_csvDown->Draw(var+">>h_TTJetsSemiLept_LF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvDown->Write("TTJetsLF_csvDown", TObject::kOverwrite);



  ////////////////////////////// DIBOSON
  cut = basecut;
  
  // NOMINAL
  TFile* f_DiBoson = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_DiBoson.root");
  TTree* t_DiBoson = (TTree*)f_DiBoson->Get("tree");
  TH1F* h_DiBoson = (TH1F*)h->Clone("h_DiBoson");
  h_DiBoson->Reset();
  h_DiBoson->Sumw2();
  t_DiBoson->Draw(var+">>h_DiBoson",   "weight"*TCut(cut.c_str()) );
  h_DiBoson->Scale(lumiScale);
  dir->cd();
  h_DiBoson->Write("DiBoson", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_DiBoson->GetNbinsX(); bin++ ){
    TH1F* h_DiBoson_b_up   = (TH1F*)h->Clone(Form("h_DiBoson_%d_Up",bin));
    TH1F* h_DiBoson_b_down = (TH1F*)h->Clone(Form("h_DiBoson_%d_Down",bin));
    bbb( h_DiBoson, h_DiBoson_b_up, h_DiBoson_b_down, bin);

    dir->cd();
    h_DiBoson_b_up  ->Write(Form("DiBoson_DiBoson%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_DiBoson_b_down->Write(Form("DiBoson_DiBoson%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }


  ////////////////////////////// SINGLE T
  
  // NOMINAL
  TFile* f_SingleT = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_SingleT.root");
  TTree* t_SingleT = (TTree*)f_SingleT->Get("tree");
  TH1F* h_SingleT = (TH1F*)h->Clone("h_SingleT");
  h_SingleT->Reset();
  h_SingleT->Sumw2();
  t_SingleT->Draw(var+">>h_SingleT",   "weight"*TCut(cut.c_str()) );
  h_SingleT->Scale(lumiScale);
  dir->cd();
  h_SingleT->Write("SingleT", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    TH1F* h_SingleT_b_up   = (TH1F*)h->Clone(Form("h_SingleT_%d_Up",bin));
    TH1F* h_SingleT_b_down = (TH1F*)h->Clone(Form("h_SingleT_%d_Down",bin));
    bbb( h_SingleT, h_SingleT_b_up, h_SingleT_b_down, bin);

    dir->cd();
    h_SingleT_b_up  ->Write(Form("SingleT_SingleT%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_SingleT_b_down->Write(Form("SingleT_SingleT%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_TTV = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_TTV.root");
  TTree* t_TTV = (TTree*)f_TTV->Get("tree");
  TH1F* h_TTV = (TH1F*)h->Clone("h_TTV");
  h_TTV->Reset();
  h_TTV->Sumw2();
  t_TTV->Draw(var+">>h_TTV",   "weight"*TCut(cut.c_str()) );
  h_TTV->Scale(lumiScale);
  dir->cd();
  h_TTV->Write("TTV", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    TH1F* h_TTV_b_up   = (TH1F*)h->Clone(Form("h_TTV_%d_Up",bin));
    TH1F* h_TTV_b_down = (TH1F*)h->Clone(Form("h_TTV_%d_Down",bin));
    bbb( h_TTV, h_TTV_b_up, h_TTV_b_down, bin);

    dir->cd();
    h_TTV_b_up  ->Write(Form("TTV_TTV%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTV_b_down->Write(Form("TTV_TTV%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_EWK = TFile::Open("MEAnalysis_"+fname+"_nominal"+version+"_EWK.root");
  TTree* t_EWK = (TTree*)f_EWK->Get("tree");
  TH1F* h_EWK = (TH1F*)h->Clone("h_EWK");
  h_EWK->Reset();
  h_EWK->Sumw2();
  t_EWK->Draw(var+">>h_EWK",   "weight"*TCut(cut.c_str()) );
  h_EWK->Scale(lumiScale);
  dir->cd();
  h_EWK->Write("EWK", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_EWK->GetNbinsX(); bin++ ){
    TH1F* h_EWK_b_up   = (TH1F*)h->Clone(Form("h_EWK_%d_Up",bin));
    TH1F* h_EWK_b_down = (TH1F*)h->Clone(Form("h_EWK_%d_Down",bin));
    bbb( h_EWK, h_EWK_b_up, h_EWK_b_down, bin);

    dir->cd();
    h_EWK_b_up  ->Write(Form("EWK_EWK%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_EWK_b_down->Write(Form("EWK_EWK%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }




  cout << "TTH125:      "  << h_TTH125->Integral()         << endl;
  cout << "TTJets (HF): "  << h_TTJetsSemiLept_HF->Integral() << endl;
  cout << "TTJets (LF): "  << h_TTJetsSemiLept_LF->Integral() << endl;
  cout << "SingleT:     "  << h_SingleT->Integral() << endl;
  cout << "TTV:         "  << h_TTV->Integral() << endl;
  cout << "EWK:         "  << h_EWK->Integral() << endl;
  cout << "DiBoson:     "  << h_DiBoson->Integral() << endl;

  float rTTH125    = h_TTH125->Integral() ;
  float rTTJets_HF = h_TTJetsSemiLept_HF->Integral();
  float rTTJets_LF = h_TTJetsSemiLept_LF->Integral();
  float rSingleT   = h_SingleT->Integral() ;
  float rTTV       = h_TTV->Integral() ;
  float rEWK       = h_EWK->Integral() ;
  float rDiBoson   = h_DiBoson->Integral() ;

  float observation = int(rTTH125 + rTTJets_HF + rTTJets_LF  + rTTV + rSingleT /*+ rEWK + rDiBoson*/ );

  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_TTH125,1.0);
  h_data->Add(h_TTJetsSemiLept_HF,1.0);
  h_data->Add(h_TTJetsSemiLept_LF,1.0);
  h_data->Add(h_TTV);
  h_data->Add(h_SingleT);
  //h_data->Add(h_EWK);
  //h_data->Add(h_DiBoson);

  h_data->Scale(int(observation)/h_data->Integral());
  dir->cd();
  h_data->Write("data_obs", TObject::kOverwrite);

  fout->Close();


  ofstream out(Form("%s.txt",(fname+"_"+category).Data()));
  out.precision(8);
  string longspace  = "              ";
  string shortspace = "          ";
  string null = "";
  string line1("imax 1"); out << line1 << endl;
  string line2("jmax *"); out << line2 << endl;
  string line3("kmax *"); out << line3 << endl;
  string line4(Form("shapes *  *    %s.root  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC", fname.Data())); out << line4 << endl;
  out<<endl;
  string line5(Form("observation %.0f", observation )); out << line5 << endl;
  out<<endl;


  // PROC AND NORM
  string line("bin                         ");
  if( rTTH125>0 )    line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_LF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTV>0 )       line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rSingleT>0 )   line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";
  if( rTTH125>0 )    line += "TTH125      ";
  if( rTTJets_HF>0 ) line += "TTJetsHF    ";
  if( rTTJets_LF>0 ) line += "TTJetsLF    ";
  if( rTTV>0 )       line += "TTV      ";
  if( rSingleT>0 )   line += "SingleT     ";
  out<<line;
  out<<endl;
  line = "process                       ";
  if( rTTH125>0 )    line += "0          ";
  if( rTTJets_HF>0 ) line += "1          ";
  if( rTTJets_LF>0 ) line += "2          ";
  if( rTTV>0 )       line += "3          ";
  if( rSingleT>0 )   line += "4          ";
  out<<line;
  out<<endl;
  line = "rate                        ";
  if( rTTH125>0 )    line += string(Form("%.4f     ", rTTH125 ));
  if( rTTJets_HF>0 ) line += string(Form("%.4f     ", rTTJets_HF ));
  if( rTTJets_LF>0 ) line += string(Form("%.4f     ", rTTJets_LF));
  if( rTTV>0 )       line += string(Form("%.4f     ", rTTV));
  if( rSingleT>0 )   line += string(Form("%.4f     ", rSingleT));
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;


  // SYSTEMATICS
  line = "lumi                  lnN  ";
  if( rTTH125>0 )    line += "1.022      ";
  if( rTTJets_HF>0 ) line += "1.022      ";
  if( rTTJets_LF>0 ) line += "1.022      ";
  if( rTTV>0 )       line += "1.022      ";
  if( rSingleT>0 )   line += "1.022      ";
  out<<line;
  out<<endl;
  line = "csv                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HF>0 ) line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;
  line = "JEC                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HF>0 ) line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    if(h_TTH125->GetBinContent(bin)>0){
      line = string(Form("TTH125%sbin%d        shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += "1.0        ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HF->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HF->GetBinContent(bin)>0){
      line = string(Form("TTJetsHF%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += "1.0        ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_LF->GetBinContent(bin)>0){
      line = string(Form("TTJetsLF%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += "1.0        ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    if(h_TTV->GetBinContent(bin)>0){
      line = string(Form("TTV%sbin%d           shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += "1.0        ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  } 
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    if(h_SingleT->GetBinContent(bin)>0){
      line = string(Form("SingleT%sbin%d       shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += "1.0        ";
      out<<line;
      out<<endl;
    }
  } 

  line = "Norm_TTbb             lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += "1.50       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTV              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += "1.20       ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_SingleT          lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += "1.20       ";
  out<<line;
  out<<endl;

  line = "QCDscale_TTH          lnN    ";
  if( rTTH125>0 )    line += "1.12       ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsHF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += "1.35       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsLF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += "1.35       ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "pdf_gg                lnN    ";
  if( rTTH125>0 )    line += "1.03       ";
  if( rTTJets_HF>0 ) line += "1.03       ";
  if( rTTJets_LF>0 ) line += "1.03       ";
  if( rTTV>0 )       line += "1.03       ";
  if( rSingleT>0 )   line += "1.03       ";
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

}



void produceAll( float LumiScale = 19.5/12.1){

  produce("SL", "type==0",  "cat1",                    0.02, 0.50, LumiScale ,   1, 8);
  produce("SL", "type==1",  "cat2",                    0.02, 0.7,  LumiScale ,   1, 8);
  produce("SL", "type==2 && flag_type2>0",  "cat3",    0.02, 0.8,  LumiScale ,   1, 8);
  produce("SL", "type==2 && flag_type2<0",  "cat4",    0.02, 0.7,  LumiScale ,   1, 8);
  produce("SL", "type==3",  "cat5",                    0.02, 0.5,  LumiScale ,   1, 8);  
  produce("DL", "type==6",  "cat6",                    0.02, 0.,   LumiScale*2 , 0, 6);
  produce("DL", "type==7",  "cat7",                    0.02, 0.,   LumiScale*2 , 0, 6);

}


// cat1  
// 0.1  => 8.28
// 0.25 => 7.93
// 0.5  => 7.95
// 0.75 => 8.03
// 1.0  => 8.09

// cat2
// 0.10 => 12.2
// 0.25 => 11.3
// 0.50 => 11.3
// 1.50 => 12.0

// cat4
// 0.25 => 8.8
// 0.50 => 8.6
// 0.75 => 9.2

// cat5
// 0.10 => 10.8
// 0.25 => 11.3
// 0.5  => 10.8
// 0.75 => 11.3

// cat6
// 0.25 => 9.3
// 0.50 => 8.8
// 0.75 => 8.8
// 1.50 => 8.8


void makePlot(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);
  
  TLegend* leg = new TLegend(0.665635,0.73951,0.866873,0.88986,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);


  float X[]        = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 };
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expXs[]  = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float expYs[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  vector<string> categories;
  categories.push_back("SL_cat1");
  categories.push_back("SL_cat2");
  categories.push_back("SL_cat3");
  categories.push_back("SL_cat4");
  categories.push_back("SL_cat5");
  categories.push_back("DL_cat6");
  categories.push_back("DL_cat7");
  categories.push_back("SL");
  categories.push_back("DL");
  categories.push_back("COMB");

  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = TFile::Open(("higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);

    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
      if(k==0) expY2sL[b] = r;
      if(k==1) expY1sL[b] = r;
      if(k==2) expY[b]    = r;
      if(k==3) expY1sH[b] = r;
      if(k==4) expY2sH[b] = r;
    }


    cout << categories[b] << ": r<" <<  expY[b] << " @ 95% CL" << endl;

    expY1sH[b] = TMath::Abs(expY1sH[b]-expY[b]);
    expY1sL[b] = TMath::Abs(expY1sL[b]-expY[b]);
    expY2sH[b] = TMath::Abs(expY2sH[b]-expY[b]);
    expY2sL[b] = TMath::Abs(expY2sL[b]-expY[b]);


  }


  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");

  TGraphAsymmErrors* expected = new TGraphAsymmErrors(10, X, expY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(10, X, expY, expXs, expXs,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma = new TGraphAsymmErrors(10, X, expY, expXs, expXs,  expY2sL, expY2sH);


  oneSigma->SetMarkerColor(kBlack);
  oneSigma->SetMarkerStyle(kOpenCircle);
  oneSigma->SetFillColor(kGreen);
  oneSigma->SetFillStyle(1001);

  twoSigma->SetMarkerColor(kBlack);
  twoSigma->SetMarkerStyle(kOpenCircle);
  twoSigma->SetFillColor(kYellow);
  twoSigma->SetFillStyle(1001);

  expected->SetMarkerColor(kBlack);
  expected->SetMarkerStyle(kOpenCircle);
  expected->SetMarkerSize(1.5);
  expected->SetLineColor(kBlack);
  expected->SetLineWidth(2);
 
  mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);

  mg->Draw("a2");
  //mg->Draw("p");
  
  expected->Draw("pSAME");
  //twoSigma->Draw("a2");
  //oneSigma->Draw("a2SAME");
  //expected->Draw("pSAME");


  TF1 *line = new TF1("line","1",0,10);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  line->Draw("SAME");

  TF1 *lineML = new TF1("lineML","2.4",0,10);
  lineML->SetLineColor(kBlue);
  lineML->SetLineStyle(kDashed);
  lineML->SetLineWidth(3);

  lineML->Draw("SAME");

  TF1 *lineTTH = new TF1("lineTTH","4.1",0,10);
  lineTTH->SetLineColor(kMagenta);
  lineTTH->SetLineStyle(kDashed);
  lineTTH->SetLineWidth(3);

  lineTTH->Draw("SAME");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->Set(10,0,10);

  //mg->GetXaxis()->SetRange(-1,11);

  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, categories[b].c_str() );
  }

  mg->GetYaxis()->SetTitleOffset(0.97);
  mg->SetMinimum(0.8);
  mg->SetMaximum( 40);
  mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");

  leg->AddEntry( lineTTH, "HIG-13-019", "L");
  leg->AddEntry( lineML,  "HIG-13-020", "L");
  leg->Draw();

  TPaveText *pt = new TPaveText(0.106811,0.155594,0.407121,0.286713,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(11);
  pt->AddText(Form("Comb: #mu < %.1f at 95%% CL", expY[categories.size()-1]))->SetTextColor(kRed);
  pt->AddText(Form("(SL: #mu < %.1f)", expY[categories.size()-3]))->SetTextColor(kBlack);
  pt->AddText(Form("(DL: #mu < %.1f)", expY[categories.size()-2]))->SetTextColor(kBlack);

  pt->Draw();

  //c1->Modified();
  //c1->Draw();
  c1->SaveAs("Limits.png");

}
