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
  hout_Down->SetBinContent(bin, bin_down);

  hout_Up->Reset();
  hout_Up->Add(hin,1.0);
  hout_Up->SetBinContent  (bin, bin_up);

  return;

}


void produce( TString fname = "SL_VType2", string cut = "type==0", TString category = "cat0", float fact = 0.5, float lumiScale = 20./12.1 ){


  TFile* fout = new TFile(fname+".root","UPDATE");
  TDirectory* dir = fout->mkdir( fname+"_"+category );

  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units", 6, 0,1 );
 
  TString var(Form("p_125_all_s/(p_125_all_s+p_125_all_b/100/%f)", fact));

  string basecut = cut;
  TString ttjets = (string(fname.Data())).find("SL")!=string::npos ? "TTJetsSemiLept" : "TTJetsFullLept";

  ////////////////////////////// TTH
  
  // NOMINAL
  TFile* f_TTH125 = TFile::Open("MEAnalysis_"+fname+"_nominal_TTH125.root");
  TTree* t_TTH125 = (TTree*)f_TTH125->Get("tree");
  TH1F* h_TTH125 = (TH1F*)h->Clone("h_TTH125");
  h_TTH125->Reset();
  h_TTH125->Sumw2();
  t_TTH125->Draw(var+">>h_TTH125",   "weight"*TCut(cut.c_str()) );
  h_TTH125->Scale(lumiScale);
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
  TFile* f_TTH125_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp_TTH125.root");
  TTree* t_TTH125_JECUp = (TTree*)f_TTH125_JECUp->Get("tree");
  TH1F* h_TTH125_JECUp = (TH1F*)h->Clone("h_TTH125_JECUp");
  h_TTH125_JECUp->Reset();
  h_TTH125_JECUp->Sumw2();
  t_TTH125_JECUp->Draw(var+">>h_TTH125_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECUp->Scale(lumiScale);
  dir->cd();
  h_TTH125_JECUp->Write("TTH125_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTH125_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown_TTH125.root");
  TTree* t_TTH125_JECDown = (TTree*)f_TTH125_JECDown->Get("tree");
  TH1F* h_TTH125_JECDown = (TH1F*)h->Clone("h_TTH125_JECDown");
  h_TTH125_JECDown->Reset();
  h_TTH125_JECDown->Sumw2();
  t_TTH125_JECDown->Draw(var+">>h_TTH125_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECDown->Scale(lumiScale);
  dir->cd();
  h_TTH125_JECDown->Write("TTH125_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTH125_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp_TTH125.root");
  TTree* t_TTH125_csvUp = (TTree*)f_TTH125_csvUp->Get("tree");
  TH1F* h_TTH125_csvUp = (TH1F*)h->Clone("h_TTH125_csvUp");
  h_TTH125_csvUp->Reset();
  h_TTH125_csvUp->Sumw2();
  t_TTH125_csvUp->Draw(var+">>h_TTH125_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvUp->Scale(lumiScale);
  dir->cd();
  h_TTH125_csvUp->Write("TTH125_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTH125_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown_TTH125.root");
  TTree* t_TTH125_csvDown = (TTree*)f_TTH125_csvDown->Get("tree");
  TH1F* h_TTH125_csvDown = (TH1F*)h->Clone("h_TTH125_csvDown");
  h_TTH125_csvDown->Reset();
  h_TTH125_csvDown->Sumw2();
  t_TTH125_csvDown->Draw(var+">>h_TTH125_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvDown->Scale(lumiScale);
  dir->cd();
  h_TTH125_csvDown->Write("TTH125_csvDown", TObject::kOverwrite);

  ////////////////////////////// TTJets: HF

  cut = basecut+"&&nSimBs>2";

  // NOMINAL
  TFile* f_TTJetsSemiLept_HF = TFile::Open("MEAnalysis_"+fname+"_nominal_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF = (TTree*)f_TTJetsSemiLept_HF->Get("tree");
  TH1F* h_TTJetsSemiLept_HF = (TH1F*)h->Clone("h_TTJetsSemiLept_HF");
  h_TTJetsSemiLept_HF->Reset();
  h_TTJetsSemiLept_HF->Sumw2();
  t_TTJetsSemiLept_HF->Draw(var+">>h_TTJetsSemiLept_HF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF->Scale(lumiScale);
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
  TFile* f_TTJetsSemiLept_HF_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECUp = (TTree*)f_TTJetsSemiLept_HF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECUp");
  h_TTJetsSemiLept_HF_JECUp->Reset();
  h_TTJetsSemiLept_HF_JECUp->Sumw2();
  t_TTJetsSemiLept_HF_JECUp->Draw(var+">>h_TTJetsSemiLept_HF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECUp->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_HF_JECUp->Write("TTJetsHF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HF_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECDown = (TTree*)f_TTJetsSemiLept_HF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECDown");
  h_TTJetsSemiLept_HF_JECDown->Reset();
  h_TTJetsSemiLept_HF_JECDown->Sumw2();
  t_TTJetsSemiLept_HF_JECDown->Draw(var+">>h_TTJetsSemiLept_HF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECDown->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_HF_JECDown->Write("TTJetsHF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HF_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvUp = (TTree*)f_TTJetsSemiLept_HF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvUp");
  h_TTJetsSemiLept_HF_csvUp->Reset();
  h_TTJetsSemiLept_HF_csvUp->Sumw2();
  t_TTJetsSemiLept_HF_csvUp->Draw(var+">>h_TTJetsSemiLept_HF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvUp->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_HF_csvUp->Write("TTJetsHF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HF_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvDown = (TTree*)f_TTJetsSemiLept_HF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvDown");
  h_TTJetsSemiLept_HF_csvDown->Reset();
  h_TTJetsSemiLept_HF_csvDown->Sumw2();
  t_TTJetsSemiLept_HF_csvDown->Draw(var+">>h_TTJetsSemiLept_HF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvDown->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_HF_csvDown->Write("TTJetsHF_csvDown", TObject::kOverwrite);


  ////////////////////////////// TTJets: LF

  cut = basecut+"&&nSimBs==2";

  // NOMINAL
  TFile* f_TTJetsSemiLept_LF = TFile::Open("MEAnalysis_"+fname+"_nominal_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF = (TTree*)f_TTJetsSemiLept_LF->Get("tree");
  TH1F* h_TTJetsSemiLept_LF = (TH1F*)h->Clone("h_TTJetsSemiLept_LF");
  h_TTJetsSemiLept_LF->Reset();
  h_TTJetsSemiLept_LF->Sumw2();
  t_TTJetsSemiLept_LF->Draw(var+">>h_TTJetsSemiLept_LF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF->Scale(lumiScale);
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
  TFile* f_TTJetsSemiLept_LF_JECUp = TFile::Open("MEAnalysis_"+fname+"_JECUp_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECUp = (TTree*)f_TTJetsSemiLept_LF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECUp");
  h_TTJetsSemiLept_LF_JECUp->Reset();
  h_TTJetsSemiLept_LF_JECUp->Sumw2();
  t_TTJetsSemiLept_LF_JECUp->Draw(var+">>h_TTJetsSemiLept_LF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECUp->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_LF_JECUp->Write("TTJetsLF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_LF_JECDown = TFile::Open("MEAnalysis_"+fname+"_JECDown_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECDown = (TTree*)f_TTJetsSemiLept_LF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECDown");
  h_TTJetsSemiLept_LF_JECDown->Reset();
  h_TTJetsSemiLept_LF_JECDown->Sumw2();
  t_TTJetsSemiLept_LF_JECDown->Draw(var+">>h_TTJetsSemiLept_LF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECDown->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_LF_JECDown->Write("TTJetsLF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_LF_csvUp = TFile::Open("MEAnalysis_"+fname+"_csvUp_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvUp = (TTree*)f_TTJetsSemiLept_LF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvUp");
  h_TTJetsSemiLept_LF_csvUp->Reset();
  h_TTJetsSemiLept_LF_csvUp->Sumw2();
  t_TTJetsSemiLept_LF_csvUp->Draw(var+">>h_TTJetsSemiLept_LF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvUp->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_LF_csvUp->Write("TTJetsLF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_LF_csvDown = TFile::Open("MEAnalysis_"+fname+"_csvDown_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvDown = (TTree*)f_TTJetsSemiLept_LF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvDown");
  h_TTJetsSemiLept_LF_csvDown->Reset();
  h_TTJetsSemiLept_LF_csvDown->Sumw2();
  t_TTJetsSemiLept_LF_csvDown->Draw(var+">>h_TTJetsSemiLept_LF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvDown->Scale(lumiScale);
  dir->cd();
  h_TTJetsSemiLept_LF_csvDown->Write("TTJetsLF_csvDown", TObject::kOverwrite);



  cout << "TTH125: "      << h_TTH125->Integral()         << endl;
  cout << "TTJets (HF): " << h_TTJetsSemiLept_HF->Integral() << endl;
  cout << "TTJets (LF): " << h_TTJetsSemiLept_LF->Integral() << endl;

  float rTTH125    = h_TTH125->Integral() ;
  float rTTJets_HF = h_TTJetsSemiLept_HF->Integral();
  float rTTJets_LF = h_TTJetsSemiLept_LF->Integral();

  float observation = int(rTTH125 + rTTJets_HF + rTTJets_LF);

  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_TTH125,1.0);
  h_data->Add(h_TTJetsSemiLept_HF,1.0);
  h_data->Add(h_TTJetsSemiLept_LF,1.0);
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
  string line6(Form("bin                %s  %s         %s          %s",       longspace.c_str(), (fname+"_"+category).Data(), (fname+"_"+category).Data(), (fname+"_"+category).Data() )); out << line6 << endl;
  string line7(Form("process            %s  TTH125                 TTJetsHF                TTJetsLF", longspace.c_str())); out << line7 << endl;
  string line8(Form("process            %s  0                      1                       2",        longspace.c_str())); out << line8 << endl;
  string line9(Form("rate               %s  %.4f                 %.4f                   %.4f",      longspace.c_str(), rTTH125, rTTJets_HF, rTTJets_LF)); out << line9 << endl;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

  out<< string("lumi                    lnN         1.022          1.022            1.022") <<endl;
  out<< string("csv                   shape         1.00           1.00             1.00") <<endl;
  out<< string("JEC                   shape         1.00           1.00             1.00") <<endl;

  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    out<< string(Form("TTH125%sbin%d            shape         1.0            -                -  ",  category.Data(), bin)) <<endl;
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HF->GetNbinsX(); bin++ ){
    out<< string(Form("TTJetsHF%sbin%d          shape       -             1.0               -  ", category.Data(), bin)) <<endl;
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    out<< string(Form("TTJetsLF%sbin%d          shape       -              -               1.0 ", category.Data(), bin)) <<endl;
  }  

  out<< string     ("QCDscale_TTbb           lnN        -            1.50              -  ") <<endl;
  out<< string     ("QCDscale_TTH            lnN      1.12            -                -  ") <<endl;
  out<< string(Form("QCDscale%s_TTJetsHF       lnN        -            1.35              -  ",  null.c_str() )) <<endl;
  out<< string(Form("QCDscale%s_TTJetsLF       lnN        -             -               1.35",  null.c_str() )) <<endl;
  out<< string     ("pdf_gg                  lnN      1.03           1.03             1.03") <<endl;

  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

}



void produceAll( float LumiScale = 20./12.1){

  produce("SL", "type==0",  "cat1",  0.5,  LumiScale );
  produce("SL", "type==1",  "cat2",  0.5,  LumiScale );
  produce("SL", "type==2 && flag_type2>0",  "cat3",  0.5,  LumiScale );
  produce("SL", "type==2 && flag_type2<0",  "cat4",  0.5,  LumiScale );
  produce("SL", "type==3",  "cat5",  0.5,  LumiScale );

  produce("DL", "type==6",  "cat6",  0.5,  LumiScale*2 );
  produce("DL", "type==7",  "cat7",  0.5,  LumiScale*2 );

}



void makePlot(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);
  
  TLegend* leg = new TLegend(0.65,0.70,0.85,0.85,NULL,"brNDC");
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
  lineML->SetLineWidth(1);

  lineML->Draw("SAME");

  TF1 *lineTTH = new TF1("lineTTH","4.1",0,10);
  lineTTH->SetLineColor(kMagenta);
  lineTTH->SetLineStyle(kDashed);
  lineTTH->SetLineWidth(1);

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
  mg->SetMaximum( 100);
  mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");

  leg->AddEntry( lineML,  "HIG-13-020", "L");
  leg->AddEntry( lineTTH, "HIG-13-019", "L");
  leg->Draw();

  c1->SaveAs("Limits.png");

}
