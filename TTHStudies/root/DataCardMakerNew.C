#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTBase.h"
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
#include "TArrayF.h"
#include "TLine.h"

typedef TMatrixT<double> TMatrixD;


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



pair<double,double> getMaxValue( TH1F* hMassProb){
 
  float est  = -99;
  float prob = -99;
  float a = -99;
  float b = -99;
  float c = -99;

  if( hMassProb->GetEntries()==0 ){
    est  =  0.;
    prob =  0.;
    return make_pair(est, prob);
  }

  
  int maxBin = hMassProb->GetMaximumBin();
  int bD = -999;
  int bC = -999;
  int bU = -999;

  int isRightMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin+k)>0 ) isRightMostBin=0;
  }
  int isLeftMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin-k)>0 ) isLeftMostBin=0;
  }

  if( (maxBin < hMassProb->GetNbinsX() && maxBin > 1) && !isRightMostBin && !isLeftMostBin ){
    //cout << "Center" << endl;
    bC = maxBin;
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(maxBin-k)>0 ){
	bD = maxBin-k;	
	////cout << "Found bin " << bD << " with " <<  hMassProb->GetBinContent(bD) << endl;
      }
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(maxBin+k)>0 ){
	bU = maxBin+k;	
	////cout << "Found bin " << bU << " with " <<  hMassProb->GetBinContent(bU) << endl;	
      }
    }   
  }
  else if( (maxBin == hMassProb->GetNbinsX()) || isRightMostBin){
    //cout << "Right" << endl;
    bU = maxBin; 
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(maxBin-k)>0 ) bC = maxBin-k;	
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(bC-k)>0 ) bD = bC-k;	
    }
  }
  else if( (maxBin == 1) || isLeftMostBin){
    //cout << "Left" << endl;
    bD = maxBin; 
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(maxBin+k)>0 ) bC = maxBin+k;	
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(bC+k)>0 ) bU = bU+k;	
    }
  }
  else{
    //cout << "Unknown" << endl;
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    //cout << "M=" << est << endl;
    return make_pair(est, prob);
  }

  if( bD==-999 || bU==-999 || bC==-999){
    //cout << "Problems" << endl;   
    //cout << bD << "," << bC << "," << bU << endl;
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    //cout << "M=" << est << endl;
    return make_pair(est, prob);
  }

  double xD =  hMassProb->GetBinCenter (bD);
  double xC =  hMassProb->GetBinCenter (bC);
  double xU =  hMassProb->GetBinCenter (bU);
  double yD =  hMassProb->GetBinContent(bD);
  double yC =  hMassProb->GetBinContent(bC);
  double yU =  hMassProb->GetBinContent(bU);
  
  TMatrixD A(3,3);
  const double elements[9] = 
    { 1, xD,  xD*xD,
      1, xC,  xC*xC,
      1, xU,  xU*xU 
    };
  A.SetMatrixArray(elements, "C");
  TMatrixD AInv(3,3);
  double det;
  AInv = A.Invert(&det);
  
  TMatrixD Y(3,1);
  const double yPos[3] = 
    { yD, yC, yU
    };
  Y.SetMatrixArray(yPos, "C");
  
  TMatrixD C(3,1);
  const double dummy[3] = 
    { 1., 1., 1.
    };
  C.SetMatrixArray(dummy,"C");
  C.Mult(AInv,Y);
  
  a = C(2,0);
  b = C(1,0);
  c = C(0,0);
  
  est  = -b/2/a ;
  prob = a*est*est + b*est + c;
  
  //cout << "xD" << xD << " => " <<  yD << endl;
  //cout << "xC" << xC << " => " <<  yC << endl;
  //cout << "xU" << xU << " => " <<  yU << endl;
  //cout << "M=" << est << endl;
  
  return make_pair(est, prob);
  
}



void produce( TString fname = "SL", 
	      string cut = "type==0", 
	      TString category = "cat0", 
	      float fact1 = 0.02, float fact2 = 0., float factbb = 0.5,
	      float lumiScale = 20./12.1,
	      int nBins =6 , int splitFirstBin=0){


  string version = (string(fname.Data())).find("SL")!=string::npos ? "_v1" : "_v1";
  //string version = "";

  cout << "Doing version " << version << " and category " << category << endl;

  string basecut = cut;
  TString ttjets = (string(fname.Data())).find("SL")!=string::npos ? "_TTJets" : "_TTJets";


  float scaleTTH      = 1.0;//4.7;
  float scaleTTJetsLF = 1.0;//3.9;
  float scaleTTJetsHF = 1.0;//3.9;

  cout << "*****************************************************" << endl;
  cout << "Attention!!! The following rescaling will be applied:" << endl;
  cout << "TTH:       " << scaleTTH << endl;
  cout << "TTJets LF: " << scaleTTJetsLF << endl;
  cout << "TTJets HF: " << scaleTTJetsHF << endl;
  cout << "*****************************************************" << endl;



  TFile* f_B = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+ttjets+".root");
  TTree* tB = (TTree*)f_B->Get("tree");
  TFile* f_S = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_TTH125.root");
  TTree* tS = (TTree*)f_S->Get("tree");

  TH1F* hInt = new TH1F("hInt","",200,-100,100);
  

  tS->Draw("log(p_125_all_s_ttbb)>>hInt", TCut((basecut+" && nSimBs>=2 && p_125_all_s_ttbb>0").c_str())); 
  float S_bb     = TMath::Exp(hInt->GetMean());
  float Int_S_bb = hInt->Integral();
  cout << "S_bb = " << S_bb <<endl;
  hInt->Reset();

  tB->Draw("log(p_125_all_b_ttbb)>>hInt", TCut((basecut+" && nMatchSimBs>=2 && nSimBs>2 && p_125_all_b_ttbb>0").c_str())); 
  float B_bb     = TMath::Exp(hInt->GetMean());
  float Int_B_bb = hInt->Integral();
  cout << "B_bb = " << B_bb <<endl;
  hInt->Reset();

  tB->Draw("log(p_125_all_b_ttjj)>>hInt", TCut((basecut+" && nMatchSimBs<2  && nSimBs==2 && p_125_all_b_ttjj>0").c_str())); 
  float B_jj     = TMath::Exp(hInt->GetMean());
  float Int_B_jj = hInt->Integral();
  cout << "B_jj = " << B_jj <<endl;
  hInt->Reset();

  cout << "Ratio S_bb/B_bb=" << S_bb/B_bb << endl;
  cout << "Ratio S_bb/B_jj=" << S_bb/B_jj << endl;
  cout << "Ratio bb/(bb+jj)=" << Int_B_bb << "/(" << Int_B_bb << "+" << Int_B_jj << ")=" << Int_B_bb/(Int_B_bb+Int_B_jj) << endl;


  vector<float> param;
  param.clear();
  if(splitFirstBin){
    param.push_back( B_bb/B_jj*factbb );
    param.push_back( 0.50 );
  }

  TFile* fout = new TFile(fname+"_New.root","UPDATE");
  TDirectory* dir =  fout->GetDirectory( fname+"_"+category); 
  if( !dir) dir = fout->mkdir( fname+"_"+category ) ;

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


  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray());
 
  TString var("");
  var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*(%f*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, (1-fact2)*S_bb/B_bb*(Int_B_bb/(Int_B_bb+Int_B_jj)), (1+fact2)*S_bb/B_jj*(Int_B_jj/(Int_B_bb+Int_B_jj)) ));

  cout << "Variable=" << string(var.Data()) << endl;

  //var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));


  ////////////////////////////// TTH
  
  // NOMINAL
  TFile* f_TTH125 = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_TTH125.root");
  TTree* t_TTH125 = (TTree*)f_TTH125->Get("tree");
  TH1F* h_TTH125 = (TH1F*)h->Clone("h_TTH125");
  h_TTH125->Reset();
  h_TTH125->Sumw2();
  draw( param, t_TTH125, var, h_TTH125, TCut(cut.c_str()));
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
  TFile* f_TTH125_JECUp = TFile::Open("MEAnalysisNew_"+fname+"_JECUp"+version+"_TTH125.root");
  TTree* t_TTH125_JECUp = (TTree*)f_TTH125_JECUp->Get("tree");
  TH1F* h_TTH125_JECUp = (TH1F*)h->Clone("h_TTH125_JECUp");
  h_TTH125_JECUp->Reset();
  h_TTH125_JECUp->Sumw2();
  draw( param, t_TTH125_JECUp, var, h_TTH125_JECUp, TCut(cut.c_str()));
  h_TTH125_JECUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECUp->Write("TTH125_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTH125_JECDown = TFile::Open("MEAnalysisNew_"+fname+"_JECDown"+version+"_TTH125.root");
  TTree* t_TTH125_JECDown = (TTree*)f_TTH125_JECDown->Get("tree");
  TH1F* h_TTH125_JECDown = (TH1F*)h->Clone("h_TTH125_JECDown");
  h_TTH125_JECDown->Reset();
  h_TTH125_JECDown->Sumw2();
  draw( param, t_TTH125_JECDown, var, h_TTH125_JECDown, TCut(cut.c_str()));
  h_TTH125_JECDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECDown->Write("TTH125_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTH125_csvUp = TFile::Open("MEAnalysisNew_"+fname+"_csvUp"+version+"_TTH125.root");
  TTree* t_TTH125_csvUp = (TTree*)f_TTH125_csvUp->Get("tree");
  TH1F* h_TTH125_csvUp = (TH1F*)h->Clone("h_TTH125_csvUp");
  h_TTH125_csvUp->Reset();
  h_TTH125_csvUp->Sumw2();
  draw( param, t_TTH125_csvUp, var, h_TTH125_csvUp, TCut(cut.c_str()));
  h_TTH125_csvUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvUp->Write("TTH125_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTH125_csvDown = TFile::Open("MEAnalysisNew_"+fname+"_csvDown"+version+"_TTH125.root");
  TTree* t_TTH125_csvDown = (TTree*)f_TTH125_csvDown->Get("tree");
  TH1F* h_TTH125_csvDown = (TH1F*)h->Clone("h_TTH125_csvDown");
  h_TTH125_csvDown->Reset();
  h_TTH125_csvDown->Sumw2();
  draw( param, t_TTH125_csvDown, var, h_TTH125_csvDown, TCut(cut.c_str()));
  h_TTH125_csvDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvDown->Write("TTH125_csvDown", TObject::kOverwrite);




  ////////////////////////////// TTJets: HFbb


  cut = basecut+"&&nSimBs>2 && nMatchSimBs>=2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HFbb = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb = (TTree*)f_TTJetsSemiLept_HFbb->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb");
  h_TTJetsSemiLept_HFbb->Reset();
  h_TTJetsSemiLept_HFbb->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb, var, h_TTJetsSemiLept_HFbb, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFbb->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb->Write("TTJetsHFbb", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFbb->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HFbb_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFbb_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HFbb_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFbb_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HFbb, h_TTJetsSemiLept_HFbb_b_up, h_TTJetsSemiLept_HFbb_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HFbb_b_up  ->Write(Form("TTJetsHFbb_TTJetsHFbb%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HFbb_b_down->Write(Form("TTJetsHFbb_TTJetsHFbb%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HFbb_JECUp = TFile::Open("MEAnalysisNew_"+fname+"_JECUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_JECUp = (TTree*)f_TTJetsSemiLept_HFbb_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_JECUp");
  h_TTJetsSemiLept_HFbb_JECUp->Reset();
  h_TTJetsSemiLept_HFbb_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_JECUp , var, h_TTJetsSemiLept_HFbb_JECUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFbb_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_JECUp->Write("TTJetsHFbb_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HFbb_JECDown = TFile::Open("MEAnalysisNew_"+fname+"_JECDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_JECDown = (TTree*)f_TTJetsSemiLept_HFbb_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_JECDown");
  h_TTJetsSemiLept_HFbb_JECDown->Reset();
  h_TTJetsSemiLept_HFbb_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_JECDown, var, h_TTJetsSemiLept_HFbb_JECDown, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFbb_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_JECDown->Write("TTJetsHFbb_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HFbb_csvUp = TFile::Open("MEAnalysisNew_"+fname+"_csvUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_csvUp = (TTree*)f_TTJetsSemiLept_HFbb_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_csvUp");
  h_TTJetsSemiLept_HFbb_csvUp->Reset();
  h_TTJetsSemiLept_HFbb_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_csvUp, var, h_TTJetsSemiLept_HFbb_csvUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFbb_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_csvUp->Write("TTJetsHFbb_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HFbb_csvDown = TFile::Open("MEAnalysisNew_"+fname+"_csvDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_csvDown = (TTree*)f_TTJetsSemiLept_HFbb_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_csvDown");
  h_TTJetsSemiLept_HFbb_csvDown->Reset();
  h_TTJetsSemiLept_HFbb_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_csvDown, var, h_TTJetsSemiLept_HFbb_csvDown, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFbb_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_csvDown->Write("TTJetsHFbb_csvDown", TObject::kOverwrite);


  ////////////////////////////// TTJets: HFb


  cut = basecut+"&&nSimBs>2 && nMatchSimBs<2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HFb = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb = (TTree*)f_TTJetsSemiLept_HFb->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb");
  h_TTJetsSemiLept_HFb->Reset();
  h_TTJetsSemiLept_HFb->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb, var, h_TTJetsSemiLept_HFb, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFb->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb->Write("TTJetsHFb", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFb->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HFb_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFb_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HFb_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFb_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HFb, h_TTJetsSemiLept_HFb_b_up, h_TTJetsSemiLept_HFb_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HFb_b_up  ->Write(Form("TTJetsHFb_TTJetsHFb%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HFb_b_down->Write(Form("TTJetsHFb_TTJetsHFb%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HFb_JECUp = TFile::Open("MEAnalysisNew_"+fname+"_JECUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_JECUp = (TTree*)f_TTJetsSemiLept_HFb_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_JECUp");
  h_TTJetsSemiLept_HFb_JECUp->Reset();
  h_TTJetsSemiLept_HFb_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_JECUp , var, h_TTJetsSemiLept_HFb_JECUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFb_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_JECUp->Write("TTJetsHFb_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HFb_JECDown = TFile::Open("MEAnalysisNew_"+fname+"_JECDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_JECDown = (TTree*)f_TTJetsSemiLept_HFb_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_JECDown");
  h_TTJetsSemiLept_HFb_JECDown->Reset();
  h_TTJetsSemiLept_HFb_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_JECDown, var, h_TTJetsSemiLept_HFb_JECDown, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFb_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_JECDown->Write("TTJetsHFb_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HFb_csvUp = TFile::Open("MEAnalysisNew_"+fname+"_csvUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_csvUp = (TTree*)f_TTJetsSemiLept_HFb_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_csvUp");
  h_TTJetsSemiLept_HFb_csvUp->Reset();
  h_TTJetsSemiLept_HFb_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_csvUp, var, h_TTJetsSemiLept_HFb_csvUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFb_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_csvUp->Write("TTJetsHFb_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HFb_csvDown = TFile::Open("MEAnalysisNew_"+fname+"_csvDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_csvDown = (TTree*)f_TTJetsSemiLept_HFb_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_csvDown");
  h_TTJetsSemiLept_HFb_csvDown->Reset();
  h_TTJetsSemiLept_HFb_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_csvDown, var, h_TTJetsSemiLept_HFb_csvDown, TCut(cut.c_str()));
  h_TTJetsSemiLept_HFb_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_csvDown->Write("TTJetsHFb_csvDown", TObject::kOverwrite);



  ////////////////////////////// TTJets: LF

  cut = basecut+"&&nSimBs==2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_LF = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF = (TTree*)f_TTJetsSemiLept_LF->Get("tree");
  TH1F* h_TTJetsSemiLept_LF = (TH1F*)h->Clone("h_TTJetsSemiLept_LF");
  h_TTJetsSemiLept_LF->Reset();
  h_TTJetsSemiLept_LF->Sumw2();
  draw( param, t_TTJetsSemiLept_LF, var, h_TTJetsSemiLept_LF, TCut(cut.c_str()));
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
  TFile* f_TTJetsSemiLept_LF_JECUp = TFile::Open("MEAnalysisNew_"+fname+"_JECUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECUp = (TTree*)f_TTJetsSemiLept_LF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECUp");
  h_TTJetsSemiLept_LF_JECUp->Reset();
  h_TTJetsSemiLept_LF_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECUp, var, h_TTJetsSemiLept_LF_JECUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_LF_JECUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECUp->Write("TTJetsLF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_LF_JECDown = TFile::Open("MEAnalysisNew_"+fname+"_JECDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECDown = (TTree*)f_TTJetsSemiLept_LF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECDown");
  h_TTJetsSemiLept_LF_JECDown->Reset();
  h_TTJetsSemiLept_LF_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECDown , var, h_TTJetsSemiLept_LF_JECDown, TCut(cut.c_str()));
  h_TTJetsSemiLept_LF_JECDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECDown->Write("TTJetsLF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_LF_csvUp = TFile::Open("MEAnalysisNew_"+fname+"_csvUp"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvUp = (TTree*)f_TTJetsSemiLept_LF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvUp");
  h_TTJetsSemiLept_LF_csvUp->Reset();
  h_TTJetsSemiLept_LF_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvUp, var, h_TTJetsSemiLept_LF_csvUp, TCut(cut.c_str()));
  h_TTJetsSemiLept_LF_csvUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvUp->Write("TTJetsLF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_LF_csvDown = TFile::Open("MEAnalysisNew_"+fname+"_csvDown"+version+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvDown = (TTree*)f_TTJetsSemiLept_LF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvDown");
  h_TTJetsSemiLept_LF_csvDown->Reset();
  h_TTJetsSemiLept_LF_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvDown, var, h_TTJetsSemiLept_LF_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_csvDown->Draw(var+">>h_TTJetsSemiLept_LF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvDown->Write("TTJetsLF_csvDown", TObject::kOverwrite);



  ////////////////////////////// DIBOSON
  cut = basecut;
  
  // NOMINAL
  TFile* f_DiBoson = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_DiBoson.root");
  TTree* t_DiBoson = (TTree*)f_DiBoson->Get("tree");
  TH1F* h_DiBoson = (TH1F*)h->Clone("h_DiBoson");
  h_DiBoson->Reset();
  h_DiBoson->Sumw2();
  draw( param, t_DiBoson , var, h_DiBoson, TCut(cut.c_str()));
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
  TFile* f_SingleT = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_SingleT.root");
  TTree* t_SingleT = (TTree*)f_SingleT->Get("tree");
  TH1F* h_SingleT = (TH1F*)h->Clone("h_SingleT");
  h_SingleT->Reset();
  h_SingleT->Sumw2();
  draw( param, t_SingleT , var, h_SingleT, TCut(cut.c_str()));
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
  TFile* f_TTV = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_TTV.root");
  TTree* t_TTV = (TTree*)f_TTV->Get("tree");
  TH1F* h_TTV = (TH1F*)h->Clone("h_TTV");
  h_TTV->Reset();
  h_TTV->Sumw2();
  draw( param, t_TTV , var, h_TTV, TCut(cut.c_str()));
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
  TFile* f_EWK = TFile::Open("MEAnalysisNew_"+fname+"_nominal"+version+"_EWK.root");
  TTree* t_EWK = (TTree*)f_EWK->Get("tree");
  TH1F* h_EWK = (TH1F*)h->Clone("h_EWK");
  h_EWK->Reset();
  h_EWK->Sumw2();
  draw( param, t_EWK, var, h_EWK, TCut(cut.c_str()));
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
  cout << "TTJets (HFbb): "<< h_TTJetsSemiLept_HFbb->Integral() << endl;
  cout << "TTJets (HFb): " << h_TTJetsSemiLept_HFb->Integral() << endl;
  cout << "TTJets (LF): "  << h_TTJetsSemiLept_LF->Integral() << endl;
  cout << "SingleT:     "  << h_SingleT->Integral() << endl;
  cout << "TTV:         "  << h_TTV->Integral() << endl;
  cout << "EWK:         "  << h_EWK->Integral() << endl;
  cout << "DiBoson:     "  << h_DiBoson->Integral() << endl;

  float rTTH125    = h_TTH125->Integral() ;
  float rTTJets_HFbb = h_TTJetsSemiLept_HFbb->Integral();
  float rTTJets_HFb  = h_TTJetsSemiLept_HFb->Integral();
  float rTTJets_LF = h_TTJetsSemiLept_LF->Integral();
  float rSingleT   = h_SingleT->Integral() ;
  float rTTV       = h_TTV->Integral() ;
  float rEWK       = h_EWK->Integral() ;
  float rDiBoson   = h_DiBoson->Integral() ;

  float observation = int(rTTH125 + rTTJets_HFbb + rTTJets_HFb + rTTJets_LF  + rTTV + rSingleT /*+ rEWK + rDiBoson*/ );

  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_TTH125,1.0);
  h_data->Add(h_TTJetsSemiLept_HFbb,1.0);
  h_data->Add(h_TTJetsSemiLept_HFb,1.0);
  h_data->Add(h_TTJetsSemiLept_LF,1.0);
  h_data->Add(h_TTV);
  h_data->Add(h_SingleT);
  //h_data->Add(h_EWK);
  //h_data->Add(h_DiBoson);

  h_data->Scale(int(observation)/h_data->Integral());
  dir->cd();
  h_data->Write("data_obs", TObject::kOverwrite);

  fout->Close();


  ofstream out(Form("%s_New.txt",(fname+"_"+category).Data()));
  out.precision(8);
  string longspace  = "              ";
  string shortspace = "          ";


  //string null(category.Data());
  string null("");


  string line1("imax 1"); out << line1 << endl;
  string line2("jmax *"); out << line2 << endl;
  string line3("kmax *"); out << line3 << endl;
  string line4(Form("shapes *  *    %s_New.root  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC", fname.Data())); out << line4 << endl;
  out<<endl;
  string line5(Form("observation %.0f", observation )); out << line5 << endl;
  out<<endl;


  // PROC AND NORM
  string line("bin                         ");
  if( rTTH125>0 )    line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HFbb>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HFb>0 )  line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_LF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTV>0 )       line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rSingleT>0 )   line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";
  if( rTTH125>0 )    line += "TTH125      ";
  if( rTTJets_HFbb>0 ) line += "TTJetsHFbb    ";
  if( rTTJets_HFb>0 )  line += "TTJetsHFb    ";

  if( rTTJets_LF>0 ) line += "TTJetsLF    ";
  if( rTTV>0 )       line += "TTV      ";
  if( rSingleT>0 )   line += "SingleT     ";
  out<<line;
  out<<endl;
  line = "process                       ";
  if( rTTH125>0 )    line += "0          ";
  if( rTTJets_HFbb>0 ) line += "1          ";
  if( rTTJets_HFb>0 )  line += "2          ";
  if( rTTJets_LF>0 ) line += "3          ";
  if( rTTV>0 )       line += "4          ";
  if( rSingleT>0 )   line += "5          ";
  out<<line;
  out<<endl;
  line = "rate                        ";
  if( rTTH125>0 )    line += string(Form("%.4f     ", rTTH125 ));
  if( rTTJets_HFbb>0 ) line += string(Form("%.4f     ", rTTJets_HFbb ));
  if( rTTJets_HFb>0 )  line += string(Form("%.4f     ", rTTJets_HFb ));
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
  if( rTTJets_HFbb>0 ) line += "1.022      ";
  if( rTTJets_HFb>0 )  line += "1.022      ";
  if( rTTJets_LF>0 ) line += "1.022      ";
  if( rTTV>0 )       line += "1.022      ";
  if( rSingleT>0 )   line += "1.022      ";
  out<<line;
  out<<endl;
  line = "csv                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HFbb>0 ) line += "1.0        ";
  if( rTTJets_HFb>0 )  line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;
  line = "JEC                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HFbb>0 ) line += "1.0        ";
  if( rTTJets_HFb>0 )  line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    if(h_TTH125->GetBinContent(bin)>0){
      line = string(Form("TTH125%sbin%d        shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += "1.0        ";
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFbb->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HFbb->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFbb%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += "1.0        ";
      if( rTTJets_HFb>0 )  line += " -        ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFb->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HFb->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFb%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += " -        ";
      if( rTTJets_HFb>0 )  line += "1.0        ";
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
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
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
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
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
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += "1.0        ";
      out<<line;
      out<<endl;
    }
  } 

  line = "Norm_TTbb             lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += "1.50       ";
  if( rTTJets_HFb>0 )  line += " -       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTb              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -       ";
  if( rTTJets_HFb>0 )  line += "1.50       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTV              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += "1.20       ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_SingleT          lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += "1.20       ";
  out<<line;
  out<<endl;

  line = "QCDscale_TTH          lnN    ";
  if( rTTH125>0 )    line += "1.12       ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsHF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += "1.35       ";
  if( rTTJets_HFb>0 )  line += "1.35       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsLF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += "1.35       ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "pdf_gg                lnN    ";
  if( rTTH125>0 )    line += "1.03       ";
  if( rTTJets_HFbb>0 ) line += "1.03       ";
  if( rTTJets_HFb>0 )  line += "1.03       ";
  if( rTTJets_LF>0 ) line += "1.03       ";
  if( rTTV>0 )       line += "1.03       ";
  if( rSingleT>0 )   line += "1.03       ";
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

}



void produceAll( float LumiScale = 19.5/12.1){

  produce("SL", "type==0",                                  "cat1", 1.0, 0.0, 0.2 , LumiScale   , 6,  1);
  produce("SL", "type==1",                                  "cat2", 1.7, 0.0, 0.5 , LumiScale   , 6,  1);
  produce("SL", "type==2 && flag_type2>0",                  "cat3", 2.2, 0.0, 0.5 , LumiScale   , 6,  1);
  produce("SL", "type==2 && flag_type2<=0",                 "cat4", 2.0, 0.0, 0.5 , LumiScale   , 6,  1);
  produce("SL", "type==3 && flag_type3>0 && p_125_all_s>0", "cat5", 5.5, 0.0, 0.5 , LumiScale   , 7,  1);
  produce("DL", "type==6",                                  "cat6", 1.5, 0.0, 0.1 , LumiScale*2 , 5,  1);

}


void massPlot( TString fname = "SL", 
	       string cut = "type==0", 
	       TString category = "cat0",
	       int plotShapes = 0,
	       float lumiScale = 20./12.1,
	       int nBins =20, float xLow=50, float xHigh=250,
	       int useBTag=0,
	       int normXSec=1,
	       float fact1 = 5,
	       int match1 = 111111, int match2 = 111111, int match3 = 11, int match4 = 11
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


  TLegend* leg = new TLegend(0.52,0.44,0.77,0.90,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  leg->SetHeader(category);

  //TF1* xsec6J = new TF1("xsec6J",Form("TMath::Landau(x,%f*7.61581e+01 ,%f*1.89245e+01)", fact1, fact2), 20, 500);
  //TF1* xsec5J = new TF1("xsec5J",Form("TMath::Landau(x,%f*7.40196e+01 ,%f*1.80142e+01)", fact1, fact2), 20, 500); 
  //TF1* xsec = 
  //( (string(category.Data())).find("cat2")!=string::npos || 
  //  (string(category.Data())).find("cat3")!=string::npos ||
  //  (string(category.Data())).find("cat4")!=string::npos ) ?  xsec5J : xsec6J;
  
  TF1* xsec = new TF1("xsec",Form("x^(-%f)", fact1 ), 20, 500);

  string version = (string(fname.Data())).find("SL")!=string::npos ? "_v1" : "_v1";
  cout << "Doing version " << version << " and category " << category << endl;

  THStack* aStack = new THStack("aStack","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events");
  TH1F* hMass     = new TH1F   ("hMass", "Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events",nBins, xLow, xHigh);
  TH1F* hTTH      = new TH1F("hTTH","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events",nBins, xLow, xHigh);
  TH1F* hTTV      = new TH1F("hTTV","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events",nBins, xLow, xHigh);
  TH1F* hTT       = new TH1F("hTT", "Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events",nBins, xLow, xHigh);
  TH1F* hTmp      = new TH1F("hTmp", "", 500,0,500);

  TH1F* hErr      = 0;

  vector<string> samples;
  samples.push_back("TTJetsBB");
  samples.push_back("TTJetsBJ");
  samples.push_back("TTJetsJJ");
  samples.push_back("TTV");
  samples.push_back("SingleT");
  samples.push_back("DiBoson");
  samples.push_back("EWK");
  samples.push_back("TTH125");
  

  for( unsigned int s = 0 ; s < samples.size(); s++){

    string sample = samples[s];

    TCut sample_cut(cut.c_str());
    
    int color;
    if( sample.find("TTJetsBB")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs>2 && nMatchSimBs>=2");
      color = 16;
    }
    if( sample.find("TTJetsBJ")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs>2 && nMatchSimBs<2");
      color = 17;
    }
    if( sample.find("TTJetsJJ")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs==2");
      color = 18;
    }

    TFile* f = TFile::Open(("MEAnalysisNew_MHscan_"+string(fname.Data())+"_nominal"+version+"_"+sample+".root").c_str());
    //TFile* f = TFile::Open("../bin/root/MEAnalysisNew.root");
    if(f==0 || f->IsZombie()){
      cout << "Missing " << sample << " file" << endl;
      continue;
    }
    else{
      cout << "Doing sample " << sample << endl;
    }

    TH1F* hMass_s = (TH1F*)hMass->Clone(("hMass_"+sample).c_str());
    hMass_s->Reset();
    hMass_s->Sumw2();

    TTree* tFull = (TTree*)f->Get("tree");
    TFile* dummy = new TFile("dummy.root","RECREATE");
    TTree* t = (TTree*)tFull->CopyTree( sample_cut );

    float p_vsMH_s[999];
    float p_tt_bb [999];
    int nPermut_s;
    int nTotInteg_s;
    int nMassPoints;
    float mH_scan[999];
    int perm_to_gen_s[999];
    float weight;

    t->SetBranchAddress("p_vsMH_s",     p_vsMH_s);
    t->SetBranchAddress("p_tt_bb",      p_tt_bb);
    t->SetBranchAddress("nPermut_s",    &nPermut_s);
    t->SetBranchAddress("nTotInteg_s",  &nTotInteg_s);
    t->SetBranchAddress("nMassPoints",  &nMassPoints);
    t->SetBranchAddress("mH_scan",      mH_scan);
    t->SetBranchAddress("weight",       &weight);
    t->SetBranchAddress("perm_to_gen_s",perm_to_gen_s);
    
    Long64_t nentries = t->GetEntries(); 
    cout << "Total entries: " << nentries << endl;

    int counter       = 0;
    int counterQuarks = 0;
    int counterMatch  = 0;

    for (Long64_t i = 0; i < nentries ; i++){

      t->GetEntry(i);
      counter++;

      hTmp->Reset();

      float perm_prob [nPermut_s];
      float perm_match[nPermut_s];      
      for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){
	perm_prob [perm_it] = 0.;
	perm_match[perm_it] = 0;
      }

      int quarks = 0;
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){
	float mH =  mH_scan[mH_it];
	for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){
	  float ME_prob = p_vsMH_s[mH_it*nPermut_s + perm_it];
	  float bb_prob = useBTag ? p_tt_bb[perm_it] : 1.0;
	  float norm    = normXSec ? 1./xsec->Eval( mH ) : 1.0;
	  int match     = perm_to_gen_s[perm_it];

	  double p =  ME_prob*bb_prob*norm;

	  perm_prob [perm_it] += p;
	  perm_match[perm_it] = match;
	  
	  if( match == match1 || match == match2 ) quarks++;

	  hTmp->Fill( mH, p );
	}
      }

      pair<double,double> bestMass = getMaxValue(hTmp);
      double mass = bestMass.first;
      //cout << "Ev." << i << " => mass=" << mass << endl;

      float maxP = 0.;
      int   maxM = 0;
      for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){
	if(perm_prob [perm_it] > maxP ){
	  maxP = perm_prob [perm_it];
	  maxM = perm_match[perm_it];
	}
      }
      
      if(quarks>0){
	counterQuarks++;
	if( maxM%100 == match3 || maxM == match4 ){
	  counterMatch++;
	  //hMass_s->Fill(mass, weight*lumiScale);
	}
	//else
	  //hMass_s->Fill(mass, weight*lumiScale);
      }

      hMass_s->Fill(mass, weight*lumiScale);
    }

    cout << "Total: " << counter << endl;
    cout << "All quarks in acceptance: " << counterQuarks << endl;
    cout << "Max permutation is correct: " << counterMatch << endl;
    cout << " ==> acc. = " << (counter>0 ? float(counterQuarks)/float(counter) : 0.) << endl;
    cout << " ==> eff. = " << (counterQuarks>0 ? float(counterMatch)/float(counterQuarks) : 0. ) << endl;

    if( hErr==0 ){
      hErr = (TH1F*)hMass_s->Clone("hErr");
      hErr->Reset();
      leg->AddEntry(hErr, "MC unc. (stat.)", "L");
    }
    if( sample.find("TTH")==string::npos ) hErr->Add( hMass_s, 1.0);

    if( sample.find("TTJets")!=string::npos ){
      hMass_s->SetFillColor(color);
      if(color==16)
	leg->AddEntry(hMass_s, "t#bar{t} + bb", "F");
      if(color==17)
	leg->AddEntry(hMass_s, "t#bar{t} + b", "F");
      if(color==18)
	leg->AddEntry(hMass_s, "t#bar{t} + jj", "F");
      hTT->Add(hMass_s, 1.0);
    }
    if( sample.find("TTH")!=string::npos ){
      hMass_s->SetFillColor(kRed);
      hTTH->Add(hMass_s, 1.0);
      leg->AddEntry(hMass_s, "t#bar{t}H", "F");
    }
    if( sample.find("TTV")!=string::npos ){
      hMass_s->SetFillColor(kBlue);
      hTTV->Add(hMass_s, 1.0);
      leg->AddEntry(hMass_s, "t#bar{t}V", "F");
    }
    if( sample.find("SingleT")!=string::npos ){
      hMass_s->SetFillColor(kMagenta);
      leg->AddEntry(hMass_s, "Single top", "F");
    }
    if( sample.find("DiBoson")!=string::npos ){
      hMass_s->SetFillColor(kYellow);
      leg->AddEntry(hMass_s, "VV", "F");
    }
    if( sample.find("EWK")!=string::npos ){
      hMass_s->SetFillColor(kGreen);
      leg->AddEntry(hMass_s, "V+jets", "F");
    }

    cout << "Adding " << hMass_s->Integral() << " weighted events to the stack" << endl;
    aStack->Add( hMass_s );

  }

  hErr->GetYaxis()->SetTitle("Events");
  hErr->GetXaxis()->SetTitle("M_{H} estimator [GeV]");
  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  hErr->GetMaximum()*1.45;
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);
  hErr->Draw("HISTE1");

  aStack->Draw("HISTSAME");

  hTTH->SetLineWidth(3);
  hTTH->SetLineColor(kRed);
  hTTH->SetLineStyle(kDashed);
  hTTH->SetFillColor(kRed);
  hTTH->SetFillStyle(3004);
  hTTH->Scale(5.0);
  hTTH->Draw("HISTSAME");
  hTTV->SetLineWidth(3);
  hTTV->SetLineColor(kBlue);
  hTTV->SetLineStyle(kDashed);
  hTTV->SetFillColor(kBlue);
  hTTV->SetFillStyle(3005);
  hTTV->Scale(5.0);
  hTTV->Draw("HISTSAME");
  hErr->Draw("HISTE1SAME");

  leg->AddEntry(hTTH, "t#bar{t}H x 5", "L");
  leg->AddEntry(hTTV, "t#bar{t}V x 5", "L");

  if(plotShapes){
    hTT->SetLineWidth(3);
    hTT->SetLineColor(kBlack);
    hTT->SetLineStyle(kSolid);
    //hTT->SetFillColor(16);
    //hTT->SetFillStyle(3005);

    hTTH->Scale(1./hTTH->Integral());
    hTTV->Scale(1./hTTV->Integral());
    hTT ->Scale(1./hTT ->Integral());

    hTTH->GetYaxis()->SetRangeUser(0., hTTH->GetMaximum()*1.45);
    hTTH->Draw("HIST");
    hTTV->Draw("HISTSAME");
    hTT->Draw("HISTSAME");
    leg->Clear();
    leg->AddEntry(hTTH, "t#bar{t}H",     "F");
    leg->AddEntry(hTTV, "t#bar{t}V",     "F");
    leg->AddEntry(hTT,  "t#bar{t}+jets", "F");
    leg->Draw();
  }
  else{
    leg->Draw();
  }

  //hTmp->Draw();



  if(0){
    if(plotShapes)
      c1->SaveAs("mass_estimator_shapes_"+fname+".png");
    else
      c1->SaveAs("mass_estimator_spectrum_"+fname+".png");
  }

  cout << "Remove dummy file" << endl;
  gSystem->Exec("rm dummy.root");
  return;
}
