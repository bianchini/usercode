#include "FWCore/FWLite/interface/AutoLibraryLoader.h"


#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"

#include "TauAnalysis/CandidateTools/interface/SVfitVMlineShapeIntegral.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

using namespace std;

void lineShape(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetGrid(0,0);
  c1->SetTicky();
  c1->SetObjectStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit(0111);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);

  //cout << (*lmLineShape)(atoi(polarization),atof(z)) << endl;

  fwlite::TFileService fs = fwlite::TFileService("shapes.root");

  cout << "Now doing for Rho L" << endl;
  TH1F* hRhoLMinus = new TH1F("hRhoLMinus","#tau_{L}^{-1} #rightarrow #rho_{L} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  TH1F* hRhoLPlus  = new TH1F("hRhoLPlus", "#tau_{R}^{-1} #rightarrow #rho_{L} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  SVfitVMlineShapeIntegral* lmLineShapeRhoL = 
    new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMrho, 
				 SVfitVMlineShapeIntegrand::kVMlongitudinalPol,
				 true);
  for(int i = 1; i <= 100; i++){
    hRhoLMinus->SetBinContent( i, (*lmLineShapeRhoL)(-1,i*0.01));
    hRhoLPlus->SetBinContent(  i, (*lmLineShapeRhoL)(+1,i*0.01));
  }
  c1->cd();
  hRhoLMinus->Draw();
  c1->SaveAs("hRhoLMinus.pdf");
  hRhoLPlus->Draw();
  c1->SaveAs("hRhoLPlus.pdf");
  delete lmLineShapeRhoL;


  cout << "Now doing for Rho T" << endl;
  TH1F* hRhoTMinus = new TH1F("hRhoTMinus","#tau_{L}^{-1} #rightarrow #rho_{T} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  TH1F* hRhoTPlus  = new TH1F("hRhoTPlus", "#tau_{R}^{-1} #rightarrow #rho_{T} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  SVfitVMlineShapeIntegral* lmLineShapeRhoT = 
    new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMrho, 
				 SVfitVMlineShapeIntegrand::kVMtransversePol,
				 true);
  for(int i = 1; i <= 100; i++){
    hRhoTMinus->SetBinContent( i, (*lmLineShapeRhoT)(-1,i*0.01));
    hRhoTPlus->SetBinContent(  i, (*lmLineShapeRhoT)(+1,i*0.01));
  }
  c1->cd();
  hRhoTMinus->Draw();
  c1->SaveAs("hRhoTMinus.pdf");
  hRhoTPlus->Draw();
  c1->SaveAs("hRhoTPlus.pdf");
  delete lmLineShapeRhoT;

  cout << "Now doing for a1 L" << endl;
  TH1F* hA1LMinus = new TH1F("hA1LMinus","#tau_{L}^{-1} #rightarrow a_{1}_{L} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  TH1F* hA1LPlus  = new TH1F("hA1LPlus", "#tau_{R}^{-1} #rightarrow a_{1}_{L} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  SVfitVMlineShapeIntegral* lmLineShapeA1L = 
    new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMa1Charged, 
				 SVfitVMlineShapeIntegrand::kVMlongitudinalPol,
				 true);
  for(int i = 1; i <= 100; i++){
    hA1LMinus->SetBinContent( i, (*lmLineShapeA1L)(-1,i*0.01));
    hA1LPlus->SetBinContent(  i, (*lmLineShapeA1L)(+1,i*0.01));
  }
  c1->cd();
  hA1LMinus->Draw();
  c1->SaveAs("hA1LMinus.pdf");
  hA1LPlus->Draw();
  c1->SaveAs("hA1LPlus.pdf");
  delete lmLineShapeA1L;

  cout << "Now doing for a1 T" << endl;
  TH1F* hA1TMinus = new TH1F("hA1TMinus","#tau_{L}^{-1} #rightarrow a_{1}_{T} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  TH1F* hA1TPlus  = new TH1F("hA1TPlus", "#tau_{R}^{-1} #rightarrow a_{1}_{T} ; x_{h} ; d#Gamma/dx_{h}",100,0,1);
  SVfitVMlineShapeIntegral* lmLineShapeA1T = 
    new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMa1Charged, 
				 SVfitVMlineShapeIntegrand::kVMtransversePol,
				 true);
  for(int i = 1; i <= 100; i++){
    hA1TMinus->SetBinContent( i, (*lmLineShapeA1T)(-1,i*0.01));
    hA1TPlus->SetBinContent(  i, (*lmLineShapeA1T)(+1,i*0.01));
  }
  c1->cd();
  hA1TMinus->Draw();
  c1->SaveAs("hA1TMinus.pdf");
  hA1TPlus->Draw();
  c1->SaveAs("hA1TPlus.pdf");
  delete lmLineShapeA1T;

 


}



int main(int argc, const char* argv[])
{

  std::cout << "plotMuTau()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  lineShape();

}
