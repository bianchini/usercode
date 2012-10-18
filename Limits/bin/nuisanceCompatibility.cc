#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TKey.h"
#include "TObjString.h"
#include "TObject.h"
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
#include "TCollection.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooLognormal.h"


using namespace RooFit;
using namespace std;


void produce(string nuisance = "lumi", int mH = 120, bool doCombined = true){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
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
  gStyle->SetTitleOffset(1.1,"y");

  TLegend* leg = new TLegend(0.70,0.72,0.85,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  vector<string> channels;
  channels.push_back("muTau");
//   channels.push_back("eTau");

  vector<string> categories;
  categories.push_back("0jet_low");
  categories.push_back("0jet_high");
  categories.push_back("boost_low");
  categories.push_back("boost_high");
  categories.push_back("vbf");
  //   categories.push_back("comb");
  
  int totalBins  = channels.size()*categories.size();
  int binCounter = 0;

  TH1F* histoB = new TH1F("histoB",Form("Nuisance O = %s, for m_{H}=%d GeV;;[(O_{fit}-O_{in}) #pm #sigma_{fit}]/#sigma_{in}",
					nuisance.c_str(), mH),totalBins, 0, totalBins);
  histoB->SetMarkerStyle(kFullCircle);
  histoB->SetMarkerSize(1.6);
  histoB->SetMarkerColor(kRed);
  histoB->SetLineColor(kRed);
  histoB->SetLineWidth(2);

  TH1F* histoS = new TH1F("histoS",Form("Nuisance O = %s, for m_{H}=%d GeV;;[(O_{fit}-O_{in}) #pm #sigma_{fit}]/#sigma_{in}",
					nuisance.c_str(), mH),totalBins, 0, totalBins);
  histoS->SetMarkerStyle(kOpenSquare);
  histoS->SetMarkerSize(1.6);
  histoS->SetMarkerColor(kBlue);
  histoS->SetFillColor(kBlue);
  histoS->SetFillStyle(3003);
  histoS->SetLineColor(kBlue);
  histoS->SetLineWidth(2);

  if(doCombined){

    gSystem->Exec( "cp test/muTauSM.root test/muTauSM_tmp.root" );
    TFile* fileTmp = new TFile("test/muTauSM_tmp.root","UPDATE");
    
    string combine = "combineCards.py ";
    int inter = 0;
    
    for(unsigned int ch = 0; ch<channels.size(); ch++){
      inter = 0;
      
      for(unsigned int ca = 0; ca<categories.size(); ca++){
	
	inter++;
	
	gSystem->Exec( Form("cp test/%s_%s_mH%d.txt %s_%s_mH%d_tmp0.txt",channels[ch].c_str(),categories[ca].c_str(),mH,
			    channels[ch].c_str(),categories[ca].c_str(),mH) );
	
	gSystem->Exec( Form("sed  's/%s/%s/g' %s_%s_mH%d_tmp0.txt > %s_%s_mH%d_tmp1.txt", nuisance.c_str(), 
			    (nuisance+"_"+channels[ch]+"_"+categories[ca]).c_str(),
			    channels[ch].c_str(),categories[ca].c_str(),mH,
			    channels[ch].c_str(),categories[ca].c_str(),mH));
	gSystem->Exec( Form("sed  's/%sSM.root/%sSM_tmp.root/g' %s_%s_mH%d_tmp1.txt > %s_%s_mH%d_tmp.txt",
			    channels[ch].c_str(),channels[ch].c_str(),
			    channels[ch].c_str(),categories[ca].c_str(),mH,
			    channels[ch].c_str(),categories[ca].c_str(),mH));
	gSystem->Exec( "rm *tmp0*txt" );  gSystem->Exec( "rm *tmp1*txt" );
	
	TString dirName( Form("%s_%s",channels[ch].c_str(), categories[ca].c_str()) );
	TDirectory* dir = (TDirectory*)fileTmp->Get( dirName.Data() );
	
	if( dir!=0 ){
	  
	  dir->cd();
	  
	  TIter nextkey(dir->GetListOfKeys());
	  TH1F *key;
	  while (key = (TH1F*)nextkey()) {
	    cout << string(key->GetName()) << endl;
	    string name(key->GetName());
	    if(  name.find(nuisance)!=string::npos && name.find(nuisance+"_"+channels[ch]+"_"+categories[ca])==string::npos ){
	      name.replace( name.find(nuisance), nuisance.size()/*+channels[ch].size()+categories[ca].size()+2*/, 
			    nuisance+"_"+channels[ch]+"_"+categories[ca] );
	      //cout << " ==> " << name << endl;
	      TH1F* h = (TH1F*)dir->FindObjectAny( key->GetName());
	      h->Write(name.c_str());
	    }
	}
	  
	}
	else{
	  cout << "Could not find directory " << string(Form("%s_%s",channels[ch].c_str(), categories[ca].c_str())) << endl;
	}

	combine = combine+string(Form(" Name%d=%s_%s_mH%d_tmp.txt",inter,channels[ch].c_str(),categories[ca].c_str(),mH));
	
      }
      
      combine = combine + " > " + string(Form("%s_comb_mH%d_tmp.txt",channels[ch].c_str(),mH));
      
      cout << combine << endl;
      
      gSystem->Exec( combine.c_str() );
      gSystem->Exec( "mv *tmp* test/" );
      
      
    }
    
    //return;
    fileTmp->Close();
  }


  float maxY = 0.;


  if(doCombined){

    for(unsigned int ch = 0; ch<channels.size(); ch++){

      for(unsigned int ca = 0; ca<categories.size(); ca++){
	
	binCounter++;
	
	gSystem->Exec( Form("combine -M MaxLikelihoodFit test/%s_comb_mH%d_tmp.txt",channels[ch].c_str(),mH) );
	TFile* file = new TFile("./mlfit.root","READ");
	if(!file || file->IsZombie()) continue;
	
	float value; float error;
	

	RooFitResult* fit_b = (RooFitResult*)file->Get("fit_b");
	if(!fit_b){
	  cout << "No fit_b available..." << endl; continue;
	}
	RooArgSet fit_b_list(fit_b->floatParsFinal());
	RooRealVar* nuisanceFitB   = 0;
	if(!fit_b_list.find(Form("%s",(nuisance+"_"+channels[ch]+"_"+categories[ca]).c_str()))){
	  cout << "Nusiance " << nuisance << " is not available..." << endl;
	  value = 0.0; error = 0.0;
	}
	else{
	  nuisanceFitB = (RooRealVar*)(&fit_b_list[Form("%s",(nuisance+"_"+channels[ch]+"_"+categories[ca]).c_str())]);
	  value = nuisanceFitB->getVal();
	  error = nuisanceFitB->getError();
	}
	histoB->SetBinContent(binCounter, value);
	histoB->SetBinError(binCounter,   error);
	histoB->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s",channels[ch].c_str(),categories[ca].c_str()));
	
	if(value+error>maxY) maxY = value+error;
	
	RooFitResult* fit_s = (RooFitResult*)file->Get("fit_s");
	if(!fit_s){
	  cout << "No fit_s available..." << endl; continue;
	}
	RooArgSet fit_s_list(fit_s->floatParsFinal());
	RooRealVar* nuisanceFitS   = 0;
	if(!fit_s_list.find(Form("%s",(nuisance+"_"+channels[ch]+"_"+categories[ca]).c_str()))){
	  cout << "Nusiance " << nuisance << " is not available..." << endl;
	  value = 0.0; error = 0.0;
	}
	else{
	  nuisanceFitS = (RooRealVar*)(&fit_s_list[Form("%s",(nuisance+"_"+channels[ch]+"_"+categories[ca]).c_str())]);
	  value = nuisanceFitS->getVal();
	  error = nuisanceFitS->getError();
	}
	histoS->SetBinContent(binCounter, value);
	histoS->SetBinError(binCounter,   error);
	histoS->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s",channels[ch].c_str(),categories[ca].c_str()));
	
	if(value+error>maxY) maxY = value+error;
	
	file->Close();
	delete file;
      }
    }

  }
  else{

    for(unsigned int ch = 0; ch<channels.size(); ch++){

      for(unsigned int ca = 0; ca<categories.size(); ca++){
	
	binCounter++;
	
	gSystem->Exec( Form("combine -M MaxLikelihoodFit test/%s_%s_mH%d.txt",channels[ch].c_str(),categories[ca].c_str(),mH) );
	TFile* file = new TFile("./mlfit.root","READ");
	if(!file || file->IsZombie()) continue;
	
	float value; float error;
	

	RooFitResult* fit_b = (RooFitResult*)file->Get("fit_b");
	if(!fit_b){
	  cout << "No fit_b available..." << endl; continue;
	}
	RooArgSet fit_b_list(fit_b->floatParsFinal());
	RooRealVar* nuisanceFitB   = 0;
	if(!fit_b_list.find(Form("%s",nuisance.c_str()))){
	  cout << "Nusiance " << nuisance << " is not available..." << endl;
	  value = 0.0; error = 0.0;
	}
	else{
	  nuisanceFitB = (RooRealVar*)(&fit_b_list[Form("%s",nuisance.c_str())]);
	  value = nuisanceFitB->getVal();
	  error = nuisanceFitB->getError();
	}
	histoB->SetBinContent(binCounter, value);
	histoB->SetBinError(binCounter,   error);
	histoB->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s",channels[ch].c_str(),categories[ca].c_str()));
	
	if(value+error>maxY) maxY = value+error;
	
	RooFitResult* fit_s = (RooFitResult*)file->Get("fit_s");
	if(!fit_s){
	  cout << "No fit_s available..." << endl; continue;
	}
	RooArgSet fit_s_list(fit_s->floatParsFinal());
	RooRealVar* nuisanceFitS   = 0;
	if(!fit_s_list.find(Form("%s",nuisance.c_str()))){
	  cout << "Nusiance " << nuisance << " is not available..." << endl;
	  value = 0.0; error = 0.0;
	}
	else{
	  nuisanceFitS = (RooRealVar*)(&fit_s_list[Form("%s",nuisance.c_str())]);
	  value = nuisanceFitS->getVal();
	  error = nuisanceFitS->getError();
	}
	histoS->SetBinContent(binCounter, value);
	histoS->SetBinError(binCounter,   error);
	histoS->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s",channels[ch].c_str(),categories[ca].c_str()));
	
	if(value+error>maxY) maxY = value+error;
	
	file->Close();
	delete file;
      }
    }
  }
  
  float MaxY = TMath::Max( maxY, float(2.0/(1.10)));

  c1->cd();
  histoB->SetAxisRange(-MaxY*1.10, +MaxY*1.10 ,"Y");
  histoB->Draw("PE");
  histoS->Draw("PSAMEE2");
  histoB->Draw("PESAME");

  leg->AddEntry(histoB,"b fit","P");
  leg->AddEntry(histoS,"s+b fit","F");
  leg->Draw();

  TF1* line = new TF1("line","0",histoB->GetXaxis()->GetXmin(),histoB->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  c1->SaveAs(Form("nuisance_%s.png", nuisance.c_str()));


  return;
}


  

void produceAll(){

  //produce("r",120);
 //  produce("lumi",120);
  produce("CMS_scale_t",120, false);
//   produce("CMS_eff_t",120);
//   produce("CMS_htt_zttNorm",120);
//   produce("CMS_scale_j",120);
//   produce("CMS_scale_met",120);
//   produce("CMS_htt_muTau_ZJetFakeTau"   ,120);
//   produce("CMS_htt_muTau_ZLeptonFakeTau",120);
//   produce("CMS_htt_eTau_ZJetFakeTau"   ,120);
//   produce("CMS_htt_eTau_ZLeptonFakeTau",120);

}




int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}
