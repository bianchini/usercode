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


void produce(string nuisance = "lumi", bool isNuisance = true,
	     int mH = 120, 
	     bool doCombined      = true, 
	     bool combineChannels = false, 
	     bool combineEnergy   = false,
	     bool forceLimits     = false){


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

  vector<string> com;
  com.push_back("7TeV");
  com.push_back("8TeV");

  vector<string> channels;  vector<string> channelsDir;
  channels.push_back("mt");   channelsDir.push_back("muTau");
  //channels.push_back("et"); channelsDir.push_back("eleTau");

  vector<string> categories;
  vector<string> categoriesIndex;
  categories.push_back("0jet_low");   categoriesIndex.push_back("0");
  categories.push_back("0jet_high");  categoriesIndex.push_back("1");
  categories.push_back("boost_low");  categoriesIndex.push_back("2");
  categories.push_back("boost_high"); categoriesIndex.push_back("3");
  categories.push_back("vbf");        categoriesIndex.push_back("5");
  
  int totalBins  = channels.size()*categories.size()*com.size();
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

  string combine = "combineCards.py ";
  int inter = 0;
  
  for(unsigned int en = 0; en<com.size(); en++){
    
    for(unsigned int ch = 0; ch<channels.size(); ch++){
      
      gSystem->Exec( Form("cp ../test/root/htt_%s.inputs-sm-%s.root ../test/htt_%s.inputs-sm-%s_tmp.root", channels[ch].c_str(), com[en].c_str(), channels[ch].c_str(), com[en].c_str()) );
      TFile* fileTmp = new TFile(Form("../test/htt_%s.inputs-sm-%s_tmp.root",channels[ch].c_str(), com[en].c_str()),"UPDATE");
      
      if(!combineChannels){
	inter = 0;
	combine = "combineCards.py ";
      }
      
      for(unsigned int ca = 0; ca<categories.size(); ca++){
	
	inter++;
	
	gSystem->Exec( Form("cp ../test/datacards/htt_%s_%s_%s-%d.txt htt_%s_%s_%s-%d_tmp0.txt",
			    channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH,
			    channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH) );	
	if(isNuisance){
	  gSystem->Exec( Form("sed  's/%s/%s/g' htt_%s_%s_%s-%d_tmp0.txt > htt_%s_%s_%s-%d_tmp1.txt", 
			      nuisance.c_str(), 
			      (nuisance+"_"+channelsDir[ch]+"_"+categories[ca]+"_"+com[en]).c_str(),
			      channels[ch].c_str(),categoriesIndex[ca].c_str(),com[en].c_str(), mH,
			      channels[ch].c_str(),categoriesIndex[ca].c_str(),com[en].c_str(), mH));
	  gSystem->Exec( Form("sed  's/htt_%s.inputs-sm-%s.root/htt_%s.inputs-sm-%s_tmp.root/g' htt_%s_%s_%s-%d_tmp1.txt > htt_%s_%s_%s-%d_tmp.txt",
			      channels[ch].c_str(), com[en].c_str(),
			      channels[ch].c_str(), com[en].c_str(),
			      channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH,
			      channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH ));
	  gSystem->Exec( "rm *tmp0*txt" );  
	  gSystem->Exec( "rm *tmp1*txt" );
	}
	else{
	  gSystem->Exec( Form("sed  's/htt_%s.inputs-sm-%s.root/htt_%s.inputs-sm-%s_tmp.root/g' htt_%s_%s_%s-%d_tmp0.txt > htt_%s_%s_%s-%d_tmp.txt",
			      channels[ch].c_str(), com[en].c_str(),
			      channels[ch].c_str(), com[en].c_str(),
			      channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH,
			      channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH ));
	  gSystem->Exec( "rm *tmp0*txt" );  
	}
	
	TString dirName( Form("%s_%s",channelsDir[ch].c_str(), categories[ca].c_str()) );
	TDirectory* dir = (TDirectory*)fileTmp->Get( dirName.Data() );
	
	if( dir!=0 ){
	  
	  dir->cd();
	  
	  TIter nextkey(dir->GetListOfKeys());
	  TH1F *key;
	  while ( (key = (TH1F*)nextkey()) ) {
	    //cout << string(key->GetName()) << endl;
	    string name(key->GetName());
	    if( isNuisance && name.find(nuisance)!=string::npos && name.find(nuisance+"_"+channelsDir[ch]+"_"+categories[ca]+"_"+com[en])==string::npos ){
	      name.replace( name.find(nuisance), nuisance.size(), 
			    nuisance+"_"+channelsDir[ch]+"_"+categories[ca]+"_"+com[en] );
	      //cout << " ==> " << name << endl;
	      TH1F* h = (TH1F*)dir->FindObjectAny( key->GetName());
	      h->Write(name.c_str());
	    }
	  }
	  
	}
	else{
	  cout << "Could not find directory " << string(Form("%s_%s",channelsDir[ch].c_str(), categories[ca].c_str())) << endl;
	}
	
	combine = combine+string(Form(" Name%d=htt_%s_%s_%s-%d_tmp.txt",inter,channels[ch].c_str(),categoriesIndex[ca].c_str(),com[en].c_str(), mH));
	
      }
      
      fileTmp->Close();
      delete fileTmp;
      
      if(!combineChannels && !combineEnergy){
	combine = combine + " > " + string(Form("htt_%s_comb_%s-%d_tmp.txt",channels[ch].c_str(),com[en].c_str(), mH));
	cout << combine << endl;
	gSystem->Exec( combine.c_str() );
	gSystem->Exec( "mv *tmp* ../test/" );
      }
      
    }
    
    if(combineChannels && !combineEnergy){
      combine = combine + " > " + string(Form("htt_comb_comb_%s-%d_tmp.txt",com[en].c_str(), mH ));
      cout << combine << endl;
      gSystem->Exec( combine.c_str() );
      gSystem->Exec( "mv *tmp* ../test/" );
    }      
  }
  
  if(combineChannels && combineEnergy){
    combine = combine + " > " + string(Form("htt_comb_comb_comb-%d_tmp.txt",mH));
    cout << combine << endl;
    gSystem->Exec( combine.c_str() );
    gSystem->Exec( "mv *tmp* ../test/" );
  }
  
  
  //return;

  float maxY = 0.;

    for(unsigned int en = 0; en<com.size(); en++){

      for(unsigned int ch = 0; ch<channels.size(); ch++){
	
	for(unsigned int ca = 0; ca<categories.size(); ca++){
	  
	  binCounter++;
	  
	  if(!doCombined)
	    gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_%s_%s_%s-%d_tmp.txt",         mH, channels[ch].c_str(),categoriesIndex[ca].c_str(),com[en].c_str(), mH) );
	  else{
	    if(!combineChannels && !combineEnergy) 
	      gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_%s_comb_%s-%d_tmp.txt",     mH, channels[ch].c_str(), com[en].c_str(), mH) );
	    else if(combineChannels && !combineEnergy)
	      gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_comb_comb_%s-%d_tmp.txt",   mH, com[en].c_str(), mH) );
	    else
	      gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_comb_comb_comb-%d_tmp.txt", mH, mH) );
	  }
	  

	  TFile* file = new TFile("./mlfit.root","READ");
	  if(!file || file->IsZombie()) continue;
	  
	  float value; float error;
	  
	  string uncorrelatedNuisanceName = isNuisance ? 
	    nuisance+"_"+channelsDir[ch]+"_"+categories[ca]+"_"+com[en] : nuisance;
	  
	  RooFitResult* fit_b = (RooFitResult*)file->Get("fit_b");
	  if(!fit_b){
	    cout << "No fit_b available..." << endl; continue;
	  }
	  RooArgSet fit_b_list(fit_b->floatParsFinal());
	  RooRealVar* nuisanceFitB   = 0;
	  if( !fit_b_list.find(Form("%s",uncorrelatedNuisanceName.c_str() ))){
	    cout << "Nuisance " << nuisance << " is not available for b fit..." << endl;
	    value = 0.0; error = 0.0;
	  }
	  else{
	    nuisanceFitB = (RooRealVar*)(&fit_b_list[uncorrelatedNuisanceName.c_str()]);
	    value = nuisanceFitB->getVal();
	    error = nuisanceFitB->getError();
	  }
	  histoB->SetBinContent(binCounter, value);
	  histoB->SetBinError(binCounter,   error);
	  histoB->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s_%s",channelsDir[ch].c_str(),categories[ca].c_str(),com[en].c_str()));
	  
	  if(value+error>maxY) maxY = value+error;
	  

	  RooFitResult* fit_s = (RooFitResult*)file->Get("fit_s");
	  if(!fit_s){
	    cout << "No fit_s available..." << endl; continue;
	  }
	  RooArgSet fit_s_list(fit_s->floatParsFinal());
	  RooRealVar* nuisanceFitS   = 0;
	  if(!fit_s_list.find(Form("%s", uncorrelatedNuisanceName.c_str() ))){
	    cout << "Nuisance " << nuisance << " is not available for s+b fit..." << endl;
	    value = 0.0; error = 0.0;
	  }
	  else{
	    nuisanceFitS = (RooRealVar*)(&fit_s_list[uncorrelatedNuisanceName.c_str()]);
	    value = nuisanceFitS->getVal();
	    error = nuisanceFitS->getError();
	  }
	  histoS->SetBinContent(binCounter, value);
	  histoS->SetBinError(binCounter,   error);
	  histoS->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s_%s",channelsDir[ch].c_str(),categories[ca].c_str(),com[en].c_str()));
	  
	  if(value+error>maxY) maxY = value+error;
	  
	  file->Close();
	  delete file;
	}
      }

    }
  
  float MaxY = TMath::Max( maxY, float(2.0/(1.10)));
  if(forceLimits) MaxY = 1.5;

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

  if(doCombined){
    if(!combineChannels && !combineEnergy) 
      c1->SaveAs(Form("nuisance_%s_mH%d_COMB_PERCHANNEL.png",  nuisance.c_str(), mH ));
    else if(combineChannels && !combineEnergy) 
      c1->SaveAs(Form("nuisance_%s_mH%d_COMB_ALLCHANNELS.png", nuisance.c_str(), mH ));
    else
      c1->SaveAs(Form("nuisance_%s_mH%d_COMB_ALLCHANNELS_ALLENERGY.png", nuisance.c_str(), mH ));
  }
  else{
    c1->SaveAs(Form("nuisance_%s_mH%d_PERCATEGORY.png", nuisance.c_str(), mH));
  }

  gSystem->Exec( "rm ../test/*tmp*" );
  return;
}


  

void produceAll(){

  //produce("CMS_scale_j_8TeV"   ,125, true, true, true);
  //produce("CMS_scale_met_8TeV" ,125, true, true, true);
  //produce("CMS_scale_met_7TeV" ,125, true, true, true);
  //produce("CMS_eff_t"          ,125, true, true, true);
  produce("CMS_scale_t",            1, 125 ,true, true, true);
  //produce("r"                    ,0, 125, true, true, true, true);

}



int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}
