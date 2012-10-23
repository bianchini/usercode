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
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"

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


void produce(string nuisance          = "lumi", 
	     bool isNuisance          = true,
	     int mH                   = 120, 
	     bool doCombined          = true,
	     vector<string>& channels = *(new  vector<string>()),
	     vector<string>& com      = *(new  vector<string>()),
	     TString title            = "",
	     bool forceLimits         = false
	     ){

  bool combineEnergy   = com.size()>1 ;
  bool combineChannels = (channels.size()>1) || combineEnergy;

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
  leg->SetHeader(title);

  //vector<string> com;
  //com.push_back("7TeV");
  //com.push_back("8TeV");
  //vector<string> channels;  
  //channels.push_back("mt"); channelsDir.push_back("muTau");
  //channels.push_back("et"); channelsDir.push_back("eleTau");
  //channels.push_back("tt"); channelsDir.push_back("tauTau");
  //channels.push_back("em"); channelsDir.push_back("emu");
  bool tauTauIsThere = false;

  for(unsigned int en = 0; en<com.size(); en++){
    cout << com[en] << endl;
  }
  vector<string> channelsDir;
  for(unsigned int ch = 0; ch<channels.size(); ch++){
    cout << channels[ch] << endl;
    if(channels[ch].find("mt")!=string::npos) channelsDir.push_back("muTau");
    if(channels[ch].find("et")!=string::npos) channelsDir.push_back("eleTau");
    if(channels[ch].find("tt")!=string::npos){
      channelsDir.push_back("tauTau");
      tauTauIsThere = true;
    }
    if(channels[ch].find("em")!=string::npos) channelsDir.push_back("emu");
  }

  vector<string> categories;
  vector<string> categoriesIndex;
  categories.push_back("0jet_low");   categoriesIndex.push_back("0");
  categories.push_back("0jet_high");  categoriesIndex.push_back("1");
  categories.push_back("boost_low");  categoriesIndex.push_back("2");
  categories.push_back("boost_high"); categoriesIndex.push_back("3");
  categories.push_back("vbf");        categoriesIndex.push_back("5");
  
  int totalBins  = !tauTauIsThere ? channels.size()*categories.size()*com.size() : (channels.size()-1)*categories.size()*com.size() + 2;
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

  float X[totalBins];
  float Y[totalBins];
  float XL[totalBins];
  float XH[totalBins];
  float YL[totalBins];
  float YH[totalBins];

  float X2[totalBins];
  float Y2[totalBins];
  float XL2[totalBins];
  float XH2[totalBins];
  float YL2[totalBins];
  float YH2[totalBins];


  string combine = "combineCards.py ";
  int inter = 0;
  
  for(unsigned int en = 0; en<com.size(); en++){
    
    for(unsigned int ch = 0; ch<channels.size(); ch++){

      if(channels[ch].find("tt")!=string::npos && com[en].find("7TeV")!=string::npos) continue;
      
      gSystem->Exec( Form("cp ../test/root/htt_%s.inputs-sm-%s.root ../test/htt_%s.inputs-sm-%s_tmp.root", channels[ch].c_str(), com[en].c_str(), channels[ch].c_str(), com[en].c_str()) );
      TFile* fileTmp = new TFile(Form("../test/htt_%s.inputs-sm-%s_tmp.root",channels[ch].c_str(), com[en].c_str()),"UPDATE");
      
      if(!combineChannels){
	inter = 0;
	combine = "combineCards.py ";
      }
      
      for(unsigned int ca = 0; ca<categories.size(); ca++){

	if(channels[ch].find("tt")!=string::npos && ca>1) continue;
	string categoryName = categories[ca];
	if(channels[ch].find("tt")!=string::npos){
	  if( ca == 0) categoryName = "boost";
	  if( ca == 1) categoryName = "vbf";
	}

	
	inter++;
	
	gSystem->Exec( Form("cp ../test/datacards/htt_%s_%s_%s-%d.txt htt_%s_%s_%s-%d_tmp0.txt",
			    channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH,
			    channels[ch].c_str(),categoriesIndex[ca].c_str(), com[en].c_str(), mH) );	
	if(isNuisance){
	  gSystem->Exec( Form("sed  's/%s/%s/g' htt_%s_%s_%s-%d_tmp0.txt > htt_%s_%s_%s-%d_tmp1.txt", 
			      nuisance.c_str(), 
			      (nuisance+"_"+channelsDir[ch]+"_"+categoryName+"_"+com[en]).c_str(),
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
       

	TString dirName( Form("%s_%s",channelsDir[ch].c_str(), categoryName.c_str()) );
	TDirectory* dir = (TDirectory*)fileTmp->Get( dirName.Data() );
	
	if( dir!=0 ){
	  
	  dir->cd();
	  
	  TIter nextkey(dir->GetListOfKeys());
	  TH1F *key;
	  while ( (key = (TH1F*)nextkey()) ) {
	    //cout << string(key->GetName()) << endl;
	    string name(key->GetName());
	    if( isNuisance && name.find(nuisance)!=string::npos && name.find(nuisance+"_"+channelsDir[ch]+"_"+categoryName+"_"+com[en])==string::npos ){
	      name.replace( name.find(nuisance), nuisance.size(), 
			    nuisance+"_"+channelsDir[ch]+"_"+categoryName+"_"+com[en] );
	      //cout << " ==> " << name << endl;
	      TH1F* h = (TH1F*)dir->FindObjectAny( key->GetName());
	      h->Write(name.c_str());
	    }
	  }
	  
	}
	else{
	  cout << "Could not find directory " << string(Form("%s_%s",channelsDir[ch].c_str(), categoryName.c_str())) << endl;
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

  if(doCombined && combineChannels && combineEnergy)
    gSystem->Exec( Form("combine -M MaxLikelihoodFit  --robustFit 1  --rMin -30 --rMax 30 -m %d ../test/htt_comb_comb_comb-%d_tmp.txt", mH, mH) );
  
  
  for(unsigned int en = 0; en<com.size(); en++){
    
    if(doCombined && combineChannels && !combineEnergy)
      gSystem->Exec( Form("combine -M MaxLikelihoodFit --robustFit 1  --rMin -30 --rMax 30 -m %d ../test/htt_comb_comb_%s-%d_tmp.txt",   mH, com[en].c_str(), mH) );
    
    for(unsigned int ch = 0; ch<channels.size(); ch++){
      
      if(channels[ch].find("tt")!=string::npos && com[en].find("7TeV")!=string::npos) continue;
      
      if(doCombined && !combineChannels && !combineEnergy)
	gSystem->Exec( Form("combine -M MaxLikelihoodFit --robustFit 1 --rMin -30 --rMax 30 -m %d ../test/htt_%s_comb_%s-%d_tmp.txt",     mH, channels[ch].c_str(), com[en].c_str(), mH) );
      
      
      for(unsigned int ca = 0; ca<categories.size(); ca++){
	
	if(channels[ch].find("tt")!=string::npos && ca>1) continue;
	string categoryName = categories[ca];
	if(channels[ch].find("tt")!=string::npos){
	  if( ca == 0) categoryName = "boost";
	  if( ca == 1) categoryName = "vbf";
	}
	
	binCounter++;
	
	if(!doCombined)
	  gSystem->Exec( Form("combine -M MaxLikelihoodFit --robustFit 1 --rMin -30 --rMax 30  -m %d ../test/htt_%s_%s_%s-%d_tmp.txt",
			      mH, channels[ch].c_str(),categoriesIndex[ca].c_str(),com[en].c_str(), mH) );
	// else{
	//if(!combineChannels && !combineEnergy) 
	//  gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_%s_comb_%s-%d_tmp.txt",     mH, channels[ch].c_str(), com[en].c_str(), mH) );
	//else if(combineChannels && !combineEnergy)
	//  gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_comb_comb_%s-%d_tmp.txt",   mH, com[en].c_str(), mH) );
	//else
	//  gSystem->Exec( Form("combine -M MaxLikelihoodFit -m %d ../test/htt_comb_comb_comb-%d_tmp.txt", mH, mH) );
	//}
	
	
	TFile* file = new TFile("./mlfit.root","READ");
	if(!file || file->IsZombie()) continue;
	
	float value; float error; float errorLo; float errorHi;
	
	string uncorrelatedNuisanceName = isNuisance ? 
	  nuisance+"_"+channelsDir[ch]+"_"+categoryName+"_"+com[en] : nuisance;
	
	RooFitResult* fit_b = (RooFitResult*)file->Get("fit_b");
	if(!fit_b){
	  cout << "No fit_b available..." << endl; continue;
	}
	RooArgSet fit_b_list(fit_b->floatParsFinal());
	RooRealVar* nuisanceFitB   = 0;
	if( !fit_b_list.find(Form("%s",uncorrelatedNuisanceName.c_str() ))){
	  cout << "Nuisance " << nuisance << " is not available for b fit..." << endl;
	  value = 0.0; error = 0.0; errorLo = 0.0; errorHi = 0.0;
	  X[binCounter-1]  = binCounter-1+0.5;
	  XL[binCounter-1] = 0.50; 
	  XH[binCounter-1] = 0.50; 
	  Y[binCounter-1]  = value;
	  YL[binCounter-1] = errorLo; 
	  YH[binCounter-1] = errorHi; 
	}
	else{
	  nuisanceFitB = (RooRealVar*)(&fit_b_list[uncorrelatedNuisanceName.c_str()]);
	  value   = nuisanceFitB->getVal();
	  error   = nuisanceFitB->getError();
	  errorLo = nuisanceFitB->getErrorLo();
	  errorHi = nuisanceFitB->getErrorHi();
	  X[binCounter-1]  = binCounter-1+0.5;
	  XL[binCounter-1] = 0.50; 
	  XH[binCounter-1] = 0.50; 
	  Y[binCounter-1]  = value;
	  YL[binCounter-1] = -errorLo; 
	  YH[binCounter-1] = errorHi; 
	}
	histoB->SetBinContent(binCounter, value);
	histoB->SetBinError(binCounter,   error);
	histoB->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s_%s",channelsDir[ch].c_str(),categoryName.c_str(),com[en].c_str()));
	
	float newMax = TMath::Max(float(fabs(value)+fabs(errorLo)),float(fabs(value)+fabs(errorHi)));
	if(newMax>maxY) maxY = newMax ;
	
	RooFitResult* fit_s = (RooFitResult*)file->Get("fit_s");
	if(!fit_s){
	  cout << "No fit_s available..." << endl; continue;
	}
	RooArgSet fit_s_list(fit_s->floatParsFinal());
	RooRealVar* nuisanceFitS   = 0;
	if(!fit_s_list.find(Form("%s", uncorrelatedNuisanceName.c_str() ))){
	  cout << "Nuisance " << nuisance << " is not available for s+b fit..." << endl;
	  value = 0.0; error = 0.0; errorLo = 0.0; errorHi = 0.0;
	  X2[binCounter-1]  = binCounter-1+0.5;
	  XL2[binCounter-1] = 0.50; 
	  XH2[binCounter-1] = 0.50; 
	  Y2[binCounter-1]  = value;
	  YL2[binCounter-1] = errorLo; 
	  YH2[binCounter-1] = errorHi; 
	}
	else{
	  nuisanceFitS = (RooRealVar*)(&fit_s_list[uncorrelatedNuisanceName.c_str()]);
	  value = nuisanceFitS->getVal();
	  error = nuisanceFitS->getError();
	  errorLo = nuisanceFitS->getErrorLo();
	  errorHi = nuisanceFitS->getErrorHi();
	  X2[binCounter-1]  = binCounter-1+0.5;
	  XL2[binCounter-1] = 0.50; 
	  XH2[binCounter-1] = 0.50; 
	  Y2[binCounter-1]  = value;
	  YL2[binCounter-1] = -errorLo; 
	  YH2[binCounter-1] = errorHi; 
	}
	histoS->SetBinContent(binCounter, value);
	histoS->SetBinError(binCounter,   error);
	histoS->GetXaxis()->SetBinLabel(binCounter,Form("%s_%s_%s",channelsDir[ch].c_str(),categoryName.c_str(),com[en].c_str()));
	
	float newMax = TMath::Max(float(fabs(value)+fabs(errorLo)),float(fabs(value)+fabs(errorHi)));
	if(newMax>maxY) maxY = newMax ;

	file->Close();
	delete file;
      }
    }
    
  }
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(histoB->GetTitle());

  TGraphAsymmErrors* background = new TGraphAsymmErrors(histoB->GetNbinsX(), X, Y, XL ,XH , YL, YH);
  background->SetMarkerStyle(kFullCircle);
  background->SetMarkerSize(1.6);
  background->SetMarkerColor(kRed);
  background->SetLineColor(kRed);
  background->SetLineWidth(2);
  TGraphAsymmErrors* signal = new TGraphAsymmErrors(histoB->GetNbinsX(), X2, Y2, XL2 ,XH2 , YL2, YH2);
  signal->SetMarkerStyle(kOpenSquare);
  signal->SetMarkerSize(1.6);
  signal->SetMarkerColor(kBlue);
  signal->SetFillStyle(3003);
  signal->SetLineColor(kBlue);  
  signal->SetFillColor(kBlue);
  signal->SetLineWidth(2);
  mg->Add(background);
  mg->Add(signal, "pe2");

  float MaxY = TMath::Max( maxY, float(2.0/(1.10)));
  if(forceLimits) MaxY = 1.5;

  c1->cd();
  mg->Draw("ap");
  //histoB->SetAxisRange(-MaxY*1.10, +MaxY*1.10 ,"Y");
  //histoB->Draw("PESAME");
  //histoS->Draw("PSAMEE2");
  //histoB->Draw("PESAME");


  gPad->Modified();
  mg->GetXaxis()->SetLimits(0,totalBins);
  mg->GetYaxis()->SetTitleOffset(0.97);
  mg->SetMinimum(-MaxY*1.10);
  mg->SetMaximum(+MaxY*1.10);

  float arrayX[totalBins+1];
  for(unsigned int i = 0; i <= (unsigned int)(totalBins) ; i++){
    if(i<(unsigned int)(totalBins)) arrayX[i] = X[i]-0.5; 
    else  arrayX[i]  = X[i-1]+0.5; 
    //cout << arrayX[i] << endl;
  }

  mg->GetXaxis()->Set(totalBins,arrayX);
  mg->GetXaxis()->SetTitle(histoB->GetXaxis()->GetTitle());
  mg->GetYaxis()->SetTitle(histoB->GetYaxis()->GetTitle());

  for(unsigned int i = 1; i <= (unsigned int)(totalBins) ; i++ )
    mg->GetXaxis()->SetBinLabel(i,histoB->GetXaxis()->GetBinLabel(i));


  leg->AddEntry(histoB,"b fit","P");
  leg->AddEntry(histoS,"s+b fit","F");
  leg->Draw();

  TF1* line = new TF1("line","0",histoB->GetXaxis()->GetXmin(),histoB->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  string pngName = "";
  for(unsigned int en = 0; en<com.size(); en++)      pngName = pngName+com[en];
  pngName = pngName+"-";
  for(unsigned int ch = 0; ch<channels.size(); ch++) pngName = pngName+channels[ch];
  
  if(doCombined){
    c1->SaveAs(Form("nuisance_%s_mH%d_COMB_%s.png",  nuisance.c_str(), mH, pngName.c_str() ));
  }
  else{
    c1->SaveAs(Form("nuisance_%s_mH%d_%s.png",  nuisance.c_str(), mH, pngName.c_str() ));
  }

  TFile* outFile = new TFile("nuisance.root","UPDATE");
  outFile->mkdir(Form("nuisance_%s_mH%d_COMB_%s",  nuisance.c_str(), mH, pngName.c_str() ));
  outFile->cd(Form("nuisance_%s_mH%d_COMB_%s",  nuisance.c_str(), mH, pngName.c_str() ));
  histoB->Write("histoB",TObject::kOverwrite);
  histoS->Write("histoS",TObject::kOverwrite);
  background->Write("graphB",TObject::kOverwrite);
  signal->Write("graphS",TObject::kOverwrite);
  outFile->Write(); 
  outFile->Close();
  delete outFile; delete histoB; delete histoS; //delete mg; delete background; delete signal;

  gSystem->Exec( "rm ../test/*tmp*" );
  return;
}


  

void produceAll(){

  vector<string> sevenTeV;
  sevenTeV.push_back("7TeV");

  vector<string> eightTeV;
  eightTeV.push_back("8TeV");
  
  vector<string> allTeV;
  allTeV.push_back("7TeV"); allTeV.push_back("8TeV");

  vector<string> mt;
  mt.push_back("mt");

  vector<string> et;
  et.push_back("et");

  vector<string> tt;
  tt.push_back("tt");

  vector<string> em;
  em.push_back("em");

  vector<string> lt;
  lt.push_back("mt"); lt.push_back("et");

  vector<string> lx;
  lx.push_back("mt"); lx.push_back("et");  lx.push_back("em"); 

  vector<string> xt;
  xt.push_back("mt"); xt.push_back("et");  xt.push_back("tt");

  vector<string> xx;
  xx.push_back("mt"); xx.push_back("et");  xx.push_back("tt");  xx.push_back("em"); 

  produce("r",0, 125   ,false, xx, allTeV,    "x+x, 7+8 TeV");
  produce("r",0, 125   ,false, mt, allTeV,    "#mu+#tau, 7+8 TeV");
  produce("r",0, 125   ,false, et, allTeV,    "e+#tau, 7+8 TeV");
  produce("r",0, 125   ,false, tt, allTeV,    "#tau+#tau, 7+8 TeV");
  produce("r",0, 125   ,false, em, allTeV,    "e+#mu, 7+8 TeV");

  produce("CMS_scale_t", 1, 125 ,true, mt, sevenTeV, "#mu+#tau, 7 TeV");
  produce("CMS_scale_t", 1, 125 ,true, mt, eightTeV, "#mu+#tau, 8 TeV");
  produce("CMS_scale_t", 1, 125 ,true, mt, allTeV,   "#mu+#tau, 7+8 TeV");

  produce("CMS_scale_t", 1, 125 ,true, et, sevenTeV, "e+#tau, 7 TeV");
  produce("CMS_scale_t", 1, 125 ,true, et, eightTeV, "e+#tau, 8 TeV");
  produce("CMS_scale_t", 1, 125 ,true, et, allTeV,   "e+#tau, 7+8 TeV");

  produce("CMS_scale_t", 1, 125 ,true, tt, sevenTeV, "#tau+#tau, 7 TeV");
  produce("CMS_scale_t", 1, 125 ,true, tt, eightTeV, "#tau+#tau, 8 TeV");
  produce("CMS_scale_t", 1, 125 ,true, tt, allTeV,   "#tau+#tau, 7+8 TeV");
  
  //produce("CMS_scale_t", 1, 125 ,true, lt, allTeV, "l+#tau, 7+8 TeV");
  //produce("CMS_scale_t", 1, 125 ,true, lx, allTeV, "l+#tau, 7+8 TeV");

  produce("CMS_scale_t", 1, 125 ,true, xt, allTeV,    "x+#tau, 7+8 TeV");

  produce("CMS_eff_t",     1, 125   ,true, xt, allTeV,    "x+#tau, 7+8 TeV");
  produce("CMS_scale_met", 1, 125   ,true, xx, allTeV,    "x+x, 7+8 TeV");
  produce("CMS_scale_jet", 1, 125   ,true, xx, allTeV,    "x+x, 7+8 TeV");
  

}



int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}
