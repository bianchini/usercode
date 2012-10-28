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
#include "TList.h"
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

#define BLIND 0

void mergeDifferentHistos(TH1F* hTarget=0, TH1F* hInput=0, float weight = 1.0,  bool verbose = false){

  
  for(unsigned int i = 1; i <= (unsigned int)(hTarget->GetNbinsX()); i++){
 
    float centralValueInput = hTarget->GetBinCenter(i);
    float widthInput        = hTarget->GetXaxis()->GetBinUpEdge(i)-hTarget->GetXaxis()->GetBinLowEdge(i);
    float minusHalfInput    = (centralValueInput-widthInput/2.) ;//+ 0.001;
    float plusHalfInput     = (centralValueInput+widthInput/2.) ;//- 0.001;

    if(verbose){
      cout << "Target Bin " << i << endl;
      cout << " ==> centralValueInput = " << centralValueInput << endl;
      cout << " ==> widthInput = " << widthInput << endl; 
      cout << " ==> minusHalfInput = " << minusHalfInput << "; plusHalfInput = " << plusHalfInput << endl; 
    }

    int binLowInput = 0; float widthBinLowInput =0. ; float binLowInputHedgeHigh =0. ; float binLowInputFracHigh=0.;
    if(minusHalfInput<=hInput->GetXaxis()->GetXmax()){
      binLowInput          = hInput->FindBin(minusHalfInput);
      widthBinLowInput     = hInput->GetXaxis()->GetBinUpEdge(binLowInput)-hInput->GetXaxis()->GetBinLowEdge(binLowInput);
      binLowInputHedgeHigh = hInput->GetXaxis()->GetBinUpEdge(binLowInput);
      binLowInputFracHigh  = (binLowInputHedgeHigh-minusHalfInput)/widthBinLowInput;
      if(binLowInputFracHigh<0 || binLowInputFracHigh>1){
	cout << "Low bin fraction is not ok... " << binLowInputFracHigh << endl; 
      }
      else
	if(verbose){
	  cout << "binLowInput = " << binLowInput << ", widthBinLowInput = " << widthBinLowInput << ", binLowInputHedgeHigh = " << binLowInputHedgeHigh <<  endl;
	}
    }
    else
      binLowInput = hInput->GetNbinsX()+1;

    int binHighInput  = 0; float widthBinHighInput=0.; float binHighInputHedgeLow=0.; float binHighInputFracLow=0.;
    if(plusHalfInput<=hInput->GetXaxis()->GetXmax()){
      binHighInput         = hInput->FindBin(plusHalfInput);
      widthBinHighInput    = hInput->GetXaxis()->GetBinUpEdge(binHighInput)-hInput->GetXaxis()->GetBinLowEdge(binHighInput);
      binHighInputHedgeLow = hInput->GetXaxis()->GetBinLowEdge(binHighInput);
      binHighInputFracLow  = (plusHalfInput-binHighInputHedgeLow)/widthBinHighInput;
      if(binHighInputFracLow<0 || binHighInputFracLow>1){
	cout << "Low High fraction is not ok... " << binHighInputFracLow << endl; 
      }
      else
	if(verbose){
	  cout << "binHighInput = " << binHighInput << ", widthBinHighInput = " << widthBinHighInput << ", binHighInputHedgeLow = " << binHighInputHedgeLow <<  endl;
	}
    }
    else
      binHighInput = hInput->GetNbinsX();
      
    float entryi = 0.; float errori2 = 0.;

    for(int j = binLowInput; j <= binHighInput; j++){

      float widthj = hInput->GetBinWidth(j);
      float inputj = hInput->GetBinContent(j);
      float errorj = hInput->GetBinError(j)*widthj;
      float entryj = inputj*widthj;

      if(j==binLowInput){
	entryi  += binLowInputFracHigh*entryj;
	errori2 += (binLowInputFracHigh*errorj)*(binLowInputFracHigh*errorj);
	if(verbose) cout << entryi << endl;
      }
      else if(j==binHighInput){
	entryi  += binHighInputFracLow*entryj;
	errori2 += (binHighInputFracLow*errorj)*(binHighInputFracLow*errorj);
	if(verbose) cout << "+" << entryi << endl;
      }
      else{
	entryi  += entryj;
	errori2 += (errorj*errorj);
	if(verbose) cout << "+" << entryj << endl;
      }
    }

    float oldContenti = hTarget->GetBinContent(i);
    float oldErrori   = hTarget->GetBinError(i);
    oldContenti *= widthInput;
    oldErrori   *= widthInput;
    float newContenti = oldContenti + entryi*weight;
    float newErrori   = TMath::Sqrt(oldErrori*oldErrori + errori2*weight*weight);
    hTarget->SetBinContent(i, newContenti/widthInput);
    hTarget->SetBinError(i,   newErrori/widthInput);

  }

  if(verbose)  cout << "#### Initial ==> " << hInput->Integral("width") << endl;
  if(verbose)  cout << "#### Final   ==> " << hTarget->Integral("width") << endl;

}

void blindHisto(TH1F* histo=0){

  int vetoBinLow  = histo->FindBin(90.);
  int vetoBinHigh = histo->FindBin(130.);

  for(int i =0; i <= histo->GetNbinsX() ; i++){
    if( i>= vetoBinLow && i<= vetoBinHigh){
      histo->SetBinContent(i, 0.0);
      histo->SetBinError(i, 0.0);
    }
  }
}


void produce(int mH                     = 120,
	     float signalRescale        = 0.20, 
	     vector<string>& channels   = *(new  vector<string>()),
	     vector<string>& com        = *(new  vector<string>()),
	     Int_t nBins_ = 80, 
	     Float_t xMin_=  0, Float_t xMax_=400,
	     bool verbose               = false
	     ){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(0);

  TPad* pad1 = new TPad("pad1","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2","",0.05,0.02,0.96,0.20);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  TCanvas *c1W = new TCanvas("c1W","",5,30,650,600);
  c1W->SetGrid(0,0);
  c1W->SetFillStyle(4000);
  c1W->SetFillColor(10);
  c1W->SetTicky();
  c1W->SetObjectStat(0);
  c1W->SetLogy(0);

  TPad* pad1W = new TPad("pad1W","",0.05,0.22,0.96,0.97);
  TPad* pad2W = new TPad("pad2W","",0.05,0.02,0.96,0.20);
  pad1W->SetFillColor(0);
  pad2W->SetFillColor(0);
  pad1W->Draw();
  pad2W->Draw();

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

  TLegend* leg = new TLegend(0.50,0.50,0.75,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);

  string title = "";

  typedef std::map<string, std::pair<vector<string>,vector<string> > >::iterator IT;
  std::map<string, std::pair<vector<string>,vector<string> > > channelToCatMap; 

  vector<string> categoriesLT; vector<string> categoriesIndexLT;
  vector<string> categoriesTT; vector<string> categoriesIndexTT;
  vector<string> categoriesMM; vector<string> categoriesIndexMM;

  categoriesLT.push_back("0jet_low");   categoriesIndexLT.push_back("0");
  categoriesLT.push_back("0jet_high");  categoriesIndexLT.push_back("1");
  categoriesLT.push_back("boost_low");  categoriesIndexLT.push_back("2");
  categoriesLT.push_back("boost_high"); categoriesIndexLT.push_back("3");
  categoriesLT.push_back("vbf");        categoriesIndexLT.push_back("5");

  categoriesTT.push_back("boost");      categoriesIndexTT.push_back("0");
  categoriesTT.push_back("vbf");        categoriesIndexTT.push_back("1");

  categoriesMM.push_back("dummy");      categoriesIndexTT.push_back("0");

  vector<string> channelsDir;
  for(unsigned int ch = 0; ch<channels.size(); ch++){
    cout << channels[ch] << endl;
    if(channels[ch].find("mt")!=string::npos){
      title = title+"#tau_{#mu}#tau_{h} ";
      channelToCatMap.insert( make_pair("mt", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("muTau");
    }
    if(channels[ch].find("et")!=string::npos){
      title = title+"#tau_{e}#tau_{h} ";
      channelToCatMap.insert( make_pair("et", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("eleTau");
    }
    if(channels[ch].find("tt")!=string::npos){
      title = title+"#tau_{h}#tau_{h} ";
      channelToCatMap.insert( make_pair("tt", make_pair(categoriesTT, categoriesIndexTT) )  );
      channelsDir.push_back("tauTau");
    }
    if(channels[ch].find("em")!=string::npos){
      title = title+"#tau_{e}#tau_{#mu} ";
      channelToCatMap.insert( make_pair("em", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("emu");
    }
    if(channels[ch].find("mm")!=string::npos){
      title = title+"#tau_{#mu}#tau_{#mu} ";
      channelToCatMap.insert( make_pair("mm", make_pair(categoriesMM, categoriesIndexMM) )  );
      channelsDir.push_back("muMu");
    }
  }

  for(unsigned int en = 0; en<com.size(); en++){
    cout << com[en] << endl;
    if(en!= (com.size()-1) ) 
      title = (title + ", " +  com[en] + "+");
    else
      title = (title + com[en]);
  }

  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open("../test/bins_coarse.txt"); 
  if(nBins_<0 &&  !is.good()){
    cout << "Bins file not found" << endl;
    return;
  }

  int nBinsFromFile = 0;
  while (is.good())     
    {
      is.getline(c,999,',');     
      if (is.good()){
	nBinsFromFile++;
      }
    }

  int nBins =  nBins_>0 ? nBins_ : nBinsFromFile-1 ;
  TArrayF bins(nBins+1);
  cout << "Making histograms with " << nBins << " bins:" << endl;

  is.close();
  is.open("../test/bins_coarse.txt");
  
  nBinsFromFile = 0;

  if(nBins_>0){
    for( ; nBinsFromFile <= nBins ; nBinsFromFile++){
      bins[nBinsFromFile] =  xMin_ + nBinsFromFile*(xMax_-xMin_)/nBins_;
    }
  }
  else{
    while (is.good())  
      {
	is.getline(c,999,',');     
	if (is.good() && nBinsFromFile<=nBins) {
	  bins[nBinsFromFile] = atof(c);
	  cout << bins[nBinsFromFile] << ", " ;
	}
	nBinsFromFile++;
      }
    cout << endl;
  }

  THStack* aStack  = new THStack("aStack","");
  THStack* aStackW = new THStack("aStackW","");

  TH1F* hData  = new TH1F( "hData" ,"Data; mass (GeV); events/ 1 GeV"                             , nBins , bins.GetArray());
  TH1F* hDataW = new TH1F( "hDataW","Data weighted; mass (GeV); S/(S+B) weighted events/ 1 GeV"   , nBins , bins.GetArray());
  TH1F* hBkg   = new TH1F( "hBkg"   ,"MC"             , nBins , bins.GetArray());
  TH1F* hBkgW  = new TH1F( "hBkgW"  ,"MC weighted"    , nBins , bins.GetArray());
  TH1F* hSgn   = new TH1F( "hSgn"  ,"Signal"          , nBins , bins.GetArray());
  TH1F* hSgnW  = new TH1F( "hSgnW" ,"Signal weightes" , nBins , bins.GetArray());
  TH1F* hErr   = new TH1F( "hErr"  ,"Signal"          , nBins , bins.GetArray());
  TH1F* hErrW  = new TH1F( "hErrW" ,"Signal weightes" , nBins , bins.GetArray());

  hBkg->SetFillColor(kBlue);
  hBkg->SetLineColor(kBlue);
  hBkg->SetLineWidth(1);
  hBkg->SetFillStyle(3003);
  hBkgW->SetFillColor(kBlue);
  hBkgW->SetFillColor(kBlue);
  hBkg->SetLineWidth(1);
  hBkgW->SetFillStyle(3003);

  hData->Sumw2();
  hData->SetLabelSize(0.04,"X");
  hData->SetTitleSize(0.06,"X");
  hData->SetLabelSize(0.04,"Y");
  hData->SetTitleSize(0.06,"Y");
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.6);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hDataW->Sumw2();
  hData->SetLabelSize(0.04,"X");
  hData->SetTitleSize(0.06,"X");
  hData->SetLabelSize(0.04,"Y");
  hData->SetTitleSize(0.06,"Y");
  hDataW->SetMarkerStyle(20);
  hDataW->SetMarkerSize(0.6);
  hDataW->SetMarkerColor(kBlack);
  hDataW->SetLineColor(kBlack);

  hSgn->SetFillColor(kRed);
  hSgn->SetLineColor(kRed);
  hSgn->SetFillStyle(1001);
  hSgn->SetLineWidth(1);
  hSgn->SetLineStyle(kSolid);
  hSgnW->SetFillColor(kRed);
  hSgnW->SetFillStyle(1001);
  hSgnW->SetLineColor(kRed);
  hSgnW->SetLineWidth(1);
  hSgnW->SetLineStyle(kSolid);

  hErr->SetMarkerSize(0);
  hErr->SetFillColor(1);
  hErr->SetFillStyle(3013);
  hErrW->SetMarkerSize(0);
  hErrW->SetFillColor(1);
  hErrW->SetFillStyle(3013);

  leg->AddEntry(hSgn,Form("SM Higgs (m_{H}=%d GeV)",mH),"F");
  leg->AddEntry(hData,"Observed","PE");
  leg->AddEntry(hBkg,"Background","F");
  leg->AddEntry(hErr,"#pm 1#sigma","F");
  leg->SetHeader(title.c_str());

  TList *listData = new TList();
  TList *listMC   = new TList();
  TList *listSgn  = new TList();
  TList *listDataW = new TList();
  TList *listMCW   = new TList();
  TList *listSgnW  = new TList();

  for(unsigned int en = 0; en<com.size(); en++){
    string iEn = com[en];

    for(unsigned int ch = 0; ch<channels.size(); ch++){

      string iCh   = channels[ch];
      string iChIn = channelsDir[ch];

      vector<string> categories      = channelToCatMap[iCh].first;
      vector<string> categoriesIndex = channelToCatMap[iCh].second;
      
      for(unsigned int cat = 0; cat< categories.size(); cat++){

	string iCat   = categories[cat];
	string iCatIn = categoriesIndex[cat];

	cout << "Now doing " << string(Form("%s_%s_rescaled_%s_.root", iChIn.c_str(),iCat.c_str(),iEn.c_str() )) << endl;
	
	TFile* fInput = new TFile(Form("../test/%s_%s_rescaled_%s_.root", iChIn.c_str(),iCat.c_str(),iEn.c_str() ),"READ");
	if(!fInput || fInput->IsZombie() ) continue;
	
	TH1F* ihSgn = 0;  TH1F* ihBkg  = 0; TH1F* ihData  = 0; TH1F* ihErr  = 0;
	TH1F* ihSgnW = 0; TH1F* ihBkgW = 0; TH1F* ihDataW = 0; TH1F* ihErrW = 0;

	TIter nextkey(gDirectory->GetListOfKeys());
	TH1F *key;
	while ( (key = (TH1F*)nextkey()) ) {
	  string name(key->GetName());
	  TH1F* h = (TH1F*)gDirectory->FindObjectAny( key->GetName());

	  if(verbose) cout << "Processing " << name << ": " << h->Integral("width") << endl;

	  if(  name.find("ggH")!=string::npos ){
	    if(verbose) cout << "... This is signal ... " << name << endl;
	    if( !ihSgn ){
	      ihSgn  = (TH1F*)h->Clone("ihSgn");
	      ihSgnW = (TH1F*)h->Clone("ihSgnW");
	      //ihSgn->Add(h);
	      //ihSgnW->Add(h);
	    }
	    else{
	      ihSgn->Add(h);
	      ihSgnW->Add(h);
	    }
	  }
	  else if(  name.find("data_obs")==string::npos &&  name.find("Ztt")!=string::npos) {
	    if(verbose) cout << "... This is background ... " << name << endl;
	    if( !ihBkg ){
	      ihBkg  = (TH1F*)h->Clone("ihBkg");
	      ihBkgW = (TH1F*)h->Clone("ihBkgW");
	      //ihBkg->Add(h);
	      //ihBkgW->Add(h);
	    }
	    else{
	      ihBkg->Add(h);
	      ihBkgW->Add(h);
	    }
	  }
	  else if( name.find("data_obs")!=string::npos){
	    if(verbose) cout << "... This is data ... " << name << endl;
	    if( !ihData ){
	      ihData  = (TH1F*)h->Clone("ihData");
	      ihDataW = (TH1F*)h->Clone("ihDataW");
	      //ihData->Add(h);
	      //ihDataW->Add(h);
	    }
	    else{
	      ihData->Add(h);
	      ihDataW->Add(h);
	    }
	  }
	  else if( name.find("errorBand")!=string::npos ){
	    if(verbose) cout << "... This is error band ... " << name << endl;
	    if( !ihErr ){
	      ihErr  = (TH1F*)h->Clone("ihErr");
	      ihErrW = (TH1F*)h->Clone("ihErrW");
	      //ihErr->Add(h);
	      //ihErrW->Add(h);
	    }
	    else{
	      ihErr->Add(h);
	      ihErrW->Add(h);
	    }
	  }
	  else{
	    if(verbose) cout << "This histo will not be considered!" << endl;
	  }
	}

	if( !ihSgn || !ihBkg || !ihData || !ihErr){
	  cout << "Could not find some histos... continue" << endl;
	  continue;
	}

	if(ihSgn->Integral() > ihBkg->Integral()){
	  ihSgn->Add(ihBkg, -1.0);
	  ihSgnW->Add(ihBkg, -1.0);
	}
	ihSgn->Scale(signalRescale);
	ihSgnW->Scale(signalRescale);


	double quantiles[] = {0.05,0.16, 0.50, 0.84, 0.95};
	double positions[] = {0.0,0.0, 0.0, 0.0, 0.0};
	ihSgn->GetQuantiles(5, positions, quantiles);
	int binLow  = ihSgn->FindBin(positions[1]);
	int binHigh = ihSgn->FindBin(positions[3]);

	//binLow  = 1;
	//binHigh = ihSgn->GetNbinsX();

	float totalSig = ihSgn->Integral(binLow,binHigh,"width");
	float totalBkg = ihBkg->Integral(binLow,binHigh,"width");

	float errorBkg2 = 0.;
	for(int i = binLow; i <= binHigh; i++){
	  errorBkg2 += (ihErr->GetBinError(i)*ihErr->GetBinWidth(i))*(ihErr->GetBinError(i)*ihErr->GetBinWidth(i));
	}

	float SoB       = totalSig/(totalBkg+totalSig);
	float SoSqrtB   = totalSig/TMath::Sqrt(totalBkg);
	float SoSqrtBdB = totalSig/TMath::Sqrt(totalBkg+errorBkg2);
	float SoBpS     = totalSig/(totalBkg+totalSig);

	cout << "Weight S/B          = " << totalSig << "/" << totalBkg+totalSig << " = "  << SoB << endl;
	cout << "Weight S/sqrt(B)    = " << totalSig << "/sqrt(" << totalBkg << ") = " << SoSqrtB << endl;
	cout << "Weight S/sqrt(B+dB) = " << totalSig << "/sqrt(" << totalBkg << "+" << errorBkg2 << ") = " << SoSqrtBdB << endl;
	//cout << "Weight S/(S+B) = "   << SoBpS << endl;

	mergeDifferentHistos(hSgn,  ihSgn,   1.0, verbose);
	mergeDifferentHistos(hSgnW, ihSgnW , SoB, verbose);
	mergeDifferentHistos(hBkg,  ihBkg,   1.0, verbose);
	mergeDifferentHistos(hBkgW, ihBkgW , SoB, verbose);
	mergeDifferentHistos(hData, ihData,  1.0, verbose);
	mergeDifferentHistos(hDataW,ihDataW, SoB, verbose);
	mergeDifferentHistos(hErr, ihErr,    1.0, verbose);
	mergeDifferentHistos(hErrW,ihErrW,   SoB, verbose);

	listData->Add(ihData);
	listDataW->Add(ihDataW);
	listMC->Add(ihBkg);
	listMCW->Add(ihBkgW);
	listSgn->Add(ihSgn);
	listSgnW->Add(ihSgnW);

      }
      
    }
     
  }

  float maxData   = hData->GetBinContent(hData->GetMaximumBin());
  float maxDataW  = hDataW->GetBinContent(hDataW->GetMaximumBin());;

  float maxRatio = 0.;
  for(unsigned int i = 1; i <= (unsigned int)(hData->GetNbinsX()); i++){
    if( TMath::Abs(hData->GetBinContent(i) - hBkg->GetBinContent(i) ) > maxRatio ) maxRatio = TMath::Abs(hData->GetBinContent(i) - hBkg->GetBinContent(i) );
  }
  cout << maxRatio << endl;
  float maxRatioW = 0.;
  for(unsigned int i = 1; i <= (unsigned int)(hDataW->GetNbinsX()); i++){
    if( TMath::Abs(hDataW->GetBinContent(i) - hBkgW->GetBinContent(i) ) > maxRatioW ) maxRatioW = TMath::Abs(hDataW->GetBinContent(i) - hBkgW->GetBinContent(i) );
  }
  cout << maxRatioW << endl;


  if(BLIND){
    blindHisto(hData);
    blindHisto(hDataW);
  }

  TH1F* hErrSub = (TH1F*)hErr->Clone("hErrSub");
  hErrSub->Add(hBkg,-1.0);

  TH1F* hErrSubW = (TH1F*)hErrW->Clone("hErrSubW");
  hErrSubW->Add(hBkgW,-1.0);

  aStack->Add(hBkg);
  aStack->Add(hSgn);

  aStackW->Add(hBkgW);
  aStackW->Add(hSgnW);

  cout << "Unweighted = " << endl;
  cout << " Sgn TOT = "   << hSgn->Integral("width") << endl;
  cout << " Bkg TOT = "   << hBkg->Integral("width") << endl;
  cout << " Data TOT = "  << hData->Integral("width") << endl;
  cout << "Weighted = "   << endl;
  cout << " SgnW TOT = "  << hSgnW->Integral("width") << endl;
  cout << " BkgW TOT = "  << hBkgW->Integral("width") << endl;
  cout << " DataW TOT = " << hDataW->Integral("width") << endl;

  c1->cd();
  pad1->cd();

  hData->SetAxisRange(0,1.2*maxData,"Y");
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hErr->Draw("e2SAME");
  hData->Draw("PSAME");
  leg->Draw();

  pad2->cd();

  TH1F* hRatio = (TH1F*)hData->Clone("hRatio");
  hRatio->Reset(); 
  hRatio->Add(hData,1.0);
  hRatio->SetTitle("");
  hRatio->SetYTitle("");
  hRatio->SetXTitle("");
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.6);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetLabelSize(0.16,"Y");
  hRatio->SetLabelSize(0.16,"X");
  hRatio->Add(hBkg,-1.0);
  if(BLIND) blindHisto(hRatio);
  for(unsigned int i = 1; i <= (unsigned int)(hRatio->GetNbinsX()); i++)
    hRatio->SetBinError(i, hData->GetBinError(i));
    

  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);

  hRatio->SetAxisRange(-1.5*maxRatio,1.5*maxRatio,"Y");
  hRatio->Draw("PE");
  hSgn->Draw("HISTSAME");
  hRatio->Draw("PESAME");
  hErrSub->Draw("e2SAME");
  line->Draw("SAME");

  c1->SaveAs("unweighted.png");


  c1W->cd();

  pad1W->cd();
  hDataW->Draw("P");
  hDataW->SetAxisRange(0,1.2*maxDataW,"Y");
  aStackW->Draw("HISTSAME");
  hErrW->Draw("e2SAME");
  hDataW->Draw("PSAME");
  leg->Draw();

  pad2W->cd();

  TH1F* hRatioW = (TH1F*)hDataW->Clone("hRatioW");
  hRatioW->SetTitle("");
  hRatioW->SetYTitle("");
  hRatioW->SetXTitle("");
  hRatioW->SetMarkerStyle(20);
  hRatioW->SetMarkerSize(0.6);
  hRatioW->SetMarkerColor(kBlack);
  hRatioW->SetLineColor(kBlack);
  hRatioW->SetLabelSize(0.16,"Y");
  hRatioW->SetLabelSize(0.16,"X");
  hRatioW->Add(hBkgW,-1.0);
  hRatioW->SetAxisRange(-1.2*maxRatioW,1.2*maxRatioW,"Y");
  if(BLIND) blindHisto(hRatioW);
  for(unsigned int i = 1; i <= (unsigned int)(hRatioW->GetNbinsX()); i++)
    hRatioW->SetBinError(i, hDataW->GetBinError(i));

  hRatioW->SetAxisRange(-1.5*maxRatioW,1.5*maxRatioW,"Y");
  hRatioW->Draw("PE");
  hSgnW->Draw("HISTSAME");
  hRatioW->Draw("PESAME");
  hErrSubW->Draw("e2SAME");
  line->Draw("SAME");

  c1W->SaveAs("weighted.png");

  delete hData; delete hDataW; delete hBkg; delete hBkgW; delete hSgn; delete hSgnW;

  return;
}


  

void produceAll( vector<string>& channels = *(new vector<string>) ,
		 vector<string>& com      = *(new vector<string>),
		 string title = ""
		 ){

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
  lx.push_back("em"); 
  lx.push_back("mt"); 
  lx.push_back("et");  

  vector<string> xt;
  xt.push_back("mt"); xt.push_back("et");  xt.push_back("tt");

  vector<string> xx;
  xx.push_back("mt"); xx.push_back("et");  xx.push_back("tt");  xx.push_back("em"); 

  produce(125, 0.20, channels , com , -1, 0, 400);

  return;

}



int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  if ( argc < 3 ) {
    std::cout << "Usage: " << argv[0] << " <channel> <period>" << std::endl;
    return 0;
  }

  vector<string> channels; vector<string> com;

  for(int i = 1; i < argc; i++){
    string input(argv[i]);
    if( input.find("mt")!=string::npos ||
	input.find("et")!=string::npos ||
	input.find("em")!=string::npos ||
	input.find("tt")!=string::npos ||
	input.find("mm")!=string::npos
	) 
      channels.push_back(input);
    else if( input.find("TeV")!=string::npos )
      com.push_back(input);
    else
      cout << input << ": invalid input" << endl;
  }

  produceAll( channels, com);

}
