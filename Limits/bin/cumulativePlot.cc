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


void produce(int mH                   = 120, 
	     vector<string>& channels = *(new  vector<string>()),
	     vector<string>& com      = *(new  vector<string>()),
	     Int_t nBins_ = 80, 
	     Float_t xMin_=  0, Float_t xMax_=400,
	     TString title            = ""
	     ){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(0);

  TCanvas *c1W = new TCanvas("c1W","",5,30,650,600);
  c1W->SetGrid(0,0);
  c1W->SetFillStyle(4000);
  c1W->SetFillColor(10);
  c1W->SetTicky();
  c1W->SetObjectStat(0);
  c1W->SetLogy(0);

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

  for(unsigned int en = 0; en<com.size(); en++){
    cout << com[en] << endl;
  }

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
      channelToCatMap.insert( make_pair("mt", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("muTau");
    }
    if(channels[ch].find("et")!=string::npos){
      channelToCatMap.insert( make_pair("et", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("eleTau");
    }
    if(channels[ch].find("tt")!=string::npos){
      channelToCatMap.insert( make_pair("tt", make_pair(categoriesTT, categoriesIndexTT) )  );
      channelsDir.push_back("tauTau");
    }
    if(channels[ch].find("em")!=string::npos){
      channelToCatMap.insert( make_pair("em", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("emu");
    }
    if(channels[ch].find("mm")!=string::npos){
      channelToCatMap.insert( make_pair("mm", make_pair(categoriesMM, categoriesIndexMM) )  );
      channelsDir.push_back("muMu");
    }
  }


  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open("../test/bins.txt"); 
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
  is.open("../test/bins.txt");
  
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


  TH1F* hData  = new TH1F( "hData" ,"Data; mass (GeV); events/GeV"            , nBins , bins.GetArray());
  TH1F* hDataW = new TH1F( "hDataW","Data weighted; mass (GeV); weighted events/GeV"   , nBins , bins.GetArray());
  TH1F* hMC    = new TH1F( "hMC"   ,"MC"              , nBins , bins.GetArray());
  TH1F* hMCW   = new TH1F( "hMCW"  ,"MC weighted"     , nBins , bins.GetArray());
  TH1F* hSgn   = new TH1F( "hSgn"  ,"Signal"          , nBins , bins.GetArray());
  TH1F* hSgnW  = new TH1F( "hSgnW" ,"Signal weightes" , nBins , bins.GetArray());

  hMC->SetFillColor(kOrange-4);
  hMC->SetFillStyle(3003);
  hMCW->SetFillColor(kOrange-4);
  hMCW->SetFillStyle(3003);

  hData->Sumw2();
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(1.2);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hDataW->Sumw2();
  hDataW->SetMarkerStyle(20);
  hDataW->SetMarkerSize(1.2);
  hDataW->SetMarkerColor(kBlack);
  hDataW->SetLineColor(kBlack);
  hSgn->SetFillColor(0);
  hSgn->SetLineColor(kBlue);
  hSgn->SetLineWidth(2);
  hSgn->SetLineStyle(kDashed);
  hSgnW->SetFillColor(0);
  hSgnW->SetLineColor(kBlue);
  hSgnW->SetLineWidth(2);
  hSgnW->SetLineStyle(kDashed);

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
	
	TFile* fInput = new TFile(Form("../test/%s_%s_rescaled_%s_.root", iChIn.c_str(),iCat.c_str(),iEn.c_str() ),"READ");
	if(!fInput || fInput->IsZombie() ) continue;
	
	TH1F* ihSgn = 0;  TH1F* ihBkg  = 0; TH1F* ihData  = 0;
	TH1F* ihSgnW = 0; TH1F* ihBkgW = 0; TH1F* ihDataW = 0;

	TIter nextkey(gDirectory->GetListOfKeys());
	TH1F *key;
	while ( (key = (TH1F*)nextkey()) ) {
	  string name(key->GetName());
	  TH1F* h = (TH1F*)gDirectory->FindObjectAny( key->GetName());

	  if(  name.find("ggH")!=string::npos ){
	    cout << "This is signal ... " << name << endl;
	    if( !ihSgn ){
	      ihSgn  = new TH1F(Form("ihSgn_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihSgn");
	      ihSgnW = new TH1F(Form("ihSgnW_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihSgn");
	      ihSgn->Add(h);
	      ihSgnW->Add(h);
	    }
	    else{
	      ihSgn->Add(h);
	      ihSgnW->Add(h);
	    }
	  }
	  else if(  name.find("data_obs")==string::npos &&  name.find("Ztt")!=string::npos) {
	    cout << "This is background ... " << name << endl;
	    if( !ihBkg ){
	      ihBkg  = new TH1F(Form("ihBkg_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihBkg");
	      ihBkgW = new TH1F(Form("ihBkgW_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihBkg");
	      ihBkg->Add(h);
	      ihBkgW->Add(h);
	    }
	    else{
	      ihBkg->Add(h);
	      ihBkgW->Add(h);
	    }
	  }
	  else if( name.find("data_obs")!=string::npos ){
	    cout << "This is data ... " << name << endl;
	    if( !ihData ){
	      ihData  = new TH1F(Form("ihData_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihData");
	      ihDataW = new TH1F(Form("ihDataW_%s_%s_%s", iChIn.c_str(),iCat.c_str(),iEn.c_str()), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()) ;//(TH1F*)h->Clone("ihData");
	      ihData->Add(h);
	      ihDataW->Add(h);
	    }
	    else{
	      ihData->Add(h);
	      ihDataW->Add(h);
	    }
	  }
	  else
	    cout << "This histo will not be considered!" << endl;
	}

	if( !ihSgn || !ihBkg || !ihData){
	  cout << "Could not find some histos... continue" << endl;
	  continue;
	}

	if(ihSgn->Integral() > ihBkg->Integral()){
	  ihSgn->Add(ihBkg, -1.0);
	  ihSgnW->Add(ihBkg, -1.0);
	}
	ihSgn->Scale(0.25);
	ihSgnW->Scale(0.25);

	cout << "Sgn  = " << ihSgn->Integral() << endl;
	cout << "Bkg  = " << ihBkg->Integral() << endl;
	cout << "Data = " << ihData->Integral() << endl;

	double quantiles[] = {0.05,0.16, 0.50, 0.84, 0.95};
	double positions[] = {0.0,0.0, 0.0, 0.0, 0.0};
	ihSgn->GetQuantiles(5, positions, quantiles);
	int binLow  = ihSgn->FindBin(positions[0]);
	int binHigh = ihSgn->FindBin(positions[4]);

	float totalSig = ihSgn->Integral(binLow,binHigh);
	float totalBkg = ihBkg->Integral(binLow,binHigh);

	float SoB     = totalSig/totalBkg;
	float SoSqrtB = totalSig/TMath::Sqrt(totalBkg);
	float SoBpS   = totalSig/(totalBkg+totalSig);

	cout << "Weight S/B = "       << SoB << endl;
	cout << "Weight S/sqrt(B) = " << SoSqrtB << endl;
	cout << "Weight S/(S+B) = "   << SoBpS << endl;

	ihDataW->Scale(SoB);
	ihSgnW->Scale(SoB);
	ihBkgW->Scale(SoB);

	listData->Add(ihData);
	listDataW->Add(ihDataW);

	listMC->Add(ihBkg);
	listMCW->Add(ihBkgW);

	listSgn->Add(ihSgn);
	listSgnW->Add(ihSgnW);

	//fInput->Close();
	//delete fInput;

      }
      
    }
     
  }
    

  hSgn->Merge(listSgn);
  hSgnW->Merge(listSgnW);
  hData->Merge(listData);
  hDataW->Merge(listDataW);
  hMC->Merge(listMC);
  hMCW->Merge(listMCW);

  cout << "Unweighted = " << endl;
  cout << "Sgn TOT = " << hSgn->Integral() << endl;
  cout << "Bkg TOT = " << hMC->Integral() << endl;
  cout << "Data TOT = " << hData->Integral() << endl;
  cout << "Weighted = " << endl;
  cout << "SgnW TOT = " << hSgnW->Integral() << endl;
  cout << "BkgW TOT = " << hMCW->Integral() << endl;
  cout << "DataW TOT = " << hDataW->Integral() << endl;

  c1->cd();
  hData->Draw("P");
  hMC->Draw("HISTSAME");
  hSgn->Draw("HISTSAME");
  hData->Draw("PSAME");
  c1->SaveAs("unweighted.png");

  c1W->cd();
  hDataW->Draw("P");
  hMCW->Draw("HISTSAME");
  hSgnW->Draw("HISTSAME");
  hDataW->Draw("PSAME");
  c1W->SaveAs("weighted.png");

  delete hData; delete hDataW; delete hMC; delete hMCW; delete hSgn; delete hSgnW;

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

  produce(125, mt, sevenTeV, -1, 0, 400, "");

  return;

}



int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}



// 	for(unsigned int i = 1; i <= (unsigned int)(hData->GetNbinsX()); i++){
// 	  float centralValue = hData->GetBinCenter(i);
// 	  //cout <<  " ==> " << centralValue << endl;

// 	  float totData = ihData->GetBinContent(ihData->FindBin(centralValue));
// 	  float totSgn  = ihSgn->GetBinContent(ihSgn->FindBin(centralValue));
// 	  float totBkg  = ihBkg->GetBinContent(ihBkg->FindBin(centralValue));

// 	  float binWidth= ihData->GetBinWidth(ihData->FindBin(centralValue));
// 	  //cout << totData << endl;cout << totSgn << endl;cout << totBkg << endl;

// 	  float oldDataContent = hData->GetBinContent(i);
// 	  float oldDataError   = hData->GetBinError(i);
// 	  hData->SetBinContent(i,  totData+oldDataContent);
// 	  hData->SetBinError(i, TMath::Sqrt( oldDataError*oldDataError + totData*totData) ); // to be checked!
// 	  hDataW->SetBinContent(i,  SoB*totData+oldDataContent);
// 	  hDataW->SetBinError(i, TMath::Sqrt( oldDataError*oldDataError + SoB*SoB*totData*totData) ); // to be checked!

// 	  float oldSgnContent = hSgn->GetBinContent(i);
// 	  float oldSgnError   = hSgn->GetBinError(i);
// 	  hSgn->SetBinContent(i,  totSgn+oldSgnContent);
// 	  hSgn->SetBinError(i, TMath::Sqrt( oldSgnError*oldSgnError + totSgn*totSgn) ); // to be checked!
// 	  hSgnW->SetBinContent(i,  SoB*totSgn+oldSgnContent);
// 	  hSgnW->SetBinError(i, TMath::Sqrt( oldSgnError*oldSgnError + SoB*SoB*totSgn*totSgn) ); // to be checked!
	  
// 	  float oldBkgContent = hMC->GetBinContent(i);
// 	  float oldBkgError   = hMC->GetBinError(i);
// 	  hMC->SetBinContent(i,  totBkg+oldBkgContent);
// 	  hMC->SetBinError(i, TMath::Sqrt( oldBkgError*oldBkgError + totBkg*totBkg) ); // to be checked!
// 	  hMCW->SetBinContent(i,  SoB*totBkg+oldBkgContent);
// 	  hMCW->SetBinError(i, TMath::Sqrt( oldBkgError*oldBkgError + SoB*SoB*totBkg*totBkg) ); // to be checked!
// 	}
