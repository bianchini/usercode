#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include "TAttMarker.h"
#include "TAxis.h"


void plot(std::string histo_ = "PatZeeAnalyzer/Zee/hDiEleMass",
	  std::string title_="di-electrons mass", 
	  std::string Xlabel_ = "mass [GeV/c^{2}]",
	  int rebin_ = 1, int firstbin_ = 0, int lastbin_ = 100, 
	  bool verbose_= false ){

  TFile *fRD = new TFile("./Data10May/PFAnalysis_DataSDGC.root","READ");
  TFile *fMC = new TFile("./MC/PFAnalysis_MC.root","READ");

  //"PatZeeAnalyzer/Zee/"+histo_;

  TH1F* hRD = (TH1F*) fRD->Get(histo_.c_str());
  TH1F* hMC = (TH1F*) fMC->Get(histo_.c_str());
  hMC->SetName("hMC");
  hRD->SetName("hRD");

  hRD->SetMarkerStyle(20);
  hRD->SetMarkerColor(1);
  hMC->SetLineColor(2);
  hMC->SetFillColor(2);
  hMC->SetFillStyle(3004);

  hRD->SetTitle(title_.c_str());
  hMC->SetTitle(title_.c_str());

  int eRD = hRD->GetEntries();
  int eMC = hMC->GetEntries();
  if(verbose_) std::cout << "Data= " << eRD << "  MC= " <<  eMC << std::endl;

  hMC->Scale( double(eRD)/double(eMC) );

  if(verbose_) std::cout
    << "Integral Data= "
    << (double)hRD->Integral() 
    << "  Integral MC= " 
    << (double)hMC->Integral() 
    << std::endl;


  TH1F* hRDnew =  (TH1F*) hRD->Rebin(rebin_,"hRDnew");
  TH1F* hMCnew =  (TH1F*) hMC->Rebin(rebin_,"hMCnew");

  TAxis* RDax = hRDnew->GetXaxis();
  TAxis* RDaxY = hRDnew->GetYaxis();
  RDax->SetTitle(Xlabel_.c_str()); 
  TAxis* MCax = hMCnew->GetXaxis();
  MCax->SetTitle(Xlabel_.c_str());  

  RDax->SetRangeUser(firstbin_,lastbin_);
  MCax->SetRangeUser(firstbin_,lastbin_);

  double RDmax = hRDnew->GetMaximum();
  double MCmax = hMCnew->GetMaximum();
  double max = std::max<double>(RDmax,MCmax);

  RDaxY->SetRangeUser(0.0,(max+TMath::Sqrt(max))*1.1);

  hRDnew->Draw("P");
  hMCnew->Draw("HISTSAME");

  if(verbose_) hMC->Print("ALL");

  TLegend *leg = new TLegend(0.4719235,0.5741525,0.8088411,0.8559322,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05);

  leg->AddEntry(hRDnew,"7 TeV");
  leg->AddEntry(hMCnew,"Minimum Bias");

  leg->SetHeader("#splitline{OS, (ch+ph^{pt>0.5})[#DeltaR<0.3]<3 GeV}{mva > -0.3, p_{T}>3 GeV}");

  leg->Draw();


}
