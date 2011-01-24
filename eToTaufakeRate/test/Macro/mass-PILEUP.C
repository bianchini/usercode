

{

  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  //c1->SetLogy(1); 

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.0275);

  TFile Z38("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root"); 
  TFile Z39("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP.root"); 

  TTree* tree38 = (TTree*)Z38.Get("etoTauMargLooseNoCracks80/fitter_tree"); 
  TTree* tree39 = (TTree*)Z39.Get("etoTauMargLooseNoCracks80/fitter_tree"); 

  TH1F* h38 = new TH1F("h38","",25,65,120);
  TH1F* h39 = new TH1F("h39","",25,65,120);
  h38->SetXTitle("mass");
  h38->SetYTitle("fraction");
  h38->SetAxisRange(0.00,1.0,"Y");

  tree38->Draw("mass>>h38","abseta>1.5 && mass>65 && mass<120 && mcTrue && tauAntiEMVA>0.5");
  tree39->Draw("mass>>h39","abseta>1.5 && mass>65 && mass<120 && mcTrue && tauAntiEMVA>0.5");

  //h38->SetBinContent(1,tree38->GetEntries("electronPreIDOutput<-0.1 && mcTrue && abseta>1.5 && mass>65 && mass<120 && mcTrue"));
  //h39->SetBinContent(1,tree39->GetEntries("electronPreIDOutput<-0.1 && mcTrue && abseta>1.5 && mass>65 && mass<120 && mcTrue"));

  h38->SetMarkerStyle(kOpenCircle);
  h39->SetMarkerStyle(kOpenSquare);
  h38->SetMarkerColor(kBlack);
  h39->SetMarkerColor(kRed);

  h38->Sumw2();
  h39->Sumw2();

  h38->DrawNormalized("P");
  h39->DrawNormalized("PSAME");

  leg->AddEntry(h38,"PYTHIA Z#rightarrow e^{+}e^{-}","P");
  leg->AddEntry(h39,"PYTHIA Z#rightarrow e^{+}e^{-} + PU","P");

  leg->Draw();
}
