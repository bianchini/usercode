
 
{
  TFile fbkg("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  
  TTree *fullTreeBkg   = (TTree*)fbkg.Get("etoTauSCMargNoCracks70/fitter_tree");

  TH1F* h = new TH1F("h","",12,60,120);
 
  h->SetXTitle("mass");
  h->SetYTitle("entries/(5 GeV)");

  fullTreeBkg->Draw("mass>>h","tauAntiEMVA<0.5 && abseta>1.5");

  h->Sumw2();

  TF1* expo = new TF1("expo","expo",60,120);

  TFitResultPtr r = h->Fit(expo,"S",0,65,120);

  //r->Print("V"); 
  std::cout << expo->GetChisquare()/ (h->GetNbinsX()-3) << std::endl;

  h->Draw("P");

  TLegend* leg = new TLegend(0.6,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->SetHeader("Shrinking Cone, |#eta|<1.5");
  leg->AddEntry(h,  "MC: passing probes","L");
  leg->AddEntry(expo,  Form("exponential: #chi^{2}/ndof=%.2f", expo->GetChisquare()/ (h->GetNbinsX()-3) ),"L");

  leg->Draw();
}
