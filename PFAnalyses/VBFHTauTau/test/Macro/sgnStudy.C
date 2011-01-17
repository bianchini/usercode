
 
{
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  
  TTree *treeHPS   = (TTree*)fsgn.Get("etoTauMargLooseNoCracks70/fitter_tree");
  TTree *treeSHC   = (TTree*)fsgn.Get("etoTauSCMargNoCracks70/fitter_tree");

  float nBins = 10;

  TH1F* h = new TH1F("h","", nBins ,15,65);
 
  h->SetXTitle("pt");
  h->SetYTitle("entries/(10 GeV)");


  for(unsigned int i = 1; i<=nBins; i++){
    float ptLow_  = 15 + 50./nBins*i;
    float ptHigh_ = 15 + 50./nBins*(i+1);
    float passing = (float)treeHPS->GetEntries(Form("pt>=%f && tauAntiEMVA>0.5 && abseta>1.5 && mass>65 && mass<120 && mcTrue",ptLow_, ptHigh_));
    float all     = (float)treeHPS->GetEntries(Form("pt>=%f && abseta>1.5 && mass>65 && mass<120  && mcTrue",ptLow_, ptHigh_));
    std::cout << all << "  " << passing << std::endl;
    h->SetBinContent(i,(float)passing/(float)all);
    h->SetBinError(i,TMath::Sqrt(passing/all*(1-passing/all)/all));
  }


  h->Draw("P");

  TLegend* leg = new TLegend(0.6,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->SetHeader("HPS, |#eta|<1.5");
  leg->AddEntry(h,  "e#rightarrow #tau_{had} fake-rate vs p_{t}","L");
  leg->Draw();

  //treeSHC->Draw("pt>>h(50,15,65)","tauAntiEMVA>0.5 && abseta<1.5 && mass>65 && mass<120 && mcTrue");
  
}
