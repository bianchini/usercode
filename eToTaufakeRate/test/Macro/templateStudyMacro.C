
{
  TFile fdata("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TFile fsoup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");

  TTree *fullTreeSgn   = (TTree*)fsgn.Get("etoTauMargLooseNoCracks70/fitter_tree");
  TTree *fullTreeSoup  = (TTree*)fsoup.Get("etoTauMargTightNoCracks60/fitter_tree");
  TTree *fullTreeData  = (TTree*)fdata.Get("etoTauMargLooseNoCracks70/fitter_tree");

  TH1F* hVBTF = new TH1F("hVBTF","",12,60,120);
  TH1F* hData = new TH1F("hData","",12,60,120);
  TH1F* hPF   = new TH1F("hPF",""  ,12,60,120);
 
  hVBTF->SetMarkerStyle(kFullTriangleUp);
  hVBTF->SetMarkerColor(kBlue);
  hVBTF->SetMarkerSize(1.0);
  hData->SetMarkerStyle(kFullCircle);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerSize(1.0);
  hPF->SetLineColor(kRed);
  hVBTF->SetXTitle("mass");
  hVBTF->SetYTitle("entries/(5 GeV)");

  fullTreeSoup->Draw("mass>>hVBTF","abseta<1.5 && leadPFChargedHadrCandTrackPt>25 && leadPFCandPt>15 && signalPFChargedHadrCands<1.5 && electronPreIDOutput<0.6");
  fullTreeSgn->Draw("mass>>hData","abseta<1.5 && leadPFChargedHadrCandTrackPt>25 && leadPFCandPt>15 && signalPFChargedHadrCands<1.5 && electronPreIDOutput<0.6");
  fullTreeSgn->Draw("mass>>hPF","mcTrue && electronPreIDOutput<0.6 && abseta<1.5");
 

  hVBTF->Sumw2();
  hData->Sumw2();
  hVBTF->DrawNormalized("P");
  hData->DrawNormalized("PSAME");
  hPF->DrawNormalized("HISTSAME");

  TLegend* leg = new TLegend(0.12,0.6,0.45,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);

  leg->SetHeader("HPS, |#eta|>1.5");
  leg->AddEntry(hPF,  "MC: passing probes (MC-truth)","L");
  leg->AddEntry(hVBTF,"MC: template from soup","P");
  leg->AddEntry(hData,"MC: template from data","P");

  leg->Draw();
}
