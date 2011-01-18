

{
  
  std::string files[6] = {"_etoTauMargLooseNoCracks70_tauAntiEMVA_abseta",
			  "_etoTauMargLooseNoCracks70_-electronPreIDOutput_abseta",
			  "_etoTauMargLooseNoCracks70_matchedID_abseta",
			  "_etoTauSCMargNoCracks80_tauAntiEMVA_abseta",
			  "_etoTauSCMargNoCracks60_-electronPreIDOutput_abseta",
			  "_etoTauSCMargNoCracks70_matchedID_abseta" }


  for(int i = 0 ; i < 6; i++){
    TFile f(("EtoTauPlots"+files[i]+".root").c_str(),"READ");
    if(f.IsZombie()) continue;
    f.cd();
    TCanvas *c0 = gDirectory->Get("effCanvas");
    c0->SaveAs( ("fitEfficiency"+files[i]+".eps").c_str() );
    f.cd("bin0.75");
    TCanvas *c1 = gDirectory->Get("fitCanvas");
    c1->SaveAs( ("fitCanvas"+files[i]+"_bin1.eps").c_str() );
    f.cd("bin2.00");
    TCanvas *c2 = gDirectory->Get("fitCanvas");
    c2->SaveAs( ("fitCanvas"+files[i]+"_bin2.eps").c_str() );
  }


  for(int i = 0 ; i < 6; i++){
    continue;
    TFile f(("EtoTauPlots"+files[i]+".root").c_str(),"READ");
    if(f.IsZombie()) continue;
    f.cd();
    TCanvas *c0 = (TCanvas *)gDirectory->Get("effCanvas");
    c0->SetName(Form("c%d",i));
    c0->Draw();
  }
  
}
