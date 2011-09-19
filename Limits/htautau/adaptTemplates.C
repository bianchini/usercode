



{

  TString channel = "mu";
  TString variable = "diTauVisMass";

  std::vector<TString> analysis;
  analysis.push_back(TString(""));
  analysis.push_back(TString("_CMS_scale_tUp"));
  analysis.push_back(TString("_CMS_scale_tDown"));
  analysis.push_back(TString("_CMS_scale_jUp"));
  analysis.push_back(TString("_CMS_scale_jDown"));
  //analysis.push_back(TString("_CMS_scale_eUp"));
  //analysis.push_back(TString("_CMS_scale_eDown"));

  TFile fin("datacards/"+channel+"TauSM_"+variable+".root","UPDATE");

  for(int k = 0 ; k < analysis.size() ; k++){

    fin.cd(channel+"Tau_SM0");
    TH1F* hQCD_new = (TH1F*)gDirectory->Get("QCD"+analysis[k]);
    TH1F* hW_new   = (TH1F*)gDirectory->Get("W"+analysis[k]);
    TH1F* hVV_new  = (TH1F*)gDirectory->Get("VV"+analysis[k]);
    TH1F* hZL_new  = (TH1F*)gDirectory->Get("ZL"+analysis[k]);
    TH1F* hZJ_new  = (TH1F*)gDirectory->Get("ZJ"+analysis[k]);
    
    
    fin.cd(channel+"Tau_SM2");
    TH1F* hQCD_old = (TH1F*)gDirectory->Get("QCD"+analysis[k]);
    TH1F* hW_old   = (TH1F*)gDirectory->Get("W"+analysis[k]);
    TH1F* hVV_old  = (TH1F*)gDirectory->Get("VV"+analysis[k]);
    TH1F* hZL_old  = (TH1F*)gDirectory->Get("ZL"+analysis[k]);
    TH1F* hZJ_old  = (TH1F*)gDirectory->Get("ZJ"+analysis[k]);
    
    ////////////////////////////////////////////////////////
    float QCD_old = hQCD_old->Integral();
    hQCD_old->Reset();
    hQCD_old->Add(hQCD_new, QCD_old/hQCD_new->Integral());
    hQCD_old->Write("QCD"+analysis[k] , TObject::kOverwrite);
    
    float W_old = hW_old->Integral();
    hW_old->Reset();
    hW_old->Add(hW_new, W_old/hW_new->Integral());
    hW_old->Write("W"+analysis[k] , TObject::kOverwrite);
    
    float VV_old = hVV_old->Integral();
    hVV_old->Reset();
    hVV_old->Add(hVV_new, VV_old/hVV_new->Integral());
    hVV_old->Write("VV"+analysis[k] , TObject::kOverwrite);
    
    float ZL_old = hZL_old->Integral();
    hZL_old->Reset();
    hZL_old->Add(hZL_new, ZL_old/hZL_new->Integral());
    hZL_old->Write("ZL"+analysis[k] , TObject::kOverwrite);
    
    float ZJ_old = hZJ_old->Integral();
    hZJ_old->Reset();
    hZJ_old->Add(hZJ_new, ZJ_old/hZJ_new->Integral());
    hZJ_old->Write("ZJ"+analysis[k] , TObject::kOverwrite);
    ///////////////////////////////////////////////////////

  }


  fin.Close();

}
