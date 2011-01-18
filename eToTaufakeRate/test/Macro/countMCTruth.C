
{
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");

  TTree *fullTreeSgn   = (TTree*)fsgn.Get("etoTauMargTightNoCracks90/fitter_tree");
  
  float all     = fullTreeSgn->GetEntries("mcTrue && abseta>1.5 && mass>65 && mass<120");
  float passing = fullTreeSgn->GetEntries("mcTrue && abseta>1.5 && mass>65 && mass<120 && tauAntiEMVA>0.5");
  float eff   = passing/all;
  float error = TMath::Sqrt(eff*(1-eff)/all);

  std::cout << eff << "+/-" << error << std::endl;

}
