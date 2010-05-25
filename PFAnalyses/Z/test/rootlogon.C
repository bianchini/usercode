void rootlogon() {

  gROOT->SetStyle("Plain");
  gStyle->SetHistMinimumZero(kTRUE);

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);

  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  //gStyle->SetTitleW(1);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);
  //gStyle->SetOptDate(3);
  gStyle->SetOptFit(0111);
  gStyle->SetOptStat(0000000);
  //gStyle->SetTextSize(10);
 
  
  gSystem->Load("libFWCoreFWLite.so");
  //gSystem->Load("libPFAnalysesPFCandidate.so");
  //gSystem->Load("libPFAnalysesRootTools.so");
  AutoLibraryLoader::enable();
  gSystem->Load("libCintex.so");
  ROOT::Cintex::Cintex::Enable();
  
}

