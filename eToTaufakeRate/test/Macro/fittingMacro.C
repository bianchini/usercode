#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooNumConvPdf.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooVoigtian.h"
#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;
using namespace std;

void testFit(float xmin=60, float xmax=120)
{

  //using namespace RooFit;

  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile f("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fbkg("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");

  // find the sgn
  fsgn.cd();
  TTree *fulltreesgn = (TTree*)fsgn.Get("etoTauSCMargNoCracks70/fitter_tree");
  float SGN = (float)fulltreesgn->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f && mcTrue",xmin,xmax));
  std::cout << "tree sgn " << SGN << std::endl;

  // find the bkg
  fbkg.cd();
  TTree *fulltreebkg = (TTree*)fbkg.Get("etoTauSCMargNoCracks70/fitter_tree");
  float BKG = (float)fulltreebkg->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree bkg " << BKG << std::endl;

  // tree for dataset
  f.cd();
  TTree *fulltree = (TTree*)f.Get("etoTauSCMargNoCracks70/fitter_tree");
  float ALL =  (float)fulltree->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree all " << ALL << std::endl;

  TFile *dummy = new TFile("dummy.root","RECREATE");
  TTree *tree = (TTree*)fulltree->CopyTree("tauAntiEMVA>0 && abseta>1.5");
  TTree *treesgn = (TTree*)fulltreesgn->CopyTree("tauAntiEMVA>0 && abseta>1.5 && mcTrue");

  RooRealVar  mass("mass","mass (GeV/c^{2})",xmin,xmax);
  RooDataSet data("data","dataset with mass",tree,mass);



  // datset
  RooDataSet data_("data_","dataset with mass",treesgn,mass);

  // crystall ball for signal only
  RooRealVar m1_("m1_","m1",91.2,88,95.);
  RooRealVar sigma_("sigma_","sigma",5.0,2.0,8.0);
  RooRealVar alfa_("alfa_","alfa",1.0,0.5,5.);
  RooRealVar n_("n_","n",1,0.,20.);
  RooCBShape signal_("signal_","crystalBall for signal_",mass,m1_,sigma_,alfa_,n_);

  // fit
  RooFitResult* signalResult = signal_.fitTo(data_,Save(kTRUE),Range(xmin,xmax));

  // crystall ball
  RooRealVar m1("m1","m1",91.2,88,95.);
  RooRealVar sigma("sigma","sigma",5.0,0.5,8.0);
  RooRealVar alfa("alfa","alfa",1.0,0.0,10.);
  RooRealVar n("n","n",1.0,0.5,20.);
  //RooCBShape signal("signal","crystalBall for signal",mass,m1,sigma,alfa,n);

  // constrained crystall ball
  RooRealVar m1c("m1c","m1",m1_.getVal());
  RooRealVar sigmac("sigmac","sigma",sigma_.getVal());
  RooRealVar alfac("alfac","alfa", alfa_.getVal());
  RooRealVar nc("nc","n", n_.getVal());
  RooCBShape signalc("signalc","crystalBall for signal constrained",mass,m1c,sigmac,alfac,nc);

  // constrained crystall ball (X) resol. model
  RooRealVar sigmar("sigmar","sigmar",2.0,0.0,8.0);
  RooRealVar meanr("meanr","meanr",0.1,0,5);
  RooGaussian resol("resol","resol",mass,meanr,sigmar);
  RooNumConvPdf signal("signal","signal",mass,signalc,resol);
  signal.setConvolutionWindow(meanr,sigmar,5);

  // voigtian
  RooRealVar mean("mean","mean",90,80,100);
  RooRealVar width("width","width",2.9);
  RooRealVar sigmaVoig("sigmaVoig","sigmaVoig",5,0.5,10);
  RooVoigtian voig("signalVoig","signalVoig",mass,mean,width,sigmaVoig);
  // bifurcated gaussina
  RooRealVar sigmaR("sigmaR","sigmaR",5,0.0,10);
  RooRealVar sigmaL("sigmaL","sigmaL",5,0.5,20);
  RooBifurGauss bifurc("bifurc","bifurc",mass,mean,sigmaL,sigmaR); 
  RooRealVar fVB("fVB","FVB",0.5,0,1);
  RooAddPdf voigPlusBifurc("voigPlusBifurc","voigPlusBifurc",RooArgList(voig,bifurc),fVB);

  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential bkg("bkg","exponential for bkg",mass,c1);
  
  RooRealVar c2("c2","c2",60,10,100);
  RooRealVar c3("c3","c3",10,0.5,50);
  RooLandau bkg2("bkg2","Landau for bkg",mass,c2,c3);

  // numbers1
  RooRealVar nsig1("nsig1","num of signal events",100,0.,1000000) ;
  RooRealVar nbkg1("nbkg1","num of bkg events",100,0.,1000000) ;
  // numbers2
  RooRealVar nsig2("nsig2","num of signal events",100,0.,1000000) ;
  RooRealVar nbkg2("nbkg2","num of bkg events",100,0.,1000000) ;

  // model 1
  RooAddPdf model1("model1","model1 for signal+bkg",RooArgList(signal,bkg),RooArgList(nsig1,nbkg1));
  // model 2
  RooAddPdf model2("model2","model2 for signal+bkg",RooArgList(voigPlusBifurc,bkg),RooArgList(nsig2,nbkg2));

  // fit
  //RooFitResult* model1Result = model1.fitTo(data,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax), SumW2Error(kTRUE) );
  RooFitResult* model2Result = model2.fitTo(data,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax), SumW2Error(kTRUE) );

  // plot
  RooPlot* frame = mass.frame();
  data.plotOn(frame);

  //model1.plotOn(frame,LineColor(kRed),Range(30,120));
  model2.plotOn(frame,LineColor(kBlue),Range(30,120));
  model2.plotOn(frame,Components(voigPlusBifurc),LineColor(kRed),LineStyle(kDashed),Range(30,120));
  model2.plotOn(frame,Components(bkg),LineColor(kGreen),LineStyle(kDashed),Range(30,120));
  
  //model1.plotOn(frame,Components(signal),LineColor(kBlue),LineStyle(kDashed),Range(30,120));
  //signalc.plotOn(frame,LineColor(kBlack),LineStyle(kDashed));

  TH1F* h1 = new TH1F("h1","",1,0,1);
  TH1F* h2 = new TH1F("h2","",1,0,1);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
  
  std::cout << "Num of signal1 MC " << SGN << " --- num of signal1 from fit " << nsig1.getVal() << std::endl;
  std::cout <<  "Num of bkg1 MC " << BKG << " --- num of bkg1 from fit " << nbkg1.getVal() << std::endl;
  std::cout << "Num of signal2 MC " << SGN << " --- num of signal2 from fit " << nsig2.getVal() << std::endl;
  std::cout <<  "Num of bkg2 MC " << BKG << " --- num of bkg2 from fit " << nbkg2.getVal() << std::endl;
  
  //std::cout << "Fit to the soup " << std::endl;
  //model1Result->Print();
  model2Result->Print();
  
  //std::cout << "Fit to the signal only " << std::endl;
  //signalResult->Print();
  
  leg->SetHeader("Zee+Z#tau#tau+We#nu +t#bar{t}; L = 235 pb^{-1}");
  leg->AddEntry(h2,"Signal+Background","l");
  leg->AddEntry(h1,"Signal component","l");

  frame->Draw();
  leg->Draw();

  c->Draw();

  delete h1;
  delete h2;
    
}

void testFitWeight(float xmin=60, float xmax=120, string fname = "testNewWriteFromPAT_soup", bool useW2=true)
{

  //using namespace RooFit;

  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TFile f(("./"+fname+".root").c_str());

  // find the sgn
  f.cd();
  TTree *fulltree = (TTree*)f.Get("etoTauSCMargNoCracks70/fitter_tree");
  float ALL = (float)fulltree->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree all " << ALL << std::endl;

  TFile *dummy = new TFile("dummy.root","RECREATE");
  TTree *tree = (TTree*)fulltree->CopyTree("tauAntiEMVA>0");

  RooRealVar  weight("weight","weight", 0,99999 );
  RooRealVar  mass("mass","mass (GeV/c^{2})",xmin,xmax);
  RooDataSet data("data","dataset with mass",RooArgSet(mass,weight), Import( *tree ), WeightVar( weight ));

  
  // crystall ball
  RooRealVar m1("m1","m1",91.2,88,95.);
  RooRealVar sigma("sigma","sigma",5.0,2.0,8.0);
  RooRealVar alfa("alfa","alfa",1.0,0.5,5.);
  RooRealVar n("n","n",1.0,0.,20.);
  RooCBShape signal("signal","crystalBall for signal",mass,m1,sigma,alfa,n);
  
  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential bkg("bkg","exponential for bkg",mass,c1);

  // numbers
  RooRealVar nsig("nsig","num of signal events",100,0.,1000000) ;
  RooRealVar nbkg("nbkg","num of bkg events",100,0.,1000000) ;

  // model
  RooAddPdf model("model","model for signal+bkg",RooArgList(signal,bkg),RooArgList(nsig,nbkg));

  // fit
  RooFitResult* modelResult = model.fitTo(data,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax), SumW2Error( useW2 ) );

  // plot
  RooPlot* frame = mass.frame();
  data.plotOn(frame);
  model.plotOn(frame,LineColor(kRed),Range(30,120));
  model.plotOn(frame,Components(bkg),LineColor(kGreen),LineStyle(kDashed),Range(30,120));
  model.plotOn(frame,Components(signal),LineColor(kBlue),LineStyle(kDashed),Range(30,120));
   
  std::cout << "Num of signal MC " << ALL << " --- num of signal from fit " << nsig.getVal() << "+/-" << nsig.getError() << std::endl;

  frame->Draw();

  c->Draw();
    
}


void testCount(float side0_low = 60, float side0_high = 70, float xmin=60, float xmax=120, float side1_low = 110, float side1_high = 120)
{

  //using namespace RooFit;
  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TFile f("./test_all.root");
  TFile fbkg("./test_bkg.root");
  TFile fsgn("./test_sgn.root");

  // find the sgn
  fsgn.cd();
  TTree *fulltreesgn = (TTree*)fsgn.Get("etoTauMargLooseNoCracks70/fitter_tree");
  float SGN = (float)fulltreesgn->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree sgn " << SGN << std::endl;

  // find the bkg
  fbkg.cd();
  TTree *fulltreebkg = (TTree*)fbkg.Get("etoTauMargLooseNoCracks70/fitter_tree");
  float BKG = (float)fulltreebkg->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree bkg " << BKG << std::endl;

  // tree for dataset
  f.cd();
  TTree *fulltree = (TTree*)f.Get("etoTauMargLooseNoCracks70/fitter_tree");
  float ALL =  (float)fulltree->GetEntries(Form("tauAntiEMVA>0 && mass>%f && mass<%f",xmin,xmax));
  std::cout << "tree all " << ALL << std::endl;

  TFile *dummy = new TFile("dummy.root","RECREATE");
  TTree *tree = (TTree*)fulltree->CopyTree("tauAntiEMVA>0");
  TTree *treesgn = (TTree*)fulltreesgn->CopyTree("tauAntiEMVA>0");

  RooRealVar  mass("mass","mass (GeV/c^{2})",xmin,xmax);
  RooDataSet data("data","dataset with mass",tree,mass);


  // extractiong fit parameters from signa only:

  // datset
  RooDataSet data_("data_","dataset with mass",treesgn,mass);

  // crystall ball for signal only
  RooRealVar m1_("m1_","m1",91.2,88,95.);
  RooRealVar sigma_("sigma_","sigma",5.0,2.0,8.0);
  RooRealVar alfa_("alfa_","alfa",1.0,0.5,5.);
  RooRealVar n_("n_","n",1,0.,20.);
  RooCBShape signal_("signal_","crystalBall for signal_",mass,m1_,sigma_,alfa_,n_);

  // fit
  RooFitResult* signalResult = signal_.fitTo(data_,Save(kTRUE),Range(xmin,xmax));



  // crystall ball
  RooRealVar m1("m1","m1",91.2,88,95.);
  RooRealVar sigma("sigma","sigma",5.0,2.0,8.0);
  RooRealVar alfa("alfa","alfa",1.0,0.5,5.);
  RooRealVar n("n","n",1.0,0.,20.);
  RooCBShape signal("signal","crystalBall for signal",mass,m1,sigma,alfa,n);

  // constrained crystall ball
  RooRealVar m1c("m1c","m1",m1_.getVal());
  RooRealVar sigmac("sigmac","sigma",sigma_.getVal());
  RooRealVar alfac("alfac","alfa", alfa_.getVal());
  RooRealVar nc("nc","n", n_.getVal());
  RooCBShape signalc("signalc","crystalBall for signal constrained",mass,m1c,sigmac,alfac,nc);

  
  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential bkg("bkg","exponential for bkg",mass,c1);
  
  RooRealVar c2("c2","c2",60,10,100);
  RooRealVar c3("c3","c3",10,0.5,50);
  RooLandau bkg2("bkg2","Landau for bkg",mass,c2,c3);

  // numbers
  RooRealVar nsig("nsig","num of signal events",100,0.,10000) ;
  RooRealVar nbkg("nbkg","num of bkg events",100,0.,10000) ;

  // model
  RooAddPdf model("model","model for signal+bkg",RooArgList(signal,bkg),RooArgList(nsig,nbkg));

  // fit
  RooFitResult* modelResult = model.fitTo(data,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax));

  // plot
  RooPlot* frame = mass.frame();
  data.plotOn(frame);
  model.plotOn(frame,LineColor(kRed),Range(30,120));
  model.plotOn(frame,Components(bkg),LineColor(kGreen),LineStyle(kDashed),Range(30,120));
  model.plotOn(frame,Components(signal),LineColor(kBlue),LineStyle(kDashed),Range(30,120));
  //signalc.plotOn(frame,LineColor(kBlack),LineStyle(kDashed));

  float peak = tree->GetEntries( Form("mass>%f && mass<%f",xmin,xmax) );
  float sidebands0 = tree->GetEntries( Form("(mass>%f && mass<%f)",side0_low,side0_high) );
  float sidebands1 = tree->GetEntries( Form("(mass>%f && mass<%f)",side1_low,side1_high) );
  float s1 = side0_high-side0_low;
  float s2 = side1_high-side1_low;
  float s3 = xmax-xmin;
  float x1 = (side0_high+side0_low)/2;
  float x2 = (side1_high+side1_low)/2;
  float x3 = (xmax+xmin)/2;

  float bkgInPeak = s3*( x3/(x1-x2)*(sidebands0/s1 - sidebands1/s2) + (sidebands1/s2*x1 - sidebands0/s1*x2)/(x1-x2));
  std::cout << "bkgInPeak " << bkgInPeak << std::endl;

  //correction for signal in sidebands:
  float N = tree->GetEntries(Form("(mass>%f && mass<%f)",xmin,xmax));
  float k1 = s3/s1*(x3-x2)/(x1-x2);
  float k2 = s3/s2*(x3-x1)/(x1-x2);
  float f1= (float)(treesgn->GetEntries( Form("(mass>%f && mass<%f)",side0_low,side0_high) ))/N;
  float f2= (float)(treesgn->GetEntries( Form("(mass>%f && mass<%f)",side1_low,side1_high) ))/N;
  float Ns = (N-k1*sidebands0+k2*sidebands1)/(1-k1*f1+k2*f2);
  float f3 = (float)(treesgn->GetEntries( Form("(mass>%f && mass<%f)",xmin,xmax) ))/(float)(treesgn->GetEntries());
  float NsAll = Ns/f3;

  std::cout << "k1,k2 " << k1 << "," << k2 << std::endl;
  std::cout << "f1,f2,f3 " << f1 << "," << f2 << "," << f3 << std::endl;
  std::cout << "N, Ns " << N << "," << Ns << std::endl;

  std::cout << "Num of signal MC " << SGN << " --- num of signal from fit " << nsig.getVal()  << std::endl;
  std::cout << " Num of signal in region [" << xmin << ","<< xmax << "]: MC " << treesgn->GetEntries(Form("(mass>%f && mass<%f)",xmin,xmax)) << " , from sidebands " <<  N-bkgInPeak << " ; from sidebands corrected " << Ns << std::endl;
  std::cout << " Num of signal all " << treesgn->GetEntries()  << " , from sidebands with correction " <<  NsAll << std::endl;
  std::cout << "Eff MC " << (float)(treesgn->GetEntries( Form("(mass>%f && mass<%f)",xmin,xmax) ))/(float)(fulltreesgn->GetEntries(Form("mass>%f && mass<%f",xmin,xmax))) << " eff sidebands corrected " << Ns/(float)(fulltree->GetEntries(Form("mass>%f && mass<%f",xmin,xmax))) << std::endl;
  
  

  frame->Draw();

  c->Draw();
    
}


void manyFit( float xmin = 60, float xmax = 120, int nexp = 1, string tnp = "etoTauMargLooseNoCracks70", string cut = "tauAntiEMVA>0"){
  
  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TH1F* h1 = new TH1F("h1","signal purity",50,0,1);
  TH1F* h1c = new TH1F("h1c","signal purity",50,0,1);
  TH1F* h1_ = new TH1F("h1_","signal purity",50,0,1);

  TFile f("./test_all.root");
  TFile fbkg("./test_bkg.root");
  TFile fsgn("./test_sgn.root");

  // find the sgn
  fsgn.cd();
  TTree *fulltreesgn = (TTree*)fsgn.Get((tnp+"/fitter_tree").c_str());
  float SGN = (float)fulltreesgn->GetEntries(Form("%s && mass>%f && mass<%f && mcTrue",cut.c_str(),xmin,xmax));
  std::cout << "tree sgn " << SGN << std::endl;
  
 
  // find the bkg
  fbkg.cd();
  TTree *fulltreebkg = (TTree*)fbkg.Get((tnp+"/fitter_tree").c_str());
  float BKG = (float)fulltreebkg->GetEntries(Form("%s && mass>%f && mass<%f",cut.c_str(),xmin,xmax));
  std::cout << "tree bkg " << BKG << std::endl;

 
  // tree for dataset
  f.cd();
  TTree *fulltree = (TTree*)f.Get((tnp+"/fitter_tree").c_str());
  float ALL =  (float)fulltree->GetEntries(Form("%s && mass>%f && mass<%f",cut.c_str(),xmin,xmax));
  std::cout << "tree all " << ALL << std::endl;

  // trees with selection
  TFile *dummy = new TFile("dummy.root","RECREATE");
  TTree *tree = (TTree*)fulltree->CopyTree( cut.c_str() );
  TTree *treesgn = (TTree*)fulltreesgn->CopyTree( ("mcTrue && "+cut).c_str() ); 
  TTree *treebkg = (TTree*)fulltreebkg->CopyTree( cut.c_str() ); 

  // mass variable
  RooRealVar  mass("mass","mass (GeV/c^{2})",xmin-15,xmax+15);
  RooDataSet data("data","dataset with mass",tree,mass);

  // constarined pdf
  // datset
  RooDataSet datasig_("datasig_","dataset with mass",treesgn,mass);
  RooDataSet databkg_("databkg_","dataset with mass",treebkg,mass);
  
   // voigtian
  RooRealVar mean_("mean_","mean",90,80,100);
  RooRealVar width_("width_","width",2.9);
  RooRealVar sigmaVoig_("sigmaVoig_","sigmaVoig",5,0.5,10);
  RooVoigtian voig_("signalVoig_","signalVoig",mass,mean_,width_,sigmaVoig_);
  // bifurcated gaussian
  RooRealVar sigmaR_("sigmaR_","sigmaR",5,0.0,10);
  RooRealVar sigmaL_("sigmaL_","sigmaL",5,0.5,30);
  RooBifurGauss bifurc_("bifurc_","bifurc",mass,mean_,sigmaL_,sigmaR_); 
  RooRealVar fVB_("fVB_","FVB",0.5,0,1);
  RooAddPdf signal_("signal_","voigPlusBifurc",RooArgList(voig_,bifurc_),fVB_);
  // fit
  RooFitResult* signalResult = signal_.fitTo(datasig_,Save(kTRUE),Range(xmin,xmax));

  // constrained exp
  RooRealVar c1_("c1_","c1_",0,-10,10.);
  RooExponential bkg_("bkg_","exponential for bkg",mass,c1_);
  // fit
  RooFitResult* bkgResult = bkg_.fitTo(databkg_,Save(kTRUE),Range(xmin,xmax));
  // model
  float MCfraction = SGN/ALL;
  RooRealVar f_("f_","fraction", MCfraction );
  RooAddPdf model_("model_","model for signal+bkg",RooArgList(signal_,bkg_),f_);
  // constrained voigtian
  RooRealVar meanc("meanc","mean",mean_.getVal());
  RooRealVar widthc("widthc","width",width_.getVal());
  RooRealVar sigmaVoigc("sigmaVoigc","sigmaVoig",sigmaVoig_.getVal());
  RooVoigtian voigc("signalVoigc","signalVoig",mass,meanc,widthc,sigmaVoigc);
  // constrained bifurcated gaussian
  RooRealVar sigmaRc("sigmaRc","sigmaR",sigmaR_.getVal());
  RooRealVar sigmaLc("sigmaLc","sigmaL",sigmaL_.getVal());
  RooBifurGauss bifurcc("bifurcc","bifurc",mass,meanc,sigmaLc,sigmaRc); 
  RooRealVar fVBc("fVBc","FVB",fVB_.getVal());
  RooAddPdf signalc("signalc","voigPlusBifurc",RooArgList(voigc,bifurcc),fVBc);

  // pdf from data
  RooKeysPdf pdfFromData("pdfFromData","non-parametric pdf from data",mass,data);
  mass.setRange(0,xmin,xmax);

  // voigtian
  RooRealVar mean("mean","mean",90,80,100);
  RooRealVar width("width","width",2.9);
  RooRealVar sigmaVoig("sigmaVoig","sigmaVoig",5,0.5,10);
  RooVoigtian voig("signalVoig","signalVoig",mass,mean,width,sigmaVoig);
  // bifurcated gaussina
  RooRealVar sigmaR("sigmaR","sigmaR",5,0.0,10);
  RooRealVar sigmaL("sigmaL","sigmaL",10,0.5,50);
  RooBifurGauss bifurc("bifurc","bifurc",mass,mean,sigmaL,sigmaR); 
  RooRealVar fVB("fVB","FVB",0.5,0,1);
  RooAddPdf signal("signal","voigPlusBifurc",RooArgList(voig,bifurc),fVB);


  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential bkg("bkg","exponential for bkg",mass,c1);

  // numbers
  RooRealVar nsig("nsig","num of signal events",100,0.,10000) ;
  RooRealVar nbkg("nbkg","num of bkg events",100,0.,10000) ;

  // model
  RooAddPdf model("model","model for signal+bkg",RooArgList(signal,bkg),RooArgList(nsig,nbkg));
  // model c
  RooAddPdf modelc("modelc","model for signal+bkg",RooArgList(signalc,bkg),RooArgList(nsig,nbkg));

  // loop
  for(int i = 0 ; i < nexp; i++){

    RooDataSet* pseudodata = 0;
    if(i == 0) pseudodata = &data;
    else pseudodata = model_.generate(mass, ALL ,Extended());
    
    RooFitResult* modelResult = model.fitTo(*pseudodata,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax));
    h1->Fill(nsig.getVal()/(nsig.getVal()+nbkg.getVal()),1);
    modelResult->Print();
  }

  // loop
  for(int i = 0 ; i < nexp; i++){

    RooDataSet* pseudodata = 0;
    if(i == 0) pseudodata = &data;
    else pseudodata = model_.generate(mass, ALL ,Extended());
    
    RooFitResult* modelResult = modelc.fitTo(*pseudodata,Extended(kTRUE),Save(kTRUE),Range(xmin,xmax));
    h1c->Fill(nsig.getVal()/(nsig.getVal()+nbkg.getVal()),1);
    modelResult->Print();
  }

  h1_->SetXTitle("fitted purity");
  h1_->SetYTitle("pseudo experiments");
  h1_->Fill( SGN/ALL,1);
  h1_->SetMarkerStyle(20);
  h1_->SetMarkerColor(kRed);
  h1_->Sumw2();
  h1_->Draw("P");

  h1->SetFillStyle(3004);
  h1->SetFillColor(kRed);
  h1c->SetFillStyle(3005);
  h1c->SetFillColor(kBlue);
  h1->DrawNormalized("HISTSAME");
  h1c->DrawNormalized("HISTSAME");

  leg->AddEntry(h1_,"MC purity","p");
  leg->AddEntry(h1,"Voigtian+BifurcatedGaussian+exp: all floating","f");
  leg->AddEntry(h1c,"Voigtian+BifurcatedGaussian+exp: fixed to signal-only fit","f");
 
  //RooPlot* frame = mass.frame();
  //model_.plotOn(frame);
  //frame->Draw();
  
  leg->Draw();

  c->Draw();
    
}


void testOnData(string tnp = "etoTauMargLooseNoCracks70", string cut = "tauAntiEMVA>0", 
		float maxEta=1.5,float minEta=0,float maxPt=9999,float minPt=0){

  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  //leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile fdata("../tagAndProbe/trees/38X/testNewWriteFromPAT_Data_runSplit.root");
  TFile fsoup("../tagAndProbe/trees/38X/testNewWriteFromPAT_soup.root");
  TFile fsgn("../tagAndProbe/trees/38X/testNewWriteFromPAT_Zee.root");

  // tree for data
  fdata.cd();
  TTree *fTdata = (TTree*)fdata.Get((tnp+"/fitter_tree").c_str()); 
  TFile *dummy1 = new TFile("dummy_data.root","RECREATE");
  TTree *cTdata = (TTree*)fTdata->CopyTree( Form("(abseta<%f) && (abseta>%f) && (pt<%f) && (pt>%f) && %s",maxEta,minEta,maxPt,minPt,cut.c_str()) );
  std::cout << Form("(abseta<%f) && (abseta>%f) && (pt<%f) && (pt>%f) && (%s)",maxEta,minEta,maxPt,minPt,cut.c_str()) << std::endl;
  std::cout << "All data = " << fTdata->GetEntries() << "  Selected data = " << cTdata->GetEntries() << std::endl;

  RooRealVar mdata("mass","mass from data (GeV/c^{2})",60,120);
  RooDataSet data("data","data set with mass",cTdata,mdata);


  // tree for soup
  //fsoup.cd();
  TTree *fTsoup = (TTree*)fsoup.Get((tnp+"/fitter_tree").c_str());
  TFile *dummy2 = new TFile("dummy_soup.root","RECREATE");
  TTree *cTsoup = (TTree*)fTsoup->CopyTree( Form("abseta<%f && abseta>%f && pt<%f && pt>%f && %s",maxEta,minEta,maxPt,minPt,cut.c_str()) );
  std::cout << "All soup = " << fTsoup->GetEntries() << "  Selected soup = " << cTsoup->GetEntries() << std::endl;
  RooRealVar  msoup("mass","mass from soup (GeV/c^{2})",60,120);
  RooDataSet soup("soup","soup set with mass",cTsoup,msoup);

  // tree for signal
  //fsgn.cd();
  TTree *fTsgn = (TTree*)fsgn.Get((tnp+"/fitter_tree").c_str());
  TFile *dummy3 = new TFile("dummy_sgn.root","RECREATE");
  TTree *cTsgn = (TTree*)fTsgn->CopyTree( Form("abseta<%f && abseta>%f && pt<%f && pt>%f && %s",maxEta,minEta,maxPt,minPt,cut.c_str()) );
  std::cout << "All signal = " << fTsgn->GetEntries() << "  Selected signal = " << cTsgn->GetEntries() << std::endl;
  RooRealVar  msgn("mass","mass from sgn (GeV/c^{2})",60,120);
  RooDataSet sgn("sgn","sgn set with mass",cTsgn,msgn);

  RooRealVar sigmar("sigmar","sigmar",2.0,0.0,8.0);
  RooRealVar meanr("meanr","meanr",0.1,0,10);
  RooGaussian resol("resol","resol",mdata,meanr,sigmar);

  // model signal shape from MC:
  RooRealVar m1_("m1_","m1",91.2,88,95.);
  RooRealVar sigma_("sigma_","sigma",5.0,2.0,8.0);
  RooRealVar alfa_("alfa_","alfa",1.0,0.5,5.);
  RooRealVar n_("n_","n",1,0.,20.);
  RooCBShape signal_("signal_","crystalBall for signal_",msgn,m1_,sigma_,alfa_,n_);

  RooNumConvPdf signal_Conv("signal_Conv","signal",msgn,signal_,resol);
  signal_Conv.setConvolutionWindow(meanr,sigmar,5);
  signal_Conv.fitTo( sgn ,Save(kTRUE),Range(60,120));

  // voigtian
  RooRealVar mean_("mean_","mean",90,80,100);
  RooRealVar width_("width_","width",2.9);
  RooRealVar sigmaVoig_("sigmaVoig_","sigmaVoig",5,0.5,10);
  RooVoigtian voig_("signalVoig_","signalVoig",msgn,mean_,width_,sigmaVoig_);
  // bifurcated gaussian
  RooRealVar sigmaR_("sigmaR_","sigmaR",5,0.0,10);
  RooRealVar sigmaL_("sigmaL_","sigmaL",5,0.5,30);
  RooBifurGauss bifurc_("bifurc_","bifurc",msgn,mean_,sigmaL_,sigmaR_); 
  RooRealVar fVB_("fVB_","FVB",0.5,0,1);
  RooAddPdf voigPlusBifurc_("voigPlusBifurc_","voigPlusBifurc",RooArgList(voig_,bifurc_),fVB_);

  // fit
  RooFitResult* signalResult = signal_.fitTo( sgn ,Save(kTRUE),Range(60,120));
  RooFitResult* voigPlusBifurcResult = voigPlusBifurc_.fitTo( sgn ,Save(kTRUE),Range(60,120));

  // constrained crystall ball
  RooRealVar m1c("m1c","m1",m1_.getVal());
  RooRealVar sigmac("sigmac","sigma",sigma_.getVal());
  RooRealVar alfac("alfac","alfa", alfa_.getVal());
  RooRealVar nc("nc","n", n_.getVal());
  RooCBShape signalc("signalc","crystalBall for signal constrained",mdata,m1c,sigmac,alfac,nc);

  // constrained crystall ball (X) resol. model
  RooRealVar m1r("m1r","m1",0.5,0,10);
  RooRealVar alfar("alfar","alfa",1.0,0.5,5.);
  RooRealVar nr("nr","n",1,0.,20.);
  RooCBShape resol2("resol2","crystalBall for signal constrained",mdata,m1r,sigmar,alfar,nr);
  RooNumConvPdf signalConv("signalConv","signal",mdata,signalc,resol);
  signalConv.setConvolutionWindow(meanr,sigmar,5);

  // constrained voigtian
  RooRealVar meanc("meanc","mean",91,85,95/*mean_.getVal()*/);
  RooRealVar widthc("widthc","width",width_.getVal());
  RooRealVar sigmaVoigc("sigmaVoigc","sigmaVoig",sigmaVoig_.getVal());
  RooVoigtian voigc("signalVoigc","signalVoig",msgn,meanc,widthc,sigmaVoigc);
  // constrained bifurcated gaussian
  RooRealVar sigmaRc("sigmaRc","sigmaR",sigmaR_.getVal());
  RooRealVar sigmaLc("sigmaLc","sigmaL",sigmaL_.getVal());
  RooBifurGauss bifurcc("bifurcc","bifurc",msgn,meanc,sigmaLc,sigmaRc); 
  RooRealVar fVBc("fVBc","FVB",fVB_.getVal());
  RooAddPdf voigPlusBifurcc("voigPlusBifurcc","voigPlusBifurc",RooArgList(voigc,bifurcc),fVBc);

  RooNumConvPdf voigPlusBifurccConv("voigPlusBifurccConv","signal",mdata,voigPlusBifurcc,resol);
  voigPlusBifurccConv.setConvolutionWindow(meanr,sigmar,5);


  // floating crystall ball
  RooRealVar m1("m1","m1",91.2,88,95.);
  RooRealVar sigma("sigma","sigma",5.0,2.0,8.0);
  RooRealVar alfa("alfa","alfa",1.0,0.5,5.);
  RooRealVar n("n","n",1.0,0.,20.);
  //RooCBShape signal("signal","crystalBall for signal",mdata,m1,sigma,alfa,n);

  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential bkg("bkg","exponential for bkg",mdata,c1);

  // numbers
  RooRealVar nsig("nsig","num of signal events",100,0.,10000) ;
  RooRealVar nbkg("nbkg","num of bkg events",100,0.,10000) ;

  // model
  RooAddPdf model("model","model for signal+bkg",RooArgList(voigPlusBifurccConv,bkg),RooArgList(nsig,nbkg));
  // model c
  RooAddPdf modelc("modelc","model for signal+bkg",RooArgList(voigPlusBifurcc,bkg),RooArgList(nsig,nbkg));

  RooFitResult* modelResults = model.fitTo( data ,Save(kTRUE),Range(60,120),Extended(kTRUE));

  RooFitResult* modelcResults = modelc.fitTo( data ,Save(kTRUE),Range(60,120),Extended(kTRUE));

  RooArgSet modelFit(modelResults->floatParsFinal());
  RooArgSet modelcFit(modelcResults->floatParsFinal());

  RooRealVar *sigFit = (RooRealVar*)(&modelFit["nsig"]);
  RooRealVar *sigcFit = (RooRealVar*)(&modelcFit["nsig"]);

  std::cout << "All floating: nsig = " << sigFit->getVal() << "+/-" << sigFit->getError() << std::endl;
  std::cout << "Constrained: nsig = " << sigcFit->getVal() << "+/-" << sigcFit->getError() << std::endl;

  c->Divide(1,2);

  c->cd(1);
  // plot
  RooPlot* frame = mdata.frame();
  data.plotOn(frame);
  model.plotOn(frame,LineColor(kRed),Range(60,120));
  modelc.plotOn(frame,LineColor(kBlue),Range(60,120));
  frame->Draw();

  TH1F* h1 = new TH1F("h1","",1,0,1);
  TH1F* h2 = new TH1F("h2","",1,0,1);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);

  leg->SetHeader("7 TeV data: L = 10.5 pb^{-1}");
  leg->AddEntry(h1,Form("all shape parameters floating: %.0f+/-%.0f",sigFit->getVal(),sigFit->getError()),"l");
  leg->AddEntry(h2,Form("shape parameters constrained to simulation: %.0f+/-%.0f",sigcFit->getVal(),sigcFit->getError()),"l");
  leg->Draw();

  c->cd(2);
  RooPlot* frameSgn = msgn.frame();
  sgn.plotOn(frameSgn);
  voigPlusBifurcc.plotOn(frameSgn,LineColor(kRed),Range(60,120));
  signal_Conv.plotOn(frameSgn,LineColor(kGreen),Range(60,120));
  frameSgn->Draw();

  c->Update();
  c->Draw();

  //delete dummy1; delete dummy2; delete dummy3;
}



void weightStudy(int nToys_ = 500, Bool_t w2Errors_ = kTRUE, bool  applyWeight_= true, bool addBkg_ = true, int seed_ = 42){


  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  RooRandom::randomGenerator()->SetSeed(seed_); 

  TFile *test = new TFile(TString::Format("toy_seed%d.root",seed_),"RECREATE");

  TTree *tRes = new TTree("tRes","my tree for the residuals");
  double mean_,sigma_,nSgn_; 
  tRes->Branch("mean", &mean_ ,"mean/D");
  tRes->Branch("sigma", &sigma_ ,"sigma/D");
  tRes->Branch("nSgn", &nSgn_ ,"nSgn/D");

  for(int i = 0; i<nToys_; i++){

    RooWorkspace *w = new RooWorkspace("w","w");
    w->factory("mass[0,20]");
    w->factory("weight[1,2]");
    w->factory("Uniform::weightPdf(weight)");
    w->factory("Gaussian::signalPdf(mass,mean[10.,5,15],sigma[1,0,10])");
    w->factory("RooExponential::backgroundPdf(mass,cF[-0.1,-10,10])");
    w->factory("numSignal[100,0,100000]");
    w->factory("numBackground[100,0,100000]");
    w->factory("SUM::model(numSignal*signalPdf, numBackground*backgroundPdf)");
    w->factory("PROD::backgroundXweightPdf(backgroundPdf,weightPdf)");

    TTree *t = new TTree("tree","my tree");
    double mass_,weight_;
    t->Branch("mass", &mass_ ,"mass/D");
    t->Branch("weight", &weight_ ,"weight/D");

    RooRealVar  *mass = w->var("mass");
    RooRealVar  *weight = w->var("weight");
    RooAbsPdf   *signalPdf = w->pdf("signalPdf");
    RooAbsPdf   *backgroundXweightPdf = w->pdf("backgroundXweightPdf");
    RooArgSet   pairObs(*mass,*weight);


    RooDataSet *signal = signalPdf->generate(*mass,50);
    for(int iSgn = 0; iSgn< signal->numEntries(); iSgn++){
      const RooArgSet *sgn = signal->get(iSgn);
      mass_ =  sgn->getRealValue("mass");
      weight_ = 1.0;
      t->Fill();
    }
    RooDataSet *background = backgroundXweightPdf->generate(pairObs,10);
    for(int iSgn = 0; iSgn<background->numEntries(); iSgn++){
      if(!addBkg_) continue;
      const RooArgSet *sgn = background->get(iSgn);
      mass_ =  sgn->getRealValue("mass");
      if( applyWeight_ ) weight_ = sgn->getRealValue("weight");
      else weight_ = 1.0;
      t->Fill();
    }

    cout << "tree with N entries " << t->GetEntries() << endl;
    cout << "... and with weight*N " << t->GetEntries("weight") << endl;
    
    RooDataSet data0("data","dataset", RooArgSet(*mass,*weight), Import(*t));
    RooDataSet data(data0.GetName(),data0.GetTitle(),&data0,*(data0.get()),0,"weight") ;
    
    RooFitResult *fitResults = w->pdf("model")->fitTo(data,Extended(kTRUE),SumW2Error( w2Errors_ ),Save(1));
    
    RooArgSet fpf(fitResults->floatParsFinal());
    
    //fitResults->Print();

    RooRealVar *fitMean = (RooRealVar*)(&fpf["mean"]);
    RooRealVar *fitSigma = (RooRealVar*)(&fpf["sigma"]);
    RooRealVar *fitSgn = (RooRealVar*)(&fpf["numSignal"]);
    mean_ = (double)fitMean->getVal();
    sigma_ = (double)fitSigma->getVal();
    nSgn_ = (double)fitSgn->getVal();
    

    tRes->Fill();
    
    if(false){
      RooPlot* frame = mass->frame(Title("weighted sample")) ;
      data.plotOn(frame);
      w->pdf("model")->plotOn(frame,LineColor(kBlue));
      w->pdf("model")->plotOn(frame,Components("signalPdf"),LineColor(kRed),LineStyle(kDashed));
      frame->Draw();
    }
    //cout << "initial tree" << endl;
    //t->Scan("mass:weight");
    //cout << "******************" << endl;
    //cout << "tree in RooDataSet0" << endl;
    //TTree* data0Tree = (TTree*) data0.tree();
    //data0Tree->Scan("mass:weight");
    //cout << "tree in RooDataSet" << endl;
    //TTree* dataTree = (TTree*) data.tree();
    //dataTree->Scan("mass:weight");
    
    delete t; delete w;
  }

  TH1F* hMean = new TH1F("hMean","pull for the mean",100,8,12);
  TH1F* hSigma = new TH1F("hSigma","pull for the sigma",150,0.0,1.5);
  TH1F* hnSgn = new TH1F("hnSgn","pull for the nSign",100,0,100);
  tRes->Draw("mean>>hMean");
  tRes->Draw("sigma>>hSigma");
  tRes->Draw("nSgn>>hnSgn");

  c->Divide(2,2);

  c->cd(1);
  hMean->Draw("HIST");
  c->cd(2);
  hSigma->Draw("HIST");
  c->cd(3);
  hnSgn->Draw("HIST");

  c->Update();
  c->Draw();

}


void testOnData2(string tnp = "etoTauMargLooseNoCracks70", string cut = "tauAntiEMVA>0", 
		float maxEta=1.5,float minEta=0,float maxPt=9999,float minPt=0){

  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  //leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");

  // tree for signal
  //fsgn.cd();
  TTree *fTsgn = (TTree*)fsgn.Get((tnp+"/fitter_tree").c_str());
  TFile *dummy3 = new TFile("dummy_sgn.root","RECREATE");
  TTree *cTsgn = (TTree*)fTsgn->CopyTree( Form("abseta<%f && abseta>%f && pt<%f && pt>%f && %s",maxEta,minEta,maxPt,minPt,cut.c_str()) );
  std::cout << "All signal = " << fTsgn->GetEntries() << "  Selected signal = " << cTsgn->GetEntries() << std::endl;
 
  RooRealVar  msgn("mass","mass from sgn (GeV/c^{2})",60,120);
  RooDataSet sgn("sgn","sgn set with mass",cTsgn,msgn);

 
  // voigtian
  RooRealVar mean_("mean_","mean",90,80,100);
  RooRealVar width_("width_","width",2.49);
  RooRealVar sigmaVoig_("sigmaVoig_","sigmaVoig",5,0.5,10);
  RooVoigtian voig_("signalVoig_","signalVoig",msgn,mean_,width_,sigmaVoig_);
  // bifurcated gaussian
  RooRealVar sigmaR_("sigmaR_","sigmaR",5,0.0,10);
  RooRealVar sigmaL_("sigmaL_","sigmaL",5,0.5,30);
  RooBifurGauss bifurc_("bifurc_","bifurc",msgn,mean_,sigmaL_,sigmaR_); 
  RooRealVar fVB_("fVB_","FVB",0.5,0,1);
  RooAddPdf signal1("signal1","voigPlusBifurc",RooArgList(voig_,bifurc_),fVB_);

  // voigtian2
  RooRealVar mean2_("mean2_","mean",90,80,100);
  RooRealVar width2_("width2_","width",2.49);
  RooRealVar sigmaVoig2_("sigmaVoig2_","sigmaVoig",5,0.5,10);
  RooVoigtian voig2_("signalVoig2_","signalVoig",msgn,mean2_,width2_,sigmaVoig2_);

  RooRealVar sigmar_("sigmar_","sigmar",2.0,0.0,8.0);
  RooRealVar meanr_("meanr_","meanr",0.1,0,10);
  //RooGaussian resol("resol","resol",msgn,meanr_,sigmar_);
  // crystall ball
  RooRealVar m1("m1","m1",0.5,0.,15);
  RooRealVar sigma("sigma","sigma",5.0,0.0,8.0);
  RooRealVar alfa("alfa","alfa",1.0,0.0,10.);
  RooRealVar n("n","n",1.0,0.5,20.);
  RooCBShape resol("resol","crystalBall for signal",msgn,m1,sigma,alfa,n);


  msgn.setBins(10000,"fft");
  RooFFTConvPdf signal2("signal2","",msgn,voig2_,resol);

  // fit
  RooFitResult* signal1Result = signal1.fitTo(sgn,Save(kTRUE),Range(60,120));
  RooFitResult* signal2Result = signal2.fitTo(sgn,Save(kTRUE),Range(60,120));

  

  c->cd();
  // plot
  RooPlot* frame = msgn.frame(Bins(30));
  sgn.plotOn(frame);
  signal1.plotOn(frame,LineColor(kRed),Range(60,120));
  signal2.plotOn(frame,LineColor(kBlue),Range(60,120), LineStyle(kDashed));
  frame->Draw();

  TH1F* h1 = new TH1F("h1","",1,0,1);
  TH1F* h2 = new TH1F("h2","",1,0,1);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h2->SetLineStyle(7);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);

  leg->AddEntry(h1,"Voigtian plus Bifurcated gaussian","l");
  leg->AddEntry(h2,"Voigtian(X)Crystall Ball","l");
  leg->Draw();

 
  c->Draw();

  signal1Result->Print();
  signal2Result->Print();

  //delete dummy1; delete dummy2; delete dummy3;
}

void testOnData3(string tnp = "etoTauMargLooseNoCracks70", string cut = "tauAntiEMVA>0", 
		float maxEta=1.5,float minEta=0,float maxPt=9999,float minPt=0){

  TCanvas *c = new TCanvas("c","Canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  //leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");

  // tree for signal
  //fsgn.cd();
  TTree *fTsgn = (TTree*)fsgn.Get((tnp+"/fitter_tree").c_str());
  TFile *dummy3 = new TFile("dummy_sgn.root","RECREATE");
  TTree *cTsgn = (TTree*)fTsgn->CopyTree( Form("abseta<%f && abseta>%f && pt<%f && pt>%f && %s  && signalPFChargedHadrCands<1.5",maxEta,minEta,maxPt,minPt,cut.c_str()) );
  std::cout << "All signal = " << fTsgn->GetEntries() << "  Selected signal = " << cTsgn->GetEntries() << std::endl;
 
  RooRealVar  msgn("mass","mass from sgn (GeV/c^{2})",60,120);
  RooDataSet sgn("sgn","sgn set with mass",cTsgn,msgn);

  // exp
  RooRealVar c1("c1","c1",0,-10,10.);
  RooExponential signal1("signal1","exponential for bkg",msgn,c1);

  RooRealVar c_0("c_0",0,-100,100);
  RooRealVar c_1("c_1",0,-100,100);
  RooRealVar c_2("c_2",0,-100,100);
  RooRealVar c_3("c_3",0,-100,100);

  RooPolynomial signal2("signal2","",msgn,RooArgList(c_0,c_1,c_2,c_3));
  // fit
  RooFitResult* signal1Result = signal1.fitTo(sgn,Save(kTRUE),Range(65,120));
  RooFitResult* signal2Result = signal2.fitTo(sgn,Save(kTRUE),Range(65,120));

  

  c->cd();
  // plot
  RooPlot* frame = msgn.frame(Bins(30));
  sgn.plotOn(frame);
  signal1.plotOn(frame,LineColor(kRed),Range(65,120));
  //signal2.plotOn(frame,LineColor(kBlue),Range(60,120), LineStyle(kDashed));
  frame->Draw();

  TH1F* h1 = new TH1F("h1","",1,0,1);
  TH1F* h2 = new TH1F("h2","",1,0,1);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h2->SetLineStyle(7);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
			      
  leg->AddEntry(h1,"exponential");
  //leg->AddEntry(h2,"polynomial");
  leg->Draw();

 
  c->Draw();

  signal1Result->Print();
  //delete dummy1; delete dummy2; delete dummy3;
}
