


{

 gSystem->Load("libRooFit") ;


 TFile f("../tagAndProbe/trees/testNewWriteFromPAT_Zee.root");

 TTree* tree = (TTree*)f.Get("etoTauSC80/fitter_tree");

 TH1F hMass = TH1F("h1Mass","",60,60,120);

 tree->Draw("mass>>h1Mass");

 //TF1 *g = new TF1("g","gaus",60,120);
 //g->SetParameters(1000,90,10);

 //hMass.FillRandom("g", 10000);


 RooRealVar  mass("mass","mass (GeV/c^{2})",60,120);

 RooRealVar m1("m1","m1",91.2,89,93.);
 RooRealVar sigma("sigma","sigma",2.3,0.5,3.);
 RooRealVar alfa("alfa","alfa",7.2.,1.,10.);
 RooRealVar n("n","n",1,0.,10.);
 RooRealVar nsig("Nsig","#signal events",5000,0.,100000) ;
 RooRealVar nbkg("nbkg","#background events",80000,0.,1000000) ;
 RooCBShape signal("crystalBall","crystalBall",mass,m1,sigma,alfa,n);
 RooGaussian gauss("gauss","gauss",mass,m1,sigma) ;

 RooRealVar zero("zero","zero",0.0);
 RooRealVar resol("resol","resol",2,0,10);
 RooGaussian smear("smear","smear",mass,zero,resol) ;

 RooRealVar m0("m0","m0",91.2);
 RooRealVar gamma("gamma","gamma",2.5);
 RooGenericPdf breit("breit","breit",
		     "gamma*gamma*mass*mass/m0/m0/m0/m0/((mass*mass - m0*m0)*(mass*mass - m0*m0) + mass*mass*mass*mass*gamma*gamma*mass*mass/m0/m0/m0/m0)",RooArgSet(mass,m0,gamma));

 RooGenericPdf interference("interference","interference",
			    "(mass*mass - m0*m0)/((mass*mass - m0*m0)*(mass*mass - m0*m0) + mass*mass*mass*mass*gamma*gamma/m0/m0)",RooArgSet(mass,m0,gamma));

 RooGenericPdf dy("dy","dy","1./mass/mass",RooArgSet(mass));


 RooRealVar n0("n0","n0",10000,-100000,1000000);
 RooRealVar n1("n1","n1",0,-1000000,1000000);
 RooRealVar n2("n2","n2",0,-1000000,1000000);

 RooAddPdf sum("sum","breit+interference+dy",RooArgList(breit,interference,dy),RooArgList(n0,n1,n2)) ;

 //mass.setBins(100,"cache") ;

 //RooFFTConvPdf conv("conv","sum (X) smear",mass,sum,smear);
 

 RooDataHist *data = new RooDataHist("data","data",mass,Import(hMass));
 
 //gauss.chi2FitTo(*data,RooLinkedList()) ;
 //signal.fitTo(*data) ;
 conv.fitTo(*data) ;

 //sum.fitTo(*data) ;

 TCanvas c1("c1","",0.,0.,800,600);
 RooPlot* mesframe = mass.frame() ;
 mesframe->SetTitle("") ;

 data->plotOn(mesframe);

 sum.plotOn(mesframe);

 //signal.plotOn(mesframe,LineColor(kBlue));
 //gauss.plotOn(mesframe,LineColor(kRed));

 mesframe->Draw();


}
