#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

#include "ratioEfficiencyElec.C"
#include "ratioEfficiencyTest.C"


Double_t myFuncRatioElec15BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();

  Float_t xx = x[0];

  return ratioEffElec15->ratio(xx, true);

}

Double_t myFuncTurnOnElec15BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();

  Float_t xx = x[0];

  return ratioEffElec15->dataEfficiency(xx, true);

}

Double_t myFuncRatioElec15EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();

  Float_t xx = x[0];

  return ratioEffElec15->ratio(xx, false);

}

Double_t myFuncTurnOnElec15EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();

  Float_t xx = x[0];

  return ratioEffElec15->dataEfficiency(xx, false);

}

Double_t myFuncRatioEle18BL(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* fitEffEle18EB = (TF1*)corrections.Get("fitEffEle18EB");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return fitEffEle18EB->Eval(xx)/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncRatioEle18EC(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* fitEffEle18EE = (TF1*)corrections.Get("fitEffEle18EE");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return fitEffEle18EE->Eval(xx)/ratioEffElec15MC->mcEfficiency(xx, false);;
}

Double_t myFuncTurnOnEle18BL(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEle18EB = (TF1*)corrections.Get("fitEffEle18EB");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffEle18EB->Eval(xx);
}

Double_t myFuncTurnOnEle18EC(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEle18EE = (TF1*)corrections.Get("fitEffEle18EE");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffEle18EE->Eval(xx);
}

Double_t myFuncRatioEle20BL(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* fitEffEle20EB = (TF1*)corrections.Get("fitEffEle20EB");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return fitEffEle20EB->Eval(xx)/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncRatioEle20EC(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* fitEffEle20EE = (TF1*)corrections.Get("fitEffEle20EE");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return fitEffEle20EE->Eval(xx)/ratioEffElec15MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnEle20BL(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEle20EB = (TF1*)corrections.Get("fitEffEle20EB");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffEle20EB->Eval(xx);
}

Double_t myFuncTurnOnEle20EC(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEle20EE = (TF1*)corrections.Get("fitEffEle20EE");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffEle20EE->Eval(xx);
}

Double_t myFuncRatioEleAllBL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();
  TFile corrections("corrections.root");
  TF1* corrEffEle18EB = (TF1*)corrections.Get("fitEffEle18EB");
  TF1* corrEffEle20EB = (TF1*)corrections.Get("fitEffEle20EB");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  float weightElec15 = 1888.;
  float weightElec18 = 1945;
  float weightElec20 = 827.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return (ratioEffElec15->dataEfficiency(xx,true) * weightElec15 +
	  corrEffEle18EB->Eval(xx) * weightElec18 +
	  corrEffEle20EB->Eval(xx) * weightElec20
	  )/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEleAllBL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();
  TFile corrections("corrections.root");
  TF1* corrEffEle18EB = (TF1*)corrections.Get("fitEffEle18EB");
  TF1* corrEffEle20EB = (TF1*)corrections.Get("fitEffEle20EB");

  float weightElec15 = 1888.;
  float weightElec18 = 1945;
  float weightElec20 = 827.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return (ratioEffElec15->dataEfficiency(xx,true) * weightElec15 +
	  corrEffEle18EB->Eval(xx) * weightElec18 +
	  corrEffEle20EB->Eval(xx) * weightElec20
	  );
}

Double_t myFuncRatioEleAllEC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();
  TFile corrections("corrections.root");
  TF1* corrEffEle18EE = (TF1*)corrections.Get("fitEffEle18EE");
  TF1* corrEffEle20EE = (TF1*)corrections.Get("fitEffEle20EE");

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  float weightElec15 = 1888.;
  float weightElec18 = 1945;
  float weightElec20 = 827.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return (ratioEffElec15->dataEfficiency(xx,false) * weightElec15 +
	  corrEffEle18EE->Eval(xx) * weightElec18 +
	  corrEffEle20EE->Eval(xx) * weightElec20
	  )/ratioEffElec15MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnEleAllEC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();
  TFile corrections("corrections.root");
  TF1* corrEffEle18EE = (TF1*)corrections.Get("fitEffEle18EE");
  TF1* corrEffEle20EE = (TF1*)corrections.Get("fitEffEle20EE");

  float weightElec15 = 1888.;
  float weightElec18 = 1945;
  float weightElec20 = 827.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return (ratioEffElec15->dataEfficiency(xx,false) * weightElec15 +
	  corrEffEle18EE->Eval(xx) * weightElec18 +
	  corrEffEle20EE->Eval(xx) * weightElec20
	  );
}





Double_t myFuncRatioElecID(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEleID = (TF1*)corrections.Get("corrEffEleID");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffEleID->Eval(xx);
}

Double_t myFuncTurnOnElecID(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEleID = (TF1*)corrections.Get("fitEffEleID");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffEleID->Eval(xx);
}

Double_t myFuncRatioElecIso(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEleIso = (TF1*)corrections.Get("corrEffEleIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffEleIso->Eval(xx);
}

Double_t myFuncTurnOnElecIso(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffEleIso = (TF1*)corrections.Get("fitEffEleIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffEleIso->Eval(xx);
}


Double_t myFuncRatioMuId(Double_t* x, Double_t *par) {
  return 0.9967;
}

Double_t myFuncTurnOnMuId(Double_t* x, Double_t *par) {
  return 1.0;
}

Double_t myFuncRatioMuIso(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffMuIso = (TF1*)corrections.Get("corrEffMuIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffMuIso->Eval(xx);
}

Double_t myFuncTurnOnMuIso(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffMuIso = (TF1*)corrections.Get("fitEffMuIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return corrEffMuIso->Eval(xx);
}


Double_t myFuncRatioMu15(Double_t* x, Double_t *par) {
  return 0.968;
}

Double_t myFuncTurnOnMu15(Double_t* x, Double_t *par) {
  return 1.0;
}


Double_t myFuncRatioMu15L114(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffMuIso = (TF1*)corrections.Get("corrEffMu15_L114");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffMuIso->Eval(xx);
}

Double_t myFuncTurnOnMu15L114(Double_t* x, Double_t *par) {

  TFile corrections("corrections.root");
  TF1* corrEffMuIso = (TF1*)corrections.Get("fitEffMu15_L114");

  Float_t xx = x[0];
  if(xx>50) xx=50.;
  return corrEffMuIso->Eval(xx);
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////

// HLT PFTauLoose10 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose10MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];

  return ratioEffTauLoose10MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose15 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose15MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.5,1.50,6.4,31,0.804);
  
  Float_t xx = x[0];

  return ratioEffTauLoose15MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose20MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(20.91,3.43,7.94,29.8,0.810);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose15 / Run2011A    / mu+tau

Double_t myFuncRatioTauLoose15MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose15MuTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose15MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
   ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);
 
  Float_t xx = x[0];

  return ratioEffTauLoose15MuTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose15 / Run2011A    / e+tau

Double_t myFuncRatioTauLoose15ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose15ElecTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose15ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);  
 
  Float_t xx = x[0];

  return ratioEffTauLoose15ElecTauRunA->turnOn(xx);

}


/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011A    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA = new ratioEfficiencyTest(19.32,1.14,2.3,1.19,1.0);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose20MuTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose20MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA = new ratioEfficiencyTest(19.32,1.14,2.3,1.19,1.0);
 

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011A    / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
 
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0); 

  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauMedium20 / Run2011A / e+tau

Double_t myFuncRatioTauMedium20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34, 3.06,1.37,1.0);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauMedium20ElecTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauMedium20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34, 3.06,1.37,1.0);
  
  Float_t xx = x[0];

  return ratioEffTauMedium20ElecTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauTight20 / Run2011B / e+tau

Double_t myFuncRatioTauTight20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);

  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauTight20ElecTauRunA->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauTight20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
   
  Float_t xx = x[0];

  return ratioEffTauTight20ElecTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////


/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011B    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose20MuTauRunB->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose20MuTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunB->turnOn(xx);

}

/////////////////////////////////////////////////


// HLT PFTauLoose20 / Run2011B / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.79,2.498,6.37,1.10,0.921);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauLoose20ElecTauRunB->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauLoose20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.79,2.498,6.37,1.10,0.921);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunB->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauMedium20 / Run2011B / e+tau

Double_t myFuncRatioTauMedium20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47,2.4184,5.161,1.11,0.875);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauMedium20ElecTauRunB->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauMedium20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47,2.4184,5.161,1.11,0.875);
  
  Float_t xx = x[0];

  return ratioEffTauMedium20ElecTauRunB->turnOn(xx);

}






/////////////////////////////////////////////////

// HLT PFTauMedium20 / Run2011A+B / e+tau

Double_t myFuncRatioTauMedium20ElecTauRunAB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return (0.14*(ratioEffTauMedium20ElecTauRunA->turnOn(xx)) + 0.86*(ratioEffTauMedium20ElecTauRunB->turnOn(xx)) )/ effMC;

}

Double_t myFuncTurnOnTauMedium20ElecTauRunAB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  Float_t xx = x[0];

  return (0.14*(ratioEffTauMedium20ElecTauRunA->turnOn(xx)) + 0.86*(ratioEffTauMedium20ElecTauRunB->turnOn(xx)) );

}


/////////////////////////////////////////////////

// HLT PFTauTight20 / Run2011B / e+tau

Double_t myFuncRatioTauTight20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunB = new ratioEfficiencyTest(19.49,2.44,4.44,1.12,0.850);
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);

  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return ratioEffTauTight20ElecTauRunB->turnOn(xx)/ effMC;

}

Double_t myFuncTurnOnTauTight20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunB = new ratioEfficiencyTest(19.49,2.44,4.44,1.12,0.850);
  
  Float_t xx = x[0];

  return ratioEffTauTight20ElecTauRunB->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTau All / e+tau

Double_t myFuncRatioTauElecTauAll(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
 
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA  = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA  = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  float weightLoose15RunA  =  168.6 ;
  float weightLoose20RunA  =  934.0 ;
  float weightTight20RunA  =  785.7 ;
  float weightMedium20RunA =  240.0 ;
  float weightMedium20RunB = 2522.0 ;
  float totalRun2011 = weightLoose15RunA+weightLoose20RunA+weightTight20RunA+weightMedium20RunA+weightMedium20RunB;

  weightLoose15RunA  /= totalRun2011;
  weightLoose20RunA  /= totalRun2011;
  weightTight20RunA  /= totalRun2011;
  weightMedium20RunA /= totalRun2011;
  weightMedium20RunB /= totalRun2011;

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);

  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose10MC->turnOn(xx);

  return (ratioEffTauLoose15ElecTauRunA->turnOn(xx) *weightLoose15RunA  +
	  ratioEffTauLoose20ElecTauRunA->turnOn(xx) *weightLoose20RunA  +
	  ratioEffTauTight20ElecTauRunA->turnOn(xx) *weightTight20RunA  +
	  ratioEffTauMedium20ElecTauRunA->turnOn(xx)*weightMedium20RunA +
	  ratioEffTauMedium20ElecTauRunB->turnOn(xx)*weightMedium20RunB
	  )
    / effMC;

}


Double_t myFuncTurnOnTauElecTauAll(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA  = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA  = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  float weightLoose15RunA  =  168.6 ;
  float weightLoose20RunA  =  934.0 ;
  float weightTight20RunA  =  785.7 ;
  float weightMedium20RunA =  240.0 ;
  float weightMedium20RunB = 2522.0 ;
  float totalRun2011 = weightLoose15RunA+weightLoose20RunA+weightTight20RunA+weightMedium20RunA+weightMedium20RunB;

  weightLoose15RunA  /= totalRun2011;
  weightLoose20RunA  /= totalRun2011;
  weightTight20RunA  /= totalRun2011;
  weightMedium20RunA /= totalRun2011;
  weightMedium20RunB /= totalRun2011;
 
  Float_t xx = x[0];

  return (ratioEffTauLoose15ElecTauRunA->turnOn(xx) *weightLoose15RunA  +
	  ratioEffTauLoose20ElecTauRunA->turnOn(xx) *weightLoose20RunA  +
	  ratioEffTauTight20ElecTauRunA->turnOn(xx) *weightTight20RunA  +
	  ratioEffTauMedium20ElecTauRunA->turnOn(xx)*weightMedium20RunA +
	  ratioEffTauMedium20ElecTauRunB->turnOn(xx)*weightMedium20RunB
	  );

}


/////////////////////////////////////////////////






void makeFile(){

  TFile* fout = new TFile("llrCorrections.root","RECREATE");

  TF1 *ratioElec15BL        = new TF1("ratioElec15BL",           myFuncRatioElec15BL,       15,400,0);
  TF1 *turnOnElec15BL       = new TF1("turnOnElec15BL",          myFuncTurnOnElec15BL,      15,400,0);
  TF1 *ratioElec15EC        = new TF1("ratioElec15EC",           myFuncRatioElec15EC,       15,400,0);
  TF1 *turnOnElec15EC       = new TF1("turnOnElec15EC",          myFuncTurnOnElec15EC,      15,400,0);

  TF1 *ratioElec18BL        = new TF1("ratioElec18BL",           myFuncRatioEle18BL,        15,400,0);
  TF1 *turnOnElec18BL       = new TF1("turnOnElec18BL",          myFuncTurnOnEle18BL,       15,400,0);
  TF1 *ratioElec18EC        = new TF1("ratioElec18EC",           myFuncRatioEle18EC,        15,400,0);
  TF1 *turnOnElec18EC       = new TF1("turnOnElec18EC",          myFuncTurnOnEle18EC,       15,400,0);

  TF1 *ratioElec20BL        = new TF1("ratioElec20BL",           myFuncRatioEle20BL,        15,400,0);
  TF1 *turnOnElec20BL       = new TF1("turnOnElec20BL",          myFuncTurnOnEle20BL,       15,400,0);
  TF1 *ratioElec20EC        = new TF1("ratioElec20EC",           myFuncRatioEle20EC,        15,400,0);
  TF1 *turnOnElec20EC       = new TF1("turnOnElec20EC",          myFuncTurnOnEle20EC,       15,400,0);

  TF1 *ratioElecAllBL       = new TF1("ratioElecAllBL",          myFuncRatioEleAllBL,       15,400,0);
  TF1 *turnOnElecAllBL      = new TF1("turnOnElecAllBL",         myFuncTurnOnEleAllBL,      15,400,0);
  TF1 *ratioElecAllEC       = new TF1("ratioElecAllEC",          myFuncRatioEleAllEC,       15,400,0);
  TF1 *turnOnElecAllEC      = new TF1("turnOnElecAllEC",         myFuncTurnOnEleAllEC,      15,400,0);

  TF1 *ratioElecID          = new TF1("ratioElecID",             myFuncRatioElecID ,        15,400,0);
  TF1 *turnOnElecID         = new TF1("turnOnElecID",            myFuncTurnOnElecID ,       15,400,0);
  TF1 *ratioElecIso         = new TF1("ratioElecIso",            myFuncRatioElecIso ,       15,400,0);
  TF1 *turnOnElecIso        = new TF1("turnOnElecIso",           myFuncTurnOnElecIso ,      15,400,0);

  TF1 *ratioMu15            = new TF1("ratioMu15",               myFuncRatioMu15,           14,400,0);
  TF1 *turnOnMu15           = new TF1("turnOnMu15",              myFuncTurnOnMu15,          14,400,0);

  TF1 *ratioMu15L114        = new TF1("ratioMu15L114",           myFuncRatioMu15L114,       14,400,0);
  TF1 *turnOnMu15L114       = new TF1("turnOnMu15L114",          myFuncTurnOnMu15L114,      14,400,0);
 
  TF1 *ratioMuId            = new TF1("ratioMuId",               myFuncRatioMuId ,          14,400,0);
  TF1 *turnOnMuId           = new TF1("turnOnMuId",              myFuncTurnOnMuId ,         14,400,0);
  TF1 *ratioMuIso           = new TF1("ratioMuIso",              myFuncRatioMuIso ,         14,400,0);
  TF1 *turnOnMuIso          = new TF1("turnOnMuIso",             myFuncTurnOnMuIso ,        14,400,0);

  TF1 *turnOnTauLoose10MuTauMC      = new TF1("turnOnTauLoose10MuTauMC",       myFuncTurnOnTauLoose10MuTauMC,      18,400,0);
  TF1 *turnOnTauLoose15MuTauMC      = new TF1("turnOnTauLoose15MuTauMC",       myFuncTurnOnTauLoose15MuTauMC,      18,400,0);
  TF1 *turnOnTauLoose20MuTauMC      = new TF1("turnOnTauLoose20MuTauMC",       myFuncTurnOnTauLoose20MuTauMC,      18,400,0);

  TF1 *ratioTauLoose15MuTauRunA      = new TF1("ratioTauLoose15MuTauRunA",     myFuncRatioTauLoose15MuTauRunA,     18,400,0);
  TF1 *turnOnTauLoose15MuTauRunA     = new TF1("turnOnTauLoose15MuTauRunA",    myFuncTurnOnTauLoose15MuTauRunA,    18,400,0);
  TF1 *ratioTauLoose15ElecTauRunA    = new TF1("ratioTauLoose15ElecTauRunA",   myFuncRatioTauLoose15ElecTauRunA ,  18,400,0);
  TF1 *turnOnTauLoose15ElecTauRunA   = new TF1("turnOnTauLoose15ElecTauRunA",  myFuncTurnOnTauLoose15ElecTauRunA,  18,400,0);
  TF1 *ratioTauLoose20MuTauRunA      = new TF1("ratioTauLoose20MuTauRunA",     myFuncRatioTauLoose20MuTauRunA,     18,400,0);
  TF1 *turnOnTauLoose20MuTauRunA     = new TF1("turnOnTauLoose20MuTauRunA",    myFuncTurnOnTauLoose20MuTauRunA,    18,400,0);
  TF1 *ratioTauLoose20ElecTauRunA    = new TF1("ratioTauLoose20ElecTauRunA",   myFuncRatioTauLoose20ElecTauRunA ,  18,400,0);
  TF1 *turnOnTauLoose20ElecTauRunA   = new TF1("turnOnTauLoose20ElecTauRunA",  myFuncTurnOnTauLoose20ElecTauRunA,  18,400,0);
  TF1 *ratioTauMedium20ElecTauRunA   = new TF1("ratioTauMedium20ElecTauRunA",  myFuncRatioTauMedium20ElecTauRunA , 18,400,0);
  TF1 *turnOnTauMedium20ElecTauRunA  = new TF1("turnOnTauMedium20ElecTauRunA", myFuncTurnOnTauMedium20ElecTauRunA, 18,400,0);
  TF1 *ratioTauTight20ElecTauRunA    = new TF1("ratioTauTight20ElecTauRunA",   myFuncRatioTauTight20ElecTauRunA ,  18,400,0);
  TF1 *turnOnTauTight20ElecTauRunA   = new TF1("turnOnTauTight20ElecTauRunA",  myFuncTurnOnTauTight20ElecTauRunA,  18,400,0);
  TF1 *ratioTauLoose20MuTauRunB      = new TF1("ratioTauLoose20MuTauRunB",     myFuncRatioTauLoose20MuTauRunB,     18,400,0);
  TF1 *turnOnTauLoose20MuTauRunB     = new TF1("turnOnTauLoose20MuTauRunB",    myFuncTurnOnTauLoose20MuTauRunB,    18,400,0);
  TF1 *ratioTauLoose20ElecTauRunB    = new TF1("ratioTauLoose20ElecTauRunB",   myFuncRatioTauLoose20ElecTauRunB ,  18,400,0);
  TF1 *turnOnTauLoose20ElecTauRunB   = new TF1("turnOnTauLoose20ElecTauRunB",  myFuncTurnOnTauLoose20ElecTauRunB,  18,400,0);
  TF1 *ratioTauMedium20ElecTauRunB   = new TF1("ratioTauMedium20ElecTauRunB",  myFuncRatioTauMedium20ElecTauRunB , 18,400,0);
  TF1 *turnOnTauMedium20ElecTauRunB  = new TF1("turnOnTauMedium20ElecTauRunB", myFuncTurnOnTauMedium20ElecTauRunB, 18,400,0);
  TF1 *ratioTauTight20ElecTauRunB    = new TF1("ratioTauTight20ElecTauRunB",   myFuncRatioTauTight20ElecTauRunB ,  18,400,0);
  TF1 *turnOnTauTight20ElecTauRunB   = new TF1("turnOnTauTight20ElecTauRunB",  myFuncTurnOnTauTight20ElecTauRunB,  18,400,0);
  TF1 *ratioTauMedium20ElecTauRunAB  = new TF1("ratioTauMedium20ElecTauRunAB", myFuncRatioTauMedium20ElecTauRunAB ,18,400,0);
  TF1 *turnOnTauMedium20ElecTauRunAB = new TF1("turnOnTauMedium20ElecTauRunAB",myFuncTurnOnTauMedium20ElecTauRunAB,18,400,0);
 
  TF1 *ratioTauElecTauAll    = new TF1("ratioTauElecTauAll", myFuncRatioTauElecTauAll ,18,400,0);
  TF1 *turnOnTauElecTauAll   = new TF1("turnOnTauElecTauAll",myFuncTurnOnTauElecTauAll,18,400,0);
 

  fout->cd();
  
  ratioElec15BL->Write();
  turnOnElec15BL->Write();
  ratioElec15EC->Write();
  turnOnElec15EC->Write();
  ratioElec18BL->Write();
  turnOnElec18BL->Write();
  ratioElec18EC->Write();
  turnOnElec18EC->Write();
  ratioElec20BL->Write();
  turnOnElec20BL->Write();
  ratioElec20EC->Write();
  turnOnElec20EC->Write();
  ratioElecAllBL->Write();
  turnOnElecAllBL->Write();
  ratioElecAllEC->Write();
  turnOnElecAllEC->Write();
  ratioElecID->Write();
  turnOnElecID->Write();
  ratioElecIso->Write();
  turnOnElecIso->Write();
  ratioMu15->Write();
  turnOnMu15->Write();
  ratioMu15L114->Write();
  turnOnMu15L114->Write();
  ratioMuId->Write();
  turnOnMuId->Write();
  ratioMuIso->Write();
  turnOnMuIso->Write();
  turnOnTauLoose10MuTauMC->Write();
  turnOnTauLoose15MuTauMC->Write();
  turnOnTauLoose20MuTauMC->Write();
  ratioTauLoose15MuTauRunA->Write();
  turnOnTauLoose15MuTauRunA->Write();
  ratioTauLoose15ElecTauRunA->Write();
  turnOnTauLoose15ElecTauRunA->Write();
  ratioTauLoose20MuTauRunA->Write();
  turnOnTauLoose20MuTauRunA->Write();
  ratioTauLoose20ElecTauRunA->Write();
  turnOnTauLoose20ElecTauRunA->Write();
  ratioTauMedium20ElecTauRunA->Write();
  turnOnTauMedium20ElecTauRunA->Write();
  ratioTauTight20ElecTauRunA->Write();
  turnOnTauTight20ElecTauRunA->Write();
  ratioTauLoose20MuTauRunB->Write();
  turnOnTauLoose20MuTauRunB->Write();
  ratioTauLoose20ElecTauRunB->Write();
  turnOnTauLoose20ElecTauRunB->Write();
  ratioTauMedium20ElecTauRunB->Write();
  turnOnTauMedium20ElecTauRunB->Write();
  ratioTauTight20ElecTauRunB->Write();
  turnOnTauTight20ElecTauRunB->Write();
  ratioTauMedium20ElecTauRunAB->Write();
  turnOnTauMedium20ElecTauRunAB->Write();
  ratioTauElecTauAll->Write();
  turnOnTauElecTauAll->Write();

  fout->Write();
  fout->Close();


}
