#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

#include "ratioEfficiencyElec.C"
#include "ratioEfficiencyTest.C"


Double_t myFuncRatioElec15BL(Double_t* x, Double_t *par) {

  //gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();

  Float_t xx = x[0];

  return ratioEffElec15->ratio(xx, true);

}

Double_t myFuncTurnOnElec15BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15 = new ratioEfficiencyElec();
  //ratioEfficiencyTest* ratioEffElec15BL = new ratioEfficiencyTest(14.87,0.3112 , 0.2210, 1.877, 0.987 );

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
  //ratioEfficiencyTest* ratioEffElec15EC = new ratioEfficiencyTest(15.66, 0.7591 , 0.4775, 2.0215, 0.99881 );
 
  Float_t xx = x[0];

  return ratioEffElec15->dataEfficiency(xx, false);

}

Double_t myFuncRatioEle18BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle18EB = new ratioEfficiencyTest(18.99, 0.0010589 , 1.2828e-05, 1.6480, 0.99604 );

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return fitEffEle18EB->turnOn(xx)/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEle18BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle18EB = new ratioEfficiencyTest(18.99,0.001058, 1.282e-05,1.648, 0.9960);

  Float_t xx = x[0];
  return fitEffEle18EB->turnOn(xx);
}


Double_t myFuncRatioEle18EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle18EC = new ratioEfficiencyTest(18.84, 1.5717e-03 , 1.6403e-06, 2.866, 0.96996 );

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return fitEffEle18EC->turnOn(xx)/ratioEffElec15MC->mcEfficiency(xx, false);;
}


Double_t myFuncTurnOnEle18EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle18EC = new ratioEfficiencyTest(18.84, 1.5717e-03 , 1.6403e-06, 2.866, 0.96996 );

  Float_t xx = x[0];
  return fitEffEle18EC->turnOn(xx);
}


Double_t myFuncRatioEle20BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EB = new ratioEfficiencyTest(20.55,0.6837,0.85557,1.459166, 1.0395732 );

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return fitEffEle20EB->turnOn(xx)/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEle20BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EB = new ratioEfficiencyTest(20.55,0.6837,0.85557,1.4591,1.0395 );

  Float_t xx = x[0];
  return fitEffEle20EB->turnOn(xx);
}

Double_t myFuncRatioEle20EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest( 23.63, -1.6077, 1.7209, 1.413, 1.1396);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return fitEffEle20EC->turnOn(xx)/ratioEffElec15MC->mcEfficiency(xx, false);
}


Double_t myFuncTurnOnEle20EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest( 23.63, -1.6077, 1.7209, 1.413, 1.1396);

  Float_t xx = x[0];
  return fitEffEle20EC->turnOn(xx);
}

Double_t myFuncRatioEleAllBL(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle15BL = new ratioEfficiencyTest(14.87,0.3112 , 0.2210, 1.877, 0.987 );
  ratioEfficiencyTest* fitEffEle18BL = new ratioEfficiencyTest(18.99,0.0010589 ,1.2828e-05, 1.6480, 0.99604 );
  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.55,0.6837,0.85557,1.4591,1.0395 );

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  float weightElec15 = 1871.;
  float weightElec18 = 1713.;
  float weightElec20 =  895.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];

  return (fitEffEle15BL->turnOn(xx) * weightElec15 +
	  fitEffEle18BL->turnOn(xx) * weightElec18 +
	  fitEffEle20BL->turnOn(xx) * weightElec20
	  )/ratioEffElec15MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEleAllBL(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle15BL = new ratioEfficiencyTest(14.87,0.3112 , 0.2210, 1.877, 0.987 );
  ratioEfficiencyTest* fitEffEle18BL = new ratioEfficiencyTest(18.99,0.0010589 , 1.2828e-05, 1.6480, 0.99604 );
  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.55,0.6837,0.85557,1.4591,1.0395 );

  float weightElec15 = 1871.;
  float weightElec18 = 1713.;
  float weightElec20 =  895.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];

  return (fitEffEle15BL->turnOn(xx) * weightElec15 +
	  fitEffEle18BL->turnOn(xx) * weightElec18 +
	  fitEffEle20BL->turnOn(xx) * weightElec20
	  );
 
}



Double_t myFuncRatioEleAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle15EC = new ratioEfficiencyTest(15.66, 0.7591 , 0.4775, 2.0215, 0.99881 );
  ratioEfficiencyTest* fitEffEle18EC = new ratioEfficiencyTest(18.84, 1.5717e-03 , 1.6403e-06, 2.866, 0.96996 );
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest( 23.63, -1.6077, 1.7209, 1.413, 1.1396);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec15MC = new ratioEfficiencyElec();

  float weightElec15 = 1871.;
  float weightElec18 = 1713.;
  float weightElec20 =  895.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];

  return (fitEffEle15EC->turnOn(xx) * weightElec15 +
	  fitEffEle18EC->turnOn(xx) * weightElec18 +
	  fitEffEle20EC->turnOn(xx) * weightElec20
	  )/ratioEffElec15MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnEleAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle15EC = new ratioEfficiencyTest( 15.66, 0.7591 , 0.4775, 2.0215, 0.99881 );
  ratioEfficiencyTest* fitEffEle18EC = new ratioEfficiencyTest( 18.84, 1.5717e-03 , 1.6403e-06, 2.866, 0.96996 );
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest( 23.63,-1.6077, 1.7209, 1.413, 1.1396);

  float weightElec15 = 1871.;
  float weightElec18 = 1713.;
  float weightElec20 =  895.;

  float total = weightElec15+weightElec18+weightElec20;
  weightElec15/=total;
  weightElec18/=total;
  weightElec20/=total;

  Float_t xx = x[0];

  return (fitEffEle15EC->turnOn(xx) * weightElec15 +
	  fitEffEle18EC->turnOn(xx) * weightElec18 +
	  fitEffEle20EC->turnOn(xx) * weightElec20
	  );
 
}






///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////




Double_t myFuncRatioElecID(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffEleID = (TF1*)corrections.Get("corrEffEleID");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffEleID->Eval(xx);
}

Double_t myFuncTurnOnElecID(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffEleID = (TF1*)corrections.Get("fitEffEleID");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffEleID->Eval(xx);
}

Double_t myFuncRatioElecIso(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffEleIso = (TF1*)corrections.Get("corrEffEleIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffEleIso->Eval(xx);
}

Double_t myFuncTurnOnElecIso(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffEleIso = (TF1*)corrections.Get("fitEffEleIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffEleIso->Eval(xx);
}




///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////



Double_t myFuncRatioMuId(Double_t* x, Double_t *par) {
  return 0.9967;
}

Double_t myFuncTurnOnMuId(Double_t* x, Double_t *par) {
  return 1.0;
}

Double_t myFuncRatioMuIso(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffMuIso = (TF1*)corrections.Get("corrEffMuIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffMuIso->Eval(xx);
}

Double_t myFuncTurnOnMuIso(Double_t* x, Double_t *par) {

  //TFile corrections("corrections.root");
  //TF1* corrEffMuIso = (TF1*)corrections.Get("fitEffMuIso");

  Float_t xx = x[0];
  if(xx>100) xx=100.;
  return 1.0;//corrEffMuIso->Eval(xx);
}






///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////











Double_t myFuncTurnOnMu12MCBL(Double_t* x, Double_t *par) {

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.93;
}

Double_t myFuncTurnOnMu12MCEC(Double_t* x, Double_t *par) {

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.84;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMu15L110BL(Double_t* x, Double_t *par) {

  TF1* turnOnMu12MCBL = new TF1("turnOnMu12MCBL",myFuncTurnOnMu12MCBL,0,400,0);

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.917/turnOnMu12MCBL->Eval(xx);

}

Double_t myFuncTurnOnMu15L110BL(Double_t* x, Double_t *par) {

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.917;
}

Double_t myFuncRatioMu15L110EC(Double_t* x, Double_t *par) {

  TF1* turnOnMu12MCEC = new TF1("turnOnMu12MCEC",myFuncTurnOnMu12MCEC,0,400,0);

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.836/turnOnMu12MCEC->Eval(xx);

}

Double_t myFuncTurnOnMu15L110EC(Double_t* x, Double_t *par) {

  Float_t xx = x[0];
  if(xx<12) 
    return 0;
  else 
    return 0.836;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMu15L114BL(Double_t* x, Double_t *par) {

  TF1* turnOnMu12MCBL = new TF1("turnOnMu12MCBL",myFuncTurnOnMu12MCBL,0,400,0);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffMuL114BL = new ratioEfficiencyTest(15.06, 0.55278, 1.3423, 1.002975, 3.3676);

  Float_t xx = x[0];
  return ratioEffMuL114BL->turnOn(xx)/turnOnMu12MCBL->Eval(xx);
}

Double_t myFuncTurnOnMu15L114BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffMuL114BL = new ratioEfficiencyTest(15.06, 0.55278, 1.3423, 1.002975, 3.3676);


  Float_t xx = x[0];
  return ratioEffMuL114BL->turnOn(xx);
}

Double_t myFuncRatioMu15L114EC(Double_t* x, Double_t *par) {

  TF1* turnOnMu12MCEC = new TF1("turnOnMu12MCEC",myFuncTurnOnMu12MCEC,0,400,0);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffMuL114EC = new ratioEfficiencyTest(15.32, 0.866114, 1.25008, 1.63711, 0.84490);

  Float_t xx = x[0];
  return ratioEffMuL114EC->turnOn(xx)/turnOnMu12MCEC->Eval(xx);
}

Double_t myFuncTurnOnMu15L114EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffMuL114EC = new ratioEfficiencyTest(15.32, 0.866114, 1.25008, 1.63711, 0.84490);

  Float_t xx = x[0];
  return ratioEffMuL114EC->turnOn(xx);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////


Double_t myFuncRatioMuAllBL(Double_t* x, Double_t *par) {

  TF1* ratioMu15L110BL = new TF1("ratioMu15L110BL", myFuncRatioMu15L110BL,0,400,0);
  TF1* ratioMu15L114BL = new TF1("ratioMu15L114BL", myFuncRatioMu15L114BL,0,400,0);
  
  float weightIsoMu15L110 = 2200.;
  float weightIsoMu15L114 = 2400.;

  float total = weightIsoMu15L110+weightIsoMu15L114;
  weightIsoMu15L110 /= total;
  weightIsoMu15L114 /= total;

  Float_t xx = x[0];

  return (ratioMu15L110BL->Eval(xx) * weightIsoMu15L110 +
	  ratioMu15L114BL->Eval(xx) * weightIsoMu15L114 );

}

Double_t myFuncTurnOnMuAllBL(Double_t* x, Double_t *par) {

  TF1* turnOnMu15L110BL = new TF1("turnOnMu15L110BL", myFuncTurnOnMu15L110BL,0,400,0);
  TF1* turnOnMu15L114BL = new TF1("turnOnMu15L114BL", myFuncTurnOnMu15L114BL,0,400,0);

  float weightIsoMu15L110 = 2200.;
  float weightIsoMu15L114 = 2400.;

  float total = weightIsoMu15L110+weightIsoMu15L114;
  weightIsoMu15L110 /=total;
  weightIsoMu15L114 /=total;

  Float_t xx = x[0];

  return (turnOnMu15L110BL->Eval(xx) * weightIsoMu15L110 +
	  turnOnMu15L114BL->Eval(xx) * weightIsoMu15L114 );
}


Double_t myFuncRatioMuAllEC(Double_t* x, Double_t *par) {

  TF1* ratioMu15L110EC = new TF1("ratioMu15L110EC", myFuncRatioMu15L110EC,0,400,0);
  TF1* ratioMu15L114EC = new TF1("ratioMu15L114EC", myFuncRatioMu15L114EC,0,400,0);
  
  float weightIsoMu15L110 = 2200.;
  float weightIsoMu15L114 = 2400.;

  float total = weightIsoMu15L110+weightIsoMu15L114;
  weightIsoMu15L110 /= total;
  weightIsoMu15L114 /= total;

  Float_t xx = x[0];

  return (ratioMu15L110EC->Eval(xx) * weightIsoMu15L110 +
	  ratioMu15L114EC->Eval(xx) * weightIsoMu15L114 );

}

Double_t myFuncTurnOnMuAllEC(Double_t* x, Double_t *par) {

  TF1* turnOnMu15L110EC = new TF1("turnOnMu15L110EC", myFuncTurnOnMu15L110EC,0,400,0);
  TF1* turnOnMu15L114EC = new TF1("turnOnMu15L114EC", myFuncTurnOnMu15L114EC,0,400,0);

  float weightIsoMu15L110 = 2200.;
  float weightIsoMu15L114 = 2400.;

  float total = weightIsoMu15L110+weightIsoMu15L114;
  weightIsoMu15L110 /=total;
  weightIsoMu15L114 /=total;

  Float_t xx = x[0];

  return (turnOnMu15L110EC->Eval(xx) * weightIsoMu15L110 +
	  turnOnMu15L114EC->Eval(xx) * weightIsoMu15L114 );
}






/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////



// HLT PFTauLoose10 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose10MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(13.06,2.60,5.03,1.11,1.0);
  
  Float_t xx = x[0];

  return ratioEffTauLoose10MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose15 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose15MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.5,1.50,6.4,31,0.804);
  ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.5,1.76,3.37,1.11,1.00);
  
  Float_t xx = x[0];

  return ratioEffTauLoose15MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / MC Summer11 / mu+tau

Double_t myFuncTurnOnTauLoose20MuTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(20.91,3.43,7.94,29.8,0.810);
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(20.89,3.25,7.99,34.5,0.811);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////

// HLT PFTauLoose15 / Run2011A    / mu+tau

Double_t myFuncRatioTauLoose15MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);  
  TF1* turnOnTauLoose10MuTauMC = new TF1("turnOnTauLoose10MuTauMC", myFuncTurnOnTauLoose10MuTauMC,0,400,0);

  Float_t xx = x[0];

  return ratioEffTauLoose15MuTauRunA->turnOn(xx)/turnOnTauLoose10MuTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose15MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);
 
  Float_t xx = x[0];

  return ratioEffTauLoose15MuTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011A    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA = new ratioEfficiencyTest(19.32,1.14,2.3,1.19,1.0);

  TF1* turnOnTauLoose10MuTauMC = new TF1("turnOnTauLoose10MuTauMC", myFuncTurnOnTauLoose10MuTauMC,0,400,0);

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunA->turnOn(xx)/ turnOnTauLoose10MuTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose20MuTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA = new ratioEfficiencyTest(19.32,1.14,2.3,1.19,1.0);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA =  new ratioEfficiencyTest(19.72,0.84,1.17,1.007,9.35 );
    //new ratioEfficiencyTest(19.72,0.84,1.17,1.01,9.35);
  //

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011B    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);
  
  TF1* turnOnTauLoose10MuTauMC = new TF1("turnOnTauLoose10MuTauMC", myFuncTurnOnTauLoose10MuTauMC,0,400,0);

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunB->turnOn(xx)/ turnOnTauLoose10MuTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose20MuTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunB->turnOn(xx);

}

/////////////////////////////////////////////////


Double_t myFuncRatioTauMuTauAll(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
 
  //IVO:
  //ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);  
  //ratioEfficiencyTest* ratioEffTauLoose20MuTauRunA = new ratioEfficiencyTest(19.32,1.14,2.3,1.19,1.0);
  //ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);

  //NOTE:
  ratioEfficiencyTest* ratioEffTauLoose10MuTau = new ratioEfficiencyTest(16.785,-0.6938,4.5728,94127.12,0.8831);
  ratioEfficiencyTest* ratioEffTauLoose15MuTau = new ratioEfficiencyTest(14.67,0.4082,0.5519,1.4477,0.9613);
  ratioEfficiencyTest* ratioEffTauLoose20MuTau = new ratioEfficiencyTest(19.19,-1.36,2.827,1.02721,1.5086);
 

  TF1* turnOnTauLoose10MuTauMC = new TF1("turnOnTauLoose10MuTauMC", myFuncTurnOnTauLoose10MuTauMC,0,400,0);

  float weightLoose10 =  169.;
  float weightLoose15 = 1968.;
  float weightLoose20 = 2460.;
  
  float total = weightLoose10+weightLoose15+weightLoose20;
  weightLoose10 /= total;
  weightLoose15 /= total;
  weightLoose20 /= total;
 
  Float_t xx = x[0];

  return (ratioEffTauLoose10MuTau->turnOn(xx) * weightLoose10 + 
	  ratioEffTauLoose15MuTau->turnOn(xx) * weightLoose15 +
	  ratioEffTauLoose20MuTau->turnOn(xx) * weightLoose20
	  )/turnOnTauLoose10MuTauMC->Eval(xx) ;

}

Double_t myFuncTurnOnTauMuTauAll(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose15MuTauRunA = new ratioEfficiencyTest(14.274,1.048,1.706,1.336,1.0);
  //ratioEfficiencyTest* ratioEffTauLoose20MuTauRunB = new ratioEfficiencyTest(18.6748,2.4703,6.59219,1.10415,0.908873);
  //ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.96,1.41,3.31,1.10,0.868);
  
  //NOTE:
  ratioEfficiencyTest* ratioEffTauLoose10MuTau = new ratioEfficiencyTest(16.785,-0.6938,4.5728,94127.12,0.8831);
  ratioEfficiencyTest* ratioEffTauLoose15MuTau = new ratioEfficiencyTest(14.67,0.4082,0.5519,1.4477,0.9613);
  ratioEfficiencyTest* ratioEffTauLoose20MuTau = new ratioEfficiencyTest(19.19,-1.36,2.827,1.02721,1.5086);
  
  float weightLoose10 =  169.;
  float weightLoose15 = 1968.;
  float weightLoose20 = 2460.;
  
  float total = weightLoose10+weightLoose15+weightLoose20;
  weightLoose10 /= total;
  weightLoose15 /= total;
  weightLoose20 /= total;
  
  Float_t xx = x[0];
  
  return (ratioEffTauLoose10MuTau->turnOn(xx) * weightLoose10 + 
	  ratioEffTauLoose15MuTau->turnOn(xx) * weightLoose15 +
	  ratioEffTauLoose20MuTau->turnOn(xx) * weightLoose20
	  );
  
 

}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////




// HLT PFTauLoose10 / MC Summer11 / elec+tau

Double_t myFuncTurnOnTauLoose10ElecTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose10MC = new ratioEfficiencyTest(12.92, 2.41,4.09,1.33,0.90);
  
  Float_t xx = x[0];

  return ratioEffTauLoose10MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose15 / MC Summer11 / elec+tau

Double_t myFuncTurnOnTauLoose15ElecTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.57,1.77,3.34,1.11,1.00);
  
  Float_t xx = x[0];

  return ratioEffTauLoose15MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / MC Summer11 / elec+tau

Double_t myFuncTurnOnTauLoose20ElecTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(20.99,3.13,7.99,34.9,0.813);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////

// HLT PFTauLoose15 / Run2011A    / e+tau

Double_t myFuncRatioTauLoose15ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);

  Float_t xx = x[0];

  return ratioEffTauLoose15ElecTauRunA->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose15ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);  
 
  Float_t xx = x[0];

  return ratioEffTauLoose15ElecTauRunA->turnOn(xx);

}



/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2011A    / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);

  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

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
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);
  
  Float_t xx = x[0];

  return ratioEffTauMedium20ElecTauRunA->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

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
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);
  
  Float_t xx = x[0];

  return ratioEffTauTight20ElecTauRunA->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauTight20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
   
  Float_t xx = x[0];

  return ratioEffTauTight20ElecTauRunA->turnOn(xx);

}

/////////////////////////////////////////////////





// HLT PFTauLoose20 / Run2011B / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.79,2.498,6.37,1.10,0.921);
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunB->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

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
  
  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);
 
  Float_t xx = x[0];

  return ratioEffTauMedium20ElecTauRunB->turnOn(xx)/ turnOnTauLoose15ElecTauMC->Eval(xx);

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
  //IVO:
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);

  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);

  Float_t xx = x[0];

  return (0.14*(ratioEffTauMedium20ElecTauRunA->turnOn(xx)) + 0.86*(ratioEffTauMedium20ElecTauRunB->turnOn(xx)) )/turnOnTauLoose15ElecTauMC->Eval(xx);

 

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
  ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.57,1.77,3.34,1.11,1.00);

  
  Float_t xx = x[0];
  float effMC = ratioEffTauLoose15MC->turnOn(xx);

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
 
  ratioEfficiencyTest* ratioEffTauLoose20ElecTau  = new ratioEfficiencyTest(19.63,-0.986,1.94,1.023979,1.911);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTau = new ratioEfficiencyTest(19.35,0.37,0.1582,3.311,0.763);
  ratioEfficiencyTest* ratioEffTauTight20ElecTau  = new ratioEfficiencyTest(19.72,0.844,1.167,1.007466,9.351);

  TF1* turnOnTauLoose15ElecTauMC = new TF1("turnOnTauLoose15ElecTauMC", myFuncTurnOnTauLoose10ElecTauMC,0,400,0);

  float weightLoose20  =  1103.;
  float weightMedium20 =  2772.;
  float weightTight20  =   786.;

  float total = weightLoose20+weightMedium20+weightTight20;
  weightLoose20  /= total;
  weightMedium20 /= total;
  weightTight20  /= total;

  Float_t xx = x[0];

  return (ratioEffTauLoose20ElecTau->turnOn(xx)  * weightLoose20 +
	  ratioEffTauMedium20ElecTau->turnOn(xx) * weightMedium20 +
	  ratioEffTauTight20ElecTau->turnOn(xx)  * weightTight20
	  )/turnOnTauLoose15ElecTauMC->Eval(xx);

  //ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA  = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  //ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  //ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA  = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
  //ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  //ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  //float weightLoose15RunA  =  168.6 ;
  //float weightLoose20RunA  =  934.0 ;
  //float weightTight20RunA  =  785.7 ;
  //float weightMedium20RunA =  240.0 ;
  //float weightMedium20RunB = 2522.0 ;
  //float totalRun2011 = weightLoose15RunA+weightLoose20RunA+weightTight20RunA+weightMedium20RunA+weightMedium20RunB;

  //weightLoose15RunA  /= totalRun2011;
  //weightLoose20RunA  /= totalRun2011;
  //weightTight20RunA  /= totalRun2011;
  //weightMedium20RunA /= totalRun2011;
  //weightMedium20RunB /= totalRun2011;

  //gSystem->Load("ratioEfficiencyTest_C.so");
  //ratioEfficiencyTest* ratioEffTauLoose15MC = new ratioEfficiencyTest(15.57,1.77,3.34,1.11,1.00);

  
  //Float_t xx = x[0];
  //float effMC = ratioEffTauLoose15MC->turnOn(xx);

  //return (ratioEffTauLoose15ElecTauRunA->turnOn(xx) *weightLoose15RunA  +
  //	  ratioEffTauLoose20ElecTauRunA->turnOn(xx) *weightLoose20RunA  +
  //  ratioEffTauTight20ElecTauRunA->turnOn(xx) *weightTight20RunA  +
  // ratioEffTauMedium20ElecTauRunA->turnOn(xx)*weightMedium20RunA +
  //  ratioEffTauMedium20ElecTauRunB->turnOn(xx)*weightMedium20RunB
  //  )
  /// effMC;

}


Double_t myFuncTurnOnTauElecTauAll(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20ElecTau  = new ratioEfficiencyTest(19.63,-0.986,1.94,1.023979,1.911);
  ratioEfficiencyTest* ratioEffTauMedium20ElecTau = new ratioEfficiencyTest(19.35,0.37,0.1582,3.311,0.763);
  ratioEfficiencyTest* ratioEffTauTight20ElecTau  = new ratioEfficiencyTest(19.72,0.844,1.167,1.007466,9.351);

  float weightLoose20  =  1103.;
  float weightMedium20 =  2772.;
  float weightTight20  =   786.;

  float total = weightLoose20+weightMedium20+weightTight20;
  weightLoose20  /= total;
  weightMedium20 /= total;
  weightTight20  /= total;

  Float_t xx = x[0];

  return (ratioEffTauLoose20ElecTau->turnOn(xx)  * weightLoose20 +
	  ratioEffTauMedium20ElecTau->turnOn(xx) * weightMedium20 +
	  ratioEffTauTight20ElecTau->turnOn(xx)  * weightTight20
	  );

  //ratioEfficiencyTest* ratioEffTauLoose15ElecTauRunA  = new ratioEfficiencyTest(14.32,1.11,1.84,1.34,1.0);
  //ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(19.38,1.05,1.99,1.25,1.0);
  //ratioEfficiencyTest* ratioEffTauTight20ElecTauRunA  = new ratioEfficiencyTest(19.49,1.16,1.82,1.10,1.04);
  //ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunA = new ratioEfficiencyTest(19.37 ,2.34,  3.06, 1.37,1.0);
  //ratioEfficiencyTest* ratioEffTauMedium20ElecTauRunB = new ratioEfficiencyTest(19.47 ,2.4184,5.161,1.11,0.875);
  
  //float weightLoose15RunA  =  168.6 ;
  //float weightLoose20RunA  =  934.0 ;
  //float weightTight20RunA  =  785.7 ;
  //float weightMedium20RunA =  240.0 ;
  //float weightMedium20RunB = 2522.0 ;
  //float totalRun2011 = weightLoose15RunA+weightLoose20RunA+weightTight20RunA+weightMedium20RunA+weightMedium20RunB;

  //weightLoose15RunA  /= totalRun2011;
  //weightLoose20RunA  /= totalRun2011;
  //weightTight20RunA  /= totalRun2011;
  //weightMedium20RunA /= totalRun2011;
  //weightMedium20RunB /= totalRun2011;
 
  //Float_t xx = x[0];

  //return (ratioEffTauLoose15ElecTauRunA->turnOn(xx) *weightLoose15RunA  +
  //  ratioEffTauLoose20ElecTauRunA->turnOn(xx) *weightLoose20RunA  +
  //  ratioEffTauTight20ElecTauRunA->turnOn(xx) *weightTight20RunA  +
  //  ratioEffTauMedium20ElecTauRunA->turnOn(xx)*weightMedium20RunA +
  //  ratioEffTauMedium20ElecTauRunB->turnOn(xx)*weightMedium20RunB
  //  );

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

  TF1 *ratioMu15L110BL      = new TF1("ratioMu15L110BL",         myFuncRatioMu15L110BL,     14,400,0);
  TF1 *turnOnMu15L110BL     = new TF1("turnOnMu15L110BL",        myFuncTurnOnMu15L110BL,    14,400,0);

  TF1 *ratioMu15L110EC      = new TF1("ratioMu15L110EC",         myFuncRatioMu15L110EC,     14,400,0);
  TF1 *turnOnMu15L110EC     = new TF1("turnOnMu15L110EC",        myFuncTurnOnMu15L110EC,    14,400,0);

  TF1 *ratioMu15L114BL      = new TF1("ratioMu15L114BL",         myFuncRatioMu15L114BL,     14,400,0);
  TF1 *turnOnMu15L114BL     = new TF1("turnOnMu15L114BL",        myFuncTurnOnMu15L114BL,    14,400,0);

  TF1 *ratioMu15L114EC      = new TF1("ratioMu15L114EC",         myFuncRatioMu15L114EC,     14,400,0);
  TF1 *turnOnMu15L114EC     = new TF1("turnOnMu15L114EC",        myFuncTurnOnMu15L114EC,    14,400,0);
 
  TF1 *ratioMuAllBL         = new TF1("ratioMuAllBL",            myFuncRatioMuAllBL,        14,400,0);
  TF1 *turnOnMuAllBL        = new TF1("turnOnMuAllBL",           myFuncTurnOnMuAllBL,       14,400,0);

  TF1 *ratioMuAllEC         = new TF1("ratioMuAllEC",            myFuncRatioMuAllEC,        14,400,0);
  TF1 *turnOnMuAllEC        = new TF1("turnOnMuAllEC",           myFuncTurnOnMuAllEC,       14,400,0);
 
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
 
  TF1 *ratioTauMuTauAll      = new TF1("ratioTauMuTauAll", myFuncRatioTauMuTauAll ,18,400,0);
  TF1 *turnOnTauMuTauAll     = new TF1("turnOnTauMuTauAll",myFuncTurnOnTauMuTauAll,18,400,0);
 

  fout->cd();
  

  ratioElec15BL->SetNpx(1600);
  turnOnElec15BL->SetNpx(1600);
  ratioElec15EC->SetNpx(1600);
  turnOnElec15EC->SetNpx(1600);
  ratioElec18BL->SetNpx(1600);
  turnOnElec18BL->SetNpx(1600);
  ratioElec18EC->SetNpx(1600);
  turnOnElec18EC->SetNpx(1600);
  ratioElec20BL->SetNpx(1600);
  turnOnElec20BL->SetNpx(1600);
  ratioElec20EC->SetNpx(1600);
  turnOnElec20EC->SetNpx(1600);
  ratioElecAllBL->SetNpx(1600);
  turnOnElecAllBL->SetNpx(1600);
  ratioElecAllEC->SetNpx(1600);
  turnOnElecAllEC->SetNpx(1600);
  ratioMu15L114BL->SetNpx(1600);
  turnOnMu15L114BL->SetNpx(1600);
  ratioMu15L114EC->SetNpx(1600);
  turnOnMu15L114EC->SetNpx(1600);
  ratioMuAllBL->SetNpx(1600);
  turnOnMuAllBL->SetNpx(1600);
  ratioMuAllEC->SetNpx(1600);
  turnOnMuAllEC->SetNpx(1600);
  ratioTauElecTauAll->SetNpx(1600);
  turnOnTauElecTauAll->SetNpx(1600);
  ratioTauMuTauAll->SetNpx(1600);
  turnOnTauMuTauAll->SetNpx(1600);

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
  ratioMu15L110BL->Write();
  turnOnMu15L110BL->Write();
  ratioMu15L110EC->Write();
  turnOnMu15L110EC->Write();
  ratioMu15L114BL->Write();
  turnOnMu15L114BL->Write();
  ratioMu15L114EC->Write();
  turnOnMu15L114EC->Write();
  ratioMuAllBL->Write();
  turnOnMuAllBL->Write();
  ratioMuAllEC->Write();
  turnOnMuAllEC->Write();
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
  ratioTauMuTauAll->Write();
  turnOnTauMuTauAll->Write();

  fout->Write();
  fout->Close();


}
