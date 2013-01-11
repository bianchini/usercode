#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

//#include "ratioEfficiencyElec.C"
#include "ratioEfficiencyTest.C"

// Weights of differents periods
#define weightRunA  0.04185 //0.067   L=810.99
#define weightRunB  0.22924 //0.367   L=4442.288743
#define weightRunC  0.35354 //0.566   L=6851.050214
#define weightRunD  0.37537 //0       L=7274
// total L(ABC)= 12104.329
// total L(D)  = 7274
// total ABCD  = 19378.329

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//ELEC ID
Double_t myFuncRatioElecIDBL(Double_t* x, Double_t *par) {  
  double ratio = 1.; 
  Float_t xx = x[0]; 
  if( xx < 24 )  
    ratio = 0.9130; 
  else if( xx >= 24 && xx< 30)  
    ratio = 0.9130; 
  else if( xx >= 30)
    ratio = 0.9567; 
  return ratio; 
}
Double_t myFuncTurnOnElecIDBL(Double_t* x, Double_t *par) {  
  double ratio = 1.;  
  Float_t xx = x[0];  
  if( xx < 24 )   
    ratio = 0.7532;  
  else if( xx >= 24 && xx< 30)   
    ratio = 0.7532;  
  else if( xx >= 30) 
    ratio = 0.8587;  
  return ratio;  
}
Double_t myFuncRatioElecIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;   
  Float_t xx = x[0];   
  if( xx < 24 )    
    ratio = 0.8509;   
  else if( xx >= 24 && xx< 30)    
    ratio = 0.8509 ;   
  else if( xx >= 30)  
    ratio = 0.9239 ;   
  return ratio;   
}
Double_t myFuncTurnOnElecIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;    
  Float_t xx = x[0];    
  if( xx < 24 )     
    ratio = 0.3631;    
  else if( xx >= 24 && xx< 30)     
    ratio = 0.3631;    
  else if( xx >= 30)   
    ratio = 0.5611 ;    
  return ratio;    
}

// Run D
Double_t myFuncRatioElecRunDIDBL(Double_t* x, Double_t *par) {  
  double ratio = 1.; 
  Float_t xx = x[0]; 
  if( xx < 24 )  
    ratio = 0.9236; 
  else if( xx >= 24 && xx< 30)  
    ratio = 0.9236; 
  else if( xx >= 30)
    ratio = 0.9497; 
  return ratio; 
}
Double_t myFuncTurnOnElecRunDIDBL(Double_t* x, Double_t *par) {  
  double ratio = 1.;  
  Float_t xx = x[0];  
  if( xx < 24 )   
    ratio = 0.7299;  
  else if( xx >= 24 && xx< 30)   
    ratio = 0.7299;  
  else if( xx >= 30) 
    ratio = 0.8510;  
  return ratio;  
}
Double_t myFuncRatioElecRunDIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;   
  Float_t xx = x[0];   
  if( xx < 24 )    
    ratio = 0.8082;   
  else if( xx >= 24 && xx< 30)    
    ratio = 0.8082;   
  else if( xx >= 30)  
    ratio = 0.9175;   
  return ratio;   
}
Double_t myFuncTurnOnElecRunDIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;    
  Float_t xx = x[0];    
  if( xx < 24 )     
    ratio = 0.3294;    
  else if( xx >= 24 && xx< 30)     
    ratio = 0.3294;    
  else if( xx >= 30)   
    ratio = 0.5491;    
  return ratio;    
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//ELECISO
Double_t myFuncRatioElecIsoBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 ) 
    ratio = 0.9602;
  else if( xx >= 24 && xx< 30) 
    ratio = 0.9602;
  else if( xx >= 30)
    ratio = 0.9858 ;
  return ratio;
}
Double_t myFuncTurnOnElecIsoBL(Double_t* x, Double_t *par) { 
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 ) 
    ratio = 0.7496;
  else if( xx >= 24 && xx< 30) 
    ratio = 0.7496;
  else if( xx >= 30)
    ratio = 0.8956;
  return ratio; 
}
Double_t myFuncRatioElecIsoEC(Double_t* x, Double_t *par) { 
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 )  
    ratio = 0.9661; 
  else if( xx >= 24 && xx< 30)  
    ratio = 0.9661; 
  else if( xx >= 30) 
    ratio = 0.9942; 
  return ratio;
}
Double_t myFuncTurnOnElecIsoEC(Double_t* x, Double_t *par) {  
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 )   
    ratio = 0.8244;  
  else if( xx >= 24 && xx< 30)   
    ratio = 0.8244;  
  else if( xx >= 30)  
    ratio = 0.9187 ;  
  return ratio;
}

// Run D
Double_t myFuncRatioElecRunDIsoBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 ) 
    ratio = 0.9440;
  else if( xx >= 24 && xx< 30) 
    ratio = 0.9440;
  else if( xx >= 30)
    ratio =  0.9764;
  return ratio;
}
Double_t myFuncTurnOnElecRunDIsoBL(Double_t* x, Double_t *par) { 
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 ) 
    ratio = 0.7302;
  else if( xx >= 24 && xx< 30) 
    ratio = 0.7302;
  else if( xx >= 30)
    ratio = 0.8823;
  return ratio; 
}
Double_t myFuncRatioElecRunDIsoEC(Double_t* x, Double_t *par) { 
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 )  
    ratio = 0.9757; 
  else if( xx >= 24 && xx< 30)  
    ratio = 0.9757; 
  else if( xx >= 30) 
    ratio = 0.9764; 
  return ratio;
}
Double_t myFuncTurnOnElecRunDIsoEC(Double_t* x, Double_t *par) {  
  double ratio = 1.;
  Float_t xx = x[0];
  if( xx < 24 )   
    ratio = 0.8304;  
  else if( xx >= 24 && xx< 30)   
    ratio = 0.8304;  
  else if( xx >= 30)  
    ratio = 0.9174;  
  return ratio;
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//ELECTRIGGER
//Run 2012A
Double_t myFuncRatioEleRunABL(Double_t* x, Double_t *par) {  
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffElecRunABL = new ratioEfficiencyTest(20.4669,1.20429,1.84954,1.38645,0.891122);
  ratioEfficiencyTest* fitEffElecBLMC   = new ratioEfficiencyTest(21.4136,0.000422,2.47314e-06,1.42487,1.00104);  
  Float_t xx = x[0];  
  return fitEffElecRunABL->turnOn(xx)/fitEffElecBLMC->turnOn(xx);  
}
Double_t myFuncTurnOnEleRunABL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* fitEffElecRunABL = new ratioEfficiencyTest(20.4669,1.20429,1.84954,1.38645,0.891122);
  Float_t xx = x[0];   
  return fitEffElecRunABL->turnOn(xx);
}
Double_t myFuncRatioEleRunAEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElecRunAEC = new ratioEfficiencyTest(21.4136,1.93922,2.43562,1.00186,51.947);
  ratioEfficiencyTest* fitEffElecECMC = new ratioEfficiencyTest(20.9985,0.002918,3.43131e-05,1.41479,1.06506);
  Float_t xx = x[0];   
  return fitEffElecRunAEC->turnOn(xx)/fitEffElecECMC->turnOn(xx);   
}
Double_t myFuncTurnOnEleRunAEC(Double_t* x, Double_t *par) {
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElecRunAEC = new ratioEfficiencyTest(21.4136,1.93922,2.43562,1.00186,51.947);
  Float_t xx = x[0];    
  return fitEffElecRunAEC->turnOn(xx);
}
//Run2012B
Double_t myFuncRatioEleRunBBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElecRunBBL = new ratioEfficiencyTest(22.8618,0.844755,1.07941,1.27956,1.07722);
  ratioEfficiencyTest* fitEffElecBLMC = new ratioEfficiencyTest(21.4136,0.000422,2.47314e-06,1.42487,1.00104);  
  Float_t xx = x[0];   
  return fitEffElecRunBBL->turnOn(xx)/fitEffElecBLMC->turnOn(xx);   
}
Double_t myFuncTurnOnEleRunBBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElecRunBBL = new ratioEfficiencyTest(22.8618,0.844755,1.07941,1.27956,1.07722);
  Float_t xx = x[0];    
  return fitEffElecRunBBL->turnOn(xx);
}
/////////
Double_t myFuncRatioEleRunBEC(Double_t* x, Double_t *par) {
  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElecRunBEC = new ratioEfficiencyTest(22.1045,1.08481,0.780119,1.91846,0.962174);
  ratioEfficiencyTest* fitEffElecECMC = new ratioEfficiencyTest(20.9985,0.002918,3.43131e-05,1.41479,1.06506);
  Float_t xx = x[0];    
  return fitEffElecRunBEC->turnOn(xx)/fitEffElecECMC->turnOn(xx);     
}
Double_t myFuncTurnOnEleRunBEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");     
  ratioEfficiencyTest* fitEffElecRunBEC = new ratioEfficiencyTest(22.1045,1.08481,0.780119,1.91846,0.962174);
  Float_t xx = x[0];     
  return fitEffElecRunBEC->turnOn(xx);
}
//Run2012C
Double_t myFuncRatioEleRunCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElecRunCBL = new ratioEfficiencyTest(22.8598,0.855666,1.02951,1.32713,1.05486);
  ratioEfficiencyTest* fitEffElecBLMC = new ratioEfficiencyTest(21.4136,0.000422,2.47314e-06,1.42487,1.00104);  
  Float_t xx = x[0];   
  return fitEffElecRunCBL->turnOn(xx)/fitEffElecBLMC->turnOn(xx);   
}
Double_t myFuncTurnOnEleRunCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffElecRunCBL = new ratioEfficiencyTest(22.8598,0.855666,1.02951,1.32713,1.05486);
  Float_t xx = x[0];    
  return fitEffElecRunCBL->turnOn(xx);
}
Double_t myFuncRatioEleRunCEC(Double_t* x, Double_t *par) {
  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElecRunCEC = new ratioEfficiencyTest(21.7643,1.45024,0.785753,3.14722,0.926788);
  ratioEfficiencyTest* fitEffElecECMC = new ratioEfficiencyTest(20.9985,0.002918,3.43131e-05,1.41479,1.06506);
  Float_t xx = x[0];    
  return fitEffElecRunCEC->turnOn(xx)/fitEffElecECMC->turnOn(xx);     
}
Double_t myFuncTurnOnEleRunCEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");     
  ratioEfficiencyTest* fitEffElecRunCEC = new ratioEfficiencyTest(21.7643,1.45024,0.785753,3.14722,0.926788);
  Float_t xx = x[0];     
  return fitEffElecRunCEC->turnOn(xx);
}

//Run2012D
Double_t myFuncRatioEleRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElecRunDBL = new ratioEfficiencyTest(23.2037, 0.947222, 1.29024, 1.09804, 1.53015 );
  ratioEfficiencyTest* fitEffElecBLMC = new ratioEfficiencyTest(21.4136,0.000422,2.47314e-06,1.42487,1.00104);  
  Float_t xx = x[0];   
  return fitEffElecRunDBL->turnOn(xx)/fitEffElecBLMC->turnOn(xx);   
}
Double_t myFuncTurnOnEleRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffElecRunDBL = new ratioEfficiencyTest(23.2037, 0.947222, 1.29024, 1.09804, 1.53015 );
  Float_t xx = x[0];    
  return fitEffElecRunDBL->turnOn(xx);
}
Double_t myFuncRatioEleRunDEC(Double_t* x, Double_t *par) {
  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElecRunDEC = new ratioEfficiencyTest(21.86, 0.979008, 0.505753, 2.2701, 0.94213 );
  ratioEfficiencyTest* fitEffElecECMC = new ratioEfficiencyTest(20.9985,0.002918,3.43131e-05,1.41479,1.06506);
  Float_t xx = x[0];    
  return fitEffElecRunDEC->turnOn(xx)/fitEffElecECMC->turnOn(xx);     
}
Double_t myFuncTurnOnEleRunDEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");     
  ratioEfficiencyTest* fitEffElecRunDEC = new ratioEfficiencyTest(21.86, 0.979008, 0.505753, 2.2701, 0.94213 );
  Float_t xx = x[0];     
  return fitEffElecRunDEC->turnOn(xx);
}

//Combined Electron 
Double_t myFuncRatioEleAllBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffElecRunABL = new ratioEfficiencyTest(20.4669,1.20429,1.84954,1.38645,0.891122);
  ratioEfficiencyTest* fitEffElecRunBBL = new ratioEfficiencyTest(22.8618,0.844755,1.07941,1.27956,1.07722);
  ratioEfficiencyTest* fitEffElecRunCBL = new ratioEfficiencyTest(22.8598,0.855666,1.02951,1.32713,1.05486);
  ratioEfficiencyTest* fitEffElecRunDBL = new ratioEfficiencyTest(23.2037, 0.947222, 1.29024, 1.09804, 1.53015 );
  ratioEfficiencyTest* fitEffElecBLMC = new ratioEfficiencyTest(21.4136,0.000422,2.47314e-06,1.42487,1.00104);  

  Float_t xx = x[0];
  return (fitEffElecRunABL->turnOn(xx) * weightRunA +
	  fitEffElecRunBBL->turnOn(xx) * weightRunB +
	  fitEffElecRunCBL->turnOn(xx) * weightRunC +
	  fitEffElecRunDBL->turnOn(xx) * weightRunD
	  )/fitEffElecBLMC->turnOn(xx);
}
Double_t myFuncTurnOnEleAllBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* fitEffElecRunABL = new ratioEfficiencyTest(20.4669,1.20429,1.84954,1.38645,0.891122);
  ratioEfficiencyTest* fitEffElecRunBBL = new ratioEfficiencyTest(22.8618,0.844755,1.07941,1.27956,1.07722);
  ratioEfficiencyTest* fitEffElecRunCBL = new ratioEfficiencyTest(22.8598,0.855666,1.02951,1.32713,1.05486);
  ratioEfficiencyTest* fitEffElecRunDBL = new ratioEfficiencyTest(23.2037, 0.947222, 1.29024, 1.09804, 1.53015 );

  Float_t xx = x[0];
  return (fitEffElecRunABL->turnOn(xx) * weightRunA +
	  fitEffElecRunBBL->turnOn(xx) * weightRunB +
	  fitEffElecRunCBL->turnOn(xx) * weightRunC +
	  fitEffElecRunDBL->turnOn(xx) * weightRunD
	  );
}
Double_t myFuncRatioEleAllEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffElecRunAEC = new ratioEfficiencyTest(21.4136,1.93922,2.43562,1.00186,51.947);
  ratioEfficiencyTest* fitEffElecRunBEC = new ratioEfficiencyTest(22.1045,1.08481,0.780119,1.91846,0.962174);
  ratioEfficiencyTest* fitEffElecRunCEC = new ratioEfficiencyTest(21.7643,1.45024,0.785753,3.14722,0.926788);
  ratioEfficiencyTest* fitEffElecRunDEC = new ratioEfficiencyTest(21.86, 0.979008, 0.505753, 2.2701, 0.94213 );
  ratioEfficiencyTest* fitEffElecECMC = new ratioEfficiencyTest(20.9985,0.002918,3.43131e-05,1.41479,1.06506);

  Float_t xx = x[0];
  return (fitEffElecRunAEC->turnOn(xx) * weightRunA +
	  fitEffElecRunBEC->turnOn(xx) * weightRunB +
	  fitEffElecRunCEC->turnOn(xx) * weightRunC +
	  fitEffElecRunDEC->turnOn(xx) * weightRunD
	  )/fitEffElecECMC->turnOn(xx);
}

Double_t myFuncTurnOnEleAllEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffElecRunAEC = new ratioEfficiencyTest(21.4136,1.93922,2.43562,1.00186,51.947);
  ratioEfficiencyTest* fitEffElecRunBEC = new ratioEfficiencyTest(22.1045,1.08481,0.780119,1.91846,0.962174);
  ratioEfficiencyTest* fitEffElecRunCEC = new ratioEfficiencyTest(21.7643,1.45024,0.785753,3.14722,0.926788);
  ratioEfficiencyTest* fitEffElecRunDEC = new ratioEfficiencyTest(21.86, 0.979008, 0.505753, 2.2701, 0.94213 );

  Float_t xx = x[0];
  return (fitEffElecRunAEC->turnOn(xx) * weightRunA +
	  fitEffElecRunBEC->turnOn(xx) * weightRunB +
	  fitEffElecRunCEC->turnOn(xx) * weightRunC +
	  fitEffElecRunDEC->turnOn(xx) * weightRunD
	  );
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//MUID
///eta < 0.8

Double_t myFuncTurnOnMuIdIso(Double_t *x, Double_t *par) {

  Double_t pt  = x[0];
  Double_t eta = par[0];
  Int_t run    = (Int_t)par[1];
  Int_t ratio  = (Int_t)par[2];
  Int_t idiso  = (Int_t)par[3]; // 0 : id ; 1 : iso

  const int nEta=3; // [0;0.8[  [0.8;1.2[  [1.2;2.1[
  const int nRun=6; // ABCD, MC-ABCD, ABC, MC-ABC, D, MC-D
  const int nPt =2; // ]-inf;30[ [30;+inf[

  // aborting cases //
  if(run<0 || run>=nRun) { cout << "Choose run>=0 and run<" << nRun << endl; return 0; }
  if(ratio<0 || ratio>1) { cout << "Choose ratio=0 or 1"            << endl; return 0; }
  if(idiso<0 || idiso>1) { cout << "Choose idiso=0 or 1"            << endl; return 0; }
  
  // define results
  //
  // ID
  //
  Double_t resID[nRun][nEta][nPt]
    = { { { 0.9505 , 0.9493} , { 0.9412 , 0.9391} , { 0.9415 , 0.9396} } ,  // ABCD
	{ { 0.9647 , 0.9631} , { 0.9586 , 0.9578} , { 0.9511 , 0.9490} } ,  // MC ABCD
	{ { 0.9524 , 0.9508} , { 0.9430 , 0.9406} , { 0.9430 , 0.9508} } ,  // ABC
	{ { 0.9649 , 0.9632} , { 0.9586 , 0.9581} , { 0.9512 , 0.9491} } ,  // MC ABC
	{ { 0.9485 , 0.9470} , { 0.9382 , 0.9370} , { 0.9392 , 0.9379} } ,  // D
	{ { 0.9644 , 0.9629} , { 0.9585 , 0.9576} , { 0.9510 , 0.9489} } }; // MC D
	//                  //                   //                         //
	//      [0;0.8[             [0.8;1.2[            [1.2;2.1[
  // ISO
  //
  Double_t resISO[nRun][nEta][nPt]
    = { { { 0.7591 , 0.9016} , { 0.8044 , 0.9218} , { 0.8491 , 0.9393} } ,  // ABCD
	{ { 0.7838 , 0.9133} , { 0.8201 , 0.9289} , { 0.8515 , 0.9382} } ,  // MC ABCD
	{ { 0.7642 , 0.9041} , { 0.8080 , 0.9239} , { 0.8497 , 0.9397} } ,  // ABC
	{ { 0.7866 , 0.9151} , { 0.8222 , 0.9361} , { 0.8531 , 0.9389} } ,  // MC ABC
	{ { 0.7493 , 0.8973} , { 0.7986 , 0.9187} , { 0.8486 , 0.9387} } ,  // D
	{ { 0.7794 , 0.9337} , { 0.8168 , 0.9270} , { 0.8491 , 0.9370} } }; // MC D
	//                  //                   //                         //
	//      [0;0.8[             [0.8;1.2[            [1.2;2.1[
	

  // choose eta index
  int iEta;
  if(abs(eta) < 0.8)      iEta = 0; // [0;0.8[
  else if(abs(eta) < 1.2) iEta = 1; // [0.8;1.2[
  else                    iEta = 2; // [1.2;2.1[
  
  // choose pt index
  int iPt;
  if(pt<30.) iPt=0;
  else       iPt=1;

  // Return relevant result //
  if(ratio==0) {
    if(idiso==0) return resID[run][iEta][iPt];
    else         return resISO[run][iEta][iPt];
  }
  else if(run==0 || run==2 || run==4) { // numerator is always data
    if(idiso==0) return resID[run+1][iEta][iPt]!=0 ? resID[run][iEta][iPt]/resID[run+1][iEta][iPt] : 0 ;
    else         return resISO[run+1][iEta][iPt]!=0 ? resISO[run][iEta][iPt]/resISO[run+1][iEta][iPt] : 0 ;
  }
  else return 0;
}

//
// MUON TRIGGER
//
Double_t myFuncTurnOnMu(Double_t *x, Double_t *par) { 

  Double_t pt  = x[0];
  Double_t eta = par[0];
  Int_t run    = (Int_t)par[1];
  Int_t ratio  = (Int_t)par[2];

  const int nEta=6; // ]-inf,-1.2[ [-1.2,-0.8[ [-0.8,0[ [0,0.8[ [0.8,1.2[ [1.2,+inf[
  const int nRun=7; // A, B, C, D, MC-old, ABCD, MC-new

  // aborting cases //
  if(run<0   || run>=nRun)   { cout << "Choose run>=0 and run<" << nRun << endl; return 0; }
  if(ratio<0 || ratio>1) { cout << "Choose ratio=0 or 1"            << endl; return 0; }

  ratioEfficiencyTest* fitEffMu[nRun][nEta];

  // 2012A
  fitEffMu[0][0] = new ratioEfficiencyTest(16.9993, 8.82202e-05, 7.91529e-08, 1.40792, 0.928102); 
  fitEffMu[0][1] = new ratioEfficiencyTest(16.9824, 0.0694986,   0.0186614,   1.66577, 0.908218); 
  fitEffMu[0][2] = new ratioEfficiencyTest(17.2736, 0.13896,     0.198452,    1.13119, 1.21897);
  fitEffMu[0][3] = new ratioEfficiencyTest(17.9605, 0.500059,    0.865294,    1.04633, 1.69027);
  fitEffMu[0][4] = new ratioEfficiencyTest(18.094,  0.607997,    0.89385,     1.36337, 0.92399); 
  fitEffMu[0][5] = new ratioEfficiencyTest(16.9805, 9.18396e-05, 2.81836e-08, 1.83783, 0.858988); 

  // 2012B
  fitEffMu[1][0] = new ratioEfficiencyTest(16.0015, 5.59745e-07, 1.3395e-07,  1.37357, 0.891284);
  fitEffMu[1][1] = new ratioEfficiencyTest(18.015,  0.0512973,   0.0603545,   1.36001, 0.907481);
  fitEffMu[1][2] = new ratioEfficiencyTest(16.4569, 0.214484,    0.302707,    1.42363, 0.982643);
  fitEffMu[1][3] = new ratioEfficiencyTest(15.9829, 0.0435624,   0.0196399,   1.71605, 0.967839);
  fitEffMu[1][4] = new ratioEfficiencyTest(17.4688, 0.0494554,   0.0628053,   1.34067, 0.904989);
  fitEffMu[1][5] = new ratioEfficiencyTest(16.0029, 4.01862e-05, 6.62491e-08, 1.42189, 0.880251);

  // 2012C
  fitEffMu[2][0] = new ratioEfficiencyTest(15.9974, 7.20337e-05,  7.72238e-08, 1.5461,  0.87064);
  fitEffMu[2][1] = new ratioEfficiencyTest(17.446,  0.760355,     1.58032,     1.0623,  1.10472);
  fitEffMu[2][2] = new ratioEfficiencyTest(15.9788, 0.044455,     0.0215911,   1.71024, 0.965673);
  fitEffMu[2][3] = new ratioEfficiencyTest(15.9762, 0.0552286,    0.0231409,   1.78576, 0.96848);
  fitEffMu[2][4] = new ratioEfficiencyTest(17.462,  0.804351,     1.62323,     1.22776, 0.900085);
  fitEffMu[2][5] = new ratioEfficiencyTest(16.0051, -4.10239e-05, 1.15509e-08, 1.82463, 0.865417);

  // 2012D ( computed for abs(eta) )
  fitEffMu[3][0] = new ratioEfficiencyTest(15.9994, 7.37077e-05, 7.21076e-08, 1.58178, 0.861339);
  fitEffMu[3][1] = new ratioEfficiencyTest(16.7041, 0.383545,    0.467605,    1.59941, 0.882451);
  fitEffMu[3][2] = new ratioEfficiencyTest(15.9852, 0.0428581,   0.0160247,   1.69952, 0.971443);
  fitEffMu[3][3] = new ratioEfficiencyTest(15.9852, 0.0428581,   0.0160247,   1.69952, 0.971443);
  fitEffMu[3][4] = new ratioEfficiencyTest(16.7041, 0.383545,    0.467605,    1.59941, 0.882451);
  fitEffMu[3][5] = new ratioEfficiencyTest(15.9994, 7.37077e-05, 7.21076e-08, 1.58178, 0.861339);
					   
  // MC-old
  fitEffMu[4][0] = new ratioEfficiencyTest(15.997,  8.73042e-05, 5.36172e-08, 1.67934, 0.871415);
  fitEffMu[4][1] = new ratioEfficiencyTest(17.3339, 0.768105,    1.31172,     1.35161, 0.942887);
  fitEffMu[4][2] = new ratioEfficiencyTest(15.959,  0.0229759,   0.00597735,  1.76124, 0.980734);
  fitEffMu[4][3] = new ratioEfficiencyTest(15.9618, 0.0587497,   0.0189749,   1.94016, 0.978294);
  fitEffMu[4][4] = new ratioEfficiencyTest(16.7859, 0.443337,    0.571078,    1.62214, 0.919211);
  fitEffMu[4][5] = new ratioEfficiencyTest(15.9974, 8.50572e-05, 5.53033e-08, 1.64714, 0.888026);
					   
  // ABCD
  fitEffMu[5][0] = new ratioEfficiencyTest(15.9825, 	7.90724e-05, 	5.49275e-08, 	1.6403, 	0.858285);
  fitEffMu[5][1] = new ratioEfficiencyTest(17.3283, 	0.707103, 	1.2047, 	1.3732, 	0.900519);
  fitEffMu[5][2] = new ratioEfficiencyTest(15.9828, 	0.0412999, 	0.0177441, 	1.66934, 	0.970097);
  fitEffMu[5][3] = new ratioEfficiencyTest(15.9802, 	0.0548775, 	0.020313, 	1.79791, 	0.968398);
  fitEffMu[5][4] = new ratioEfficiencyTest(16.8396, 	0.458636, 	0.633185, 	1.5706, 	0.8848);
  fitEffMu[5][5] = new ratioEfficiencyTest(15.9987, 	8.94398e-05, 	5.18549e-08, 	1.8342, 	0.854625);
					   
  // MC-new
  fitEffMu[6][0] = new ratioEfficiencyTest(16.0051, 	2.45144e-05, 	4.3335e-09, 	1.66134, 	0.87045);
  fitEffMu[6][1] = new ratioEfficiencyTest(17.3135, 	0.747636, 	1.21803, 	1.40611, 	0.934983);
  fitEffMu[6][2] = new ratioEfficiencyTest(15.9556, 	0.0236127, 	0.00589832, 	1.75409, 	0.981338);
  fitEffMu[6][3] = new ratioEfficiencyTest(15.9289, 	0.0271317, 	0.00448573, 	1.92101, 	0.978625);
  fitEffMu[6][4] = new ratioEfficiencyTest(16.5678, 	0.328333, 	0.354533, 	1.67085, 	0.916992);
  fitEffMu[6][5] = new ratioEfficiencyTest(15.997,	7.90069e-05, 	4.40036e-08, 	1.66272, 	0.884502);
					   
  // choose iEta //
  int iEta;
  if(eta < -1.2)      iEta = 0;
  else if(eta < -0.8) iEta = 1;
  else if(eta < 0)    iEta = 2;
  else if(eta < 0.8)  iEta = 3;
  else if(eta < 1.2)  iEta = 4;
  else                iEta = 5;


  // return relevant result //
  if(ratio==0) return fitEffMu[run][iEta]->turnOn(pt) ;

  else if(run>=0 && run<=3) 
    return fitEffMu[4][iEta]->turnOn(pt)!=0 ? fitEffMu[run][iEta]->turnOn(pt) / fitEffMu[4][iEta]->turnOn(pt) : 0;

  else if(run==5)
    return fitEffMu[6][iEta]->turnOn(pt)!=0 ? fitEffMu[run][iEta]->turnOn(pt) / fitEffMu[6][iEta]->turnOn(pt) : 0;

  else return 0;

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
//TAU TRIGGER
//
/////////////////////////////////////////////////////////
// HLT PFTauLoose20 / MC Summer12 53X / elec+tau
//
Double_t myFuncTurnOnTauElecTauMCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMCBL = new ratioEfficiencyTest(18.40815138, 1.53235636, 3.55989632, 1.74542709, 0.90118450);
  Float_t xx = x[0];
  return ratioEffTauMCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauElecTauMCEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauMCEC = new ratioEfficiencyTest(18.29028052, 1.56239255, 11.03605631, 155.89290151, 0.85683995);
  Float_t xx = x[0]; 
  return ratioEffTauMCEC->turnOn(xx); 
} 
/////////////////////////////////////////////////
// HLT PFTau All / e+tau
Double_t myFuncTurnOnTauElecTauBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauBL = new ratioEfficiencyTest(18.43442868,2.08967536,3.27357845,6.96327309,0.85564484);
  Float_t xx = x[0];
  return ratioEffTauElecTauBL->turnOn(xx);
}

Double_t myFuncRatioTauElecTauBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauBL = new ratioEfficiencyTest(18.43442868,2.08967536,3.27357845,6.96327309,0.85564484);
  ratioEfficiencyTest* ratioEffTauMCBL = new ratioEfficiencyTest(18.40815138, 1.53235636, 3.55989632, 1.74542709, 0.90118450);
  Float_t xx = x[0];
  return ratioEffTauElecTauBL->turnOn(xx)/ratioEffTauMCBL->turnOn(xx);
}

Double_t myFuncTurnOnTauElecTauEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauEC = new ratioEfficiencyTest(18.16839440,1.86184564,4.39116712,1.01410741,1.39240481);
  Float_t xx = x[0];
  return ratioEffTauElecTauEC->turnOn(xx);
}

Double_t myFuncRatioTauElecTauEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauEC = new ratioEfficiencyTest(18.16839440,1.86184564,4.39116712,1.01410741,1.39240481);
  ratioEfficiencyTest* ratioEffTauMCEC = new ratioEfficiencyTest(18.29028052, 1.56239255, 11.03605631, 155.89290151, 0.85683995);
  Float_t xx = x[0];
  return ratioEffTauElecTauEC->turnOn(xx)/ratioEffTauMCEC->turnOn(xx);
}

// Run 2012D
Double_t myFuncTurnOnTauElecTauRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauRunDBL = new ratioEfficiencyTest(18.73, 0.374578, 0.136068, 5.789410, 0.8638);
  Float_t xx = x[0];
  return ratioEffTauElecTauRunDBL->turnOn(xx);
}

Double_t myFuncRatioTauElecTauRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauRunDBL = new ratioEfficiencyTest(18.73, 0.374578, 0.136068, 5.789410, 0.8638);
  ratioEfficiencyTest* ratioEffTauMCBL = new ratioEfficiencyTest(19.22, 0.204905, 0.175676, 2.644803, 0.8974);
  Float_t xx = x[0];
  return ratioEffTauElecTauRunDBL->turnOn(xx)/ratioEffTauMCBL->turnOn(xx);
}

Double_t myFuncTurnOnTauElecTauRunDEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauRunDEC = new ratioEfficiencyTest(19.32, 0.146243, 0.123579, 3.126114, 0.8313);
  Float_t xx = x[0];
  return ratioEffTauElecTauRunDEC->turnOn(xx);
}

Double_t myFuncRatioTauElecTauRunDEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauElecTauRunDEC = new ratioEfficiencyTest(19.32, 0.146243, 0.123579, 3.126114, 0.8313);
  ratioEfficiencyTest* ratioEffTauMCEC = new ratioEfficiencyTest(18.62, 0.037935, 0.002134, 95.090919, 0.8515);
  Float_t xx = x[0];
  return ratioEffTauElecTauRunDEC->turnOn(xx)/ratioEffTauMCEC->turnOn(xx);
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// HLT PFTauLoose20 / MC Summer12 53X / mu+tau
//
Double_t myFuncTurnOnTauMuTauMCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauMCBL = new ratioEfficiencyTest(18.80484409, 0.19082817, 0.19983010, 1.81979820, 0.93270649);
  Float_t xx = x[0];
  return ratioEffTauMuTauMCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauMuTauMCEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauMCEC = new ratioEfficiencyTest(18.25975478, 1.32745225, 1.70380810, 149.18410074, 0.87377770);
  Float_t xx = x[0]; 
  return ratioEffTauMuTauMCEC->turnOn(xx); 
} 
/////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012A+B+C     / mu+tau
Double_t myFuncRatioTauMuTauRunABCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauRunABCBL   = new ratioEfficiencyTest(18.50940288, 1.62285299, 2.73232995, 1.79135412, 0.91481432);
  ratioEfficiencyTest* ratioEffTauMuTauMCBL = new ratioEfficiencyTest(18.80484409, 0.19082817, 0.19983010, 1.81979820, 0.93270649);
  Float_t xx = x[0];
  return ratioEffTauMuTauRunABCBL->turnOn(xx)/ratioEffTauMuTauMCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauMuTauRunABCBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauRunABCBL = new ratioEfficiencyTest(18.50940288, 1.62285299, 2.73232995, 1.79135412, 0.91481432);
  Float_t xx = x[0];
  return ratioEffTauMuTauRunABCBL->turnOn(xx);
}
Double_t myFuncRatioTauMuTauRunABCEC(Double_t* x, Double_t *par) {    
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauMuTauRunABCEC = new ratioEfficiencyTest(18.45678784, 0.68697618, 0.57008697, 3.73470825, 0.84747211);
  ratioEfficiencyTest* ratioEffTauMuTauMCEC = new ratioEfficiencyTest(18.25975478, 1.32745225, 1.70380810, 149.18410074, 0.87377770);
  Float_t xx = x[0]; 
  return ratioEffTauMuTauRunABCEC->turnOn(xx)/ratioEffTauMuTauMCEC->turnOn(xx); 
} 
Double_t myFuncTurnOnTauMuTauRunABCEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauMuTauRunABCEC = new ratioEfficiencyTest(18.45678784, 0.68697618, 0.57008697, 3.73470825, 0.84747211);
  Float_t xx = x[0];  
  return ratioEffTauMuTauRunABCEC->turnOn(xx);
}

// HLT PFTauLoose20 / Run2012D     / mu+tau
Double_t myFuncRatioTauMuTauRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauBL   = new ratioEfficiencyTest(18.5793, 0.2706, 0.1356, 2.4432, 0.9190); // Run D only
  ratioEfficiencyTest* ratioEffTauMuTauMCBL = new ratioEfficiencyTest(18.80484409, 0.19082817, 0.19983010, 1.81979820, 0.93270649);
  Float_t xx = x[0];
  return ratioEffTauMuTauBL->turnOn(xx)/ratioEffTauMuTauMCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauMuTauRunDBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauBL = new ratioEfficiencyTest(18.5793, 0.2706, 0.1356, 2.4432, 0.9190); // Run D only
  Float_t xx = x[0];
  return ratioEffTauMuTauBL->turnOn(xx);
}
Double_t myFuncRatioTauMuTauRunDEC(Double_t* x, Double_t *par) {    
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauMuTauEC   = new ratioEfficiencyTest(18.8394, 0.1354, 0.1648, 1.3407, 0.9696); // Run D only
  ratioEfficiencyTest* ratioEffTauMuTauMCEC = new ratioEfficiencyTest(18.25975478, 1.32745225, 1.70380810, 149.18410074, 0.87377770);
  Float_t xx = x[0]; 
  return ratioEffTauMuTauEC->turnOn(xx)/ratioEffTauMuTauMCEC->turnOn(xx); 
} 
Double_t myFuncTurnOnTauMuTauRunDEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauMuTauEC = new ratioEfficiencyTest(18.8394, 0.1354, 0.1648, 1.3407, 0.9696); // Run D only
  Float_t xx = x[0];  
  return ratioEffTauMuTauEC->turnOn(xx);
}

// Tau trigger in mu+tau channel : 2012ABCD all together
Double_t myFuncRatioTauMuTauBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauBL   = new ratioEfficiencyTest(18.53484769, 1.55026068, 2.78382728, 1.40048225, 0.94356699);
  ratioEfficiencyTest* ratioEffTauMuTauMCBL = new ratioEfficiencyTest(18.90334548, 0.08929792, 0.10788389, 1.47675914, 0.95116761);
  Float_t xx = x[0];
  return ratioEffTauMuTauBL->turnOn(xx)/ratioEffTauMuTauMCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauMuTauBL(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauMuTauBL = new ratioEfficiencyTest(18.53484769, 1.55026068, 2.78382728, 1.40048225, 0.94356699);
  Float_t xx = x[0];
  return ratioEffTauMuTauBL->turnOn(xx);
}
Double_t myFuncRatioTauMuTauEC(Double_t* x, Double_t *par) {    
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauMuTauEC   = new ratioEfficiencyTest(18.28798198, 0.95234835, 0.75147523, 7.78047104, 0.84161836);
  ratioEfficiencyTest* ratioEffTauMuTauMCEC = new ratioEfficiencyTest(18.65016918, 0.40513575, 0.38648072, 3.63287413, 0.87121168);
  Float_t xx = x[0]; 
  return ratioEffTauMuTauEC->turnOn(xx)/ratioEffTauMuTauMCEC->turnOn(xx); 
} 
Double_t myFuncTurnOnTauMuTauEC(Double_t* x, Double_t *par) { 
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauMuTauEC = new ratioEfficiencyTest(18.28798198, 0.95234835, 0.75147523, 7.78047104, 0.84161836);
  Float_t xx = x[0];  
  return ratioEffTauMuTauEC->turnOn(xx);
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void makeFile(){

  TFile* fout = new TFile("llrCorrections_Moriond.root","RECREATE");

  TF1 *ratioElecIDBL        = new TF1("ratioElecIDBL",           myFuncRatioElecIDBL ,      14,800,0);
  TF1 *turnOnElecIDBL       = new TF1("turnOnElecIDBL",          myFuncTurnOnElecIDBL ,     14,800,0);
  TF1 *ratioElecIDEC        = new TF1("ratioElecIDEC",           myFuncRatioElecIDEC ,      14,800,0);
  TF1 *turnOnElecIDEC       = new TF1("turnOnElecIDEC",          myFuncTurnOnElecIDEC ,     14,800,0);

  TF1 *ratioElecIsoBL       = new TF1("ratioElecIsoBL",          myFuncRatioElecIsoBL ,     14,800,0);
  TF1 *turnOnElecIsoBL      = new TF1("turnOnElecIsoBL",         myFuncTurnOnElecIsoBL ,    14,800,0);
  TF1 *ratioElecIsoEC       = new TF1("ratioElecIsoEC",          myFuncRatioElecIsoEC ,     14,800,0);
  TF1 *turnOnElecIsoEC      = new TF1("turnOnElecIsoEC",         myFuncTurnOnElecIsoEC ,    14,800,0);

  TF1 *ratioElecRunDIDBL        = new TF1("ratioElecRunDIDBL",           myFuncRatioElecRunDIDBL ,      14,800,0);
  TF1 *turnOnElecRunDIDBL       = new TF1("turnOnElecRunDIDBL",          myFuncTurnOnElecRunDIDBL ,     14,800,0);
  TF1 *ratioElecRunDIDEC        = new TF1("ratioElecRunDIDEC",           myFuncRatioElecRunDIDEC ,      14,800,0);
  TF1 *turnOnElecRunDIDEC       = new TF1("turnOnElecRunDIDEC",          myFuncTurnOnElecRunDIDEC ,     14,800,0);

  TF1 *ratioElecRunDIsoBL       = new TF1("ratioElecRunDIsoBL",          myFuncRatioElecRunDIsoBL ,     14,800,0);
  TF1 *turnOnElecRunDIsoBL      = new TF1("turnOnElecRunDIsoBL",         myFuncTurnOnElecRunDIsoBL ,    14,800,0);
  TF1 *ratioElecRunDIsoEC       = new TF1("ratioElecRunDIsoEC",          myFuncRatioElecRunDIsoEC ,     14,800,0);
  TF1 *turnOnElecRunDIsoEC      = new TF1("turnOnElecRunDIsoEC",         myFuncTurnOnElecRunDIsoEC ,    14,800,0);

  TF1 *ratioElecRunABL      = new TF1("ratioElecRunABL",         myFuncRatioEleRunABL,      24,800,0);
  TF1 *turnOnElecRunABL     = new TF1("turnOnElecRunABL",        myFuncTurnOnEleRunABL,     24,800,0);
  TF1 *ratioElecRunAEC      = new TF1("ratioElecRunAEC",         myFuncRatioEleRunAEC,      24,800,0);
  TF1 *turnOnElecRunAEC     = new TF1("turnOnElecRunAEC",        myFuncTurnOnEleRunAEC,     24,800,0);

  TF1 *ratioElecRunBBL      = new TF1("ratioElecRunBBL",         myFuncRatioEleRunBBL,      24,800,0);
  TF1 *turnOnElecRunBBL     = new TF1("turnOnElecRunBBL",        myFuncTurnOnEleRunBBL,     24,800,0);
  TF1 *ratioElecRunBEC      = new TF1("ratioElecRunBEC",         myFuncRatioEleRunBEC,      24,800,0);
  TF1 *turnOnElecRunBEC     = new TF1("turnOnElecRunBEC",        myFuncTurnOnEleRunBEC,     24,800,0);

  TF1 *ratioElecRunCBL      = new TF1("ratioElecRunCBL",         myFuncRatioEleRunCBL,      24,800,0);
  TF1 *turnOnElecRunCBL     = new TF1("turnOnElecRunCBL",        myFuncTurnOnEleRunCBL,     24,800,0);
  TF1 *ratioElecRunCEC      = new TF1("ratioElecRunCEC",         myFuncRatioEleRunCEC,      24,800,0);
  TF1 *turnOnElecRunCEC     = new TF1("turnOnElecRunCEC",        myFuncTurnOnEleRunCEC,     24,800,0);

  TF1 *ratioElecRunDBL      = new TF1("ratioElecRunDBL",         myFuncRatioEleRunDBL,      24,800,0);
  TF1 *turnOnElecRunDBL     = new TF1("turnOnElecRunDBL",        myFuncTurnOnEleRunDBL,     24,800,0);
  TF1 *ratioElecRunDEC      = new TF1("ratioElecRunDEC",         myFuncRatioEleRunDEC,      24,800,0);
  TF1 *turnOnElecRunDEC     = new TF1("turnOnElecRunDEC",        myFuncTurnOnEleRunDEC,     24,800,0);

  TF1 *ratioElecAllBL       = new TF1("ratioElecAllBL",          myFuncRatioEleAllBL,       24,800,0);
  TF1 *turnOnElecAllBL      = new TF1("turnOnElecAllBL",         myFuncTurnOnEleAllBL,      24,800,0);
  TF1 *ratioElecAllEC       = new TF1("ratioElecAllEC",          myFuncRatioEleAllEC,       24,800,0);
  TF1 *turnOnElecAllEC      = new TF1("turnOnElecAllEC",         myFuncTurnOnEleAllEC,      24,800,0);

  // MUON //
  TF1 *turnOnMuIdIso = new TF1("turnOnMuIdIso", myFuncTurnOnMuIdIso, 0,800,4); // don't forget to SetParameters !
  TF1 *turnOnMu      = new TF1("turnOnMu",      myFuncTurnOnMu,      0,800,3); // don't forget to SetParameters !
  //////////

  TF1 *turnOnTauElecTauMCBL      = new TF1("turnOnTauElecTauMCBL",   myFuncTurnOnTauElecTauMCBL, 20,800,0);
  TF1 *turnOnTauElecTauMCEC      = new TF1("turnOnTauElecTauMCEC",   myFuncTurnOnTauElecTauMCEC, 20,800,0);

  TF1 *ratioTauElecTauBL    = new TF1("ratioTauElecTauBL", myFuncRatioTauElecTauBL, 20,800,0);
  TF1 *turnOnTauElecTauBL   = new TF1("turnOnTauElecTauBL", myFuncTurnOnTauElecTauBL, 20,800,0);
  TF1 *ratioTauElecTauEC    = new TF1("ratioTauElecTauEC", myFuncRatioTauElecTauEC, 20,800,0);
  TF1 *turnOnTauElecTauEC   = new TF1("turnOnTauElecTauEC", myFuncTurnOnTauElecTauEC, 20,800,0);

  TF1 *ratioTauElecTauRunDBL    = new TF1("ratioTauElecTauRunDBL", myFuncRatioTauElecTauRunDBL, 20,800,0);
  TF1 *turnOnTauElecTauRunDBL   = new TF1("turnOnTauElecTauRunDBL", myFuncTurnOnTauElecTauRunDBL, 20,800,0);
  TF1 *ratioTauElecTauRunDEC    = new TF1("ratioTauElecTauRunDEC", myFuncRatioTauElecTauRunDEC, 20,800,0);
  TF1 *turnOnTauElecTauRunDEC   = new TF1("turnOnTauElecTauRunDEC", myFuncTurnOnTauElecTauRunDEC, 20,800,0);

  TF1 *turnOnTauMuTauMCBL      = new TF1("turnOnTauMuTauMCBL",   myFuncTurnOnTauMuTauMCBL,      20,800,0);
  TF1 *turnOnTauMuTauMCEC      = new TF1("turnOnTauMuTauMCEC",   myFuncTurnOnTauMuTauMCEC,	20,800,0); 

  TF1 *ratioTauMuTauBL      = new TF1("ratioTauMuTauBL", myFuncRatioTauMuTauBL ,20,800,0);
  TF1 *turnOnTauMuTauBL     = new TF1("turnOnTauMuTauBL",myFuncTurnOnTauMuTauBL,20,800,0);
  TF1 *ratioTauMuTauEC      = new TF1("ratioTauMuTauEC", myFuncRatioTauMuTauEC ,20,800,0); 
  TF1 *turnOnTauMuTauEC     = new TF1("turnOnTauMuTauEC",myFuncTurnOnTauMuTauEC,20,800,0);
  
  TF1 *ratioTauMuTauRunABCBL  = new TF1("ratioTauMuTauRunABCBL", myFuncRatioTauMuTauRunABCBL ,20,800,0);
  TF1 *turnOnTauMuTauRunABCBL = new TF1("turnOnTauMuTauRunABCBL",myFuncTurnOnTauMuTauRunABCBL,20,800,0);
  TF1 *ratioTauMuTauRunABCEC  = new TF1("ratioTauMuTauRunABCEC", myFuncRatioTauMuTauRunABCEC ,20,800,0); 
  TF1 *turnOnTauMuTauRunABCEC = new TF1("turnOnTauMuTauRunABCEC",myFuncTurnOnTauMuTauRunABCEC,20,800,0);

  TF1 *ratioTauMuTauRunDBL  = new TF1("ratioTauMuTauRunDBL", myFuncRatioTauMuTauRunDBL ,20,800,0);
  TF1 *turnOnTauMuTauRunDBL = new TF1("turnOnTauMuTauRunDBL",myFuncTurnOnTauMuTauRunDBL,20,800,0);
  TF1 *ratioTauMuTauRunDEC  = new TF1("ratioTauMuTauRunDEC", myFuncRatioTauMuTauRunDEC ,20,800,0); 
  TF1 *turnOnTauMuTauRunDEC = new TF1("turnOnTauMuTauRunDEC",myFuncTurnOnTauMuTauRunDEC,20,800,0);

  fout->cd();

  ratioElecIDBL->SetNpx(25600);
  turnOnElecIDBL->SetNpx(25600);
  ratioElecIDEC->SetNpx(25600);
  turnOnElecIDEC->SetNpx(25600);

  ratioElecIsoBL->SetNpx(25600);
  turnOnElecIsoBL->SetNpx(25600);
  ratioElecIsoEC->SetNpx(25600);
  turnOnElecIsoEC->SetNpx(25600);

  ratioElecRunABL->SetNpx(25600);
  turnOnElecRunABL->SetNpx(25600);
  ratioElecRunAEC->SetNpx(25600);
  turnOnElecRunAEC->SetNpx(25600);
  
  ratioElecRunBBL->SetNpx(25600);
  turnOnElecRunBBL->SetNpx(25600);
  ratioElecRunBEC->SetNpx(25600);
  turnOnElecRunBEC->SetNpx(25600);
  
  ratioElecRunCBL->SetNpx(25600);
  turnOnElecRunCBL->SetNpx(25600);
  ratioElecRunCEC->SetNpx(25600);
  turnOnElecRunCEC->SetNpx(25600);
  
  ratioElecRunDBL->SetNpx(25600);
  turnOnElecRunDBL->SetNpx(25600);
  ratioElecRunDEC->SetNpx(25600);
  turnOnElecRunDEC->SetNpx(25600);
  
  ratioElecAllBL->SetNpx(25600);
  turnOnElecAllBL->SetNpx(25600);
  ratioElecAllEC->SetNpx(25600);
  turnOnElecAllEC->SetNpx(25600);

  //
  turnOnMu->SetNpx(25600);
  turnOnMuIdIso->SetNpx(25600);
  //

  turnOnTauElecTauMCBL->SetNpx(25600);
  turnOnTauElecTauMCEC->SetNpx(25600);

  ratioTauElecTauBL->SetNpx(25600);
  turnOnTauElecTauBL->SetNpx(25600);
  ratioTauElecTauEC->SetNpx(25600);
  turnOnTauElecTauEC->SetNpx(25600);

  ratioTauElecTauRunDBL->SetNpx(25600);
  turnOnTauElecTauRunDBL->SetNpx(25600);
  ratioTauElecTauRunDEC->SetNpx(25600);
  turnOnTauElecTauRunDEC->SetNpx(25600);

  turnOnTauMuTauMCBL->SetNpx(25600);
  turnOnTauMuTauMCEC->SetNpx(25600);

  ratioTauMuTauBL->SetNpx(25600);
  turnOnTauMuTauBL->SetNpx(25600);
  ratioTauMuTauEC->SetNpx(25600); 
  turnOnTauMuTauEC->SetNpx(25600);

  ratioTauMuTauRunABCBL->SetNpx(25600);
  turnOnTauMuTauRunABCBL->SetNpx(25600);
  ratioTauMuTauRunABCEC->SetNpx(25600); 
  turnOnTauMuTauRunABCEC->SetNpx(25600);

  ratioTauMuTauRunDBL->SetNpx(25600);
  turnOnTauMuTauRunDBL->SetNpx(25600);
  ratioTauMuTauRunDEC->SetNpx(25600); 
  turnOnTauMuTauRunDEC->SetNpx(25600);

  //==============================
  ratioElecIDBL->Write();
  turnOnElecIDBL->Write();
  ratioElecIDEC->Write();
  turnOnElecIDEC->Write();

  ratioElecIsoBL->Write();
  turnOnElecIsoBL->Write();
  ratioElecIsoEC->Write();
  turnOnElecIsoEC->Write();

  ratioElecRunABL->Write();
  turnOnElecRunABL->Write();
  ratioElecRunAEC->Write();
  turnOnElecRunAEC->Write();
  
  ratioElecRunBBL->Write();
  turnOnElecRunBBL->Write();
  ratioElecRunBEC->Write();
  turnOnElecRunBEC->Write();
  
  ratioElecRunCBL->Write();
  turnOnElecRunCBL->Write();
  ratioElecRunCEC->Write();
  turnOnElecRunCEC->Write();
  
  ratioElecRunDBL->Write();
  turnOnElecRunDBL->Write();
  ratioElecRunDEC->Write();
  turnOnElecRunDEC->Write();
  
  ratioElecAllBL->Write();
  turnOnElecAllBL->Write();
  ratioElecAllEC->Write();
  turnOnElecAllEC->Write();

  //
  turnOnMu->Write();
  turnOnMuIdIso->Write();
  //

  turnOnTauElecTauMCBL->Write();
  turnOnTauElecTauMCEC->Write();

  ratioTauElecTauBL->Write();
  turnOnTauElecTauBL->Write();
  ratioTauElecTauEC->Write();
  turnOnTauElecTauEC->Write();

  ratioTauElecTauRunDBL->Write();
  turnOnTauElecTauRunDBL->Write();
  ratioTauElecTauRunDEC->Write();
  turnOnTauElecTauRunDEC->Write();

  turnOnTauMuTauMCBL->Write();
  turnOnTauMuTauMCEC->Write();

  ratioTauMuTauBL->Write();
  turnOnTauMuTauBL->Write();
  ratioTauMuTauEC->Write(); 
  turnOnTauMuTauEC->Write();

  ratioTauMuTauRunABCBL->Write();
  turnOnTauMuTauRunABCBL->Write();
  ratioTauMuTauRunABCEC->Write(); 
  turnOnTauMuTauRunABCEC->Write();

  ratioTauMuTauRunDBL->Write();
  turnOnTauMuTauRunDBL->Write();
  ratioTauMuTauRunDEC->Write(); 
  turnOnTauMuTauRunDEC->Write();

  fout->Write();
  fout->Close();
}
