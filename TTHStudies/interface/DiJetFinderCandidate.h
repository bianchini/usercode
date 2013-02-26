#ifndef DIJETCANDIDATEFINDER_H
#define DIJETCANDIDATEFINDER_H






class DiJetCandidateFinder{

 public:
  
  DiJetCandidateFinder(bool verbose, float ptCutCount, float etaCutCount, float ptCutRecover, float etaCutRecover, std::map<string, TF1*> func1D, std::map<string, TF2*> func2D):  
    verbose_(verbose), ptCutCount_(ptCutCount), etaCutCount_(etaCutCount), ptCutRecover_(ptCutRecover), etaCutRecover_(etaCutRecover),
    eventFlag_(-99),numJetCount_(-99), numJetRecov_(-99), error_(0) {
    
    for(std::map<string, TF1*>::iterator it = func1D.begin(); it!=func1D.end(); it++){
      func1D_[it->first] = it->second;
    }
    for(std::map<string, TF2*>::iterator it = func2D.begin(); it!=func2D.end(); it++){
      func2D_[it->first] = it->second;
    }
    
  }

  void run(int, std::vector<JetByPt>, std::vector<LV>);
  float jetCountMult();
  float jetRecoverMult();
  int eventFlag();
  int errorFlag();
    
 private:

  typedef struct
  {
    int index;
    float pt;
    float eta;
    float phi;
    float mass;
    float csv;
    float topB;
    float topW;
    float atopB;
    float atopW;
    float higgsB;
    float flavor;
    float unc;
    void reset(){
      index = -99; pt = -99; eta = -99; phi = -99; mass = -99; csv = -99; topB = -99; topW = -99;  atopW = -99;  atopB = -99;
      higgsB = -99; flavor = -99; unc = -99;
    }
    void set(int index_,  float pt_, float eta_, float phi_, float mass_,float csv_, float topB_, float topW_, float atopB_,
	     float atopW_, float higgsB_, float flavor_, float unc_){
      index = index_; 
      pt = pt_; 
      eta = eta_; 
      phi = phi_; 
      mass = mass_; 
      csv = csv_; 
      topB = topB_; 
      topW = topW_;  
      atopW = atopW_;  
      atopB = atopB_;
      higgsB = higgsB_; 
      flavor = flavor_; 
      unc = unc_;
    }
  } JetByPt;
  

  struct sorterByPt {
    bool operator() (float i,float j) const { return (i>j);}
  };

  struct sorterByJetPt {
    bool operator() (JetByPt i, JetByPt j) const { return (i.pt>j.pt);}
  };
    
  void runLJ(std::vector<JetByPt>, std::vector<LV>);
  void runLL(std::vector<JetByPt>, std::vector<LV>);

  bool findMatchLJ(int, std::vector<JetByPt>&, std::vector<LV>, std::map<unsigned int, int>&);

  LV buildNeutrino(LV, LV);

  bool verbose_;
  float ptCutCount_;
  float etaCutCount_;
  float ptCutRecover_;
  float etaCutRecover_;

  int eventFlag_;
  int numJetCount_;
  int numJetRecov_;

  int error_;

  std::map<float, unsigned int , sorterByPt> jMapPt_;
  std::map<float, unsigned int , sorterByPt> lMapPt_;

  std::map<string, TF1*> func1D_;
  std::map<string, TF2*> func2D_;

};

#endif
