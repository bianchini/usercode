#include "PFAnalyses/PFCandidate/interface/TreePFCandidate.h"

using namespace PF; 


void Particle::setIsoDeposit( IsolationKeys key, const IsoDeposit& isoDep ) {
  isoDeposits_[key] = isoDep;
}

void Particle::setIsoDeposits(const IsoDeposits& isoDeps){
  isoDeposits_ = isoDeps; 
}





