#ifndef JET_HORNS
#define JET_HORNS

#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

using namespace std;

bool Jet_inHorns(RVecF CleanJet_pt,
		 RVecF CleanJet_eta,
		 bool minReq = false){

  bool lowpt = false;
  bool inhorns = false; 
  
  for (unsigned int i = 0; i<CleanJet_pt.size(); i++){
    if (minReq)
      lowpt = (CleanJet_pt[i] > 30.0 && CleanJet_pt[i] < 50.0);
    else
      lowpt = (CleanJet_pt[i] < 50.0);
    inhorns = (fabs(CleanJet_eta[i])>2.6 && fabs(CleanJet_eta[i])<3.1);
    if (lowpt && inhorns)
      return false;
  } 
  return true;
}

#endif
