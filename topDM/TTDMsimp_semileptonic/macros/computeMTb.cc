#include <vector>
#include <cmath>

// Compute MTb
double computeMTb(
    const ROOT::VecOps::RVec<float> &CleanJet_pt,
    const ROOT::VecOps::RVec<float> &CleanJet_eta,
    const ROOT::VecOps::RVec<float> &CleanJet_phi,
    const ROOT::VecOps::RVec<float> &CleanJet_btag,
    float btagWP,
    float MET_pt,
    float MET_phi
) {
    int bjet_idx = -1;
    double max_btag = -1.0;
    
    for (size_t i = 0; i < CleanJet_pt.size(); i++) {
        // Apply kinematic cuts first
        if (CleanJet_pt[i] < 30.0) continue;
        if (std::abs(CleanJet_eta[i]) > 2.5) continue;
    
        // Now consider b-tag score
        if (CleanJet_btag[i] > btagWP && CleanJet_btag[i] > max_btag) {
            max_btag = CleanJet_btag[i];
            bjet_idx = i;
        }
    }

    if(bjet_idx < 0) return 0.0;

    double pt = CleanJet_pt[bjet_idx];
    double phi = CleanJet_phi[bjet_idx];

    double mtb = sqrt(2.0 * pt * MET_pt * (1.0 - cos(DeltaPhi(phi, MET_phi))));
    return mtb;
}
