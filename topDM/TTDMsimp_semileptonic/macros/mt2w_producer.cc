#ifndef MT2W_PRODUCER_CC
#define MT2W_PRODUCER_CC

#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/macros/mt2w_bisect.h"
#include <vector>
#include <cmath>
#include "ROOT/RVec.hxx"

/***********************************************************************/
/*                                                                     */
/*              Finding MT2W                                           */
/*              Reference:  arXiv:1203.4813 [hep-ph]                   */
/*              Authors: Yang Bai, Hsin-Chia Cheng,                    */
/*                       Jason Gallicchio, Jiayin Gu                   */
/*              Based on MT2 by: Hsin-Chia Cheng, Zhenyu Han           */
/*              May 8, 2012, v1.00a                                    */
/*                                                                     */
/***********************************************************************/

using namespace ROOT::VecOps;

double mt2w_producer(float lep_pt,
                  float lep_eta,
                  float lep_phi,
                  const RVec<float>& jet_pt,
                  const RVec<float>& jet_eta,
                  const RVec<float>& jet_phi,
                  const RVec<float>& jet_btag,
                  float met_pt,
                  float met_phi,
                  float btagWP)
    {
        double mt2w = -1.0;
	    
	double lep_px = lep_pt * cos(lep_phi);
        double lep_py = lep_pt * sin(lep_phi);
        double lep_pz = lep_pt * sinh(lep_eta);
        double lep_E  = sqrt(lep_pt*lep_pt + lep_pz*lep_pz);

        std::vector<int> b_indices;
        std::vector<int> nonb_indices;
        for (size_t i = 0; i < jet_pt.size(); i++) {
            // Apply basic pT cut, and eta for btagged jets
            if (jet_pt[i] < 30.0) continue;
        
            if (jet_btag[i] > btagWP) {
                if (std::abs(jet_eta[i]) < 2.5)
                    b_indices.push_back(i);
            } else {
                nonb_indices.push_back(i);
            }
        }

        std::vector<std::pair<int,int>> jet_pairs;
        if (b_indices.size() >= 2) {
            for (size_t i = 0; i < b_indices.size(); i++)
                for (size_t j = i+1; j < b_indices.size(); j++)
                    jet_pairs.emplace_back(b_indices[i], b_indices[j]);
        } else if (b_indices.size() == 1) {
            int b = b_indices[0];
            size_t max_nonb = std::min((size_t)3, nonb_indices.size());
            for (size_t i = 0; i < max_nonb; i++)
                jet_pairs.emplace_back(b, nonb_indices[i]);
        } else {
            mt2w = -1.0;
            return mt2w;
        }

        double mt2w_min = 1e9;

        for (auto &pair : jet_pairs) {
            int j1 = pair.first;
            int j2 = pair.second;

            double b1_px = jet_pt[j1] * cos(jet_phi[j1]);
            double b1_py = jet_pt[j1] * sin(jet_phi[j1]);
            double b1_pz = jet_pt[j1] * sinh(jet_eta[j1]);
            double b1_E  = jet_pt[j1] * cosh(jet_eta[j1]);

            double b2_px = jet_pt[j2] * cos(jet_phi[j2]);
            double b2_py = jet_pt[j2] * sin(jet_phi[j2]);
            double b2_pz = jet_pt[j2] * sinh(jet_eta[j2]);
            double b2_E  = jet_pt[j2] * cosh(jet_eta[j2]);

            double pmiss[3] = {0., met_pt*cos(met_phi), met_pt*sin(met_phi)};
            double pl[4]  = {lep_E, lep_px, lep_py, lep_pz};
            double pb1[4] = {b1_E, b1_px, b1_py, b1_pz};
            double pb2[4] = {b2_E, b2_px, b2_py, b2_pz};

            mt2w_bisect::mt2w mt2w_event;
            mt2w_event.set_momenta(pl, pb1, pb2, pmiss);

            double mt2w_val = mt2w_event.get_mt2w();
            if (mt2w_val < mt2w_min) mt2w_min = mt2w_val;
        }

        mt2w = (mt2w_min > 1e8) ? -1.0 : mt2w_min;

	return mt2w;
}

#endif
