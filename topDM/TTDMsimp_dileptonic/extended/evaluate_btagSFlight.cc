#ifndef BTAGSFLIGHT
#define BTAGSFLIGHT

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "correction.h"

using namespace ROOT;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

class btagSFlight {
public:
    btagSFlight(TString eff_map, const std::string year, TString algo_extension = "");
    ~btagSFlight();

    RVecF operator()(const RVecF& CleanJet_pt, const RVecF& CleanJet_eta, const RVecI& CleanJet_jetIdx, unsigned int nCleanJet,
                     const RVecI& Jet_hadronFlavour, const RVecF& Jet_btag, const std::string WP, const RVec<std::string>& systematic) {

        RVecF results(systematic.size(), 1.0f);

        const bool isUParTAK4 = (year_ == "2024_Summer24" || year_ == "2025_Summer24");
        const std::string correctionKey = isUParTAK4 ? "UParTAK4_light" : "particleNet_light";
        const std::string wpKey = isUParTAK4 ? "UParTAK4_wp_values" : "particleNet_wp_values";

        auto cset_btag = cset->at(correctionKey);
        auto cset_wps = cset->at(wpKey);
        const double wpCut = cset_wps->evaluate({WP});

        const std::size_t nJets = std::min({static_cast<std::size_t>(nCleanJet), CleanJet_pt.size(), CleanJet_eta.size(), CleanJet_jetIdx.size()});

        for (std::size_t i = 0; i < systematic.size(); ++i) {
            float btag_sf = 1.0f;

            for (std::size_t iJ = 0; iJ < nJets; ++iJ) {
                if (CleanJet_pt[iJ] <= 30.0f || std::abs(CleanJet_eta[iJ]) >= 2.5f) continue;

                const int jetIdx = CleanJet_jetIdx[iJ];
                if (jetIdx < 0 || jetIdx >= static_cast<int>(Jet_btag.size()) || jetIdx >= static_cast<int>(Jet_hadronFlavour.size())) continue;

                if (Jet_hadronFlavour[jetIdx] != 0) continue;

                const double pt = CleanJet_pt[iJ];
                const double absEta = std::abs(static_cast<double>(CleanJet_eta[iJ]));
                const double sf = cset_btag->evaluate({systematic[i], WP, 0, absEta, pt});

                if (Jet_btag[jetIdx] > wpCut) {
                    btag_sf *= sf;
                } else {
                    const float eff = getEff(CleanJet_pt[iJ], CleanJet_eta[iJ], 0);
                    if (eff >= 1.0f) continue;
                    btag_sf *= (1.0 - eff * sf) / (1.0 - eff);
                }
            }

            results[i] = btag_sf;
        }

        return results;
    }

private:
    TH2* h_ljet_eff = nullptr;

    std::string year_;
    std::unique_ptr<CorrectionSet> cset;

    float getEff(float pt, float eta, int flavour) const;
    static int findClampedBin(const TAxis* axis, double value);
};

btagSFlight::btagSFlight(TString eff_map, const std::string year, TString algo_extension) : year_(year) {
    const std::string home = "/afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/processor/data/jsonpog-integration/POG/BTV/" + year;
    cset = CorrectionSet::from_file(home + "/btagging.json.gz");

    std::unique_ptr<TFile> reff(TFile::Open(eff_map, "READ"));
    if (!reff || reff->IsZombie()) throw std::runtime_error("btagSFlight: cannot open efficiency map " + std::string(eff_map.Data()));

    const TString lHistName = "ljet" + algo_extension + "_eff";
    auto* hl = dynamic_cast<TH2*>(reff->Get(lHistName));

    if (!hl) throw std::runtime_error("btagSFlight: missing efficiency histogram " + std::string(lHistName.Data()));

    h_ljet_eff = static_cast<TH2*>(hl->Clone());
    h_ljet_eff->SetDirectory(nullptr);
}

int btagSFlight::findClampedBin(const TAxis* axis, double value) {
    return std::max(1, std::min(axis->FindBin(value), axis->GetNbins()));
}

float btagSFlight::getEff(float pt, float eta, int flavour) const {
    if (flavour != 0) return 1.0f;

    const int xbin = findClampedBin(h_ljet_eff->GetXaxis(), pt);
    const int ybin = findClampedBin(h_ljet_eff->GetYaxis(), eta);

    return h_ljet_eff->GetBinContent(xbin, ybin);
}

btagSFlight::~btagSFlight() {
    delete h_ljet_eff;
}

#endif
