#ifndef BTAGSF
#define BTAGSF

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

class btagSFbc {
public:
    btagSFbc(TString eff_map, const std::string year, TString algo_extension = "");
    ~btagSFbc();

    RVecF operator()(const RVecF& CleanJet_pt, const RVecF& CleanJet_eta, const RVecI& CleanJet_jetIdx, unsigned int nCleanJet,
                     const RVecI& Jet_hadronFlavour, const RVecF& Jet_btag, const std::string WP, const RVec<std::string>& systematic) {

        RVecF results(systematic.size(), 1.0f);

        const bool isUParTAK4 = (year_ == "2024_Summer24" || year_ == "2025_Summer24");
        const std::string correctionKey = isUParTAK4 ? "UParTAK4_comb" : "particleNet_comb";
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

                const int flavour = Jet_hadronFlavour[jetIdx];
                if (flavour != 4 && flavour != 5) continue;

                const double pt = CleanJet_pt[iJ];
                const double absEta = std::abs(static_cast<double>(CleanJet_eta[iJ]));
                const double sf = cset_btag->evaluate({systematic[i], WP, flavour, absEta, pt});

                if (Jet_btag[jetIdx] > wpCut) {
                    btag_sf *= sf;
                } else {
                    const float eff = getEff(CleanJet_pt[iJ], CleanJet_eta[iJ], flavour);
                    if (eff >= 1.0f) continue;
                    btag_sf *= (1.0 - eff * sf) / (1.0 - eff);
                }
            }

            results[i] = btag_sf;
        }

        return results;
    }

private:
    TH2* h_bjet_eff = nullptr;
    TH2* h_cjet_eff = nullptr;

    std::string year_;
    std::unique_ptr<CorrectionSet> cset;

    float getEff(float pt, float eta, int flavour) const;
    static int findClampedBin(const TAxis* axis, double value);
};

btagSFbc::btagSFbc(TString eff_map, const std::string year, TString algo_extension) : year_(year) {
    const std::string home = "/afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/processor/data/jsonpog-integration/POG/BTV/" + year;
    cset = CorrectionSet::from_file(home + "/btagging.json.gz");

    std::unique_ptr<TFile> reff(TFile::Open(eff_map, "READ"));
    if (!reff || reff->IsZombie()) throw std::runtime_error("btagSFbc: cannot open efficiency map " + std::string(eff_map.Data()));

    const TString bHistName = "bjet" + algo_extension + "_eff";
    const TString cHistName = "cjet" + algo_extension + "_eff";

    auto* hb = dynamic_cast<TH2*>(reff->Get(bHistName));
    auto* hc = dynamic_cast<TH2*>(reff->Get(cHistName));

    if (!hb || !hc) throw std::runtime_error("btagSFbc: missing efficiency histogram(s) " + std::string(bHistName.Data()) + " / " + std::string(cHistName.Data()));

    h_bjet_eff = static_cast<TH2*>(hb->Clone());
    h_cjet_eff = static_cast<TH2*>(hc->Clone());

    h_bjet_eff->SetDirectory(nullptr);
    h_cjet_eff->SetDirectory(nullptr);
}

int btagSFbc::findClampedBin(const TAxis* axis, double value) {
    return std::max(1, std::min(axis->FindBin(value), axis->GetNbins()));
}

float btagSFbc::getEff(float pt, float eta, int flavour) const {
    const TH2* hist = nullptr;

    if (flavour == 5) hist = h_bjet_eff;
    else if (flavour == 4) hist = h_cjet_eff;
    else return 1.0f;

    const int xbin = findClampedBin(hist->GetXaxis(), pt);
    const int ybin = findClampedBin(hist->GetYaxis(), eta);

    return hist->GetBinContent(xbin, ybin);
}

btagSFbc::~btagSFbc() {
    delete h_bjet_eff;
    delete h_cjet_eff;
}

#endif
