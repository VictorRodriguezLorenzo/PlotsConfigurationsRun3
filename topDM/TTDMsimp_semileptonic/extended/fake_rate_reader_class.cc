#ifndef FAKERATEREADER_1L
#define FAKERATEREADER_1L

#include "TSystem.h"
#include <iostream>
#include <vector>
#include <tuple>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH2F.h"
#include "TFile.h"
#include <map>
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace std;

typedef std::map<std::string,std::map<std::string,std::string>> map_dict;

class fake_rate_reader_1l {
public:
    fake_rate_reader_1l(TString year, TString ele_WP, TString muon_WP, TString kind, TString electron_tight_charge);

    std::tuple<double,double> GetRate(TH2F* fake_rate_histo, double pt, double eta, double lepton_pt_max);
    float GetFR_1l(double pt, double eta, double pdg, double isTight, TH2F* fake_rate_ele_, TH2F* fake_rate_muon_, TString stat);

    TString year_;
    TString electron_tight_charge_;
    TString ele_WP_;
    TString muon_WP_;
    TString kind_;

    RVecI Lepton_pdgId;
    RVecF Lepton_pt;
    RVecF Lepton_eta;
    RVecI Lepton_isTightMuon_cut_Tight_HWWW;
    RVecI Lepton_isTightElectron_mvaFall17V2Iso_WP90;
    RVecF CleanJet_pt;
    int nCleanJet;

    TH2F* fake_rate_muon_10_;
    TH2F* fake_rate_muon_15_;
    TH2F* fake_rate_muon_20_;
    TH2F* fake_rate_muon_25_;
    TH2F* fake_rate_muon_30_;
    TH2F* fake_rate_muon_35_;
    TH2F* fake_rate_muon_45_;
    TH2F* fake_rate_ele_25_;
    TH2F* fake_rate_ele_35_;
    TH2F* fake_rate_ele_45_;
    TH2F* prompt_rate_muon_;
    TH2F* prompt_rate_ele_;

    int isTight_;

    float operator()(RVecI Lepton_pdgId,
                     RVecF Lepton_pt,
                     RVecF Lepton_eta,
                     RVecI Lepton_isTightMuon_cut_Tight_HWWW,
                     RVecI Lepton_isTightElectron_mvaFall17V2Iso_WP90,
                     RVecF CleanJet_pt,
                     int nCleanJet) {

	double SF = 1.;
        isTight_ = 0;
        this->Lepton_pdgId = Lepton_pdgId;
        this->Lepton_pt = Lepton_pt;
        this->Lepton_eta = Lepton_eta;
        this->Lepton_isTightMuon_cut_Tight_HWWW = Lepton_isTightMuon_cut_Tight_HWWW;
        this->Lepton_isTightElectron_mvaFall17V2Iso_WP90 = Lepton_isTightElectron_mvaFall17V2Iso_WP90;
        this->CleanJet_pt = CleanJet_pt;
        this->nCleanJet = nCleanJet;

        // Determine tightness
        if (TMath::Abs(Lepton_pdgId[0]) == 11)
            if (Lepton_isTightElectron_mvaFall17V2Iso_WP90[0] == 1)
                isTight_ = 1;
        if (TMath::Abs(Lepton_pdgId[0]) == 13)
            if (Lepton_isTightMuon_cut_Tight_HWWW[0] == 1)
                isTight_ = 1;

        // Compute per-jet-bin weights
        float w0j = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_,
                             fake_rate_ele_35_, fake_rate_muon_20_, "Nominal");
        float w1j = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_,
                             fake_rate_ele_35_, fake_rate_muon_25_, "Nominal");
        float w2j = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_,
                             fake_rate_ele_35_, fake_rate_muon_35_, "Nominal");

        // Assign based on jet bins
        float fakeWeight = w0j*( nCleanJet == 0 || CleanJet_pt[0] <  30) +
                           w1j*((nCleanJet == 1 && CleanJet_pt[0] >= 30) ||
                                (nCleanJet >  1 && CleanJet_pt[0] >= 30 && CleanJet_pt[1] < 30)) +
                           w2j*( nCleanJet >  1 && CleanJet_pt[1] >= 30);

        if (kind_ == "nominal") return fakeWeight;

        // -------------------------------
        // Electron systematics
        // -------------------------------
        if (kind_ == "EleUp") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_45_, fake_rate_muon_20_, "Nominal");
            return num/fakeWeight;
        }
        if (kind_ == "EleDown") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_25_, fake_rate_muon_20_, "Nominal");
            return num/fakeWeight;
        }
        if (kind_ == "StatEleUp") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_20_, "ElUp");
            return num/fakeWeight;
        }
        if (kind_ == "StatEleDown") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_20_, "ElDown");
            return num/fakeWeight;
        }

        // -------------------------------
        // Muon systematics
        // -------------------------------
        if (kind_ == "MuUp") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_30_, "Nominal");
            return num/fakeWeight;
        }
        if (kind_ == "MuDown") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_15_, "Nominal");
            return num/fakeWeight;
        }
        if (kind_ == "StatMuUp") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_20_, "MuUp");
            return num/fakeWeight;
        }
        if (kind_ == "StatMuDown") {
            float num = GetFR_1l(Lepton_pt[0], Lepton_eta[0], Lepton_pdgId[0], isTight_, fake_rate_ele_35_, fake_rate_muon_20_, "MuDown");
            return num/fakeWeight;
        }

        return -1;
    }

private:
    void loadSF2D(std::string filename);
    std::tuple<double,double> GetSF(double pt_in, double eta_in);
};

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
fake_rate_reader_1l::fake_rate_reader_1l(TString year, TString ele_WP, TString muon_WP, TString kind, TString electron_tight_charge) {

    cout << "Year:             " << year                  << endl;
    cout << "Ele WP:           " << ele_WP                << endl;
    cout << "Muon WP:          " << muon_WP               << endl;
    cout << "Kind:             " << kind                  << endl;
    cout << "Ele tight charge: " << electron_tight_charge << endl;
    year_ = year;
    ele_WP_ = ele_WP;
    muon_WP_ = muon_WP;
    kind_ = kind;
    electron_tight_charge_ = electron_tight_charge;

    TString mkShapesRDF_base = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/";

    // Muon FR
    fake_rate_muon_10_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet10.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_15_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet15.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_20_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet20.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_25_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet25.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_30_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet30.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_35_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet35.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_muon_45_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonFR_jet45.root")->Get("FR_pT_eta_EWKcorr"));

    // Electron FR
    fake_rate_ele_25_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cutBased_LooseID_tthMVA_" + ele_WP + "/EleFR_jet25.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_ele_35_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cutBased_LooseID_tthMVA_" + ele_WP + "/EleFR_jet35.root")->Get("FR_pT_eta_EWKcorr"));
    fake_rate_ele_45_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cutBased_LooseID_tthMVA_" + ele_WP + "/EleFR_jet45.root")->Get("FR_pT_eta_EWKcorr"));

    // Prompt rates
    prompt_rate_muon_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cut_TightID_pfIsoLoose_HWW_tthmva_" + muon_WP + "/MuonPR.root")->Get("h_Muon_signal_pt_eta_bin"));
    prompt_rate_ele_ = (TH2F*) (TFile::Open(mkShapesRDF_base + "/data/fakerate/" + year + "/cutBased_LooseID_tthMVA_" + ele_WP + "/ElePR.root")->Get("h_Ele_signal_pt_eta_bin"));
}

// ------------------------------------------------------------
// Get FR from histogram
// ------------------------------------------------------------
std::tuple<double,double> fake_rate_reader_1l::GetRate(TH2F* fake_rate_histo, double pt, double eta, double lepton_pt_max){
    double aeta = abs(eta);
    int nbinsx = fake_rate_histo->GetNbinsX();
    if (lepton_pt_max <= 0.) lepton_pt_max = fake_rate_histo->GetXaxis()->GetBinCenter(nbinsx);
    double rate_value = fake_rate_histo->GetBinContent(fake_rate_histo->FindBin(min(pt, lepton_pt_max), aeta));
    double rate_error = fake_rate_histo->GetBinError(fake_rate_histo->FindBin(min(pt, lepton_pt_max), aeta));
    if (rate_value < 0) return std::make_tuple(0.0,0.0);
    return std::make_tuple(rate_value,rate_error);
}

// ------------------------------------------------------------
// Compute 1-lepton fake weight
// ------------------------------------------------------------
float fake_rate_reader_1l::GetFR_1l(double pt, double eta, double pdg, double isTight, TH2F* fake_rate_ele_, TH2F* fake_rate_muon_, TString stat){
    double p=1., f=0., pE=0., fE=0.;
    if (abs(pdg)==11){
        std::tie(p,pE) = GetRate(prompt_rate_ele_, pt, eta, -999.);
        std::tie(f,fE) = GetRate(fake_rate_ele_, pt, eta, 35.);
        if (stat=="ElUp") f+=fE;
        if (stat=="ElDown") f-=fE;
    } else if (abs(pdg)==13){
        std::tie(p,pE) = GetRate(prompt_rate_muon_, pt, eta, -999.);
        std::tie(f,fE) = GetRate(fake_rate_muon_, pt, eta, 35.);
        if (stat=="MuUp") f+=fE;
        if (stat=="MuDown") f-=fE;
    }

    double prompt_prob = 1., fake_prob = 0.;
    if (isTight==1){
        prompt_prob = p*(1.-f)/(p-f);
        fake_prob = f*(1.-p)/(p-f);
    } else{
        prompt_prob = p*f/(p-f);
        fake_prob = f*p/(p-f);
    }
    float weight = -1.;
    if (isTight==1) weight = -fake_prob;
    else weight =  fake_prob;
    return weight;
}

#endif
