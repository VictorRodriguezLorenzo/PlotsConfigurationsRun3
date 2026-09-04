#include <iostream>
#include <string>
#include <map>
#include <filesystem>
#include <cmath>
#include <stdexcept>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

struct CampaignConfig {
    std::string sample;
    std::string algo;
    std::map<std::string, double> wp;
};

std::map<std::string, CampaignConfig> getCampaigns() {
    const std::string base = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/";

    return {
        {"Full2022v12",     {base+"Summer22_130x_nAODv12_Full2022v12/MCl2loose2022v12__MCCorr2022v12JetScaling__l2tight/nanoLatino_TTTo2L2Nu__part*.root",                     "PNetB",      {{"loose",0.0470},{"medium",0.2450},{"tight",0.6734},{"xtight",0.7862},{"xxtight",0.9610}}}},
        {"Full2022EEv12",   {base+"Summer22EE_130x_nAODv12_Full2022v12/MCl2loose2022EEv12__MCCorr2022EEv12JetScaling__l2tight/nanoLatino_TTTo2L2Nu__part*.root",             "PNetB",      {{"loose",0.0499},{"medium",0.2605},{"tight",0.6915},{"xtight",0.8033},{"xxtight",0.9664}}}},
        {"Full2023v12",     {base+"Summer23_130x_nAODv12_Full2023v12/MCl2loose2023v12__MCCorr2023v12JetScaling__l2tight/nanoLatino_TTTo2L2Nu__part*.root",                     "PNetB",      {{"loose",0.0358},{"medium",0.1919}}}},
        {"Full2023BPixv12", {base+"Summer23BPix_130x_nAODv12_Full2023BPixv12/MCl2loose2023BPixv12__MCCorr2023BPixv12JetScaling__l2tight/nanoLatino_TTTo2L2Nu__part*.root",     "PNetB",      {{"loose",0.0359},{"medium",0.1919}}}},
        {"Full2024v15",     {base+"Summer24_150x_nAODv15_Full2024v15/MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight/nanoLatino_TTTo2L2Nu__part*.root",               "UParTAK4B",  {{"loose",0.0246},{"medium",0.1272},{"tight",0.4648},{"xtight",0.6298},{"xxtight",0.9739}}}}
    };
}


void bTagEff(std::string campaign, std::string WP="medium") {

    TH1::SetDefaultSumw2(true);

    auto campaigns = getCampaigns();
    if (!campaigns.count(campaign)) throw std::runtime_error("Unknown campaign: "+campaign);

    const auto &cfg = campaigns.at(campaign);
    if (!cfg.wp.count(WP)) throw std::runtime_error("WP '"+WP+"' not defined for campaign '"+campaign+"'");

    const std::string algo = cfg.algo;
    const double wp = cfg.wp.at(WP);

    TString samples = cfg.sample.c_str();
    TString btag = ("Jet_btag"+algo).c_str();
    TString fname = Form("bTagEff_%s_ttbar_%s_%s.root", campaign.c_str(), algo.c_str(), WP.c_str());

    std::filesystem::create_directories("efficiencies");

    std::cout << "\nCampaign: " << campaign
              << "\nAlgorithm: " << algo
              << "\nWP: " << WP
              << "\nThreshold: " << wp
              << "\nInput: " << samples
              << "\nOutput: " << fname << "\n" << std::endl;

    TChain *Events = new TChain("Events");
    Events->Add(samples);

    const Long64_t entries = Events->GetEntries();
    if (entries == 0) throw std::runtime_error("No events found for "+campaign);

    std::cout << "Entries: " << entries << std::endl;

    // Branches
    Events->SetBranchStatus("*",0);
    Events->SetBranchStatus("XSWeight",1);
    Events->SetBranchStatus("nCleanJet",1);
    Events->SetBranchStatus("CleanJet_pt",1);
    Events->SetBranchStatus("CleanJet_eta",1);
    Events->SetBranchStatus("CleanJet_jetIdx",1);
    Events->SetBranchStatus("Jet_hadronFlavour",1);
    Events->SetBranchStatus(btag,1);
    Events->SetBranchStatus("nLepton",1);
    Events->SetBranchStatus("Lepton_pt",1);
    Events->SetBranchStatus("Lepton_eta",1);
    Events->SetBranchStatus("Lepton_pdgId",1);
    Events->SetBranchStatus("mll",1);

    // Variables
    double XSWeight;
    Int_t nCleanJet, nLepton;
    float CleanJet_pt[100], CleanJet_eta[100], Jet_btag[100];
    ULong64_t CleanJet_jetIdx[100];
    char Jet_hadronFlavour[100];
    float Lepton_pt[100], Lepton_eta[100];
    int Lepton_pdgId[100];
    Double_t mll;

    // Addresses
    Events->SetBranchAddress("XSWeight",&XSWeight);
    Events->SetBranchAddress("nCleanJet",&nCleanJet);
    Events->SetBranchAddress("CleanJet_pt",CleanJet_pt);
    Events->SetBranchAddress("CleanJet_eta",CleanJet_eta);
    Events->SetBranchAddress("CleanJet_jetIdx",CleanJet_jetIdx);
    Events->SetBranchAddress("Jet_hadronFlavour",Jet_hadronFlavour);
    Events->SetBranchAddress(btag,Jet_btag);
    Events->SetBranchAddress("nLepton",&nLepton);
    Events->SetBranchAddress("Lepton_pt",Lepton_pt);
    Events->SetBranchAddress("Lepton_eta",Lepton_eta);
    Events->SetBranchAddress("Lepton_pdgId",Lepton_pdgId);
    Events->SetBranchAddress("mll",&mll);

    // Binning
    Float_t ptbins[]  = {20,30,50,70,100,140,200,300,600,1000};
    Float_t etabins[] = {-2.5,-1.479,1.479,2.5};

    TH2F *bjet_den = new TH2F("bjet_den","bjet_den",9,ptbins,3,etabins);
    TH2F *bjet_num = new TH2F("bjet_num","bjet_num",9,ptbins,3,etabins);
    TH2F *cjet_den = new TH2F("cjet_den","cjet_den",9,ptbins,3,etabins);
    TH2F *cjet_num = new TH2F("cjet_num","cjet_num",9,ptbins,3,etabins);
    TH2F *ljet_den = new TH2F("ljet_den","ljet_den",9,ptbins,3,etabins);
    TH2F *ljet_num = new TH2F("ljet_num","ljet_num",9,ptbins,3,etabins);

    for (auto h : {bjet_den,bjet_num,cjet_den,cjet_num,ljet_den,ljet_num}) {
        h->GetXaxis()->SetTitle("p_{T} [GeV]");
        h->GetYaxis()->SetTitle("#eta");
    }

    // Event loop
    for (Long64_t i=0; i<entries; ++i) {

        Events->GetEntry(i);

        if (i%1000000==0) std::cout << "Processing " << i << "/" << entries << " (" << 100.*i/entries << "%)" << std::endl;

        if (nLepton<2) continue;

        const int id0 = std::abs(Lepton_pdgId[0]);
        const int id1 = std::abs(Lepton_pdgId[1]);

        if (!((id0==11 || id0==13) && (id1==11 || id1==13))) continue;
        if (Lepton_pt[0]<=25. || Lepton_pt[1]<=20.) continue;
        if (!((id0==11 && std::abs(Lepton_eta[0])<2.5) || (id0==13 && std::abs(Lepton_eta[0])<2.4))) continue;
        if (!((id1==11 && std::abs(Lepton_eta[1])<2.5) || (id1==13 && std::abs(Lepton_eta[1])<2.4))) continue;
        if (mll<=20.) continue;

        for (int jet=0; jet<nCleanJet; ++jet) {

            if (CleanJet_pt[jet]<=30. || std::abs(CleanJet_eta[jet])>=2.5) continue;

            const auto jetIdx = CleanJet_jetIdx[jet];
            const int flavour = std::abs(static_cast<int>(Jet_hadronFlavour[jetIdx]));
            const bool tagged = Jet_btag[jetIdx]>wp;

            if (flavour==5) {
                bjet_den->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
                if (tagged) bjet_num->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
            }
            else if (flavour==4) {
                cjet_den->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
                if (tagged) cjet_num->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
            }
            else if (flavour==0) {
                ljet_den->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
                if (tagged) ljet_num->Fill(CleanJet_pt[jet],CleanJet_eta[jet],XSWeight);
            }
        }
    }

    // Efficiencies
    TFile *outfile = new TFile(fname,"RECREATE");

    TH2F *bjet_eff = (TH2F*)bjet_num->Clone("bjet_eff"); bjet_eff->Divide(bjet_den);
    TH2F *cjet_eff = (TH2F*)cjet_num->Clone("cjet_eff"); cjet_eff->Divide(cjet_den);
    TH2F *ljet_eff = (TH2F*)ljet_num->Clone("ljet_eff"); ljet_eff->Divide(ljet_den);

    bjet_den->Write(); bjet_num->Write(); bjet_eff->Write();
    cjet_den->Write(); cjet_num->Write(); cjet_eff->Write();
    ljet_den->Write(); ljet_num->Write(); ljet_eff->Write();

    // Plots
    gStyle->SetPaintTextFormat("2.3f");
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c","c",850,700);
    c->SetLogx();

    TString prefix = Form("efficiencies/%s_%s_%s",campaign.c_str(),algo.c_str(),WP.c_str());

    bjet_eff->Draw("COLZ E TEXT"); c->SaveAs(prefix+"_bjet.png");
    cjet_eff->Draw("COLZ E TEXT"); c->SaveAs(prefix+"_cjet.png");
    ljet_eff->Draw("COLZ E TEXT"); c->SaveAs(prefix+"_lightjet.png");

    outfile->Close();

    std::cout << "Written: " << fname << std::endl;
}
