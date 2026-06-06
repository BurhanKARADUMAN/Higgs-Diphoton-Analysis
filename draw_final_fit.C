#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMVA/Reader.h>
#include <iostream>

void draw_final_fit() {
    TMVA::Reader *tmvaReader = new TMVA::Reader("!Color:!Silent");
    float pt1, eta1, phi1, pt2, eta2, phi2;
    tmvaReader->AddVariable("pt1", &pt1);
    tmvaReader->AddVariable("eta1", &eta1);
    tmvaReader->AddVariable("phi1", &phi1);
    tmvaReader->AddVariable("pt2", &pt2);
    tmvaReader->AddVariable("eta2", &eta2);
    tmvaReader->AddVariable("phi2", &phi2);
    tmvaReader->BookMVA("BDT", "dataset/weights/TMVAClassification_BDT.weights.xml");

    TChain *chain = new TChain("analysis");
    chain->Add("/hepdata/Atlas/GamGam/FullData/*.root");

    TTreeReader reader(chain);
    TTreeReaderValue<Int_t> photon_n(reader, "photon_n");
    TTreeReaderArray<float> photon_pt(reader, "photon_pt");
    TTreeReaderArray<float> photon_eta(reader, "photon_eta");
    TTreeReaderArray<float> photon_phi(reader, "photon_phi");

    double optimal_cut = -0.16; 

    TH1F *h_opt = new TH1F("h_opt", "Optimized Di-photon Invariant Mass; m_{\\gamma\\gamma} [GeV]; Events / 1 GeV", 60, 100, 160);

    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            TLorentzVector g1, g2;
            g1.SetPtEtaPhiM(photon_pt[0], photon_eta[0], photon_phi[0], 0.0);
            g2.SetPtEtaPhiM(photon_pt[1], photon_eta[1], photon_phi[1], 0.0);
            
            pt1 = photon_pt[0]; eta1 = photon_eta[0]; phi1 = photon_phi[0];
            pt2 = photon_pt[1]; eta2 = photon_eta[1]; phi2 = photon_phi[1];
            
            if (tmvaReader->EvaluateMVA("BDT") > optimal_cut) {
                h_opt->Fill((g1 + g2).M());
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TF1 *fitFunc = new TF1("fitFunc", "pol2(0) + gaus(3)", 100, 160);
    fitFunc->SetParameter(0, 500); 
    fitFunc->SetParameter(1, -2);
    fitFunc->SetParameter(3, 100);
    fitFunc->SetParLimits(3, 0, 500);
    fitFunc->SetParameter(4, 125);
    fitFunc->SetParLimits(4, 122, 128);
    fitFunc->SetParameter(5, 1.5);
    fitFunc->SetParLimits(5, 0.5, 4.0);

    h_opt->SetMarkerStyle(20);
    h_opt->SetLineColor(kBlack);
    h_opt->Fit(fitFunc, "S");
    h_opt->Draw("E1");
    fitFunc->Draw("SAME");
    c1->SaveAs("higgs_mass_spectrum_optimized.pdf");
    std::cout << "2. Grafik guncellendi: higgs_mass_spectrum_optimized.pdf" << std::endl;
}
