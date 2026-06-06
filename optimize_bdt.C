#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMVA/Reader.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void optimize_bdt() {
    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
    float pt1, eta1, phi1, pt2, eta2, phi2;
    
    reader->AddVariable("pt1", &pt1);
    reader->AddVariable("eta1", &eta1);
    reader->AddVariable("phi1", &phi1);
    reader->AddVariable("pt2", &pt2);
    reader->AddVariable("eta2", &eta2);
    reader->AddVariable("phi2", &phi2);
    reader->BookMVA("BDT", "dataset/weights/TMVAClassification_BDT.weights.xml");

    TChain *chain = new TChain("analysis");
    chain->Add("/hepdata/Atlas/GamGam/FullData/*.root");

    TTreeReader treeReader(chain);
    TTreeReaderValue<Int_t> photon_n(treeReader, "photon_n");
    TTreeReaderArray<float> photon_pt(treeReader, "photon_pt");
    TTreeReaderArray<float> photon_eta(treeReader, "photon_eta");
    TTreeReaderArray<float> photon_phi(treeReader, "photon_phi");

    struct EventData { double mass; double bdt_score; };
    std::vector<EventData> events;

    std::cout << "Tum veriler okunuyor ve BDT skorlari hesaplaniyor..." << std::endl;
    
    while (treeReader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            TLorentzVector gamma1, gamma2;
            
            gamma1.SetPtEtaPhiM(photon_pt[0], photon_eta[0], photon_phi[0], 0.0);
            gamma2.SetPtEtaPhiM(photon_pt[1], photon_eta[1], photon_phi[1], 0.0);
            
            pt1 = photon_pt[0]; eta1 = photon_eta[0]; phi1 = photon_phi[0];
            pt2 = photon_pt[1]; eta2 = photon_eta[1]; phi2 = photon_phi[1];
            
            double m_yy = (gamma1 + gamma2).M();
            double bdt_score = reader->EvaluateMVA("BDT");
            
            events.push_back({m_yy, bdt_score});
        }
    }

    double best_cut = -999.0;
    double max_significance = -1.0;
    double best_S = 0.0, best_B = 0.0;

    std::cout << "Optimizasyon dongusu basliyor..." << std::endl;
    
    for (double cut = -0.2; cut <= 0.4; cut += 0.04) {
        
        TH1F *h_temp = new TH1F(Form("h_%.2f", cut), "Temp", 60, 100, 160);
        
        for (const auto& ev : events) {
            if (ev.bdt_score > cut) {
                h_temp->Fill(ev.mass);
            }
        }

        TF1 *fitFunc = new TF1(Form("fit_%.2f", cut), "pol2(0) + gaus(3)", 100, 160);
        
        fitFunc->SetParameter(0, 500);
        fitFunc->SetParameter(1, -1);
        fitFunc->SetParameter(3, 20); 
        fitFunc->SetParLimits(3, 0, 500); 
        fitFunc->SetParameter(4, 125); 
        fitFunc->SetParLimits(4, 122, 128);
        fitFunc->SetParameter(5, 1.5); 
        fitFunc->SetParLimits(5, 0.5, 4.0);
        
        h_temp->Fit(fitFunc, "Q0"); 

        TF1 *bgFunc = new TF1("bg", "pol2(0)", 100, 160);
        bgFunc->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
        TF1 *sigFunc = new TF1("sig", "gaus(0)", 100, 160);
        sigFunc->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5));

        double binW = h_temp->GetBinWidth(1);
        double S = sigFunc->Integral(120, 130) / binW;
        double B = bgFunc->Integral(120, 130) / binW;
        
        if (B <= 0 || S <= 0) { 
            delete h_temp; 
            continue; 
        }

        double significance = S / TMath::Sqrt(B);
        std::cout << "BDT Cut: " << cut << " | S: " << S << " | B: " << B << " | Signif: " << significance << std::endl;

        if (significance > max_significance) {
            max_significance = significance;
            best_cut = cut;
            best_S = S;
            best_B = B;
        }
        
        delete h_temp;
    }

    std::cout << "\n==========================================" << std::endl;
    std::cout << "OPTIMAL BDT KESIMI (TUM VERILER ILE) BULUNDU: > " << best_cut << std::endl;
    std::cout << "Maksimum Anlamlilik (S/sqrt(B)): " << max_significance << std::endl;
    std::cout << "Optimal Sinyal (S): " << best_S << std::endl;
    std::cout << "Optimal Arkaplan (B): " << best_B << std::endl;
    std::cout << "==========================================\n" << std::endl;
}
