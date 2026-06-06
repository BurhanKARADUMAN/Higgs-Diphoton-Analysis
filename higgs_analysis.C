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
#include <TMath.h>
#include <iostream>

// Körleme Fonksiyonu (121-129 arası kapalı)
Double_t exclude_sig(Double_t *x, Double_t *par) {
    if (x[0] > 121.0 && x[0] < 129.0) {
        TF1::RejectPoint();
        return 0;
    }
    // Daha katı bir pol3 (Yılan gibi kıvrılmaz, köprü kurar)
    return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
}

void higgs_analysis() {
    TChain *chain = new TChain("analysis");
    chain->Add("/hepdata/Atlas/GamGam/FullData/*.root");

    TTreeReader reader(chain);
    TTreeReaderValue<Int_t>  photon_n  (reader, "photon_n");
    TTreeReaderArray<float>  photon_pt (reader, "photon_pt");
    TTreeReaderArray<float>  photon_eta(reader, "photon_eta");
    TTreeReaderArray<float>  photon_phi(reader, "photon_phi");

    TH1F *h_myy = new TH1F("h_myy", "Di-photon Invariant Mass;m_{#gamma#gamma} [GeV];Events / 1 GeV", 60, 100, 160);

    double scale = 1.0;
    bool scale_determined = false;
    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            scale = (photon_pt[0] > 1000.0) ? 1000.0 : 1.0;
            scale_determined = true; break;
        }
    }
    if (!scale_determined) return;
    reader.Restart();

    std::cout << "Reading baseline events..." << std::endl;
    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            double pt1_gev = photon_pt[0] / scale;
            double pt2_gev = photon_pt[1] / scale;
            if (pt1_gev > 45.0 && pt2_gev > 45.0) {
                TLorentzVector g1, g2;
                g1.SetPtEtaPhiM(pt1_gev, photon_eta[0], photon_phi[0], 0.0);
                g2.SetPtEtaPhiM(pt2_gev, photon_eta[1], photon_phi[1], 0.0);
                h_myy->Fill((g1 + g2).M());
            }
        }
    }

    // ADIM 1: Sadece yan bantlara bakarak Arkaplanı öğren (110-150 arası daha stabildir)
    TF1 *bg_pre = new TF1("bg_pre", exclude_sig, 110, 150, 4);
    h_myy->Fit(bg_pre, "Q0 R");

    // ADIM 2: Asıl Fit. Pol3 + Gaus. 
    TF1 *fitFunc = new TF1("fitFunc", "pol3(0) + gaus(4)", 110, 150);
    
    // FİZİK HİLESİ: Arkaplanı kilitliyoruz ki sinyali yutamasın!
    for(int i=0; i<4; i++) {
        fitFunc->FixParameter(i, bg_pre->GetParameter(i));
    }
    
    // Gaus'un Higgs'i bulmasına izin ver
    fitFunc->SetParameter(4, 1000.0);        
    fitFunc->SetParLimits(4, 0.0, 50000.0);

    fitFunc->FixParameter(5, 125.0); // Higgs'in yeri KİLİTLENDİ
    fitFunc->FixParameter(6, 1.5);   // Dağılım genişliği KİLİTLENDİ
    // Kilitli arkaplanla fiti ateşle
    h_myy->Fit(fitFunc, "R S");

    double bw = h_myy->GetBinWidth(1);
    TF1 *gaus_only = new TF1("gaus_only", "gaus", 123, 127);
    gaus_only->SetParameters(fitFunc->GetParameter(4), fitFunc->GetParameter(5), fitFunc->GetParameter(6));
    double S0 = gaus_only->Integral(123, 127) / bw;

    TF1 *bg_only = new TF1("bg_only", "pol3", 123, 127);
    bg_only->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2), fitFunc->GetParameter(3));
    double B0 = bg_only->Integral(123, 127) / bw;

    std::cout << "\nS0 (signal in 123-127 GeV)     = " << S0 << std::endl;
    std::cout << "B0 (background in 123-127 GeV) = " << B0 << std::endl;
    if (B0 > 0) std::cout << "Significance S/sqrt(B)         = " << S0 / TMath::Sqrt(B0) << " sigma" << std::endl;

    TCanvas *c1 = new TCanvas("c1", "Baseline Higgs", 800, 600);
    h_myy->SetMarkerStyle(20); h_myy->Draw("E1");
    fitFunc->SetLineColor(kRed); fitFunc->Draw("SAME");
    bg_only->SetRange(100, 160); bg_only->SetLineColor(kBlue); bg_only->SetLineStyle(2); bg_only->Draw("SAME");
    
    c1->SaveAs("higgs_mass_spectrum_fitted.pdf");

    delete gaus_only; delete bg_only; delete fitFunc; delete bg_pre; delete h_myy; delete chain;
}
