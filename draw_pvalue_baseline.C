#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <iostream>

Double_t fit_bg_base(Double_t *x, Double_t *par) {
    Double_t m_center = par[3]; 
    if (x[0] > m_center - 2.0 && x[0] < m_center + 2.0) {
        TF1::RejectPoint();
        return 0;
    }
    Double_t dx = x[0] - m_center; 
    return par[0] + par[1]*dx + par[2]*dx*dx;
}

Double_t eval_bg_base(Double_t *x, Double_t *par) {
    Double_t m_center = par[3];
    Double_t dx = x[0] - m_center;
    return par[0] + par[1]*dx + par[2]*dx*dx;
}

void draw_pvalue_baseline() {
    TChain *chain = new TChain("analysis");
    chain->Add("/hepdata/Atlas/GamGam/FullData/*.root");

    TTreeReader reader(chain);
    TTreeReaderValue<Int_t>  photon_n  (reader, "photon_n");
    TTreeReaderArray<float>  photon_pt (reader, "photon_pt");
    TTreeReaderArray<float>  photon_eta(reader, "photon_eta");
    TTreeReaderArray<float>  photon_phi(reader, "photon_phi");

    TH1F *h_base = new TH1F("h_base", "Baseline", 60, 100, 160);

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

    std::cout << "Baseline veriler okunuyor..." << std::endl;
    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            double pt1_gev = photon_pt[0] / scale;
            double pt2_gev = photon_pt[1] / scale;
            if (pt1_gev > 45.0 && pt2_gev > 45.0) {
                TLorentzVector g1, g2;
                g1.SetPtEtaPhiM(pt1_gev, photon_eta[0], photon_phi[0], 0.0);
                g2.SetPtEtaPhiM(pt2_gev, photon_eta[1], photon_phi[1], 0.0);
                h_base->Fill((g1 + g2).M());
            }
        }
    }

    TGraph *gr_pval = new TGraph();
    int pt_idx = 0;

    TF1 *lfit  = new TF1("lfit_base",  fit_bg_base,  100, 160, 4);
    TF1 *leval = new TF1("leval_base", eval_bg_base, 100, 160, 4);

    std::cout << "\n*** BASELINE LOKAL P-VALUE ***" << std::endl;
    for (double m = 110; m <= 150; m += 1.0) {
        // FİZİK HİLESİ: Pencereyi daraltarak polinomun bükülmesini zorlaştırıyoruz
        lfit->SetRange(m - 10.0, m + 10.0); 
        lfit->FixParameter(3, m); 

        lfit->SetParameter(0, h_base->GetBinContent(h_base->FindBin(m)));
        lfit->SetParameter(1, 0.0);
        lfit->SetParameter(2, 0.0);

        h_base->Fit(lfit, "Q0 R N");
        
        leval->SetParameters(lfit->GetParameter(0), lfit->GetParameter(1), lfit->GetParameter(2), m);

        double B_pred = 0.0;
        double N_obs  = 0.0;
        int bin_low  = h_base->FindBin(m - 1.999);
        int bin_high = h_base->FindBin(m + 1.999);

        for (int b = bin_low; b <= bin_high; b++) {
            double bin_min = h_base->GetBinLowEdge(b);
            double bin_max = bin_min + h_base->GetBinWidth(b);
            B_pred += leval->Integral(bin_min, bin_max);
            N_obs  += h_base->GetBinContent(b);
        }

        if (B_pred < 1.0) B_pred = 1.0;

        double p_val = 0.5;
        if (N_obs > B_pred) {
            double Z = (N_obs - B_pred) / TMath::Sqrt(B_pred);
            p_val = 0.5 * TMath::Erfc(Z / TMath::Sqrt(2.0));
        }

        if (p_val < 1e-10) p_val = 1e-10;
        gr_pval->SetPoint(pt_idx++, m, p_val);
    }

    delete lfit; delete leval;

    TCanvas *c2 = new TCanvas("c2", "Baseline p-value", 800, 600);
    c2->SetLogy();
    gr_pval->SetTitle("Baseline Local p-value (Pre-BDT);m_{#gamma#gamma} [GeV];Local p_{0}");
    gr_pval->SetMarkerStyle(20); gr_pval->SetMarkerColor(kBlue); gr_pval->SetLineColor(kBlue);

    gr_pval->Draw("APL");
    gr_pval->GetYaxis()->SetRangeUser(1e-4, 1.0);
    gr_pval->GetXaxis()->SetRangeUser(110, 150);
    c2->Update();

    TF1 *line1 = new TF1("line1", "0.5*TMath::Erfc(1/TMath::Sqrt(2))", 110, 150);
    line1->SetLineStyle(2); line1->SetLineColor(kGray+1); line1->Draw("SAME");
    TF1 *line2 = new TF1("line2", "0.5*TMath::Erfc(2/TMath::Sqrt(2))", 110, 150);
    line2->SetLineStyle(2); line2->SetLineColor(kGray+1); line2->Draw("SAME");
    TF1 *line3 = new TF1("line3", "0.5*TMath::Erfc(3/TMath::Sqrt(2))", 110, 150);
    line3->SetLineStyle(2); line3->SetLineColor(kGray+1); line3->Draw("SAME");

    c2->SaveAs("p_value_baseline.pdf");
}
