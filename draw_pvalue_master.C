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
#include <TMVA/Reader.h>
#include <TMath.h>
#include <TSystem.h>
#include <iostream>

// FUNCTION 1: For fitting only (with blind region)
Double_t fit_bg(Double_t *x, Double_t *par) {
    Double_t m_center = par[3];
    if (x[0] > m_center - 2.0 && x[0] < m_center + 2.0) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// FUNCTION 2: For integration only (pure polynomial, no RejectPoint)
Double_t eval_bg(Double_t *x, Double_t *par) {
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

void draw_pvalue_master() {

    // Weight file safety check
    if (gSystem->AccessPathName("dataset/weights/TMVAClassification_BDT.weights.xml")) {
        std::cerr << "CRITICAL ERROR: BDT Weight file not found!" << std::endl;
        return;
    }

    TMVA::Reader *tmvaReader = new TMVA::Reader("!Color:!Silent");
    float pt1, eta1, phi1, pt2, eta2, phi2;
    // IMPORTANT: These names must exactly match what was used during training
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
    TTreeReaderValue<Int_t>    photon_n  (reader, "photon_n");
    TTreeReaderArray<float>    photon_pt (reader, "photon_pt");
    TTreeReaderArray<float>    photon_eta(reader, "photon_eta");
    TTreeReaderArray<float>    photon_phi(reader, "photon_phi");

    double optimal_cut = -0.16;
    TH1F *h_opt = new TH1F("h_opt", "Optimized Di-photon Invariant Mass", 60, 100, 160);

    // --- MeV/GeV unit detection ---
    // Read the first valid event to determine scale, BEFORE the main loop
    double scale = 1.0;
    bool scale_determined = false;
    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {
            scale = (photon_pt[0] > 1000.0) ? 1000.0 : 1.0;
            std::cout << "DEBUG: First photon pT = " << photon_pt[0]
                      << " -> Scale factor = " << scale
                      << " (Unit: " << (scale > 1.0 ? "MeV" : "GeV") << ")" << std::endl;
            scale_determined = true;
            break;
        }
    }
    if (!scale_determined) {
        std::cerr << "ERROR: Could not determine unit scale. No valid events found." << std::endl;
        delete tmvaReader;
        delete chain;
        return;
    }

    // CRITICAL NOTE: The TMVA BDT must have been trained with the same unit as we feed here.
    // We always feed GeV (after scaling). If training used MeV, retrain with GeV values.
    std::cout << "WARNING: Ensure BDT was trained with pT in GeV (post-scale). "
              << "If trained in MeV, results will be incorrect." << std::endl;

    // Reset reader to beginning for the main loop
    reader.Restart();

    std::cout << "Reading data and applying BDT..." << std::endl;
    while (reader.Next()) {
        if (*photon_n >= 2 && photon_pt.GetSize() >= 2) {

            double pt1_gev = photon_pt[0] / scale;
            double pt2_gev = photon_pt[1] / scale;

            // Kinematic pre-selection (same as baseline — required for fair comparison)
            if (pt1_gev > 45.0 && pt2_gev > 45.0) {

                TLorentzVector g1, g2;
                g1.SetPtEtaPhiM(pt1_gev, photon_eta[0], photon_phi[0], 0.0);
                g2.SetPtEtaPhiM(pt2_gev, photon_eta[1], photon_phi[1], 0.0);

                // Feed GeV values into TMVA (must match training scale)
                pt1  = pt1_gev;  eta1 = photon_eta[0]; phi1 = photon_phi[0];
                pt2  = pt2_gev;  eta2 = photon_eta[1]; phi2 = photon_phi[1];

                if (tmvaReader->EvaluateMVA("BDT") > optimal_cut) {
                    h_opt->Fill((g1 + g2).M());
                }
            }
        }
    }

    std::cout << "Total entries after BDT cut: " << h_opt->GetEntries() << std::endl;

    // --- Local p-value calculation ---
    TGraph *gr_pval = new TGraph();
    int pt_idx = 0;

    TF1 *lfit  = new TF1("lfit",  fit_bg,  100, 160, 4); // blinded fit function
    TF1 *leval = new TF1("leval", eval_bg, 100, 160, 3); // pure polynomial for integration

    std::cout << "\n*** COMPUTING LOCAL P-VALUES ***" << std::endl;
    for (double m = 110; m <= 150; m += 1.0) {

        lfit->SetRange(m - 12.0, m + 12.0);
        lfit->FixParameter(3, m);
        lfit->SetParameter(0, h_opt->GetBinContent(h_opt->FindBin(m)));
        lfit->SetParameter(1, -10.0);
        lfit->SetParameter(2, 0.0);

        // Fit with quality check via NDF
        h_opt->Fit(lfit, "Q0 R N");
        if (lfit->GetNDF() <= 0) {
            std::cerr << "WARNING: Fit invalid (NDF<=0) at m=" << m << ". Assigning p=0.5." << std::endl;
            gr_pval->SetPoint(pt_idx++, m, 0.5);
            continue;
        }

        // Transfer fitted parameters to the pure evaluation function
        leval->SetParameters(lfit->GetParameter(0),
                             lfit->GetParameter(1),
                             lfit->GetParameter(2));

        double B_pred = 0.0;
        double N_obs  = 0.0;

        int bin_low  = h_opt->FindBin(m - 1.999);
        int bin_high = h_opt->FindBin(m + 1.999);

        for (int b = bin_low; b <= bin_high; b++) {
            double bin_min = h_opt->GetBinLowEdge(b);
            double bin_max = bin_min + h_opt->GetBinWidth(b);
            B_pred += leval->Integral(bin_min, bin_max); // correct bin-width-aware integration
            N_obs  += h_opt->GetBinContent(b);
        }

        // Floor to prevent sqrt(0) instability
        if (B_pred < 1.0) B_pred = 1.0;

        double p_val = 0.5;
        if (N_obs > B_pred) {
            double Z = (N_obs - B_pred) / TMath::Sqrt(B_pred);
            p_val = 0.5 * TMath::Erfc(Z / TMath::Sqrt(2.0));
        }

        if (p_val < 1e-10) p_val = 1e-10;
        gr_pval->SetPoint(pt_idx++, m, p_val);

        std::cout << "  m=" << m << " GeV | N_obs=" << N_obs
                  << " | B_pred=" << B_pred
                  << " | p0=" << p_val << std::endl;
    }

    delete lfit;
    delete leval;
    delete tmvaReader;

    // --- Plotting ---
    TCanvas *c2 = new TCanvas("c2", "Local p-value (BDT Optimized)", 800, 600);
    c2->SetLogy();

    gr_pval->SetTitle("Local p-value (BDT Optimized);m_{#gamma#gamma} [GeV];Local p_{0}");
    gr_pval->SetMarkerStyle(20);
    gr_pval->SetMarkerColor(kRed);
    gr_pval->SetLineColor(kRed);
    gr_pval->Draw("APL");

    // Axis range MUST be set after Draw() for TGraph
    gr_pval->GetYaxis()->SetRangeUser(1e-4, 1.0);
    gr_pval->GetXaxis()->SetRangeUser(110, 150);
    c2->Update(); // force axis refresh after range change

    // Reference sigma lines
    TF1 *line1 = new TF1("line1", "0.5*TMath::Erfc(1/TMath::Sqrt(2))", 110, 150);
    line1->SetLineStyle(2); line1->SetLineColor(kGray+1); line1->Draw("SAME");

    TF1 *line2 = new TF1("line2", "0.5*TMath::Erfc(2/TMath::Sqrt(2))", 110, 150);
    line2->SetLineStyle(2); line2->SetLineColor(kGray+1); line2->Draw("SAME");

    TF1 *line3 = new TF1("line3", "0.5*TMath::Erfc(3/TMath::Sqrt(2))", 110, 150);
    line3->SetLineStyle(2); line3->SetLineColor(kGray+1); line3->Draw("SAME");

    c2->SaveAs("p_value_master.pdf");
    std::cout << "\nSaved: p_value_master.pdf" << std::endl;

    // Cleanup
    delete h_opt;
    delete chain;
    delete gr_pval;
}
