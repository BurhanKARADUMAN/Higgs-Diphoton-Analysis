// higgs_diphoton_analysis.C
// Complete rewrite using TTreeReader for RVec compatibility

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <vector>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>

// ============================================================
// CONSTANTS
// ============================================================
const double MYY_LO = 100., MYY_HI = 160.;
const double SR_LO  = 120., SR_HI  = 130.;
const double SB1_LO = 105., SB1_HI = 120.;
const double SB2_LO = 130., SB2_HI = 160.;
const int    NBINS  = 60;
const double BIN_W  = (MYY_HI - MYY_LO) / NBINS;

// ============================================================
// ATLAS STYLE
// ============================================================
void SetAtlasStyle() {
    TStyle* s = new TStyle("ATLAS","ATLAS style");
    s->SetOptStat(0); s->SetOptTitle(0); s->SetOptFit(0);
    s->SetFrameBorderMode(0); s->SetCanvasBorderMode(0);
    s->SetPadColor(0); s->SetCanvasColor(0);
    s->SetPadTopMargin(0.06);  s->SetPadRightMargin(0.08);
    s->SetPadBottomMargin(0.16); s->SetPadLeftMargin(0.16);
    s->SetLabelSize(0.04,"xyz"); s->SetTitleSize(0.05,"xyz");
    s->SetMarkerStyle(20); s->SetMarkerSize(0.9);
    s->SetHistLineWidth(2); s->SetErrorX(0.0001);
    gROOT->SetStyle("ATLAS"); gROOT->ForceStyle();
}

// ============================================================
// FIT FUNCTIONS
// ============================================================
double SigPlusBkg(double* x, double* p) {
    return p[0]*TMath::Exp(p[1]*x[0])
         + p[2]*TMath::Gaus(x[0], p[3], p[4], false);
}
double BkgOnly(double* x, double* p) {
    return p[0]*TMath::Exp(p[1]*x[0]);
}
double LocalBkg_fit(double* x, double* p) {
    if (x[0] > p[3]-2.0 && x[0] < p[3]+2.0) {
        TF1::RejectPoint(); return 0;
    }
    return p[0] + p[1]*x[0] + p[2]*x[0]*x[0];
}
double LocalBkg_eval(double* x, double* p) {
    return p[0] + p[1]*x[0] + p[2]*x[0]*x[0];
}

// ============================================================
// EXTRACT S AND B FROM GLOBAL FIT
// ============================================================
void ExtractSB(TF1* f, double& S, double& B, double& Z) {
    TF1 fb("_fb_", BkgOnly, MYY_LO, MYY_HI, 2);
    fb.SetParameters(f->GetParameter(0), f->GetParameter(1));
    B = fb.Integral(SR_LO, SR_HI) / BIN_W;

    TF1 fs("_fs_","[0]*TMath::Gaus(x,[1],[2],false)",
            MYY_LO, MYY_HI);
    fs.SetParameters(f->GetParameter(2),
                     f->GetParameter(3),
                     f->GetParameter(4));
    S = std::max(0.0, fs.Integral(SR_LO, SR_HI) / BIN_W);
    Z = (B > 0) ? S/TMath::Sqrt(B) : 0;
}

// ============================================================
// SLIDING WINDOW P-VALUE
// ============================================================
TGraph* SlidingWindowPvalue(TH1F* h, const char* name,
                             double window=12.0,
                             double blind=2.0) {
    TGraph* gr = new TGraph();
    gr->SetName(name);
    int idx = 0;

    TF1* lfit  = new TF1(Form("lfit_%s",name),
                          LocalBkg_fit,  MYY_LO, MYY_HI, 4);
    TF1* leval = new TF1(Form("leval_%s",name),
                          LocalBkg_eval, MYY_LO, MYY_HI, 3);

    for (double m = 110; m <= 150; m += 1.0) {
        lfit->SetRange(m-window, m+window);
        lfit->FixParameter(3, m);
        double c_lo = h->GetBinContent(h->FindBin(m-window*0.6));
        double c_hi = h->GetBinContent(h->FindBin(m+window*0.6));
        lfit->SetParameter(0, h->GetBinContent(h->FindBin(m)));
        lfit->SetParameter(1, (c_hi-c_lo)/(window*1.2));
        lfit->SetParameter(2, 0.0);

        h->Fit(lfit, "Q0 R N");
        if (lfit->GetNDF() <= 0) {
            lfit->SetParameter(1, -50.0);
            h->Fit(lfit, "Q0 R N M");
        }
        if (lfit->GetNDF() <= 0) {
            gr->SetPoint(idx++, m, 0.5); continue;
        }

        leval->SetParameters(lfit->GetParameter(0),
                             lfit->GetParameter(1),
                             lfit->GetParameter(2));

        double B=0, N=0;
        int bl = h->FindBin(m-blind+0.001);
        int bh = h->FindBin(m+blind-0.001);
        for (int b=bl; b<=bh; b++) {
            double blo = h->GetBinLowEdge(b);
            B += leval->Integral(blo, blo+h->GetBinWidth(b));
            N += h->GetBinContent(b);
        }

        if (B<=0 || B>1e9) { gr->SetPoint(idx++,m,0.5); continue; }
        if (B<1.0) B=1.0;

        double p=0.5;
        if (N>B) {
            double Z=(N-B)/TMath::Sqrt(B);
            p=0.5*TMath::Erfc(Z/TMath::Sqrt(2.0));
        }
        if (p<1e-12) p=1e-12;
        gr->SetPoint(idx++, m, p);
        printf("  m=%.0f N=%.0f B=%.1f Z=%.2f p=%.2e\n",
               m,N,B,(N>B?(N-B)/TMath::Sqrt(B):0.),p);
    }
    delete lfit; delete leval;
    return gr;
}

// ============================================================
// DRAW P-VALUE GRAPH
// ============================================================
void DrawPvalGraph(TGraph* gr, const char* label,
                   const char* fname, int color=kRed+1) {
    TCanvas* c = new TCanvas(fname,fname,800,600);
    c->SetLogy();
    c->SetLeftMargin(0.14); c->SetRightMargin(0.10);
    gr->SetLineColor(color); gr->SetLineWidth(2);
    gr->SetMarkerColor(color); gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetTitle(";m_{#gamma#gamma} [GeV];Local p_{0}");
    gr->Draw("APL");
    gr->GetXaxis()->SetRangeUser(110,150);
    gr->GetYaxis()->SetRangeUser(1e-12,2.0);
    c->Update();

    const double sigLev[]={1,2,3,4,5};
    const int    sigCol[]={kGray+1,kGray+2,kOrange+1,kOrange+2,kRed};
    TLatex lt; lt.SetTextSize(0.030);
    for (int k=0;k<5;k++) {
        double p=0.5*TMath::Erfc(sigLev[k]/TMath::Sqrt2());
        TLine* l=new TLine(110,p,150,p);
        l->SetLineColor(sigCol[k]); l->SetLineStyle(2);
        l->SetLineWidth(1); l->Draw();
        lt.SetTextColor(sigCol[k]);
        lt.DrawLatex(150.3,p,Form("%d#sigma",(int)sigLev[k]));
    }
    TLatex tit; tit.SetNDC(); tit.SetTextSize(0.034);
    tit.DrawLatex(0.16,0.88,label);
    c->Update();
    c->SaveAs(Form("plots/%s.pdf",fname));
    c->SaveAs(Form("plots/%s.png",fname));
    printf("[PLOT] %s saved.\n",fname);
    delete c;
}

// ============================================================
// BDT VARIABLES STRUCT
// ============================================================
struct BDTVars {
    float pt1n, pt2n, aeta1, aeta2, adeta,
          adphi, costh, iso1, iso2, ptiso1, ptiso2;
};

// ============================================================
// CHAIN BUILDER
// ============================================================
TChain* MakeChain(const char* dir, const char* label) {
    std::vector<std::string> names =
        {"mini","analysis","GamGam","CollectionTree"};
    for (auto& n : names) {
        TChain* ch = new TChain(n.c_str());
        ch->Add(Form("%s*.root",dir));
        if (ch->GetEntries()>0) {
            printf("[INFO] %-12s %lld entries (tree: %s)\n",
                   label, ch->GetEntries(), n.c_str());
            return ch;
        }
        delete ch;
    }
    printf("[WARN] %-12s 0 entries\n", label);
    return new TChain("mini");
}

// ============================================================
// PRESELECTION USING TTreeReader
// Fills histogram and optionally a BDT vars vector
// isMC=true: reads MC weight branches
// fillVars=true: fills vOut for BDT training/scoring
// sbOnly=true: only fills sideband events (for bkg tree)
// srOnly=true: only fills signal region events (for sig tree)
// ============================================================
struct EvScore { double myy; float score; };

void ProcessChain(TChain* chain, bool isMC,
                  TH1F* hOut,           // fill invariant mass (can be nullptr)
                  TTree* tOut,          // fill output tree (can be nullptr)
                  std::vector<EvScore>* evOut, // fill scored events (nullptr if not scoring)
                  TMVA::Reader* scorer, // TMVA reader (nullptr if not scoring)
                  bool sbOnly,          // sideband only (for bkg tree)
                  bool srOnly,          // signal region only (for sig tree)
                  // Tree branches for output tree:
                  float* t_pt1n, float* t_pt2n,
                  float* t_ae1,  float* t_ae2,
                  float* t_ade,  float* t_adp,
                  float* t_cs,   float* t_i1,
                  float* t_i2,   float* t_pi1,
                  float* t_pi2,  float* t_w,
                  double* t_myy,
                  // TMVA scorer variables:
                  float* r_pt1n, float* r_pt2n,
                  float* r_ae1,  float* r_ae2,
                  float* r_ade,  float* r_adp,
                  float* r_cs,   float* r_i1,
                  float* r_i2,   float* r_pi1,
                  float* r_pi2,
                  double bestT=0.)
{
    TTreeReader reader(chain);

    // --- Kinematic branches (always present) ---
    TTreeReaderArray<float> photon_pt (reader, "photon_pt");
    TTreeReaderArray<float> photon_eta(reader, "photon_eta");
    TTreeReaderArray<float> photon_phi(reader, "photon_phi");

    // --- Optional branches (check existence first) ---
    bool hasE       = chain->GetBranch("photon_E")        != nullptr;
    bool hasTight   = chain->GetBranch("photon_isTightID") != nullptr;
    bool hasEtcone  = chain->GetBranch("photon_etcone20") != nullptr;
    bool hasPtcone  = chain->GetBranch("photon_ptcone30") != nullptr;
    bool hasTrig    = chain->GetBranch("trigP")           != nullptr;
    bool hasMCW     = chain->GetBranch("mcWeight")        != nullptr && isMC;
    bool hasSFPh    = chain->GetBranch("scaleFactor_PHOTON") != nullptr && isMC;
    bool hasSFTrig  = chain->GetBranch("scaleFactor_PhotonTRIGGER") != nullptr && isMC;

    // Use unique_ptr-style optional readers
    // We use raw pointers and check hasXxx before access
    TTreeReaderValue<Bool_t>* rTrig   = hasTrig   ?
        new TTreeReaderValue<Bool_t>(reader,"trigP") : nullptr;
    TTreeReaderArray<float>*  rE      = hasE       ?
        new TTreeReaderArray<float>(reader,"photon_E") : nullptr;
    TTreeReaderArray<bool>*   rTight  = hasTight   ?
        new TTreeReaderArray<bool>(reader,"photon_isTightID") : nullptr;
    TTreeReaderArray<float>*  rEtc    = hasEtcone  ?
        new TTreeReaderArray<float>(reader,"photon_etcone20") : nullptr;
    TTreeReaderArray<float>*  rPtc    = hasPtcone  ?
        new TTreeReaderArray<float>(reader,"photon_ptcone30") : nullptr;
    TTreeReaderValue<float>*  rMCW    = hasMCW     ?
        new TTreeReaderValue<float>(reader,"mcWeight") : nullptr;
    TTreeReaderValue<float>*  rSFPh   = hasSFPh    ?
        new TTreeReaderValue<float>(reader,"scaleFactor_PHOTON") : nullptr;
    TTreeReaderValue<float>*  rSFTr   = hasSFTrig  ?
        new TTreeReaderValue<float>(reader,"scaleFactor_PhotonTRIGGER") : nullptr;

    Long64_t total = chain->GetEntries();
    Long64_t proc  = 0;

    while (reader.Next()) {
        proc++;
        if (proc % 1000000 == 0)
            printf("  %.0f / %.0f M\r",
                   proc/1e6, total/1e6);

        // Require at least 2 photons
        if (photon_pt.GetSize() < 2) continue;

        // Trigger requirement
        if (rTrig && !(**rTrig)) continue;

        // TightID
        if (rTight && rTight->GetSize() >= 2)
            if (!(*rTight)[0] || !(*rTight)[1]) continue;

        double pt1  = photon_pt[0]/1000.;
        double pt2  = photon_pt[1]/1000.;
        double eta1 = photon_eta[0];
        double eta2 = photon_eta[1];
        double phi1 = photon_phi[0];
        double phi2 = photon_phi[1];

        // pT threshold
        if (pt1 < 25. || pt2 < 25.) continue;

        // Crack veto
        auto crack=[](double e){
            return (fabs(e)>2.37)||(fabs(e)>1.37&&fabs(e)<1.52);
        };
        if (crack(eta1)||crack(eta2)) continue;

        // Energy
        double E1 = (rE && rE->GetSize()>=2) ?
                    (*rE)[0]/1000. : pt1*TMath::CosH(eta1);
        double E2 = (rE && rE->GetSize()>=2) ?
                    (*rE)[1]/1000. : pt2*TMath::CosH(eta2);

        // Isolation
        float iso1=0,iso2=0,ptiso1=0,ptiso2=0;
        if (rEtc && rEtc->GetSize()>=2) {
            if ((*rEtc)[0]/1000.>4.||(*rEtc)[1]/1000.>4.) continue;
            iso1=(*rEtc)[0]/1000.; iso2=(*rEtc)[1]/1000.;
        }
        if (rPtc && rPtc->GetSize()>=2) {
            ptiso1=(*rPtc)[0]/1000.; ptiso2=(*rPtc)[1]/1000.;
        }

        // 4-vectors
        TLorentzVector p1,p2;
        p1.SetPtEtaPhiE(pt1,eta1,phi1,E1);
        p2.SetPtEtaPhiE(pt2,eta2,phi2,E2);
        TLorentzVector pgg=p1+p2;
        double myy=pgg.M();

        // Mass window + relative pT cuts
        if (myy<MYY_LO||myy>MYY_HI) continue;
        if (pt1/myy<0.35||pt2/myy<0.25) continue;

        // BDT variables
        float pt1n = pt1/myy, pt2n = pt2/myy;
        float ae1  = fabs(eta1), ae2=fabs(eta2);
        float ade  = fabs(eta1-eta2);
        float adp  = fabs(TVector2::Phi_mpi_pi(phi1-phi2));
        double denom = myy*TMath::Sqrt(myy*myy+pgg.Pt()*pgg.Pt());
        float costh = (denom>0) ?
            fabs((p1.Pz()*(E2-p2.Pz())-p2.Pz()*(E1-p1.Pz()))/denom) : 0.f;

        // Event weight
        float wt = 1.0f;
        if (hasMCW  && rMCW)  wt *= **rMCW;
        if (hasSFPh && rSFPh) wt *= **rSFPh;
        if (hasSFTrig&& rSFTr) wt *= **rSFTr;

        // Region flags
        bool inSR = (myy>SR_LO && myy<SR_HI);
        bool inSB = (myy>SB1_LO&&myy<SB1_HI)||
                    (myy>SB2_LO&&myy<SB2_HI);

        // Fill histogram
        if (hOut && !sbOnly) hOut->Fill(myy);

        // Fill output tree
        if (tOut) {
            bool doFill = (!sbOnly && !srOnly) ||
                          (sbOnly && inSB)     ||
                          (srOnly && inSR);
            if (doFill) {
                *t_pt1n=pt1n; *t_pt2n=pt2n;
                *t_ae1=ae1;   *t_ae2=ae2;
                *t_ade=ade;   *t_adp=adp;
                *t_cs=costh;
                *t_i1=iso1;   *t_i2=iso2;
                *t_pi1=ptiso1;*t_pi2=ptiso2;
                *t_myy=myy;   *t_w=wt;
                tOut->Fill();
            }
        }

        // Score and collect events
        if (scorer && evOut) {
            *r_pt1n=pt1n; *r_pt2n=pt2n;
            *r_ae1=ae1;   *r_ae2=ae2;
            *r_ade=ade;   *r_adp=adp;
            *r_cs=costh;
            *r_i1=iso1;   *r_i2=iso2;
            *r_pi1=ptiso1;*r_pi2=ptiso2;
            float sc = scorer->EvaluateMVA("BDT");
            evOut->push_back({myy, sc});
        }
    }
    printf("\n");

    // Cleanup optional readers
    delete rTrig; delete rE; delete rTight;
    delete rEtc;  delete rPtc;
    delete rMCW;  delete rSFPh; delete rSFTr;
}

// ============================================================
// MAIN ANALYSIS
// ============================================================
void higgs_diphoton_analysis(
    const char* dataDir     = "/hepdata/Atlas/GamGam/Data/",
    const char* mcDir       = "/hepdata/Atlas/GamGam/MC/",
    const char* fullDataDir = "/hepdata/Atlas/GamGam/FullData/")
{
    SetAtlasStyle();
    mkdir("plots",       0755);
    mkdir("tmva_output", 0755);

    // --------------------------------------------------------
    // STEP 1: Open chains
    // --------------------------------------------------------
    printf("\n========== STEP 1: OPENING DATA ==========\n");
    TChain* chData     = MakeChain(dataDir,     "Data");
    TChain* chMC       = MakeChain(mcDir,       "MC");
    TChain* chFullData = MakeChain(fullDataDir, "FullData");

    bool useFullData = (chFullData && chFullData->GetEntries()>0);
    TChain* chBest   = useFullData ? chFullData : chData;
    printf("[INFO] Best chain for scoring: %lld entries\n",
           chBest->GetEntries());

    if (!chData||chData->GetEntries()==0||
        !chMC  ||chMC->GetEntries()==0) {
        printf("[ERROR] Data or MC empty.\n"); return;
    }

    // --------------------------------------------------------
    // Output tree branches (shared across steps)
    // --------------------------------------------------------
    float t_pt1n,t_pt2n,t_ae1,t_ae2,t_ade,t_adp,
          t_cs,t_i1,t_i2,t_pi1,t_pi2,t_w;
    double t_myy;

    // TMVA reader variables
    float r_pt1n,r_pt2n,r_ae1,r_ae2,r_ade,r_adp,
          r_cs,r_i1,r_i2,r_pi1,r_pi2;

    // --------------------------------------------------------
    // STEP 2: BASELINE SPECTRUM
    // --------------------------------------------------------
    printf("\n========== STEP 2: BASELINE ==========\n");

    TH1F* hBase = new TH1F("hBase","Baseline",NBINS,MYY_LO,MYY_HI);
    hBase->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hBase->GetYaxis()->SetTitle("Events / 1 GeV");
    hBase->SetMarkerStyle(20); hBase->SetMarkerSize(0.8);

    ProcessChain(chData, false,
                 hBase,  nullptr, nullptr, nullptr,
                 false, false,
                 &t_pt1n,&t_pt2n,&t_ae1,&t_ae2,
                 &t_ade,&t_adp,&t_cs,&t_i1,&t_i2,
                 &t_pi1,&t_pi2,&t_w,&t_myy,
                 &r_pt1n,&r_pt2n,&r_ae1,&r_ae2,
                 &r_ade,&r_adp,&r_cs,&r_i1,&r_i2,
                 &r_pi1,&r_pi2);

    printf("[STEP 2] Baseline entries: %d\n",(int)hBase->GetEntries());

    // Global S+B fit
    TF1* fSpBBase = new TF1("fSpBBase",SigPlusBkg,MYY_LO,MYY_HI,5);
    fSpBBase->SetParameters(hBase->GetMaximum(),-0.03,10.,125.09,1.8);
    fSpBBase->SetParLimits(2,0.,hBase->GetMaximum()*5.);
    fSpBBase->SetParLimits(3,122.,128.);
    fSpBBase->SetParLimits(4,0.5,4.0);
    fSpBBase->SetLineColor(kRed); fSpBBase->SetLineWidth(2);
    hBase->Fit(fSpBBase,"RQN");

    TF1* fBkgBase = new TF1("fBkgBase",BkgOnly,MYY_LO,MYY_HI,2);
    fBkgBase->SetParameters(fSpBBase->GetParameter(0),
                             fSpBBase->GetParameter(1));
    fBkgBase->SetLineColor(kBlue+1);
    fBkgBase->SetLineStyle(2); fBkgBase->SetLineWidth(2);

    double S0=0,B0=0,sigB0=0;
    ExtractSB(fSpBBase,S0,B0,sigB0);
    printf("[STEP 2] S0=%.1f B0=%.1f S0/sqrt(B0)=%.2f\n",S0,B0,sigB0);

    // Plot
    TCanvas* c1=new TCanvas("c1","Baseline",800,600);
    hBase->Draw("E");
    fSpBBase->Draw("SAME"); fBkgBase->Draw("SAME");
    TLegend* lg1=new TLegend(0.55,0.65,0.93,0.90);
    lg1->SetBorderSize(0); lg1->SetFillStyle(0);
    lg1->AddEntry(hBase,   "Data",    "lpe");
    lg1->AddEntry(fSpBBase,"Sig+Bkg", "l");
    lg1->AddEntry(fBkgBase,"Bkg only","l");
    lg1->Draw();
    TLatex tx1; tx1.SetNDC(); tx1.SetTextSize(0.034);
    tx1.DrawLatex(0.18,0.87,
        "ATLAS Open Data  #sqrt{s} = 13 TeV");
    tx1.DrawLatex(0.18,0.82,
        Form("S_{0}=%.0f  B_{0}=%.0f  S_{0}/#sqrt{B_{0}}=%.2f",
             S0,B0,sigB0));
    c1->SaveAs("plots/Fig1_baseline_spectrum.pdf");
    c1->SaveAs("plots/Fig1_baseline_spectrum.png");
    printf("[STEP 2] Fig1 saved.\n");

    // Baseline p-value
    printf("[STEP 2] Computing baseline p-value...\n");
    TGraph* grPvalBase=SlidingWindowPvalue(hBase,"grPvalBase");
    grPvalBase->SetLineColor(kBlue+1);
    grPvalBase->SetMarkerColor(kBlue+1);
    DrawPvalGraph(grPvalBase,
        "Baseline local p-value (Pre-BDT)",
        "Fig7_pvalue_baseline", kBlue+1);

    // --------------------------------------------------------
    // STEP 3: GENERATE ML DATASETS
    // --------------------------------------------------------
    printf("\n========== STEP 3: ML DATASETS ==========\n");

    // Signal tree from MC
    TFile* fSigOut=new TFile("tmva_output/signal.root","RECREATE");
    TTree* tSig=new TTree("sig","signal");
    auto BookTree=[&](TTree* t){
        t->Branch("pt1_norm",&t_pt1n,"pt1_norm/F");
        t->Branch("pt2_norm",&t_pt2n,"pt2_norm/F");
        t->Branch("aeta1",  &t_ae1, "aeta1/F");
        t->Branch("aeta2",  &t_ae2, "aeta2/F");
        t->Branch("adeta",  &t_ade, "adeta/F");
        t->Branch("adphi",  &t_adp, "adphi/F");
        t->Branch("costh",  &t_cs,  "costh/F");
        t->Branch("iso1",   &t_i1,  "iso1/F");
        t->Branch("iso2",   &t_i2,  "iso2/F");
        t->Branch("ptiso1", &t_pi1, "ptiso1/F");
        t->Branch("ptiso2", &t_pi2, "ptiso2/F");
        t->Branch("myy",    &t_myy, "myy/D");
        t->Branch("weight", &t_w,   "weight/F");
    };
    fSigOut->cd(); BookTree(tSig);

    ProcessChain(chMC, true,
                 nullptr, tSig, nullptr, nullptr,
                 false, false,  // all MC events (no region filter)
                 &t_pt1n,&t_pt2n,&t_ae1,&t_ae2,
                 &t_ade,&t_adp,&t_cs,&t_i1,&t_i2,
                 &t_pi1,&t_pi2,&t_w,&t_myy,
                 &r_pt1n,&r_pt2n,&r_ae1,&r_ae2,
                 &r_ade,&r_adp,&r_cs,&r_i1,&r_i2,
                 &r_pi1,&r_pi2);

    Long64_t nSig=tSig->GetEntries();
    tSig->Write(); fSigOut->Close(); gROOT->cd();
    printf("[STEP 3] Signal events: %lld\n",nSig);

    // Background tree from sidebands
    TFile* fBkgOut=new TFile("tmva_output/background.root","RECREATE");
    TTree* tBkg=new TTree("bkg","background");
    fBkgOut->cd(); BookTree(tBkg);

    // Use FullData sidebands if available, else small Data
    ProcessChain(chBest, false,
                 nullptr, tBkg, nullptr, nullptr,
                 true, false,  // sbOnly=true
                 &t_pt1n,&t_pt2n,&t_ae1,&t_ae2,
                 &t_ade,&t_adp,&t_cs,&t_i1,&t_i2,
                 &t_pi1,&t_pi2,&t_w,&t_myy,
                 &r_pt1n,&r_pt2n,&r_ae1,&r_ae2,
                 &r_ade,&r_adp,&r_cs,&r_i1,&r_i2,
                 &r_pi1,&r_pi2);

    Long64_t nBkg=tBkg->GetEntries();
    tBkg->Write(); fBkgOut->Close(); gROOT->cd();
    printf("[STEP 3] Background events: %lld\n",nBkg);

    // --------------------------------------------------------
    // STEP 4: TRAIN BDT
    // --------------------------------------------------------
    printf("\n========== STEP 4: BDT TRAINING ==========\n");

    TMVA::Tools::Instance();
    TFile* fTMVA=new TFile("tmva_output/TMVA.root","RECREATE");
    TMVA::Factory* factory=new TMVA::Factory(
        "TMVAClassif",fTMVA,
        "!V:!Silent:Color:DrawProgressBar:"
        "Transformations=I:AnalysisType=Classification");
    TMVA::DataLoader* loader=
        new TMVA::DataLoader("tmva_output/dataset");

    TFile* fS2=TFile::Open("tmva_output/signal.root");
    TFile* fB2=TFile::Open("tmva_output/background.root");
    loader->AddSignalTree    ((TTree*)fS2->Get("sig"),1.0);
    loader->AddBackgroundTree((TTree*)fB2->Get("bkg"),1.0);
    loader->SetSignalWeightExpression("weight");

    loader->AddVariable("pt1_norm","","",'F');
    loader->AddVariable("pt2_norm","","",'F');
    loader->AddVariable("aeta1",   "","'F");
    loader->AddVariable("aeta2",   "","'F");
    loader->AddVariable("adeta",   "","'F");
    loader->AddVariable("adphi",   "","'F");
    loader->AddVariable("costh",   "","'F");
    loader->AddVariable("iso1",    "","'F");
    loader->AddVariable("iso2",    "","'F");
    loader->AddVariable("ptiso1",  "","'F");
    loader->AddVariable("ptiso2",  "","'F");

    Long64_t nTrS=std::min((Long64_t)60000,nSig*7/10);
    Long64_t nTeS=std::min((Long64_t)20000,nSig*3/10);
    Long64_t nTrB=std::min((Long64_t)60000,nBkg*7/10);
    Long64_t nTeB=std::min((Long64_t)20000,nBkg*3/10);

    loader->PrepareTrainingAndTestTree("",
        Form("nTrain_Signal=%lld:nTest_Signal=%lld:"
             "nTrain_Background=%lld:nTest_Background=%lld:"
             "SplitMode=Random:SplitSeed=42:"
             "NormMode=NumEvents:!V",
             nTrS,nTeS,nTrB,nTeB));

    factory->BookMethod(loader,TMVA::Types::kBDT,"BDT",
        "!H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=4:"
        "BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:"
        "BaggedSampleFraction=0.8:nCuts=30:DoBoostMonitor");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    fTMVA->Close();
    printf("[STEP 4] BDT training complete.\n");

    // --------------------------------------------------------
    // STEP 5: APPLY BDT + OPTIMIZE
    // --------------------------------------------------------
    printf("\n========== STEP 5: BDT SCORING ==========\n");

    TMVA::Reader* tmvaReader=new TMVA::Reader("!Color:!Silent");
    tmvaReader->AddVariable("pt1_norm",&r_pt1n);
    tmvaReader->AddVariable("pt2_norm",&r_pt2n);
    tmvaReader->AddVariable("aeta1",   &r_ae1);
    tmvaReader->AddVariable("aeta2",   &r_ae2);
    tmvaReader->AddVariable("adeta",   &r_ade);
    tmvaReader->AddVariable("adphi",   &r_adp);
    tmvaReader->AddVariable("costh",   &r_cs);
    tmvaReader->AddVariable("iso1",    &r_i1);
    tmvaReader->AddVariable("iso2",    &r_i2);
    tmvaReader->AddVariable("ptiso1",  &r_pi1);
    tmvaReader->AddVariable("ptiso2",  &r_pi2);
    tmvaReader->BookMVA("BDT",
        "tmva_output/dataset/weights/TMVAClassif_BDT.weights.xml");

    std::vector<EvScore> evts;
    evts.reserve(5000000);

    ProcessChain(chBest, false,
                 nullptr, nullptr, &evts, tmvaReader,
                 false, false,
                 &t_pt1n,&t_pt2n,&t_ae1,&t_ae2,
                 &t_ade,&t_adp,&t_cs,&t_i1,&t_i2,
                 &t_pi1,&t_pi2,&t_w,&t_myy,
                 &r_pt1n,&r_pt2n,&r_ae1,&r_ae2,
                 &r_ade,&r_adp,&r_cs,&r_i1,&r_i2,
                 &r_pi1,&r_pi2);

    printf("[STEP 5] Events scored: %zu\n",evts.size());
    if (evts.empty()) { printf("[ERROR] No events.\n"); return; }

    // BDT threshold scan
    float sMin=1.,sMax=-1.;
    for (auto& e:evts) {
        sMin=std::min(sMin,e.score);
        sMax=std::max(sMax,e.score);
    }
    const int NSCAN=40;
    double scanT[NSCAN],scanZ[NSCAN];
    double bestT=0.,bestZ=0.,bestS=0.,bestB=0.;
    double tLo=sMin+0.05*(sMax-sMin);
    double tHi=sMin+0.95*(sMax-sMin);

    printf("[STEP 5] Scanning thresholds...\n");
    for (int k=0;k<NSCAN;k++) {
        double tcut=tLo+(tHi-tLo)*k/(NSCAN-1);
        scanT[k]=tcut;
        TH1F hTmp(Form("hSc_%d",k),"",NBINS,MYY_LO,MYY_HI);
        for (auto& e:evts) if(e.score>tcut) hTmp.Fill(e.myy);
        if (hTmp.GetEntries()<50) {scanZ[k]=0.;continue;}
        TF1 fTmp("fTmpSB",SigPlusBkg,MYY_LO,MYY_HI,5);
        fTmp.SetParameters(hTmp.GetMaximum()*0.9,-0.04,
                            hTmp.GetMaximum()*0.05,125.09,2.0);
        fTmp.SetParLimits(2,0.,hTmp.GetMaximum()*5.);
        fTmp.SetParLimits(3,115.,135.);
        fTmp.SetParLimits(4,0.5,6.0);
        hTmp.Fit(&fTmp,"RQBN");
        double S_k=0,B_k=0,z_k=0;
        ExtractSB(&fTmp,S_k,B_k,z_k);
        scanZ[k]=z_k;
        if(z_k>bestZ){bestZ=z_k;bestT=tcut;bestS=S_k;bestB=B_k;}
        printf("  t=%.3f S=%.1f B=%.1f Z=%.2f\n",tcut,S_k,B_k,z_k);
    }
    printf("[STEP 5] Optimal: t=%.4f S=%.1f B=%.1f Z=%.2f\n",
           bestT,bestS,bestB,bestZ);

    // Optimization plot
    TCanvas* c2=new TCanvas("c2","Opt",800,600);
    TGraph* grOpt=new TGraph(NSCAN,scanT,scanZ);
    grOpt->SetTitle(";BDT Cut;S/#sqrt{B}");
    grOpt->SetMarkerStyle(20); grOpt->SetMarkerColor(kBlue+1);
    grOpt->SetLineColor(kBlue+1); grOpt->SetLineWidth(2);
    grOpt->Draw("APL");
    TLine* lB=new TLine(bestT,0,bestT,bestZ*1.15);
    lB->SetLineColor(kRed); lB->SetLineStyle(2);
    lB->SetLineWidth(2); lB->Draw();
    TLatex tx2; tx2.SetNDC(); tx2.SetTextSize(0.034);
    tx2.DrawLatex(0.18,0.88,"BDT Cut Optimisation");
    tx2.SetTextColor(kRed);
    tx2.DrawLatex(0.18,0.83,
        Form("t_{BDT}=%.3f  S/#sqrt{B}=%.2f",bestT,bestZ));
    c2->SaveAs("plots/Fig2_bdt_optimisation.pdf");
    c2->SaveAs("plots/Fig2_bdt_optimisation.png");

    // --------------------------------------------------------
    // STEP 6: BDT SPECTRUM + FIT
    // --------------------------------------------------------
    printf("\n========== STEP 6: BDT SPECTRUM ==========\n");

    TH1F* hBDT=new TH1F("hBDT","BDT",NBINS,MYY_LO,MYY_HI);
    hBDT->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hBDT->GetYaxis()->SetTitle("Events / 1 GeV");
    hBDT->SetMarkerStyle(20); hBDT->SetMarkerSize(0.8);

    TH1F* hSAll=new TH1F("hSAll","",100,-1,1);
    TH1F* hSAcc=new TH1F("hSAcc","",100,-1,1);
    TH1F* hSRej=new TH1F("hSRej","",100,-1,1);
    for (auto& e:evts) {
        hSAll->Fill(e.score);
        if(e.score>bestT){hBDT->Fill(e.myy);hSAcc->Fill(e.score);}
        else hSRej->Fill(e.score);
    }
    printf("[STEP 6] Events after BDT: %d\n",(int)hBDT->GetEntries());

    TF1* fSpBBDT=new TF1("fSpBBDT",SigPlusBkg,MYY_LO,MYY_HI,5);
    fSpBBDT->SetParameters(hBDT->GetMaximum()*0.9,-0.04,
                            hBDT->GetMaximum()*0.05,125.09,2.0);
    fSpBBDT->SetParLimits(2,0.,hBDT->GetMaximum()*5.);
    fSpBBDT->SetParLimits(3,115.,135.);
    fSpBBDT->SetParLimits(4,0.5,6.0);
    fSpBBDT->SetLineColor(kRed); fSpBBDT->SetLineWidth(2);
    hBDT->Fit(fSpBBDT,"RQ");

    TF1* fBkgBDT=new TF1("fBkgBDT",BkgOnly,MYY_LO,MYY_HI,2);
    fBkgBDT->SetParameters(fSpBBDT->GetParameter(0),
                            fSpBBDT->GetParameter(1));
    fBkgBDT->SetLineColor(kBlue+1);
    fBkgBDT->SetLineStyle(2); fBkgBDT->SetLineWidth(2);

    double S_bdt=0,B_bdt=0,Z_bdt=0;
    ExtractSB(fSpBBDT,S_bdt,B_bdt,Z_bdt);
    double mu_fit=fSpBBDT->GetParameter(3);
    double sg_fit=fSpBBDT->GetParameter(4);
    printf("[STEP 6] S=%.1f B=%.1f Z=%.2f mu=%.2f sg=%.2f\n",
           S_bdt,B_bdt,Z_bdt,mu_fit,sg_fit);

    // BDT spectrum plot
    TCanvas* c4=new TCanvas("c4","BDT",800,600);
    hBDT->Draw("E");
    fSpBBDT->Draw("SAME"); fBkgBDT->Draw("SAME");
    TBox* srb=new TBox(SR_LO,0,SR_HI,hBDT->GetMaximum()*0.55);
    srb->SetFillColorAlpha(kGreen+1,0.18); srb->Draw();
    hBDT->Draw("E SAME");
    TLegend* lg4=new TLegend(0.55,0.62,0.93,0.88);
    lg4->SetBorderSize(0); lg4->SetFillStyle(0);
    lg4->AddEntry(hBDT,   "Data",    "lpe");
    lg4->AddEntry(fSpBBDT,"Sig+Bkg", "l");
    lg4->AddEntry(fBkgBDT,"Bkg only","l");
    lg4->Draw();
    TLatex tx4; tx4.SetNDC(); tx4.SetTextSize(0.032);
    tx4.DrawLatex(0.18,0.88,
        "ATLAS Open Data  #sqrt{s} = 13 TeV");
    tx4.DrawLatex(0.18,0.83,
        Form("BDT > %.3f   N=%d",bestT,(int)hBDT->GetEntries()));
    tx4.SetTextColor(kRed+1);
    tx4.DrawLatex(0.18,0.78,
        Form("#mu=%.2f GeV  #sigma=%.2f GeV",mu_fit,sg_fit));
    tx4.SetTextColor(kBlack);
    tx4.DrawLatex(0.18,0.73,
        Form("S=%.0f  B=%.0f  S/#sqrt{B}=%.2f",
             S_bdt,B_bdt,Z_bdt));
    c4->SaveAs("plots/Fig4_bdt_spectrum.pdf");
    c4->SaveAs("plots/Fig4_bdt_spectrum.png");

    // Residual plot
    TH1F* hResid=(TH1F*)hBDT->Clone("hResid");
    hResid->Reset();
    hResid->GetYaxis()->SetTitle("Data #minus Background");
    for (int ib=1;ib<=hBDT->GetNbinsX();ib++) {
        double x=hBDT->GetBinCenter(ib);
        hResid->SetBinContent(ib,
            hBDT->GetBinContent(ib)-fBkgBDT->Eval(x));
        hResid->SetBinError(ib,hBDT->GetBinError(ib));
    }
    TCanvas* c6=new TCanvas("c6","Resid",800,600);
    c6->SetGridy();
    hResid->SetMarkerStyle(20); hResid->SetMarkerSize(0.8);
    hResid->Draw("E"); c6->Update();
    TBox* srb6=new TBox(SR_LO,c6->GetUymin(),
                         SR_HI,c6->GetUymax());
    srb6->SetFillColorAlpha(kGreen+1,0.18); srb6->Draw();
    TLine* lz=new TLine(MYY_LO,0,MYY_HI,0);
    lz->SetLineColor(kRed); lz->SetLineStyle(2);
    lz->SetLineWidth(2); lz->Draw();
    hResid->Draw("E SAME");
    TLatex tx6; tx6.SetNDC(); tx6.SetTextSize(0.034);
    tx6.DrawLatex(0.18,0.88,"Background-Subtracted Residual");
    c6->SaveAs("plots/Fig6_residuals.pdf");
    c6->SaveAs("plots/Fig6_residuals.png");

    // --------------------------------------------------------
    // STEP 7: BDT P-VALUE
    // --------------------------------------------------------
    printf("\n========== STEP 7: BDT P-VALUE ==========\n");
    TGraph* grPvalBDT=SlidingWindowPvalue(hBDT,"grPvalBDT");
    DrawPvalGraph(grPvalBDT,
        Form("BDT local p-value (t_{BDT} > %.3f)",bestT),
        "Fig8_pvalue_bdt", kRed+1);

    // --------------------------------------------------------
    // STEP 8: PERFORMANCE TABLE
    // --------------------------------------------------------
    printf("\n========== STEP 8: PERFORMANCE ==========\n");

    double eps_base=1.0;
    double P_base  =(S0+B0>0)?S0/(S0+B0):0.;
    double eP_base = P_base;

    double eps_bdt=(S0>0)?S_bdt/S0:0.;
    double P_bdt  =(S_bdt+B_bdt>0)?S_bdt/(S_bdt+B_bdt):0.;
    double eP_bdt =(S0>0&&S_bdt+B_bdt>0)?
                   (S_bdt*S_bdt)/(S0*(S_bdt+B_bdt)):0.;

    printf("\n╔══════════════════════════════════════════════════════════╗\n");
    printf(  "║                   PERFORMANCE TABLE                      ║\n");
    printf(  "╠══════════════════════════════════════════════════════════╣\n");
    printf(  "║ Method    S       B        S/√B   ε      P       εP     ║\n");
    printf(  "╠══════════════════════════════════════════════════════════╣\n");
    printf(  "║ Baseline %6.1f %8.1f %5.2f %5.3f %6.4f %6.4f  ║\n",
             S0,B0,sigB0,eps_base,P_base,eP_base);
    printf(  "║ BDT      %6.1f %8.1f %5.2f %5.3f %6.4f %6.4f  ║\n",
             S_bdt,B_bdt,Z_bdt,eps_bdt,P_bdt,eP_bdt);
    printf(  "╠══════════════════════════════════════════════════════════╣\n");
    printf(  "║ Purity improvement: %.1fx                                ║\n",
             (P_base>0?P_bdt/P_base:0.));
    printf(  "║ εP improvement:     %.1fx                                ║\n",
             (eP_base>0?eP_bdt/eP_base:0.));
    printf(  "╚══════════════════════════════════════════════════════════╝\n\n");

    printf("Fitted Higgs mass: %.2f ± %.2f GeV\n",mu_fit,sg_fit);
    printf("BDT significance:  %.2f sigma\n\n", Z_bdt);

    // Cleanup
    delete tmvaReader;
    fS2->Close(); fB2->Close();
    printf("[DONE] All plots saved to plots/\n");
}
