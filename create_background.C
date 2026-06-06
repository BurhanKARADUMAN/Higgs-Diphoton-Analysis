#include <TFile.h>

#include <TTree.h>

#include <TLorentzVector.h>

#include <iostream>

#include <vector>



void create_background() {


    TFile *fin = TFile::Open("/hepdata/Atlas/GamGam/FullData/*.root");

    TTree *tin = (TTree*)fin->Get("analysis");



    std::vector<float> *photon_pt = 0;

    std::vector<float> *photon_eta = 0;

    std::vector<float> *photon_phi = 0;

    std::vector<float> *photon_e = 0;

    Int_t photon_n = 0;



    tin->SetBranchAddress("photon_n", &photon_n);

    tin->SetBranchAddress("photon_pt", &photon_pt);

    tin->SetBranchAddress("photon_eta", &photon_eta);

    tin->SetBranchAddress("photon_phi", &photon_phi);

    tin->SetBranchAddress("photon_e", &photon_e);



    // 2. Yeni oluşturulacak background.root dosyasını ve ağacını tanımla

    TFile *fout = new TFile("background.root", "RECREATE");

    TTree *tout = new TTree("mlTree", "Background data for ML");



    // ML için düz değişkenler (Flat variables)

    float pt1, eta1, phi1, e1;

    float pt2, eta2, phi2, e2;

    float m_yy;



    // Yeni ağacın dallarını oluştur

    tout->Branch("pt1", &pt1, "pt1/F");

    tout->Branch("eta1", &eta1, "eta1/F");

    tout->Branch("phi1", &phi1, "phi1/F");

    tout->Branch("e1", &e1, "e1/F");

    tout->Branch("pt2", &pt2, "pt2/F");

    tout->Branch("eta2", &eta2, "eta2/F");

    tout->Branch("phi2", &phi2, "phi2/F");

    tout->Branch("e2", &e2, "e2/F");

    tout->Branch("m_yy", &m_yy, "m_yy/F");



    // 3. Veri üzerinde dön ve sadece yan bantları seç

    Long64_t nentries = tin->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {

        tin->GetEntry(i);



        if (photon_n >= 2) {

            TLorentzVector gamma1, gamma2, higgs;

            gamma1.SetPtEtaPhiE(photon_pt->at(0), photon_eta->at(0), photon_phi->at(0), photon_e->at(0));

            gamma2.SetPtEtaPhiE(photon_pt->at(1), photon_eta->at(1), photon_phi->at(1), photon_e->at(1));

            higgs = gamma1 + gamma2;

            m_yy = higgs.M();



            // YAN BANT (SIDEBAND) KOŞULU: (105-120) VEYA (130-160) arasındaysa

            if ((m_yy > 105 && m_yy < 120) || (m_yy > 130 && m_yy < 160)) {

                

                // Değerleri yeni düz değişkenlere ata

                pt1 = photon_pt->at(0); eta1 = photon_eta->at(0); phi1 = photon_phi->at(0); e1 = photon_e->at(0);

                pt2 = photon_pt->at(1); eta2 = photon_eta->at(1); phi2 = photon_phi->at(1); e2 = photon_e->at(1);

                

                // Ağaca ekle

                tout->Fill();

            }

        }

    }



    // 4. Dosyayı kaydet ve kapat

    fout->Write();

    fout->Close();

    fin->Close();

    

    std::cout << "background.root dosyasi basariyla olusturuldu!" << std::endl;

}
