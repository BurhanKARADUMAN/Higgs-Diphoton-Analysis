#include <TFile.h>

#include <TTree.h>

#include <TChain.h>

#include <TLorentzVector.h>

#include <iostream>

#include <vector>



void create_signal() {


    TChain *tin = new TChain("mini");
    tin->Add("mc_*H125*.root");

    std::vector<float> *photon_pt = 0;

    std::vector<float> *photon_eta = 0;

    std::vector<float> *photon_phi = 0;


    UInt_t photon_n = 0;


    tin->SetBranchAddress("photon_n", &photon_n);

    tin->SetBranchAddress("photon_pt", &photon_pt);

    tin->SetBranchAddress("photon_eta", &photon_eta);

    tin->SetBranchAddress("photon_phi", &photon_phi);





    TFile *fout = new TFile("signal.root", "RECREATE");

    TTree *tout = new TTree("mlTree", "Signal data for ML");



    float pt1, eta1, phi1, e1;

    float pt2, eta2, phi2, e2;

    float m_yy;



    tout->Branch("pt1", &pt1, "pt1/F");

    tout->Branch("eta1", &eta1, "eta1/F");

    tout->Branch("phi1", &phi1, "phi1/F");

    tout->Branch("e1", &e1, "e1/F");

    tout->Branch("pt2", &pt2, "pt2/F");

    tout->Branch("eta2", &eta2, "eta2/F");

    tout->Branch("phi2", &phi2, "phi2/F");

    tout->Branch("e2", &e2, "e2/F");

    tout->Branch("m_yy", &m_yy, "m_yy/F");




    Long64_t nentries = tin->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {

        tin->GetEntry(i);




        if (photon_n >= 2) {

            TLorentzVector gamma1, gamma2, higgs;

            gamma1.SetPtEtaPhiM(photon_pt->at(0), photon_eta->at(0), photon_phi->at(0), 0.0);
            gamma2.SetPtEtaPhiM(photon_pt->at(1), photon_eta->at(1), photon_phi->at(1), 0.0);

            higgs = gamma1 + gamma2;

            m_yy = higgs.M();




           pt1 = photon_pt->at(0); eta1 = photon_eta->at(0); phi1 = photon_phi->at(0); e1 = gamma1.E();
           pt2 = photon_pt->at(1); eta2 = photon_eta->at(1); phi2 = photon_phi->at(1); e2 = gamma2.E();
            

            tout->Fill();

        }

    }




    fout->Write();

    fout->Close();

    

    std::cout << "signal.root dosyasi basariyla olusturuldu!" << std::endl;

}
