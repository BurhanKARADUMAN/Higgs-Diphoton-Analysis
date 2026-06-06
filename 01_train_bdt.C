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
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
const double MYY_LO=100., MYY_HI=160., SR_LO=120., SR_HI=130., SB1_LO=105., SB1_HI=120., SB2_LO=130., SB2_HI=160.;
const int NBINS=60;
void SetAtlasStyle() {
    TStyle* s = new TStyle("ATLAS","ATLAS style"); s->SetOptStat(0); s->SetOptTitle(0); s->SetOptFit(0);
    s->SetPadTopMargin(0.06); s->SetPadRightMargin(0.05); s->SetPadBottomMargin(0.16); s->SetPadLeftMargin(0.16);
    s->SetMarkerStyle(20); s->SetMarkerSize(0.9); s->SetHistLineWidth(2); gROOT->SetStyle("ATLAS"); gROOT->ForceStyle();
}
double SigPlusBkg(double* x, double* p) { return p[0]*TMath::Exp(p[1]*x[0]) + p[2]*TMath::Gaus(x[0],p[3],p[4],false); }
double BkgOnly(double* x, double* p) { return p[0]*TMath::Exp(p[1]*x[0]); }
void ExtractSB(TF1* fSpB, double binWidth, double& S, double& B, double& sigSqrtB) {
    TF1 fb("_fb", BkgOnly, MYY_LO, MYY_HI, 2); fb.SetParameters(fSpB->GetParameter(0), fSpB->GetParameter(1)); B = fb.Integral(SR_LO, SR_HI) / binWidth;
    TF1 fs("_fs","[0]*TMath::Gaus(x,[1],[2],false)", MYY_LO, MYY_HI); fs.SetParameters(fSpB->GetParameter(2), fSpB->GetParameter(3), fSpB->GetParameter(4));
    S = std::max(0.0, fs.Integral(SR_LO, SR_HI)/binWidth); sigSqrtB = (B>0) ? S/TMath::Sqrt(B) : 0;
}
struct BDTVars { float pt1n, pt2n, aeta1, aeta2, adeta, adphi, costh, iso1, iso2; };
struct EventData {
    Bool_t trigP = false;
    std::vector<float> *pt=new std::vector<float>(), *eta=new std::vector<float>(), *phi=new std::vector<float>(), *E=new std::vector<float>(), *etcone20=new std::vector<float>(), *ptcone30=new std::vector<float>();
    std::vector<bool>* isTightID = new std::vector<bool>(); float mcWeight=1.0, scaleFactor_PHOTON=1.0, scaleFactor_PhotonTRIGGER=1.0;
    void Connect(TChain* ch, bool isMC) {
        ch->SetBranchStatus("*",0); auto bindVar = [&](const char* n, auto addr){ if(ch->GetBranch(n)){ ch->SetBranchStatus(n,1); ch->SetBranchAddress(n,addr); } };
        bindVar("trigP",&trigP); bindVar("photon_pt",&pt); bindVar("photon_eta",&eta); bindVar("photon_phi",&phi); bindVar("photon_E",&E);
        bindVar("photon_isTightID",&isTightID); bindVar("photon_etcone20",&etcone20); bindVar("photon_ptcone30",&ptcone30);
        if(isMC){ bindVar("mcWeight",&mcWeight); bindVar("scaleFactor_PHOTON",&scaleFactor_PHOTON); bindVar("scaleFactor_PhotonTRIGGER",&scaleFactor_PhotonTRIGGER); }
    }
};
void train_bdt(const char* dataDir="/hepdata/Atlas/GamGam/Data/", const char* mcDir="/hepdata/Atlas/GamGam/MC/") {
    SetAtlasStyle(); mkdir("plots",0755); mkdir("tmva_output",0755);
    auto MakeChain = [](const char* dir, const char* name) {
        std::vector<std::string> treeNames = {"mini", "analysis", "GamGam", "CollectionTree"};
        for(const auto& tname : treeNames){ TChain* ch=new TChain(tname.c_str()); ch->Add(Form("%s*.root",dir)); if(ch->GetEntries()>0){ printf("[INFO] %s: %lld\n",name,ch->GetEntries()); return ch; } delete ch; }
        return new TChain("mini"); 
    };
    TChain* chData = MakeChain(dataDir,"Data"); TChain* chMC = MakeChain(mcDir,"MC");
    if(chData->GetEntries()==0 || chMC->GetEntries()==0){ printf("[ERROR] Empty tree!\n"); return; }
    EventData dData; dData.Connect(chData,false); EventData dMC; dMC.Connect(chMC,true);
    auto Preselect = [&](BDTVars& v, EventData& evt) -> double {
        v={}; if(!evt.trigP || !evt.pt || evt.pt->size()<2) return -1.;
        if(evt.isTightID && evt.isTightID->size()>=2) if(!(*evt.isTightID)[0] || !(*evt.isTightID)[1]) return -1.;
        double pt1=(*evt.pt)[0]/1000., pt2=(*evt.pt)[1]/1000., eta1=(*evt.eta)[0], eta2=(*evt.eta)[1], phi1=(*evt.phi)[0], phi2=(*evt.phi)[1];
        double E1=(evt.E && evt.E->size()>=2)?(*evt.E)[0]/1000.:pt1*TMath::CosH(eta1), E2=(evt.E && evt.E->size()>=2)?(*evt.E)[1]/1000.:pt2*TMath::CosH(eta2);
        if(pt1<25. || pt2<25. || fabs(eta1)>2.37 || (fabs(eta1)>1.37&&fabs(eta1)<1.52) || fabs(eta2)>2.37 || (fabs(eta2)>1.37&&fabs(eta2)<1.52)) return -1.;
        if(evt.etcone20 && evt.etcone20->size()>=2){ if((*evt.etcone20)[0]/1000.>4. || (*evt.etcone20)[1]/1000.>4.) return -1.; v.iso1=(*evt.etcone20)[0]/1000.; v.iso2=(*evt.etcone20)[1]/1000.; } else { v.iso1=0.f; v.iso2=0.f; }
        TLorentzVector p1v, p2v; p1v.SetPtEtaPhiE(pt1,eta1,phi1,E1); p2v.SetPtEtaPhiE(pt2,eta2,phi2,E2); TLorentzVector pgg=p1v+p2v; double myy=pgg.M();
        if(myy<MYY_LO || myy>MYY_HI || pt1/myy<0.35 || pt2/myy<0.25) return -1.;
        v.pt1n=pt1/myy; v.pt2n=pt2/myy; v.aeta1=fabs(eta1); v.aeta2=fabs(eta2); v.adeta=fabs(eta1-eta2); v.adphi=fabs(TVector2::Phi_mpi_pi(phi1-phi2));
        double denom=myy*TMath::Sqrt(myy*myy+pgg.Pt()*pgg.Pt()); v.costh=(denom>0)?fabs((p1v.Pz()*(E2-p2v.Pz())-p2v.Pz()*(E1-p1v.Pz()))/denom):0.f; return myy;
    };
    TFile* fSigOut=new TFile("tmva_output/signal.root","RECREATE"); TTree* tSig=new TTree("sig","signal");
    TFile* fBkgOut=new TFile("tmva_output/background.root","RECREATE"); TTree* tBkg=new TTree("bkg","background");
    float t_pt1n, t_pt2n, t_ae1, t_ae2, t_ade, t_adp, t_cs, t_i1, t_i2, t_w; double t_myy;
    auto BookTree = [&](TTree* t){ t->Branch("pt1_norm",&t_pt1n,"pt1_norm/F"); t->Branch("pt2_norm",&t_pt2n,"pt2_norm/F"); t->Branch("aeta1",&t_ae1,"aeta1/F"); t->Branch("aeta2",&t_ae2,"aeta2/F"); t->Branch("adeta",&t_ade,"adeta/F"); t->Branch("adphi",&t_adp,"adphi/F"); t->Branch("costh",&t_cs,"costh/F"); t->Branch("iso1",&t_i1,"iso1/F"); t->Branch("iso2",&t_i2,"iso2/F"); t->Branch("myy",&t_myy,"myy/D"); t->Branch("weight",&t_w,"weight/F"); };
    fSigOut->cd(); BookTree(tSig); fBkgOut->cd(); BookTree(tBkg); BDTVars vv;
    for(Long64_t i=0; i<chMC->GetEntries(); i++){ chMC->GetEntry(i); double myy=Preselect(vv,dMC); if(myy<0) continue; t_pt1n=vv.pt1n; t_pt2n=vv.pt2n; t_ae1=vv.aeta1; t_ae2=vv.aeta2; t_ade=vv.adeta; t_adp=vv.adphi; t_cs=vv.costh; t_i1=vv.iso1; t_i2=vv.iso2; t_myy=myy; t_w=dMC.mcWeight*dMC.scaleFactor_PHOTON*dMC.scaleFactor_PhotonTRIGGER; tSig->Fill(); }
    tSig->Write(); fSigOut->Close(); gROOT->cd();
    for(Long64_t i=0; i<chData->GetEntries(); i++){ chData->GetEntry(i); double myy=Preselect(vv,dData); if(myy<0) continue; if((myy>SB1_LO && myy<SB1_HI)||(myy>SB2_LO && myy<SB2_HI)){ t_pt1n=vv.pt1n; t_pt2n=vv.pt2n; t_ae1=vv.aeta1; t_ae2=vv.aeta2; t_ade=vv.adeta; t_adp=vv.adphi; t_cs=vv.costh; t_i1=vv.iso1; t_i2=vv.iso2; t_myy=myy; t_w=1.f; tBkg->Fill(); } }
    tBkg->Write(); fBkgOut->Close(); gROOT->cd();
    TMVA::Tools::Instance(); TFile* fTMVA=new TFile("tmva_output/TMVA.root","RECREATE"); TMVA::Factory* factory=new TMVA::Factory("TMVAClassif",fTMVA,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification"); TMVA::DataLoader* loader=new TMVA::DataLoader("tmva_output/dataset");
    TFile* fS2=TFile::Open("tmva_output/signal.root"); TFile* fB2=TFile::Open("tmva_output/background.root"); loader->AddSignalTree((TTree*)fS2->Get("sig"),1.0); loader->AddBackgroundTree((TTree*)fB2->Get("bkg"),1.0); loader->SetSignalWeightExpression("weight");
    loader->AddVariable("pt1_norm","","",'F'); loader->AddVariable("pt2_norm","","",'F'); loader->AddVariable("aeta1","","",'F'); loader->AddVariable("aeta2","","",'F'); loader->AddVariable("adeta","","",'F'); loader->AddVariable("adphi","","",'F'); loader->AddVariable("costh","","",'F'); loader->AddVariable("iso1","","",'F'); loader->AddVariable("iso2","","",'F');
    loader->PrepareTrainingAndTestTree("","nTrain_Signal=50000:nTest_Signal=50000:nTrain_Background=20000:nTest_Background=20000:SplitMode=Random:SplitSeed=42:NormMode=NumEvents:!V");
    factory->BookMethod(loader,TMVA::Types::kBDT,"BDT","!H:!V:NTrees=200:MinNodeSize=5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.8:nCuts=20:DoBoostMonitor");
    factory->TrainAllMethods(); factory->TestAllMethods(); factory->EvaluateAllMethods(); fTMVA->Close(); printf("[INFO] TMVA Trained.\n");
    TMVA::Reader* reader=new TMVA::Reader("!Color:!Silent"); float r_pt1n,r_pt2n,r_ae1,r_ae2,r_ade,r_adp,r_cs,r_i1,r_i2;
    reader->AddVariable("pt1_norm",&r_pt1n); reader->AddVariable("pt2_norm",&r_pt2n); reader->AddVariable("aeta1",&r_ae1); reader->AddVariable("aeta2",&r_ae2); reader->AddVariable("adeta",&r_ade); reader->AddVariable("adphi",&r_adp); reader->AddVariable("costh",&r_cs); reader->AddVariable("iso1",&r_i1); reader->AddVariable("iso2",&r_i2);
    reader->BookMVA("BDT","tmva_output/dataset/weights/TMVAClassif_BDT.weights.xml");
    struct Ev{ double myy; float score; }; std::vector<Ev> evts; evts.reserve(300000);
    for(Long64_t i=0; i<chData->GetEntries(); i++){ chData->GetEntry(i); double myy=Preselect(vv,dData); if(myy<0) continue; r_pt1n=vv.pt1n; r_pt2n=vv.pt2n; r_ae1=vv.aeta1; r_ae2=vv.aeta2; r_ade=vv.adeta; r_adp=vv.adphi; r_cs=vv.costh; r_i1=vv.iso1; r_i2=vv.iso2; evts.push_back({myy,(float)reader->EvaluateMVA("BDT")}); }
    float sMin=1., sMax=-1.; for(auto& e:evts){ sMin=std::min(sMin,e.score); sMax=std::max(sMax,e.score); }
    double bestT=0., bestZ=0., bestS=0., bestB=0., tLo=sMin+0.10*(sMax-sMin), tHi=sMin+0.90*(sMax-sMin); const int NSCAN=25; double binW=(MYY_HI-MYY_LO)/NBINS;
    for(int k=0; k<NSCAN; k++){ double tcut=tLo+(tHi-tLo)*k/(NSCAN-1); TH1F hTmp(Form("hT_%d",k),"",NBINS,MYY_LO,MYY_HI); for(auto& e:evts) if(e.score>tcut) hTmp.Fill(e.myy);
        if(hTmp.GetEntries()<20) continue; TF1 fTmp("fTmp",SigPlusBkg,MYY_LO,MYY_HI,5); fTmp.SetParameters(hTmp.GetMaximum()*0.9,-0.04,hTmp.GetMaximum()*0.05,125.09,2.0); fTmp.SetParLimits(2,0.,hTmp.GetMaximum()*2.); fTmp.SetParLimits(3,115.,135.); fTmp.SetParLimits(4,0.5,6.0); hTmp.Fit(&fTmp,"RQBN");
        double S_k=0., B_k=0., z_k=0.; ExtractSB(&fTmp,binW,S_k,B_k,z_k); if(z_k>bestZ){ bestZ=z_k; bestT=tcut; } }
    std::ofstream out("optimal_cut.txt"); out << bestT; out.close(); printf("[SUCCESS] Optimal BDT Cut: %.4f saved to optimal_cut.txt\n", bestT);
}
