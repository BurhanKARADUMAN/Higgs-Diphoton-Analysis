# -------------------------------------------------------------------------------------
# skeleton_atlas_H2llll.py
#
# Example PyROOT to read/analyse Atlas data files
# for analyzing H -> 4 leptons events
# -------------------------------------------------------------------------------------

import ROOT
import math

# --- Enable multi-threading ---
ROOT.EnableImplicitMT()

MonteCarlo = False

# --- Chain setup ---
chain = ROOT.TChain("mini")
if MonteCarlo:
   print("Running ATLAS Monte Carlo Simulation...")
   fileout = ROOT.TFile("m4lepMC.root", "recreate")
   chain.Add("/hepdata/Atlas/4lep/MC/mc_363490.llll.4lep.root")
else:
   print("Running ATLAS Recorded Data...")
   fileout = ROOT.TFile("m4lepData.root", "recreate")
   chain.Add("/hepdata/Atlas/4lep/Data/*.root")

# --- Histogram ---
h_m4lep = ROOT.TH1F("h_m4lep",";m_{llll} [GeV]; Entries",36, 80, 170)

# --- Event loop ---
fraction = 1
nentries = int(fraction * chain.GetEntries())
print(f"Nevents: {nentries}")

for ev in range(nentries):
    chain.GetEntry(ev)

    if ev % 1000 == 0:
        pct = math.ceil(100.0 * ev / nentries)
        print(f"{ev} / {nentries} [{pct}%]", end="\r", flush=True)

    # Read branches as Python lists
    pt       = list(chain.lep_pt)
    eta      = list(chain.lep_eta)
    phi      = list(chain.lep_phi)
    e        = list(chain.lep_E)

    # Build four-vectors
    p1 = ROOT.TLorentzVector()
    p2 = ROOT.TLorentzVector()
    p3 = ROOT.TLorentzVector()
    p4 = ROOT.TLorentzVector()
    p1.SetPtEtaPhiE(pt[0], eta[0], phi[0], e[0])
    p2.SetPtEtaPhiE(pt[1], eta[1], phi[1], e[1])
    p3.SetPtEtaPhiE(pt[2], eta[2], phi[2], e[2])
    p4.SetPtEtaPhiE(pt[3], eta[3], phi[3], e[3])

    fourlep = p1 + p2 + p3 + p4
    m = fourlep.M() * 1e-3 # convert of GeV
    h_m4lep.Fill(m)

# --- Save and plot ---
fileout.Write()
fileout.Close()
print("\nDone.")
