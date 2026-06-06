import ROOT

file1 =ROOT.TFile("m4lepData.root")
file2 =ROOT.TFile("m4lepMC.root")

h1 = file1.Get("h_m4lep")
h2 = file2.Get("h_m4lep")

can = ROOT.TCanvas("c1","c1",800,600)
h1.Rebin(2)
h2.Rebin(2)

h1.SetMarkerSize(1)
h1.SetMarkerStyle(20)

h2.Scale(h1.Integral()/h2.Integral())

h2.Draw("histo")
h1.Draw("same,error")

input("press enter to quit")
