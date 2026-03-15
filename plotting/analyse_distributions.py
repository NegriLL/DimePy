#!/usr/bin/env python3

import ROOT
from pathlib import Path
from numpy import arctanh, linalg

# Calculate pseudorapidity
def pseudorapidity(px, py, pz):
    P = linalg.norm([px, py, pz])
    return arctanh(pz/P)

def main():
    # Open the ROOT file
    file_path = Path(__file__).parent.parent / "data" / "exrec.root"
    root_file = ROOT.TFile(str(file_path.resolve()))

    # Making sure tree and files exist
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open file {file_path}")
        return
    
    tree = root_file.Get("particles")

    if not tree:
        print("Error: Could not find tree in file")
        root_file.Close()
        return
    
    # Cut parameters
    proton_py_min = 0.18
    proton_py_max = 0.68
    
    # Create histograms
    xmin = 1.8
    xmax = 3
    bin_size = 0.02 # in GeV
    bins = int((xmax - xmin) / bin_size)
    # K distributions
    hist_K_mass_dist = ROOT.TH1F("hist_K_mass_dist", "K* Generated Mass Distribution", bins, xmin, xmax)
    hist_Kbar_mass_dist = ROOT.TH1F("hist_Kbar_mass_dist", "K*bar Generated Mass Distribution", bins, xmin, xmax)
    # p distributions K
    hist_Kpx_dist =  ROOT.TH1F("hist_Kpx_dist", "K* px", bins, xmin, xmax)
    hist_Kpy_dist =  ROOT.TH1F("hist_Kpy_dist", "K* py", bins, xmin, xmax)
    hist_Kpz_dist =  ROOT.TH1F("hist_Kpz_dist", "K* pz", bins, xmin, xmax)
    hist_Ke_dist =  ROOT.TH1F("hist_Ke_dist", "K* pe", bins, xmin, xmax)

    hist_Kbarpx_dist =  ROOT.TH1F("hist_Kbarpx_dist", "K*bar px", bins, xmin, xmax)
    hist_Kbarpy_dist =  ROOT.TH1F("hist_Kbarpy_dist", "K*bar py", bins, xmin, xmax)
    hist_Kbarpz_dist =  ROOT.TH1F("hist_Kbarpz_dist", "K*bar pz", bins, xmin, xmax)
    hist_Kbare_dist =  ROOT.TH1F("hist_Kbare_dist", "K*bar pe", bins, xmin, xmax)
    # eta rapidity distis
    hist_Keta_dist = ROOT.TH1F("hist_Keta_dist", "K* Pseudorapidity Distribution", bins, xmin, xmax)
    hist_Kbareta_dist = ROOT.TH1F("hist_Kbareta_dist", "K*bar Pseudorapidity Distribution", bins, xmin, xmax)

    # Now same but with cuts
    hist_K_mass_dist_cut =  ROOT.TH1F("hist_K_mass_dist_cut", "K* Generated Mass Distribution (py cut)", bins, xmin, xmax)
    hist_Kbar_mass_dist_cut = ROOT.TH1F("hist_Kbar_mass_dist_cut", "K*bar Generated Mass Distribution (py cut)", bins, xmin, xmax)

    hist_Kpx_dist_cut =  ROOT.TH1F("hist_Kpx_dist_cut", "K* px", bins, xmin, xmax)
    hist_Kpy_dist_cut =  ROOT.TH1F("hist_Kpy_dist_cut", "K* py", bins, xmin, xmax)
    hist_Kpz_dist_cut =  ROOT.TH1F("hist_Kpz_dist_cut", "K* pz", bins, xmin, xmax)
    hist_Ke_dist_cut =  ROOT.TH1F("hist_Ke_dist_cut", "K* pe", bins, xmin, xmax)

    hist_Kbarpx_dist_cut =  ROOT.TH1F("hist_Kbarpx_dist_cut", "K*bar px", bins, xmin, xmax)
    hist_Kbarpy_dist_cut =  ROOT.TH1F("hist_Kbarpy_dist_cut", "K*bar py", bins, xmin, xmax)
    hist_Kbarpz_dist_cut =  ROOT.TH1F("hist_Kbarpz_dist_cut", "K*bar pz", bins, xmin, xmax)
    hist_Kbare_dist_cut =  ROOT.TH1F("hist_Kbare_dist_cut", "K*bar pe", bins, xmin, xmax)

    hist_Keta_dist_cut = ROOT.TH1F("hist_Keta_dist_cut", "K* Pseudorapidity Distribution", bins, xmin, xmax)
    hist_Kbareta_dist_cut = ROOT.TH1F("hist_Kbareta_dist_cut", "K*bar Pseudorapidity Distribution", bins, xmin, xmax)

    # Iterate through entries
    for entry in tree:
        # First two particles are K* and K*bar in this generation
        # Fill K* first
        Km = entry.produced_m[0]
        Kpx = entry.produced_px[0]
        Kpy = entry.produced_py[0]
        Kpz = entry.produced_pz[0]
        Ke = entry.produced_e[0]
        Keta = pseudorapidity(Kpx, Kpy, Kpz)

        hist_K_mass_dist.Fill(Km)
        hist_Kpx_dist.Fill(Kpx)
        hist_Kpy_dist.Fill(Kpy)
        hist_Kpz_dist.Fill(Kpz)
        hist_Ke_dist.Fill(Ke)
        hist_Keta_dist.Fill(Keta)

        # Now K*bar
        Kbarm = entry.produced_m[1]
        Kbarpx = entry.produced_px[1]
        Kbarpy = entry.produced_py[1]
        Kbarpz = entry.produced_pz[1]
        Kbare = entry.produced_e[1]
        Kbareta = pseudorapidity(Kbarpx, Kbarpy, Kbarpz)

        hist_Kbar_mass_dist.Fill(Kbarm)
        hist_Kbarpx_dist.Fill(Kbarpx)
        hist_Kbarpy_dist.Fill(Kbarpy)
        hist_Kbarpz_dist.Fill(Kbarpz)
        hist_Kbare_dist.Fill(Kbare)
        hist_Kbareta_dist.Fill(Kbareta)

        # Now with cuts
        if ((proton_py_min < abs(entry.p1_out_py[0]) < proton_py_max) and (proton_py_min < abs(entry.p2_out_py[0]) < proton_py_max)):
            hist_K_mass_dist_cut.Fill(Km)
            hist_Kpx_dist_cut.Fill(Kpx)
            hist_Kpy_dist_cut.Fill(Kpy)
            hist_Kpz_dist_cut.Fill(Kpz)
            hist_Ke_dist_cut.Fill(Ke)
            hist_Keta_dist_cut.Fill(Keta)
            hist_Kbar_mass_dist_cut.Fill(Kbarm)
            hist_Kbarpx_dist_cut.Fill(Kbarpx)
            hist_Kbarpy_dist_cut.Fill(Kbarpy)
            hist_Kbarpz_dist_cut.Fill(Kbarpz)
            hist_Kbare_dist_cut.Fill(Kbare)
            hist_Kbareta_dist_cut.Fill(Kbareta)


    print("Done!")

    root_file.Close()

    # Mass distributions
    c1 = ROOT.TCanvas("c1", "Mass Distributions", 800, 600)
    c1.Divide(1, 2)
    c1.cd(1)
    hist_K_mass_dist.SetLineColor(ROOT.kRed)
    hist_K_mass_dist.Draw()
    hist_Kbar_mass_dist.SetLineColor(ROOT.kBlue)
    hist_Kbar_mass_dist.Draw("same")
    c1.cd(2)
    hist_K_mass_dist_cut.SetLineColor(ROOT.kRed)
    hist_K_mass_dist_cut.Draw()
    hist_Kbar_mass_dist_cut.SetLineColor(ROOT.kBlue)
    hist_Kbar_mass_dist_cut.Draw("same")

    # Momentum distributions
    c2 = ROOT.TCanvas("c2", "Momentum Distributions", 800, 600)
    c2.Divide(2, 2)
    c2.cd(1)
    hist_Kpx_dist.SetLineColor(ROOT.kRed)
    hist_Kpx_dist.Draw()
    hist_Kbarpx_dist.SetLineColor(ROOT.kBlue)
    hist_Kbarpx_dist.Draw("same")
    hist_Kpx_dist_cut.SetLineColor(ROOT.kGreen)
    hist_Kpx_dist_cut.Draw("same")
    hist_Kbarpx_dist_cut.SetLineColor(ROOT.kMagenta)
    hist_Kbarpx_dist_cut.Draw("same")
    c2.cd(2)
    hist_Kpy_dist.SetLineColor(ROOT.kRed)
    hist_Kpy_dist.Draw()
    hist_Kbarpy_dist.SetLineColor(ROOT.kBlue)
    hist_Kbarpy_dist.Draw("same")
    hist_Kpy_dist_cut.SetLineColor(ROOT.kGreen)
    hist_Kpy_dist_cut.Draw("same")
    hist_Kbarpy_dist_cut.SetLineColor(ROOT.kMagenta)
    hist_Kbarpy_dist_cut.Draw("same")
    c2.cd(3)
    hist_Kpz_dist.SetLineColor(ROOT.kRed)
    hist_Kpz_dist.Draw()
    hist_Kbarpz_dist.SetLineColor(ROOT.kBlue)
    hist_Kbarpz_dist.Draw("same")
    hist_Kpz_dist_cut.SetLineColor(ROOT.kGreen)
    hist_Kpz_dist_cut.Draw("same")
    hist_Kbarpz_dist_cut.SetLineColor(ROOT.kMagenta)
    hist_Kbarpz_dist_cut.Draw("same")
    c2.cd(4)
    hist_Ke_dist.SetLineColor(ROOT.kRed)
    hist_Ke_dist.Draw()
    hist_Kbare_dist.SetLineColor(ROOT.kBlue)
    hist_Kbare_dist.Draw("same")
    hist_Ke_dist_cut.SetLineColor(ROOT.kGreen)
    hist_Ke_dist_cut.Draw("same")
    hist_Kbare_dist_cut.SetLineColor(ROOT.kMagenta)
    hist_Kbare_dist_cut.Draw("same")

    # Pseudorapidities
    c3 = ROOT.TCanvas("c3", "Pseudorapidities", 800, 600)
    hist_Keta_dist.SetLineColor(ROOT.kRed)
    hist_Keta_dist.Draw()
    hist_Kbareta_dist.SetLineColor(ROOT.kBlue)
    hist_Kbareta_dist.Draw("same")
    hist_Keta_dist_cut.SetLineColor(ROOT.kGreen)
    hist_Keta_dist_cut.Draw("same")
    hist_Kbareta_dist_cut.SetLineColor(ROOT.kMagenta)
    hist_Kbareta_dist_cut.Draw("same")

    # Save canvases
    save_path = Path(__file__).parent.parent / "plots"
    c1.SaveAs(save_path / "K_mass_dist.png")
    c1.SaveAs(save_path / "K_mass_dist.root")
    c2.SaveAs(save_path / "K_momentum_dist.png")
    c2.SaveAs(save_path / "K_momentum_dist.root")
    c3.SaveAs(save_path / "K_pseudorapidity_dist.png")
    c3.SaveAs(save_path / "K_pseudorapidity_dist.root")

    # Keep windows open
    input("Press Enter to exit")



if __name__ == "__main__":
    main()