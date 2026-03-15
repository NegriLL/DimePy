#!/usr/bin/env python3

import ROOT
from pathlib import Path

def main():
    # Open the ROOT file
    file_path = Path(__file__).parent.parent / "data" / "exrec.root"
    root_file = ROOT.TFile(str(file_path.resolve()))



if __name__ == "__main__":
    main()