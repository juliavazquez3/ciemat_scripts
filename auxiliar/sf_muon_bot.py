import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

list_names = ["iso","iso_log","z"]

histFile = TFile.Open("/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/final_ratios_muons_all.root","READ")
hist_z = histFile.Get("ratio_bot1plusbot2_nocut_z") 
hist_iso = histFile.Get("ratio_bot1plusbot2_nocut_iso") 
hist_iso_log = histFile.Get("ratio_bot1plusbot2_nocut_iso_log") 

gInterpreter.Declare("""
   Double_t xbins[22] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,20.};
""")

hist_iso_new = hist_iso.Rebin(21,"hist_iso_new",xbins)

myratiofile = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/test_sf.root', 'RECREATE')
myratiofile.WriteObject(hist_iso_new,"test_sf")
myratiofile.Close()
