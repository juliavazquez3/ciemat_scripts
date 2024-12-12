import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--etabin", type=string, default="one",
                    help="etabin to process")
args = parser.parse_args()

if args.etabin == "none": eta_bin = "none"
elif args.etabin == "one": eta_bin = "one"
elif args.etabin == "two": eta_bin = "two"
elif args.etabin == "three": eta_bin = "three"
elif args.etabin == "four": eta_bin = "four"
elif args.etabin == "five": eta_bin = "five"
else: raise NameError('Incorrect channel')
print(eta_bin)

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

datayears = ["all"]
channelslist = ["one_chi","two_chi","one","two"]
list_names = ["iso_abs"]

names_channel = {}
names_channel["one_chi"] = "muon_bot1_";names_channel["two_chi"] = "muon_bot2_";
names_channel["one"] = "muon_bot1_";names_channel["two"] = "muon_bot2_";

histFile = {}
hist_num = {}
hist_den = {}

for chan in channelslist:
   histFile[chan] = {}
   hist_num[chan] = {}
   hist_den[chan] = {}
   for data_op in datayears:
      histFile[chan][data_op] = TFile.Open("/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/etabins/ratios_muons_"+str(data_op)+"_"+str(chan)+"_"+str(eta_bin)+".root","READ")
      hist_num[chan][data_op] = {}
      hist_den[chan][data_op] = {}
      for name in list_names:
         hist_num[chan][data_op][name] = histFile[chan][data_op].Get("ratio_num_"+str(names_channel[chan])+str(name)) 
         hist_den[chan][data_op][name] = histFile[chan][data_op].Get("ratio_den_"+str(names_channel[chan])+str(name)) 

#print(hist_num[chan][data_op].keys())

gInterpreter.Declare("""
   Double_t xbins[23] = {0.,2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,37.5,40.,45.,50.,60.,70.,80.,100.};
""")

myratiofile = {}
myratiofilerebin = {}
hist_bot1cut = {}
hist_bot2cut = {}
hist_bot1nocut = {}
hist_bot2nocut = {}
hist_bot1bot2cut = {}
hist_bot1bot2nocut = {}
hist_bot1cutbot2nocut = {}
hist_bot1nocutbot2cut = {}
hist_bot1nocutbot2nocut = {}
hist_bot1bot2cut_aux = {}
hist_bot1bot2nocut_aux = {}
hist_bot1cutbot2nocut_aux = {}
hist_bot1nocutbot2cut_aux = {}
hist_bot1nocutbot2nocut_aux = {}
for data_op in datayears:
   myratiofile[data_op] = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/etabins/final_ratios_muons_'+str(data_op)+'_'+str(eta_bin)+'.root', 'RECREATE')
   myratiofilerebin[data_op] = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/etabins/final_ratios_muons_'+str(data_op)+'_'+str(eta_bin)+'_rebin.root', 'RECREATE')
   hist_bot1cut[data_op] = {}
   hist_bot2cut[data_op] = {}
   hist_bot1nocut[data_op] = {}
   hist_bot2nocut[data_op] = {}
   hist_bot1bot2cut[data_op] = {}
   hist_bot1bot2nocut[data_op] = {}
   hist_bot1cutbot2nocut[data_op] = {}
   hist_bot1nocutbot2cut[data_op] = {}
   hist_bot1nocutbot2nocut[data_op] = {}
   hist_bot1bot2cut_aux[data_op] = {}
   hist_bot1bot2nocut_aux[data_op] = {}
   hist_bot1cutbot2nocut_aux[data_op] = {}
   hist_bot1nocutbot2cut_aux[data_op] = {}
   hist_bot1nocutbot2nocut_aux[data_op] = {}
   for name in list_names:
      #print(name)
      hist_bot1cut[data_op][name] = hist_num["one_chi"][data_op][name].Clone("ratio_bot1_chi2cut_"+str(name))
      hist_bot1cut[data_op][name].Divide(hist_den["one_chi"][data_op][name])
      hist_bot2cut[data_op][name] = hist_num["two_chi"][data_op][name].Clone("ratio_bot2_chi2cut_"+str(name))
      hist_bot2cut[data_op][name].Divide(hist_den["two_chi"][data_op][name])
      hist_bot1nocut[data_op][name] = hist_num["one"][data_op][name].Clone("ratio_bot1_nocut_"+str(name))
      hist_bot1nocut[data_op][name].Divide(hist_den["one"][data_op][name])
      hist_bot2nocut[data_op][name] = hist_num["two"][data_op][name].Clone("ratio_bot2_nocut_"+str(name))
      hist_bot2nocut[data_op][name].Divide(hist_den["two"][data_op][name])

      myratiofile[data_op].WriteObject(hist_bot1cut[data_op][name],"ratio_bot1_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot2cut[data_op][name],"ratio_bot2_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1nocut[data_op][name],"ratio_bot1_nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot2nocut[data_op][name],"ratio_bot2_nocut_"+str(name))

      hist_bot1bot2cut[data_op][name] = hist_num["one_chi"][data_op][name].Clone("ratio_num_bot1plusbot2_chi2cut_"+str(name))
      hist_bot1bot2cut_aux[data_op][name] = hist_den["one_chi"][data_op][name].Clone("ratio_den_bot1plusbot2_chi2cut_"+str(name))
      hist_bot1bot2cut[data_op][name].Add(hist_num["two_chi"][data_op][name])
      hist_bot1bot2cut_aux[data_op][name].Add(hist_den["two_chi"][data_op][name])
      hist_bot1bot2cut[data_op][name].Divide(hist_bot1bot2cut_aux[data_op][name])

      hist_bot1bot2nocut[data_op][name] = hist_num["one"][data_op][name].Clone("ratio_num_bot1plusbot2_nocut_"+str(name))
      hist_bot1bot2nocut_aux[data_op][name] = hist_den["one"][data_op][name].Clone("ratio_den_bot1plusbot2_nocut_"+str(name))
      hist_bot1bot2nocut[data_op][name].Add(hist_num["two"][data_op][name])
      hist_bot1bot2nocut_aux[data_op][name].Add(hist_den["two"][data_op][name])
      hist_bot1bot2nocut[data_op][name].Divide(hist_bot1bot2nocut_aux[data_op][name])

      hist_bot1cutbot2nocut[data_op][name] = hist_num["one_chi"][data_op][name].Clone("ratio_num_bot1cut_bot2nocut_"+str(name))
      hist_bot1cutbot2nocut_aux[data_op][name] = hist_den["one_chi"][data_op][name].Clone("ratio_den_bot1cut_bot2nocut_"+str(name))
      hist_bot1cutbot2nocut[data_op][name].Add(hist_num["two"][data_op][name])
      hist_bot1cutbot2nocut_aux[data_op][name].Add(hist_den["two"][data_op][name])
      hist_bot1cutbot2nocut[data_op][name].Divide(hist_bot1cutbot2nocut_aux[data_op][name])

      hist_bot1nocutbot2cut[data_op][name] = hist_num["one"][data_op][name].Clone("ratio_num_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut_aux[data_op][name] = hist_den["one"][data_op][name].Clone("ratio_den_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut[data_op][name].Add(hist_num["two_chi"][data_op][name])
      hist_bot1nocutbot2cut_aux[data_op][name].Add(hist_den["two_chi"][data_op][name])
      hist_bot1nocutbot2cut[data_op][name].Divide(hist_bot1nocutbot2cut_aux[data_op][name])

      hist_bot1nocutbot2cut[data_op][name] = hist_num["one"][data_op][name].Clone("ratio_num_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut_aux[data_op][name] = hist_den["one"][data_op][name].Clone("ratio_den_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut[data_op][name].Add(hist_num["two_chi"][data_op][name])
      hist_bot1nocutbot2cut_aux[data_op][name].Add(hist_den["two_chi"][data_op][name])
      hist_bot1nocutbot2cut[data_op][name].Divide(hist_bot1nocutbot2cut_aux[data_op][name])

      hist_bot1nocutbot2nocut[data_op][name] = hist_num["one"][data_op][name].Clone("ratio_num_bot1nocut_bot2nocut_"+str(name))
      hist_bot1nocutbot2nocut_aux[data_op][name] = hist_den["one"][data_op][name].Clone("ratio_den_bot1nocut_bot2nocut_"+str(name))
      hist_bot1nocutbot2nocut[data_op][name].Add(hist_num["two"][data_op][name])
      hist_bot1nocutbot2nocut_aux[data_op][name].Add(hist_den["two"][data_op][name])
      hist_bot1nocutbot2nocut[data_op][name].Divide(hist_bot1nocutbot2nocut_aux[data_op][name])

      myratiofile[data_op].WriteObject(hist_bot1bot2cut[data_op][name],"ratio_bot1plusbot2_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1bot2nocut[data_op][name],"ratio_bot1plusbot2_nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1cutbot2nocut[data_op][name],"ratio_bot1cut_bot2nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1nocutbot2cut[data_op][name],"ratio_bot1nocut_bot2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1nocutbot2nocut[data_op][name],"ratio_bot1nocut_bot2nocut_"+str(name))

      hist_1_rebin = hist_num["one"][data_op][name].Clone("ratio_rebin_1_"+str(name))
      hist_2_rebin = hist_den["one"][data_op][name].Clone("ratio_rebin_2_"+str(name))
      hist_3_rebin = hist_num["two"][data_op][name].Clone("ratio_rebin_3_"+str(name))
      hist_4_rebin = hist_den["two"][data_op][name].Clone("ratio_rebin_4_"+str(name))

      hist_1_rebin = hist_1_rebin.Rebin(22,"hist_test_bin_1_new",xbins)
      hist_2_rebin = hist_2_rebin.Rebin(22,"hist_test_bin_1_new",xbins)
      hist_3_rebin = hist_3_rebin.Rebin(22,"hist_test_bin_1_new",xbins)
      hist_4_rebin = hist_4_rebin.Rebin(22,"hist_test_bin_1_new",xbins)

      hist_1_rebin.Add(hist_3_rebin)
      hist_2_rebin.Add(hist_4_rebin)
      hist_1_rebin.Divide(hist_2_rebin)
      myratiofilerebin[data_op].WriteObject(hist_1_rebin,"ratio_bot1bot2_nochi_rebin")

   myratiofile[data_op].Close()
   myratiofilerebin[data_op].Close()
