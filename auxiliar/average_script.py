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
parser.add_argument("--rebin", type=string, default="v33",
                    help="Select year of process to run")
args = parser.parse_args()

if args.rebin == "og": nbins = 21
elif args.rebin == "v1": nbins = 22
elif args.rebin == "v2": nbins = 23
elif args.rebin == "v3": nbins = 24
elif args.rebin == "v33": nbins = 22
elif args.rebin == "v4": nbins = 25
elif args.rebin == "v5": nbins = 30

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

datayears = ["all"]
channelslist = ["jet_bot1_chi","jet_bot2_chi","jet_bot1_nochi","jet_bot2_nochi"]
list_names = ["iso","iso_log","z"]

names_channel = {}
names_channel["jet_bot1_chi"] = "muon_bot1_";names_channel["jet_bot2_chi"] = "muon_bot2_";
names_channel["jet_bot1_nochi"] = "muon_bot1_";names_channel["jet_bot2_nochi"] = "muon_bot2_";

histFile = {}
hist_num = {}
hist_den = {}

for chan in channelslist:
   histFile[chan] = {}
   hist_num[chan] = {}
   hist_den[chan] = {}
   for data_op in datayears:
      histFile[chan][data_op] = TFile.Open("/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/ratios_muons_"+str(data_op)+"_"+str(chan)+".root","READ")
      hist_num[chan][data_op] = {}
      hist_den[chan][data_op] = {}
      for name in list_names:
         hist_num[chan][data_op][name] = histFile[chan][data_op].Get("ratio_num_"+str(names_channel[chan])+str(name)) 
         hist_den[chan][data_op][name] = histFile[chan][data_op].Get("ratio_den_"+str(names_channel[chan])+str(name)) 

gInterpreter.Declare("""
   Double_t xbins[22] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,20.};
   Double_t xbinsv1[23] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,15.,20.};
   Double_t xbinsv2[24] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,12.5,15.,20.};
   Double_t xbinsv3[25] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,12.5,15.,17.5,20.};
   Double_t xbinsv33[23] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.5,15.,17.5,20.};
   Double_t xbinsv4[26] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,12.,14.,16.,18.,20.};
   Double_t xbinsv5[31] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.};
   Double_t xbinsz[25] = {0.,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.};
""")

myratiofile = {}
hist_bot1cut = {}
hist_bot2cut = {}
hist_bot1nocut = {}
hist_bot2nocut = {}
hist_bot1bot2cut = {}
hist_bot1bot2nocut = {}
hist_bot1cutbot2nocut = {}
hist_bot1nocutbot2cut = {}
hist_bot1bot2cut_aux = {}
hist_bot1bot2nocut_aux = {}
hist_bot1cutbot2nocut_aux = {}
hist_bot1nocutbot2cut_aux = {}
for data_op in datayears:
   myratiofile[data_op] = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/final_ratios_muons_'+str(data_op)+'.root', 'RECREATE')
   hist_bot1cut[data_op] = {}
   hist_bot2cut[data_op] = {}
   hist_bot1nocut[data_op] = {}
   hist_bot2nocut[data_op] = {}
   hist_bot1bot2cut[data_op] = {}
   hist_bot1bot2nocut[data_op] = {}
   hist_bot1cutbot2nocut[data_op] = {}
   hist_bot1nocutbot2cut[data_op] = {}
   hist_bot1bot2cut_aux[data_op] = {}
   hist_bot1bot2nocut_aux[data_op] = {}
   hist_bot1cutbot2nocut_aux[data_op] = {}
   hist_bot1nocutbot2cut_aux[data_op] = {}
   for name in list_names:
      #print(name)
      hist_bot1cut[data_op][name] = hist_num["jet_bot1_chi"][data_op][name].Clone("ratio_bot1_chi2cut_"+str(name))
      hist_bot1cut[data_op][name].Divide(hist_den["jet_bot1_chi"][data_op][name])
      hist_bot2cut[data_op][name] = hist_num["jet_bot2_chi"][data_op][name].Clone("ratio_bot2_chi2cut_"+str(name))
      hist_bot2cut[data_op][name].Divide(hist_den["jet_bot2_chi"][data_op][name])
      hist_bot1nocut[data_op][name] = hist_num["jet_bot1_nochi"][data_op][name].Clone("ratio_bot1_nocut_"+str(name))
      hist_bot1nocut[data_op][name].Divide(hist_den["jet_bot1_nochi"][data_op][name])
      hist_bot2nocut[data_op][name] = hist_num["jet_bot2_nochi"][data_op][name].Clone("ratio_bot2_nocut_"+str(name))
      hist_bot2nocut[data_op][name].Divide(hist_den["jet_bot2_nochi"][data_op][name])

      myratiofile[data_op].WriteObject(hist_bot1cut[data_op][name],"ratio_bot1_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot2cut[data_op][name],"ratio_bot2_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1nocut[data_op][name],"ratio_bot1_nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot2nocut[data_op][name],"ratio_bot2_nocut_"+str(name))

      hist_bot1bot2cut[data_op][name] = hist_num["jet_bot1_chi"][data_op][name].Clone("ratio_num_bot1plusbot2_chi2cut_"+str(name))
      hist_bot1bot2cut_aux[data_op][name] = hist_den["jet_bot1_chi"][data_op][name].Clone("ratio_den_bot1plusbot2_chi2cut_"+str(name))
      hist_bot1bot2cut[data_op][name].Add(hist_num["jet_bot2_chi"][data_op][name])
      hist_bot1bot2cut_aux[data_op][name].Add(hist_den["jet_bot2_chi"][data_op][name])
      hist_bot1bot2cut[data_op][name].Divide(hist_bot1bot2cut_aux[data_op][name])

      hist_bot1bot2nocut[data_op][name] = hist_num["jet_bot1_nochi"][data_op][name].Clone("ratio_num_bot1plusbot2_nocut_"+str(name))
      hist_bot1bot2nocut_aux[data_op][name] = hist_den["jet_bot1_nochi"][data_op][name].Clone("ratio_den_bot1plusbot2_nocut_"+str(name))
      hist_bot1bot2nocut[data_op][name].Add(hist_num["jet_bot2_nochi"][data_op][name])
      hist_bot1bot2nocut_aux[data_op][name].Add(hist_den["jet_bot2_nochi"][data_op][name])
      hist_bot1bot2nocut[data_op][name].Divide(hist_bot1bot2nocut_aux[data_op][name])

      hist_bot1cutbot2nocut[data_op][name] = hist_num["jet_bot1_chi"][data_op][name].Clone("ratio_num_bot1cut_bot2nocut_"+str(name))
      hist_bot1cutbot2nocut_aux[data_op][name] = hist_den["jet_bot1_chi"][data_op][name].Clone("ratio_den_bot1cut_bot2nocut_"+str(name))
      hist_bot1cutbot2nocut[data_op][name].Add(hist_num["jet_bot2_nochi"][data_op][name])
      hist_bot1cutbot2nocut_aux[data_op][name].Add(hist_den["jet_bot2_nochi"][data_op][name])
      hist_bot1cutbot2nocut[data_op][name].Divide(hist_bot1cutbot2nocut_aux[data_op][name])

      hist_bot1nocutbot2cut[data_op][name] = hist_num["jet_bot1_nochi"][data_op][name].Clone("ratio_num_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut_aux[data_op][name] = hist_den["jet_bot1_nochi"][data_op][name].Clone("ratio_den_bot1nocut_bot2cut_"+str(name))
      hist_bot1nocutbot2cut[data_op][name].Add(hist_num["jet_bot2_chi"][data_op][name])
      hist_bot1nocutbot2cut_aux[data_op][name].Add(hist_den["jet_bot2_chi"][data_op][name])
      hist_bot1nocutbot2cut[data_op][name].Divide(hist_bot1nocutbot2cut_aux[data_op][name])

      myratiofile[data_op].WriteObject(hist_bot1bot2cut[data_op][name],"ratio_bot1plusbot2_chi2cut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1bot2nocut[data_op][name],"ratio_bot1plusbot2_nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1cutbot2nocut[data_op][name],"ratio_bot1cut_bot2nocut_"+str(name))
      myratiofile[data_op].WriteObject(hist_bot1nocutbot2cut[data_op][name],"ratio_bot1nocut_bot2cut_"+str(name))

   if args.rebin == "og": x_bins = xbins
   elif args.rebin == "v1": x_bins = xbinsv1
   elif args.rebin == "v2": x_bins = xbinsv2
   elif args.rebin == "v3": x_bins = xbinsv3
   elif args.rebin == "v33": x_bins = xbinsv33
   elif args.rebin == "v4": x_bins = xbinsv4
   elif args.rebin == "v5": x_bins = xbinsv5

   hist_test_bin_1 = hist_num["jet_bot1_nochi"][data_op]["iso"].Clone("ratio_num_bot1_nochi_iso")
   hist_test_bin_2 = hist_den["jet_bot1_nochi"][data_op]["iso"].Clone("ratio_den_bot1_nochi_iso")
   hist_test_bin_3 = hist_num["jet_bot2_nochi"][data_op]["iso"].Clone("ratio_num_bot2_nochi_iso")
   hist_test_bin_4 = hist_den["jet_bot2_nochi"][data_op]["iso"].Clone("ratio_den_bot2_nochi_iso")
   hist_test_bin_1_new = hist_test_bin_1.Rebin(nbins,"hist_test_bin_1_new",x_bins)
   hist_test_bin_2_new = hist_test_bin_2.Rebin(nbins,"hist_test_bin_2_new",x_bins)
   hist_test_bin_3_new = hist_test_bin_3.Rebin(nbins,"hist_test_bin_3_new",x_bins)
   hist_test_bin_4_new = hist_test_bin_4.Rebin(nbins,"hist_test_bin_4_new",x_bins)
   hist_test_bin_1_new.Add(hist_test_bin_3_new)
   hist_test_bin_2_new.Add(hist_test_bin_4_new)
   hist_test_bin_1_new.Divide(hist_test_bin_2_new)
   myratiofile[data_op].WriteObject(hist_test_bin_1_new,"ratio_bot1bot2_nochi_iso_rebin")

   hist_rebin_1 = hist_num["jet_bot1_chi"][data_op]["iso"].Clone("ratio_num_bot1_chi_iso")
   hist_rebin_2 = hist_den["jet_bot1_chi"][data_op]["iso"].Clone("ratio_den_bot1_chi_iso")
   hist_rebin_3 = hist_num["jet_bot2_chi"][data_op]["iso"].Clone("ratio_num_bot2_chi_iso")
   hist_rebin_4 = hist_den["jet_bot2_chi"][data_op]["iso"].Clone("ratio_den_bot2_chi_iso")
   hist_rebin_1_new = hist_rebin_1.Rebin(nbins,"hist_rebin_1_new",x_bins)
   hist_rebin_2_new = hist_rebin_2.Rebin(nbins,"hist_rebin_2_new",x_bins)
   hist_rebin_3_new = hist_rebin_3.Rebin(nbins,"hist_rebin_3_new",x_bins)
   hist_rebin_4_new = hist_rebin_4.Rebin(nbins,"hist_rebin_4_new",x_bins)
   hist_rebin_1_new.Add(hist_rebin_3_new)
   hist_rebin_2_new.Add(hist_rebin_4_new)
   hist_rebin_1_new.Divide(hist_rebin_2_new)
   myratiofile[data_op].WriteObject(hist_rebin_1_new,"ratio_bot1bot2_chi_iso_rebin")

   hist_bin_1 = hist_num["jet_bot1_nochi"][data_op]["iso"].Clone("ratio_num_bot1_nochi_iso")
   hist_bin_2 = hist_den["jet_bot1_nochi"][data_op]["iso"].Clone("ratio_den_bot1_nochi_iso")
   hist_bin_3 = hist_num["jet_bot2_chi"][data_op]["iso"].Clone("ratio_num_bot2_chi_iso")
   hist_bin_4 = hist_den["jet_bot2_chi"][data_op]["iso"].Clone("ratio_den_bot2_chi_iso")
   hist_bin_1_new = hist_bin_1.Rebin(nbins,"hist_bin_1_new",x_bins)
   hist_bin_2_new = hist_bin_2.Rebin(nbins,"hist_bin_2_new",x_bins)
   hist_bin_3_new = hist_bin_3.Rebin(nbins,"hist_bin_3_new",x_bins)
   hist_bin_4_new = hist_bin_4.Rebin(nbins,"hist_bin_4_new",x_bins)
   hist_bin_1_new.Add(hist_bin_3_new)
   hist_bin_2_new.Add(hist_bin_4_new)
   hist_bin_1_new.Divide(hist_bin_2_new)
   myratiofile[data_op].WriteObject(hist_bin_1_new,"ratio_bot1nocut_bot2cut_iso_rebin")

   hist_z_bin_1 = hist_num["jet_bot1_nochi"][data_op]["z"].Clone("ratio_num_bot1_nochi_z")
   hist_z_bin_2 = hist_den["jet_bot1_nochi"][data_op]["z"].Clone("ratio_den_bot1_nochi_z")
   hist_z_bin_3 = hist_num["jet_bot2_chi"][data_op]["z"].Clone("ratio_num_bot2_chi_z")
   hist_z_bin_4 = hist_den["jet_bot2_chi"][data_op]["z"].Clone("ratio_den_bot2_chi_z")
   hist_z_bin_1_new = hist_z_bin_1.Rebin(24,"hist_z_bin_1_new",xbinsz)
   hist_z_bin_2_new = hist_z_bin_2.Rebin(24,"hist_z_bin_2_new",xbinsz)
   hist_z_bin_3_new = hist_z_bin_3.Rebin(24,"hist_z_bin_3_new",xbinsz)
   hist_z_bin_4_new = hist_z_bin_4.Rebin(24,"hist_z_bin_4_new",xbinsz)
   hist_z_bin_1_new.Add(hist_z_bin_3_new)
   hist_z_bin_2_new.Add(hist_z_bin_4_new)
   hist_z_bin_1_new.Divide(hist_z_bin_2_new)
   myratiofile[data_op].WriteObject(hist_z_bin_1_new,"ratio_bot1nocut_bot2cut_z_rebin")

   myratiofile[data_op].Close()
