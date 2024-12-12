import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir
import CMS_lumi, tdrstyle

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

#set the tdr style
#tdrstyle.setTDRStyle()
gStyle.SetPadGridX(False)
gStyle.SetPadGridY(False)
#gStyle.SetGridStyle(3)

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
#CMS_lumi.extraText = ""

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--postfit", action="store_true", default=False,
                    help="Plots normalization of signal")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalising test")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--wcs", action="store_true", default=False,
                    help="classification for ttbar and st")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--folder", type=string, default="",
                    help="folder to grab files from")
parser.add_argument("--charmtag", type=string, default="no",
                    help="Perform ssos substraction")
parser.add_argument("--presel", type=string, default="btagMM_chitest",
                    help="preselection")
parser.add_argument("--nosyst", action="store_true", default=False,
                    help="systematics inclusion")
parser.add_argument("--sumEM", action="store_true", default=False,
                    help="union of M and E channels")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

if args.charmtag == "full": channel = "sl_full"
elif args.charmtag == "ss": channel = "sl_ss"
elif args.charmtag == "os": channel = "sl_os"
elif args.charmtag == "ssos": channel = "sl_ssos"
elif args.charmtag == "csv": channel = "csv"
elif args.charmtag == "no": channel = "wqq"
elif args.charmtag == "antisl": channel = "antisl"
else: raise NameError('Incorrect channel')
print(channel)

if args.presel == "nobtag": chan = "nobtag"
elif args.presel == "btagMM": chan = "btagMM"
elif args.presel == "lepton50": chan = "lepton50"
elif args.presel == "btagMM_chitest": chan = "btagMM_chitest"
elif args.presel == "lepton50_chitest": chan = "lepton50_chitest"
else: raise NameError('Incorrect channel')
print(chan)

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
###### if you want to specify the folder in /nfs or any other location
#### script_dir = ''

######### folder name
term1 = ""
if args.wcs:
   term1 = "hists_wcs_"+str(chan)
else:
   term1 = "hists_jetflav_"+str(chan)

path1 = os.path.join(script_dir, term1)

if "sl" in channel:
   termm = str(channel)
   folder_path = os.path.join(path1,termm)
else:
   folder_path = path1

sl_channel = ["sl_full","sl_ss","sl_os","sl_ssos"]

plotdir = script_dir

##### some settings ######

iPeriod = 0
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

c_rat = 1.5
c_rat2 = 0.5
nrebin = 2; nrat = 11;
mrk_size = 6; 
############## titY_off2 = 0.6; titY_size2=0.07;
titY_off2 = 0.9; titY_size2=0.055;
titY_off = titY_off2/1.857; titY_size=titY_size2*1.857;
############## tit_size = titY_size/0.99; titX_off = titY_off2*1.4;
tit_size = titY_size/0.99; titX_off = titY_off2*1.1;
titY_lab=0.055; titY2_lab=titY_lab*1.857;titX_lab=titY_lab*1.857;
if not args.nodata:
   leg1 = 0.48; leg2 = 0.7; leg3 = 0.94; leg4 = 0.9;
else:
   leg1 = 0.4; leg2 = 0.76; leg3 = 0.88; leg4 = 0.89;

if args.nodata: titY_size = 0.07;titY_off=0.5;
if not args.ratio: 
    titY_size2=0.04; titY_off2 = 1.2;
    titX_size= 0.05; titX_off = 1.1;
    titY_lab=0.04; titY2_lab=titY_lab;titX_lab=titY_lab;

aux_ratio = 0.977; aux_rwqq = 0.999; aux_conratio = 1.022;

syst_val = {}

syst_val["btagMM_chitest"] = 0.065; syst_val["btagMM_chitest_antisl"] = 0.065; syst_val["btagMM_chitest_sl"] = 0.05;
syst_val["btagMM_chitest_slssos"] = 0.066; syst_val["btagMM_chitest_slos"] = 0.05; syst_val["btagMM_chitest_slss"] = 0.05;

##########################

norm_factorM = {}
norm_factorE = {}

norm_factorM["all"] = {}; norm_factorM["2016"] = {}; norm_factorM["2016B"] = {}; norm_factorM["2017"] = {}; norm_factorM["2018"] = {};
norm_factorE["all"] = {}; norm_factorE["2016"] = {}; norm_factorE["2016B"] = {}; norm_factorE["2017"] = {}; norm_factorE["2018"] = {};

### btagMM smeared, chitest
norm_factorM["all"]["wqq"] = 0.915;
norm_factorE["all"]["wqq"] = 0.894;

norm_factorM["2016"]["wqq"] = 0.89;norm_factorE["2016"]["wqq"] = 0.87;
norm_factorM["2016B"]["wqq"] = 0.92;norm_factorE["2016B"]["wqq"] = 0.94;
norm_factorM["2017"]["wqq"] = 0.92;norm_factorE["2017"]["wqq"] = 0.88;
norm_factorM["2018"]["wqq"] = 0.91;norm_factorE["2018"]["wqq"] = 0.90;

### btagMM smeared, chitest, SL
norm_factorM["all"]["sl_full"] = 0.92; norm_factorM["all"]["sl_ss"] = 0.92;
norm_factorM["all"]["sl_os"] = 0.92; norm_factorM["all"]["sl_ssos"] = 0.92;
norm_factorE["all"]["sl_full"] = 0.90; norm_factorE["all"]["sl_ss"] = 0.90;
norm_factorE["all"]["sl_os"] = 0.90; norm_factorE["all"]["sl_ssos"] = 0.90;

norm_factorM["all"]["antisl"] = 0.92; norm_factorE["all"]["antisl"] = 0.9;

########### pT jet not corrected
#norm_factorM["all"]["btagMM_chitest"] = 0.930;
#norm_factorE["all"]["btagMM_chitest"] = 0.909;
#norm_factorM["2016"]["btagMM_chitest"] = 0.91;norm_factorE["2016"]["btagMM_chitest"] = 0.89;
#norm_factorM["2016B"]["btagMM_chitest"] = 0.93;norm_factorE["2016B"]["btagMM_chitest"] = 0.95;
#norm_factorM["2017"]["btagMM_chitest"] = 0.94;norm_factorE["2017"]["btagMM_chitest"] = 0.89;
#norm_factorM["2018"]["btagMM_chitest"] = 0.93;norm_factorE["2018"]["btagMM_chitest"] = 0.92;

observable_names = ["InvM_2jets","nJetGood","jet_1_pt", "jet_2_pt", "jet_1_eta", "jet_1_nmu", "jet_2_eta", "jet_2_mass", "jet_2_qgl","jet_2_nmu","jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_bot_closer", "InvM_bot_farther",
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt_aux", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "pT_Wlep",
   "jet_1_btag", "jet_2_btag", "deltaphi_MET_lep", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber",
   "jet_1_cvltag_csv", "jet_2_cvltag_csv", "jet_1_cvltag", "jet_2_cvltag","InvM30","InvM31","InvMl0","InvMl1","chi2_test0","chi2_test1",
   "InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","chi2_test_good","chi2_test_bad","jet_max_cvltag","jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag",
   "jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
   "InvM_2jets_thick","InvM_2jets_short","bot1_muons","bot2_muons","muon_bot1_eta","muon_bot2_eta","muon_bot1_pt","muon_bot2_pt",
   "jet_bot1_tracks","jet_bot2_tracks", "jet_1_pt_long", "jet_2_pt_long", "jet_bot1_pt_long", "jet_bot2_pt_long",
   "lepton_pt_long","transverse_mass_long","MET_pt_aux_long"]


observable_names = ["InvM_2jets","nJetGood","jet_1_pt", "jet_2_pt", "jet_1_eta", "jet_2_eta","lepton_pt","lepton_eta","lepton_eta_thick",
     "deltaR_jet1_jet2", "MET_pt_aux","transverse_mass","pT_Wlep","deltaphi_MET_lep","jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
     "InvM3_good","InvMl_good","chi2_test_good","jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
     "InvM_2jets_short","MET_my_sig", "jet_1_pt_long", "jet_2_pt_long", "jet_bot1_pt_long", "jet_bot2_pt_long",
     "lepton_pt_long","transverse_mass_long","MET_pt_aux_long","jet_bot1_btag","jet_bot2_btag","jet_1_btagnumber","jet_2_btagnumber"]


if (channel in sl_channel): observable_names = ["nJetGood","jet_1_pt", "jet_2_pt", "jet_1_eta", 
     "jet_2_eta","lepton_pt","lepton_eta_thick","jet_bot1_btag","jet_bot2_btag",
     "deltaR_jet1_jet2", "MET_pt_aux","transverse_mass","deltaphi_MET_lep",
     "InvM3_good","InvMl_good","chi2_test_good","jet_bot1_eta_thick","jet_bot2_eta_thick",
     "InvM_2jets_short","MET_my_sig", "jet_bot1_pt_long", "jet_bot2_pt_long",
     "muon_jet_pt","muon_jet_eta",
     "muon_jet_xy","muon_jet_dz","muon_jet_z_short","muon_jet_iso_abs"]

cosmetic_names = {}
axisY_label_name = {}

for name in observable_names:
    cosmetic_names[name] = name
    axisY_label_name[name] = "Events"

cosmetic_names["InvM_2jets_short"] = "m_{jj}^{W} [GeV]";cosmetic_names["nJetGood"] = "Number of jets";
cosmetic_names["InvM_2jets"] = "m_{jj}^{W} [GeV]"; cosmetic_names["InvM3_good"] = "m_{jjb}^{t} [GeV]"; cosmetic_names["InvMl_good"] = "m_{lb} [GeV]"; 
cosmetic_names["jet_1_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_2_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_bot1_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_bot2_pt"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_1_pt_long"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_2_pt_long"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_bot1_pt_long"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_bot2_pt_long"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_1_eta"] = "\eta^{ jet}";cosmetic_names["jet_2_eta"] = "\eta^{ jet}";cosmetic_names["jet_bot1_eta"] = "\eta^{ jet}";
cosmetic_names["jet_bot2_eta"] = "\eta^{ jet}";
cosmetic_names["jet_max_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_min_pt"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_max_eta"] = "\eta^{ jet}";cosmetic_names["jet_min_eta"] = "\eta^{ jet}";
cosmetic_names["jet_1_flavourP"] = "Parton flavour";cosmetic_names["jet_2_flavourP"] = "Parton flavour";
cosmetic_names["jet_bot1_flavourP"] = "Parton flavour";cosmetic_names["jet_bot2_flavourP"] = "Parton flavour";
cosmetic_names["lepton_pt"] = "p_{T}^{l} [GeV]";cosmetic_names["lepton_eta"] = "\eta^{l}"; cosmetic_names["lepton_eta_thick"] = "\eta^{l}";
cosmetic_names["chi2_test_good"] = "\chi^{2}";cosmetic_names["transverse_mass"] = "m_{T} [GeV]"; cosmetic_names["MET_pt_aux"] = "MET p_{T} [GeV]";
cosmetic_names["jet_bot1_btag"] = "DeepJet (b)"; cosmetic_names["jet_bot2_btag"] = "DeepJet (b)";
cosmetic_names["jet_bot1_btagnumber"] = "DeepJet WP (b)"; cosmetic_names["jet_bot2_btagnumber"] = "DeepJet WP (b)";
cosmetic_names["jet_1_btagnumber"] = "DeepJet WP (b)"; cosmetic_names["jet_2_btagnumber"] = "DeepJet WP (b)";
cosmetic_names["deltaR_jet1_jet2"] = "\Delta R (jet 1, jet 2)"; cosmetic_names["deltaphi_jet1_jet2"] = "\Delta \phi (jet 1, jet 2)";
cosmetic_names["deltaR_jet1_jet2"] = "\Delta R (jet 1, jet 2)"; cosmetic_names["deltaphi_MET_lep"] = "\Delta \phi (l, MET)";
cosmetic_names["pT_Wlep"] = "W^{l} p_{T} [GeV]";
cosmetic_names["muon_jet_pt"] = "p_{T}^{#mu} [GeV]"; cosmetic_names["muon_jet_eta"] = "#eta^{#mu}"; cosmetic_names["muon_jet_iso_abs"] = "I^{#mu} [GeV]";
cosmetic_names["muon_jet_z_short"] = "p_{T}^{#mu}/p_{T}^{jet}";cosmetic_names["muon_jet_z"] = "p_{T}^{#mu}/p_{T}^{jet}";
cosmetic_names["jet_1_eta_thick"] = "\eta^{ jet}";cosmetic_names["jet_2_eta_thick"] = "\eta^{ jet}";
cosmetic_names["jet_bot1_eta_thick"] = "\eta^{ jet}";cosmetic_names["jet_bot2_eta_thick"] = "\eta^{ jet}";
cosmetic_names["muon_jet_iso"] = "I_{PF}/p_{T}^{#mu}";
cosmetic_names["muon_jet_pt_short"] = "p_{T}^{#mu} [GeV]"; cosmetic_names["muon_jet_eta_short"] = "#eta^{#mu}";
cosmetic_names["lepton_pt_long"] = "p_{T}^{l} [GeV]";cosmetic_names["transverse_mass_long"] = "m_{T} [GeV]";cosmetic_names["MET_pt_aux_long"] = "MET p_{T} [GeV]";
cosmetic_names["jet_1_flavourP"] = "Jet Flavour";cosmetic_names["jet_2_flavourP"] = "Jet Flavour";cosmetic_names["jet_bot1_flavourP"] = "Jet Flavour";cosmetic_names["jet_bot2_flavourP"] = "Jet Flavour";

### Ptmiss, MT, pt y eta de los bjets, btag discriminant cvl for the b-jets, pt y eta
### de los dos jets del W, deltaR(jet1W, jet2W), number of b-tag jets, deltaphi (lepton, MET), pt
### leptonic W

not_rebin = ["nJetGood","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_pt",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP",
      "tau_discr_jetbot2","muon_jet_iso","jet_muon_btagnumber","jet_notmuon_btagnumber",
      "jet_1_nmu","jet_2_nmu","jet_muon_nmu","jet_notmuon_nmu","ttsl_lepflav","ttdl_lepflav",
      "muon_jet_pt","muon_jet_eta","muon_jet_pt_short","muon_jet_eta_short","muon_jet_z_short",
      "jet_1_eta_thick","jet_2_eta_thick"]

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

if args.wcs:
   samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz"]

samples_d = ["2016","2016B","2017","2018"]

## lumi info
lumi = {}
xsecs = {}
nevents = {}

for data_op in samples_d:
        files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+data_op+".json"))
        lumi[data_op] = {}
        for p in samples:
                #num_events = files[p]["events"] # Number of events
                num_files = files[p]["files"] # Number of files
                luminosity = files[p]["lumi"] # Luminosity
                #print(files[p]["type"])
                lumi[data_op][p] = luminosity
        lumi[data_op]["ttbar_sl_charm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_nocharm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_light"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottom"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_else"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_charmgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottomgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["st_1_charm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_charm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_charm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_charm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_nocharm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_nocharm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_nocharm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_nocharm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_else"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_else"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_else"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_else"] = lumi[data_op]["st_4"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_charm", "ttbar_sl_charmgluon",
        "ttbar_sl_bottom", "ttbar_sl_bottomgluon", "ttbar_sl_else", "ttbar_sl_light", "ttbar_sl_nocharm",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

for s in listsampl:
   lumi["2016B"][s] = lumi["2016"][s]

lumi_d = {}
lumi_d["2016"] = 19.5
lumi_d["2016B"] = 16.8
lumi_d["2017"] = 41.5
lumi_d["2018"] = 59.8

#######################################################
########### Start of plot creation ####################
#######################################################

if args.year == "all": datayears = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): datayears = [str(args.year)]
else: raise NameError('Incorrect year')

histFile = {}
histFileDM = {}
histFileDE = {}

for name in observable_names:
  if name in ["muon_jet_iso","deltaphi_MET_lep"]: 
     nrebin = 5;
  else:
     nrebin = 2;
  ## Open hists files
  filePath = folder_path
  term = "hist_wqqfromJF_"
  end_term = ".root"
  ## mc files
  histFile[name] = {}
  histFileDM[name] = {}
  histFileDE[name] = {}
  for data_op in datayears:
       if isfile(filePath + term+"MC_"+data_op+"_"+name+end_term):
           histFile[name][data_op] = TFile.Open(filePath + term+"MC_"+data_op+"_"+name+end_term,"READ")
       # data files
       if not args.nodata: histFileDM[name][data_op] = TFile.Open(filePath + term+"dataM_"+data_op+"_"+name+".root","READ")
       if not args.nodata: histFileDE[name][data_op] = TFile.Open(filePath + term+"dataE_"+data_op+"_"+name+".root","READ")
  #print(data_op)
  #print(histFile[data_op].keys())
  #print(histFileDM[data_op].keys())

  if args.png: c1 = TCanvas("c1","",1200,1000)
  else: c1 = TCanvas("c1","",600,400)

  if args.ratio:
    ## In case of ratio plot
    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.095)
        p.SetRightMargin(0.04)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    upper_pad.SetTopMargin(0.075)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.25)

    if not args.linear: upper_pad.SetLogy()
    upper_pad.Draw()
    lower_pad.Draw()

  samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4","ttbar_sl_else","ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

  if args.wcs: 
    samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
            "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
            "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
            "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
            "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

  ## HISTS
  samples_foryear = {}
  hist_nom_M = {}
  hist_nom_E = {}
  hist_syst_M = {}
  hist_syst_E = {}
  #list_syst = ["btaglightup","btaglightdown","btagheavyup","btagheavydown","seclepup","seclepdown"]
  list_syst = ["ctagup","ctagdown"]
  for syst in list_syst:
      hist_syst_M[syst] = {}
      hist_syst_E[syst] = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    hist_nom_M[data_op] = {}
    hist_nom_E[data_op] = {}
    for syst in list_syst:
       hist_syst_M[syst][data_op] = {}
       hist_syst_E[syst][data_op] = {}
    data_term = data_op
    #print(data_op+"M_"+name+"_M")
    #print(histFile[data_op][name].ls())
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      elif s[0:2] == "st":
         s_term = s[0:4]+data_term+s[4:] 
      else:
         s_term = s+data_term
      hist_nom_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M")
      hist_nom_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E")
      if False:
         for syst in list_syst:
             hist_syst_M[syst][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_"+str(syst))
             hist_syst_E[syst][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_"+str(syst))
    if not args.nodata: histdata_M[data_op] = histFileDM[name][data_op].Get("data"+data_op+"M_"+name+"_M")
    if not args.nodata: histdata_E[data_op] = histFileDE[name][data_op].Get("data"+data_op+"E_"+name+"_E")

  samples_st = {}
  samples_wjets = {}
  samples_zjets = {}
  ## Scaling to lumi
  #print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      #print(data_op)
      hist_nom_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_nom_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if False:
         for syst in list_syst:
           hist_syst_M[syst][data_op][s].Scale(lumi_data/lumi[data_op][s])
           hist_syst_E[syst][data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(channel)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(channel)])
      if False:
         for syst in list_syst:
             if args.norm: hist_syst_M[syst][data_op][s].Scale(norm_factorM[str(args.year)][str(channel)])
             if args.norm: hist_syst_E[syst][data_op][s].Scale(norm_factorE[str(args.year)][str(channel)])
    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    if args.wcs:
       list_st_aux = ["_charm","_nocharm","_else"]
    else:
       list_st_aux = [""]
    for apx in list_st_aux:
      hist_nom_M[data_op]["st"+apx] = hist_nom_M[data_op][list_st[0]+apx]
      hist_nom_E[data_op]["st"+apx] = hist_nom_E[data_op][list_st[0]+apx]
      if False:
         for syst in list_syst:
              hist_syst_M[syst][data_op]["st"+apx] = hist_syst_M[syst][data_op][list_st[0]+apx]
              hist_syst_E[syst][data_op]["st"+apx] = hist_syst_E[syst][data_op][list_st[0]+apx]
      for l in list_st[1:]:
        hist_nom_M[data_op]["st"+apx].Add(hist_nom_M[data_op][l+apx])
        hist_nom_E[data_op]["st"+apx].Add(hist_nom_E[data_op][l+apx])
        if False:
           for syst in list_syst:
              hist_syst_M[syst][data_op]["st"+apx].Add(hist_syst_M[syst][data_op][l+apx])
              hist_syst_E[syst][data_op]["st"+apx].Add(hist_syst_E[syst][data_op][l+apx])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    if False:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["wjets"] = hist_syst_M[syst][data_op][list_wjets[0]]
            hist_syst_E[syst][data_op]["wjets"] = hist_syst_E[syst][data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      if False:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["wjets"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["wjets"].Add(hist_syst_E[syst][data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    if False:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["zjets"] = hist_syst_M[syst][data_op][list_zjets[0]]
            hist_syst_E[syst][data_op]["zjets"] = hist_syst_E[syst][data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      if False:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["zjets"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["zjets"].Add(hist_syst_E[syst][data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    if False:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["vv"] = hist_syst_M[syst][data_op][list_vv[0]]
            hist_syst_E[syst][data_op]["vv"] = hist_syst_E[syst][data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      if False:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["vv"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["vv"].Add(hist_syst_E[syst][data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]
  if args.wcs: samples = ["ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom_M = {}
  histT_nom_E = {}
  histT_syst_M = {}
  histT_syst_E = {}
  for syst in list_syst:
      histT_syst_M[syst] = {}
      histT_syst_E[syst] = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       if False:
          for syst in list_syst:
              histT_syst_M[syst][s] =  hist_syst_M[syst][datayears[0]][s]
              histT_syst_E[syst][s] =  hist_syst_E[syst][datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          if False:
             for syst in list_syst:
                 histT_syst_M[syst][s].Add(hist_syst_M[syst][d][s])
                 histT_syst_E[syst][s].Add(hist_syst_E[syst][d][s])
       if (channel in sl_channel) and (name not in not_rebin):
          histT_nom_M[s].Rebin(nrebin)
          histT_nom_E[s].Rebin(nrebin)
          if False:
             for syst in list_syst:
                 histT_syst_M[syst][s].Rebin(nrebin)
                 histT_syst_E[syst][s].Rebin(nrebin)

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])
    if (channel in sl_channel) and (name not in not_rebin): 
       histD_M.Rebin(nrebin)
       histD_E.Rebin(nrebin)
    if ("_pt" in name or ("transverse" in name)) and name != "muon_jet_pt":
       histD_M.GetXaxis().SetRange(1,histD_M.GetNbinsX() + 1)
       histD_E.GetXaxis().SetRange(1,histD_E.GetNbinsX() + 1)

  if False:
    histT_sT_syst_M = {}
    histT_sT_syst_E = {}
    ### Para los histogramas de sistematica up y down sumamos todas las contribuciones
    for syst in list_syst:
        histT_sT_syst_M[syst] = histT_syst_M[syst][samples[0]]
        histT_sT_syst_E[syst] = histT_syst_E[syst][samples[0]]
        for s in samples[1:]:
            histT_sT_syst_M[syst].Add(histT_syst_M[syst][s])
            histT_sT_syst_E[syst].Add(histT_syst_E[syst][s])

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ttbar_dl"] = (222,90,106)
  colors["wjets"] = (155,152,204)
  colors["ttbar_sl_bottomgluon"] = (236,164,207)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_light"] = (120,154,86)
  colors["ttbar_sl_nocharm"] = (120,154,86)
  colors["ttbar_sl_bottom"] = (201,79,152)
  colors["ttbar_sl_else"] = (242,193,121)
  colors["ttbar_sl_charmgluon"] = (222,212,74)
  colors["vv"] = (255,180,85)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (113,209,223)
  colors["st"] = (153,51,255)
  colors["st_charm"] = (102,0,204)
  colors["st_nocharm"] = (198,101,222)
  colors["st_else"] = (207,176,235)

  ##### para los daltonicos de los cojones
  colors["ttbar_sl_charm"] = (152,228,230)
  colors["ttbar_sl_nocharm"] = (255,178,102)
  colors["st_charm"] = (131,45,182)
  colors["st_nocharm"] = (192,192,192)
  colors["st_else"] = (166,114,101)
  colors["ttbar_dl"] = (243,89,0)
  colors["zjets"] = (198,201,114)
  colors["wjets"] = (108,121,152)
  colors["vv"] = (0,128,255)
  colors["ttbar_dh"] = (35,192,30)
  colors["vjets_plus"] = (0,128,255)

  #colors["ttbar_sl_charm"] = (204,255,153)
  #colors["ttbar_sl_nocharm"] = (120,154,86)
  #colors["st_charm"] = (102,0,204)
  #colors["st_nocharm"] = (198,101,222)
  #colors["st_else"] = (207,176,235)
  #colors["ttbar_dl"] = (222,90,106)
  #colors["zjets"] = (113,209,223)
  #colors["wjets"] = (155,152,204)
  #colors["vv"] = (255,180,85)
  #colors["vjets_plus"] = (155,152,204)

  ############# Axis Y name deff
  a = histT_nom_M[samples[0]].GetXaxis().GetXmin()
  b = histT_nom_M[samples[0]].GetXaxis().GetXmax()
  c = histT_nom_M[samples[0]].GetNbinsX()
  aux_n = round((b-a)/c,2)
  if (("pt" in name or "pT" in name) or ("Inv" in name or "mass" in name) or "bs" in name) and "eta" not in name:
     axisY_label_name[name] = "Events /" + str(aux_n)+" GeV"
  else:
     axisY_label_name[name] = "Events /" + str(aux_n)

  if args.stack:
    ymax_M = 0
    ymin_M = 0
    ymax_E = 0
    ymin_E = 0
    for s in samples:
      histT_nom_M[s].SetLineWidth(1)
      histT_nom_M[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_M[s].GetYaxis().SetMaxDigits(3)
      histT_nom_M[s].GetYaxis().SetTitle(axisY_label_name[name])
      histT_nom_M[s].GetXaxis().SetTitle(name)
      histT_nom_M[s].GetYaxis().SetTitleOffset(titY_off2)
      histT_nom_M[s].GetYaxis().SetTitleSize(titY_size2)
      if ("_pt" in name or ("transverse" in name))  and name != "muon_jet_pt": histT_nom_M[s].GetXaxis().SetRange(1, histT_nom_M[s].GetNbinsX() + 1);
      histT_nom_E[s].SetLineWidth(1)
      histT_nom_E[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_E[s].GetYaxis().SetMaxDigits(3)
      histT_nom_E[s].GetYaxis().SetTitle(axisY_label_name[name])
      histT_nom_E[s].GetXaxis().SetTitle(name)
      histT_nom_E[s].GetYaxis().SetTitleOffset(titY_off2)
      histT_nom_E[s].GetYaxis().SetTitleSize(titY_size2)
      if ("_pt" in name or ("transverse" in name))  and name != "muon_jet_pt": histT_nom_E[s].GetXaxis().SetRange(1, histT_nom_E[s].GetNbinsX() + 1);

      y = histT_nom_M[s].GetMaximum()
      ym = histT_nom_M[s].GetMinimum()
      if y>ymax_M: ymax_M=y
      if ymin_M>ym: ymin_M=ym
      if not args.nodata:
        y = histD_M.GetMaximum()
        if y>ymax_M: ymax_M=y

      y = histT_nom_E[s].GetMaximum()
      ym = histT_nom_E[s].GetMinimum()
      if y>ymax_E: ymax_E=y
      if ymin_E>ym: ymin_E=ym
      if not args.nodata:
        y = histD_E.GetMaximum()
        if y>ymax_E: ymax_E=y

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if not args.linear: 
         histT_nom_M["vv"].SetMaximum(3*ymax_M)
      else:
         histT_nom_M["vv"].SetMaximum(1.3*ymax_M)
         histT_nom_M["vv"].SetMinimum(1.3*ymin_M)

      #histT_E["vv"].SetMinimum(1000.)
      #histT_E["vv"].SetMaximum(3*ymax_E)
      if not args.linear: 
         histT_nom_E["vv"].SetMaximum(3*ymax_M)
      else:
         histT_nom_E["vv"].SetMaximum(1.3*ymax_E)
         histT_nom_E["vv"].SetMinimum(1.3*ymin_E)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
    if args.wcs: samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

    ############### Backgrounds reorganisation
    if args.wcs: 
       histT_nom_M["vjets_plus"] = histT_nom_M["vv"]
       histT_nom_E["vjets_plus"] = histT_nom_E["vv"]
       for s in ["ttbar_dh","zjets","wjets","st_else","st_nocharm"]:
           histT_nom_M["vjets_plus"].Add(histT_nom_M[s])
           histT_nom_E["vjets_plus"].Add(histT_nom_E[s])
       samples = ["vjets_plus","ttbar_dl","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
    ###############################################

    if args.ratio: upper_pad.cd()

    if args.sumEM:
       histT_nom_T = {}
       for s in samples:
         histT_nom_T[s] = histT_nom_M[s]
         histT_nom_T[s].Add(histT_nom_E[s])
       #stack_T = ROOT.THStack("hs", ";;Events")
       stack_T = ROOT.THStack()
       if args.postfit:
          for s in ["st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]:
             histT_nom_T[s].Scale(aux_rwqq)
          histT_nom_T["ttbar_sl_charm"].Scale(aux_ratio)
          histT_nom_T["ttbar_sl_nocharm"].Scale(aux_conratio)
       for s in samples:
         if ("jet" in name) and ("pt" in name) and (not "muon_jet" in name): 
            if (not "long" in name):
               if (channel in sl_channel) and (name not in not_rebin):
                  histT_nom_T[s].SetBinContent(1,0.)
               else:
                  histT_nom_T[s].SetBinContent(5,0.)
            else:
               if (channel in sl_channel) and (name not in not_rebin):
                  histT_nom_T[s].SetBinContent(1,0.)
               else:
                  histT_nom_T[s].SetBinContent(2,0.)
         if "pt_long" in name:
            if (channel in sl_channel) and (name not in not_rebin):
               histT_nom_T[s].SetBinContent(1,0.)
            else:
               histT_nom_T[s].SetBinContent(2,0.)
         stack_T.Add(histT_nom_T[s])
       last_T = stack_T.GetStack().Last()
       if ("_pt" in name or ("transverse" in name))  and name != "muon_jet_pt": last_T.GetXaxis().SetRange(1,last_T.GetNbinsX()+1);

    else:
       stack_M = ROOT.THStack()
       stack_E = ROOT.THStack()
       if args.postfit:
          for s in ["st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]:
             histT_nom_M[s].Scale(aux_rwqq)
             histT_nom_E[s].Scale(aux_rwqq)
          histT_nom_M["ttbar_sl_charm"].Scale(aux_ratio)
          histT_nom_E["ttbar_sl_nocharm"].Scale(aux_conratio)
          histT_nom_M["ttbar_sl_charm"].Scale(aux_ratio)
          histT_nom_E["ttbar_sl_nocharm"].Scale(aux_conratio)
       for s in samples:
         if ("jet" in name) and ("pt" in name) and (not "muon_jet" in name):
            if (not "long" in name):
               if (channel in sl_channel) and (name not in not_rebin):
                  histT_nom_M[s].SetBinContent(1,0.)
                  histT_nom_E[s].SetBinContent(1,0.)
       	       else:
                  histT_nom_M[s].SetBinContent(5,0.)
                  histT_nom_E[s].SetBinContent(5,0.)
            else:
               if (channel in sl_channel) and (name not in not_rebin):
                  histT_nom_M[s].SetBinContent(1,0.)
                  histT_nom_E[s].SetBinContent(1,0.)
               else:
                  histT_nom_M[s].SetBinContent(2,0.)
                  histT_nom_E[s].SetBinContent(2,0.)
         stack_M.Add(histT_nom_M[s])
         stack_E.Add(histT_nom_E[s])

       last_M = stack_M.GetStack().Last()
       last_E = stack_E.GetStack().Last()
       if ("_pt" in name or ("transverse" in name))  and name != "muon_jet_pt": 
              last_M.GetXaxis().SetRange(1,last_M.GetNbinsX()+1);
              last_E.GetXaxis().SetRange(1,last_E.GetNbinsX()+1);

    if args.sumEM:
      if not args.nosyst:
        graph_err_T = TGraphAsymmErrors();
        ratio_graph_err_T = TGraphAsymmErrors();

        n_aux = last_T.GetNbinsX()+1
        for bin in range(n_aux):
          graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),last_T.GetBinContent(bin+1));
          graph_err_T.SetPointEYhigh(bin,last_T.GetBinContent(bin+1)*syst_val[channel]);
          graph_err_T.SetPointEYlow(bin,last_T.GetBinContent(bin+1)*syst_val[channel]);
 
          ratio_graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),1.);
          if (last_T.GetBinContent(bin+1)>0.001):
             ratio_graph_err_T.SetPointEYhigh(bin,syst_val[channel]);
             ratio_graph_err_T.SetPointEYlow(bin,syst_val[channel]);
          else:
             ratio_graph_err_T.SetPointEYhigh(bin,0.);
             ratio_graph_err_T.SetPointEYlow(bin,0.);

        graph_err_T.SetFillColorAlpha(14,1);
        graph_err_T.SetMarkerColor(kGray);
        graph_err_T.SetMarkerSize(0);
        graph_err_T.SetFillStyle(3002);
        ratio_graph_err_T.SetFillColorAlpha(15,1);
        ratio_graph_err_T.SetMarkerColor(kGray);
        ratio_graph_err_T.SetMarkerSize(0);
        ratio_graph_err_T.SetFillStyle(3001);
    else:
      if not args.nosyst:
        graph_err_M = TGraphAsymmErrors();
        graph_err_E = TGraphAsymmErrors();
        ratio_graph_err_M = TGraphAsymmErrors();
        ratio_graph_err_E = TGraphAsymmErrors();

        n_aux = last_M.GetNbinsX()+1
        for bin in range(n_aux):
          graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),last_M.GetBinContent(bin+1));
          graph_err_M.SetPointEYhigh(bin,last_M.GetBinContent(bin+1)*syst_val[channel]);
          graph_err_M.SetPointEYlow(bin,last_M.GetBinContent(bin+1)*syst_val[channel]);

          graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),last_E.GetBinContent(bin+1));
          graph_err_E.SetPointEYhigh(bin,last_E.GetBinContent(bin+1)*syst_val[channel]);
          graph_err_E.SetPointEYlow(bin,last_E.GetBinContent(bin+1)*syst_val[channel]);

          ratio_graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),1.);
          if (last_M.GetBinContent(bin+1)>0.01):
             ratio_graph_err_M.SetPointEYhigh(bin,syst_val[channel]);
             ratio_graph_err_M.SetPointEYlow(bin,syst_val[channel]);
          else:
             ratio_graph_err_M.SetPointEYhigh(bin,0.);
             ratio_graph_err_M.SetPointEYlow(bin,0.);

          ratio_graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),1.);
          if (last_E.GetBinContent(bin+1)>0.01):
             ratio_graph_err_E.SetPointEYhigh(bin,syst_val[channel]);
             ratio_graph_err_E.SetPointEYlow(bin,syst_val[channel]);
          else:
             ratio_graph_err_E.SetPointEYhigh(bin,0.);
             ratio_graph_err_E.SetPointEYlow(bin,0.);

        graph_err_M.SetFillColorAlpha(14,1);
        graph_err_M.SetMarkerColor(kGray);
        graph_err_M.SetMarkerSize(0);
        graph_err_M.SetFillStyle(3002);
        graph_err_E.SetFillColorAlpha(14,1);
        graph_err_E.SetMarkerColor(kGray);
        graph_err_E.SetMarkerSize(0);
        graph_err_E.SetFillStyle(3002);
        ratio_graph_err_M.SetFillColorAlpha(15,1);
        ratio_graph_err_M.SetMarkerColor(kGray);
        ratio_graph_err_M.SetMarkerSize(0);
        ratio_graph_err_M.SetFillStyle(3001);
        ratio_graph_err_E.SetFillColorAlpha(15,1);
        ratio_graph_err_E.SetMarkerColor(kGray);
        ratio_graph_err_E.SetMarkerSize(0);
        ratio_graph_err_E.SetFillStyle(3001);

    if args.sumEM:
      #####################################
      ########## sumEM channel ############
      #####################################
      ymax_T = ymax_E+ymax_M
      y_T = stack_T.GetMaximum()
      if y_T>ymax_T: ymax_T=y_T
      ymin_T = stack_T.GetMinimum()
      #stack_M.SetMinimum(1000.)
      #stack_M.SetMaximum(3*ymax_T)
      if not args.linear: 
         stack_T.SetMaximum(3*ymax_T)
      else:
         stack_T.SetMaximum(1.25*ymax_T)
         stack_T.SetMinimum(1.)
         
      if not args.nodata:
        histD_T = histD_M
        histD_T.Add(histD_E)
        if ("_pt" in name or ("transverse" in name)) and name != "muon_jet_pt": histD_T.GetXaxis().SetRange(1,histD_T.GetNbinsX() + 1)

      stack_T.Draw("HIST")
      stack_T.GetYaxis().SetTitle(axisY_label_name[name])
      if not args.ratio: 
         stack_T.GetXaxis().SetTitle(cosmetic_names[name])
         stack_T.GetYaxis().SetLabelSize(titY_size2)
      if ("_pt" in name or ("transverse" in name))   and name != "muon_jet_pt": stack_T.GetXaxis().SetRange(1,last_T.GetNbinsX() + 1)
      stack_T.GetYaxis().SetTitleSize(titY_size2)
      stack_T.GetYaxis().SetMaxDigits(3)
      stack_T.GetXaxis().SetLabelSize(titX_lab)
      stack_T.GetYaxis().SetLabelSize(titY_lab)
      stack_T.GetYaxis().SetLabelFont(42)
      stack_T.GetYaxis().SetTitleFont(42)
      if not args.ratio:
         stack_T.GetXaxis().SetLabelFont(42)
         stack_T.GetXaxis().SetTitleFont(42)
      stack_T.GetYaxis().SetTitleOffset(titY_off2)
      if not args.nosyst: 
         graph_err_T.Draw("SAME 2")
         if ("_pt" in name or ("transverse" in name))   and name != "muon_jet_pt": graph_err_T.GetXaxis().SetRange(1,last_T.GetNbinsX()+1) 
      if not args.nodata:
        histD_T.SetMarkerStyle(20)
        if channel == "sl_ssos":
           histD_T.SetMarkerSize(1.5)
        else:
           histD_T.SetMarkerSize(1)
        histD_T.SetLineWidth(1)
        histD_T.SetLineColor(ROOT.kBlack)
        histD_T.Draw("EX0 SAME")

      if args.ratio and not args.nodata:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        #frame.Draw()

      if args.ratio:
        lower_pad.cd()
        lower_pad.SetGridy();
        if not args.nodata:
          ratio = histD_T.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(20)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/Pred")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(titX_lab)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetLabelFont(42)
          ratio.GetYaxis().SetLabelFont(42)
          ratio.GetYaxis().SetNdivisions(1008)
          ratio.GetXaxis().SetTitleFont(42)
          ratio.GetYaxis().SetTitleFont(42)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(titY2_lab)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle(False)
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          if args.wcs:
             hTotal = histT_nom_T["vjets_plus"].Clone('hTotal')
          else:
             hTotal = histT_nom_T["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_T[s])
          ratio.Divide(hTotal)
          ratio.Draw("epx0")
          if not args.nosyst: ratio_graph_err_T.Draw("same 2")
          ratio.Draw("epx0 same")
        else:
          if not args.nosyst:
             ratio = histT_nom_T["vjets_plus"].Clone('ratio')
             ratio.SetLineColor(kBlack)
             ratio.SetMarkerStyle(6)
             ratio.SetTitle("")
             ratio.SetMinimum(c_rat2)
             ratio.SetMaximum(c_rat)
             ratio.GetYaxis().SetTitle("Ctag syst ratio")
             ratio.GetXaxis().SetTitle(cosmetic_names[name])
             ratio.GetXaxis().SetLabelSize(0.08)
             ratio.GetXaxis().SetTitleSize(tit_size)
             ratio.GetXaxis().SetTitleOffset(titX_off)
             ratio.GetYaxis().SetLabelSize(0.05)
             ratio.GetYaxis().SetTitleSize(titY_size)
             ratio.GetYaxis().CenterTitle(False)
             ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
             ratio.GetYaxis().SetTitleOffset(titY_off)
             # Set up plot for markers and errors
             ratio.Sumw2()
             ratio.SetStats(0)
             ratio.Divide(ratio)
             ratio.Draw("hist p")
             ratio_graph_err_T.Draw("same 2")

      if "flavourP" in str(name) and args.wcs:
         print(name)
         for s in ["ttbar_dl","ttbar_sl_nocharm","ttbar_sl_charm"]:
           print(s)
           qaux1 = histT_nom_T[s].Integral()
           print("Total events are "+str(qaux1))
           for i in [2,3,11,12,28]:
              print("Content of bin "+str(i)+" is "+str(histT_nom_T[s].GetBinContent(i)))
              print("Percentage over total is "+str(histT_nom_T[s].GetBinContent(i)/qaux1))

      if name == "jet_max_cvltag" and (not args.nodata):
          print("charm tag")
          print("Integral of data is "+str(histD_T.Integral()))
          print("Integral of MC is "+str(hTotal.Integral()))
          for s in samples:
              print("Integral of "+str(s)+" MC is "+str(histT_nom_T[s].Integral()))
          print("Integral of data from 0.3 is "+str(histD_T.Integral(16,51)))
          print("Integral of MC from 0.3 is "+str(hTotal.Integral(16,51)))
          for s in samples:
              print("Integral of "+str(s)+" MC from 0.3 is "+str(histT_nom_T[s].Integral(16,51)))

      if (name == "InvM_2jets" or name == "InvM3_good") and (not args.nodata):
        print("Integral of data is "+str(histD_T.Integral()))
        print("Integral of MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_T.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_T[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_T[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_T[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_T["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["st_charm"],"Single top cq","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["vjets_plus"],"V+jets and other","f")
        leg.AddEntry(histT_nom_T["ttbar_dl"],"t#bar{t} dileptonic","f")
        #leg.AddEntry(histT_nom_T["st_nocharm"],"Single top tW uq","f")
        #leg.AddEntry(histT_nom_T["st_else"],"Single top s/t channel","f")
        #leg.AddEntry(histT_nom_T["zjets"],"Z+jets","f")
        #leg.AddEntry(histT_nom_T["wjets"],"W+jets","f")
        if not args.nosyst: leg.AddEntry(graph_err_T, "Syst. Unc." ,"f")
        if args.stack and not args.nodata: leg.AddEntry(histD_T, "Data" ,"ep")
      else:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["ttbar_dl"],"Dileptonic t#bar{t}","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["st"],"Single top","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_T["zjets"],"Z+jets","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_T["wjets"],"W+jets","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        if args.stack and not args.nodata: leg.AddEntry(histD_T, "Data" ,"lep")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      #leg.AddEntry(histT_nom_T["ttbar_dh"],"Hadronic t#bar{t}","f")
      leg.SetTextFont(42);
      leg.Draw()
      termp= "totalHT_wqq"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: 
         c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    else:
      #####################################
      ############ M channel ##############
      #####################################
      y_M = stack_M.GetMaximum()
      if y_M>ymax_M: ymax_M=y_M
      #stack_M.SetMinimum(1000.)
      #stack_M.SetMaximum(3*ymax_M)
      if args.linear: stack_M.SetMaximum(1.3*ymax_M)
      if args.linear: stack_M.SetMinimum(1.)
      y_E = stack_E.GetMaximum()
      if y_E>ymax_E: ymax_E=y_E
      #stack_E.SetMinimum(1000.)
      #stack_E.SetMaximum(3*ymax_E)
      if args.linear: stack_E.SetMaximum(1.3*ymax_E)
      if args.linear: stack_M.SetMinimum(1.)

      stack_M.Draw("HIST")
      stack_M.GetYaxis().SetTitle(axisY_label_name[name])
      stack_M.GetYaxis().SetMaxDigits(3)
      stack_M.GetYaxis().SetTitleSize(titY_size2)
      stack_M.GetYaxis().SetTitleOffset(titY_off2)
      stack_M.GetXaxis().SetLabelSize(titX_lab)
      stack_M.GetYaxis().SetLabelSize(titY_lab)
      stack_M.GetYaxis().SetLabelFont(42)
      stack_M.GetYaxis().SetTitleFont(42)
      if not args.ratio:
         stack_M.GetXaxis().SetLabelFont(42)
         stack_M.GetXaxis().SetTitleFont(42)
      if ("_pt" in name or ("transverse" in name))   and name != "muon_jet_pt": stack_M.GetXaxis().SetRange(1,last_M.GetNbinsX() + 1)
      if not args.ratio:
         stack_M.GetXaxis().SetTitle(cosmetic_names[name])
         stack_M.GetYaxis().SetLabelSize(titY_size2)
      if not args.nosyst: 
         graph_err_M.Draw("SAME 2")
         if ("_pt" in name or ("transverse" in name))  and name != "muon_jet_pt": graph_err_M.GetXaxis().SetRange(1,last_M.GetNbinsX()+1) 
      if not args.nodata:
        histD_M.SetMarkerStyle(20)
        if channel == "sl_ssos":
           histD_M.SetMarkerSize(1.5)
        else:
           histD_M.SetMarkerSize(1)
        histD_M.SetLineWidth(1)
        histD_M.SetLineColor(ROOT.kBlack)
        histD_M.Draw("EX0 SAME")

      if args.ratio:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        lower_pad.cd()
        lower_pad.SetGridy();
        if not args.nodata:
          ratio = histD_M.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(20)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/Pred")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(titX_lab)
          ratio.GetXaxis().SetLabelFont(42)
          ratio.GetYaxis().SetLabelFont(42)
          ratio.GetYaxis().SetNdivisions(1008)
          ratio.GetXaxis().SetTitleFont(42)
          ratio.GetYaxis().SetTitleFont(42)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(titY2_lab)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle()
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          if args.wcs:
             hTotal = histT_nom_M["vjets_plus"].Clone('hTotal')
          else:
             hTotal = histT_nom_M["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_M[s])
          ratio.Divide(hTotal)
          ratio.Draw("epx0")
          if not args.nosyst: ratio_graph_err_M.Draw("same 2")
          ratio.Draw("epx0 same")

      if (name == "InvM_2jets" or name == "InvM3_good") and (not args.nodata):
        print("Integral of M data is "+str(histD_M.Integral()))
        print("Integral of M MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_M.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_M[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_M[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_M[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_M["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["st_charm"],"Single top cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_M["vjets_plus"],"V+jets and other","f")
        leg.AddEntry(histT_nom_M["ttbar_dl"],"t#bar{t} dileptonic","f")
        #leg.AddEntry(histT_nom_M["st_nocharm"],"Single top tW uq","f")
        #leg.AddEntry(histT_nom_M["st_else"],"Single top s/t channel","f")
        #leg.AddEntry(histT_nom_M["zjets"],"Z+jets","f")
        #leg.AddEntry(histT_nom_M["wjets"],"W+jets","f")
        if not args.nosyst: leg.AddEntry(graph_err_M, "Syst. Unc." ,"f")
        if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"ep")
      else:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_M["st"],"Single top","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
        leg.AddEntry(histT_nom_M["zjets"],"Z+jets","f")
        leg.AddEntry(histT_nom_M["wjets"],"W+jets","f")
        leg.AddEntry(histT_nom_M["ttbar_dl"],"Dileptonic t#bar{t}","f")
        #leg.AddEntry(histT_nom_M["ttbar_dh"],"Hadronic t#bar{t}","f")
        if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"lep")
      leg.SetTextFont(42);
      leg.Draw()
      termp= "totalHT_wqq_M"
      if args.ratio: 
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else: 
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      #####################################
      ############ E channel ##############
      #####################################
      stack_E.Draw("HIST")
      if not args.ratio:
         stack_E.GetXaxis().SetTitle(cosmetic_names[name])
         stack_E.GetYaxis().SetLabelSize(titY_size2)
      stack_E.GetYaxis().SetTitle(axisY_label_name[name])
      stack_E.GetYaxis().SetMaxDigits(3)
      stack_E.GetYaxis().SetTitleSize(titY_size2)
      stack_E.GetYaxis().SetTitleOffset(titY_off2)
      stack_E.GetXaxis().SetLabelSize(titX_lab)
      stack_E.GetYaxis().SetLabelSize(titY_lab)
      stack_E.GetYaxis().SetLabelFont(42)
      stack_E.GetYaxis().SetTitleFont(42)
      if not args.ratio:
         stack_E.GetXaxis().SetLabelFont(42)
         stack_E.GetXaxis().SetTitleFont(42)
      if ("_pt" in name or ("transverse" in name))   and name != "muon_jet_pt": stack_E.GetXaxis().SetRange(1,last_E.GetNbinsX() + 1)
      if not args.nosyst: 
         graph_err_E.Draw("SAME 2")
         if ("_pt" in name or ("transverse" in name))   and name != "muon_jet_pt": graph_err_E.GetXaxis().SetRange(1,last_E.GetNbinsX()+1) 
      if not args.nodata:
        histD_E.SetMarkerStyle(20)
        if channel == "sl_ssos":
           histD_E.SetMarkerSize(1.5)
        else:
           histD_E.SetMarkerSize(1)
        histD_E.SetLineWidth(1)
        histD_E.SetLineColor(ROOT.kBlack)
        histD_E.Draw("EX0 SAME")

      if args.ratio and not args.nodata:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        lower_pad.cd()
        lower_pad.SetGridy();
        if not args.nodata:
          ratio = histD_E.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(20)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/Pred")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(titX_lab)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetLabelFont(42)
          ratio.GetYaxis().SetLabelFont(42)
          ratio.GetYaxis().SetNdivisions(1008)
          ratio.GetXaxis().SetTitleFont(42)
          ratio.GetYaxis().SetTitleFont(42)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(titY2_lab)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle()
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          if args.wcs:
             hTotal = histT_nom_E["vjets_plus"].Clone('hTotal')
          else:
             hTotal = histT_nom_E["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_E[s])
          ratio.Divide(hTotal)
          ratio.Draw("epx0")
          if not args.nosyst: ratio_graph_err_E.Draw("same 2")
          ratio.Draw("epx0 same")

      if (name == "InvM_2jets" or name == "InvM3_good") and (not args.nodata):
        print("Integral of E data is "+str(histD_E.Integral()))
        print("Integral of E MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_E.Integral()/hTotal.Integral()))
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_E[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_E[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_E[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_E["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["st_charm"],"Single top cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_E["vjets_plus"],"V+jets and other","f")
        leg.AddEntry(histT_nom_E["ttbar_dl"],"t#bar{t} dileptonic","f")
        #leg.AddEntry(histT_nom_E["st_nocharm"],"Single top tW uq","f")
        #leg.AddEntry(histT_nom_E["st_else"],"Single top s/t channel","f")
        #leg.AddEntry(histT_nom_E["zjets"],"Z+jets","f")
        #leg.AddEntry(histT_nom_E["wjets"],"W+jets","f")
        if not args.nosyst: leg.AddEntry(graph_err_E, "Syst. Unc." ,"f")
        if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"ep")
      else:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_E["st"],"Single top","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
        leg.AddEntry(histT_nom_E["zjets"],"Z+jets","f")
        leg.AddEntry(histT_nom_E["wjets"],"W+jets","f")
        leg.AddEntry(histT_nom_E["ttbar_dl"],"Dileptonic t#bar{t}","f")
        if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"lep")
      leg.SetTextFont(42);
      leg.Draw()
      termp= "totalHT_wqq_E"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    for data_op in datayears:
                        histFile[name][data_op].Close()
                        if not args.nodata:
                          histFileDM[name][data_op].Close()
                          histFileDE[name][data_op].Close()


