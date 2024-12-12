######################################                                        
######################################                                        
#######      GOOD  VERSION     #######
######################################                                     
######################################                                     

print('SL CHANNEL')

## Modified version, meant to be run locally, not sent to condor 

## Selection for noth MC and data samples, created to distinguish between years
## Different histogram files are produced for each situation 
## It includes an option for SSOS substraction

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir

import json
import argparse

sys.path.append('/nfs/cms/vazqueze/ttbaranalisis/last_corrections/')

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

ROOT.EnableImplicitMT()

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--process", type=string, default="WW",
                    help="Select type of process to run")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--type", type=string, default="data",
                    help="Selec type of data to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if args.process == "allMCH": proc = ["higgs_tocs_m80","higgs_tocs_m90","higgs_tocs_m100","higgs_tocs_m120","higgs_tocs_m140","higgs_tocs_m150",
           "higgs_tocs_m155","higgs_tocs_m160"]
elif args.process == "allMC": proc = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zz","wz","st_1","st_2","st_3","st_4"]
elif (args.process == "higgs_tocs_m80" or args.process == "higgs_tocs_m90" or args.process == "higgs_tocs_m100" or args.process == "higgs_tocs_m120"
        or args.process == "higgs_tocs_m140" or args.process == "higgs_tocs_m150"
        or args.process == "higgs_tocs_m155" or args.process == "higgs_tocs_m160"
        or args.process == "ww" or args.process == "wjets_1" or args.process == "wjets_2" or args.process == "wjets_3" or args.process == "wjets_4" 
        or args.process == "ttbar_sl" or args.process == "wjets_5" or args.process == "wjets_6" or args.process == "wjets_7" or args.process == "wjets_8"
        or args.process == "zjets_1" or args.process == "zjets_2" or args.process == "zjets_3" or args.process == "zjets_4"
       	or args.process == "zjets_5" or args.process == "zjets_6" or args.process == "zjets_7" or args.process == "zjets_8"
        or args.process == "DY" or args.process == "WZ"or args.process == "ZZ" or args.process == "ST1"
        or args.process == "ST2" or args.process == "ST3" or args.process == "ST4"  or args.process == "DYjets"  or args.process == "ttbar_dl"
        or args.process == "ttbar_dh" or args.process == "M" or args.process == "E"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if args.year == "all": years = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = [str(args.year)]
else: raise NameError('Incorrect year')

if args.type == "data": mode = "data"
elif args.type == "mc": mode = "mc"
else: raise NameError('Incorrect type')

samples = []

for p in proc:
  for y in years:
    if mode == "mc": samples.append(p+y)
    else: samples.append(y+p)

samples_test = samples[:]

print("the samples treated are",samples)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

if mode == "mc": files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+years[0]+".json"))
else: files = json.load(open("/nfs/cms/vazqueze/higgssearch/datainfo.json"))

processes = files.keys()

df = {}
xsecs = {}
sumws = {}
archives = {}
event_test = {}

for s in samples:
  archives[s]=[]

names_proc = {"M":{"2016":["2016M1","2016M2","2016M3","2016M4","2016M5","2016M6"],"2016B":["2016BM1","2016BM2","2016BM3"],
   "2017":["2017M1","2017M2","2017M3","2017M4","2017M5"],"2018":["2018M1","2018M2","2018M3","2018M4"]},
   "E":{"2016":["2016E1","2016E2","2016E3","2016E4","2016E5","2016E6"],"2016B":["2016BE1","2016BE2","2016BE3"],
   "2017":["2017E1","2017E2","2017E3","2017E4","2017E5"],"2018":["2018E1","2018E2","2018E3","2018E4"]}}

files_data = {"M":{"2016":230,"2016B":157,"2017":436,"2018":393},
   "E":{"2016":335,"2016B":156,"2017":245,"2018":738}}

term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/"

for p in proc:
    # Construct the dataframes
    folder = term1+"myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if args.type == "mc":
        folder = folder.format(sample = p) # Sample name
        num_files = files[p]["files"] # Number of files
        list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
        if (num_files == len(list_files)):
           for f in list_files:
               file_r = TFile(join(folder, f))
               if file_r.GetListOfKeys().Contains("Events"):
                  archives[p+years[0]].append(join(folder,f))
    else:
        if p == "M":
          folder1 = term1+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_mu")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["M"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     archives[years[0]+"M"].append(join(folder,fil))
        if p == "E":
          folder1 = term1+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_el")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["E"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     archives[years[0]+"E"].append(join(folder,fil))

for s in samples:
  if args.notfull: archives[s]=archives[s][int(args.list[0]):int(args.list[1])]
  df[s] = ROOT.RDataFrame("Events",set(archives[s]))
  print("Number of files for",s,len(archives[s]))

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

cuts_btag["2016"]=[0.0614, 0.3093, 0.7221]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=26
muon_pt["2016B"]=26
muon_pt["2017"]=29
muon_pt["2018"]=26

el_pt = {}

el_pt["2016"]=30
el_pt["2016B"]=30
el_pt["2017"]=35
el_pt["2018"]=35

## Triggers per year

muon_trig = {}

muon_trig["2016"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2016B"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2017"]="HLT_IsoMu27"
muon_trig["2018"]="HLT_IsoMu24"

el_trig = {}

el_trig["2016"]="HLT_Ele27_WPTight_Gsf"
el_trig["2016B"]="HLT_Ele27_WPTight_Gsf"
el_trig["2017"]="HLT_Ele32_WPTight_Gsf"
el_trig["2018"]="HLT_Ele32_WPTight_Gsf"

### MET Filters per year

met_filter = {}

met_filter["2016"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2016B"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2017"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_BadPFMuonFilter")
met_filter["2018"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")

###################################################
################   DEFINITIONS   ##################
###################################################

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi, Vint puID, Vint jetID) {
            vector<int> vb;
	    bool cond = false;
	    bool cond1 = false;
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>25. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>25. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };

      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
	    auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
	    auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
	    auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
      auto InvariantM3(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1, const float pt2, const float eta2, const float phi2, const float mass2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);
            auto e2 = std::sqrt(x2*x2+y2*y2+z2*z2+mass2*mass2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

######### bottom jets ordering

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      template <typename T>
      vector<size_t> sort_indexes(const vector<T> &v) {

          // initialize original index locations
          vector<size_t> idx(v.size());
          iota(idx.begin(), idx.end(), 0);

          // sort indexes based on comparing values in v
          // using std::stable_sort instead of std::sort
          // to avoid unnecessary index re-orderings
          // when v contains elements of equal values 
          stable_sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

          return idx;
      };
      auto bottomjets(UInt_t njetgood, Vint jetgood, UInt_t njet, Vfloat jet_btag) {
          vector<float> vb;
          vector<int> fin;
          if (njetgood>0) {
                for (unsigned int i=0; i<njetgood; ++i){
                    vb.push_back(jet_btag[jetgood[i]]);
                }
          }
          for (auto i: sort_indexes(vb)){
                fin.push_back(i);
          }
          return fin;
      };
""")

## Eleccion de los jets que no sean bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInds(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, Vint jetbotind, Vfloat qgl) {
            vector<int> vb;
            bool condb = true;
            int ind1 = -1;
            int ind2 = -1;
            float ptJ{-10.};
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        //if(pt[good[j]]>ptJ && condb && qgl[good[j]]>0.25){
                        if(pt[good[j]]>ptJ && condb){
                                ind1 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind1 > -1) {
                vb.push_back(ind1);
            }
            ptJ = -10.;
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        //if(pt[good[j]]>ptJ && condb && good[j] != ind1 && qgl[good[j]]>0.25){
                        if(pt[good[j]]>ptJ && condb && good[j] != ind1){
                                ind2 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind2 > -1) {
                vb.push_back(ind2);
            }
            return vb;
      };
""")

######### Masa invariante con los bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMassBot(Vint jetbot, Vint jetmuon, Vint jetgood, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, UInt_t jet_notmuon, Vfloat Jet_mass) {
            int ind = -1;
            vector<float> vb;
            float dR1 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));           
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }
            return vb;
      };
""")

## pT component calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto ptot = plep1+plep2+phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py(); auto ptotX = ptot.Px(); auto ptotY = ptot.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = std::sqrt(ptotX*ptotX + ptotY*ptotY);
            return ptotmod/(plepmod+phadmod);
      };
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector plep1;
            float deltaRlepjet1;
       	    float deltaRlepjet2;
       	    float deltaPhilepjet1;
            float deltaPhilepjet2;
            float deltaEtalepjet1;
            float deltaEtalepjet2;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[0]],mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
       	       	  deltaRlepjet2	= ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
       	       	  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
       	       	  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
       	       	  deltaRlepjet1	= ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
       	       	  deltaPhilepjet1 = fabs(el_phi[el_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(el_phi[el_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(fabs(plep2.Phi()-phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(fabs(plep.Phi()-phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(fabs(plep1.Phi()-phad.Phi()));
            vb.push_back(fabs(plep.Eta()-phad.Eta()));
            vb.push_back(fabs(plep1.Eta()-phad.Eta()));
            vb.push_back(fabs(plep.Pt()-phad.Pt()));
            vb.push_back(fabs(plep1.Pt()-phad.Pt()));
            vb.push_back(phad2.Pt()*std::sin(phad1.Phi()-phad2.Phi()));
            vb.push_back(phad.Pt()/(phad1.Pt()+phad2.Pt()));
            vb.push_back(plep.Pt());
            vb.push_back(deltaRlepjet1);
            vb.push_back(deltaRlepjet2);
            vb.push_back(deltaPhilepjet1);
            vb.push_back(deltaPhilepjet2);
            vb.push_back(deltaEtalepjet1);
            vb.push_back(deltaEtalepjet2);
            return vb;
      };
""")

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto btag_eff(Vint good, Vfloat btag, Vfloat pt, Vfloat eta, Vint flav, UInt_t flav_num, Float_t btag_limit) {
            int ind1 = 0;
            int ind2 = 0;
            vector<int> vb;
            for (unsigned int j=0; j<good.size(); ++j){
                        if(fabs(flav[good[j]]) == flav_num){
                                ind1++;
                                if(btag[good[j]] >= btag_limit) {
                                   ind2++;
                                }
                        }
            }
            vb.push_back(ind1);
            vb.push_back(ind2);
            return vb;
      };
      auto btag_eff_pt(Vint good, Vfloat btag, Vfloat pt, Vfloat eta, Vint flav, UInt_t flav_num, Float_t btag_limit, Float_t pt_limit_down, Float_t pt_limit_up) {
            int ind1 = 0;
            int ind2 = 0;
            vector<int> vb;
            for (unsigned int j=0; j<good.size(); ++j){
                        if(fabs(flav[good[j]]) == flav_num && pt[good[j]] > pt_limit_down && pt[good[j]] < pt_limit_up){
                                ind1++;
                                if(btag[good[j]] >= btag_limit) {
                                   ind2++;
                                }
                        }
            }
            vb.push_back(ind1);
            vb.push_back(ind2);
            return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

if mode == "data":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_nom','Jet_pt')
           df[s] = df[s].Define('Jet_mass_nom','Jet_mass')
           df[s] = df[s].Define('MET_smeared_phi','MET_phi')
           df[s] = df[s].Define('MET_smeared_pt','MET_pt')

df_muon = {}
df_electron = {}
df_test = {}

###########    Extra definitions

for s in samples:
        df[s] = df[s].Define('JetQInd','JetInds(nJet, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, JetBotInd, Jet_qgl)')
        ########### Filtering and further definitions
        df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>=4')
        df[s] = df[s].Filter('JetQInd.size() > 1')
        df[s] = df[s].Define('JetnotMuonInd','JetQInd[1]')
        ### hists definitions
        df[s] = df[s].Define('jet_1_pt','Jet_pt_nom[JetQInd[0]]')
        df[s] = df[s].Define('jet_1_nmu','Jet_nMuons[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_pt','nJetGood>1 ? Jet_pt_nom[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_1_eta','Jet_eta[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_2_mass','nJetGood>1 ? Jet_mass_nom[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_2_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_2_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
        df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Jet_pt_nom[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],0.) : 0')
        df[s] = df[s].Define('deltaR_jet1_jet2','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetQInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetQInd[0]])  : 10')
        df[s] = df[s].Define('deltaphi_jet1_jet2','fabs(Jet_phi[JetQInd[0]]-Jet_phi[JetnotMuonInd])')
        df[s] = df[s].Define('deltaeta_jet1_jet2','fabs(Jet_eta[JetQInd[0]]-Jet_eta[JetnotMuonInd])')
        df[s] = df[s].Define('deltapt_jet1_jet2','fabs(Jet_pt_nom[JetQInd[0]]-Jet_pt_nom[JetnotMuonInd])')
        df[s] = df[s].Define('tracks_jet1','Jet_nConstituents[JetQInd[0]]')
        df[s] = df[s].Define('tracks_jet2','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
        df[s] = df[s].Define('EMN_jet1','Jet_neEmEF[JetQInd[0]]')
        df[s] = df[s].Define('EMC_jet1','Jet_chEmEF[JetQInd[0]]')
        df[s] = df[s].Define('EMtotal_jet1','Jet_chEmEF[JetQInd[0]]+Jet_neEmEF[JetQInd[0]]')
        df[s] = df[s].Define('jet_bot1_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]]')
        df[s] = df[s].Define('jet_bot2_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]]')
        df[s] = df[s].Define('jet_1_btag','Jet_btagDeepFlavB[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_btag','Jet_btagDeepFlavB[JetnotMuonInd]')

############ Gen level definitions

if mode == "mc":
	for s in samples:
		df[s] = df[s].Define('jet_1_flavourH','Jet_hadronFlavour[JetQInd[0]]')
		df[s] = df[s].Define('jet_2_flavourH','Jet_hadronFlavour[JetnotMuonInd]')
		df[s] = df[s].Define('jet_1_flavourP','Jet_partonFlavour[JetQInd[0]]')
		df[s] = df[s].Define('jet_2_flavourP','Jet_partonFlavour[JetnotMuonInd]')
		df[s] = df[s].Define('jet_bot1_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[0]]]')
		df[s] = df[s].Define('jet_bot2_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[1]]]')
		df[s] = df[s].Define('btag_sf_aux1','(fabs(jet_bot1_flavourP) == 4 || fabs(jet_bot1_flavourP) == 5) ? btag_MED_sf[JetGoodInd[JetBotInd[0]]] : btag_MED_incl_sf[JetGoodInd[JetBotInd[0]]]')
		df[s] = df[s].Define('btag_sf_aux2','(fabs(jet_bot2_flavourP) == 4 || fabs(jet_bot2_flavourP) == 5) ? btag_LOO_sf[JetGoodInd[JetBotInd[1]]] : btag_LOO_incl_sf[JetGoodInd[JetBotInd[1]]]')
		df[s] = df[s].Define('btag_sf','btag_sf_aux1*btag_sf_aux2')

############################################################
#### Distinguishing between same sign and opposite sign ####
############################################################

#################### Final definitions and filters ###############################

for s in samples:
        df[s] = df[s].Define('weight_aux','1.')
        df[s] = df[s].Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetQInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetQInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetQInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetQInd[0]])')
        df[s] = df[s].Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df[s] = df[s].Define('InvM_muon_jet','nMuonGood>0 ? InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : 50.')
        ### Filters
        #df[s] = df[s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
        #df[s] = df[s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][0]))
        #df[s] = df[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
        df[s] = df[s].Filter('transverse_mass > 50')
        df[s] = df[s].Filter('InvM_2jets > 30. && InvM_2jets < 300.')

for s in samples:
       df[s] = df[s].Define('Aux_number_bottom_medium','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+')')
       df[s] = df[s].Define('Aux_number_bottom_loose','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+')')
       df[s] = df[s].Define('Aux_number_charm_medium','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+')')
       df[s] = df[s].Define('Aux_number_charm_loose','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+')')
       df[s] = df[s].Define('Aux_number_light_medium','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+')')
       df[s] = df[s].Define('Aux_number_light_loose','btag_eff(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+')')

       df[s] = df[s].Define('weight_bottom_medium_tot','Aux_number_bottom_medium[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg','Aux_number_bottom_medium[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot','Aux_number_bottom_loose[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg','Aux_number_bottom_loose[1]')
       df[s] = df[s].Define('weight_charm_medium_tot','Aux_number_charm_medium[0]')
       df[s] = df[s].Define('weight_charm_medium_btg','Aux_number_charm_medium[1]')
       df[s] = df[s].Define('weight_charm_loose_tot','Aux_number_charm_loose[0]')
       df[s] = df[s].Define('weight_charm_loose_btg','Aux_number_charm_loose[1]')
       df[s] = df[s].Define('weight_light_medium_tot','Aux_number_light_medium[0]')
       df[s] = df[s].Define('weight_light_medium_btg','Aux_number_light_medium[1]')
       df[s] = df[s].Define('weight_light_loose_tot','Aux_number_light_loose[0]')
       df[s] = df[s].Define('weight_light_loose_btg','Aux_number_light_loose[1]')

       df[s] = df[s].Define('eff_bottom_medium','weight_bottom_medium_tot > 0 ? weight_bottom_medium_btg/weight_bottom_medium_tot : 0.')
       df[s] = df[s].Define('eff_bottom_loose','weight_bottom_loose_tot > 0 ? weight_bottom_loose_btg/weight_bottom_loose_tot : 0.')
       df[s] = df[s].Define('eff_charm_medium','weight_charm_medium_tot > 0 ? weight_charm_medium_btg/weight_charm_medium_tot : 0.')
       df[s] = df[s].Define('eff_charm_loose','weight_charm_loose_tot > 0 ? weight_charm_loose_btg/weight_charm_loose_tot : 0.')
       df[s] = df[s].Define('eff_light_medium','weight_light_medium_tot > 0 ? weight_light_medium_btg/weight_light_medium_tot : 0.')
       df[s] = df[s].Define('eff_light_loose','weight_light_loose_tot > 0 ? weight_light_loose_btg/weight_light_loose_tot : 0.')

       df[s] = df[s].Define('Aux_numberjets_total','1.')
       df[s] = df[s].Define('Aux_numberjets_btagged','1.')

##############################################
### Calculating eff taking pT into account ###
##############################################

for s in samples:
       df[s] = df[s].Define('Aux_number_bottom_medium_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 20, 25)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 20, 25)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 20, 25)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 20, 25)')
       df[s] = df[s].Define('Aux_number_light_medium_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 20, 25)')
       df[s] = df[s].Define('Aux_number_light_loose_pt2025','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 20, 25)')

       df[s] = df[s].Define('Aux_number_bottom_medium_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 25, 30)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 25, 30)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 25, 30)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 25, 30)')
       df[s] = df[s].Define('Aux_number_light_medium_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 25, 30)')
       df[s] = df[s].Define('Aux_number_light_loose_pt2530','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 25, 30)')

       df[s] = df[s].Define('Aux_number_bottom_medium_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 30, 40)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 30, 40)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 30, 40)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 30, 40)')
       df[s] = df[s].Define('Aux_number_light_medium_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 30, 40)')
       df[s] = df[s].Define('Aux_number_light_loose_pt3040','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 30, 40)')

       df[s] = df[s].Define('Aux_number_bottom_medium_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 40, 50)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 40, 50)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 40, 50)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 40, 50)')
       df[s] = df[s].Define('Aux_number_light_medium_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 40, 50)')
       df[s] = df[s].Define('Aux_number_light_loose_pt4050','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 40, 50)')

       df[s] = df[s].Define('Aux_number_bottom_medium_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 50, 60)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 50, 60)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 50, 60)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 50, 60)')
       df[s] = df[s].Define('Aux_number_light_medium_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 50, 60)')
       df[s] = df[s].Define('Aux_number_light_loose_pt5060','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 50, 60)')

       df[s] = df[s].Define('Aux_number_bottom_medium_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][1])+', 60, 300)')
       df[s] = df[s].Define('Aux_number_bottom_loose_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 5,'+str(cuts_btag[years[0]][0])+', 60, 300)')
       df[s] = df[s].Define('Aux_number_charm_medium_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][1])+', 60, 300)')
       df[s] = df[s].Define('Aux_number_charm_loose_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 4,'+str(cuts_btag[years[0]][0])+', 60, 300)')
       df[s] = df[s].Define('Aux_number_light_medium_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][1])+', 60, 300)')
       df[s] = df[s].Define('Aux_number_light_loose_pt60','btag_eff_pt(JetGoodInd, Jet_btagDeepFlavB, Jet_pt, Jet_eta, Jet_hadronFlavour, 0,'+str(cuts_btag[years[0]][0])+', 60, 300)')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt2025','Aux_number_bottom_medium_pt2025[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt2025','Aux_number_bottom_medium_pt2025[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt2025','Aux_number_bottom_loose_pt2025[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt2025','Aux_number_bottom_loose_pt2025[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt2025','Aux_number_charm_medium_pt2025[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt2025','Aux_number_charm_medium_pt2025[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt2025','Aux_number_charm_loose_pt2025[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt2025','Aux_number_charm_loose_pt2025[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt2025','Aux_number_light_medium_pt2025[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt2025','Aux_number_light_medium_pt2025[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt2025','Aux_number_light_loose_pt2025[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt2025','Aux_number_light_loose_pt2025[1]')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt2530','Aux_number_bottom_medium_pt2530[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt2530','Aux_number_bottom_medium_pt2530[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt2530','Aux_number_bottom_loose_pt2530[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt2530','Aux_number_bottom_loose_pt2530[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt2530','Aux_number_charm_medium_pt2530[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt2530','Aux_number_charm_medium_pt2530[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt2530','Aux_number_charm_loose_pt2530[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt2530','Aux_number_charm_loose_pt2530[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt2530','Aux_number_light_medium_pt2530[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt2530','Aux_number_light_medium_pt2530[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt2530','Aux_number_light_loose_pt2530[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt2530','Aux_number_light_loose_pt2530[1]')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt3040','Aux_number_bottom_medium_pt3040[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt3040','Aux_number_bottom_medium_pt3040[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt3040','Aux_number_bottom_loose_pt3040[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt3040','Aux_number_bottom_loose_pt3040[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt3040','Aux_number_charm_medium_pt3040[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt3040','Aux_number_charm_medium_pt3040[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt3040','Aux_number_charm_loose_pt3040[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt3040','Aux_number_charm_loose_pt3040[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt3040','Aux_number_light_medium_pt3040[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt3040','Aux_number_light_medium_pt3040[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt3040','Aux_number_light_loose_pt3040[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt3040','Aux_number_light_loose_pt3040[1]')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt4050','Aux_number_bottom_medium_pt4050[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt4050','Aux_number_bottom_medium_pt4050[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt4050','Aux_number_bottom_loose_pt4050[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt4050','Aux_number_bottom_loose_pt4050[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt4050','Aux_number_charm_medium_pt4050[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt4050','Aux_number_charm_medium_pt4050[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt4050','Aux_number_charm_loose_pt4050[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt4050','Aux_number_charm_loose_pt4050[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt4050','Aux_number_light_medium_pt4050[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt4050','Aux_number_light_medium_pt4050[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt4050','Aux_number_light_loose_pt4050[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt4050','Aux_number_light_loose_pt4050[1]')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt5060','Aux_number_bottom_medium_pt5060[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt5060','Aux_number_bottom_medium_pt5060[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt5060','Aux_number_bottom_loose_pt5060[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt5060','Aux_number_bottom_loose_pt5060[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt5060','Aux_number_charm_medium_pt5060[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt5060','Aux_number_charm_medium_pt5060[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt5060','Aux_number_charm_loose_pt5060[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt5060','Aux_number_charm_loose_pt5060[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt5060','Aux_number_light_medium_pt5060[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt5060','Aux_number_light_medium_pt5060[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt5060','Aux_number_light_loose_pt5060[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt5060','Aux_number_light_loose_pt5060[1]')

       df[s] = df[s].Define('weight_bottom_medium_tot_pt60','Aux_number_bottom_medium_pt60[0]')
       df[s] = df[s].Define('weight_bottom_medium_btg_pt60','Aux_number_bottom_medium_pt60[1]')
       df[s] = df[s].Define('weight_bottom_loose_tot_pt60','Aux_number_bottom_loose_pt60[0]')
       df[s] = df[s].Define('weight_bottom_loose_btg_pt60','Aux_number_bottom_loose_pt60[1]')
       df[s] = df[s].Define('weight_charm_medium_tot_pt60','Aux_number_charm_medium_pt60[0]')
       df[s] = df[s].Define('weight_charm_medium_btg_pt60','Aux_number_charm_medium_pt60[1]')
       df[s] = df[s].Define('weight_charm_loose_tot_pt60','Aux_number_charm_loose_pt60[0]')
       df[s] = df[s].Define('weight_charm_loose_btg_pt60','Aux_number_charm_loose_pt60[1]')
       df[s] = df[s].Define('weight_light_medium_tot_pt60','Aux_number_light_medium_pt60[0]')
       df[s] = df[s].Define('weight_light_medium_btg_pt60','Aux_number_light_medium_pt60[1]')
       df[s] = df[s].Define('weight_light_loose_tot_pt60','Aux_number_light_loose_pt60[0]')
       df[s] = df[s].Define('weight_light_loose_btg_pt60','Aux_number_light_loose_pt60[1]')

############################################################
####################     HISTS    ##########################
############################################################

hist_botflav_mediumbtag_jets_tot = {}
hist_botflav_mediumbtag_jets_btg = {}
hist_botflav_loosebtag_jets_tot = {}
hist_botflav_loosebtag_jets_btg = {}
hist_charmflav_mediumbtag_jets_tot = {}
hist_charmflav_mediumbtag_jets_btg = {}
hist_charmflav_loosebtag_jets_tot = {}
hist_charmflav_loosebtag_jets_btg = {}
hist_lightflav_mediumbtag_jets_tot = {}
hist_lightflav_mediumbtag_jets_btg = {}
hist_lightflav_loosebtag_jets_tot = {}
hist_lightflav_loosebtag_jets_btg = {}
hist_entries = {}

hist_botflav_mediumbtag_jets_tot_pt2025 = {}
hist_botflav_mediumbtag_jets_btg_pt2025 = {}
hist_botflav_loosebtag_jets_tot_pt2025 = {}
hist_botflav_loosebtag_jets_btg_pt2025 = {}
hist_charmflav_mediumbtag_jets_tot_pt2025 = {}
hist_charmflav_mediumbtag_jets_btg_pt2025 = {}
hist_charmflav_loosebtag_jets_tot_pt2025 = {}
hist_charmflav_loosebtag_jets_btg_pt2025 = {}
hist_lightflav_mediumbtag_jets_tot_pt2025 = {}
hist_lightflav_mediumbtag_jets_btg_pt2025 = {}
hist_lightflav_loosebtag_jets_tot_pt2025 = {}
hist_lightflav_loosebtag_jets_btg_pt2025 = {}

hist_botflav_mediumbtag_jets_tot_pt2530 = {}
hist_botflav_mediumbtag_jets_btg_pt2530 = {}
hist_botflav_loosebtag_jets_tot_pt2530 = {}
hist_botflav_loosebtag_jets_btg_pt2530 = {}
hist_charmflav_mediumbtag_jets_tot_pt2530 = {}
hist_charmflav_mediumbtag_jets_btg_pt2530 = {}
hist_charmflav_loosebtag_jets_tot_pt2530 = {}
hist_charmflav_loosebtag_jets_btg_pt2530 = {}
hist_lightflav_mediumbtag_jets_tot_pt2530 = {}
hist_lightflav_mediumbtag_jets_btg_pt2530 = {}
hist_lightflav_loosebtag_jets_tot_pt2530 = {}
hist_lightflav_loosebtag_jets_btg_pt2530 = {}

hist_botflav_mediumbtag_jets_tot_pt3040 = {}
hist_botflav_mediumbtag_jets_btg_pt3040 = {}
hist_botflav_loosebtag_jets_tot_pt3040 = {}
hist_botflav_loosebtag_jets_btg_pt3040 = {}
hist_charmflav_mediumbtag_jets_tot_pt3040 = {}
hist_charmflav_mediumbtag_jets_btg_pt3040 = {}
hist_charmflav_loosebtag_jets_tot_pt3040 = {}
hist_charmflav_loosebtag_jets_btg_pt3040 = {}
hist_lightflav_mediumbtag_jets_tot_pt3040 = {}
hist_lightflav_mediumbtag_jets_btg_pt3040 = {}
hist_lightflav_loosebtag_jets_tot_pt3040 = {}
hist_lightflav_loosebtag_jets_btg_pt3040 = {}

hist_botflav_mediumbtag_jets_tot_pt4050 = {}
hist_botflav_mediumbtag_jets_btg_pt4050 = {}
hist_botflav_loosebtag_jets_tot_pt4050 = {}
hist_botflav_loosebtag_jets_btg_pt4050 = {}
hist_charmflav_mediumbtag_jets_tot_pt4050 = {}
hist_charmflav_mediumbtag_jets_btg_pt4050 = {}
hist_charmflav_loosebtag_jets_tot_pt4050 = {}
hist_charmflav_loosebtag_jets_btg_pt4050 = {}
hist_lightflav_mediumbtag_jets_tot_pt4050 = {}
hist_lightflav_mediumbtag_jets_btg_pt4050 = {}
hist_lightflav_loosebtag_jets_tot_pt4050 = {}
hist_lightflav_loosebtag_jets_btg_pt4050 = {}

hist_botflav_mediumbtag_jets_tot_pt5060 = {}
hist_botflav_mediumbtag_jets_btg_pt5060 = {}
hist_botflav_loosebtag_jets_tot_pt5060 = {}
hist_botflav_loosebtag_jets_btg_pt5060 = {}
hist_charmflav_mediumbtag_jets_tot_pt5060 = {}
hist_charmflav_mediumbtag_jets_btg_pt5060 = {}
hist_charmflav_loosebtag_jets_tot_pt5060 = {}
hist_charmflav_loosebtag_jets_btg_pt5060 = {}
hist_lightflav_mediumbtag_jets_tot_pt5060 = {}
hist_lightflav_mediumbtag_jets_btg_pt5060 = {}
hist_lightflav_loosebtag_jets_tot_pt5060 = {}
hist_lightflav_loosebtag_jets_btg_pt5060 = {}

hist_botflav_mediumbtag_jets_tot_pt60 = {}
hist_botflav_mediumbtag_jets_btg_pt60 = {}
hist_botflav_loosebtag_jets_tot_pt60 = {}
hist_botflav_loosebtag_jets_btg_pt60 = {}
hist_charmflav_mediumbtag_jets_tot_pt60 = {}
hist_charmflav_mediumbtag_jets_btg_pt60 = {}
hist_charmflav_loosebtag_jets_tot_pt60 = {}
hist_charmflav_loosebtag_jets_btg_pt60 = {}
hist_lightflav_mediumbtag_jets_tot_pt60 = {}
hist_lightflav_mediumbtag_jets_btg_pt60 = {}
hist_lightflav_loosebtag_jets_tot_pt60 = {}
hist_lightflav_loosebtag_jets_btg_pt60 = {}

for s in samples:

  hist_botflav_mediumbtag_jets_tot[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot")
  hist_botflav_mediumbtag_jets_btg[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg")
  hist_botflav_loosebtag_jets_tot[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot")
  hist_botflav_loosebtag_jets_btg[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg")

  hist_charmflav_mediumbtag_jets_tot[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot")
  hist_charmflav_mediumbtag_jets_btg[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg")
  hist_charmflav_loosebtag_jets_tot[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot")
  hist_charmflav_loosebtag_jets_btg[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg")

  hist_lightflav_mediumbtag_jets_tot[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot")
  hist_lightflav_mediumbtag_jets_btg[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg")
  hist_lightflav_loosebtag_jets_tot[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot")
  hist_lightflav_loosebtag_jets_btg[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg")

  hist_botflav_mediumbtag_jets_tot_pt2025[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt2025")
  hist_botflav_mediumbtag_jets_btg_pt2025[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt2025")
  hist_botflav_loosebtag_jets_tot_pt2025[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt2025")
  hist_botflav_loosebtag_jets_btg_pt2025[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt2025")

  hist_charmflav_mediumbtag_jets_tot_pt2025[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt2025")
  hist_charmflav_mediumbtag_jets_btg_pt2025[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt2025")
  hist_charmflav_loosebtag_jets_tot_pt2025[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt2025")
  hist_charmflav_loosebtag_jets_btg_pt2025[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt2025")

  hist_lightflav_mediumbtag_jets_tot_pt2025[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt2025")
  hist_lightflav_mediumbtag_jets_btg_pt2025[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt2025")
  hist_lightflav_loosebtag_jets_tot_pt2025[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt2025","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt2025")
  hist_lightflav_loosebtag_jets_btg_pt2025[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt2025","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt2025")

  hist_botflav_mediumbtag_jets_tot_pt2530[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt2530")
  hist_botflav_mediumbtag_jets_btg_pt2530[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt2530")
  hist_botflav_loosebtag_jets_tot_pt2530[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt2530")
  hist_botflav_loosebtag_jets_btg_pt2530[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt2530")

  hist_charmflav_mediumbtag_jets_tot_pt2530[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt2530")
  hist_charmflav_mediumbtag_jets_btg_pt2530[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt2530")
  hist_charmflav_loosebtag_jets_tot_pt2530[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt2530")
  hist_charmflav_loosebtag_jets_btg_pt2530[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt2530")

  hist_lightflav_mediumbtag_jets_tot_pt2530[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt2530")
  hist_lightflav_mediumbtag_jets_btg_pt2530[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt2530")
  hist_lightflav_loosebtag_jets_tot_pt2530[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt2530","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt2530")
  hist_lightflav_loosebtag_jets_btg_pt2530[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt2530","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt2530")

  hist_botflav_mediumbtag_jets_tot_pt3040[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt3040")
  hist_botflav_mediumbtag_jets_btg_pt3040[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt3040")
  hist_botflav_loosebtag_jets_tot_pt3040[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt3040")
  hist_botflav_loosebtag_jets_btg_pt3040[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt3040")

  hist_charmflav_mediumbtag_jets_tot_pt3040[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt3040")
  hist_charmflav_mediumbtag_jets_btg_pt3040[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt3040")
  hist_charmflav_loosebtag_jets_tot_pt3040[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt3040")
  hist_charmflav_loosebtag_jets_btg_pt3040[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt3040")

  hist_lightflav_mediumbtag_jets_tot_pt3040[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt3040")
  hist_lightflav_mediumbtag_jets_btg_pt3040[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt3040")
  hist_lightflav_loosebtag_jets_tot_pt3040[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt3040","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt3040")
  hist_lightflav_loosebtag_jets_btg_pt3040[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt3040","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt3040")

  hist_botflav_mediumbtag_jets_tot_pt4050[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt4050")
  hist_botflav_mediumbtag_jets_btg_pt4050[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt4050")
  hist_botflav_loosebtag_jets_tot_pt4050[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt4050")
  hist_botflav_loosebtag_jets_btg_pt4050[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt4050")

  hist_charmflav_mediumbtag_jets_tot_pt4050[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt4050")
  hist_charmflav_mediumbtag_jets_btg_pt4050[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt4050")
  hist_charmflav_loosebtag_jets_tot_pt4050[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt4050")
  hist_charmflav_loosebtag_jets_btg_pt4050[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt4050")

  hist_lightflav_mediumbtag_jets_tot_pt4050[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt4050")
  hist_lightflav_mediumbtag_jets_btg_pt4050[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt4050")
  hist_lightflav_loosebtag_jets_tot_pt4050[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt4050","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt4050")
  hist_lightflav_loosebtag_jets_btg_pt4050[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt4050","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt4050")

  hist_botflav_mediumbtag_jets_tot_pt5060[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt5060")
  hist_botflav_mediumbtag_jets_btg_pt5060[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt5060")
  hist_botflav_loosebtag_jets_tot_pt5060[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt5060")
  hist_botflav_loosebtag_jets_btg_pt5060[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt5060")

  hist_charmflav_mediumbtag_jets_tot_pt5060[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt5060")
  hist_charmflav_mediumbtag_jets_btg_pt5060[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt5060")
  hist_charmflav_loosebtag_jets_tot_pt5060[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt5060")
  hist_charmflav_loosebtag_jets_btg_pt5060[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt5060")

  hist_lightflav_mediumbtag_jets_tot_pt5060[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt5060")
  hist_lightflav_mediumbtag_jets_btg_pt5060[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt5060")
  hist_lightflav_loosebtag_jets_tot_pt5060[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt5060","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt5060")
  hist_lightflav_loosebtag_jets_btg_pt5060[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt5060","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt5060")

  hist_botflav_mediumbtag_jets_tot_pt60[s] = df[s].Histo1D(("botflav_mediumbtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_bottom_medium_tot_pt60")
  hist_botflav_mediumbtag_jets_btg_pt60[s] = df[s].Histo1D(("botflav_mediumbtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_bottom_medium_btg_pt60")
  hist_botflav_loosebtag_jets_tot_pt60[s] = df[s].Histo1D(("botflav_loosebtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_bottom_loose_tot_pt60")
  hist_botflav_loosebtag_jets_btg_pt60[s] = df[s].Histo1D(("botflav_loosebtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_bottom_loose_btg_pt60")

  hist_charmflav_mediumbtag_jets_tot_pt60[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_charm_medium_tot_pt60")
  hist_charmflav_mediumbtag_jets_btg_pt60[s] = df[s].Histo1D(("charmflav_mediumbtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_charm_medium_btg_pt60")
  hist_charmflav_loosebtag_jets_tot_pt60[s] = df[s].Histo1D(("charmflav_loosebtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_charm_loose_tot_pt60")
  hist_charmflav_loosebtag_jets_btg_pt60[s] = df[s].Histo1D(("charmflav_loosebtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_charm_loose_btg_pt60")

  hist_lightflav_mediumbtag_jets_tot_pt60[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_light_medium_tot_pt60")
  hist_lightflav_mediumbtag_jets_btg_pt60[s] = df[s].Histo1D(("lightflav_mediumbtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_light_medium_btg_pt60")
  hist_lightflav_loosebtag_jets_tot_pt60[s] = df[s].Histo1D(("lightflav_loosebtag_jets_tot_pt60","",5,0,5),"Aux_numberjets_total","weight_light_loose_tot_pt60")
  hist_lightflav_loosebtag_jets_btg_pt60[s] = df[s].Histo1D(("lightflav_loosebtag_jets_btg_pt60","",5,0,5),"Aux_numberjets_btagged","weight_light_loose_btg_pt60")

  hist_entries[s] = df[s].Histo1D(("entries_mine","",5,0,5),"Aux_numberjets_btagged") 

#############################
####     DATA SAVING     ####
#############################

for s in samples:
                path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/btag_eff_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_botflav_mediumbtag_jets_tot[s].Write()
                hist_botflav_mediumbtag_jets_btg[s].Write()
                hist_botflav_loosebtag_jets_tot[s].Write()
                hist_botflav_loosebtag_jets_btg[s].Write()
                hist_charmflav_mediumbtag_jets_tot[s].Write()
                hist_charmflav_mediumbtag_jets_btg[s].Write()
                hist_charmflav_loosebtag_jets_tot[s].Write()
                hist_charmflav_loosebtag_jets_btg[s].Write()
                hist_lightflav_mediumbtag_jets_tot[s].Write()
                hist_lightflav_mediumbtag_jets_btg[s].Write()
                hist_lightflav_loosebtag_jets_tot[s].Write()
                hist_lightflav_loosebtag_jets_btg[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt2025[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt2025[s].Write()
                hist_botflav_loosebtag_jets_tot_pt2025[s].Write()
                hist_botflav_loosebtag_jets_btg_pt2025[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt2025[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt2025[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt2025[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt2025[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt2025[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt2025[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt2025[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt2025[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt2530[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt2530[s].Write()
                hist_botflav_loosebtag_jets_tot_pt2530[s].Write()
                hist_botflav_loosebtag_jets_btg_pt2530[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt2530[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt2530[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt2530[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt2530[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt2530[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt2530[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt2530[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt2530[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt3040[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt3040[s].Write()
                hist_botflav_loosebtag_jets_tot_pt3040[s].Write()
                hist_botflav_loosebtag_jets_btg_pt3040[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt3040[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt3040[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt3040[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt3040[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt3040[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt3040[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt3040[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt3040[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt4050[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt4050[s].Write()
                hist_botflav_loosebtag_jets_tot_pt4050[s].Write()
                hist_botflav_loosebtag_jets_btg_pt4050[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt4050[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt4050[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt4050[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt4050[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt4050[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt4050[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt4050[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt4050[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt5060[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt5060[s].Write()
                hist_botflav_loosebtag_jets_tot_pt5060[s].Write()
                hist_botflav_loosebtag_jets_btg_pt5060[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt5060[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt5060[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt5060[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt5060[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt5060[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt5060[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt5060[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt5060[s].Write()
                hist_botflav_mediumbtag_jets_tot_pt60[s].Write()
                hist_botflav_mediumbtag_jets_btg_pt60[s].Write()
                hist_botflav_loosebtag_jets_tot_pt60[s].Write()
                hist_botflav_loosebtag_jets_btg_pt60[s].Write()
                hist_charmflav_mediumbtag_jets_tot_pt60[s].Write()
                hist_charmflav_mediumbtag_jets_btg_pt60[s].Write()
                hist_charmflav_loosebtag_jets_tot_pt60[s].Write()
                hist_charmflav_loosebtag_jets_btg_pt60[s].Write()
                hist_lightflav_mediumbtag_jets_tot_pt60[s].Write()
                hist_lightflav_mediumbtag_jets_btg_pt60[s].Write()
                hist_lightflav_loosebtag_jets_tot_pt60[s].Write()
                hist_lightflav_loosebtag_jets_btg_pt60[s].Write()

                myfile.Close()

for s in samples:
   print("Sample is "+s)

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg[s].Integral()/hist_botflav_mediumbtag_jets_tot[s].Integral(), hist_botflav_loosebtag_jets_btg[s].Integral()/hist_botflav_loosebtag_jets_tot[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg[s].Integral()/hist_charmflav_mediumbtag_jets_tot[s].Integral(), hist_charmflav_loosebtag_jets_btg[s].Integral()/hist_charmflav_loosebtag_jets_tot[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg[s].Integral()/hist_lightflav_mediumbtag_jets_tot[s].Integral(), hist_lightflav_loosebtag_jets_btg[s].Integral()/hist_lightflav_loosebtag_jets_tot[s].Integral()))

   print("%s entries" %(hist_entries[s].Integral()))
   print("%s entries from df" %(df[s].Count().GetValue()))

   print("--------------------------------------------------")

   print("pT between 20 and 25")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt2025[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt2025[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt2025[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt2025[s].Integral(), hist_botflav_loosebtag_jets_btg_pt2025[s].Integral()/hist_botflav_loosebtag_jets_tot_pt2025[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt2025[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt2025[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt2025[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt2025[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt2025[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt2025[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt2025[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt2025[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt2025[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt2025[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt2025[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt2025[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt2025[s].Integral()))

   print("--------------------------------------------------")

   print("pT between 25 and 30")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt2530[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt2530[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt2530[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt2530[s].Integral(), hist_botflav_loosebtag_jets_btg_pt2530[s].Integral()/hist_botflav_loosebtag_jets_tot_pt2530[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt2530[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt2530[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt2530[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt2530[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt2530[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt2530[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt2530[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt2530[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt2530[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt2530[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt2530[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt2530[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt2530[s].Integral()))

   print("--------------------------------------------------")

   print("pT between 30 and 40")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt3040[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt3040[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt3040[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt3040[s].Integral(), hist_botflav_loosebtag_jets_btg_pt3040[s].Integral()/hist_botflav_loosebtag_jets_tot_pt3040[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt3040[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt3040[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt3040[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt3040[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt3040[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt3040[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt3040[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt3040[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt3040[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt3040[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt3040[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt3040[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt3040[s].Integral()))

   print("--------------------------------------------------")

   print("pT between 40 and 50")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt4050[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt4050[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt4050[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt4050[s].Integral(), hist_botflav_loosebtag_jets_btg_pt4050[s].Integral()/hist_botflav_loosebtag_jets_tot_pt4050[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt4050[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt4050[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt4050[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt4050[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt4050[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt4050[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt4050[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt4050[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt4050[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt4050[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt4050[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt4050[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt4050[s].Integral()))

   print("--------------------------------------------------")

   print("pT between 50 and 60")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt5060[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt5060[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt5060[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt5060[s].Integral(), hist_botflav_loosebtag_jets_btg_pt5060[s].Integral()/hist_botflav_loosebtag_jets_tot_pt5060[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt5060[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt5060[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt5060[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt5060[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt5060[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt5060[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt5060[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt5060[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt5060[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt5060[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt5060[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt5060[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt5060[s].Integral()))

   print("--------------------------------------------------")

   print("pT over 60")

   print("%s jets in total with bottom flav" %(hist_botflav_mediumbtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag medium  with bottom flav" %(hist_botflav_mediumbtag_jets_btg_pt60[s].Integral()))
   print("%s jets in total with bottom flav" %(hist_botflav_loosebtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag loose with bottom flav" %(hist_botflav_loosebtag_jets_btg_pt60[s].Integral()))

   print("Efficiency for bottom jets is %s for medium WP and %s for loose WP" %(hist_botflav_mediumbtag_jets_btg_pt60[s].Integral()/hist_botflav_mediumbtag_jets_tot_pt60[s].Integral(), hist_botflav_loosebtag_jets_btg_pt60[s].Integral()/hist_botflav_loosebtag_jets_tot_pt60[s].Integral()))

   print("%s jets in total with charm flav" %(hist_charmflav_mediumbtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag medium  with charm flav" %(hist_charmflav_mediumbtag_jets_btg_pt60[s].Integral()))
   print("%s jets in total with charm flav" %(hist_charmflav_loosebtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag loose with charm flav" %(hist_charmflav_loosebtag_jets_btg_pt60[s].Integral()))

   print("Efficiency for charm jets is %s for medium WP and %s for loose WP" %(hist_charmflav_mediumbtag_jets_btg_pt60[s].Integral()/hist_charmflav_mediumbtag_jets_tot_pt60[s].Integral(), hist_charmflav_loosebtag_jets_btg_pt60[s].Integral()/hist_charmflav_loosebtag_jets_tot_pt60[s].Integral()))

   print("%s jets in total with light flav" %(hist_lightflav_mediumbtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag medium  with light flav" %(hist_lightflav_mediumbtag_jets_btg_pt60[s].Integral()))
   print("%s jets in total with light flav" %(hist_lightflav_loosebtag_jets_tot_pt60[s].Integral()))
   print("%s jets passing btag loose with light flav" %(hist_lightflav_loosebtag_jets_btg_pt60[s].Integral()))

   print("Efficiency for light jets is %s for medium WP and %s for loose WP" %(hist_lightflav_mediumbtag_jets_btg_pt60[s].Integral()/hist_lightflav_mediumbtag_jets_tot_pt60[s].Integral(), hist_lightflav_loosebtag_jets_btg_pt60[s].Integral()/hist_lightflav_loosebtag_jets_tot_pt60[s].Integral()))

   print("--------------------------------------------------")

print('Ended succesfully')

