print('Trigger scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

trigger_sf_files_2016 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers_schemaV2.json"
trigger_sf_files_2016B = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers_schemaV2.json"
trigger_sf_files_2017 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers_schemaV2.json"
trigger_sf_files_2018 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers_schemaV2.json"
fName_2016 = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
fName_2016B = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
fName_2017 = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"
fName_2018 = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"

file2017 = json.load(open(trigger_sf_files_2017))
file_32017 = [item for item in file2017["corrections"] if item["name"] == fName_2017 and item["inputs"][0]["name"]=="abseta"][0]
pt_edges2017 = [29.,30.,40.,50.,60.,120.]

file2016 = json.load(open(trigger_sf_files_2016))
file_32016 = [item for item in file2016["corrections"] if item["name"] == fName_2016 and item["inputs"][0]["name"]=="abseta"][0]
pt_edges2016 = [26.,30.,40.,50.,60.,120.]

file2016B = json.load(open(trigger_sf_files_2016B))
file_32016B = [item for item in file2016B["corrections"] if item["name"] == fName_2016B and item["inputs"][0]["name"]=="abseta"][0]
pt_edges2016B = [26.,30.,40.,50.,60.,120.]

file2018 = json.load(open(trigger_sf_files_2018))
file_32018 = [item for item in file2018["corrections"] if item["name"] == fName_2018 and item["inputs"][0]["name"]=="abseta"][0]
pt_edges2018 = [26.,30.,40.,50.,60.,120.]

abseta_edges = [0,1,2,3]

dict_to_cpp2016 = "{"
for i,iabseta in enumerate(file_32016["data"]["content"]):
   dict_to_cpp2016 += "{%s, {" % abseta_edges[i]
   for j,ipt in enumerate(iabseta["content"]):
   # if i == 0:
   # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
   # break
     dict_to_cpp2016 += "{%s, %s}%s" % (pt_edges2016[j], ipt["content"][5]["value"],
         (", " if j < len(iabseta["content"]) - 1 else ""))
   dict_to_cpp2016 += "}}%s" % (", " if i < len(file_32016["data"]["content"]) - 1 else "")
dict_to_cpp2016 += "}"

dict_to_cpp2016B = "{"
for i,iabseta in enumerate(file_32016B["data"]["content"]):
   dict_to_cpp2016B += "{%s, {" % abseta_edges[i]
   for j,ipt in enumerate(iabseta["content"]):
   # if i == 0:
   # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
   # break
     dict_to_cpp2016B += "{%s, %s}%s" % (pt_edges2016B[j], ipt["content"][5]["value"],
         (", " if j < len(iabseta["content"]) - 1 else ""))
   dict_to_cpp2016B += "}}%s" % (", " if i < len(file_32016B["data"]["content"]) - 1 else "")
dict_to_cpp2016B += "}"

dict_to_cpp2017 = "{"
for i,iabseta in enumerate(file_32017["data"]["content"]):
   dict_to_cpp2017 += "{%s, {" % abseta_edges[i]
   for j,ipt in enumerate(iabseta["content"]):
   # if i == 0:
   # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
   # break
     dict_to_cpp2017 += "{%s, %s}%s" % (pt_edges2017[j], ipt["content"][5]["value"],
         (", " if j < len(iabseta["content"]) - 1 else ""))
   dict_to_cpp2017 += "}}%s" % (", " if i < len(file_32017["data"]["content"]) - 1 else "")
dict_to_cpp2017 += "}"

dict_to_cpp2018 = "{"
for i,iabseta in enumerate(file_32018["data"]["content"]):
   dict_to_cpp2018 += "{%s, {" % abseta_edges[i]
   for j,ipt in enumerate(iabseta["content"]):
   # if i == 0:
   # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
   # break
     dict_to_cpp2018 += "{%s, %s}%s" % (pt_edges2018[j], ipt["content"][5]["value"],
         (", " if j < len(iabseta["content"]) - 1 else ""))
   dict_to_cpp2018 += "}}%s" % (", " if i < len(file_32018["data"]["content"]) - 1 else "")
dict_to_cpp2018 += "}"

print("2016")
print(dict_to_cpp2016)
print("2016B")
print(dict_to_cpp2016B)
print("2017")
print(dict_to_cpp2017)
print("2018")
print(dict_to_cpp2018)

# Funciones para los sf del trigger
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sf_2016(Vfloat pt, Vfloat eta, float quant) {
     std::map <int, std::vector<std::vector<float>>> abseta_pt_trig_sf2016 = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it;
     vector<float> vb;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (fabs(eta[i]) <= 0.9) {
         auto ptranges = abseta_pt_trig_sf2016[0];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 0.9 && fabs(eta[i]) <= 1.2) {
         auto ptranges = abseta_pt_trig_sf2016[1];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 1.2 && fabs(eta[i]) <= 2.1) {
         auto ptranges = abseta_pt_trig_sf2016[2];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 2.1 && fabs(eta[i]) <= 2.4) {
         auto ptranges = abseta_pt_trig_sf2016[3];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else {
         vb.push_back(1.);
       }
     }
     return vb;
   }
""" % dict_to_cpp2016)

gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sf_2016B(Vfloat pt, Vfloat eta, float quant) {
     std::map <int, std::vector<std::vector<float>>> abseta_pt_trig_sf2016B = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it;
     vector<float> vb;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (fabs(eta[i]) <= 0.9) {
         auto ptranges = abseta_pt_trig_sf2016B[0];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 0.9 && fabs(eta[i]) <= 1.2) {
         auto ptranges = abseta_pt_trig_sf2016B[1];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 1.2 && fabs(eta[i]) <= 2.1) {
         auto ptranges = abseta_pt_trig_sf2016B[2];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 2.1 && fabs(eta[i]) <= 2.4) {
         auto ptranges = abseta_pt_trig_sf2016B[3];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else {
         vb.push_back(1.);
       }
     }
     return vb;
   }
""" % dict_to_cpp2016B)

gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sf_2017(Vfloat pt, Vfloat eta, float quant) {
     std::map <int, std::vector<std::vector<float>>> abseta_pt_trig_sf2017 = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it;
     vector<float> vb;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (fabs(eta[i]) <= 0.9) {
         auto ptranges = abseta_pt_trig_sf2017[0];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 0.9 && fabs(eta[i]) <= 1.2) {
         auto ptranges = abseta_pt_trig_sf2017[1];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 1.2 && fabs(eta[i]) <= 2.1) {
         auto ptranges = abseta_pt_trig_sf2017[2];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 2.1 && fabs(eta[i]) <= 2.4) {
         auto ptranges = abseta_pt_trig_sf2017[3];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else {
         vb.push_back(1.);
       }
     }
     return vb;
   }
""" % dict_to_cpp2017)

gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sf_2018(Vfloat pt, Vfloat eta, float quant) {
     std::map <int, std::vector<std::vector<float>>> abseta_pt_trig_sf2018 = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it;
     vector<float> vb;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (fabs(eta[i]) <= 0.9) {
         auto ptranges = abseta_pt_trig_sf2018[0];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 0.9 && fabs(eta[i]) <= 1.2) {
         auto ptranges = abseta_pt_trig_sf2018[1];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 1.2 && fabs(eta[i]) <= 2.1) {
         auto ptranges = abseta_pt_trig_sf2018[2];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 2.1 && fabs(eta[i]) <= 2.4) {
         auto ptranges = abseta_pt_trig_sf2018[3];
         if (pt[i] > quant && pt[i] <= 30.) {
           vb.push_back(ptranges[0][1]);
         } else if (pt[i] >= 30. && pt[i] < 40.) {
           vb.push_back(ptranges[1][1]);
         } else if (pt[i] >= 40. && pt[i] < 50.) {
           vb.push_back(ptranges[2][1]);
         } else if (pt[i] >= 50. && pt[i] < 60.) {
           vb.push_back(ptranges[3][1]);
         } else if (pt[i] >= 60. && pt[i] < 120.) {
           vb.push_back(ptranges[4][1]);
         } else if (pt[i] >= 120. && pt[i] < 200.) {
           vb.push_back(ptranges[5][1]);
         } else {
           vb.push_back(1.);
         }
       } else {
         vb.push_back(1.);
       }
     }
     return vb;
   }
""" % dict_to_cpp2018)

