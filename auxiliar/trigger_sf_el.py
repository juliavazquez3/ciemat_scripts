print('Trigger scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

trigger_sf_files_2016 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/trig_2016pretVFP.json"
trigger_sf_files_2016B = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/trig_2016postVFP.json"
trigger_sf_files_2017 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/trig_2017.json"
trigger_sf_files_2018 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/trig_2018.json"
fName_2016 = "SF_values_2016preVFP"
fName_2016B = "SF_values_2016postVFP"
fName_2017 = "SF_values_2017"
fName_2018 = "SF_values_2018"

file2017 = json.load(open(trigger_sf_files_2017))
file_32017 = file2017[fName_2017]

file2016 = json.load(open(trigger_sf_files_2016))
file_32016 = file2016[fName_2016]

file2016B = json.load(open(trigger_sf_files_2016B))
file_32016B = file2016B[fName_2016B]

file2018 = json.load(open(trigger_sf_files_2018))
file_32018 = file2018[fName_2018]

abseta_edges = ["abseta:[-2.5,-2]","abseta:[-2,-1.566]","abseta:[-1.566,-1.4442]","abseta:[-1.442,-0.8]","abseta:[-0.8,0]",
           "abseta:[0,0.8]","abseta:[0.8,1.4442]","abseta:[1.4442,1.566]","abseta:[1.566,2]","abseta:[2,2.5]"]
pt_edges1 = ["pt:[29,35]","pt:[35,50]","pt:[50,100]","pt:[100,200]","pt:[200,500]"]
pt_edges2 = ["pt:[34,50]","pt:[50,100]","pt:[100,200]","pt:[200,500]"]

pt_edges1_aux = [0,1,2,3,4]
pt_edges2_aux = [0,1,2,3]

dict_to_cpp2016 = "{"
for i,ipt in enumerate(pt_edges1):
  dict_to_cpp2016 += "{%s, {" % str(pt_edges1_aux[i])
  for j,iabseta in enumerate(abseta_edges):
  # if i == 0:
  # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
  # break
    dict_to_cpp2016 += "{%s,%s} %s" % (file_32016[ipt][iabseta][0],file_32016[ipt][iabseta][1],
        (", " if j < len(file_32016[ipt]) - 1 else ""))
  dict_to_cpp2016 += "}}%s" % (", " if i < len(file_32016) - 1 else "")
dict_to_cpp2016 += "}"

dict_to_cpp2016B = "{"
for i,ipt in enumerate(pt_edges1):
  dict_to_cpp2016B += "{%s, {" % str(pt_edges1_aux[i])
  for j,iabseta in enumerate(abseta_edges):
  # if i == 0:
  # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
  # break
    dict_to_cpp2016B += "{%s,%s} %s" % (file_32016B[ipt][iabseta][0],file_32016B[ipt][iabseta][1],
        (", " if j < len(file_32016B[ipt]) - 1 else ""))
  dict_to_cpp2016B += "}}%s" % (", " if i < len(file_32016B) - 1 else "")
dict_to_cpp2016B += "}"

dict_to_cpp2017 = "{"
for i,ipt in enumerate(pt_edges2):
  dict_to_cpp2017 += "{%s, {" % str(pt_edges2_aux[i])
  for j,iabseta in enumerate(abseta_edges):
  # if i == 0:
  # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
  # break
    dict_to_cpp2017 += "{%s,%s} %s" % (file_32017[ipt][iabseta][0],file_32017[ipt][iabseta][1],
        (", " if j < len(file_32017[ipt]) - 1 else ""))
  dict_to_cpp2017 += "}}%s" % (", " if i < len(file_32017) - 1 else "")
dict_to_cpp2017 += "}"

dict_to_cpp2018 = "{"
for i,ipt in enumerate(pt_edges2):
  dict_to_cpp2018 += "{%s, {" % str(pt_edges2_aux[i])
  for j,iabseta in enumerate(abseta_edges):
  # if i == 0:
  # dict_to_cpp += "{%s, %s} " % (value[0], value[1])
  # break
    dict_to_cpp2018 += "{%s,%s} %s" % (file_32018[ipt][iabseta][0],file_32018[ipt][iabseta][1],
        (", " if j < len(file_32018[ipt]) - 1 else ""))
  dict_to_cpp2018 += "}}%s" % (", " if i < len(file_32018) - 1 else "")
dict_to_cpp2018 += "}"

# Funciones para los sf del trigger
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sfel_2016(Vfloat pt, Vfloat eta, const string syst) {
     std::map <int, std::vector<std::vector<float>>> pt_eta_trig_sf2016 = %s;
     std::map <int, std::vector<float>>::iterator it;
     vector<float> vb;
     vector<float> vb_err;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (pt[i]>=29. && pt[i]<=35.) {
         auto ptranges = pt_eta_trig_sf2016[0];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=35. && pt[i]<=50.) {
         auto ptranges = pt_eta_trig_sf2016[1];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=50. && pt[i]<=100.) {
         auto ptranges = pt_eta_trig_sf2016[2];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=100. && pt[i]<=200.) {
         auto ptranges = pt_eta_trig_sf2016[3];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=200. && pt[i]<=500.) {
         auto ptranges = pt_eta_trig_sf2016[4];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else {
         vb.push_back(1.);
         vb_err.push_back(0.);
       }
     }
     if (syst == "nom") {
       return vb;
     } else {
       return vb_err;
     }
   }
""" % dict_to_cpp2016)

# Funciones para los sf del trigger
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sfel_2016B(Vfloat pt, Vfloat eta, const string syst) {
     std::map <int, std::vector<std::vector<float>>> pt_eta_trig_sf2016 = %s;
     std::map <int, std::vector<float>>::iterator it;
     vector<float> vb;
     vector<float> vb_err;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (pt[i]>=29. && pt[i]<=35.) {
         auto ptranges = pt_eta_trig_sf2016[0];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=35. && pt[i]<=50.) {
         auto ptranges = pt_eta_trig_sf2016[1];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=50. && pt[i]<=100.) {
         auto ptranges = pt_eta_trig_sf2016[2];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=100. && pt[i]<=200.) {
         auto ptranges = pt_eta_trig_sf2016[3];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=200. && pt[i]<=500.) {
         auto ptranges = pt_eta_trig_sf2016[4];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else {
         vb.push_back(1.);
         vb_err.push_back(0.);
       }
     }
     if (syst == "nom") {
       return vb;
     } else {
       return vb_err;
     }
   }
""" % dict_to_cpp2016B)

# Funciones para los sf del trigger
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sfel_2017(Vfloat pt, Vfloat eta, const string syst) {
     std::map <int, std::vector<std::vector<float>>> pt_eta_trig_sf2016 = %s;
     std::map <int, std::vector<float>>::iterator it;
     vector<float> vb;
     vector<float> vb_err;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (pt[i]>=34. && pt[i]<=50.) {
         auto ptranges = pt_eta_trig_sf2016[0];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=50. && pt[i]<=100.) {
         auto ptranges = pt_eta_trig_sf2016[1];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=100. && pt[i]<=200.) {
         auto ptranges = pt_eta_trig_sf2016[2];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=200. && pt[i]<=500.) {
         auto ptranges = pt_eta_trig_sf2016[3];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else {
         vb.push_back(1.);
         vb_err.push_back(0.);
       }
     }
     if (syst == "nom") {
       return vb;
     } else {
       return vb_err;
     }
   }
""" % dict_to_cpp2017)

# Funciones para los sf del trigger
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto trigger_sfel_2018(Vfloat pt, Vfloat eta, const string syst) {
     std::map <int, std::vector<std::vector<float>>> pt_eta_trig_sf2016 = %s;
     std::map <int, std::vector<float>>::iterator it;
     vector<float> vb;
     vector<float> vb_err;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (pt[i]>=34. && pt[i]<=50.) {
         auto ptranges = pt_eta_trig_sf2016[0];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=50. && pt[i]<=100.) {
         auto ptranges = pt_eta_trig_sf2016[1];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=100. && pt[i]<=200.) {
         auto ptranges = pt_eta_trig_sf2016[2];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else if (pt[i]>=200. && pt[i]<=500.) {
         auto ptranges = pt_eta_trig_sf2016[3];
         if (eta[i]>= -2.5 && eta[i]<= -2.0) {
           vb.push_back(ptranges[0][0]);
           vb_err.push_back(ptranges[0][1]);
         } else if (eta[i]>= -2.0 && eta[i]<= -1.566) {
           vb.push_back(ptranges[1][0]);
           vb_err.push_back(ptranges[1][1]);
         } else if (eta[i]>= -1.566 && eta[i]<= -1.4442) {
           vb.push_back(ptranges[2][0]);
           vb_err.push_back(ptranges[2][1]);
         } else if (eta[i]>= -1.4442 && eta[i]<= -0.8) {
           vb.push_back(ptranges[3][0]);
           vb_err.push_back(ptranges[3][1]);
         } else if (eta[i]>= -0.8 && eta[i]<= 0.0) {
           vb.push_back(ptranges[4][0]);
           vb_err.push_back(ptranges[4][1]);
         } else if (eta[i]>= 0.0 && eta[i]<= 0.8) {
           vb.push_back(ptranges[5][0]);
           vb_err.push_back(ptranges[5][1]);
         } else if (eta[i]>= 0.8 && eta[i]<= 1.4442) {
           vb.push_back(ptranges[6][0]);
           vb_err.push_back(ptranges[6][1]);
         } else if (eta[i]>= 1.4442 && eta[i]<= 1.566) {
           vb.push_back(ptranges[7][0]);
           vb_err.push_back(ptranges[7][1]);
         } else if (eta[i]>= 1.566 && eta[i]<= 2.0) {
           vb.push_back(ptranges[8][0]);
           vb_err.push_back(ptranges[8][1]);
         } else if (eta[i]>= 2.0 && eta[i]<= 2.5) {
           vb.push_back(ptranges[9][0]);
           vb_err.push_back(ptranges[9][1]);
         } else {
           vb.push_back(1.);
           vb_err.push_back(0.);
         }
       } else {
         vb.push_back(1.);
         vb_err.push_back(0.);
       }
     }
     if (syst == "nom") {
       return vb;
     } else {
       return vb_err;
     }
   }
""" % dict_to_cpp2018)

####################################################################################################################################################################################

class trigger_el_sf():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2017":
               df = df.Define('trigger_sf_el_aux','trigger_sfel_2017(Electron_pt, Electron_eta, "nom")')
               df = df.Define('trigger_sf_el_aux_err','trigger_sfel_2017(Electron_pt, Electron_eta, "err")')
           elif self.year == "2016":
               df = df.Define('trigger_sf_el_aux','trigger_sfel_2016(Electron_pt, Electron_eta, "nom")')
               df = df.Define('trigger_sf_el_aux_err','trigger_sfel_2016(Electron_pt, Electron_eta, "err")')
           elif self.year == "2016B":
               df = df.Define('trigger_sf_el_aux','trigger_sfel_2016B(Electron_pt, Electron_eta, "nom")')
               df = df.Define('trigger_sf_el_aux_err','trigger_sfel_2016B(Electron_pt, Electron_eta, "err")')
           elif self.year == "2018":
               df = df.Define('trigger_sf_el_aux','trigger_sfel_2018(Electron_pt, Electron_eta, "nom")')
               df = df.Define('trigger_sf_el_aux_err','trigger_sfel_2018(Electron_pt, Electron_eta, "err")')

        variables = ['trigger_sf_el_aux', 'trigger_sf_el_aux_err']

        branches = variables

        return df

def trigger_el_sfRDF(**kwargs):
    """
    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: jetVarRDF
            path: Base.Modules.smearing
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                isUL: self.dataset.has_tag('ul')
                ¿¿ proc: self.dataset.process ??
    """
    return lambda: trigger_el_sf(**kwargs)

