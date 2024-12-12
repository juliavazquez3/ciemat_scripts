print('Low pt muons scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

muon_files_2016 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/NUM_DisplacedID_DEN_dSAMuons_abseta_pt_2016.json"
muon_files_2017 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/NUM_DisplacedID_DEN_dSAMuons_abseta_pt_2017.json"
muon_files_2018 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/NUM_DisplacedID_DEN_dSAMuons_abseta_pt_2018.json"
fName = "NUM_DisplacedID_DEN_dSAMuons"

file2016 = json.load(open(muon_files_2016))
file2017 = json.load(open(muon_files_2017))
file2018 = json.load(open(muon_files_2018))

table_2016 = file2016[fName]["abseta_pt"]
table_2017 = file2017[fName]["abseta_pt"]
table_2018 = file2018[fName]["abseta_pt"]

abseta_edges = ["abseta:[0.0,0.9]","abseta:[0.9,1.2]","abseta:[1.2,2.1]","abseta:[2.1,2.4]"]
pt_edges = ["pt:[3,4]","pt:[4,5]","pt:[5,6]","pt:[6,7]","pt:[7,8]","pt:[8,9]","pt:[9,10]","pt:[10,30]"]

dict_to_cpp2016 = "{"
for i,iabseta in enumerate(abseta_edges):
   dict_to_cpp2016 += "{"
   for j,ipt in enumerate(pt_edges):
     dict_to_cpp2016 += "%s %s" % (table_2016[iabseta][ipt]["value"],
         (", " if j < len(pt_edges) - 1 else ""))
   dict_to_cpp2016 += "}%s" % (", " if i < len(abseta_edges) - 1 else "")
dict_to_cpp2016 += "}"

dict_to_cpp2017 = "{"
for i,iabseta in enumerate(abseta_edges):
   dict_to_cpp2017 += "{"
   for j,ipt in enumerate(pt_edges):
     dict_to_cpp2017 += "%s %s" % (table_2017[iabseta][ipt]["value"],
         (", " if j < len(pt_edges) - 1 else ""))
   dict_to_cpp2017 += "}%s" % (", " if i < len(abseta_edges) - 1 else "")
dict_to_cpp2017 += "}"

dict_to_cpp2018 = "{"
for i,iabseta in enumerate(abseta_edges):
   dict_to_cpp2018 += "{"
   for j,ipt in enumerate(pt_edges):
     dict_to_cpp2018 += "%s %s" % (table_2018[iabseta][ipt]["value"],
         (", " if j < len(pt_edges) - 1 else ""))
   dict_to_cpp2018 += "}%s" % (", " if i < len(abseta_edges) - 1 else "")
dict_to_cpp2018 += "}"

new_dict = "{{0,"+dict_to_cpp2016+"},{1,"+dict_to_cpp2017+"},{2,"+dict_to_cpp2018+"}}"

gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto muon_lowpt(Vfloat pt, Vfloat eta, const int yearT) {
     std::map <int, std::vector<std::vector<float>>> first_tab = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it;
     auto final_table = first_tab[yearT];
     vector<float> vb;
     for (unsigned int i=0; i<eta.size(); ++i) {
       if (fabs(eta[i]) <= 0.9) {
         auto ptranges = final_table[0];
         if (pt[i] <= 4. && pt[i] > 3.) {
           vb.push_back(ptranges[0]);
         } else if (pt[i] >= 4. && pt[i] < 5.) {
           vb.push_back(ptranges[1]);
         } else if (pt[i] >= 5. && pt[i] < 6.) {
           vb.push_back(ptranges[2]);
         } else if (pt[i] >= 6. && pt[i] < 7.) {
           vb.push_back(ptranges[3]);
         } else if (pt[i] >= 7. && pt[i] < 8.) {
           vb.push_back(ptranges[4]);
         } else if (pt[i] >= 8. && pt[i] < 9.) {
           vb.push_back(ptranges[5]);
         } else if (pt[i] >= 9. && pt[i] < 10.) {
           vb.push_back(ptranges[6]);
         } else if (pt[i] >= 10. && pt[i] < 30.) {
           vb.push_back(ptranges[7]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 0.9 && fabs(eta[i]) <= 1.2) {
         auto ptranges = final_table[1];
         if (pt[i] <= 4. && pt[i] > 3.) {
           vb.push_back(ptranges[0]);
         } else if (pt[i] >= 4. && pt[i] < 5.) {
           vb.push_back(ptranges[1]);
         } else if (pt[i] >= 5. && pt[i] < 6.) {
           vb.push_back(ptranges[2]);
         } else if (pt[i] >= 6. && pt[i] < 7.) {
           vb.push_back(ptranges[3]);
         } else if (pt[i] >= 7. && pt[i] < 8.) {
           vb.push_back(ptranges[4]);
         } else if (pt[i] >= 8. && pt[i] < 9.) {
           vb.push_back(ptranges[5]);
         } else if (pt[i] >= 9. && pt[i] < 10.) {
           vb.push_back(ptranges[6]);
         } else if (pt[i] >= 10. && pt[i] < 30.) {
           vb.push_back(ptranges[7]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 1.2 && fabs(eta[i]) <= 2.1) {
         auto ptranges = final_table[2];
         if (pt[i] <= 4. && pt[i] > 3.) {
           vb.push_back(ptranges[0]);
         } else if (pt[i] >= 4. && pt[i] < 5.) {
           vb.push_back(ptranges[1]);
         } else if (pt[i] >= 5. && pt[i] < 6.) {
           vb.push_back(ptranges[2]);
         } else if (pt[i] >= 6. && pt[i] < 7.) {
           vb.push_back(ptranges[3]);
         } else if (pt[i] >= 7. && pt[i] < 8.) {
           vb.push_back(ptranges[4]);
         } else if (pt[i] >= 8. && pt[i] < 9.) {
           vb.push_back(ptranges[5]);
         } else if (pt[i] >= 9. && pt[i] < 10.) {
           vb.push_back(ptranges[6]);
         } else if (pt[i] >= 10. && pt[i] < 30.) {
           vb.push_back(ptranges[7]);
         } else {
           vb.push_back(1.);
         }
       } else if (fabs(eta[i]) > 2.1 && fabs(eta[i]) <= 2.4) {
         auto ptranges = final_table[3];
         if (pt[i] <= 4. && pt[i] > 3.) {
           vb.push_back(ptranges[0]);
         } else if (pt[i] >= 4. && pt[i] < 5.) {
           vb.push_back(ptranges[1]);
         } else if (pt[i] >= 5. && pt[i] < 6.) {
           vb.push_back(ptranges[2]);
         } else if (pt[i] >= 6. && pt[i] < 7.) {
           vb.push_back(ptranges[3]);
         } else if (pt[i] >= 7. && pt[i] < 8.) {
           vb.push_back(ptranges[4]);
         } else if (pt[i] >= 8. && pt[i] < 9.) {
           vb.push_back(ptranges[5]);
         } else if (pt[i] >= 9. && pt[i] < 10.) {
           vb.push_back(ptranges[6]);
         } else if (pt[i] >= 10. && pt[i] < 30.) {
           vb.push_back(ptranges[7]);
         } else {
           vb.push_back(1.);
         }
       } else {
         vb.push_back(1.);
       }
     }
     return vb;
   }
""" % new_dict)

####################################################################################################################################################################################

class displaced_mu_id():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('displaced_muon_low_id_sf','muon_lowpt(Muon_pt, Muon_eta, 0)')
           elif self.year == "2016B":
               df = df.Define('displaced_muon_low_id_sf','muon_lowpt(Muon_pt, Muon_eta, 0)')
           elif self.year == "2017":
               df = df.Define('displaced_muon_low_id_sf','muon_lowpt(Muon_pt, Muon_eta, 1)')
           elif self.year == "2018":
               df = df.Define('displaced_muon_low_id_sf','muon_lowpt(Muon_pt, Muon_eta, 2)')

        variables = ['displaced_muon_low_id_sf']

        branches = variables

        return df

def displaced_mu_idRDF(**kwargs):
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
    return lambda: displaced_mu_id(**kwargs)

