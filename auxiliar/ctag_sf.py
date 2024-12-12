print('CTAG scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join
import pandas as pd

import json
import argparse

filecsv = {}

yearsL = ["2016pre","2016post","2017","2018"]

filecsv["2016pre"] = pd.read_csv('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/ctagging_wp_deepJet2016.csv')
filecsv["2016post"] = pd.read_csv('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/ctagging_wp_deepJet2016B.csv')
filecsv["2017"] = pd.read_csv('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/ctagging_wp_deepJet2017.csv')
filecsv["2018"] = pd.read_csv('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/ctagging_wp_deepJet2018.csv')

for ye in yearsL:
  filecsv[ye] = filecsv[ye][filecsv[ye]["OperatingPoint"] == "T"]

flavL = [0,4,5]

ptbins = {}
ptbinsX = {}

ptbins20174 = [20.0,30.0,40.0,50.0,110.0,210.0]
ptbins20175 = [30.0,50.0,70.0,100.0,140.0,200.0]
ptbins20170 = [20.0,1000.0]
ptbins2016pre4 = [25.0,30.0,50.0,70.0,100.0,130.0,160.0,250.0]
ptbins2016pre5 = [30.0,50.0,70.0,100.0,140.0,200.0]
ptbins2016pre0 = [20.0,1000.0]
ptbins2016post4 = [25.0,30.0,50.0,70.0,100.0,130.0,160.0,250.0]
ptbins2016post5 = [30.0,50.0,70.0,100.0,140.0,200.0]
ptbins2016post0 = [20.0,1000.0]
ptbins20184 = [25.0,30.0,50.0,70.0,100.0,130.0,160.0,250.0]
ptbins20185 = [30.0,50.0,70.0,100.0,140.0,200.0]
ptbins20180 = [20.0,1000.0]

ptbins["0"] = ptbins20180[:-1]; ptbins["4"] = ptbins20184[:-1]; ptbins["5"] = ptbins20185[:-1];
ptbinsX["0"] = ptbins20170[:-1]; ptbinsX["4"] = ptbins20174[:-1]; ptbinsX["5"] = ptbins20175[:-1];

dict_to_cpp2018_tag = "{"
for k,keta in enumerate([4,5]):
  dict_to_cpp2018_tag += "{"
  for i,ptval in enumerate(ptbins[str(keta)]):
     file_aux = filecsv["2018"][(filecsv["2018"]["jetFlavor"] == keta) & (filecsv["2018"]["ptMin"] == ptval)];
     indcent = file_aux[file_aux["sysType"]=="central"].index[0]
     indup = file_aux[file_aux["sysType"]=="up"].index[0]
     inddown = file_aux[file_aux["sysType"]=="down"].index[0]
     dict_to_cpp2018_tag += "{%s,%s,%s} %s" % (file_aux[file_aux["sysType"]=="central"].at[indcent,"formula"],file_aux[file_aux["sysType"]=="up"].at[indup,"formula"],file_aux[file_aux["sysType"]=="down"].at[inddown,"formula"],
      ("," if i< len(ptbins[str(keta)])-1 else "}"))
  dict_to_cpp2018_tag +="%s" % (", " if k < len([4,5]) - 1 else "}")

dict_to_cpp2017_tag = "{"
for k,keta in enumerate([4,5]):
  dict_to_cpp2017_tag += "{"
  for i,ptval in enumerate(ptbinsX[str(keta)]):
     file_aux = filecsv["2017"][(filecsv["2017"]["jetFlavor"] == keta) & (filecsv["2017"]["ptMin"] == ptval)];
     indcent = file_aux[file_aux["sysType"]=="central"].index[0]
     indup = file_aux[file_aux["sysType"]=="up"].index[0]
     inddown = file_aux[file_aux["sysType"]=="down"].index[0]
     dict_to_cpp2017_tag += "{%s,%s,%s} %s" % (file_aux[file_aux["sysType"]=="central"].at[indcent,"formula"],file_aux[file_aux["sysType"]=="up"].at[indup,"formula"],file_aux[file_aux["sysType"]=="down"].at[inddown,"formula"],
      ("," if i< len(ptbinsX[str(keta)])-1 else "}"))
  dict_to_cpp2017_tag +="%s" % (", " if k < len([4,5]) - 1 else "}")

dict_to_cpp2016_tag = "{"
for k,keta in enumerate([4,5]):
  dict_to_cpp2016_tag += "{"
  for i,ptval in enumerate(ptbins[str(keta)]):
     file_aux = filecsv["2016pre"][(filecsv["2016pre"]["jetFlavor"] == keta) & (filecsv["2016pre"]["ptMin"] == ptval)];
     indcent = file_aux[file_aux["sysType"]=="central"].index[0]
     indup = file_aux[file_aux["sysType"]=="up"].index[0]
     inddown = file_aux[file_aux["sysType"]=="down"].index[0]
     dict_to_cpp2016_tag += "{%s,%s,%s} %s" % (file_aux[file_aux["sysType"]=="central"].at[indcent,"formula"],file_aux[file_aux["sysType"]=="up"].at[indup,"formula"],file_aux[file_aux["sysType"]=="down"].at[inddown,"formula"],
      ("," if i< len(ptbins[str(keta)])-1 else "}"))
  dict_to_cpp2016_tag +="%s" % (", " if k < len([4,5]) - 1 else "}")

dict_to_cpp2016B_tag = "{"
for k,keta in enumerate([4,5]):
  dict_to_cpp2016B_tag += "{"
  for i,ptval in enumerate(ptbins[str(keta)]):
     file_aux = filecsv["2016post"][(filecsv["2016post"]["jetFlavor"] == keta) & (filecsv["2016post"]["ptMin"] == ptval)];
     indcent = file_aux[file_aux["sysType"]=="central"].index[0]
     indup = file_aux[file_aux["sysType"]=="up"].index[0]
     inddown = file_aux[file_aux["sysType"]=="down"].index[0]
     dict_to_cpp2016B_tag += "{%s,%s,%s} %s" % (file_aux[file_aux["sysType"]=="central"].at[indcent,"formula"],file_aux[file_aux["sysType"]=="up"].at[indup,"formula"],file_aux[file_aux["sysType"]=="down"].at[inddown,"formula"],
      ("," if i< len(ptbins[str(keta)])-1 else "}"))
  dict_to_cpp2016B_tag +="%s" % (", " if k < len([4,5]) - 1 else "}")

mistag1 = "{{0,{{1.24249,0.000710784,-4.7407e-7,3.50092},{0.158654,0.000101332,-1.32546e-7,0.0}}},"
mistag2 = "{1,{{0.412919,0.000974572,-5.40298e-7,7.53779},{0.210975,2.4705e-5,-2.79805e-8,0.0}}},"
mistag3 = "{2,{{1.17856,-0.00124081,7.9488e-7,-10.5641},{0.193428,0.000128445,-1.75534e-7,0.0}}},"
mistag4 = "{3,{{1.37717,-0.00107803,7.1253e-7,-8.46103},{0.165593,6.43651e-5,1.04852e-7,0.0}}}}"

dict_mistag = mistag1+mistag2+mistag3+mistag4

#### Comentamos aquí cómo se usaria
#### """ % eff_btag ) ## al final
###     std::map <int, std::vector<std::vector<std::vector<float>>>> eff_btag_all = %s;
###     std::map <int, std::vector<std::vector<std::vector<float>>>>::iterator it; ## al principio

dict_to_cpp = "{{0,"+dict_to_cpp2016_tag+"},{1,"+dict_to_cpp2016B_tag+"},{2,"+dict_to_cpp2017_tag+"},{3,"+dict_to_cpp2018_tag+"}}"

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_pol5(const float aa, const float bb, const float cc, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx;
            return res;
   };
   auto total_ctag_weight(const int jet_ind, Vint flav, Vfloat jetpt, const int yearT) {
     std::map <int, std::vector<std::vector<std::vector<float>>>> ctag_sf = %s;
     std::map <int, std::vector<std::vector<std::vector<float>>>>::iterator it;
     float weight_mc = 1.;
     auto flavranges = ctag_sf[yearT];
     if (fabs(flav[jet_ind]) == 4) {
        if (yearT == 2){
                 if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 30.) {
       	       	     weight_mc = flavranges[0][0][0];
                 } else if (jetpt[jet_ind] > 30. && jetpt[jet_ind] <= 40.) {
       	       	     weight_mc = flavranges[0][1][0];
                 } else if (jetpt[jet_ind] > 40. && jetpt[jet_ind] <= 50.) {
                     weight_mc = flavranges[0][2][0];
                 } else if (jetpt[jet_ind] > 50. && jetpt[jet_ind] <= 110.) {
                     weight_mc = flavranges[0][3][0];
                 } else if (jetpt[jet_ind] > 110. && jetpt[jet_ind] <= 210.) {
                     weight_mc = flavranges[0][4][0];
                 } else {
                     weight_mc = 1.;
                 }
        } else	{
                 if (jetpt[jet_ind] > 25. && jetpt[jet_ind] <= 30.) {
                     weight_mc = flavranges[0][0][0];
                 } else if (jetpt[jet_ind] > 30. && jetpt[jet_ind] <= 50.) {
                     weight_mc = flavranges[0][1][0];
                 } else if (jetpt[jet_ind] > 50. && jetpt[jet_ind] <= 70.) {
                     weight_mc = flavranges[0][2][0];
                 } else if (jetpt[jet_ind] > 70. && jetpt[jet_ind] <= 100.) {
                     weight_mc = flavranges[0][3][0];
                 } else if (jetpt[jet_ind] > 100. && jetpt[jet_ind] <= 130.) {
                     weight_mc = flavranges[0][4][0];
                 } else if (jetpt[jet_ind] > 130. && jetpt[jet_ind] <= 160.) {
                     weight_mc = flavranges[0][5][0];
                 } else if (jetpt[jet_ind] > 160. && jetpt[jet_ind] <= 250.) {
                     weight_mc = flavranges[0][6][0];
                 } else {
                     weight_mc = 1.;
                 } 
        }
     } else if (fabs(flav[jet_ind]) == 5) {
                 if (jetpt[jet_ind] > 30. && jetpt[jet_ind] <= 50.) {
                     weight_mc = flavranges[1][0][0];
                 } else if (jetpt[jet_ind] > 50. && jetpt[jet_ind] <= 70.) {
                     weight_mc = flavranges[1][1][0];
                 } else if (jetpt[jet_ind] > 70. && jetpt[jet_ind] <= 100.) {
                     weight_mc = flavranges[1][2][0];
                 } else if (jetpt[jet_ind] > 100. && jetpt[jet_ind] <= 140.) {
                     weight_mc = flavranges[1][3][0];
                 } else if (jetpt[jet_ind] > 140. && jetpt[jet_ind] <= 200.) {
                     weight_mc = flavranges[1][4][0];
                 } else {
                     weight_mc = 1.;
                 }
     }
     return weight_mc;
   };
""" % dict_to_cpp)

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_pol6(const float aa, const float bb, const float cc, const float dd, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx + dd/xx;
            return res;
   };
   auto total_cmistag_weight(const int jet_ind, Vint flav, Vfloat jetpt, const int yearT) {
     std::map <int, std::vector<std::vector<float>>> cmistag_sf = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it2;
     float weight_mc = 1.;
     auto flavranges = cmistag_sf[yearT];
     if (fabs(flav[jet_ind]) != 4 && fabs(flav[jet_ind]) != 5) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 1000.) {
            weight_mc = myfunc_pol6(flavranges[0][0], flavranges[0][1], flavranges[0][2], flavranges[0][3], jetpt[jet_ind]);
        } else {
            weight_mc = 1.;
        }
     }
     return weight_mc;
   }
""" % dict_mistag)

####################################################################################################################################################################################

class ctag_weights_tot():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('ctag_weight','total_ctag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 0)')
               df = df.Define('cmistag_weight','total_cmistag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 0)')
           elif self.year == "2016B":
               df = df.Define('ctag_weight','total_ctag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 1)')
               df = df.Define('cmistag_weight','total_cmistag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 1)')
           elif self.year == "2017":
               df = df.Define('ctag_weight','total_ctag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 2)')
               df = df.Define('cmistag_weight','total_cmistag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 2)')
           elif self.year == "2018":
               df = df.Define('ctag_weight','total_ctag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 3)')
               df = df.Define('cmistag_weight','total_cmistag_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 3)')

        variables = ['ctag_weight', 'cmistag_weight']

        branches = variables

        return df

def ctag_weights_totRDF(**kwargs):
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
    return lambda: ctag_weights_tot(**kwargs)

