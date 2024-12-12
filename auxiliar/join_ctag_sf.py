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

## mistag1 is 2016 then 2016B,2017 and 2018
## first is flavour 0, second 4 and third 5
## values correspond to first coeff [0], second [1] and third [3] in equation [0]+[1]*exp([2]*x)

effmc1 = "{{0,{{4.13182e-02,-2.37595e-02,-2.36221e-02},{3.97975e-01,-3.75733e-01,-3.30244e-02},{4.09319e-01,-2.57558e-01,-1.51208e-02}}},"
effmc2 = "{1,{{4.17324e-02,-2.35894e-02,-2.32517e-02},{3.99713e-01,-3.73863e-01,-3.28992e-02},{4.12717e-01,-2.55883e-01,-1.51056e-02}}},"
effmc3 = "{2,{{2.70404e-01,-2.59587e-01,-7.16768e-05},{3.88679e-01,-4.04733e-01,-3.96979e-02},{3.48276e-01,-2.47078e-01,-2.50963e-02}}},"
effmc4 = "{3,{{2.81567e-01,-2.48132e-01,-2.02721e-04},{4.72195e-01,-4.11363e-01,-3.98227e-02},{4.33088e-01,-2.61989e-01,-1.94343e-02}}}}"

dict_effmc = effmc1+effmc2+effmc3+effmc4

#### Comentamos aquí cómo se usaria
#### """ % eff_btag ) ## al final
###     std::map <int, std::vector<std::vector<std::vector<float>>>> eff_btag_all = %s;
###     std::map <int, std::vector<std::vector<std::vector<float>>>>::iterator it; ## al principio

dict_flav45 = "{{0,"+dict_to_cpp2016_tag+"},{1,"+dict_to_cpp2016B_tag+"},{2,"+dict_to_cpp2017_tag+"},{3,"+dict_to_cpp2018_tag+"}}"

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_pol55(const float aa, const float bb, const float cc, const float dd, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx + dd/xx;
            return res;
   };
   auto myfunc_pol66(const float aa, const float bb, const float cc, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx;
            return res;
   };
   auto myfunc_expjoin(const float aa, const float bb, const float cc, const float xx) {
            float res;
            res = aa + bb*std::exp(cc*xx);
            return res;
   };
   auto joined_ctag_weight(Vint jet_ind, Vint flav, Vfloat jetpt, const int yearT, Vfloat cvl, Vfloat cvb, const string syst) {
     std::map <int, std::vector<std::vector<std::vector<float>>>> ctag_flav45_sf = %s;
     std::map <int, std::vector<std::vector<float>>> ctag_flav0_sf = %s;
     std::map <int, std::vector<std::vector<float>>> effmc_map = %s;
     std::vector<float> jet_sf;
     jet_sf.push_back(1.);
     jet_sf.push_back(1.);
     std::vector<float> jet_eff;
     jet_eff.push_back(0.5);
     jet_eff.push_back(0.5);
     auto flavranges = ctag_flav45_sf[yearT];
     auto flavrangeslight = ctag_flav0_sf[yearT];
     auto effmcranges = effmc_map[yearT];
     for (unsigned int i=0; i<2; ++i) {
       if (fabs(flav[jet_ind[i]]) == 4) {
          if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 200.) {
            jet_eff[i] = myfunc_expjoin(effmcranges[1][0], effmcranges[1][1], effmcranges[1][2], jetpt[jet_ind[i]]);
          } else {
            jet_eff[i] = myfunc_expjoin(effmcranges[1][0], effmcranges[1][1], effmcranges[1][2], 200.);
          }
          if (yearT == 2){
            if (syst == "nom"){
                 if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 30.) {
       	       	     jet_sf[i] = flavranges[0][0][0];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 40.) {
       	       	     jet_sf[i] = flavranges[0][1][0];
                 } else if (jetpt[jet_ind[i]] > 40. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][2][0];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 110.) {
                     jet_sf[i] = flavranges[0][3][0];
                 } else if (jetpt[jet_ind[i]] > 110. && jetpt[jet_ind[i]] <= 210.) {
                     jet_sf[i] = flavranges[0][4][0];
                 } else {
                     jet_sf[i] = 1.;
                 }
            } else if (syst == "up") {
                 if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 30.) {
                     jet_sf[i] = flavranges[0][0][1];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 40.) {
                     jet_sf[i] = flavranges[0][1][1];
                 } else if (jetpt[jet_ind[i]] > 40. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][2][1];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 110.) {
                     jet_sf[i] = flavranges[0][3][1];
                 } else if (jetpt[jet_ind[i]] > 110. && jetpt[jet_ind[i]] <= 210.) {
                     jet_sf[i] = flavranges[0][4][1];
                 } else {
                     jet_sf[i] = 1.;
                 }
            } else if (syst == "down") {
                 if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 30.) {
                     jet_sf[i] = flavranges[0][0][2];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 40.) {
                     jet_sf[i] = flavranges[0][1][2];
                 } else if (jetpt[jet_ind[i]] > 40. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][2][2];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 110.) {
                     jet_sf[i] = flavranges[0][3][2];
                 } else if (jetpt[jet_ind[i]] > 110. && jetpt[jet_ind[i]] <= 210.) {
                     jet_sf[i] = flavranges[0][4][2];
                 } else {
                     jet_sf[i] = 1.;
                 }
            }
          } else {
            if (syst == "nom") {
                 if (jetpt[jet_ind[i]] > 25. && jetpt[jet_ind[i]] <= 30.) {
                     jet_sf[i] = flavranges[0][0][0];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][1][0];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
                     jet_sf[i] = flavranges[0][2][0];
                 } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
                     jet_sf[i] = flavranges[0][3][0];
                 } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 130.) {
                     jet_sf[i] = flavranges[0][4][0];
                 } else if (jetpt[jet_ind[i]] > 130. && jetpt[jet_ind[i]] <= 160.) {
                     jet_sf[i] = flavranges[0][5][0];
                 } else if (jetpt[jet_ind[i]] > 160. && jetpt[jet_ind[i]] <= 250.) {
                     jet_sf[i] = flavranges[0][6][0];
                 } else {
                     jet_sf[i] = 1.;
                 } 
            } else if (syst == "up") {
                 if (jetpt[jet_ind[i]] > 25. && jetpt[jet_ind[i]] <= 30.) {
                     jet_sf[i] = flavranges[0][0][1];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][1][1];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
                     jet_sf[i] = flavranges[0][2][1];
                 } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
                     jet_sf[i] = flavranges[0][3][1];
                 } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 130.) {
                     jet_sf[i] = flavranges[0][4][1];
                 } else if (jetpt[jet_ind[i]] > 130. && jetpt[jet_ind[i]] <= 160.) {
                     jet_sf[i] = flavranges[0][5][1];
                 } else if (jetpt[jet_ind[i]] > 160. && jetpt[jet_ind[i]] <= 250.) {
                     jet_sf[i] = flavranges[0][6][1];
                 } else {
                     jet_sf[i] = 1.;
                 }
            } else if (syst == "down") {
                 if (jetpt[jet_ind[i]] > 25. && jetpt[jet_ind[i]] <= 30.) {
                     jet_sf[i] = flavranges[0][0][2];
                 } else if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
                     jet_sf[i] = flavranges[0][1][2];
                 } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
                     jet_sf[i] = flavranges[0][2][2];
                 } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
                     jet_sf[i] = flavranges[0][3][2];
                 } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 130.) {
                     jet_sf[i] = flavranges[0][4][2];
                 } else if (jetpt[jet_ind[i]] > 130. && jetpt[jet_ind[i]] <= 160.) {
                     jet_sf[i] = flavranges[0][5][2];
                 } else if (jetpt[jet_ind[i]] > 160. && jetpt[jet_ind[i]] <= 250.) {
                     jet_sf[i] = flavranges[0][6][2];
                 } else {
                     jet_sf[i] = 1.;
                 }
            }
          }
       } else if (fabs(flav[jet_ind[i]]) == 5) {
          if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 200.) {
            jet_eff[i] = myfunc_expjoin(effmcranges[2][0], effmcranges[2][1], effmcranges[2][2], jetpt[jet_ind[i]]);
          } else {
            jet_eff[i] = myfunc_expjoin(effmcranges[2][0], effmcranges[2][1], effmcranges[2][2], 200.);
          }
          if (syst == "nom"){
            if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
               jet_sf[i] = flavranges[1][0][0];
            } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
               jet_sf[i] = flavranges[1][1][0];
            } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
               jet_sf[i] = flavranges[1][2][0];
            } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 140.) {
               jet_sf[i] = flavranges[1][3][0];
            } else if (jetpt[jet_ind[i]] > 140. && jetpt[jet_ind[i]] <= 200.) {
               jet_sf[i] = flavranges[1][4][0];
            } else {
               jet_sf[i] = 1.;
            }
          } else if (syst == "up"){
            if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
               jet_sf[i] = flavranges[1][0][1];
            } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
               jet_sf[i] = flavranges[1][1][1];
            } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
               jet_sf[i] = flavranges[1][2][1];
            } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 140.) {
               jet_sf[i] = flavranges[1][3][1];
            } else if (jetpt[jet_ind[i]] > 140. && jetpt[jet_ind[i]] <= 200.) {
               jet_sf[i] = flavranges[1][4][1];
            } else {
               jet_sf[i] = 1.;
            }
          } else if (syst == "down"){
            if (jetpt[jet_ind[i]] > 30. && jetpt[jet_ind[i]] <= 50.) {
               jet_sf[i] = flavranges[1][0][2];
            } else if (jetpt[jet_ind[i]] > 50. && jetpt[jet_ind[i]] <= 70.) {
               jet_sf[i] = flavranges[1][1][2];
            } else if (jetpt[jet_ind[i]] > 70. && jetpt[jet_ind[i]] <= 100.) {
               jet_sf[i] = flavranges[1][2][2];
            } else if (jetpt[jet_ind[i]] > 100. && jetpt[jet_ind[i]] <= 140.) {
               jet_sf[i] = flavranges[1][3][2];
            } else if (jetpt[jet_ind[i]] > 140. && jetpt[jet_ind[i]] <= 200.) {
               jet_sf[i] = flavranges[1][4][2];
            } else {
               jet_sf[i] = 1.;
            }
          }
       } else {
          if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 200.) {
            jet_eff[i] = myfunc_expjoin(effmcranges[0][0], effmcranges[0][1], effmcranges[0][2], jetpt[jet_ind[i]]);
          } else {
            jet_eff[i] = myfunc_expjoin(effmcranges[0][0], effmcranges[0][1], effmcranges[0][2], 200.);
          }
          if (syst == "nom"){
            if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 1000.) {
               jet_sf[i] = myfunc_pol55(flavrangeslight[0][0], flavrangeslight[0][1], flavrangeslight[0][2], flavrangeslight[0][3], jetpt[jet_ind[i]]);
            } else {
               jet_sf[i] = 1.;
            }
          } else if (syst == "up") {
            if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 1000.) {
               jet_sf[i] = myfunc_pol55(flavrangeslight[0][0], flavrangeslight[0][1], flavrangeslight[0][2], flavrangeslight[0][3], jetpt[jet_ind[i]])*(1+myfunc_pol66(flavrangeslight[1][0], flavrangeslight[1][1], flavrangeslight[1][2], jetpt[jet_ind[i]]));
            } else {
               jet_sf[i] = 1.;
            }
          } else if (syst == "down") {
            if (jetpt[jet_ind[i]] > 20. && jetpt[jet_ind[i]] <= 1000.) {
               jet_sf[i] = myfunc_pol55(flavrangeslight[0][0], flavrangeslight[0][1], flavrangeslight[0][2], flavrangeslight[0][3], jetpt[jet_ind[i]])*(1-myfunc_pol66(flavrangeslight[1][0], flavrangeslight[1][1], flavrangeslight[1][2], jetpt[jet_ind[i]]));
            } else {
               jet_sf[i] = 1.;
            }
          }
       }
     }
     std::vector<float> term_sf;
     term_sf.push_back(1.);
     term_sf.push_back(1.);
     term_sf.push_back(1.);
     term_sf.push_back(1.);
     float weight_mc = 1.;
     std::vector<std::vector<float>> cut_ctag = {{0.270, 0.256},{0.269, 0.247},{0.520, 0.050},{0.289,0.267}};
     for (unsigned int i=0; i<2; ++i) {
       if (cvl[jet_ind[i]]>cut_ctag[yearT][0] && cvb[jet_ind[i]]>cut_ctag[yearT][1]) {
          term_sf[i] = jet_sf[i]*jet_eff[i]; 
          term_sf[i+2] = jet_eff[i]; 
       } else {
          term_sf[i] = 1-jet_sf[i]*jet_eff[i];
          term_sf[i+2] = 1-jet_eff[i];
       }
     }
     weight_mc = fabs(term_sf[2]*term_sf[3])>0 ? ((term_sf[0]*term_sf[1])/(term_sf[2]*term_sf[3])) : 1.;
     return weight_mc;
   };
""" % (dict_flav45,dict_mistag,dict_effmc))


####################################################################################################################################################################################

class ctag_joined():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('ctag_join_weight','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 0, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"nom")')
               df = df.Define('ctag_join_weight_up','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 0, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"up")')
               df = df.Define('ctag_join_weight_down','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 0, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"down")')
           elif self.year == "2016B":
               df = df.Define('ctag_join_weight','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 1, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"nom")')
               df = df.Define('ctag_join_weight_up','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 1, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"up")')
               df = df.Define('ctag_join_weight_down','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 1, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"down")')
           elif self.year == "2017":
               df = df.Define('ctag_join_weight','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 2, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"nom")')
               df = df.Define('ctag_join_weight_up','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 2, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"up")')
               df = df.Define('ctag_join_weight_down','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 2, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"down")')
           elif self.year == "2018":
               df = df.Define('ctag_join_weight','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 3, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"nom")')
               df = df.Define('ctag_join_weight_up','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 3, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"up")')
               df = df.Define('ctag_join_weight_down','joined_ctag_weight(JetXInd, Jet_hadronFlavour, Jet_pt_aux, 3, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB,"down")')

        variables = ['ctag_join_weight','ctag_join_weight_up','ctag_join_weight_down']

        branches = variables

        return df

def ctag_joinedRDF(**kwargs):
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
    return lambda: ctag_joined(**kwargs)

