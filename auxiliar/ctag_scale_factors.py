print('BTAG scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

ctag_2016_M_charm = "{1.000049,1.001356,1.050858,1.080369,1.003036,0.932686,1.160132}"
ctag_2016_M_bottom = "{1.00765764803,1.07777849339,1.01052439536,0.995782637473,0.950639166205}"
ctag_2016_M_light = "{1.19674,0.00011874,-1.89974e-07,-2.35069}"
ctag_2016_M = "{"+ctag_2016_M_charm+","+ctag_2016_M_bottom+","+ctag_2016_M_charm+"}"

ctag_2017_M_charm = "{0.898456,0.855116,0.876741,0.851636,0.855525}"
ctag_2017_M_bottom = "{1.0578398855,1.1336078783,1.12709741266,1.03942437216,1.0525262067}"
ctag_2017_M_light = "{1.09044,1.0214e-05,-1.01308e-07,0.519089}"
ctag_2017_M = "{"+ctag_2017_M_charm+","+ctag_2017_M_bottom+","+ctag_2017_M_charm+"}"

ctag_2018_M_charm = "{0.901666,0.927927,0.956623,0.989703,0.973638,0.945285,1.01199}"
ctag_2018_M_bottom = "{1.00390653345,0.920789734322,1.0486769655,1.06848197619,1.05542359778}"
ctag_2018_M_light = "{1.20542,-0.000267544,1.21665e-07,-1.91998}"
ctag_2018_M = "{"+ctag_2018_M_charm+","+ctag_2018_M_bottom+","+ctag_2018_M_charm+"}"

ctag_M_all= "{{0,"+ctag_2016_M+"},{1,"+ctag_2017_M+"},{2,"+ctag_2018_M+"}}"

pt_edges = "{{0,{25.,30.,50.,70.,100.,130.,160.,250.}},{1,{20.,30.,40.,50.,110.,210.}},{2,{25.,30.,50.,70.,100.,130.,160.,250.}}}"

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_pol2_inv(const float aa, const float bb, const float cc, const float dd, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx + dd/xx;
            return res;
   };
   auto ctag_scale_factors(Vint flav, Vfloat jetpt, const int yearT) {
     std::map <int, std::vector<std::vector<float>>> ctag_all = %s;
     std::map <int, std::vector<std::vector<std::vector<float>>>>::iterator it;
     auto flavranges = ctag_all[yearT];
     auto pt_terms = pt_all[yearT];
     vector<float> vb;
     std::vector<std::vector<float>> cut_ctag = {{0.181,0.228},{0.144,0.290},{0.1532,0.363}};
     for (unsigned int i=0; i<jetpt.size(); ++i) {
           if (fabs(flav[i]) == 4) {
              if (yearT == 1) {
                 if (jetpt[i] >= 20. && jetpt[i] < 30.) {
                    vb.push_back(flavranges[0][0]);
                 } else if (jetpt[i] >= 30. && jetpt[i] < 40.) {
                    vb.push_back(flavranges[0][1]);
                 } else if (jetpt[i] >= 40. && jetpt[i] < 50.) {
                    vb.push_back(flavranges[0][2]);
                 } else if (jetpt[i] >= 50. && jetpt[i] < 110.) {
                    vb.push_back(flavranges[0][3]);
                 } else if (jetpt[i] >= 110. && jetpt[i] < 210.) {
                    vb.push_back(flavranges[0][4]);
                 } else {
                    vb.push_back(1.);
              } else {
                 if (jetpt[i] >= 25. && jetpt[i] < 30.) {
                    vb.push_back(flavranges[0][0]);
                 } else if (jetpt[i] >= 30. && jetpt[i] < 50.) {
                    vb.push_back(flavranges[0][1]);
                 } else if (jetpt[i] >= 50. && jetpt[i] < 70.) {
                    vb.push_back(flavranges[0][2]);
                 } else if (jetpt[i] >= 70. && jetpt[i] < 100.) {
                    vb.push_back(flavranges[0][3]);
                 } else if (jetpt[i] >= 100. && jetpt[i] < 130.) {
                    vb.push_back(flavranges[0][4]);
                 } else if (jetpt[i] >= 130. && jetpt[i] < 160.) {
                    vb.push_back(flavranges[0][5]);
                 } else if (jetpt[i] >= 160. && jetpt[i] < 250.) {
                    vb.push_back(flavranges[0][6]);
                 } else {
                    vb.push_back(1.);
              }
           } else if (fabs(flav[i]) == 5) {
                 if (jetpt[i] >= 30. && jetpt[i] < 50.) {
                    vb.push_back(flavranges[1][0]);
                 } else if (jetpt[i] >= 50. && jetpt[i] < 70.) {
                    vb.push_back(flavranges[1][1]);
                 } else if (jetpt[i] >= 70. && jetpt[i] < 100.) {
                    vb.push_back(flavranges[1][2]);
                 } else if (jetpt[i] >= 100. && jetpt[i] < 140.) {
                    vb.push_back(flavranges[1][3]);
                 } else if (jetpt[i] >= 140. && jetpt[i] < 200.) {
                    vb.push_back(flavranges[1][4]);
                 } else {
                    vb.push_back(1.);
           } else {
               if (jetpt[i] > 0. && jetpt[i] < 1000.) {
                  vb.push_back(myfunc_pol2_inv(flavranges[2][0], flavranges[2][1], flavranges[2][2], flavranges[2][3], jetpt[i]));
               } else {
                  vb.push_back(1.);
               }
           }
     }
     return vb;
   }
""" % ctag_M_all )

####################################################################################################################################################################################

class ctag_sf():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('Aux_ctag_sf','ctag_scale_factors(Jet_partonFlavour, Jet_pt_aux, 0)')
           elif self.year == "2016B":
               df = df.Define('Aux_ctag_sf','ctag_scale_factors(Jet_partonFlavour, Jet_pt_aux, 0)')
           elif self.year == "2017":
               df = df.Define('Aux_ctag_sf','ctag_scale_factors(Jet_partonFlavour, Jet_pt_aux, 1)')
           elif self.year == "2018":
               df = df.Define('Aux_ctag_sf','ctag_scale_factors(Jet_partonFlavour, Jet_pt_aux, 2)')

        variables = ['Aux_ctag_sf']

        branches = variables

        return df

def ctag_sfRDF(**kwargs):
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
    return lambda: ctag_sf(**kwargs)

