print('CTAG scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join
import pandas as pd

import json
import argparse

## mistag1 is 2016 then 2016B,2017 and 2018
## first is flavour 0, second 4 and third 5
## values correspond to first coeff [0], second [1] and third [3] in equation [0]+[1]*exp([2]*x)

mistag1 = "{{0,{{4.13182e-02,-2.37595e-02,-2.36221e-02},{3.97975e-01,-3.75733e-01,-3.30244e-02},{4.09319e-01,-2.57558e-01,-1.51208e-02}}},"
mistag2 = "{1,{{4.17324e-02,-2.35894e-02,-2.32517e-02},{3.99713e-01,-3.73863e-01,-3.28992e-02},{4.12717e-01,-2.55883e-01,-1.51056e-02}}},"
mistag3 = "{2,{{2.70404e-01,-2.59587e-01,-7.16768e-05},{3.88679e-01,-4.04733e-01,-3.96979e-02},{3.48276e-01,-2.47078e-01,-2.50963e-02}}},"
mistag4 = "{3,{{2.81567e-01,-2.48132e-01,-2.02721e-04},{4.72195e-01,-4.11363e-01,-3.98227e-02},{4.33088e-01,-2.61989e-01,-1.94343e-02}}}}"

dict_mistag = mistag1+mistag2+mistag3+mistag4

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_exp11(const float aa, const float bb, const float cc, const float xx) {
            float res;
            res = aa + bb*std::exp(cc*xx);
            return res;
   };
   auto total_veto_weight(const int jet_ind, Vint flav, Vfloat jetpt, const int yearT, const float ctagsf, const float cmistagsf) {
     std::map <int, std::vector<std::vector<float>>> dict_sf = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it2;
     float eff_mc = 1.;
     float weight_mc = 1.;
     auto flavranges = dict_sf[yearT];
     if (fabs(flav[jet_ind]) != 4 && fabs(flav[jet_ind]) != 5) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[0][0], flavranges[0][1], flavranges[0][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[0][0], flavranges[0][1], flavranges[0][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-cmistagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     } else if (fabs(flav[jet_ind]) == 4) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[1][0], flavranges[1][1], flavranges[1][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[1][0], flavranges[1][1], flavranges[1][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-ctagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     } else if (fabs(flav[jet_ind]) == 5) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[2][0], flavranges[2][1], flavranges[2][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[2][0], flavranges[2][1], flavranges[2][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-ctagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     }
     return weight_mc;
   }
""" % dict_mistag)

####################################################################################################################################################################################

class vetoctag_weights_tot():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 0, ctag_weight, cmistag_weight)')
           elif self.year == "2016B":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 1, ctag_weight, cmistag_weight)')
           elif self.year == "2017":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 2, ctag_weight, cmistag_weight)')
           elif self.year == "2018":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 3, ctag_weight, cmistag_weight)')

        variables = ['vetoctag_weight']

        branches = variables

        return df

def vetoctag_weights_totRDF(**kwargs):
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
    return lambda: vetoctag_weights_tot(**kwargs)

