print('BTAG scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

eff_btag = "{{0,{{0.755221, 0.88767}, {0.12435, 0.38731},{0.00873,0.08478}}}, {1,{{0.78891, 0.90842}, {0.16213, 0.44549},{0.01253,0.11345}}}, {2,{{0.81954, 0.93174}, {0.16239, 0.49414},{0.01165,0.11567}}}, {3,{{0.82241, 0.93137}, {0.17038, 0.49652},{0.01212,0.11764}}}}"
#### Comentamos aquí cómo se usaria
#### """ % eff_btag ) ## al final
###     std::map <int, std::vector<std::vector<float>>> eff_btag_all = %s;
###     std::map <int, std::vector<std::vector<float>>>::iterator it; ## al principio
###      auto flavranges = eff_btag_all[yearT];
####  and then one can access the values with
#### flavranges[0][0]

bterm_2016 = "{{0,{{{8.96457e-01,-1.52849+02,6.63671e-03},{-3.26046,-5.67045,7.28429, -6.70161,3.67505}}, {{4.61519e-01, -2.43475e-03,1.87056e-05,-4.15439e-08},{1.68619e-01,-1.50050e-03, 1.17150e-05,-2.64113e-08}}, {{-5.04361e-02,3.56263e+01,1.23237e+04,-3.56235e+03,1.26185e-01},{-7.21077, -1.11074e+03,1.77001e-03,4.76605e-03,7.21835}}}},"
bterm_2016B = "{1,{{{9.01335e-01,-1.52945e+02,6.74166e-03},{6.00690,-5.58941e+01,2.30254e-02,-8.09513e-03,-5.20650}}, {{4.76616e-01, -2.54351e-03,1.90570e-05,-4.19702e-08},{1.76326e-01,-1.56466e-03,1.18777e-05,-2.64430e-08}}, {{-2.66795,-1.68069e+02,6.53654e-03,5.12816e-03,2.73776},{-2.70087,-2.30732e2,7.06126e-03,5.12817e-03,2.70903}}}},"
bterm_2017 = "{2, {{{9.33198e-01,-1.70363e+02,6.54397e-03},{5.14848,5.15611e-01,3.79144e-01,-7.55046e-01,-1.82755}}, {{5.67350e-01,-2.35569e-03,1.52796e-05,-2.86500e-08},{2.45945e-01,-2.71556e-03,2.11635e-05,-4.93887e-08}}, {{-7.52950e-02,3.40907e+01,1.13523e-01, -2.11876e-02,1.57648e-01},{ -6.93023e-03,3.55232e+01,1.44266e-01,-2.07173e-02,1.61089e-02}}}}, "
bterm_2018 = "{3,{{{-9.33433e-01,-1.86375e+02,-6.08555e-03},{6.53872,1.55142,1.20742, -3.72436,-1.44919}}, {{5.74930e-01,-2.56256e-03,1.74877e-05,-3.50985e-08},{2.59213e-01,-2.86323e-03,2.20375e-05,-5.10648e-08}}, {{-2.31672e+01,6.34169,2.51682e+02,-9.72341e+01,2.32381e+01},{-2.31984e+01,5.22932,2.72228e+02,-8.94228e+01,2.32068e+01}}}}}"

eff_btag_pt = bterm_2016+bterm_2016B+bterm_2017+bterm_2018

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_pol3(const float aa, const float bb, const float cc, const float dd, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx + dd*xx*xx*xx;
            return res;
   };
   auto myfunc_pol4(const float aa, const float bb, const float cc, const float dd, const float ee, const float xx) {
            float res;
            res = aa + bb*xx + cc*xx*xx + dd*xx*xx*xx + ee*xx*xx*xx*xx;
            return res;
   };
   auto myfunc_erf1(const float aa, const float bb, const float cc, const float xx) {
            float res;
       	    res	= aa*ROOT::Math::erf(cc*(xx-bb));
            return res;
   };
   auto myfunc_erf2(const float aa, const float bb, const float cc, const float dd, const float ee, const float xx) {
            float res;
       	    res	= aa*ROOT::Math::erf(cc*(xx-bb)/(1-dd*xx))+ee;
            return res;
   };
   auto total_btag_weight(Vint good, Vint flav, Vfloat btag, Vfloat jetpt, const int yearT, Vfloat sf_med_heavy, Vfloat sf_med_light, Vfloat sf_loo_heavy, Vfloat sf_loo_light) {
     std::map <int, std::vector<std::vector<std::vector<float>>>> eff_btag_all = %s;
     std::map <int, std::vector<std::vector<std::vector<float>>>>::iterator it;
     float weight_mc_complx = 1.;
     float weight_data_complx = 1.;
     float weight_mc_simpl = 1.;
     float weight_data_simpl = 1.;
     float eff_loose_evaluated;
     float eff_medium_evaluated;
     auto flavranges = eff_btag_all[yearT];
     std::vector<std::vector<float>> cut_btag = {{0.0614, 0.3093},{0.0480, 0.2489},{0.0532, 0.3040},{0.0490,0.2783}};
     for (unsigned int i=0; i<good.size(); ++i) {
           if (fabs(flav[good[i]]) == 5) {
                 if (jetpt[good[i]] > 200.) {
                     eff_loose_evaluated = myfunc_erf1(flavranges[0][0][0], flavranges[0][0][1], flavranges[0][0][2], 200.);
       	       	     eff_medium_evaluated = myfunc_erf2(flavranges[0][1][0], flavranges[0][1][1], flavranges[0][1][2], flavranges[0][1][3], flavranges[0][1][4], 200.);
                 } else {
                     eff_loose_evaluated = myfunc_erf1(flavranges[0][0][0], flavranges[0][0][1], flavranges[0][0][2], jetpt[good[i]]);
       	       	     eff_medium_evaluated = myfunc_erf2(flavranges[0][1][0], flavranges[0][1][1], flavranges[0][1][2], flavranges[0][1][3], flavranges[0][1][4], jetpt[good[i]]);
                 }
                 if (btag[good[i]] > cut_btag[yearT][1]) {
                     weight_mc_complx = weight_mc_complx*eff_medium_evaluated;
                     weight_data_complx = weight_data_complx*eff_medium_evaluated*sf_med_heavy[good[i]];
                     weight_mc_simpl = weight_mc_simpl*eff_medium_evaluated;
                     weight_data_simpl = weight_data_simpl*eff_medium_evaluated*sf_med_heavy[good[i]];
                 } else if (btag[good[i]] <= cut_btag[yearT][1] && btag[good[i]] > cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(eff_loose_evaluated - eff_medium_evaluated);
                     weight_data_complx = weight_data_complx*(eff_loose_evaluated*sf_loo_heavy[good[i]] - eff_medium_evaluated*sf_med_heavy[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_heavy[good[i]]);
                 } else if (btag[good[i]] <= cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(1 - eff_loose_evaluated);
                     weight_data_complx = weight_data_complx*( 1- eff_loose_evaluated*sf_loo_heavy[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_heavy[good[i]]);
                 }
           } else if (fabs(flav[good[i]]) == 4) {
                 if (jetpt[good[i]] > 200.) {
                    eff_loose_evaluated = myfunc_pol3(flavranges[1][0][0], flavranges[1][0][1], flavranges[1][0][2], flavranges[1][0][3], 200.);
    	       	    eff_medium_evaluated = myfunc_pol3(flavranges[1][1][0], flavranges[1][1][1], flavranges[1][1][2], flavranges[1][1][3], 200.);
                 } else {
                    eff_loose_evaluated = myfunc_pol3(flavranges[1][0][0], flavranges[1][0][1], flavranges[1][0][2], flavranges[1][0][3], jetpt[good[i]]);
    	       	    eff_medium_evaluated = myfunc_pol3(flavranges[1][1][0], flavranges[1][1][1], flavranges[1][1][2], flavranges[1][1][3], jetpt[good[i]]);
                 }
                 if (btag[good[i]] > cut_btag[yearT][1]) {
                     weight_mc_complx = weight_mc_complx*eff_medium_evaluated;
                     weight_data_complx = weight_data_complx*eff_medium_evaluated*sf_med_heavy[good[i]];
                     weight_mc_simpl = weight_mc_simpl*eff_medium_evaluated;
                     weight_data_simpl = weight_data_simpl*eff_medium_evaluated*sf_med_heavy[good[i]];
                 } else if (btag[good[i]] <= cut_btag[yearT][1] && btag[good[i]] > cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(eff_loose_evaluated - eff_medium_evaluated);
                     weight_data_complx = weight_data_complx*(eff_loose_evaluated*sf_loo_heavy[good[i]] - eff_medium_evaluated*sf_med_heavy[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_heavy[good[i]]);
                 } else if (btag[good[i]] <= cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(1 - eff_loose_evaluated);
                     weight_data_complx = weight_data_complx*( 1- eff_loose_evaluated*sf_loo_heavy[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_heavy[good[i]]);
                 }
           } else {
                 if (jetpt[good[i]] > 140.) {
   	       	    eff_loose_evaluated = myfunc_erf2(flavranges[2][0][0], flavranges[2][0][1], flavranges[2][0][2], flavranges[2][0][3], flavranges[2][0][4], 140.);
       	       	    eff_medium_evaluated = myfunc_erf2(flavranges[2][1][0], flavranges[2][1][1], flavranges[2][1][2], flavranges[2][1][3], flavranges[2][1][4], 140.);
                 } else {
   	       	    eff_loose_evaluated = myfunc_erf2(flavranges[2][0][0], flavranges[2][0][1], flavranges[2][0][2], flavranges[2][0][3], flavranges[2][0][4], jetpt[good[i]]);
       	       	    eff_medium_evaluated = myfunc_erf2(flavranges[2][1][0], flavranges[2][1][1], flavranges[2][1][2], flavranges[2][1][3], flavranges[2][1][4], jetpt[good[i]]);
                 }
                 if (btag[good[i]] > cut_btag[yearT][1]) {
                     weight_mc_complx = weight_mc_complx*eff_medium_evaluated;
                     weight_data_complx = weight_data_complx*eff_medium_evaluated*sf_med_light[good[i]];
                     weight_mc_simpl = weight_mc_simpl*eff_medium_evaluated;
                     weight_data_simpl = weight_data_simpl*eff_medium_evaluated*sf_med_light[good[i]];
                 } else if (btag[good[i]] <= cut_btag[yearT][1] && btag[good[i]] > cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(eff_loose_evaluated - eff_medium_evaluated);
                     weight_data_complx = weight_data_complx*(eff_loose_evaluated*sf_loo_light[good[i]] - eff_medium_evaluated*sf_med_light[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_light[good[i]]);
                 } else if (btag[good[i]] <= cut_btag[yearT][0]) {
                     weight_mc_complx = weight_mc_complx*(1 - eff_loose_evaluated);
                     weight_data_complx = weight_data_complx*( 1- eff_loose_evaluated*sf_loo_light[good[i]]);
                     weight_mc_simpl = weight_mc_simpl*(1 - eff_medium_evaluated);
                     weight_data_simpl = weight_data_simpl*(1 - eff_medium_evaluated*sf_med_light[good[i]]);
                 }
           }
     }
     vector<float> vb;
     vb.push_back(weight_mc_complx);
     vb.push_back(weight_data_complx);
     vb.push_back(weight_mc_simpl);
     vb.push_back(weight_data_simpl);
     return vb;
   }
""" % eff_btag_pt )

####################################################################################################################################################################################

class btag_weights_tot():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('Aux_btag_weight','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 0, btag_MED_sf, btag_MED_incl_sf, btag_LOO_sf, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 0, btag_MED_sf, btag_MED_incl_sf_up, btag_LOO_sf, btag_LOO_incl_sf_up)')
               df = df.Define('Aux_btag_weight_heavy_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 0, btag_MED_sf_up, btag_MED_incl_sf, btag_LOO_sf_up, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 0, btag_MED_sf, btag_MED_incl_sf_down, btag_LOO_sf, btag_LOO_incl_sf_down)')
               df = df.Define('Aux_btag_weight_heavy_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 0, btag_MED_sf_down, btag_MED_incl_sf, btag_LOO_sf_down, btag_LOO_incl_sf)')
           elif self.year == "2016B":
               df = df.Define('Aux_btag_weight','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 1, btag_MED_sf, btag_MED_incl_sf, btag_LOO_sf, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 1, btag_MED_sf, btag_MED_incl_sf_up, btag_LOO_sf, btag_LOO_incl_sf_up)')
               df = df.Define('Aux_btag_weight_heavy_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 1, btag_MED_sf_up, btag_MED_incl_sf, btag_LOO_sf_up, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 1, btag_MED_sf, btag_MED_incl_sf_down, btag_LOO_sf, btag_LOO_incl_sf_down)')
               df = df.Define('Aux_btag_weight_heavy_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 1, btag_MED_sf_down, btag_MED_incl_sf, btag_LOO_sf_down, btag_LOO_incl_sf)')
           elif self.year == "2017":
               df = df.Define('Aux_btag_weight','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 2, btag_MED_sf, btag_MED_incl_sf, btag_LOO_sf, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 2, btag_MED_sf, btag_MED_incl_sf_up, btag_LOO_sf, btag_LOO_incl_sf_up)')
               df = df.Define('Aux_btag_weight_heavy_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 2, btag_MED_sf_up, btag_MED_incl_sf, btag_LOO_sf_up, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 2, btag_MED_sf, btag_MED_incl_sf_down, btag_LOO_sf, btag_LOO_incl_sf_down)')
               df = df.Define('Aux_btag_weight_heavy_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 2, btag_MED_sf_down, btag_MED_incl_sf, btag_LOO_sf_down, btag_LOO_incl_sf)')
           elif self.year == "2018":
               df = df.Define('Aux_btag_weight','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 3, btag_MED_sf, btag_MED_incl_sf, btag_LOO_sf, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 3, btag_MED_sf, btag_MED_incl_sf_up, btag_LOO_sf, btag_LOO_incl_sf_up)')
               df = df.Define('Aux_btag_weight_heavy_up','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 3, btag_MED_sf_up, btag_MED_incl_sf, btag_LOO_sf_up, btag_LOO_incl_sf)')
               df = df.Define('Aux_btag_weight_light_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 3, btag_MED_sf, btag_MED_incl_sf_down, btag_LOO_sf, btag_LOO_incl_sf_down)')
               df = df.Define('Aux_btag_weight_heavy_down','total_btag_weight(JetGoodInd, Jet_hadronFlavour, Jet_btagDeepFlavB, Jet_pt_aux, 3, btag_MED_sf_down, btag_MED_incl_sf, btag_LOO_sf_down, btag_LOO_incl_sf)')

        variables = ['Aux_btag_weight', 'Aux_btag_weight_light_up', 'Aux_btag_weight_heavy_up', 'Aux_btag_weight_light_down', 'Aux_btag_weight_heavy_down']

        branches = variables

        return df

def btag_weights_totRDF(**kwargs):
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
    return lambda: btag_weights_tot(**kwargs)

