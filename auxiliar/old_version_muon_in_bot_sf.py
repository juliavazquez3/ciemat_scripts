print('Muon in bottom jet scale factor')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

#############################################
###### Branching fractions corrections ######
#############################################

#
gInterpreter.Declare("""

""")

gInterpreter.Declare("""
      #include <bitset>
      #include <string>
      #include <iostream>
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
      auto sf_muon_obt_iso(UInt_t nmu, Vfloat muon_iso){
            vector<float> vb;
            float sf_val = 1.;
            for (unsigned int i=0; i<nmu; ++i) {
                if (muon_iso[i] > 0. && muon_iso[i] < 0.5) {
                   sf_val =1.1160817;
                } else if (muon_iso[i] > 0.5 && muon_iso[i] < 1.) {
                   sf_val =1.0438291;
                } else if (muon_iso[i] > 1. && muon_iso[i] < 1.5) {
                   sf_val =0.97615421;
                } else if (muon_iso[i] > 1.5 && muon_iso[i] < 2.) {
                   sf_val =0.95700548;
                } else if (muon_iso[i] > 2. && muon_iso[i] < 2.5) {
                   sf_val =0.93035452;
                } else if (muon_iso[i] > 2.5 && muon_iso[i] < 3.) {
                   sf_val =0.90793874;
                } else if (muon_iso[i] > 3. && muon_iso[i] < 3.5) {
                   sf_val =0.90069424;
                } else if (muon_iso[i] > 3.5 && muon_iso[i] < 4.) {
                   sf_val =0.88897392;
                } else if (muon_iso[i] > 4. && muon_iso[i] < 4.5) {
                   sf_val =0.87461712;
                } else if (muon_iso[i] > 4.5 && muon_iso[i] < 5.) {
                   sf_val =0.86002409;
                } else if (muon_iso[i] > 5. && muon_iso[i] < 5.5) {
                   sf_val =0.85386769;
                } else if (muon_iso[i] > 5.5 && muon_iso[i] < 6.) {
                   sf_val =0.85858902;
                } else if (muon_iso[i] > 6. && muon_iso[i] < 6.5) {
                   sf_val =0.84793999;
                } else if (muon_iso[i] > 6.5 && muon_iso[i] < 7.) {
                   sf_val =0.84042996;
                } else if (muon_iso[i] > 7. && muon_iso[i] < 7.5) {
                   sf_val =0.83703991;
                } else if (muon_iso[i] > 7.5 && muon_iso[i] < 8.) {
                   sf_val =0.83748975;
                } else if (muon_iso[i] > 8. && muon_iso[i] < 9.) {
                   sf_val =0.84238776;
                } else if (muon_iso[i] > 9. && muon_iso[i] < 10.) {
                   sf_val =0.82480928;
                } else if (muon_iso[i] > 10. && muon_iso[i] < 12.5) {
                   sf_val =0.82456387;
                } else if (muon_iso[i] > 12.5 && muon_iso[i] < 15.) {
                   sf_val =0.81624343;
                } else if (muon_iso[i] > 15. && muon_iso[i] < 17.5) {
                   sf_val =0.81598700;
                } else if (muon_iso[i] > 17.5) {
                   sf_val =0.82430490;
                }
                vb.push_back(sf_val);
            }
            return vb;
      }
      auto third_grade_pol(const float a, const float b, const float c, const float d, const float x) {
            float vb;
            vb = a+b*x+c*x*x+d*x*x*x;
            return vb;
      }
      auto sf_muon_test_bots(const float muon_isoabs, const float muon_eta, const string syst) {
            vector<float> vb;
            float sf_val = 1.;
            if (syst == "bot1") {
              if (muon_isoabs<100.) {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.1421,-9.3994e-3,1.00973e-04,-3.92594e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.16587,-1.15339e-02,1.43742e-04,-6.0e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.11705,-9.08506e-03,9.62117e-05,-3.6e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.18271, -1.54494e-02,2.05526e-04,-9.0e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.11089, -1.14153e-02,1.14017e-04,-4.0e-07,muon_isoabs);
                  } else {
                     sf_val =1.;
                  }
              } else {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.1421,-9.3994e-3,1.00973e-04,-3.92594e-07,100.);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.16587,-1.15339e-02,1.43742e-04,-6.0e-07,100.);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.11705,-9.08506e-03,9.62117e-05,-3.6e-07,100.);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.18271, -1.54494e-02,2.05526e-04,-9.0e-07,100.);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.11089, -1.14153e-02,1.14017e-04,-4.0e-07,100.);
                  } else {
                     sf_val =1.;
                  }
              }
            } else {
              if (muon_isoabs<100.) {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.14600,-1.04757e-02,1.22985e-04,-4.93673e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.14660,-1.09570e-02,1.37406e-04,-5.82058e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.1284,-1.04855e-02,1.16279e-04,-4.0e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.17682,-1.38722e-02,1.70299e-04,-7.0e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.13303,-1.35750e-02,1.64187e-04,-6.5e-07,muon_isoabs);
                  } else {
                     sf_val =1.;
                  }
              } else {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.14600,-1.04757e-02,1.22985e-04,-4.93673e-07,100.);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.14660,-1.09570e-02,1.37406e-04,-5.82058e-07,100.);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.1284,-1.04855e-02,1.16279e-04,-4.0e-07,100.);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.17682,-1.38722e-02,1.70299e-04,-7.0e-07,100.);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.13303,-1.35750e-02,1.64187e-04,-6.5e-07,100.);
                  } else {
                     sf_val =1.;
                  }
              }
            }
            vb.push_back(sf_val);
            return vb;
      }
      auto sf_muon_obt_iso_abs(const float muon_isoabs, const float muon_eta, const string syst) {
            vector<float> vb;
            float sf_val = 1.;
            if (syst == "nom") {
              if (muon_isoabs<100.) {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.14408,-9.91869e-03,1.11677e-04,-4.42003e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.14855,-1.05787e-02,1.26e-04,-5.0e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.12236,-9.75158e-03,1.05912e-04,-3.79540e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.17323,-1.40088e-02,1.73e-04,-7.05705e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.11646,-1.19960e-02,1.3e-04,-4.78773e-07,muon_isoabs);
                  } else {
                     sf_val =1.;
                  }
              } else {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.14408,-9.91869e-03,1.11677e-04,-4.42003e-07,100.);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.14855,-1.05787e-02,1.26e-04,-5.0e-07,100.);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.12236,-9.75158e-03,1.05912e-04,-3.79540e-07,100.);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.17323,-1.40088e-02,1.73e-04,-7.05705e-07,100.);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.11646,-1.19960e-02,1.3e-04,-4.78773e-07,100.);
                  } else {
                     sf_val =1.;
                  }
              }
            } else if (syst == "up") {
              if (muon_isoabs<100.) {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.15802,-8.84418e-03,1.35648e-04,-2.84616e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.16225,-9.50077e-03,1.5044e-04,-3.37632e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.1395,-8.35395e-03,1.384e-04,-1.60025e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.18542,-1.29562e-02,1.98404e-04,-5.29722e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.13137,-1.05678e-02,1.67053e-04,-2.08976e-07,muon_isoabs);
                  } else {
                     sf_val =1.;
                  }
              } else {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.15802,-8.84418e-03,1.35648e-04,-2.84616e-07,100.);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.16225,-9.50077e-03,1.5044e-04,-3.37632e-07,100.);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.1395,-8.35395e-03,1.384e-04,-1.60025e-07,100.);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.18542,-1.29562e-02,1.98404e-04,-5.29722e-07,100.);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.13137,-1.05678e-02,1.67053e-04,-2.08976e-07,100.);
                  } else {
                     sf_val =1.;
                  } 
              }
            } else if (syst == "down") {
              if (muon_isoabs<100.) {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.13014,-1.09932e-02,8.77062e-05,-5.9939e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.13485,-1.16566e-02,1.0156e-04,-6.62368e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.10522,-1.11501e-02,7.34239e-05,-5.99055e-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.16104,-1.50614e-02,1.47596e-04,-8.81688-07,muon_isoabs);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.10155,-1.34242e-02,9.29466e-05,-7.4857-07,muon_isoabs);
                  } else {
                     sf_val =1.;
                  }
              } else {
                  if (fabs(muon_eta)<0.45) {
                     sf_val =third_grade_pol(1.13014,-1.09932e-02,8.77062e-05,-5.9939e-07,100.);
                  } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                     sf_val =third_grade_pol(1.13485,-1.16566e-02,1.0156e-04,-6.62368e-07,100.);
                  } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                     sf_val =third_grade_pol(1.10522,-1.11501e-02,7.34239e-05,-5.99055e-07,100.);
                  } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                     sf_val =third_grade_pol(1.16104,-1.50614e-02,1.47596e-04,-8.81688-07,100.);
                  } else if (fabs(muon_eta)>=1.8 && fabs(muon_eta)<2.4) {
                     sf_val =third_grade_pol(1.10155,-1.34242e-02,9.29466e-05,-7.4857-07,100.);
                  } else {
                     sf_val =1.;
                  } 
              }
            }
            vb.push_back(sf_val);
            return vb;
      }
      auto sf_muon_obt_z(UInt_t nmu, const float muon_z){
            vector<float> vb;
            float sf_val = 1.;
            if (muon_z > 0. && muon_z < 0.025) {
                   sf_val =1.0072096;
            } else if (muon_z > 0.025 && muon_z < 0.05) {
                   sf_val =0.90479816;
            } else if (muon_z > 0.05 && muon_z < 0.075) {
                   sf_val =0.88403459;
            } else if (muon_z > 0.075 && muon_z < 0.1) {
                   sf_val =0.88231193;
            } else if (muon_z > 0.1 && muon_z < 0.125) {
                   sf_val =0.88855841;
            } else if (muon_z > 0.125 && muon_z < 0.15) {
                   sf_val =0.89619713;
            } else if (muon_z > 0.15 && muon_z < 0.175) {
                   sf_val =0.90461392;
            } else if (muon_z > 0.175 && muon_z < 0.2) {
                   sf_val =0.91561947;
            } else if (muon_z > 0.2 && muon_z < 0.225) {
                   sf_val =0.92735799;
            } else if (muon_z > 0.225 && muon_z < 0.25) {
                   sf_val =0.94194931;
            } else if (muon_z > 0.25 && muon_z < 0.275) {
                   sf_val =0.96199020;
            } else if (muon_z > 0.275 && muon_z < 0.3) {
                   sf_val =0.96312640;
            } else if (muon_z > 0.3 && muon_z < 0.325) {
                   sf_val =0.94698137;
       	    } else if (muon_z > 0.325 && muon_z < 0.35) {
                   sf_val =0.96263222;
            } else if (muon_z > 0.35 && muon_z < 0.375) {
                   sf_val =0.97310410;
            } else if (muon_z > 0.375 && muon_z < 0.4) {
                   sf_val =0.95317820;
            } else if (muon_z > 0.4 && muon_z < 0.45) {
                   sf_val =0.93077223;
       	    } else if (muon_z > 0.45 && muon_z < 0.5) {
                   sf_val =0.90063700;
       	    } else if (muon_z > 0.5 && muon_z < 0.55) {
                   sf_val =0.84658799;
       	    } else if (muon_z > 0.55 && muon_z < 0.6) {
                   sf_val =0.79590258;
       	    } else if (muon_z > 0.6 && muon_z < 0.7) {
                   sf_val =0.69489125;
       	    } else if (muon_z > 0.7 && muon_z < 0.8) {
                   sf_val =0.62335481;
       	    } else if (muon_z > 0.8 && muon_z < 0.9) {
                   sf_val =0.54019142;
       	    } else if (muon_z > 0.9) {
                   sf_val =0.37773328;
       	    }
            vb.push_back(sf_val); 
            return vb;
      }
""")

####################################################################################################################################

class muon_frombot_sf():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")
        self.canal = kwargs.pop("canal")

    def run(self, df):

        if self.isMC:
           if self.canal == "selection_sl":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_jet_iso_abs, muon_jet_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso_abs_up','sf_muon_obt_iso_abs(muon_jet_iso_abs, muon_jet_eta, "up")')
              df = df.Define('muon_from_bot_sf_iso_abs_down','sf_muon_obt_iso_abs(muon_jet_iso_abs, muon_jet_eta, "down")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot1','sf_muon_test_bots(muon_jet_iso_abs, muon_jet_eta, "bot1")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot2','sf_muon_test_bots(muon_jet_iso_abs, muon_jet_eta, "bot2")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_jet_z)')
           if self.canal == "bot1":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_bot1_iso_abs, muon_bot1_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot1','sf_muon_test_bots(muon_bot1_iso_abs, muon_bot1_eta, "bot1")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot2','sf_muon_test_bots(muon_bot1_iso_abs, muon_bot1_eta, "bot2")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_bot1_z)')
           if self.canal == "bot2":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_bot2_iso_abs, muon_bot2_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot1','sf_muon_test_bots(muon_bot2_iso_abs, muon_bot2_eta, "bot1")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot2','sf_muon_test_bots(muon_bot2_iso_abs, muon_bot2_eta, "bot2")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_bot2_z)')
           if self.canal == "both":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot1','sf_muon_test_bots(muon_both_iso_abs, muon_both_eta, "bot1")')
              df = df.Define('muon_from_bot_sf_iso_abs_bot2','sf_muon_test_bots(muon_both_iso_abs, muon_both_eta, "bot2")')
              df = df.Define('muon_from_bot_sf_iso_abs_up','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "up")')
              df = df.Define('muon_from_bot_sf_iso_abs_down','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "down")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_both_z)')

        variables = ['muon_from_bot_sf_iso','muon_from_bot_sf_z','muon_from_bot_sf_iso_abs_bot1','muon_from_bot_sf_iso_abs_bot2']

        branches = variables

        return df

def muon_frombot_sfRDF(**kwargs):

    return lambda: muon_frombot_sf(**kwargs)

