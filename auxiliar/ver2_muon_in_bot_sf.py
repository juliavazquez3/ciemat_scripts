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

histfile = TFile.Open("higgssearch/fromJF/ratios_muons_botjet/hists_confidenceintervals/confidence_intervals_etabins_all.root","READ")

h1 = histfile.Get("etabinone")
h2 = histfile.Get("etabintwo")
h3 = histfile.Get("etabinthree")
h4 = histfile.Get("etabinfour")
h5 = histfile.Get("etabinfive")

hist_d = [h1,h2,h3,h4,h5]
limi = h1.GetNbinsX()

isoabs_edges = [2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,37.5,40.,45.,50.,60.,70.,80.,100.]
etabins = [0,1,2,3,4]

dict_to_cpp = "{"
for k,keta in enumerate(etabins):
  dict_to_cpp += "{%s, {" % etabins[k]
  for i,isoed in enumerate(isoabs_edges):
     dict_to_cpp += "{%s, %s, %s} %s" % (hist_d[i].GetBinContent(i+1), (hist_d[i].GetBinContent(i+1))+(hist_d[i].GetBinError(i+1)), (hist_d[i].GetBinContent(i+1))-(hist_d[i].GetBinError(i+1)),
      ("," if i< len(isoabs_edges)-1 else "}"))
  dict_to_cpp +="}%s" % (", " if k < len(etabins) - 1 else "}")

gInterpreter.Declare("""
      #include <bitset>
      #include <string>
      #include <iostream>
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
      std::map<int, std::vector<std::vector<float>>> m = %s;
      vector<float> arr = {0.,2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,37.5,40.,45.,50.,60.,70.,80.};
      auto sf_muon_obt_iso_abs_ver2(const float muon_isoabs, const float muon_eta, const string syst) {
            int binx = 0;
            for (unsigned int i=0; i<arr.size(); ++i) {
               if (muon_isoabs >= arr[i]) {
                  binx = i+1;
               }
            }
            #include <vector>
            vector<float> vb;
            float sf_val = 1.;
            std::vector<std::vector<float>>> mvals;
            if (fabs(muon_eta)<0.45) {
               mvals = m[0];
            } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
               mvals = m[1];
            } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
               mvals = m[2];
            } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
               mvals = m[3];
            } else {
               mvals = m[4];
            }
            if (syst == "nom") {
                sf_val = mvals[binx][0];
            } else if (syst == "up") {
                sf_val = mvals[binx][1];
            } else if (syst == "down") {
                sf_val = mvals[binx][2];
            }
            vb.push_back(sf_val);
            return vb;
      }
""" %dict_to_cpp)


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
      TFile *file1 = new TFile("/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/hists_confidenceintervals/confidence_intervals_etabins_all.root");
      TH1F * h1 = (TH1F*)file1->Get("etabinone");
      TH1F * h2 = (TH1F*)file1->Get("etabintwo");
      TH1F * h3 = (TH1F*)file1->Get("etabinthree");
      TH1F * h4 = (TH1F*)file1->Get("etabinfour");
      TH1F * h5 = (TH1F*)file1->Get("etabinfive");
      TH1F * hcl;
      auto sf_muon_obt_iso_abs(const float muon_isoabs, const float muon_eta, const string syst) {
            #include <vector>
            vector<float> vb;
            float sf_val = 1.;
            if (fabs(muon_eta)<0.45) {
                hcl = (TH1F*)(h1->Clone("hcl"));
            } else if (fabs(muon_eta)>=0.45 && fabs(muon_eta)<0.9) {
                hcl = (TH1F*)(h2->Clone("hcl"));
            } else if (fabs(muon_eta)>=0.9 && fabs(muon_eta)<1.2) {
                hcl = (TH1F*)(h3->Clone("hcl"));
            } else if (fabs(muon_eta)>=1.2 && fabs(muon_eta)<1.8) {
                hcl = (TH1F*)(h4->Clone("hcl"));
            } else {
                hcl = (TH1F*)(h5->Clone("hcl"));
            }
            int limi = hcl->GetNbinsX();
            int binx = std::max(1, std::min(hcl->GetNbinsX(), hcl->GetXaxis()->FindBin(muon_isoabs)));
            if (binx<1) binx = 1;
            if (binx>22) binx = 22;
            if (syst == "nom") {
                 if (muon_isoabs<100.) {
                    sf_val = hcl->GetBinContent(binx);
                 } else {
                    sf_val = hcl->GetBinContent(limi);
                 }
            } else if (syst == "up") { 
                 if (muon_isoabs<100.) {
                    sf_val = (hcl->GetBinContent(binx))+(hcl->GetBinError(binx));
                 } else {
                    sf_val = (hcl->GetBinContent(limi))+(hcl->GetBinError(limi));
                 }
            } else if (syst == "down") {
                 if (muon_isoabs<100.) {
                    sf_val = (hcl->GetBinContent(binx))-(hcl->GetBinError(binx));
                 } else {
                    sf_val = (hcl->GetBinContent(limi))-(hcl->GetBinError(limi));
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
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_jet_z)')
           if self.canal == "bot1":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_bot1_iso_abs, muon_bot1_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_bot1_z)')
           if self.canal == "bot2":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_bot2_iso_abs, muon_bot2_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_bot2_z)')
           if self.canal == "both":
              df = df.Define('muon_from_bot_sf_iso_abs','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "nom")')
              df = df.Define('muon_from_bot_sf_iso_abs_up','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "up")')
              df = df.Define('muon_from_bot_sf_iso_abs_down','sf_muon_obt_iso_abs(muon_both_iso_abs, muon_both_eta, "down")')
              df = df.Define('muon_from_bot_sf_iso','sf_muon_obt_iso(nMuon, Muon_pfRelIso04_all)')
              df = df.Define('muon_from_bot_sf_z','sf_muon_obt_z(nMuon, muon_both_z)')

        variables = ['muon_from_bot_sf_iso','muon_from_bot_sf_z']

        branches = variables

        return df

def muon_frombot_sfRDF(**kwargs):

    return lambda: muon_frombot_sf(**kwargs)

