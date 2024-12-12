print('Trigger scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

#############################################
###### Branching fractions corrections ######
#############################################

## Funciones para los sf del trigger
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto trigger_sf_2016(const string& fFile, const string& fName, Vfloat pt, Vfloat eta, float quant) {
        TH2F *h;
        TFile my_file_sf(fFile.c_str(),"READ");
        h = (TH2F*)my_file_sf.Get(fName.c_str()); 
        vector<float> vb;
        for (unsigned int i=0; i<pt.size(); ++i) {
              if (pt[i] > quant && pt[i] <= 25.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,1)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,1)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,1)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,1)));
                  } else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 25. && pt[i] <= 30.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,2)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,2)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,2)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,2)));
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 30. && pt[i] <= 40.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,3)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,3)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,3)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,3)));
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 40. && pt[i] <= 50.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,4)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,4)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,4)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,4)));
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 50. && pt[i] <= 60.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,5)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,5)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,5)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,5)));
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 60. && pt[i] <= 120.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(1,6)));
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(h->GetBinContent(h->GetBin(2,6)));
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(h->GetBinContent(h->GetBin(3,6)));
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(h->GetBinContent(h->GetBin(4,6)));
                  }  else {
                     vb.push_back(1.);
                  }
              } else {
                    vb.push_back(1.);
              }
        }
        my_file_sf.Close();
        return vb;
      }
""")

####################################################################################################################################################################################

class trigger_mu_sf():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        trigger_sf_files_2016 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Runpre2016_UL_HIPM_SingleMuonTriggers.root"
        trigger_sf_files_2016B = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Runpost2016_UL_SingleMuonTriggers.root"
        trigger_sf_files_2017 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root"
        trigger_sf_files_2018 = "/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root"
        fName_2016 = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt"
        fName_2016B = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt"
        fName_2017 = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt"
        fName_2018 = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt"
        if self.isMC:
           if self.year == "2017":
               df = df.Define('trigger_sf_mu_aux','trigger_sf_2016("/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root","NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt", Muon_pt, Muon_eta, 29.0)')
           elif self.year == "2016":
               df = df.Define('trigger_sf_mu_aux','trigger_sf_2016("/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Runpre2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt", Muon_pt, Muon_eta, 26.0)')
           elif self.year == "2016B":
               df = df.Define('trigger_sf_mu_aux','trigger_sf_2016("/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Runpost2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt", Muon_pt, Muon_eta, 26.0)')
           elif self.year == "2018":
               df = df.Define('trigger_sf_mu_aux','trigger_sf_2016("/nfs/cms/vazqueze/ttbaranalisis/trigger_jsons/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt", Muon_pt, Muon_eta, 26.0)')

        variables = ['trigger_sf_mu_aux']

        branches = variables

        return df

def trigger_mu_sfRDF(**kwargs):
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
    return lambda: trigger_mu_sf(**kwargs)

