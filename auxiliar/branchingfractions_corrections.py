print('D hadrons corrections')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

#############################################
###### Branching fractions corrections ######
#############################################

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;

// Function to print the
// index of an element
// 1: Full hadronic 2: One W hadronic, the other leptonic (SEMI) 3: Full leptonic
      auto getIndex(UInt_t nGenPart,Vint v, int K){
            int ind = 2;
            while (v[ind]!=K) {
                ++ind;
            }
            return ind;
      }
      auto some_features_sl(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetmuonind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi){
            vector<int> vb;
            bool cond_D = false;
            bool cond1 = false;
            bool isitok = false;
            int muon_ind = -1;
            int mother_ind = -1;
            float mother_pt = 0;
            float mother_eta = 0;
            int muon_char = 0;
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetmuonind[0]],gen_eta[mother[i]],jet_phi[jetmuonind[0]],gen_phi[mother[i]]) < 0.4;
                cond_D = (fabs(pdg[mother[i]])== 431 || fabs(pdg[mother[i]]) == 411 || fabs(pdg[mother[i]]) == 421 || fabs(pdg[mother[i]]) == 4122);
                if (fabs(pdg[i])==13 && cond_D && cond1) {
                          isitok = true;
                          muon_ind = i;
                          mother_ind = mother[i];
                          muon_char = pdg[i]/fabs(pdg[i]);
                          mother_pt = gen_pt[mother[i]];
                          mother_eta = gen_eta[mother[i]];
                }
                vb.push_back(isitok);
                vb.push_back(muon_ind);
                vb.push_back(mother_ind);
                vb.push_back(muon_char);
            }
            return isitok;
      }
      auto some_features_sl2(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetmuonind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi){
            vector<float> vb;
            bool cond_D = false;
            bool cond1 = false;
            bool isitok = false;
            int muon_ind = -1;
            int mother_ind = -1;
            float mother_pt = 0;
            float mother_eta = 0;
            int muon_char = 0;
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetmuonind[0]],gen_eta[mother[i]],jet_phi[jetmuonind[0]],gen_phi[mother[i]]) < 0.4;
                cond_D = (fabs(pdg[mother[i]])== 431 || fabs(pdg[mother[i]]) == 411 || fabs(pdg[mother[i]]) == 421 || fabs(pdg[mother[i]]) == 4122);
                if (fabs(pdg[i])==13 && cond_D && cond1) {
                          isitok = true;
                          muon_ind = i;
                          mother_ind = mother[i];
                          muon_char = pdg[i]/fabs(pdg[i]);
                          mother_pt = gen_pt[mother[i]];
                          mother_eta = gen_eta[mother[i]];
                }
                vb.push_back(mother_pt);
                vb.push_back(mother_eta);
            }
            return vb;
      }
      auto branchingfracrat_factor_sl(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, const bool some_feat, Vfloat some2, Vint jetmuonind, Vfloat jet_eta, Vfloat jet_phi) {
            float br_weight = 1.;
            float frag_weight = 0.992127;
            vector<float> vb;
            for (unsigned int i=0; i<nGenPart; ++i) {
              if (ROOT::VecOps::DeltaR(jet_eta[jetmuonind[0]],gen_eta[i],jet_phi[jetmuonind[0]],gen_phi[i]) < 0.4) {
                 if (fabs(pdg[i]) == 411) {
                        if (some_feat && gen_pt[i] == some2[0] && gen_eta[i] == some2[1]) {
                             br_weight = 1.0283; // D+>muX    
                        } else {
                             br_weight = 0.995356;
                        }
                        frag_weight = 0.8289; // c->D+
                 } else if (fabs(pdg[i]) == 421) {
                        if (some_feat && gen_pt[i] == some2[0] && gen_eta[i] == some2[1]) {
                             br_weight = 0.9927; // D0>muX
                        } else {
                             br_weight = 1.000411; 
                        }
                        frag_weight = 1.0778; // c->D0
                 } else if (fabs(pdg[i]) == 431) {
                        if (some_feat && gen_pt[i] == some2[0] && gen_eta[i] == some2[1]) {
                             br_weight = 0.8233; // Ds>muX
                        } else {
                             br_weight = 1.011687;
                        }
                        frag_weight = 0.8279; // c->Ds
                 } else if (fabs(pdg[i]) == 4122) {
                        if (some_feat && gen_pt[i] == some2[0] && gen_eta[i] == some2[1]) {
                             br_weight = 0.7629; //
                        } else {
                             br_weight = 1.009146;
                        }
                        frag_weight = 1.7435; // c->lambdac
                 }
              }
            }
            //if (!some_feat[0]) frag_weight = 0.992127;
            vb.push_back(br_weight);
            vb.push_back(frag_weight);
            return vb;
      };
      auto branchingfracrat_factor_sv1(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetsvind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi) {
            float br_weight = 1.;
            float frag_weight = 1.;
            bool cond1 = false;
            bool cond2 = false;
            bool cond3 = false;
            vector<bool> vb;
            bool tiened0 = false;
            bool d0tienepi0 = false;
            bool d0tienek = false;
            bool d0tienepip = false;
            bool d0tienepip2 = false;
            int d0signpip = 0;
            int d0igpip = -1;
            int d0igpip2  = -1;
            float d0pt{0.};
            float d0eta{0.};
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetsvind[0]],gen_eta[mother[i]],jet_phi[jetsvind[0]],gen_phi[mother[i]]) < 0.4;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==421 && cond1) {
                       tiened0 = true;
                       d0tienepip = true;
                       d0signpip = pdg[i]/fabs(pdg[i]);
                       d0igpip = i;
                       d0pt = gen_pt[mother[i]];
                       d0eta = gen_eta[mother[i]];
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = d0signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==421 && d0tienepip && i != d0igpip && cond2 && cond3) {
                       d0tienepip2 = true;
                       d0igpip2 = i;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==321 && fabs(pdg[mother[i]])==421 && d0tienepip2 && cond3) {
                       d0tienek = true;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==421 && d0tienek && cond3 && i != d0igpip && i != d0igpip2) {
                       d0tienepi0 = true;
                }
            }
            vb.push_back(tiened0);
            vb.push_back(d0pt);
            vb.push_back(d0eta);
            vb.push_back(d0tienepip);
            vb.push_back(d0tienepip2);
            vb.push_back(d0tienek);
            vb.push_back(d0tienepi0);
            return vb;
      };
      auto branchingfracrat_factor_sv2(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetsvind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi) {
            float br_weight = 1.;
            float frag_weight = 1.;
            bool cond1 = false;
            bool cond2 = false;
            bool cond3 = false;
            vector<bool> vb;
            bool tienedp = false;
            bool tienek = false;
            bool tieneele = false;
            bool tienepi0 = false;
            bool tienepip = false;
            bool tienepip2 = false;
            int signpip = 0;
            int igpip = -1;
            int igpip2  = -1;
            float d0pt{0.};
            float d0eta{0.};
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetsvind[0]],gen_eta[mother[i]],jet_phi[jetsvind[0]],gen_phi[mother[i]]) < 0.4;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==411 && cond1) {
                       tienedp = true;
                       tienepip = true;
                       signpip = pdg[i]/fabs(pdg[i]);
                       igpip = i;
                       d0pt = gen_pt[mother[i]];
                       d0eta = gen_eta[mother[i]];
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==411 && tienepip && i != igpip && cond2 && cond3) {
                       tienepip2 = true;
                       igpip2 = i;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==321 && fabs(pdg[mother[i]])==411 && tienepip && cond3) {
                       tienek = true;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==11 && fabs(pdg[mother[i]])==411 && tienepip && cond3) {
                       tieneele = true;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==111 && fabs(pdg[mother[i]])==411 && tienepip && tienepip2 && tienek && cond3) {
                       tienepi0 = true;
                }
            }
            vb.push_back(tienedp);
            vb.push_back(d0pt);
            vb.push_back(d0eta);
            vb.push_back(tienepip);
            vb.push_back(tienepip2);
            vb.push_back(tienek);
            vb.push_back(tieneele);
            vb.push_back(tienepi0);
            return vb;
      };
      auto branchingfracrat_factor_sv3(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetsvind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi) {
            float br_weight = 1.;
            float frag_weight = 1.;
            bool cond1 = false;
            bool cond2 = false;
            bool cond3 = false;
            vector<bool> vb;
            bool tieneds = false;
            bool tienepi0 = false;
            bool tienek = false;
            bool tienepip = false;
            bool tienepip2 = false;
            int signpip = 0;
            int igpip = -1;
            int igpip2  = -1;
            float d0pt{0.};
            float d0eta{0.};
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetsvind[0]],gen_eta[mother[i]],jet_phi[jetsvind[0]],gen_phi[mother[i]]) < 0.4;
                if (fabs(pdg[i])==321 && fabs(pdg[mother[i]])==431 && cond1) {
                       tieneds = true;
                       tienepip = true;
                       signpip = pdg[i]/fabs(pdg[i]);
                       igpip = i;
                       d0pt = gen_pt[mother[i]];
                       d0eta = gen_eta[mother[i]];
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==321 && fabs(pdg[mother[i]])==431 && tienepip && i != igpip && cond2 && cond3) {
                       tienepip2 = true;
                       igpip2 = i;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = pdg[mother[i]]*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==431 && tienepip2 && cond3 && cond2) {
                       tienek = true;
                }
            }
            vb.push_back(tieneds);
            vb.push_back(d0pt);
            vb.push_back(d0eta);
            vb.push_back(tienepip);
            vb.push_back(tienepip2);
            vb.push_back(tienek);
            return vb;
      };
      auto branchingfracrat_factor_sv4(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vint jetsvind, Vfloat muon_charge, Vfloat jet_eta, Vfloat jet_phi) {
            float br_weight = 1.;
            float frag_weight = 1.;
            bool cond1 = false;
            bool cond2 = false;
            bool cond3 = false;
            vector<bool> vb;
            bool tienelambda = false;
            bool tienepip = false;
            bool tienek = false;
            bool tieneproton = false;
            int signpip = 0;
            int igpip = -1;
            int igpip2  = -1;
            float d0pt{0.};
            float d0eta{0.};
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond1 = ROOT::VecOps::DeltaR(jet_eta[jetsvind[0]],gen_eta[mother[i]],jet_phi[jetsvind[0]],gen_phi[mother[i]]) < 0.4;
                if (fabs(pdg[i])==211 && fabs(pdg[mother[i]])==4122 && cond1) {
                       tienelambda = true;
                       tienepip = true;
                       signpip = pdg[i]/fabs(pdg[i]);
                       igpip = i;
                       d0pt = gen_pt[mother[i]];
                       d0eta = gen_eta[mother[i]];
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond2 = signpip*pdg[i] > 0;
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==2212 && fabs(pdg[mother[i]])==4122 && tienepip && i != igpip && cond2 && cond3) {
                       tieneproton = true;
                       igpip2 = i;
                }
            }
            for (unsigned int i=0; i<nGenPart; ++i) {
                cond3 = gen_eta[mother[i]] == d0eta && gen_pt[mother[i]] == d0pt;
                if (fabs(pdg[i])==321 && fabs(pdg[mother[i]])==4122 && tienepip && tieneproton && cond3) {
                       tienek = true;
                }
            }
            vb.push_back(tienelambda);
            vb.push_back(d0pt);
            vb.push_back(d0eta);
            vb.push_back(tienepip);
            vb.push_back(tieneproton);
            vb.push_back(tienek);
            return vb;
      };
      auto factors_obtention_sv(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vbool aux1, Vbool aux2, Vbool aux3, Vbool aux4){
            vector<float> vb;
            bool cond_D = false;
            bool cond1 = false;
            bool isitok = false;
            int muon_ind = -1;
            int mother_ind = -1;
            int muon_char = 0;
            float Br_weight = 1.;
            float Frag_weight = 1.;
            for (unsigned int i=0; i<nGenPart; ++i) {
                 if (fabs(pdg[i]) == 421 && aux1[0]) {
                         if (aux1[6]) {  
                                Br_weight =  1.11216;
                         } else {
                                Br_weight = 0.99868;
                         }
	                 Frag_weight = 1.0778; // c->D0     
                 } else if (fabs(pdg[i]) == 411 && aux2[0]) {
                         if (aux2[7]) {
                                Br_weight = 5.2; // D+>Kpipipi0
                         } else if (aux2[5] && aux2[4]) {  
                                Br_weight = 1.019565; // D+>Kpipi
                         } else if (aux2[6]) {
                                Br_weight = 14.888889;
                         } else {
                                Br_weight = 0.871384;
                         }
                         Frag_weight =  0.8289; // c->D+     
                 } else if (fabs(pdg[i]) == 431 && aux3[0]) {
                         if (aux3[5]) {  
                                Br_weight = 4.9; // Ds>KKpi 5.2 in PDG and 1.1 in MC
                         } else { 
                                Br_weight = 0.9350; // 0.935847 para W -> ELECTRON and  0.934959  para W MUON weight on Ds events that are not Kpipi or Kpipipi0 
                         }
                         Frag_weight =  0.8279; // c->Ds     
                 } else if (fabs(pdg[i]) == 4122 && aux4[0]) {
                         if (aux4[5]) {
                                Br_weight = 2.94; // lambdac >p K pi 6.28 in PDG and 2.1 in MC
                         } else {
                                Br_weight = 0.955469; //
                         }
                         Frag_weight = 1.7435; // c->lambdac                 }
                 }
            }
            vb.push_back(Br_weight);
            vb.push_back(Frag_weight);
            return vb;
      }
""")

#####################################################################################################################################333


class branchingfractions_SL_cor():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):

        if self.isMC:
           df = df.Define('Br_cor_sl_aux1','some_features_sl(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetXInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sl_aux2','some_features_sl2(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetXInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sl_aux3','branchingfracrat_factor_sl(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, Br_cor_sl_aux1, Br_cor_sl_aux2, JetXInd, Jet_eta, Jet_phi)')
           df = df.Define('Frag_weight_sl','Br_cor_sl_aux3[1]')
           df = df.Define('Br_weight_sl','Br_cor_sl_aux3[0]')

        variables = ['Br_cor_sl_aux1','Br_cor_sl_aux2','Br_cor_sl_aux3','Frag_weight_sl','Br_weight_sl']

        branches = variables

        return df

def branchingfractions_SL_corRDF(**kwargs):
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
    return lambda: branchingfractions_SL_cor(**kwargs)

class branchingfractions_SV_cor():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):

        if self.isMC:
           df = df.Define('Br_cor_sv_aux1','branchingfracrat_factor_sv1(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetSVInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sv_aux2','branchingfracrat_factor_sv2(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetSVInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sv_aux3','branchingfracrat_factor_sv3(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetSVInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sv_aux4','branchingfracrat_factor_sv4(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, JetSVInd, Muon_charge, Jet_eta, Jet_phi)')
           df = df.Define('Br_cor_sv_aux5','factors_obtention_sv(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, Br_cor_sv_aux1, Br_cor_sv_aux2, Br_cor_sv_aux3, Br_cor_sv_aux4)')
           df = df.Define('Frag_weight_sv','Br_cor_sv_aux5[1]')
           df = df.Define('Br_weight_sv','Br_cor_sv_aux5[0]')

        variables = ['Br_cor_sv_aux5','Frag_weight_sv','Br_weight_sv']

        branches = variables

        return df

def branchingfractions_SV_corRDF(**kwargs):
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
    return lambda: branchingfractions_SV_cor(**kwargs)

