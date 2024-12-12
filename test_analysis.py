###########################################
############## RDataFrame #################
###########################################

#### Para mas informacion visitar: https://root.cern/doc/master/classROOT_1_1RDataFrame.html

#### Importamos las librerias necesarias

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir
import json
import argparse

#### Activamos la ejecucion en paralelo

ROOT.EnableImplicitMT()

#### Usamos un fichero de prueba para leerlo y crear un dataset

testFile_path = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/ttbar_sl/cat_base/prod_test/data_0.root"
df = ROOT.RDataFrame("Events",testFile_path)

#### Con el dataframe creado podemos definir nuevas variables con la función Define(a,b) donde a y b son strings, siendo a el nombre de la nueva variable
#### y b es una línea de código C++ con el contenido de la variable

df = df.Define('var0','1.') ### Esto crea una nueva variable float que es 1. para todos los eventos de la ntupla
df = df.Define('var1','var0<0. ? var0*var0 : var0*var0*-1.') ### Se pueden usar las operaciones de C++

#### En vez de actualizar la ntupla df se puede crear otra distinta segun lo que interese

df_aux = df.Define('var2','var0*var1')

#### Para manipular las variables que ya están en la ntupla se puede saber sus nombres con el siguiente comando:

str = df.GetColumnNames()
## print(str)

#### Si queremos definir variables más sofisticadas podemos declarar código C++ previo y luego usarlo
#### Por ejemplo, si queremos el pT del muon mas energetico de cada evento con ciertas condiciones

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muon_maxpt(Vfloat muonpt, Vfloat muoneta, Vbool muonid) {
            float ptM{-10.};
            for (unsigned int i=0; i<muonpt.size(); ++i){
               if (muonpt[i]>ptM && fabs(muoneta[i])<2.4 && muonid[i]) {
                  ptM = muonpt[i];
       	       }
            }
            return ptM;
      };
""")

df = df.Define('max_muon_pt','muon_maxpt(Muon_pt, Muon_eta, Muon_tightId)')

#### Ademas podemos filtrar eventos

df_muon = df.Filter('max_muon_pt > 10.') # Nos quedamos con los eventos en los que hay al menos un muon con 10 GeV de momento

#### Podemos definir histogramas

histo1 = df_muon.Histo1D(("Hist_muon_pt", "Muon pt;p_T^{#mu} (GeV);N_{Events}", 100, 5, 55), "max_muon_pt") ## histograma de 100 bines de 5 a 55 GeV

#### IMPORTANTE ####

#### Hasta ahora nada de lo que hemos hecho colapsa el loop, todo son "lazy actions" para trigerear el loop hace falta o guardar una ntupla o un histograma
#### Pero es importante escribir primero todas las lazy actions y después ya las acciones que hagan que el loop sobre los datos comience 

myfile = TFile('hist_file_test.root', 'RECREATE')
histo1.Write()
myfile.Close()
