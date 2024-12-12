#####################################################################
## Scripts para correr el análisis en las máquinas gaes del CIEMAT ##
#####################################################################

El entorno que uso normalmente es:

source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

_________________________________________________________________________________________________________________________________________________

Ya se puede correr entonces con el comando:

python test_selection_analysis.py --process="M" --year="2016" --type="data" --presel="btagMM_chitest" --charmtag="no"  --wcs

en caso de querer desanclarlo de la terminal 

nohup python test_selection_analysis.py --process="M" --year="2016" --type="data" --presel="btagMM_chitest" --charmtag="no"  --wcs &> fich4.log &

Un ejemplo de MC sería

python test_selection_analysis.py --process="ww" --year="2017" --type="MC" --presel="btagMM_chitest" --charmtag="no"  --wcs

__________________________________________________________________________________________________________________________________________________

Los argumentos corresponden al tipo de proceso, agno, si es datos o MC, el tipo de preseleccion, el tipo de charmtag y si la clasificacion es por sabor de jets o 
como lo hacemos nosotros. Dentro del script estan todas las opciones posibles y se entiende mejor lo que hace si se busca dentro. 

El proceso del script es primero leer los datos, los abre como un RDataFrame, luego se definen las funciones C++ que necesitaremos, luego con df.Define() 
se aplican las funciones y con df.Filter() se hacen los cortes. Luego se procede a crear los histogramas y guardarlos.

Es importante mantener el orden para que RDataFrame funcione correctamente
