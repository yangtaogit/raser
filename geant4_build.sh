#!/usr/bin/env bash

# Main driver to build GEANT4 as user inside raser 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2021-06-28 Mon 12:47]

    
echo "-----build geant4-----" 
mkdir -p $HOME/geant4/src
mkdir -p $HOME/geant4/build 
mkdir -p $HOME/geant4/install
mkdir -p $HOME/geant4/data 
# wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p02.tar.gz
# tar xzf geant4.10.07.p02.tar.gz -C $HOME/geant4/src 
cd $HOME/geant4/build 
cmake -DCMAKE_INSTALL_PREFIX=$HOME/geant4/install  -DGEANT4_INSTALL_DATA=ON -DGEANT4_INSTALL_DATADIR=$HOME/geant4/data -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_USE_QT=ON  -DGEANT4_BUILD_TLS_MODEL=global-dynamic  -DGEANT4_USE_PYTHON=ON ../src/geant4.10.07.p02  
make -j24 
make install 

