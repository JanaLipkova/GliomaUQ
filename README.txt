=========================================================

               Glioma UQ + MRAG solver

All rights reserverd: Jana Lipkova (jana.lipkova@tum.de)
==========================================================

Installation:
1) Install libraries in /lib/
2) Modify the Makefile (see notes below) 
3) If using brain anatomy in /source/Anatomy/ set up path to the anatomy inside the solver
4) run, see /source/makefile/Example

Below are notes for setting up the solver on Linux, Mac and cluster (SLURM)
Followed  by compilation notes and example (same for all enviroments)


************************
     Linux:
************************

---------------------------
1) Get the compiler
---------------------------

Check if you have gcc combiler (gcc47 or gcc49 are preferd but shoudl work with all versions starting from gcc4.2)
gcc --version

If not install with:
sudo apt-get install gcc

--------------------------------------
2) install libraries
--------------------------------------
- unpack libraries in /lib/ folder
- install with calling: make clean, make

---------------------------------------
3) set up the enviroment
---------------------------------------
Setup and make file with extension Kraken are for Linux enviroment:

i  ) set the path to the libraries see: source/make/setupGlioma_Kraken.sh
ii ) export path to libraries see make.kraken
iii) Makefile calls corresponding make.* depending on the hostname, modify it so it calls your make.*


**********************************************
 		MAC OX
**********************************************
--------------------------------------
1) get the compiler
--------------------------------------
Check your gcc compiler:
gcc --version

i) if gcc is not installed, install i.e. with MacPorts(https://www.macports.org/install.php)
ii) if you want to switch version of gcc compiler, swich it as follows:

sudo port select --list gcc
sudo port select --set gcc mp-gcc44
sudo port select --set gcc mp-gcc47

 --------------------------------------
2) install libraries
---------------------------------------
used libraries in /lib/ 
- unpack libraries in /lib/ folder
- install with calling: make clean, make


---------------------------------------
3) set up the enviroment 
---------------------------------------
Setup and make file with extension jana are for MAC OS:

i  ) set the path to the libraries see: source/make/setupGlioma_Jana.sh
ii ) export path to libraries see make.jana
iii) Makefile calls corresponding make.* depending on the hostname, modify it so it calls your make.*


**********************************************
          SLURM CLUSTER
**********************************************
--------------------------------------
1) get the compiler
--------------------------------------
module load gcc/4.8

--------------------------------------
2) install libraries
---------------------------------------
Follow the instructions for the Linux enviroment above


---------------------------------------
3) set up the enviroment
---------------------------------------
Setup and make file with extension lrz are for LRZ(SLURM) cluster, use mpp2 for MUC2 and mpp3 for MUC3

i  ) set the path to the libraries see: source/make/setupGlioma_mpp*.sh
ii ) export path to libraries see make.lrz
iii) Makefile calls corresponding make.* depending on the hostname, modify it so it calls your make.*

**********************************************
**********************************************
          Compilation + Example
**********************************************
**********************************************
--------------------------------------
1) Makefile + compilation
--------------------------------------
i) set up the enviroment:
source setupGlioma_*.sh

ii) The source/makefile/Makefile is used for compilation. It calls file make.* depending on the hostname. Modfy the Makefile so it calls your make.*

iii) compile with
make clean && make -j 4

-> creates executable called brain

iv) run as
./brain

--------------------------------------
2) Example
--------------------------------------
i)  copy the executable brain into /source/makefile/Example/
ii) attached script is for running:
    - on linux and mac computer see: runHGG_kraken.sh
    - on lrz see runHGG_lrz.sh (submit as sbatch runHGG_lrz.sh)


