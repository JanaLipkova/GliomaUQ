/*
 *  mainSTD_TST_L4.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 11/19/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

//pls visit http://chaton.inf.ethz.ch/wiki/index.php/MRAG:_Multiple_Resolution_Adapted_Grids for more info
//CVS version has this comment

#include <iostream>





#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"

//#include "MRAGcore/MRAG_STDTestL3_LevelsetReinitialization.h"
#include "MRAGtests/MRAG_STDTestL4_LevelsetAdvection.h"

int main (int argc, char **  argv) 
{
	MRAG::Environment::setup();
	
	MRAG::STDTestL4_LevelsetAdvection<MRAG::Wavelets_Interp4thOrder >::runTests(argc, argv, true);
	
    return 0;
}

