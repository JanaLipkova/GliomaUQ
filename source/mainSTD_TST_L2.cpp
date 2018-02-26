/*
 *  mainSTD_TST_L2.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

//pls visit http://chaton.inf.ethz.ch/wiki/index.php/MRAG:_Multiple_Resolution_Adapted_Grids for more info
//CVS version has this comment


#include <iostream>

#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAG_STDTestL2_Generator.h"
#include "MRAGcore/MRAGBlock.h"
//#include "MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGWavelets_AverageInterp3rdOrder.h"

int main (int argc, char **  argv) 
{
	MRAG::Environment::setup();
	//MRAG::MRAG_STDTestL2_Generator::run(argc, argv);
	//MRAG::MRAG_STDTestL2_Refinement< MRAG::Wavelets_Interp4thOrder, MRAG::Block<float, 4,4,1> >::runTests();
	//MRAG::MRAG_STDTestL2_Refinement< MRAG::Wavelets_AverageInterp5thOrder, MRAG::Block<float, 8,8,1> >::runTests();
	MRAG::MRAG_STDTestL2_Refinement< MRAG::Wavelets_AverageInterp3rdOrder, MRAG::Block<float, 8,8,1> >::runTests();
    return 0;
}

