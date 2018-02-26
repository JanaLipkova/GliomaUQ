/*
 *  mainSTD_TST_L3.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

//pls visit http://chaton.inf.ethz.ch/wiki/index.php/MRAG:_Multiple_Resolution_Adapted_Grids for more info
//CVS version has this comment


#include <iostream>

#if _MRAG_OS == _MRAG_OS_APPLE
#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#include "GL/glew.h"
#include "GL/glut.h"
#define _USE_MATH_DEFINES
#endif

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_Haar.h"

//#include "MRAGcore/MRAG_STDTestL3_LevelsetReinitialization.h"
#include "MRAGcore/MRAG_STDTestL3_Diffusion.h"
//#include "MRAGcore/MRAG_STDTestL3_IO.h"


struct TestBlock: MRAG::Block<float, 10,10,10>
{
	typedef float ElementType;
	TestBlock(ElementType e = ElementType()): MRAG::Block<float, 10, 10, 10>(e){}
	
	//User needs to implement these functions for his blocks:
	//If you have enough disk-space for everything: this will do:
	void serialize(FILE* outputBinary)
	{
		//binary_stream.write(reinterpret_cast<char*> (this),sizeof(TestBlock));
		fwrite(reinterpret_cast<char*> (this),sizeof(TestBlock),1,outputBinary);
	}
	
	void deserialize(FILE* inputBinary )
	{
		//binary_stream.read(reinterpret_cast<char*> (this),sizeof(TestBlock));
		fread(reinterpret_cast<char*> (this),sizeof(TestBlock),1,inputBinary);

	}
	
};

int main (int argc, char **  argv) 
{
	
	MRAG::Environment::setup();
	
//	MRAG::STDTestL3_LevelsetReinitialization<MRAG::Wavelets_Interp4thOrder >::runTests(argc, argv, true);
	MRAG::STDTestL3_Diffusion<MRAG::Wavelets_Interp4thOrder >::runTests(argc, argv, false);
	//MRAG::MRAG_STDTestL3_IO<MRAG::Wavelets_Interp4thOrder, TestBlock >::runTests();
	
    return 0;
}

