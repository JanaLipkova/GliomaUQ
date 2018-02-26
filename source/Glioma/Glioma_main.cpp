/*
 *  mainGlioma.cpp
 *  Untitled
 *
 *  Created by Lipkova on 9/14/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


// Brain stuffs
#include "Glioma_HG_UQ.h"
#include "Glioma_ProcessSyntheticData.h"
#include "Glioma_HG_ProcessPatientData.h"
#include "Glioma_HG_Propagation.h"
#include "Glioma_HG_Visualizations.h"
#include "Glioma_HG_UQ_MALA.h"
//#include "Glioma_HG_Recurrance.h"

//Bone stuffs
//#include "Glioma_Bone_UQ.h"
//#include "Glioma_Bone_BMD_UQ.h"
//#include "Glioma_Bone_ProcessSyntheticData.h"
//#include "Glioma_Bone_ProcessPatientData.h"

// Necrosis
#include "Glioma_Necrosis.h"

// Radiotherpy
#include "Glioma_HG_AddUniformMargin.h"

// VP Visualization
#include "Glioma_dat2VP.h"

// Turing patterns
#include "Glioma_TuringPatterns.h"

#include "Test.h"

#include <iostream>
#include <xmmintrin.h>

using namespace std;
using namespace MRAG;

int main(int argc,const char ** argv)
{
	std::cout << std::endl << "MRAG Launched" << std::endl << std::endl;
	
	ArgumentParser parser(argc, argv);
	Environment::setup(max(1, parser("-nthreads").asInt()));

	cout <<" Precission: sizeof(Real)="<< sizeof(Real) << endl;
	
	Glioma * s = NULL;
	
	printf("INPUT IS %s\n", parser("-anatomy").asString().data());
	
      if(parser("-anatomy").asString() == "HGG_UQ")
      s = new Glioma_HG_UQ(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "HGG_UQ_MALA")
          s = new Glioma_HG_UQ_MALA(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "recurrence")
          s = new Glioma_HG_Recurrance(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "syntheticData")
      s = new Glioma_ProcessSyntheticData(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "patientData")
      s = new Glioma_HG_ProcessPatientData(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "propagation")
          s = new Glioma_HG_Propagation(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "visualization")
          s = new Glioma_HG_Visualizations(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "simpleBone")
          s = new Glioma_Bone_UQ(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "bmdBone")
          s = new Glioma_Bone_BMD_UQ(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "syntheticBoneData")
          s = new Glioma_Bone_ProcessSyntheticData(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "patientBoneData")
          s = new Glioma_Bone_ProcessPatientData(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "necrosis")
          s = new Glioma_Necrosis(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "addMargin")
          s = new Glioma_HG_AddUniformMargin(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "VPvisualisation")
          s = new Glioma_dat2VP(argc, (const char **)argv);
      else if(parser("-anatomy").asString() == "Turing")
          s = new Glioma_TuringPatterns(argc, (const char **)argv);
      else
      s = new Test(argc, (const char **)argv);

	
	tbb::tick_count t1,t0;
	{
		t0=tbb::tick_count::now();
		s->run();
		t1=tbb::tick_count::now();
	}
	
	printf("we spent: %2.2f \n",(t1-t0).seconds());	
	delete s;
	
	std::cout << std::endl << "MRAG Terminated" << std::endl;
}



