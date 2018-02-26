/*
 *  mainSTD_TST_L1.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#include <iostream>
#undef min
#undef max
#include "MRAGcore/MRAGEnvironment.h"
//#include "MRAGcore/MRAG_STDTestL1_Generator.h"

#include "MRAG_STDTestL1_BoundaryInfo.h"

#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGWavelets_Haar.h"
#include "MRAGcore/MRAGHuffmanEncoder.h"

int main (int argc, char **  argv) 
{
	MRAG::Environment::setup();
	//MRAG_STDTestL1_Generator::run(argc, argv);
	//Encoder<
	/*MRAG::BitStream bs;
	bs.setup();
	bs.append_bits(23, 5);
	MRAG::BitStream bs2;
	bs2.setup();
	bs2.append_bits(5,3);
	bs.append_bits(bs2);*/
//	bs.append_bits(0xFF000000, 32);
//	bs.append_bits(0xFF000000, 32);
	//bs.printBits();
	//bs.get_bits(0,8);
	
//	bs.append_bits(5, 3);
	/*MRAG::HuffmanEncoder<char> encoder;
	vector<char> input;
	input.push_back('0');
	input.push_back('0');
	input.push_back('0');
	input.push_back('0');
	input.push_back('1');
	input.push_back('2');
	input.push_back('3');
	encoder.encode(input, 4);
	vector<char> output;
	encoder.decode(output);
	
	for(int i=0; i<output.size(); i++)
		printf("ouput %d: %c\n", i, output[i]);*/
	MRAG::MRAG_STDTestL1_BoundaryInfo< MRAG::Wavelets_AverageInterp5thOrder,  MRAG::Block< float,16,16,1> >::runTests(argc, argv, true);
	
    return 0;
}


