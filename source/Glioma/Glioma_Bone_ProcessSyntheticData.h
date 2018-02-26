//
//  Glioma_Bone_ProcessSyntheticData.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 10/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __Glioma_Bone_ProcessSyntheticData__
#define __Glioma_Bone_ProcessSyntheticData__

#pragma once
#include "Glioma_Types.h"


class Glioma_Bone_ProcessSyntheticData: public Glioma
{
private:
    Grid<W, B>								* grid;
    BlockProcessing							blockProcessing;
    Refiner_SpaceExtension					*refiner;
    Compressor								*compressor;
    BlockFWT<W, B, RD_Projector_Wavelets>	blockfwt;			// refinment based on single channel
    SpaceTimeSorter							stSorter;
    Profiler								profiler;
    ArgumentParser							parser;
    IO_VTK< W, B, RD_Projector_VTK >		vtk;
    BlockLab< B >							lab;
    bool                                    isDone;
    bool                                    bAdaptivity;
    bool                                    bVerbose;
    
    Real                                    maxPhi, maxPET;
    
    static void _ic(Grid<W,B>& grid);   // read in tumour and cut values below threahsold
    
    void        _getPETstatistic();
    void        _addNoise();
    void        _normalizePETsignal();
    void        _dumpUQoutput();
    void		_dump(int counter);
    
    
    
public:
    Glioma_Bone_ProcessSyntheticData(int argc, const char ** argv);
    ~Glioma_Bone_ProcessSyntheticData();
    void run();
    
};


#endif /* defined(__Glioma_Bone_ProcessSyntheticData__) */
