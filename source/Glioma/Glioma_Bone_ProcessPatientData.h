//
//  Glioma_Bone_ProcessPatientData.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 15/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __Glioma_Bone_ProcessPatientData__
#define __Glioma_Bone_ProcessPatientData__


#pragma once
#include "Glioma_Types.h"

class Glioma_Bone_ProcessPatientData: public Glioma
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
    Real                                    maxPET;
    
    static void _icPatientData(Grid<W,B>& grid);
    
    void        _getPETstatistic();
    void        _threasholdPETsignal();
    void        _normalisePETsignal();
    void		_dump(int counter);
    void        _dumpOutput();
    
public:
    Glioma_Bone_ProcessPatientData(int argc, const char ** argv);
    ~Glioma_Bone_ProcessPatientData();
    void run();
};


#endif /* defined(Glioma_HG_ProcessPatientData__) */
