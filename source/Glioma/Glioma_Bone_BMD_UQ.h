//
//  Glioma_Bone_BMD_UQ.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 09/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __GliomaBrutusXcode__Glioma_Bone_BMD_UQ__
#define __GliomaBrutusXcode__Glioma_Bone_BMD_UQ__

#pragma once
#include "Glioma_Types.h"


class Glioma_Bone_BMD_UQ: public Glioma
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
    int										numberOfIterations;
    double                                  whenToWrite;
    double                                  whenToWriteOffset;
    bool                                    isDone;
    bool                                    bAdaptivity;
    bool                                    bVerbose;
    
    static void _ic_Vertebra(Grid<W,B>& grid);
    static void _ic_Patient_Vertebrea(Grid<W,B>& grid);
    static void _readInTumorPosition(vector<Real>& tumorIC);
    
    
    void        _reactionDiffusionStep_BMD(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dscale, const Real rho, double dt);
    void		_dump(int counter);
    void        _dumpUQoutput(Grid<W,B>& grid);
    
    
public:
    Glioma_Bone_BMD_UQ(int argc, const char ** argv);
    ~Glioma_Bone_BMD_UQ();
    void run();
    
};

#endif /* defined(__GliomaBrutusXcode__Glioma_Bone_BMD_UQ__) */
