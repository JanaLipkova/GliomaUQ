//
//  Glioma_TuringPatterns.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 07/11/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

#ifndef __Glioma_TuringPatterns__
#define __Glioma_TuringPatterns__

#pragma once
#include "Glioma_Types.h"


class Glioma_TuringPatterns: public Glioma
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
    int                                     ICtype;
    
    static void _ic_Schnackenberg(Grid<W,B>& grid);
    static void _ic_GrayScott(Grid<W,B>& grid);
    
    void        _TuringReactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, Real dt, Real Dpsi, Real Dphi, Real k1, Real k2, Real k3 = 0, Real k4 = 0);
    void		_dump(int counter);
    
public:
    Glioma_TuringPatterns(int argc, const char ** argv);
    ~Glioma_TuringPatterns();
    void run();
    
};

#endif /* defined(__Glioma_TuringPatterns__) */
