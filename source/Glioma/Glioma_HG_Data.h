//
//  Glioma_HG_Data.h
//  GliomaXcode
//
//  Created by Lipkova on 18/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef __GliomaXcode__Glioma_HG_Data__
#define __GliomaXcode__Glioma_HG_Data__

#pragma once
#include "Glioma_Types.h"


class Glioma_HG_Data: public Glioma
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
    
    Real           maxPhi, maxPET;
    Real           ucT1, ucT2;
    
    static void _icTumor(Grid<W,B>& grid);
    
    void        _setScalingFactor();
    void        _rescalePETsignal();
    void        _generateBinaryData();
    void        _generatePETdata();
    void        _detectSegmentationsBC(BoundaryInfo* boundaryInfo, const int nParallelGranularity);
    void        _getAvPETatSegmentationsBC();
    void        _computeError(Grid<W,B>& grid);
    void        _pickRandomPoints(int N);

    
    void		_dump(int counter);
    void        _dumpOutput();

    
public:
    Glioma_HG_Data(int argc, const char ** argv);
    ~Glioma_HG_Data();
    void run();
    
};
#endif /* defined(__GliomaXcode__Glioma_HG_Data__) */
