//
//  Glioma_HG_ProcessPatientData.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 13/11/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef Glioma_HG_ProcessPatientData__
#define Glioma_HG_ProcessPatientData__

#pragma once
#include "Glioma_Types.h"

class Glioma_HG_ProcessPatientData: public Glioma
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
    int                                     pID;
    
    static void _icPatientData(Grid<W,B>& grid, int pID);
    
    void        _getPETstatistic();
    void        _normalisePETsignal();
    void		_dump(int counter);
    void        _dumpOutput();
    void        _dumpSubBrainPoints(Grid<W,B>& grid );
    void        _readInTumorPosition(vector<Real>& tumorIC );

public:
    Glioma_HG_ProcessPatientData(int argc, const char ** argv);
    ~Glioma_HG_ProcessPatientData();
    void run();
};


#endif /* defined(Glioma_HG_ProcessPatientData__) */
