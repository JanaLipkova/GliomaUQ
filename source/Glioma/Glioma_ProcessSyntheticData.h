//
//  Glioma_ProcessSyntheticData.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 12/11/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef Glioma_ProcessSyntheticData__
#define Glioma_ProcessSyntheticData__

#pragma once
#include "Glioma_Types.h"


class Glioma_ProcessSyntheticData: public Glioma
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
    void        _normalizePETsignal();
    void        _addNoise();
    void        _dumpOutput();
    void		_dump(int counter);
    void        _readInTumorPosition(vector<Real>& tumorIC );
    void        _dumpSubBrainPoints(Grid<W,B>& grid );


    
public:
    Glioma_ProcessSyntheticData(int argc, const char ** argv);
    ~Glioma_ProcessSyntheticData();
    void run();
    
};


#endif /* defined(Glioma_ProcessSyntheticData__) */
