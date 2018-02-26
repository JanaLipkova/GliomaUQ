//
//  Glioma_HG_PatientData.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 28/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef __GliomaBrutusXcode__Glioma_HG_PatientData__
#define __GliomaBrutusXcode__Glioma_HG_PatientData__

#pragma once
#include "Glioma_Types.h"

class Glioma_HG_PatientData: public Glioma
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
    
    static void _icPatient1(Grid<W,B>& grid);
    static void _icPatient7(Grid<W,B>& grid);
    static void _icPatient22(Grid<W,B>& grid);

    
    void        _setScalingFactor();
    void        _rescalePETsignal();
    void        _detectSegmentationsBC(BoundaryInfo* boundaryInfo, const int nParallelGranularity);
    void        _getAvPETatSegmentationsBC();
    
    
    void		_dump(int counter);
    void        _dumpOutput();
    
    
public:
    Glioma_HG_PatientData(int argc, const char ** argv);
    ~Glioma_HG_PatientData();
    void run();
};

#endif /* defined(__GliomaBrutusXcode__Glioma_HG_PatientData__) */
