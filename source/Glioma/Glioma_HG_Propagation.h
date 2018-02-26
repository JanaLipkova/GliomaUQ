//
//  Glioma_HG_Propagation.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 11/01/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef Glioma_HG_Propagation__
#define Glioma_HG_Propagation__

#pragma once
#include "Glioma_Types.h"


class Glioma_HG_Propagation: public Glioma
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
    bool                                    bVerbose;
    int                                     pID;

    
    static void _ic_SubjectBrainPropagation(Grid<W,B>& grid);
    static void _ic_PatientCasePropagation(Grid<W,B>& grid, int pID);

    static void	_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold );
    static void _selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize);
    static void _readInTumorPosition(vector<Real>& tumorIC);
    static int  _readInNumberOfSamples();

    void        _computeStatistics(const int nParallelGranularity, const int nSamples);
    void		_dump(int counter);
    void        _dumpBinaryOutput(Grid<W,B>& grid, int nSamples);

    
public:
    Glioma_HG_Propagation(int argc, const char ** argv);
    ~Glioma_HG_Propagation();
    void run();
    
};

#endif /* defined(Glioma_HG_Propagation__) */
