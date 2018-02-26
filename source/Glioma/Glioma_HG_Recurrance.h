//
//  Glioma_HG_Recurrance.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 29/11/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __GliomaBrutusXcode__Glioma_HG_Recurrance__
#define __GliomaBrutusXcode__Glioma_HG_Recurrance__

#pragma once
#include "Glioma_Types.h"

class Glioma_HG_Recurrance: public Glioma
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
    
    
    // IC
    static void _ic_SubjectBrain(Grid<W,B>& grid);
    static void _ic_AtlasBrain(Grid<W,B>& grid);
    static void _ic_PatientCase(Grid<W,B>& grid);
    
    
    // IC helpers
    static void _readInTumorPosition(vector<Real>& tumorIC);
    static void	_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold );
    static void _selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize);
    
    
    // Operators
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt);
    
    // other helpers
    void        _removeTumour();
    void		_dump(int counter);
    void        _dumpTumourBinary(int counter);

    
public:
    Glioma_HG_Recurrance(int argc, const char ** argv);
    ~Glioma_HG_Recurrance();
    void run();
    
};

#endif /* defined(__GliomaBrutusXcode__Glioma_HG_Recurrance__) */
