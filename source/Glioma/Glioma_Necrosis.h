//
//  Glioma_Necrosis.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 15/03/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __Glioma_Necrosis__
#define __Glioma_Necrosis__


#pragma once
#include "Glioma_Types.h"


class Glioma_Necrosis: public Glioma
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
    
    
    static void _ic_SubjectBrain(Grid<W,B>& grid);
    static void _ic_PatientCase(Grid<W,B>& grid);
    
    static void	_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold );
    static void _selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize);
    static void _readInTumorPosition(vector<Real>& tumorIC);
    
    void        _reactionDiffusionNecrosisStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, const Real gamma, double dt);
    void		_dump(int counter);
    
    
public:
    Glioma_Necrosis(int argc, const char ** argv);
    ~Glioma_Necrosis();
    void run();
    
};


#endif /* defined(__Glioma_Necrosis__) */
