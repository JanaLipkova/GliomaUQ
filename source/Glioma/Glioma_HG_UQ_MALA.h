//
//  Glioma_HG_UQ_MALA.h
//  
//
//  Created by Lipkova on 29/07/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
// ---------------------------------------------------
/*   Description:
     - solve Glioma model for UQ (same as Glioma_HG_UQ.c solver)
     - in addition also compute derivatives of model output w.r.t model parameters
     - derivatives are used for Metropolis Adjusted Langevin Algorithm (MALA)
 */

#ifndef ____Glioma_HG_UQ_MALA__
#define ____Glioma_HG_UQ_MALA__

#pragma once
#include "Glioma_Types.h"

class Glioma_HG_UQ_MALA: public Glioma
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
    
    static void _ic_SubjectBrain(Grid<W,B>& grid);
    static void _ic_PatientCase(Grid<W,B>& grid);
    
    static void	_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold );
    static void _selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize);
    static void _readInTumorPosition(vector<Real>& tumorIC);
    
    void        _reactionDiffusionMALAStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt);
    void		_dump(int counter);
    void        _dumpUQoutput(Grid<W,B>& grid);
    
    
public:
    Glioma_HG_UQ_MALA(int argc, const char ** argv);
    ~Glioma_HG_UQ_MALA();
    void run();
    
};

#endif /* defined(____Glioma_HG_UQ_MALA__) */



