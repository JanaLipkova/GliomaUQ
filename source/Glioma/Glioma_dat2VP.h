//
//  Glioma_dat2VP.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 19/03/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef __Glioma_dat2VP__
#define __Glioma_dat2VP__



#pragma once
#include "Glioma_Types.h"
#include "MRAGio/DumpScalarToVP.h"


class Glioma_dat2VP: public Glioma
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

    
    static void _ic(Grid<W,B>& grid, int pID);
    void        _normalisePET(Grid<W,B>& grid);

    void		_dumpVTK(int counter);
    void		_dumpVP(int counter);
    void        _copyFields(Grid<W,B>& grid,int filedID);
    void        _combineFields(Grid<W,B>& grid);


    
    
public:
    Glioma_dat2VP(int argc, const char ** argv);
    ~Glioma_dat2VP();
    void run();
    
};

#endif /* defined(__Glioma_dat2VP__) */
