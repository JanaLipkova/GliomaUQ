//
//  Glioma_HG_Visualizations.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 21/01/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef Glioma_HG_Visualizations__
#define Glioma_HG_Visualizations__


#pragma once
#include "Glioma_Types.h"


class Glioma_HG_Visualizations: public Glioma
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
    void		_dump(int counter);
    void        _dumpBinaryOutput(Grid<W,B>& grid);

    
public:
    Glioma_HG_Visualizations(int argc, const char ** argv);
    ~Glioma_HG_Visualizations();
    void run();
    
};

#endif /* defined(Glioma_HG_Visualizations__) */
