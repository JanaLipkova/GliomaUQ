//
//  Glioma_TuringPatterns.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 07/11/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

#include "Glioma_TuringPatterns.h"

// need biger stencil for the refinment !!!
static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_TuringPatterns::Glioma_TuringPatterns(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////        Turing Reaction Diffusion Systems     ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    
    if(bVerbose) printf("suggested commands:\n");
    if(bVerbose) printf("mv test test_t%d_b%d_w%s\n", nThreads, blockSize, "w");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    bAdaptivity = parser("-adaptive").asBool();
    ICtype = parser("-IC").asInt();
    
    switch (ICtype)
    {
        case 0: _ic_Schnackenberg(*grid);
            break;
            
        case 1: _ic_GrayScott(*grid);
            break;
            
        default:  printf("Selected IC type is not supported");
            break;
    }
    
    if(parser("-dumpIC").asBool(1))
    _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
    
}

Glioma_TuringPatterns::~Glioma_TuringPatterns()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
//Literatue use U and V, here for convience we use G and W
// initialise with steady state solution with added white noise
void Glioma_TuringPatterns::_ic_Schnackenberg(Grid<W,B>& grid)
{
    
    Real noiseLevel1 = 0.01;
    Real noiseLevel2 = 0.01;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        srand48(i);  // set seed for rn generator
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    block(ix,iy,iz).psi= 2.   + noiseLevel1 * drand48();
                    block(ix,iy,iz).phi = 0.75 + noiseLevel2 * drand48();
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


void Glioma_TuringPatterns::_ic_GrayScott(Grid<W,B>& grid)
{
    Real noiseLevel1 = 0.005;
    Real noiseLevel2 = 0.0025;
    Real radius      = 0.2;
    Real center[3]   = {0.5, 0.5, 0.5};
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block        = grid.getBlockCollection()[info.blockID];
        
        srand48(i);  // set seed for rn generator
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real x[3];
                    info.pos(x, ix, iy, iz);
                    
                    const Real p[3]  = { x[0] - center[0], x[1] - center[1], x[2] - center[2]};
                    const Real dist  = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                    
                    if ( dist <= radius){
                        block(ix,iy,iz).psi = 0.5 + noiseLevel1 * drand48();
                        block(ix,iy,iz).phi = 0.25 + noiseLevel2 * drand48();
                    }
                    else{
                        block(ix,iy,iz).psi = 1.;
                        block(ix,iy,iz).phi = 0.;
                    }
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}





#pragma mark ReactionDiffusion

void  Glioma_TuringPatterns:: _TuringReactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, Real dt, Real Dpsi, Real Dphi, Real k1, Real k2, Real k3, Real k4)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    if(ICtype==0){
        SchnackenbergReactionDiffusionOperator  rhs_sch(Dpsi,Dphi,k1,k2,k3,k4);
        TuringTimeUpdate                        update_sch(dt);
        
        blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs_sch);
        BlockProcessing::process(vInfo, collecton, update_sch, nParallelGranularity);
    }
    else{
        GrayScottReactionDiffusionOperator      rhs_gs(Dpsi,Dphi,k1,k2);
        TuringTimeUpdate                        update_gs(dt);
        
        std::cout<<"In GS, Dpsi = "<<Dpsi <<" Dphi = "<< Dphi << std::endl;
        
        blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs_gs);
        BlockProcessing::process(vInfo, collecton, update_gs, nParallelGranularity);
    }


}



#pragma mark DumpingOutput
void Glioma_TuringPatterns:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"Turing_Data%04d", counter);
        
        IO_VTKNative3D<W,B, 2,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
    
}


void Glioma_TuringPatterns::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    Real Dpsi = parser("-Dpsi").asDouble();
    Real Dphi = parser("-Dphi").asDouble();
    Real Tend = parser("-Tend").asDouble(1000.);
    
    // Parameters, set k1 - k, k2=F
    Real k1 = parser("-k1").asDouble(1.);
    Real k2 = parser("-k2").asDouble(1.);
    Real k3 = parser("-k3").asDouble(2.);
    Real k4 = parser("-k4").asDouble(3.);
    
    
    double t			= 0.0;
    int iCounter        = 1;
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.99 * h*h / ( 2.* _DIM * max(Dpsi, Dphi) );
    if(bVerbose)  printf("Dg=%e, Dw=%e, dt= %f, h=%f\n", Dpsi, Dphi, dt,h);
    
    
    
    while (t <= Tend)
    {
        if(ICtype==0)
            _TuringReactionDiffusionStep(boundaryInfo, nParallelGranularity, Dpsi, Dphi, dt, k1, k2, k3, k4);
        else
            _TuringReactionDiffusionStep(boundaryInfo, nParallelGranularity, Dpsi, Dphi, dt, k1, k2);
        
        
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( numberOfIterations >= ((double)(whenToWrite)) )
        {
            if(bAdaptivity)
            {
                Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            _dump(iCounter);
            iCounter++;
            whenToWrite = whenToWrite + whenToWriteOffset;
            if(bVerbose) printf("Dumping data at time t=%f\n", t);
            
        }
    }
    

    // Refine final state & dump for UQ Likelihood
    if(bAdaptivity)
        Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);

    _dump(iCounter);
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
