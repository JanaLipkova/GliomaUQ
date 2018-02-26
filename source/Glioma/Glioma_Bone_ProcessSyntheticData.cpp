//
//  Glioma_Bone_ProcessSyntheticData.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 10/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_Bone_ProcessSyntheticData.h"


static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};


Glioma_Bone_ProcessSyntheticData::Glioma_Bone_ProcessSyntheticData(int argc, const char ** argv): parser(argc, argv)
{
    
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    printf("//////       High Grade Glioma UQ  - GENERATION OF SYNTHETIC DATA        ///////\n");
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    
    printf("suggested commands:\n");
    printf("mv test test_t%d_b%d_w%s\n", nThreads, blockSize, "w");
    printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    bVerbose = parser("-verbose").asBool();
    
    _ic(*grid);
    _dump(0);
    
    isDone              = false;
}

Glioma_Bone_ProcessSyntheticData::~Glioma_Bone_ProcessSyntheticData()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
// 1) Read in tumour concentration
// 2) threashold it
// 3) add noise
void Glioma_Bone_ProcessSyntheticData:: _ic(Grid<W,B>& grid)
{
    
    MatrixD3D tumor("UQ_data.dat");
    int tumor_size = (int) (tumor.getSizeX() * tumor.getSizeY() * tumor.getSizeZ());
    printf("tumor_size=%i\n",tumor_size);
    
    Real uc = 0.2;  // threashold

    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    Real petSignal = tumor(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2]);
                    block(ix,iy,iz).phi  = ( petSignal >= uc ) ? petSignal : 0.;
                    
                }
        
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



// compute normalizing factor & center of mass
void Glioma_Bone_ProcessSyntheticData::_getPETstatistic()
{
    maxPET = 0.;
    float mass = 0.;
    float cx = 0.;
    float cy = 0.;
    float cz = 0.;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    maxPET = max(maxPET, block(ix,iy,iz).phi);
                    mass += block(ix,iy,iz).phi;
                    
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    cx += x[0] * block(ix,iy,iz).phi;
                    cy += x[1] * block(ix,iy,iz).phi;
                    cz += x[2] * block(ix,iy,iz).phi;
                }
        
    }
    
    cx = cx / mass;
    cy = cy / mass;
    cz = cz / mass;
    
    printf("maxPET=%f \n", maxPET);
    printf("Center of mass: cx=%f, cy=%f, cz=%f \n", cx,cy,cz);
    
    Real mx,my,mz;
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    if(maxPET == block(ix,iy,iz).phi)
                    {
                        double x[3];
                        info.pos(x, ix, iy, iz);
                        mx = x[0];
                        my = x[1];
                        mz = x[2];
                    }
                }
        
    }
    
    printf("Max of PET: mx=%f, my=%f, mz=%f \n", mx,my,mz);
    
    FILE * pFile;
    float buffer[3] = { cx , cy , cz };
    pFile = fopen ("TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);
    fclose (pFile);
    
}

// scale PET signal is [0,1]
void Glioma_Bone_ProcessSyntheticData::_normalizePETsignal()
{
    double iscale = 1./maxPET;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).phi = block(ix,iy,iz).phi * iscale;
    }
}

// add noise
void Glioma_Bone_ProcessSyntheticData::_addNoise()
{
    
    double xi_1, xi_2;
    int n = 0;
    double alpha = 0.02;  // 5% noise
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        srand48(i);
        
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    // Box-Muller
                    if(  (n % 2) == 0 )
                    {
                        double u1 = drand48();
                        double u2 = drand48();
                        
                        xi_1 = sqrt( - 2.0*log( u1 ) ) * cos( 2.0*M_PI*u2 );
                        xi_2 = sqrt( - 2.0*log( u1 ) ) * sin( 2.0*M_PI*u2 );
                    }
                    
                    
                    if( block(ix,iy,iz).phi > 0.  )
                    {
                        block(ix,iy,iz).phi = max( 0., block(ix,iy,iz).phi + alpha * xi_1);
                        xi_1 = xi_2;
                        n++;
                    }
                    
                    
                }
    }
}


#pragma mark DumpingOutput
void Glioma_Bone_ProcessSyntheticData:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%dD_Patient42_Data%04d",_DIM, counter);
        
        IO_VTKNative3D<W,B, 1,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_Bone_ProcessSyntheticData::_dumpUQoutput()
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    MatrixD3D PET(gpd,gpd,gpd);
    
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        double h = info.h[0];
        
        
        if(h == hf)
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                        PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).phi;
            
        }
        else
        {
            for(int iz=0; iz<blockSize; iz++)
                for(int iy=0; iy<blockSize; iy++)
                    for(int ix=0; ix<blockSize; ix++)
                        PET(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
            
        }
    }
    
    char filename[256];
    sprintf(filename,"PET_data.dat");
    PET.dump(filename);
    
}





void Glioma_Bone_ProcessSyntheticData::run()
{
    _dump(0);

    _getPETstatistic();
    _addNoise();
    _getPETstatistic();
    _normalizePETsignal();
    _getPETstatistic();
    
    _dump(1);
    _dumpUQoutput();
    
        
    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}

