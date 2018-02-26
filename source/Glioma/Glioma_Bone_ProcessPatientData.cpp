//
//  Glioma_Bone_ProcessPatientData.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 15/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_Bone_ProcessPatientData.h"

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};


Glioma_Bone_ProcessPatientData::Glioma_Bone_ProcessPatientData(int argc, const char ** argv): parser(argc, argv)
{
    
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    printf("//////       High Grade Glioma UQ  - GENERATION OF PATIENT DATA        ///////\n");
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
    _icPatientData(*grid);
    
    _dump(0);
    isDone              = false;
}

Glioma_Bone_ProcessPatientData::~Glioma_Bone_ProcessPatientData()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* 1) Read in all patient data
 2) Rescaled into desired simulation domain
 3) Integrate time out of PET signals
 4) Restric PET signal only to T1 and T2 region */
void Glioma_Bone_ProcessPatientData:: _icPatientData(Grid<W,B>& grid)
{
#if !defined(Patient6) && !defined(Patient10) && !defined(Patient1)
#define Patient6
#define V2
#endif
    
#if defined(Patient1) && defined(V1)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V1/ \n");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V1/P1_CT_V1.dat");
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V1/P1_PET_V1.dat");
#endif
#if  defined(Patient1) && defined(V2)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V2/ \n");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V2/P1_CT_V2.dat");
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V2/P1_PET_V2.dat");
#endif
#if  defined(Patient1) && defined(V3)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V3/ \n");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V3/P1_CT_V3.dat");
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V3/P1_PET_V3.dat");
#endif
    
    
#if defined(Patient6) && defined(V2)
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V2/P6_CT_V2.dat");
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V2/P6_PET_V2.dat");
#endif

 
#if defined(Patient6) && defined(V3)
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V3/P6_PET_V3.dat");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V3/P6_CT_V3.dat");
#endif
    
#if  defined(Patient6) && defined(V4)
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V4/P6_PET_V4.dat");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V4/P6_CT_V4.dat");
#endif
    
#if  defined(Patient6) && defined(V5)
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V5/P6_PET_V5.dat");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V5/P6_CT_V5.dat");
#endif
    
    
#if defined(Patient10) && defined(V1)
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V1/P10_PET_V1.dat");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V1/P10_CT_V1.dat");
#endif
    
#if defined(Patient10) && defined(V2)
    MatrixD3D PET("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V2/P10_PET_V2.dat");
    MatrixD3D CT("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V2/P10_CT_V2.dat");
#endif
    
    int boneSizeX = (int) PET.getSizeX();
    int boneSizeY = (int) PET.getSizeY();
    int boneSizeZ = (int) PET.getSizeZ();
    printf("Data sizeX=%i, sizeY=%i, sizeZ=%i\n", boneSizeX, boneSizeY, boneSizeZ);

    int boneSizeMax = max(boneSizeX, max(boneSizeY,boneSizeZ));
    L    = boneSizeMax * 0.1 * 1.52344 * 0.5;   // voxel spacing 1.52344 mm, (*0.5 due to interpolation), convert from mm to cm  //
    
    
    double boneHx = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHy = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHz = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
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
                    
                    int mappedBoneX = (int)floor( (x[0]) / boneHx  );
                    int mappedBoneY = (int)floor( (x[1]) / boneHy  );
                    int mappedBoneZ = (int)floor( (x[2]) / boneHz  );
                    
                    // aspect ratio correction
                    mappedBoneX -= (int) ( (boneSizeMax - boneSizeX) * 0.5);
                    mappedBoneY -= (int) ( (boneSizeMax - boneSizeY) * 0.5);
                    mappedBoneZ -= (int) ( (boneSizeMax - boneSizeZ) * 0.5);
                    
                    if ( (mappedBoneX >= 0 && mappedBoneX < boneSizeX) && (mappedBoneY >= 0 && mappedBoneY < boneSizeY) && (mappedBoneZ >= 0 && mappedBoneZ < boneSizeZ) )
                    {
                        Real anatomy = CT(mappedBoneX,mappedBoneY,mappedBoneZ);
                        Real petSingal = PET(mappedBoneX,mappedBoneY,mappedBoneZ);
                        block(ix,iy,iz).phi  =  (anatomy >1) ? petSingal : 0.;   // restrict PET only to bone anaotmy
                        block(ix,iy,iz).tmp = (anatomy >1) ? anatomy : 0.;
                        
                    }
                }
        grid.getBlockCollection().release(info.blockID);

    }
    
    
}



void Glioma_Bone_ProcessPatientData::_getPETstatistic()
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
    
    printf("maxPET= %f \n", maxPET);
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

void Glioma_Bone_ProcessPatientData::_threasholdPETsignal()
{
    double uc = parser("-uc").asDouble(3.);
    printf("Threasholded with uc = %f \n",uc);
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).phi = (block(ix,iy,iz).phi > uc ) ? block(ix,iy,iz).phi : 0. ;
        
    }
}

void Glioma_Bone_ProcessPatientData::_normalisePETsignal()
{
    
    double scale = parser("-scale").asDouble(0.);
    double factor = max( (double)maxPET, scale);
    
    printf("Scale factor of PET :%f \n",factor);


    double iscale = 1./factor;
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).phi = (block(ix,iy,iz).phi > 0.1) ? block(ix,iy,iz).phi * iscale : 0.;
        
    }
}




#pragma mark DumpingOutput
void Glioma_Bone_ProcessPatientData:: _dump(int counter)
{
    
    if (parser("-vtk").asBool())
    {
        if(bVerbose) printf("dumping data \n");
        char filename[256];
        sprintf(filename,"%dD_Patient42_Data%04d",_DIM, counter);
        
        
        IO_VTKNative3D<W,B, 2,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_Bone_ProcessPatientData::_dumpOutput()
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    MatrixD3D PET(gpd,gpd,gpd);
    MatrixD3D CT(gpd,gpd,gpd);

    
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
                    {
                        PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).phi;
                        
                        Real bone = (block(ix,iy,iz).tmp > 50. ) ? 1. : 0.;
                        CT(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = bone;
                    }
        }
        else
        {
            for(int iz=0; iz<blockSize; iz++)
                for(int iy=0; iy<blockSize; iy++)
                    for(int ix=0; ix<blockSize; ix++)
                    {
                        PET(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                        CT(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                    }
            
            
        }
    }
    
    char filename[256];
    sprintf(filename,"PET_data.dat");
    PET.dump(filename);
    
    sprintf(filename,"Bone_data.dat");
    CT.dump(filename);
}


void Glioma_Bone_ProcessPatientData::run()
{
    
    _threasholdPETsignal();
    _getPETstatistic();
    _normalisePETsignal();
    _dump(1);
    _dumpOutput();
    
    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}
