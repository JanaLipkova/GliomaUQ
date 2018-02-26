//
//  Glioma_HG_PatientData.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 28/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_HG_PatientData.h"

static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_HG_PatientData::Glioma_HG_PatientData(int argc, const char ** argv): parser(argc, argv)
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
    
    bAdaptivity = parser("-adaptive").asBool();
    bVerbose = parser("-verbose").asBool();
    
    int ICtype = parser("-IC").asInt();
    
    
    switch (ICtype) {
        case 1:
        {
            _icPatient1(*grid);
            
            if (bAdaptivity)
            {
                Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _icPatient1);
                Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler,_icPatient1);
            }
        }
            break;
            
        case 7:
        {
            _icPatient7(*grid);
            
            if (bAdaptivity)
            {
                Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _icPatient7);
                Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler,_icPatient7);
            }
        }
            break;
            
        case 22:
        {
            _icPatient22(*grid);
            
            if (bAdaptivity)
            {
                Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _icPatient22);
                Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler,_icPatient22);
            }
        }
            break;
            
        default:
            break;
    }
    
    _dump(0);
    isDone              = false;
}

Glioma_HG_PatientData::~Glioma_HG_PatientData()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions

void Glioma_HG_PatientData:: _icPatient1(Grid<W,B>& grid)
{
//#ifdef BRUTUS
    MatrixD3D GM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/GM_1.dat");
    MatrixD3D WM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/WM_1.dat");
    MatrixD3D CSF("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/CSF_1.dat");
    MatrixD3D T1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/T1_1.dat");
    MatrixD3D T2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/T2_1.dat");
    MatrixD3D PET1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/FET1_1.dat");
    MatrixD3D PET2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/FET1_2.dat");
    MatrixD3D PET3("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/FET1_3.dat");
//#else
//    MatrixD3D GM("/Users/lipkova/WORK/Glioma/sourcestochasticF/Anatomy/Patient01/GM_1.dat");
//#endif
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    L = 16.2; //[cm]
    
    int GM_size = (int) (GM.getSizeX() * GM.getSizeY() * GM.getSizeZ());
    printf("GM_size=%i, sizeX=%i, sizeY=%i, sizeZ=%i\n",GM_size,(int)brainSizeX,(int)brainSizeY,(int)brainSizeZ);
    
    double brainHx = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. Z for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio

    
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
                    
                    int mappedBrainX = (int)floor( (x[0]) / brainHx  );
                    int mappedBrainY = (int)floor( (x[1]) / brainHy  );
                    int mappedBrainZ = (int)floor( (x[2]) / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= 8;
                    mappedBrainY -= 48;
                    
                    if ( ( mappedBrainX < 0 || mappedBrainX >= brainSizeX ) || ( mappedBrainY < 0 || mappedBrainY >= brainSizeY ) )
                    {
                    
                    }
                    else
                    {
                        
                        Real PGt     = GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PWt     = WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Pcsf    = CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT1     = T1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT2     = T2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet1   = PET1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet2   = PET2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet3   = PET3(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        
                        // Anatomy
                        if( (PGt + PWt) > 0.005)
                        {
                            block(ix,iy,iz).p_g     = PGt / (PGt + PWt);
                            block(ix,iy,iz).p_w     = PWt / (PGt + PWt);
                            block(ix,iy,iz).p_csf   = 0.;
                        }
                        else
                        {
                            block(ix,iy,iz).p_g     = 0.;
                            block(ix,iy,iz).p_w     = 0.;
                            block(ix,iy,iz).p_csf   = (Pcsf>0.005) ? 1. : 0.;
                        }
                        
                        // Binary Segmentations
                        block(ix,iy,iz).ux = PT1;
                        block(ix,iy,iz).uy = PT2;
                        
                        // PET - Integrate out time and restrict to T1 u T2 region
                        if ( (PT1+PT2) > 0.9)
                            block(ix,iy,iz).chi = Ppet1 + Ppet2 + Ppet3;
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


void Glioma_HG_PatientData:: _icPatient7(Grid<W,B>& grid)
{
    //#ifdef BRUTUS
    MatrixD3D GM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/GM_7.dat");
    MatrixD3D WM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/WM_7.dat");
    MatrixD3D CSF("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/CSF_7.dat");
    MatrixD3D T1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/T1_7.dat");
    MatrixD3D T2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/T2_7.dat");
    MatrixD3D PET1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/FET7_1.dat");
    MatrixD3D PET2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/FET7_2.dat");
    MatrixD3D PET3("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient07/FET7_3.dat");
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    L = 15.3; //[cm]
    
    int GM_size = (int) (GM.getSizeX() * GM.getSizeY() * GM.getSizeZ());
    printf("GM_size=%i, sizeX=%i, sizeY=%i, sizeZ=%i\n",GM_size,(int)brainSizeX,(int)brainSizeY,(int)brainSizeZ);
    
    double brainHx = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. Z for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio
    
    
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
                    
                    int mappedBrainX = (int)floor( (x[0]) / brainHx  );
                    int mappedBrainY = (int)floor( (x[1]) / brainHy  );
                    int mappedBrainZ = (int)floor( (x[2]) / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= 8;
                    mappedBrainY -= 48;
                    
                    if ( ( mappedBrainX < 0 || mappedBrainX >= brainSizeX ) || ( mappedBrainY < 0 || mappedBrainY >= brainSizeY ) )
                    {
                        
                    }
                    else
                    {
                        
                        Real PGt     = GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PWt     = WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Pcsf    = CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT1     = T1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT2     = T2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet1   = PET1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet2   = PET2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet3   = PET3(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        
                        // Anatomy
                        if( (PGt + PWt) > 0.005)
                        {
                            block(ix,iy,iz).p_g     = PGt / (PGt + PWt);
                            block(ix,iy,iz).p_w     = PWt / (PGt + PWt);
                            block(ix,iy,iz).p_csf   = 0.;
                        }
                        else
                        {
                            block(ix,iy,iz).p_g     = 0.;
                            block(ix,iy,iz).p_w     = 0.;
                            block(ix,iy,iz).p_csf   = (Pcsf>0.005) ? 1. : 0.;
                        }
                        
                        // Binary Segmentations
                        block(ix,iy,iz).ux = PT1;
                        block(ix,iy,iz).uy = PT2;
                        
                        // PET - Integrate out time and restrict to T1 u T2 region
                        if ( (PT1+PT2) > 0.9)
                            block(ix,iy,iz).chi = Ppet1 + Ppet2 + Ppet3;
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


void Glioma_HG_PatientData:: _icPatient22(Grid<W,B>& grid)
{
    //#ifdef BRUTUS
    MatrixD3D GM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/GM_22.dat");
    MatrixD3D WM("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/WM_22.dat");
    MatrixD3D CSF("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/CSF_22.dat");
    MatrixD3D T1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/T1_22.dat");
    MatrixD3D T2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/T2_22.dat");
    MatrixD3D PET1("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/FET22_1.dat");
    MatrixD3D PET2("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/FET22_2.dat");
    MatrixD3D PET3("/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient22/FET22_3.dat");
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    L = 15.0; //[cm]
    
    int GM_size = (int) (GM.getSizeX() * GM.getSizeY() * GM.getSizeZ());
    printf("GM_size=%i, sizeX=%i, sizeY=%i, sizeZ=%i\n",GM_size,(int)brainSizeX,(int)brainSizeY,(int)brainSizeZ);
    
    double brainHx = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. Z for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeZ-1)); // should be w.r.t. z for correct aspect ratio
    
    
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
                    
                    int mappedBrainX = (int)floor( (x[0]) / brainHx  );
                    int mappedBrainY = (int)floor( (x[1]) / brainHy  );
                    int mappedBrainZ = (int)floor( (x[2]) / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= 8;
                    mappedBrainY -= 48;
                    
                    if ( ( mappedBrainX < 0 || mappedBrainX >= brainSizeX ) || ( mappedBrainY < 0 || mappedBrainY >= brainSizeY ) )
                    {
                        
                    }
                    else
                    {
                        
                        Real PGt     = GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PWt     = WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Pcsf    = CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT1     = T1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real PT2     = T2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet1   = PET1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet2   = PET2(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Real Ppet3   = PET3(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        
                        // Anatomy
                        if( (PGt + PWt) > 0.005)
                        {
                            block(ix,iy,iz).p_g     = PGt / (PGt + PWt);
                            block(ix,iy,iz).p_w     = PWt / (PGt + PWt);
                            block(ix,iy,iz).p_csf   = 0.;
                        }
                        else
                        {
                            block(ix,iy,iz).p_g     = 0.;
                            block(ix,iy,iz).p_w     = 0.;
                            block(ix,iy,iz).p_csf   = (Pcsf>0.005) ? 1. : 0.;
                        }
                        
                        // Binary Segmentations
                        block(ix,iy,iz).ux = PT1;
                        block(ix,iy,iz).uy = PT2;
                        
                        // PET - Integrate out time and restrict to T1 u T2 region
                        if ( (PT1+PT2) > 0.9)
                            block(ix,iy,iz).chi = Ppet1 + Ppet2 + Ppet3;
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



void Glioma_HG_PatientData::_setScalingFactor()
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
                    maxPET = max(maxPET, block(ix,iy,iz).chi);
                    mass += block(ix,iy,iz).chi;
                    
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    cx += x[0] * block(ix,iy,iz).chi;
                    cy += x[1] * block(ix,iy,iz).chi;
                    cz += x[2] * block(ix,iy,iz).chi;
                }
        
    }
    
    cx = cx / mass;
    cy = cy / mass;
    cz = cz / mass;
    
    printf("maxPET= %f \n", maxPET);
    printf("Center of mass: cx=%f, cy=%f, cz=%f \n", cx,cy,cz);
    
    FILE * pFile;
    float buffer[3] = { cx , cy , cz };
    pFile = fopen ("HGG_TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);
    fclose (pFile);
}

// so PET signal is [0,1]
void Glioma_HG_PatientData::_rescalePETsignal()
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
                block(ix,iy,iz).chi = block(ix,iy,iz).chi * iscale;
    }
}

void Glioma_HG_PatientData::_detectSegmentationsBC(BoundaryInfo* boundaryInfo, const int nParallelGranularity)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_SegmentationBC_Operator<_DIM>  rhs;
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}

void Glioma_HG_PatientData::_getAvPETatSegmentationsBC()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    double T1bcPETav = 0.;
    double T2bcPETav = 0.;
    
    int nT1 = 0;
    int nT2 = 0;
        
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
    
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        if ( block(ix,iy,iz).t1bc == 1. )
                        {
                            T1bcPETav += block(ix,iy,iz).chi;
                            nT1++;
                        }
                        
                        if ( block(ix,iy,iz).t2bc == 1. )
                        {
                            T2bcPETav += block(ix,iy,iz).chi;
                            nT2++;
                        }
                    }
        
    }
    
    
    T1bcPETav = T1bcPETav / nT1;
    T2bcPETav = T2bcPETav / nT2;
    
    printf("T1bcPETav=%f, T2bcPETav=%f,\n", T1bcPETav,T2bcPETav );
}


#pragma mark DumpingOutput
void Glioma_HG_PatientData:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%dD_Patient42_Data%04d",_DIM, counter);
        
        if( _DIM == 2)
        {
            IO_VTKNative<W,B, 2,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
        else
        {
            IO_VTKNative3D<W,B, 11,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_HG_PatientData::_dumpOutput()
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    MatrixD3D PET(gpd,gpd,gpd);
    MatrixD3D T1(gpd,gpd,gpd);
    MatrixD3D T2(gpd,gpd,gpd);
    
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
                        T1(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).ux;
                        T2(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).uy;
                        PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).chi;
                    }
        }
        else
        {
            for(int iz=0; iz<blockSize; iz++)
                for(int iy=0; iy<blockSize; iy++)
                    for(int ix=0; ix<blockSize; ix++)
                    {
                        T1(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                        T2(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                        PET(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                    }
        }
    }
    
    char filename1[256];
    sprintf(filename1,"PET_data.dat");
    PET.dump(filename1);
    
    sprintf(filename1,"T1_data.dat");
    T1.dump(filename1);
    
    sprintf(filename1,"T2_data.dat");
    T2.dump(filename1);
}


void Glioma_HG_PatientData::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    _setScalingFactor();
    _rescalePETsignal();
    _detectSegmentationsBC(boundaryInfo, nParallelGranularity);
    _getAvPETatSegmentationsBC();
    _dump(1);
        
    _dumpOutput();
    
    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}
