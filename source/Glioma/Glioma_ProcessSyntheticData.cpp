//
//  Glioma_ProcessSyntheticData.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 12/11/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_ProcessSyntheticData.h"


static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_ProcessSyntheticData::Glioma_ProcessSyntheticData(int argc, const char ** argv): parser(argc, argv)
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
    
    bAdaptivity = parser("-adaptive").asBool();
    bVerbose = parser("-verbose").asBool();
    
    _ic(*grid);
    
    if (bAdaptivity)
    {
        Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _ic);
        Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler,_ic);
    }
    
    _dump(0);
    isDone              = false;
}

Glioma_ProcessSyntheticData::~Glioma_ProcessSyntheticData()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
// Patient Brain anatomy
void Glioma_ProcessSyntheticData:: _ic(Grid<W,B>& grid)
{
    MatrixD3D tumor("HGG_data.dat");

    int tumor_size = (int) (tumor.getSizeX() * tumor.getSizeY() * tumor.getSizeZ());
    printf("tumor_size=%i\n",tumor_size);
    
    
//    // Synthetic 1
//    Real ucT2 = 0.25;
//    Real ucT1 = 0.7;
    
    // Synthetic 2
    Real ucT2 = 0.3;
    Real ucT1 = 0.8;
    
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
                    block(ix,iy,iz).phi  = ( petSignal >= ucT2 ) ? petSignal : 0.;
                    block(ix,iy,iz).t1bc = ( petSignal >= ucT1 ) ? 1. : 0. ;
                    block(ix,iy,iz).t2bc = ( petSignal >= ucT2 ) ? 1. : 0. ;
                    
                }
        
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



// compute normalizing factor & center of mass
void Glioma_ProcessSyntheticData::_getPETstatistic()
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
    
    FILE * pFile;
    float buffer[3] = { cx , cy , cz };
    pFile = fopen ("HGG_TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);
    fclose (pFile);
    
}

// scale PET signal is [0,1]
void Glioma_ProcessSyntheticData::_normalizePETsignal()
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
void Glioma_ProcessSyntheticData::_addNoise()
{
    double xi_1, xi_2;
//    double alpha = 0.05;  // 5% noise
    double alpha = 0.02;  // 5% noise

    int n = 0;
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
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
void Glioma_ProcessSyntheticData:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%dD_Patient42_Data%04d",_DIM, counter);
        
        IO_VTKNative3D<W,B, 7,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_ProcessSyntheticData::_dumpOutput()
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
                        PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).phi;
                        T1(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).t1bc;
                        T2(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).t2bc;
                    }

            
        }
        else
        {
            for(int iz=0; iz<blockSize; iz++)
                for(int iy=0; iy<blockSize; iy++)
                    for(int ix=0; ix<blockSize; ix++)
                    {
                        PET(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                        T1(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                        T2(ix + blockSize * info.index[0], iy + blockSize * info.index[1], iz + blockSize * info.index[2] ) = 0.;
                    }
        }
    }
    
    char filename[256];
    sprintf(filename,"PET_data.dat");
    PET.dump(filename);
    
    sprintf(filename,"T1_data.dat");
    T1.dump(filename);
    
    sprintf(filename,"T2_data.dat");
    T2.dump(filename);
    
}



void Glioma_ProcessSyntheticData::_readInTumorPosition(vector<Real>& tumorIC )
{
    typedef float dataType;
    
    FILE* fp;
    fp = fopen("HGG_TumorIC.bin", "rb");
    if (fp == NULL) {fputs ("File error", stderr); exit (1);}
    
    // obtain file size
    fseek (fp , 0 , SEEK_END);
    long int size = ftell (fp);
    rewind (fp);
    
    
    // allocate memory to contain the whole file:
    dataType * buffer;
    buffer = (dataType*) malloc (sizeof(dataType)*size);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
    
    
    // copy the file into the buffer:
    size_t result;
    result = fread (buffer, 1, size, fp);
    if (result != size) {fputs ("Reading error",stderr); exit (3);}
    
    
    for (int i = 0; i < 3; ++i)
        tumorIC[i] = (Real)buffer[i];
    
    free(buffer);
    fclose (fp);
}

void Glioma_ProcessSyntheticData::_dumpSubBrainPoints(Grid<W,B>& grid )
{
    if (bAdaptivity)
    {
        printf("Aborting ... dumpBrainPoints needs uniform grid");
        abort();
    }
    
    int points   = 0;
    vector<Real> cm(3);
    _readInTumorPosition(cm);
    printf("cm: cm(0) = %f, cm(1) =%f, cm(2)=%f \n", cm[0], cm[1], cm[2]);
    
    float radius = parser("-radius").asDouble(0.21);  //made it max_tumor radius + 2cm margin
    
    vector<double>     brain;
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++ )
            for(int iy=0; iy<B::sizeY; iy++ )
                for(int ix=0; ix<B::sizeX; ix++ )
                {
                    Real x[3];
                    info.pos(x, ix, iy, iz);
                    
                    const Real p[3] = {x[0] - cm[0], x[1] - cm[1], x[2] - cm[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                    
                    Real tissue = block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if(( tissue > 0 )&&(dist <= radius))
                    {
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        brain.push_back(gix);
                        brain.push_back(giy);
                        brain.push_back(giz);
                        
                        points++;
                    }
                }
    }
    
    
    printf("points=%i \n", points);
    std::cout<<"points="<<points<<std::endl;
    
    MatrixD2D out(points,3);
    for (int i = 0; i<points; i++)
    {
        out(i,0) = brain[i*3    ];
        out(i,1) = brain[i*3 + 1];
        out(i,2) = brain[i*3 + 2];
    }
    
    char filename2[256];
    sprintf(filename2,"SphereR_%f.dat",radius);
    out.dump(filename2);
    
}



void Glioma_ProcessSyntheticData::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    

        _getPETstatistic();
        _normalizePETsignal();
        _dump(1);
        
        _addNoise();
        _dump(2);
        
        _getPETstatistic();
        
        _normalizePETsignal();
        _dump(3);
        
        _dumpOutput();
    
    
//    _dumpSubBrainPoints(*grid);  //<-- don't use this one !!! or fixed it first

    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}

