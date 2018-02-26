//
//  Glioma_HG_Data.cpp
//  GliomaXcode
//
//  Created by Lipkova on 18/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_HG_Data.h"

static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_HG_Data::Glioma_HG_Data(int argc, const char ** argv): parser(argc, argv)
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
    
    _icTumor(*grid);
    
    if (bAdaptivity)
    {
        Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _icTumor);
        Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler,_icTumor);
    }
    
    _dump(0);
    isDone              = false;
    ucT1                = 0.7;
    ucT2                = 0.2;
}

Glioma_HG_Data::~Glioma_HG_Data()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
// Patient Brain anatomy - subject04 from BrainWeb database
void Glioma_HG_Data:: _icTumor(Grid<W,B>& grid)
{
#ifdef BRUTUS
    MatrixD3D tumor("HGG_data.dat");
#else
    MatrixD3D tumor("/Users/lipkova/WORK/Glioma/sourcestochasticF/makefile/HGG_data.dat");
#endif
    int tumor_size = (int) (tumor.getSizeX() * tumor.getSizeY() * tumor.getSizeZ());
    printf("tumor_size=%i\n",tumor_size);
    
    double xi_1, xi_2;
    double eps = 0.02;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        int n = 0;
        
        if(_DIM==2)
        {
            
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy);
                    int middle = (int) (0.5*tumor.getSizeZ()) ;
                    block(ix,iy).phi = tumor(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1],middle );
                    
                    if(  (n % 2) == 0 )
                    {
                        double u1 = drand48();
                        double u2 = drand48();
                        
                        xi_1 = sqrt( - 2.0*log( u1 ) ) * cos( 2.0*M_PI*u2 );
                        xi_2 = sqrt( - 2.0*log( u1 ) ) * sin( 2.0*M_PI*u2 );
                    }
                    
                    if( block(ix,iy).phi > 0. )
                    {
                        block(ix,iy).chi = max( 0., block(ix,iy).phi + eps * xi_1);
                        xi_1 = xi_2;
                        
                        n++;
                    }
                    
                }
        }
        else
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        double x[3];
                        info.pos(x, ix, iy, iz);
                        
                        block(ix,iy,iz).phi = tumor(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2]);
                        
                        if(  (n % 2) == 0 )
                        {
                            double u1 = drand48();
                            double u2 = drand48();
                            
                            xi_1 = sqrt( - 2.0*log( u1 ) ) * cos( 2.0*M_PI*u2 );
                            xi_2 = sqrt( - 2.0*log( u1 ) ) * sin( 2.0*M_PI*u2 );
                        }
                        
                        
                        if( block(ix,iy,iz).phi > 0.  )
                        {
                            block(ix,iy,iz).chi = max( 0., block(ix,iy,iz).phi + eps * xi_1);
                            xi_1 = xi_2;
                            
                            n++;
                        }
                        
                    }
        }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}

void Glioma_HG_Data::_setScalingFactor()
{
    maxPET = 0.;
    maxPhi = 0.;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        if(_DIM==2)
        {
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    maxPET = max(maxPET, block(ix,iy).chi);
                    maxPhi = max(maxPhi, block(ix,iy).phi);
                    
                }
        }
        else
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        maxPET = max(maxPET, block(ix,iy,iz).chi);
                        maxPhi = max(maxPhi, block(ix,iy,iz).phi);
                        
                    }
        }
    }
    
    printf("maxPET= %f, maxPhi=%f \n", maxPET,maxPhi);
    
}

// so PET signal is [0,1]
void Glioma_HG_Data::_rescalePETsignal()
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


void Glioma_HG_Data::_generateBinaryData()
{
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        if(_DIM==2)
        {
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    block(ix,iy).ux = (block(ix,iy).phi < ucT1) ? 0. : 1.;
                    block(ix,iy).uy = (block(ix,iy).phi < ucT2) ? 0. : 1.;
                }
            
        }
        else
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        block(ix,iy,iz).ux = (block(ix,iy,iz).phi < ucT1) ? 0. : 1.;
                        block(ix,iy,iz).uy = (block(ix,iy,iz).phi < ucT2) ? 0. : 1.;
                    }
        }
    }
}

void Glioma_HG_Data::_generatePETdata()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        if(_DIM==2)
        {
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy).uz = ( (block(ix,iy).ux + block(ix,iy).uy)>0) ? block(ix,iy).chi : 0.;
        }
        else
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                        block(ix,iy,iz).uz = ( (block(ix,iy,iz).ux + block(ix,iy,iz).uy)>0) ? block(ix,iy,iz).chi : 0.;
            
        }
    }
}

void Glioma_HG_Data::_getAvPETatSegmentationsBC()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    double T1bcPETav = 0.;
    double T2bcPETav = 0.;
    
    double T1bcav = 0.;
    double T2bcav = 0.;
    
    
    int nT1 = 0;
    int nT2 = 0;
    
    double h = 1./(blocksPerDimension * blockSize);
    double eps = 2. * h;
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        if(_DIM==2)
        {
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    if ( block(ix,iy).t1bc == 1. )
                    {
                        T1bcPETav   += block(ix,iy).chi;
                        T1bcav      += block(ix,iy).phi;
                        
                        block(ix,iy).p_w = 100.;
                        
                        nT1++;
                    }
                    
                    if ( block(ix,iy).t2bc == 1. )
                    {
                        T2bcPETav   += block(ix,iy).chi;
                        T2bcav      += block(ix,iy).phi;
                        
                        nT2++;
                    }
                }
            
        }
        else
        {
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        if ( block(ix,iy,iz).t1bc == 1. )
                        {
                            T1bcPETav += block(ix,iy,iz).chi;
                            T1bcav    += block(ix,iy,iz).phi;
                            block(ix,iy,iz).p_w = 100.;
                            
                            nT1++;
                        }
                        
                        if ( block(ix,iy,iz).t2bc == 1. )
                        {
                            T2bcPETav += block(ix,iy,iz).chi;
                            T2bcav    += block(ix,iy,iz).phi;
                            
                            nT2++;
                        }
                    }
        }
    }
    
    
    T1bcPETav = T1bcPETav / nT1;
    T2bcPETav = T2bcPETav / nT2;
    
    T2bcav = T2bcav/nT2;
    T1bcav = T1bcav/nT1;
    
    
    printf("T1bcPETav=%f, T2bcPETav=%f, ucT1av = %f, ucT2av=%f \n", T1bcPETav,T2bcPETav,T1bcPETav*maxPhi, T2bcPETav*maxPhi );
    printf("T1bcav=%f, T2bcav=%f, nT1=%i, nT2 =%i  \n",T1bcav, T2bcav, nT1, nT2 );
    
}

void Glioma_HG_Data::_detectSegmentationsBC(BoundaryInfo* boundaryInfo, const int nParallelGranularity)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_SegmentationBC_Operator<_DIM>  rhs;
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}

void Glioma_HG_Data::_computeError(Grid<W,B>& grid)
{
    Real L1 = 0.0;
    Real L2 = 0.0;
    Real LI = 0.0;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for reduction (+:L1,L2,LI)
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real dv = vInfo[i].h[0]*vInfo[i].h[1]*vInfo[i].h[2];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    const Real err = fabs(block(ix,iy,iz).uz * maxPhi - block(ix,iy,iz).phi);
                    L1 += err * dv;
                    L2 += err * err * dv;
                    LI = std::max(LI,err);
                }
    }
    
    L2 =std::sqrt(L2);
    
    const int res = blockSize * blocksPerDimension;
    
    printf("========= PET DATA ERRORS %d ========\n",res);
    printf("L1, L2, LI: %e %e %e\n",L1,L2,LI);
    printf("========= END PET DATA ===========\n");
}

#pragma mark DumpingOutput
void Glioma_HG_Data:: _dump(int counter)
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
            IO_VTKNative3D<W,B, 8,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_HG_Data::_dumpOutput()
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
                        PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).uz;
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


// pick N points at random from tumour and N from background
void Glioma_HG_Data::_pickRandomPoints(int N)
{
    float cm[3] = {0.315, 0.65, 0.5};
    float radius = 0.2;
    
    vector<int> tumor;
    vector<int> background;
    
    
    
    // pick area around tumour
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
                    
                    const Real p[3] = {x[0] - cm[0], x[1] - cm[1], x[2] - cm[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    
                    if(dist <= radius)
                    {
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        if( block(ix,iy,iz).chi > 0 )
                        {
                            tumor.push_back(gix);
                            tumor.push_back(giy);
                            tumor.push_back(giz);
                        }
                        else
                        {
                            background.push_back(gix);
                            background.push_back(giy);
                            background.push_back(giz);
                        }
                    }
                }
    }
    
    
    // pick 50 independt points form tumour and 50 from background
    
    MatrixD2D out(2*N,3);
    int Nt = tumor.size() / 3.;
    int Nb = background.size() / 3.;
    
    for (int i = 0; i< N; i++)
    {
        int ut = round( 0. + Nt * drand48());
        int ub = round( 0. + Nb * drand48());
        
        out(i,0) =  tumor[ut*3    ];
        out(i,1) =  tumor[ut*3 + 1];
        out(i,2) =  tumor[ut*3 + 2];
        
        out(i + N,0) =  background[ub*3    ];
        out(i + N,1) =  background[ub*3 + 1];
        out(i + N,2) =  background[ub*3 + 2];
        
    }
    
    char filename1[256];
    sprintf(filename1,"Points.dat");
    out.dump(filename1);
}


void Glioma_HG_Data::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    _setScalingFactor();
    _rescalePETsignal();
    _setScalingFactor();
    _generateBinaryData();
    _detectSegmentationsBC(boundaryInfo, nParallelGranularity);
    _getAvPETatSegmentationsBC();
    _generatePETdata();
    _computeError(*grid);
    _pickRandomPoints(100);
    _dump(1);
    
    
    _dumpOutput();
    
    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}
