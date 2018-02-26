//
//  Glioma_Bone_BMD_UQ.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 09/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_Bone_BMD_UQ.h"

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_Bone_BMD_UQ::Glioma_Bone_BMD_UQ(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////                  Bone UQ                     ////////////////\n");
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
    
    int ICtype = parser("-IC").asInt(0);
    
    switch (ICtype) {
            case 0:
            _ic_Vertebra(*grid);
            break;
            case 1:
            _ic_Patient_Vertebrea(*grid);
            break;
            
        default:
            break;
    }

        _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
    
}

Glioma_Bone_BMD_UQ:: ~Glioma_Bone_BMD_UQ()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
void Glioma_Bone_BMD_UQ::_ic_Vertebra(Grid<W,B>& grid)
{
#if !defined(Template16) && !defined(Template23)
#define Template16
#endif
    
#ifdef Template16
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/vertebra_16.dat \n");
    MatrixD3D BMD(  "/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/vertebra_16.dat");
#endif
    
#ifdef Template23

#ifdef Healthy
    printf("healthy:\n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/vertebra_23.dat");
#endif
    
#ifdef Age_40
    printf("Age 40:\n");
    MatrixD3D BMD(  "/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/40_vertebra_23.dat");
#endif
    
#ifdef Age_50
    printf("Age 50:\n");
    MatrixD3D BMD(  "/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/50_vertebra_23.dat");
#endif
    
#ifdef Age_60
    printf("Age 60:\n");
    MatrixD3D BMD(  "/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/60_vertebra_23.dat");
#endif
    
#ifdef Age_70
    printf("Age 70:\n");
    MatrixD3D BMD(  "/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Template/70_vertebra_23.dat");
#endif
#endif
    
    
    int boneSizeX = (int) BMD.getSizeX();
    int boneSizeY = (int) BMD.getSizeY();
    int boneSizeZ = (int) BMD.getSizeZ();
    
    int boneSizeMax = max(boneSizeX, max(boneSizeY,boneSizeZ));
    L    = boneSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 14.6 cm
    
    printf("boneSizeX=%i, boneSizeY=%i, boneSizeZ= %i \n", boneSizeX, boneSizeY, boneSizeZ);
    
    double boneHx = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHy = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHz = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(3);
    _readInTumorPosition(tumor_ic);
    
    const Real tumorRadius = 0.005;
    const Real smooth_sup  = 2.;		// suppor of smoothening, over how many gp to smooth
    Real minBMD = 50;
    
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real h = vInfo[0].h[0];
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBoneX = (int)floor( x[0] / boneHx  );
                    int mappedBoneY = (int)floor( x[1] / boneHy  );
                    int mappedBoneZ = (int)floor( x[2] / boneHz  );
                    
                    // aspect ratio correction
                    mappedBoneX -= (int) ( (boneSizeMax - boneSizeX) * 0.5);
                    mappedBoneY -= (int) ( (boneSizeMax - boneSizeY) * 0.5);
                    mappedBoneZ -= (int) ( (boneSizeMax - boneSizeZ) * 0.5);
                    
                    
                    if ( (mappedBoneX >= 0 && mappedBoneX < boneSizeX) && (mappedBoneY >= 0 && mappedBoneY < boneSizeY) && (mappedBoneZ >= 0 && mappedBoneZ < boneSizeZ) )
                    {
                        // anatomy
                        Real bmd = BMD(mappedBoneX,mappedBoneY,mappedBoneZ);
                        
                        block(ix,iy,iz).bmd  = (bmd > minBMD) ? bmd : 0.;   // remove background signal, BMD [HU]
                        
                        // tumor
                        const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                        const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                        const Real psi = (dist - tumorRadius)*iw;
                        
                        if ((psi < -1)&& ( bmd > minBMD ))		// we are in tumor
                            block(ix,iy,iz).phi = 1.0;
                        else if(( (-1 <= psi) && (psi <= 1) )&& ( bmd > minBMD ))
                            block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                        else
                            block(ix,iy,iz).phi = 0.0;
                        
                        
                    }
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}

void Glioma_Bone_BMD_UQ::_ic_Patient_Vertebrea(Grid<W,B>& grid)
{
#if !defined(Patient10) && !defined(Patient6) && !defined(Patient1)
#define Patient6
#define V5
#endif
    
#if defined(Patient1) && defined(V1)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V1/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V1/P1_CT_V1.dat");
#endif
#if  defined(Patient1) && defined(V2)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V2/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V2/P1_CT_V2.dat");
#endif
#if  defined(Patient1) && defined(V3)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V3/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient1/V3/P1_CT_V3.dat");
#endif
 
#if defined(Patient6) && defined(V2)
     printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V2/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V2/P6_CT_V2.dat");
#endif
    
#if  defined(Patient6) && defined(V3)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V3/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V3/P6_CT_V3.dat");
#endif
    
#if defined(Patient6) && defined(V4)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V4/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V4/P6_CT_V4.dat");
#endif
    
#if defined(Patient6) && defined(V5)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V5/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient6/V5/P6_CT_V5.dat");
#endif
    
#if defined(Patient10) && defined(V1)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V1/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V1/P10_CT_V1.dat");
#endif
    
#if  defined(Patient10) && defined(V2)
    printf("Reading anatomy from: /cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V2/ \n");
    MatrixD3D BMD("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Bone/Patient10/V2/P10_CT_V2.dat");
#endif
    
    
    int boneSizeX = (int) BMD.getSizeX();
    int boneSizeY = (int) BMD.getSizeY();
    int boneSizeZ = (int) BMD.getSizeZ();
    
    int boneSizeMax = max(boneSizeX, max(boneSizeY,boneSizeZ));
    L    = boneSizeMax * 0.1 * 1.52344 * 0.5;   // voxel spacing 1.52344 mm (*0.5 due to interpolation), convert from mm to cm  //
    
    printf("boneSizeX=%i, boneSizeY=%i, boneSizeZ= %i \n", boneSizeX, boneSizeY, boneSizeZ);
    
    double boneHx = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHy = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double boneHz = 1.0 / ((double)(boneSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(3);
    _readInTumorPosition(tumor_ic);
    
    const Real tumorRadius = 0.005;
    const Real smooth_sup  = 2.;		// suppor of smoothening, over how many gp to smooth
    Real minBMD = 100;
    
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real h = vInfo[0].h[0];
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
        for(int iy=0; iy<B::sizeY; iy++)
        for(int ix=0; ix<B::sizeX; ix++)
        {
            double x[3];
            info.pos(x, ix, iy, iz);
            
            /* Anatomy */
            int mappedBoneX = (int)floor( x[0] / boneHx  );
            int mappedBoneY = (int)floor( x[1] / boneHy  );
            int mappedBoneZ = (int)floor( x[2] / boneHz  );
            
            // aspect ratio correction
            mappedBoneX -= (int) ( (boneSizeMax - boneSizeX) * 0.5);
            mappedBoneY -= (int) ( (boneSizeMax - boneSizeY) * 0.5);
            mappedBoneZ -= (int) ( (boneSizeMax - boneSizeZ) * 0.5);
            
            
            if ( (mappedBoneX >= 0 && mappedBoneX < boneSizeX) && (mappedBoneY >= 0 && mappedBoneY < boneSizeY) && (mappedBoneZ >= 0 && mappedBoneZ < boneSizeZ) )
            {
                // anatomy
                Real bmd = BMD(mappedBoneX,mappedBoneY,mappedBoneZ);
                block(ix,iy,iz).tmp  = bmd;
                block(ix,iy,iz).chi  = (bmd >= minBMD) ? 1. : 0.;     // domain charact. function
                
                // tumour diffusion D = d/psi^4, psi = (bmd*0.01)
                Real tmp   = bmd * 0.01;
                Real tmpA = tmp;

#ifdef bmd2
                tmpA  = tmp * tmp;
#endif
#ifdef bmd3
                tmpA  = tmp * tmp * tmp;
#endif
#ifdef bmd4
                tmpA  = tmp * tmp * tmp * tmp;
#endif
                Real ItmpA = (tmp > 0.) ? (1./tmpA) : 0. ;
                
                block(ix,iy,iz).bmd  = ItmpA * block(ix,iy,iz).chi;  // diffusion term restricted to bone region
                
                // tumor
                const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                const Real psi = (dist - tumorRadius)*iw;
                
                if ((psi < -1)&& ( bmd > minBMD ))		// we are in tumor
                block(ix,iy,iz).phi = 1.0;
                else if(( (-1 <= psi) && (psi <= 1) )&& ( bmd > minBMD ))
                block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                else
                block(ix,iy,iz).phi = 0.0;
                
            }
            
        }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}




void Glioma_Bone_BMD_UQ::_readInTumorPosition(vector<Real>& tumorIC )
{
    typedef float dataType;
    
    FILE* fp;
    fp = fopen("TumorIC.bin", "rb");
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

#pragma mark ReactionDiffusion
void Glioma_Bone_BMD_UQ::_reactionDiffusionStep_BMD(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dscale, const Real rho, double dt)
{

    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_BMD_ReactionDiffusion_Operator<_DIM>  rhs(Dscale, rho);
    UpdateTumor                          <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}


#pragma mark DumpingOutput
void Glioma_Bone_BMD_UQ:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        
#ifdef Template16
        sprintf(filename,"BMD_Spine_v16_Data%04d", counter);
#endif
        
#ifdef Template23
        
#ifdef Healthy
        sprintf(filename,"BMD_Spine_v23_Data%04d", counter);
#endif
#ifdef Age_40
        sprintf(filename,"BMD_Spine_40_v23_Data%04d", counter);
#endif
#ifdef Age_50
        sprintf(filename,"BMD_Spine_50_v23_Data%04d", counter);
#endif
#ifdef Age_60
        sprintf(filename,"BMD_Spine_60_v23_Data%04d", counter);
#endif 
#ifdef Age_70
        sprintf(filename,"BMD_Spine_70_v23_Data%04d", counter);
#endif
#endif
        
#if defined(Patient6) && defined(V2)
        sprintf(filename,"Patient6_V2_Data%04d", counter);
#endif
#if defined(Patient6) && defined(V3)
        sprintf(filename,"Patient6_V3_Data%04d", counter);
#endif
#if defined(Patient6) && defined(V4)
        sprintf(filename,"Patient6_V4_Data%04d", counter);
#endif
#if defined(Patient6) && defined(V5)
        sprintf(filename,"Patient6_V5_Data%04d", counter);
#endif
#if defined(Patient10) && defined(V1)
        sprintf(filename,"Patient10_V1_Data%04d", counter);
#endif
#if defined(Patient10) && defined(V2)
        sprintf(filename,"Patient10_V2_Data%04d", counter);
#endif
        IO_VTKNative3D<W,B, 4,0 > vtkdumper;
        vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
        
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_Bone_BMD_UQ::_dumpUQoutput(Grid<W,B>& grid)
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    if(bVerbose) printf("bpd=%i, bs=%i, hf=%f,\n",blocksPerDimension,blockSize,hf);
    
    MatrixD3D tumor(gpd,gpd,gpd);
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        double h = info.h[0];
        
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    //mapped coordinates
                    int mx = (int)floor( (x[0]) / hf  );
                    int my = (int)floor( (x[1]) / hf  );
                    int mz = (int)floor( (x[2]) / hf  );
                    
                    
                    if(h==hf)
                        tumor(mx,my,mz) = block(ix,iy,iz).phi;
                    else if(h == 2.*hf)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else if (h == 3.*hf)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                }
        
    }
    
    char filename2[256];
    sprintf(filename2,"UQ_data.dat");
    tumor.dump(filename2);
    
}



void Glioma_Bone_BMD_UQ::run()
{
    
    bool bProfiler = 0;
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* read in patien specific parameters*/
    ifstream mydata("Bone_InputParameters.txt");
    Real Dscale, maxD;
    double tend;
    
    if (mydata.is_open())
    {
        mydata >> Dscale;
        mydata >> rho;
        mydata >> tend;
        mydata.close();
    }
    
    /*rescale*/
    Real minBMD = 100;
    Real minPsi = (minBMD * 0.01);
    Real minPsiA = minPsi;
    
#ifdef bmd2
    minPsiA = minPsi * minPsi;
#endif
#ifdef bmd3
    minPsiA = minPsi * minPsi * minPsi;
#endif
#ifdef bmd4
    minPsiA = minPsi * minPsi * minPsi * minPsi;
#endif
    
    Dscale= Dscale/(L*L);
    maxD = Dscale / minPsiA;   // (  = c/min(bmd) min BMD 100 HU );
    
    double t			= 0.0;
    int iCounter        = 1;
    
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.99 * h*h / ( 2.* _DIM * maxD);
    if(bVerbose)  printf("Dscale=%e, maxD=%e, dt= %f, rho=%f , h=%f\n", Dscale, maxD, dt,rho,h);
    
    
    
    while (t <= tend)
    {
        if(bProfiler) profiler.getAgent("RD_Step").start();
        _reactionDiffusionStep_BMD(boundaryInfo, nParallelGranularity, Dscale, rho, dt);
        if(bProfiler) profiler.getAgent("RD_Step").stop();
        
        
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( t >= ((double)(whenToWrite)) )
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
    
    
    if(bProfiler) profiler.getAgent("UQ_output").start();
    _dumpUQoutput(*grid);
    if(bProfiler) profiler.getAgent("UQ_output").stop();
    
    _dump(100);
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}