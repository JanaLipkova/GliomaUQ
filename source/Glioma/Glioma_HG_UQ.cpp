//
//  Glioma_HG_UQ.cpp
//  GliomaXcode
//
//  Created by Lipkova on 10/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_HG_UQ.h"

// need biger stencil for the refinment !!!
static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_HG_UQ::Glioma_HG_UQ(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////             High Grade Glioma UQ             ////////////////\n");
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
    int ICtype = 0;
    ICtype = parser("-IC").asInt(1);
    
    pID =  parser("-pID").asInt();
    _ic_PatientCase(*grid, pID);
                    
    
    _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
    
}

Glioma_HG_UQ::~Glioma_HG_UQ()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions

// Patient 01 data
// 1) read in anatomies - rescaled to [0,1]^3
// 2) read in tumor center of mass + initialize tumor around
// 3) set length of brain
void Glioma_HG_UQ::_ic_PatientCase(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef BRUTUS
    sprintf(dataFolder,"/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/");
#elif defined(KRAKEN)
    sprintf(dataFolder,"/home/jana/Work/GliomaAdvance/source/Anatmoy/");
#elif defined(PLURIPOTENT)
    sprintf(dataFolder,"/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/");
#elif defined(LRZ_CLUSTER)
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/UQ_Section/source/Anatmoy/HGG/");
#else
    sprintf(dataFolder,"../../Anatmoy/");
#endif
    
    sprintf(patientFolder, "%sPatient%02d/P%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s", patientFolder);
    
    sprintf(anatomy, "%s_GM.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_WM.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_CSF.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%s_T1.dat", patientFolder);
    MatrixD3D T1(anatomy);
    sprintf(anatomy, "%s_T2.dat", patientFolder);
    MatrixD3D T2(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    
    printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(_DIM);
    _readInTumorPosition(tumor_ic);
    

    const Real tumorRadius = 0.005;
    const Real smooth_sup  = 2.;		// suppor of smoothening, over how many gp to smooth
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        //        const Real h = vInfo[0].h[0];
        
        const Real h = 1./128;
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    Real PGt, PWt, Pcsf, PT1, PT2;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt = 0.;
                        PWt = 0.;
                        Pcsf = 0.;
                        PT1 = 0.;
                        PT2 = 0.;
                    }
                    else
                    {
                        PGt     =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt     =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf    = CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PT1     =  T1(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PT2     =  T2(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    
                    // remove fluid below tumour that would be normally pushed away by growing tumour, not for P12, here is CSF was manaully reinforce to maintain hemisphere sepratation (csf was already removed from tumour regions, however at the interesection this step could destroy the hemisphere separations, therefore we omit it in this case
                    Real MRIsignal = PT1 + PT2;
                    Pcsf  = ( MRIsignal > 0.) ? 0. : Pcsf;
                    
#if defined(Patient20)
                    if ( (MRIsignal > 0.) && (PGt + PWt < 0.1) )
                        PGt = 1.;
#endif
                    
                    // Anatomy
                    double all = PWt + PGt + Pcsf;
                    
                    
                    
#ifdef BRAIN_MASK
                    
                    block(ix,iy,iz).p_w = 1.;
                    
#else
                    if(all > 0)
                    {
                        // normalize
                        PGt    = PGt  / all;
                        PWt    = PWt  / all;
                        Pcsf   = Pcsf / all;
                        
                        Pcsf = ( Pcsf > 0.1 ) ? 1. : Pcsf;  // threasholding to ensure hemisphere separations
                        block(ix,iy,iz).p_csf = Pcsf;
                        
                        if(Pcsf  < 1.)
                        {
                            block(ix,iy,iz).p_csf  = Pcsf / (Pcsf + PWt + PGt);
                            block(ix,iy,iz).p_w    = PWt  / (Pcsf + PWt + PGt);
                            block(ix,iy,iz).p_g    = PGt  / (Pcsf + PWt + PGt);
                        }
                    }
                    
#endif
                    
                    // tumor
                    const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    const Real psi = (dist - tumorRadius)*iw;
                    
                    if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                        block(ix,iy,iz).phi = 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                        block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi = 0.0;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



void Glioma_HG_UQ::_readInTumorPosition(vector<Real>& tumorIC )
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
    
    
    for (int i = 0; i < _DIM; ++i)
    {
        tumorIC[i] = (Real)buffer[i];
        // std::cout<<"IC["<<i<<"]="<<tumorIC[i]<<std::endl;
    }
    
    free(buffer);
    fclose (fp);
}

#pragma mark ReactionDiffusion
void Glioma_HG_UQ::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    UpdateTumor                     <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}


#pragma mark DumpingOutput
void Glioma_HG_UQ:: _dump(int counter)
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
            IO_VTKNative3D<W,B, 9,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_HG_UQ::_dumpUQoutput(Grid<W,B>& grid)
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
    sprintf(filename2,"HGG_data.dat");
    tumor.dump(filename2);
    
}

void Glioma_HG_UQ::_dumpBrainPoints(Grid<W,B>& grid )
{
    printf("InBrain\n");
    
    if (bAdaptivity)
    {
        printf("Aborting ... dumpBrainPoints needs uniform grid");
        abort();
    }
    
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    int s=0;
    
    vector<double>     brain;
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real tissue = block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    if( tissue > 0 )
                    {
                        
                        double x[3];
                        info.pos(x, ix, iy, iz);
                        
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        brain.push_back(gix);
                        brain.push_back(giy);
                        brain.push_back(giz);
                        
                        s++;
                    }
                }
    }
    
    
    MatrixD2D out(s,3);
    for (int i = 0; i<s; i++)
    {
        out(i,0) = brain[i*3    ];
        out(i,1) = brain[i*3 + 1];
        out(i,2) = brain[i*3 + 2];
    }
    
    
    char filename2[256];
    sprintf(filename2,"Shit.dat");
    out.dump(filename2);
    
}



void Glioma_HG_UQ::_dumpSubBrainPoints(Grid<W,B>& grid )
{
    if (bAdaptivity)
    {
        printf("Aborting ... dumpBrainPoints needs uniform grid");
        abort();
    }
    
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    int points   = 0;
    float cm[3]  = {0.6, 0.45, 0.65};
    float radius = 0.21;  // 0.15 max tumour radius, 2.17 cm safty margin-> 0.1 L*
    
    
    vector<double>     brain;
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++ )
                for(int ix=0; ix<B::sizeX; ix++ )
                {
                    double x[3];
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
    
    
    printf("points=%d \n", points);
    
    MatrixD2D out(points,3);
    for (int i = 0; i<points; i++)
    {
        out(i,0) = brain[i*3    ];
        out(i,1) = brain[i*3 + 1];
        out(i,2) = brain[i*3 + 2];
    }
    printf("radius = %f \n", radius);
    cout<<"radius ="<< radius<<endl;
    char filename2[256];
    sprintf(filename2,"Sphere.dat");
    out.dump(filename2);
    
}



void Glioma_HG_UQ::run()
{
    
    bool bProfiler = 0;
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* read in patien specific parameters*/
    ifstream mydata("HGG_InputParameters.txt");
    Real Dg, Dw, rho;
    double tend;
    
    if (mydata.is_open())
    {
        mydata >> Dw;
        mydata >> rho;
        mydata >> tend;
        mydata.close();
    }
    
    /*rescale*/
    Dw = Dw/(L*L);
    Dg = 0.1*Dw;
    
    double t			= 0.0;
    int iCounter        = 1;
    
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.99 * h*h / ( 2.* _DIM * max(Dw, Dg) );
    if(bVerbose)  printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho,h);
    

    
    while (t <= tend)
    {
        if(bProfiler) profiler.getAgent("RD_Step").start();
        _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
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
            
           // _dump(iCounter);
            iCounter++;
            whenToWrite = whenToWrite + whenToWriteOffset;
            if(bVerbose) printf("Dumping data at time t=%f\n", t);
            
        }
    }
    
    
    // _dumpSubBrainPoints(*grid);
    
    // Refine final state & dump for UQ Likelihood
    if(bAdaptivity)
        Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
    
    
    if(bProfiler) profiler.getAgent("UQ_output").start();
    _dumpUQoutput(*grid);
    if(bProfiler) profiler.getAgent("UQ_output").stop();
    
    _dump(1000);
    
    if(bVerbose) profiler.printSummary();
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
