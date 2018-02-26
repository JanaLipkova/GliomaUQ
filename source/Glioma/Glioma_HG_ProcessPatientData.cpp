//
//  Glioma_HG_ProcessPatientData.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 13/11/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_HG_ProcessPatientData.h"

static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_HG_ProcessPatientData::Glioma_HG_ProcessPatientData(int argc, const char ** argv): parser(argc, argv)
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
    
    pID =  parser("-pID").asInt();
    _icPatientData(*grid, pID);
    
    _dump(0);
    isDone              = false;
}

Glioma_HG_ProcessPatientData::~Glioma_HG_ProcessPatientData()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* 1) Read in all patient data
 2) Rescaled into desired simulation domain
 3) Integrate time out of PET signals
 4) Restric PET signal only to T1 and T2 region */
void Glioma_HG_ProcessPatientData:: _icPatientData(Grid<W,B>& grid, int pID)
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
    sprintf(anatomy, "%s_FET1.dat", patientFolder);
    MatrixD3D PET1(anatomy);
    sprintf(anatomy, "%s_FET2.dat", patientFolder);
    MatrixD3D PET2(anatomy);
    sprintf(anatomy, "%s_FET3.dat", patientFolder);
    MatrixD3D PET3(anatomy);
    
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L                = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm
    
    printf("Data sizeX=%i, sizeY=%i, sizeZ=%i\n", brainSizeX, brainSizeY, brainSizeZ);
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    
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
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )
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
                        
                        
                        
                        // Integrate time out of PET signal, and restrict PET into T1 u T2 region
                        Real PETsignal = Ppet1 + Ppet2 + Ppet3;
                        Real MRIsignal = PT1 + PT2;
                        block(ix,iy,iz).phi = (MRIsignal > 0.9) ? PETsignal : 0.;
                        block(ix,iy,iz).t1bc = (PT1 > 0.5) ? 1. : 0.; // T1 segm, 1 - active T1, 2 - necrotic core
                        block(ix,iy,iz).t2bc = (block(ix,iy,iz).t1bc > 0.01) ? 1. : PT2;
                        
                        // remove fluid that would be pushed away by tumour
                        Pcsf  = ( block(ix,iy,iz).phi > 0.) ? 0. : Pcsf;

                        double all = PWt + PGt + Pcsf;
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
                        
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


void Glioma_HG_ProcessPatientData::_getPETstatistic()
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
    
    FILE * pFile;
    float buffer[3] = { cx , cy , cz };
    pFile = fopen ("HGG_TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);
    fclose (pFile);
}



void Glioma_HG_ProcessPatientData::_normalisePETsignal()
{
    double iscale = 1./maxPET;
    double PETthr = parser("-PETthr").asDouble(0.); //threashold to exclude necrotic regions
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    block(ix,iy,iz).phi = block(ix,iy,iz).phi * iscale;
                    block(ix,iy,iz).phi = ( block(ix,iy,iz).phi < PETthr) ? 0. : block(ix,iy,iz).phi;
                    
                }
    }
}



#pragma mark DumpingOutput
void Glioma_HG_ProcessPatientData:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%dD_Patient42_Data%04d",_DIM, counter);
        
        
        IO_VTKNative3D<W,B, 9,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        
    }
    
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_HG_ProcessPatientData::_dumpOutput()
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



void Glioma_HG_ProcessPatientData::_readInTumorPosition(vector<Real>& tumorIC )
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
    
    
    
    printf("in read sec cm: cm[0]=%f, cm[1]=%f, cm[2]=%f \n",tumorIC[0],tumorIC[1],tumorIC[2]);
    
    free(buffer);
    fclose (fp);
}

void Glioma_HG_ProcessPatientData::_dumpSubBrainPoints(Grid<W,B>& grid )
{
    if (bAdaptivity)
    {
        printf("Aborting ... dumpBrainPoints needs uniform grid");
        abort();
    }
    
    int points   = 0;
    vector<Real> cm(3);
    _readInTumorPosition(cm);
    
    
    printf("cm: cm[0]=%f, cm[1]=%f, cm[2]=%f \n", cm[0],cm[1],cm[2]);
    
    float radius = parser("-radius").asDouble();  //made it max_tumor radius + 2cm margin
    
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
    
    
    printf("points=%d \n", points);
    
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



void Glioma_HG_ProcessPatientData::run()
{
    _getPETstatistic();
    _normalisePETsignal();
    _dump(1);
    _dumpOutput();
    _dumpSubBrainPoints(*grid);
    
    profiler.printSummary();
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}
