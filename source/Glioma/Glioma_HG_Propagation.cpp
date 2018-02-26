//
//  Glioma_HG_Propagation.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 11/01/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_HG_Propagation.h"


// need biger stencil for the refinment !!!
static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_HG_Propagation::Glioma_HG_Propagation(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////             PROPAGATION TOOL                 ////////////////\n");
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
    
    
    int ICtype = parser("-IC").asInt(0);
    
    switch (ICtype)
    {
        case 0:
            _ic_SubjectBrainPropagation(*grid);
            break;
        case 1:
            pID =  parser("-pID").asInt();
            _ic_PatientCasePropagation(*grid, pID);
            break;
            
        default:
            break;
    }
    
    isDone              = false;
}

Glioma_HG_Propagation::~Glioma_HG_Propagation()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
// Patient Brain anatomy - subject04 from BrainWeb database
void Glioma_HG_Propagation:: _ic_SubjectBrainPropagation(Grid<W,B>& grid)
{
    
#pragma mark Initialization
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    /* 2) Read in Propagation Results  */
    /* 2.1. Mean */
    int nSamples = _readInNumberOfSamples();
    char filename[256];

    for (int is = 0; is < nSamples; is++)
    {
        
#ifdef BRUTUS
        sprintf(filename, "/cluster/scratch_xp/public/lipkovaj/TMCMC/SpaceSearching/SynthethicAll_All/SynthethicAll_All_4K_Lustre_C/Propagations/tmpdir.0.0.0.%d/HGG_data.dat", is);
#elif defined(LRZ_CLUSTER)
        sprintf(filename, "/gpfs/scratch/txh01/di49zin/HGG_UQ/SyntheticBig_P103/SyntheticBig_P103_ALL_4K_newPrior/tmpdir.0.0.0.%d/HGG_data.dat", is);

        
        
//        sprintf(filename, "/naslx/ptmp/10/di49zin/HGG/Synthetic/Synthetic_MRIonly_4K_HR/tmpdir.0.0.0.%d/HGG_data.dat", is);
#endif
        
//        printf("Reading data %s \n", filename);
        
        MatrixD3D data(filename);
        
        int dataSizeX = (int) data.getSizeX();
        int dataSizeY = (int) data.getSizeY();
        int dataSizeZ = (int) data.getSizeZ();
        int dataSizeMax = max(dataSizeX, max(dataSizeY,dataSizeZ));
        
        assert(dataSizeX == dataSizeY);
        assert(dataSizeY == dataSizeZ);
        
        double dataHx = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHy = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHz = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        
        
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
                        
                        int mappedDataX = (int)floor( x[0] / dataHx  );
                        int mappedDataY = (int)floor( x[1] / dataHy  );
                        int mappedDataZ = (int)floor( x[2] / dataHz  );
                        
                        
                        if ( (mappedDataX >= 0 && mappedDataX < dataSizeX) && (mappedDataY >= 0 && mappedDataY < dataSizeY) && (mappedDataZ >= 0 && mappedDataZ < dataSizeZ) )
                            block(ix,iy,iz).mean += data(mappedDataX,mappedDataY,mappedDataZ);
                        
                    }
            
            grid.getBlockCollection().release(info.blockID);
            
        }
    }
    
    /* 2.2 var = (u - mean)^2 */
    Real inSamples = 1./nSamples;
    
    for (int is = 0; is < nSamples; is++)
    {
#ifdef BRUTUS
        sprintf(filename, "/cluster/scratch_xp/public/lipkovaj/TMCMC/SpaceSearching/SynthethicAll_All/SynthethicAll_All_4K_Lustre_C/Propagations/tmpdir.0.0.0.%d/HGG_data.dat", is);
#elif defined(LRZ_CLUSTER)
        
        sprintf(filename, "/gpfs/scratch/txh01/di49zin/HGG_UQ/SyntheticBig_P103/SyntheticBig_P103_ALL_4K_newPrior/tmpdir.0.0.0.%d/HGG_data.dat", is);
        
//        sprintf(filename, "/naslx/ptmp/10/di49zin/HGG/Synthetic/Synthetic_MRIonly_4K_HR/tmpdir.0.0.0.%d/HGG_data.dat", is);
#endif
//        printf("Reading data %s \n", filename);

        MatrixD3D data(filename);
        
        int dataSizeX = (int) data.getSizeX();
        int dataSizeY = (int) data.getSizeY();
        int dataSizeZ = (int) data.getSizeZ();
        int dataSizeMax = max(dataSizeX, max(dataSizeY,dataSizeZ));
        
        assert(dataSizeX == dataSizeY);
        assert(dataSizeY == dataSizeZ);
        
        double dataHx = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHy = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHz = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        
        
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
                        
                        int mappedDataX = (int)floor( x[0] / dataHx  );
                        int mappedDataY = (int)floor( x[1] / dataHy  );
                        int mappedDataZ = (int)floor( x[2] / dataHz  );
                        
                        
                        if ( (mappedDataX > 0 && mappedDataX <= dataSizeX) && (mappedDataY > 0 && mappedDataY <= dataSizeY) && (mappedDataZ > 0 && mappedDataZ <= dataSizeZ) )
                        {
                            Real tmp1 = data(mappedDataX,mappedDataY,mappedDataZ) - inSamples * block(ix,iy,iz).mean;
                            block(ix,iy,iz).var += tmp1 * tmp1;
                        }
                    }
            
            grid.getBlockCollection().release(info.blockID);
            
        }
    }
    
}

void Glioma_HG_Propagation::_ic_PatientCasePropagation(Grid<W,B>& grid,int pID)
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
                    
                    // remove fluid below tumour that would be normally pushed away by growing tumour
                    Real MRIsignal = PT1 + PT2;
                    Pcsf  = ( MRIsignal > 0.) ? 0. : Pcsf;
                    
                    // Anatomy
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
        
        grid.getBlockCollection().release(info.blockID);
        
    }
    
    
    
    /* 2) Read in Propagation Results  */
    int nSamples = _readInNumberOfSamples();
    char filename[256];

#ifdef BRUTUS
    sprintf(dataFolder,"/cluster/scratch_xp/public/lipkovaj/TMCMC/PatientCases/");
#elif defined(KRAKEN)
    sprintf(dataFolder,"/home/jana/Work/GliomaAdvance/source/Anatmoy/");
#elif defined(PLURIPOTENT)
    sprintf(dataFolder,"/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/");
#elif defined(LRZ_CLUSTER)
    
    sprintf(dataFolder,"/gpfs/scratch/txh01/di49zin/HGG_UQ/");  // mpp2
//    sprintf(dataFolder,"/naslx/ptmp/10/di49zin/HGG/");        // mpp1
#else
    sprintf(dataFolder,"../../Anatmoy/");
#endif
    
    
#ifdef Synthetic
    sprintf(patientFolder, "%sSyntheticBig_P103/SyntheticBig_P103_ALL_4K_newPrior_LR",dataFolder);
    //    sprintf(patientFolder, "%sSyntheticBig_P103/SyntheticBig_P103_ALL_4K_newPrior",dataFolder);
#else
//    sprintf(patientFolder, "%sPatient%02d/P%02d_MRIonly_4K",dataFolder,pID,pID);
    sprintf(patientFolder, "%sPatient%02d/Patient%02d_All_4K_Ti_step_3_mpp1",dataFolder,pID,pID);

#endif

    
    /* 2.1. Mean = 1/N sum u(ix,iy,iz,is),  1st part of variance: var = 1/N sum u(ix,iy,iz,is)^2 */
    for (int is = 0; is < nSamples; is++)
    {
        sprintf(filename, "%s/tmpdir.0.0.0.%d/HGG_data.dat",patientFolder,is);
        MatrixD3D data(filename);
        
        int dataSizeX = (int) data.getSizeX();
        int dataSizeY = (int) data.getSizeY();
        int dataSizeZ = (int) data.getSizeZ();
        int dataSizeMax = max(dataSizeX, max(dataSizeY,dataSizeZ));
        
        assert(dataSizeX == dataSizeY);
        assert(dataSizeY == dataSizeZ);
        
        double dataHx = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHy = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHz = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        
        
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
                        
                        int mappedDataX = (int)floor( x[0] / dataHx  );
                        int mappedDataY = (int)floor( x[1] / dataHy  );
                        int mappedDataZ = (int)floor( x[2] / dataHz  );
                        
                        
                        if ( (mappedDataX > 0 && mappedDataX <= dataSizeX) && (mappedDataY > 0 && mappedDataY <= dataSizeY) && (mappedDataZ > 0 && mappedDataZ <= dataSizeZ) )
                            block(ix,iy,iz).mean += data(mappedDataX,mappedDataY,mappedDataZ);
                        
                    }
            
            grid.getBlockCollection().release(info.blockID);
            
        }
    }
    
    /* 2.2 2nd part of variance */
    Real inSamples = 1./nSamples;
    
    for (int is = 0; is < nSamples; is++)
    {
        sprintf(filename, "%s/tmpdir.0.0.0.%d/HGG_data.dat",patientFolder,is);
        MatrixD3D data(filename);
        
        int dataSizeX = (int) data.getSizeX();
        int dataSizeY = (int) data.getSizeY();
        int dataSizeZ = (int) data.getSizeZ();
        int dataSizeMax = max(dataSizeX, max(dataSizeY,dataSizeZ));
        
        
//        printf("DataSizeX = %i, DataSizeY = %i, DataSizeZ = %i \n", dataSizeX, dataSizeY, dataSizeZ);
        
        assert(dataSizeX == dataSizeY);
        assert(dataSizeY == dataSizeZ);
        
        double dataHx = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHy = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        double dataHz = 1.0 / ((double)(dataSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
        
        
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
                        
                        int mappedDataX = (int)floor( x[0] / dataHx  );
                        int mappedDataY = (int)floor( x[1] / dataHy  );
                        int mappedDataZ = (int)floor( x[2] / dataHz  );
                        
                        
                        if ( (mappedDataX > 0 && mappedDataX <= dataSizeX) && (mappedDataY > 0 && mappedDataY <= dataSizeY) && (mappedDataZ > 0 && mappedDataZ <= dataSizeZ) )
                        {
                            Real tmp1 = data(mappedDataX,mappedDataY,mappedDataZ) - inSamples * block(ix,iy,iz).mean;
                            block(ix,iy,iz).var += tmp1 * tmp1;
                        }
                    }
            
            grid.getBlockCollection().release(info.blockID);
            
        }
    }
    
    
}




void Glioma_HG_Propagation::_selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize)
{
    FILE * fp;
    
    
#pragma mark GRAY_MATTER
    // cout << "Reading in Gray Matter: " << endl;
    
#ifdef BRUTUS
    fp = fopen("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_gm_v.rawb", "rb");
#elif defined(LRZ_CLUSTER)
    fp = fopen("/home/hpc/txh01/di49zin/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_gm_v.rawb", "rb");
#else
    fp = fopen("../Anatmoy/Subject42/subject42_gm_v.rawb", "rb");
#endif
    
    _readInBrainWebAnatomy(GreyTissueData,fp, dataSize, 0 );
    fclose (fp);
    
    
#pragma mark WHITE_MATTER
    //cout << "Reading in White Matter: " << endl;
    
#ifdef BRUTUS
    fp = fopen("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_wm_v.rawb", "rb");
#elif defined(LRZ_CLUSTER)
    fp = fopen("/home/hpc/txh01/di49zin/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_wm_v.rawb", "rb");
#else
    fp = fopen("../Anatmoy/Subject42/subject42_wm_v.rawb", "rb");
#endif
    
    _readInBrainWebAnatomy(WhiteTissueData,fp, dataSize, 0 );
    fclose (fp);
    
#pragma mark VESSELS
    // cout << "Reading in Vessels: " << endl;
    
#ifdef BRUTUS
    fp = fopen("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_vessels.rawb", "rb");
#elif defined(LRZ_CLUSTER)
    fp = fopen("/home/hpc/txh01/di49zin/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_vessels.rawb", "rb");
#else
    fp = fopen("../Anatmoy/Subject42/subject42_vessels.rawb", "rb");
#endif
    
    _readInBrainWebAnatomy(CsfData,fp, dataSize, 0 );
    fclose (fp);
    
    
#pragma mark CSF
    // cout << "Reading in CSF:  " << endl;
    
#ifdef BRUTUS
    fp = fopen("/cluster/home/mavt/lipkovaj/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_csf_v.rawb", "rb");
#elif defined(LRZ_CLUSTER)
    fp = fopen("/home/hpc/txh01/di49zin/GliomaAdvance/UQ_Section/source/Anatmoy/Subject42/subject42_csf_v.rawb", "rb");
#else
    fp = fopen("../Anatmoy/Subject42/subject42_csf_v.rawb", "rb");
#endif
    
    _readInBrainWebAnatomy(CsfData,fp, dataSize, 128 );
    fclose (fp);
}

void Glioma_HG_Propagation::_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold )
{
    typedef uint8_t brainWebType;
    
    // obtain file size
    if (fp == NULL) {fputs ("File error", stderr); exit (1);}
    fseek (fp , 0 , SEEK_END);
    long int size = ftell (fp);
    rewind (fp);
    
    assert(size == DataSize );
    
    // allocate memory to contain the whole file:
    brainWebType * buffer;
    buffer = (brainWebType*) malloc (sizeof(brainWebType)*size);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
    
    // copy the file into the buffer:
    size_t result;
    result = fread (buffer, 1, size, fp);
    if (result != size) {fputs ("Reading error",stderr); exit (3);}
    
    for (int i = 0; i < size; ++i)
    {
        if (tissue[i] >= threshold)
            tissue[i] = (float)buffer[i];
    }
    
    free(buffer);
}

void Glioma_HG_Propagation::_readInTumorPosition(vector<Real>& tumorIC )
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
        tumorIC[i] = (Real)buffer[i];
    
    
    free(buffer);
    fclose (fp);
}


int Glioma_HG_Propagation::_readInNumberOfSamples()
{
    typedef float dataType;
    
    FILE* fp;
    fp = fopen("HGG_Nsamples.bin", "rb");
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
    
    int nSamples = (int)buffer[0];
    
    free(buffer);
    fclose (fp);
    
    return nSamples;
}



void Glioma_HG_Propagation::_computeStatistics(const int nParallelGranularity, const int nSamples)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    PropagationStatistics  statistics(nSamples);
    BlockProcessing::process(vInfo, collecton, statistics, nParallelGranularity);
}




#pragma mark DumpingOutput
void Glioma_HG_Propagation:: _dump(int counter)
{
    
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"Propagation%04i", counter);
        
        IO_VTKNative3D<W,B, 4,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
    
}

void Glioma_HG_Propagation::_dumpBinaryOutput(Grid<W,B>& grid, int nSamples)
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    if(bVerbose) printf("bpd=%i, bs=%i, hf=%f, nSamples=%i \n",blocksPerDimension,blockSize,hf,nSamples);
    
    MatrixD3D meanOutput(gpd,gpd,gpd);
    MatrixD3D varOutput(gpd,gpd,gpd);
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp paralle for
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
                    
                    //mapped coordinates
                    int gx = (int)floor( (x[0]) / hf  );
                    int gy = (int)floor( (x[1]) / hf  );
                    int gz = (int)floor( (x[2]) / hf  );
                    
                    meanOutput(gx,gy,gz) = block(ix,iy,iz).mean;
                    varOutput(gx,gy,gz)  = block(ix,iy,iz).var;
                    
                }
        
    }
    
    char filename[256];
    sprintf(filename,"PropagationMean_%i.dat",nSamples);
    meanOutput.dump(filename);
    
    sprintf(filename,"PropagationVar_%i.dat",nSamples);
    varOutput.dump(filename);
    
}

void Glioma_HG_Propagation::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    
    int nSamples = _readInNumberOfSamples();
    _computeStatistics(nParallelGranularity,nSamples);
    
    _dump(nSamples);
    _dumpBinaryOutput(*grid, nSamples);
    
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
