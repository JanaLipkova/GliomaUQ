//
//  Glioma_dat2VP.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 19/03/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_dat2VP.h"



// need biger stencil for the refinment !!!
static int maxStencil[2][3] = {
    -2, -2, -2,
    +3, +3, +3
};

Glioma_dat2VP::Glioma_dat2VP(int argc, const char ** argv): parser(argc, argv)
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
    
    pID =  parser("-pID").asInt();
    _ic(*grid, pID);
    
    isDone              = false;
}

Glioma_dat2VP::~Glioma_dat2VP()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
void Glioma_dat2VP::_ic(Grid<W,B>& grid, int pID)
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
    sprintf(anatomy, "%s_MAP.dat", patientFolder);
    MatrixD3D MAP(anatomy);
    
    /* Read in binary matrices */
    
    // Tumour
    int brainSizeX = (int) MAP.getSizeX();
    int brainSizeY = (int) MAP.getSizeY();
    int brainSizeZ = (int) MAP.getSizeZ();
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L = brainSizeMax * 0.1; // 25.6
    
    printf("L=%f \n", L);
    std::cout<<"MAP: "<<brainSizeX<<","<<brainSizeY<<","<<brainSizeZ<<std::endl;
    std::cout<<"T1: "<<(int) T1.getSizeX() <<","<<(int) T1.getSizeY() <<","<<(int) T1.getSizeZ() <<std::endl;
    std::cout<<"PET1: "<<(int) PET1.getSizeX() <<","<<(int) PET1.getSizeY() <<","<<(int) PET1.getSizeZ() <<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    // Segmentations
    int segmSizeX = (int) T1.getSizeX();
    int segmSizeY = (int) T1.getSizeY();
    int segmSizeZ = (int) T1.getSizeZ();
    int segmSizeMax = max(segmSizeX, max(segmSizeY,segmSizeZ));
    
    double segmHx = 1.0 / ((double)(segmSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double segmHy = 1.0 / ((double)(segmSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double segmHz = 1.0 / ((double)(segmSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    vector<Real> blockMaxPET;
    
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
                    
                    
                    // switch  x and y,  mirror x and y (P10, P16, P18)
//                    Real tmp = x[1];
//                    x[1] = 1.-x[0];
//                    x[0] = 1.-tmp;
                    
                    // mirror y: direction for VP (so frontal lobe is at the top) P01,P07,P20,P22,P31
                    x[1] = 1. - x[1];
                    
                    
                    
                   
                    

                    
                    /* 1) Tumour data */
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    
                    if ( (mappedBrainX > 0 && mappedBrainX < brainSizeX-1) && (mappedBrainY > 0 && mappedBrainY-1 < brainSizeY-1) && (mappedBrainZ > 0 && mappedBrainZ < brainSizeZ-1) )
                        block(ix,iy,iz).phi = MAP(mappedBrainX,mappedBrainY,mappedBrainZ);

                    
                  
                    /* 2) Segmenation data */
                    int mappedSegmX = (int)floor( x[0] / segmHx  );
                    int mappedSegmY = (int)floor( x[1] / segmHy  );
                    int mappedSegmZ = (int)floor( x[2] / segmHz  );
                    
                    // aspect ratio correction
                    mappedSegmX -= (int) ( (segmSizeMax - segmSizeX) * 0.5);
                    mappedSegmY -= (int) ( (segmSizeMax - segmSizeY) * 0.5);
                    mappedSegmZ -= (int) ( (segmSizeMax - segmSizeZ) * 0.5);
                    
                    
                    if ( (mappedSegmX > 0 && mappedSegmX < segmSizeX-1) && (mappedSegmY > 0 && mappedSegmY-1 < segmSizeY-1) && (mappedSegmZ > 0 && mappedSegmZ < segmSizeZ-1) )
                    {
                        Real PT1     = T1(mappedSegmX,mappedSegmY,mappedSegmZ);
                        Real PT2     = T2(mappedSegmX,mappedSegmY,mappedSegmZ);
                        Real Ppet1   = PET1(mappedSegmX,mappedSegmY,mappedSegmZ);
                        Real Ppet2   = PET2(mappedSegmX,mappedSegmY,mappedSegmZ);
                        Real Ppet3   = PET3(mappedSegmX,mappedSegmY,mappedSegmZ);
                        Real pMAP    = MAP(mappedSegmX,mappedSegmY,mappedSegmZ);
                        
                        // Integrate time out of PET signal, and restrict PET into T1 u T2 region
                        Real MRIsignal = PT1 + PT2;

                        Real PETsignal = Ppet1 + Ppet2 + Ppet3;
                        block(ix,iy,iz).psi = (MRIsignal > 0.9) ? PETsignal : 0.;
                        block(ix,iy,iz).p_g = GM(mappedSegmX,mappedSegmY,mappedSegmZ);
                        block(ix,iy,iz).p_w = WM(mappedSegmX,mappedSegmY,mappedSegmZ);
                        
                        
                    }
                    

                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


void Glioma_dat2VP:: _normalisePET(Grid<W,B>& grid)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    vector<Real> blockMaxPET(vInfo.size());  // fill with zeros
    
    
    // 1) get local scaling factor
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        Real localMaxPET = 0.;
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    localMaxPET = max( localMaxPET, block(ix,iy,iz).psi);
        
        blockMaxPET[i] = localMaxPET;
    }
    
    // 2) get global scaling factor
    Real maxPET = 0.;
    
    for(int i=0; i<vInfo.size(); i++)
        maxPET = max(maxPET,blockMaxPET[i] );
    
    printf("PET scaling factor: \%f \n", maxPET);
    
    // 3) normalise
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).psi =  block(ix,iy,iz).psi/maxPET;
        
    }
    
}

// 0 - tumour
// 1 - pet
// 2 - anatomy
// 3 - tumour + pet
// 4 - tumour + pet + anatomy

void Glioma_dat2VP:: _copyFields(Grid<W,B>& grid, int filedID)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    switch (filedID) {
                        case 0:
                            block(ix,iy,iz).vp = block(ix,iy,iz).phi;
                            break;
                            
                        case 1:
                            block(ix,iy,iz).vp = block(ix,iy,iz).psi;
                            break;
                            
                        case 2:
                            block(ix,iy,iz).vp = (block(ix,iy,iz).p_g + block(ix,iy,iz).p_w > 1.) ? 1. : block(ix,iy,iz).p_g + block(ix,iy,iz).p_w ;
                            break;
                            
                        case 3:
                        {
                            // threashold + combine: tum = [0.001, 0.3], pet = [0.3, 1]
                            Real tum = (block(ix,iy,iz).phi >= 0.3) ? 0.3  : block(ix,iy,iz).phi;
                            Real pet = (block(ix,iy,iz).psi <  0.3) ? 0.   : block(ix,iy,iz).psi;
                            
                            block(ix,iy,iz).vp = (tum + pet > 1.) ? 1. : tum + pet;
                        }
                            break;
                        case 4:
                        {
                            // threashold + combine: anat = 0.2, tum = [0.201, 0.5], pet = [0.5, 1]
                            Real anatomy = (block(ix,iy,iz).p_w + block(ix,iy,iz).p_g > 0.1) ? 1. : 0.;  // anatomy mask
                            Real tum = (block(ix,iy,iz).phi >= 0.3) ? 0.3  : block(ix,iy,iz).phi;
                            Real pet = (block(ix,iy,iz).psi <  0.3) ? 0.   : block(ix,iy,iz).psi;
                            
                            Real out = 0.2 * anatomy + tum + pet;
                            block(ix,iy,iz).vp = ( out > 1.) ? 1. : out;
                        }
                            break;
                        case 5:
                        {
                            // threashold + combine: anat = 0.2, tum = [0.201, 0.5], pet = [0.5, 1]
                            Real anatomy = (block(ix,iy,iz).p_w + block(ix,iy,iz).p_g > 0.1) ? 1. : 0.;  // anatomy mask
                            Real tum = (block(ix,iy,iz).phi >= 0.2) ? 0.2  : block(ix,iy,iz).phi;
                            Real pet = (block(ix,iy,iz).psi <  0.2) ? 0.   : block(ix,iy,iz).psi;
                            
                            Real out = 0.01 * anatomy + 2. * tum + pet;
                            block(ix,iy,iz).vp = ( out > 1.) ? 1. : out;
                        }
                            break;
                            
                            
                        default:
                            break;
                    }
                }
    }
    
}


void Glioma_dat2VP:: _combineFields(Grid<W,B>& grid)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
         
                    Real tum = (block(ix,iy,iz).phi >= 0.2) ? 0.2  : block(ix,iy,iz).phi;
                    Real pet = (block(ix,iy,iz).psi <  0.2) ? 0.   : block(ix,iy,iz).psi;
                    
                    block(ix,iy,iz).vp = (tum + pet > 1.) ? 1. : tum + pet;
                }

    }
}



#pragma mark DumpingOutput
void Glioma_dat2VP:: _dumpVTK(int counter)
{
    if(bVerbose) printf("dumping VTK data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"Visualisation%04i", counter);
        
        IO_VTKNative3D<W,B, 5,0 > vtkdumper;
        vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}

void Glioma_dat2VP:: _dumpVP(int counter)
{
    if(bVerbose) printf("dumping VP data \n");
    if (parser("-vp").asBool())
    {
        char filename[256];
        sprintf(filename,"Visualisation%04i", counter);
        
       DumpScalarToVP < BlockLab<B>, W, B>	vpdumper(grid);
        vpdumper.Write(filename);
    }


    
    
}

void Glioma_dat2VP::run()
{
    
    _normalisePET(*grid);
    
//    int inputID = 0;
//    _dumpVTK(inputID);  // tumor
//    _dumpVP(inputID);
//    inputID++;
//    
//    _copyFields(*grid, 3); // tum + pet
//    _dumpVTK(inputID);
//    _dumpVP(inputID);
//    inputID++;
    
    int inputID = 0;
    _copyFields(*grid, 5); // tum + anatomy
    _dumpVTK(inputID);
    _dumpVP(inputID);
    
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
