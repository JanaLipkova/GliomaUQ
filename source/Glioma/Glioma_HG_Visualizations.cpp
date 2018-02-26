//
//  Glioma_HG_Visualizations.cpp
//  GliomaBrutusXcode
//
//  Created by Lipkova on 21/01/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#include "Glioma_HG_Visualizations.h"

// need biger stencil for the refinment !!!
static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_HG_Visualizations::Glioma_HG_Visualizations(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("/////////Volumes/FileStorage/GLIOMA/DataForPaper/10/VOIs_pat10//////////             PROPAGATION TOOL                 ////////////////\n");
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

Glioma_HG_Visualizations::~Glioma_HG_Visualizations()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
void Glioma_HG_Visualizations::_ic(Grid<W,B>& grid, int pID)
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
    
    sprintf(anatomy, "%s_MRI_T1w.dat", patientFolder);
    MatrixD3D T1w(anatomy);
    sprintf(anatomy, "%s_MRI_FLAIR.dat", patientFolder);
    MatrixD3D FLAIR(anatomy);
    sprintf(anatomy, "%s_FET1.dat", patientFolder);
    MatrixD3D PET1(anatomy);
    sprintf(anatomy, "%s_FET2.dat", patientFolder);
    MatrixD3D PET2(anatomy);
    sprintf(anatomy, "%s_FET3.dat", patientFolder);
    MatrixD3D PET3(anatomy);
    sprintf(anatomy, "%s_T1.dat", patientFolder);
    MatrixD3D T1s(anatomy);
    sprintf(anatomy, "%s_T2.dat", patientFolder);
    MatrixD3D T2s(anatomy);
    sprintf(anatomy, "%s_distMap.dat", patientFolder);
    MatrixD3D RT(anatomy);
    sprintf(anatomy, "%s_R_FLAIR_reg.dat", patientFolder);
    MatrixD3D recFLAIR(anatomy);
    sprintf(anatomy, "%s_R_T1w_reg.dat", patientFolder);
    MatrixD3D recT1w(anatomy);
    sprintf(anatomy, "%s_R_FLAIR_tumour.dat", patientFolder);
    MatrixD3D recT2s(anatomy);
    sprintf(anatomy, "%s_R_T1w_tumour.dat", patientFolder);
    MatrixD3D recT1s(anatomy);
    

    
    int brainSizeX = (int) T1w.getSizeX();
    int brainSizeY = (int) T1w.getSizeY();
    int brainSizeZ = (int) T1w.getSizeZ();
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L = brainSizeMax * 0.1; // 25.6
    
    printf("L=%f \n", L);
    std::cout<<"T1w:      "<<brainSizeX                 <<","<<brainSizeY                   <<","<<brainSizeZ                   <<std::endl;
    std::cout<<"FLAIR:    "<<(int) FLAIR.getSizeX()     <<","<<(int) FLAIR.getSizeY()       <<","<<(int) FLAIR.getSizeZ()       <<std::endl;
    std::cout<<"PET1:     "<<(int) PET1.getSizeX()      <<","<<(int) PET1.getSizeY()        <<","<<(int) PET1.getSizeZ()        <<std::endl;
    std::cout<<"PET2:     "<<(int) PET2.getSizeX()      <<","<<(int) PET2.getSizeY()        <<","<<(int) PET2.getSizeZ()        <<std::endl;
    std::cout<<"PET3:     "<<(int) PET3.getSizeX()      <<","<<(int) PET3.getSizeY()        <<","<<(int) PET3.getSizeZ()        <<std::endl;
    std::cout<<"T1s:      "<<(int) T1s.getSizeX()       <<","<<(int) T1s.getSizeY()         <<","<<(int) T1s.getSizeZ()         <<std::endl;
    std::cout<<"T2s:      "<<(int) T2s.getSizeX()       <<","<<(int) T2s.getSizeY()         <<","<<(int) T2s.getSizeZ()         <<std::endl;
    std::cout<<"RT:       "<<(int) RT.getSizeX()        <<","<<(int) RT.getSizeY()          <<","<<(int) RT.getSizeZ()          <<std::endl;
    std::cout<<"recFLAIR: "<<(int) recFLAIR.getSizeX()  <<","<<(int) recFLAIR.getSizeY()    <<","<<(int) recFLAIR.getSizeZ()    <<std::endl;
    std::cout<<"recT1w:   "<<(int) recT1w.getSizeX()    <<","<<(int) recT1w.getSizeY()      <<","<<(int) recT1w.getSizeZ()      <<std::endl;
    std::cout<<"recT2s:   "<<(int) recT2s.getSizeX()    <<","<<(int) recT2s.getSizeY()      <<","<<(int) recT2s.getSizeZ()      <<std::endl;
    std::cout<<"recT1s:   "<<(int) recT1s.getSizeX()        <<","<<(int) recT1s.getSizeY()      <<","<<(int) recT1s.getSizeZ()      <<std::endl;

    
    
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
                    
                    
                    if ( (mappedBrainX > 0 && mappedBrainX < brainSizeX-1) && (mappedBrainY > 0 && mappedBrainY-1 < brainSizeY-1) && (mappedBrainZ > 0 && mappedBrainZ < brainSizeZ-1) )
                    {
                        
                        block(ix,iy,iz).p_w     = T1w(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).p_g     = FLAIR(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).p_csf   = PET1(mappedBrainX,mappedBrainY,mappedBrainZ) + PET2(mappedBrainX,mappedBrainY,mappedBrainZ) + PET3(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).t1bc    = T1s(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).t2bc    = T2s(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).eps     = RT(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).mu      = recFLAIR(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).kappa   = recT1w(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).chi     = recT2s(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).psi     = recT1s(mappedBrainX,mappedBrainY,mappedBrainZ);

                    }
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
    
    
}




#pragma mark DumpingOutput
void Glioma_HG_Visualizations:: _dump(int counter)
{
    if(bVerbose) printf("dumping data \n");
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"Visualisation%04i", counter);
        
        IO_VTKNative3D<W,B, 10,0 > vtkdumper;
        vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
    }
    
}

void Glioma_HG_Visualizations::_dumpBinaryOutput(Grid<W,B>& grid)
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    
    if(bVerbose) printf("bpd=%i, bs=%i, hf=%f  \n",blocksPerDimension,blockSize,hf);
    
    MatrixD3D flair(gpd,gpd,gpd);
    MatrixD3D t1(gpd,gpd,gpd);
    MatrixD3D pet(gpd,gpd,gpd);

    
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
                    
                    flair(gx,gy,gz) = block(ix,iy,iz).p_g;
                    t1(gx,gy,gz)    = block(ix,iy,iz).p_w;
                    pet(gx,gy,gz)   = block(ix,iy,iz).p_csf;
                }
        
    }
    
    char filename[256];
    sprintf(filename,"FLAIR.dat");
    flair.dump(filename);
    
    sprintf(filename,"T1.dat");
    t1.dump(filename);
    
    sprintf(filename,"PET.dat");
    pet.dump(filename);
    
}



void Glioma_HG_Visualizations::run()
{
    int size = blockSize * blocksPerDimension;
    _dump(size);
    //_dumpBinaryOutput(*grid);
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}