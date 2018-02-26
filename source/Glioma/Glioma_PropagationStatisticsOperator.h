//
//  Glioma_HG_PropagationStatisticsOperator.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 13/01/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef Glioma_PropagationStatisticsOperator_h
#define Glioma_PropagationStatisticsOperator_h



/* 
 STATE OF VARIABLES AT INPUT TIME:
 block(ix,iy,iz).mean = sum_is u(ix,iy,iz,is)  -> needs to be rescaled by nSamples
 block(ix,iy,iz).var  = sum_is u(ix,iy,iz,is)^2 
 block(ix,iy,iz).tmp  = sum_is (u(ix,iy,iz,is) - mean)
 
 OUTPUT:
 block(ix,iy,iz).mean /= nSamples
 block(ix,iy,iz).var = (block(ix,iy,ix).var + block(ix,iy,iz).tmp) /nSamples 
 block(ix,iy,iz).tmp /= nSamples */
struct PropagationStatistics
{
    const int nSamples;
    
    PropagationStatistics(int nSamples_):nSamples(nSamples_)
    { }
    
    PropagationStatistics(const PropagationStatistics& copy): nSamples(copy.nSamples)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        Real inSamples = 1./nSamples;
        
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix,iy,iz).mean *= inSamples;  // 1st moment - mean
                    o(ix,iy,iz).var  *= inSamples;   // var: (phi - mean)^2
                    
                }
        
        
    }
};

#endif
