//
//  Glioma_SegmentationBC_Operator.h
//  GliomaXcode
//
//  Created by Lipkova on 21/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef GliomaXcode_Glioma_SegmentationBC_Operator_h
#define GliomaXcode_Glioma_SegmentationBC_Operator_h


template<int nDim = 3>
struct Glioma_SegmentationBC_Operator
{
    int stencil_start[3];
    int stencil_end[3];
    
    Glioma_SegmentationBC_Operator()
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    Glioma_SegmentationBC_Operator(const Glioma_SegmentationBC_Operator& copy)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        Real phiB, phiF, phiS, phiN, phiW, phiE;
        
        if(nDim == 2)
        {

            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    
                    // binnary T1 segm. in ux channel
                    phiS = lab(ix,iy-1).ux;
                    phiN = lab(ix,iy+1).ux;
                    phiW = lab(ix-1,iy).ux;
                    phiE = lab(ix+1,iy).ux;
                    
                    if ( ( (phiS+phiN) == 1.) || ( (phiW + phiE) == 1.  ) )
                        o(ix,iy).t1bc = 1.;
                    
                    // binnary T2 segm. in ux channel
                    phiS = lab(ix,iy-1).uy;
                    phiN = lab(ix,iy+1).uy;
                    phiW = lab(ix-1,iy).uy;
                    phiE = lab(ix+1,iy).uy;
                    
                    if ( ( (phiS+phiN) == 1.) || ( (phiW + phiE) == 1.  ) )
                        o(ix,iy).t2bc = 1.;
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        // binnary T1 segm. in ux channel
                        phiB = lab(ix,  iy,  iz-1).ux;
                        phiF = lab(ix,  iy,  iz+1).ux;
                        phiS = lab(ix,  iy-1,iz  ).ux;
                        phiN = lab(ix,  iy+1,iz  ).ux;
                        phiW = lab(ix-1,iy,  iz  ).ux;
                        phiE = lab(ix+1,iy,  iz  ).ux;
                        
                        if ( ( (phiS+phiN) == 1.) || ( (phiW + phiE) == 1.  ) || ( (phiF + phiB) == 1. ) )
                           o(ix,iy,iz).t1bc = 1.;
                        
                        
                        // binnary T2 segm. in uy channel
                        phiB = lab(ix,  iy,  iz-1).uy;
                        phiF = lab(ix,  iy,  iz+1).uy;
                        phiS = lab(ix,  iy-1,iz  ).uy;
                        phiN = lab(ix,  iy+1,iz  ).uy;
                        phiW = lab(ix-1,iy,  iz  ).uy;
                        phiE = lab(ix+1,iy,  iz  ).uy;
                        
                        if ( ( (phiS+phiN) == 1.) || ( (phiW + phiE) == 1.  ) || ( (phiF + phiB) == 1. ) )
                            o(ix,iy,iz).t2bc = 1.;
                        
                    }
        }
    }
};


#endif
