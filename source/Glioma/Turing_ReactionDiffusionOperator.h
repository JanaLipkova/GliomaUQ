//
//  Turing_ReactionDiffusionOperator.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 07/11/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

#ifndef _Turing_ReactionDiffusionOperator_h
#define _Turing_ReactionDiffusionOperator_h


/* Schnackenberg system:
    2U + V -> 3U   with k1
         0 -> U    with k2
         U -> 0    with k3
         0 -> V    with k4
 
 @u/@t = Du Δu + k1 * u^2 * v + k2 - k3 * u
 @v/@t = Dv Δv - k1 * u^2 * v + k4

NOTE: u = psi, v = phi
*/
struct SchnackenbergReactionDiffusionOperator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dpsi, Dphi, k1, k2, k3, k4;
    
    SchnackenbergReactionDiffusionOperator(const Real Dpsi_, const Real Dphi_, const Real k1_, const Real k2_, const Real k3_, const Real k4_): Dpsi(Dpsi_), Dphi(Dphi_), k1(k1_), k2(k2_), k3(k3_), k4(k4_)
    {
        stencil_start[0] = stencil_start[1]= stencil_start[2] = -1;
        stencil_end[0]   = stencil_end[1] =  stencil_end[2]   = +2;
    }
    
    SchnackenbergReactionDiffusionOperator(const SchnackenbergReactionDiffusionOperator& copy): Dpsi(copy.Dpsi), Dphi(copy.Dphi), k1(copy.k1), k2(copy.k2), k3(copy.k3), k4(copy.k4)
    {
        stencil_start[0] = stencil_start[1]= stencil_start[2] = -1;
        stencil_end[0]   = stencil_end[1]  = stencil_end[2]   = +2;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		= info.h[0];
        double ih2		= 1./(h*h);
        
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    
                    Real laplacePsi = ih2 * (lab(ix  ,iy  ,iz-1).psi +
                                             lab(ix  ,iy  ,iz+1).psi +
                                             lab(ix  ,iy-1,iz  ).psi +
                                             lab(ix  ,iy+1,iz  ).psi +
                                             lab(ix-1,iy  ,iz  ).psi +
                                             lab(ix+1,iy  ,iz  ).psi -
                                             6.* lab(ix  ,iy  ,iz  ).psi );
                    

                    Real laplacePhi = ih2 * (lab(ix  ,iy  ,iz-1).phi +
                                             lab(ix  ,iy  ,iz+1).phi +
                                             lab(ix  ,iy-1,iz  ).phi +
                                             lab(ix  ,iy+1,iz  ).phi +
                                             lab(ix-1,iy  ,iz  ).phi +
                                             lab(ix+1,iy  ,iz  ).phi -
                                             6.* lab(ix  ,iy  ,iz  ).phi );
                    
                    
                    o(ix,iy,iz).dpsidt = Dpsi * laplacePsi + k1 * lab(ix,iy,iz).psi * lab(ix,iy,iz).psi * lab(ix,iy,iz).phi + k2 - k3 * lab(ix,iy,iz).psi;
                    o(ix,iy,iz).dphidt = Dphi * laplacePhi - k1 * lab(ix,iy,iz).psi * lab(ix,iy,iz).psi * lab(ix,iy,iz).phi + k4;
                    
                }
    }
    
};



/* Gray-Scott system:
 U + 2V -> 3V
 V -> 0
 
 @u/@t = Du Δu - u^2 * v + F(1-u)
 @v/@t = Dv Δv + u^2 * v - (F+k)v
 
 NOTE: u = psi, v = phi
 */
struct GrayScottReactionDiffusionOperator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dpsi, Dphi, F, k;
    
    GrayScottReactionDiffusionOperator(const Real Dpsi_, const Real Dphi_, const Real F_, const Real k_): Dpsi(Dpsi_), Dphi(Dphi_), F(F_), k(k_)
    {
        stencil_start[0] = stencil_start[1]= stencil_start[2] = -1;
        stencil_end[0]   = stencil_end[1] =  stencil_end[2]   = +2;
    }
    
    GrayScottReactionDiffusionOperator(const GrayScottReactionDiffusionOperator& copy): Dpsi(copy.Dpsi), Dphi(copy.Dphi), F(copy.F), k(copy.k)
    {
        stencil_start[0] = stencil_start[1]= stencil_start[2] = -1;
        stencil_end[0]   = stencil_end[1]  = stencil_end[2]   = +2;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		= info.h[0];
        double ih2		= 1./(h*h);
        
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    
                    Real laplacePsi = ih2 * (lab(ix  ,iy  ,iz-1).psi +
                                             lab(ix  ,iy  ,iz+1).psi +
                                             lab(ix  ,iy-1,iz  ).psi +
                                             lab(ix  ,iy+1,iz  ).psi +
                                             lab(ix-1,iy  ,iz  ).psi +
                                             lab(ix+1,iy  ,iz  ).psi -
                                             6.* lab(ix  ,iy  ,iz  ).psi );
                    
                    
                    Real laplacePhi = ih2 * (lab(ix  ,iy  ,iz-1).phi +
                                             lab(ix  ,iy  ,iz+1).phi +
                                             lab(ix  ,iy-1,iz  ).phi +
                                             lab(ix  ,iy+1,iz  ).phi +
                                             lab(ix-1,iy  ,iz  ).phi +
                                             lab(ix+1,iy  ,iz  ).phi -
                                             6.* lab(ix  ,iy  ,iz  ).phi );
                    
                    
                    o(ix,iy,iz).dpsidt = Dpsi * laplacePsi - lab(ix,iy,iz).psi * lab(ix,iy,iz).psi * lab(ix,iy,iz).phi + F * (1. - lab(ix,iy,iz).psi );
                    o(ix,iy,iz).dphidt = Dphi * laplacePhi + lab(ix,iy,iz).psi * lab(ix,iy,iz).psi * lab(ix,iy,iz).phi - (F + k) * lab(ix,iy,iz).phi ;
                    
                }
    }
    
};




struct TuringTimeUpdate
{
    double dt;
    
    TuringTimeUpdate(double dt_):dt(dt_)
    { }
    
    TuringTimeUpdate(const TuringTimeUpdate& copy):dt(copy.dt)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix,iy,iz).psi += dt * o(ix,iy,iz).dpsidt;
                    o(ix,iy,iz).phi += dt * o(ix,iy,iz).dphidt;
                }
        
    }
};


#endif
