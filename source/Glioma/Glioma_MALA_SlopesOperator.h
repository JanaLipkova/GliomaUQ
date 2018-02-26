//
//  Glioma_MALA_SlopesOperator.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 29/07/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef _Glioma_MALA_SlopesOperator_h
#define _Glioma_MALA_SlopesOperator_h


template<int nDim = 3>
struct Glioma_MALAslopesRHS
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dw, Dg;
    
    Glioma_MALAslopesRHS(const Real Dw_, const Real Dg_): Dw(Dw_), Dg(Dg_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    Glioma_MALAslopesRHS(const Glioma_MALAslopesRHS& copy): Dw(copy.Dw), Dg(copy.Dg)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		= info.h[0];
        double ih2		= 1./(h*h);
        double ratioD   = Dg/Dw;
        double tt[6];   // Tissue tensor, Di = Dw * Ti = Dw * (p_w + ratioD * p_g)
        
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    tt[0] = 0.0; tt[1] = 0.0; tt[2] = 0.0; tt[3] = 0.0;
                    
                    // check if we are in the brain domain
                    if ( (lab(ix, iy).p_w > 0.0) || (lab(ix, iy).p_g > 0.0) )
                    {
                        // Harmonic averages of piecewise constant diffusion coefficients
                        // Bernstein 2005 is wrong: factor of 2 is missing
                        // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                        // so we directly obtain diffusion zero for gp out of domain :)))
                        
                        tt[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w + lab(ix, iy).p_g * ratioD) ) + (1.0 / (lab(ix-1, iy).p_w + lab(ix-1, iy).p_g * ratioD) ) ) );
                        tt[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w + lab(ix, iy).p_g * ratioD) ) + (1.0 / (lab(ix+1, iy).p_w + lab(ix+1, iy).p_g * ratioD) ) ) );
                        tt[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w + lab(ix, iy).p_g * ratioD) ) + (1.0 / (lab(ix, iy-1).p_w + lab(ix, iy-1).p_g * ratioD) ) ) );
                        tt[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w + lab(ix, iy).p_g * ratioD) ) + (1.0 / (lab(ix, iy+1).p_w + lab(ix, iy+1).p_g * ratioD) ) ) );
                    
                    // Neumann no flux boundary condition, 2nd order using ghosts
                    // need correction by factor of 2
                    // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                    if ( tt[0] == 0 ){ tt[1] *= 2.0; }
                    if ( tt[1] == 0 ){ tt[0] *= 2.0; }
                    if ( tt[2] == 0 ){ tt[3] *= 2.0; }
                    if ( tt[3] == 0 ){ tt[2] *= 2.0; }
                }
                    
                    
                    
                    /* reaction-diffusion fluxes for slopes */
                    
                    // 1) sRho: ∇ ( D ∇ sRho) + phi * (1 - phi)
                    double diffusionFluxIn_sRho  = ih2 * Dw * (tt[0] * lab(ix-1, iy).sRho +
                                                               tt[1] * lab(ix+1, iy).sRho +
                                                               tt[2] * lab(ix, iy-1).sRho +
                                                               tt[3] * lab(ix, iy+1).sRho );
                    
                    double diffusionFluxOut_sRho = - ih2 * Dw * (tt[0] + tt[1] + tt[2] + tt[3]) * lab(ix, iy).sRho;
                    double reactionFlux_sRho     = lab(ix,iy).phi * ( 1. - lab(ix,iy).phi );
                    
                    o(ix, iy).dsRhodt =   diffusionFluxOut_sRho + diffusionFluxIn_sRho + reactionFlux_sRho ;
                    
                    
                    // 2) sD:  ∇ ( D ∇ sD) +  ∇ ( T ∇ phi)
                    double diffusionFluxIn_sD    = ih2 * Dw * (tt[0] * lab(ix-1, iy).sD +
                                                               tt[1] * lab(ix+1, iy).sD +
                                                               tt[2] * lab(ix, iy-1).sD +
                                                               tt[3] * lab(ix, iy+1).sD );
                    
                    double diffusionFluxOut_sD   = - ih2 * Dw * (tt[0] + tt[1] + tt[2] + tt[3] ) * lab(ix, iy).sD;
                    
                    double reactionFluxIn_sD     = ih2  * (tt[0] * lab(ix-1, iy).phi +
                                                           tt[1] * lab(ix+1, iy).phi +
                                                           tt[2] * lab(ix, iy-1).phi +
                                                           tt[3] * lab(ix, iy+1).phi );
                    
                    double reactionFluxOut_sD   = - ih2 *  (tt[0] + tt[1] + tt[2] + tt[3] ) * lab(ix, iy).phi;
                    
                    o(ix, iy).dsDdt =   diffusionFluxOut_sD + diffusionFluxIn_sD + reactionFluxIn_sD + reactionFluxOut_sD ;
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        tt[0] = 0.0; tt[1] = 0.0; tt[2] = 0.0; tt[3] = 0.0; tt[4] = 0.0; tt[5] = 0.0;
                        
                        // check if we are in the brain domain
                        if ( (lab(ix, iy, iz).p_w > 0.0) || (lab(ix, iy, iz).p_g > 0.0) )
                        {
                            // Harmonic averages of piecewise constant diffusion coefficients
                            // Bernstein 2005 is wrong: factor of 2 is missing
                            // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                            // so we directly obtain diffusion zero for gp out of domain :)))
                            
                            tt[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix-1, iy, iz).p_w + lab(ix-1, iy, iz).p_g*ratioD) ) ) );
                            tt[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix+1, iy, iz).p_w + lab(ix+1, iy, iz).p_g*ratioD) ) ) );
                            tt[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix, iy-1, iz).p_w + lab(ix, iy-1, iz).p_g*ratioD) ) ) );
                            tt[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix, iy+1, iz).p_w + lab(ix, iy+1, iz).p_g*ratioD) ) ) );
                            tt[4] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix, iy, iz-1).p_w + lab(ix, iy, iz-1).p_g*ratioD) ) ) );
                            tt[5] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w + lab(ix, iy, iz).p_g*ratioD) ) + (1.0 / (lab(ix, iy, iz+1).p_w + lab(ix, iy, iz+1).p_g*ratioD) ) ) );
                            
                            
                            // Neumann no flux boundary condition, 2nd order using ghosts
                            // need correction by factor of 2
                            // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                            if ( tt[0] == 0 ){ tt[1] *= 2.0; }
                            if ( tt[1] == 0 ){ tt[0] *= 2.0; }
                            if ( tt[2] == 0 ){ tt[3] *= 2.0; }
                            if ( tt[3] == 0 ){ tt[2] *= 2.0; }
                            if ( tt[4] == 0 ){ tt[5] *= 2.0; }
                            if ( tt[5] == 0 ){ tt[4] *= 2.0; }
                        }
            
                        // reaction-diffusion fluxes for slopes
                        
                        // 1) sRho: ∇ ( D ∇ sRho) + phi * (1 - phi)
                        double diffusionFluxIn_sRho  = ih2 * Dw * (tt[0] * lab(ix-1, iy, iz).sRho +
                                                                   tt[1] * lab(ix+1, iy, iz).sRho +
                                                                   tt[2] * lab(ix, iy-1, iz).sRho +
                                                                   tt[3] * lab(ix, iy+1, iz).sRho +
                                                                   tt[4] * lab(ix, iy, iz-1).sRho +
                                                                   tt[5] * lab(ix, iy, iz+1).sRho   );
                        
                        double diffusionFluxOut_sRho = - ih2 * Dw * (tt[0] + tt[1] + tt[2] + tt[3] + tt[4] + tt[5]) * lab(ix, iy, iz).sRho;
                        double reactionFlux_sRho     = lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi );
                        
                        o(ix, iy, iz).dsRhodt =   diffusionFluxOut_sRho + diffusionFluxIn_sRho + reactionFlux_sRho ;


                        
                        // 2) sD:  ∇ ( D ∇ sD) +  ∇ ( T ∇ phi)
                        double diffusionFluxIn_sD    = ih2 * Dw * (tt[0] * lab(ix-1, iy, iz).sD +
                                                                   tt[1] * lab(ix+1, iy, iz).sD +
                                                                   tt[2] * lab(ix, iy-1, iz).sD +
                                                                   tt[3] * lab(ix, iy+1, iz).sD +
                                                                   tt[4] * lab(ix, iy, iz-1).sD +
                                                                   tt[5] * lab(ix, iy, iz+1).sD   );
                        
                        double diffusionFluxOut_sD   = - ih2 * Dw * (tt[0] + tt[1] + tt[2] + tt[3] + tt[4] + tt[5]) * lab(ix, iy, iz).sD;
                        
                        double reactionFluxIn_sD     = ih2  * (tt[0] * lab(ix-1, iy, iz).phi +
                                                               tt[1] * lab(ix+1, iy, iz).phi +
                                                               tt[2] * lab(ix, iy-1, iz).phi +
                                                               tt[3] * lab(ix, iy+1, iz).phi +
                                                               tt[4] * lab(ix, iy, iz-1).phi +
                                                               tt[5] * lab(ix, iy, iz+1).phi   );
                        
                        double reactionFluxOut_sD   = - ih2 *  (tt[0] + tt[1] + tt[2] + tt[3] + tt[4] + tt[5]) * lab(ix, iy, iz).phi;
                        
                        o(ix, iy, iz).dsDdt =   diffusionFluxOut_sD + diffusionFluxIn_sD + reactionFluxIn_sD + reactionFluxOut_sD ;


                    }
        }
    }
};




template<int nDim = 3>
struct UpdateMALAslopes
{
    double dt;
    
    UpdateMALAslopes(double dt_):dt(dt_)
    { }
    
    UpdateMALAslopes(const UpdateMALAslopes& copy):dt(copy.dt)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix, iy).sRho  += dt * o(ix, iy).dsRhodt;
                    o(ix, iy).sD    += dt * o(ix, iy).dsDdt;
                }
            
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        o(ix, iy, iz).sRho  += dt * o(ix, iy, iz).dsRhodt;
                        o(ix, iy, iz).sD    += dt * o(ix, iy, iz).dsDdt;
                    }
        }
        
    }
};



#endif


