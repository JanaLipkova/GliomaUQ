//
//  Glioma_BMD_ReactionDiffusion_Operator.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 07/02/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef GliomaBrutusXcode_Glioma_BMD_ReactionDiffusion_Operator_h
#define GliomaBrutusXcode_Glioma_BMD_ReactionDiffusion_Operator_h

template<int nDim = 3>
struct Glioma_BMD_ReactionDiffusion_Operator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dscale, rho;
    
    Glioma_BMD_ReactionDiffusion_Operator(const Real Dscale_, const Real rho_): Dscale(Dscale_), rho(rho_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    Glioma_BMD_ReactionDiffusion_Operator(const Glioma_BMD_ReactionDiffusion_Operator& copy): Dscale(copy.Dscale), rho(copy.rho)
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
        double df[7];
        double tau = 0.5;
        
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0; df[4] = 0.0; df[5] = 0.0; df[6] = 0.0;
               
                    // check if we are in the bone domain
                    if ( ( lab(ix, iy, iz).chi > tau) )
                    {
                        // lab().bmd = 1. / (bmd*0.01)^4)
                        df[0]  =  Dscale * lab(ix  ,iy  ,iz-1).bmd ;
                        df[1]  =  Dscale * lab(ix  ,iy  ,iz+1).bmd ;
                        df[2]  =  Dscale * lab(ix  ,iy-1,iz  ).bmd ;
                        df[3]  =  Dscale * lab(ix  ,iy+1,iz  ).bmd ;
                        df[4]  =  Dscale * lab(ix-1,iy  ,iz  ).bmd ;
                        df[5]  =  Dscale * lab(ix+1,iy  ,iz  ).bmd ;
                        df[6]  =  Dscale * lab(ix  ,iy  ,iz  ).bmd ;
                        
                        // Harmonic averages of piecewise constant diffusion coefficients
                        // Bernstein 2005 is wrong: factor of 2 is missing
                        // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                        // so we directly obtain diffusion zero for gp out of domain :)))
                        
                        df[0] = 2. / ( (1./df[6]) + (1./ df[0] )  ) ;
                        df[1] = 2. / ( (1./df[6]) + (1./ df[1] )  ) ;
                        df[2] = 2. / ( (1./df[6]) + (1./ df[2] )  ) ;
                        df[3] = 2. / ( (1./df[6]) + (1./ df[3] )  ) ;
                        df[4] = 2. / ( (1./df[6]) + (1./ df[4] )  ) ;
                        df[5] = 2. / ( (1./df[6]) + (1./ df[5] )  ) ;
                        
                        // Neumann no flux boundary condition, 2nd order using ghosts
                        // need correction by factor of 2
                        // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                        // if whole bone is eaten, then D will be zero in region with tumour -> but we shouldn't apply BC here
                        if ( (df[0] == 0)  ){ df[1] *= 2.0; }
                        if ( (df[1] == 0)  ){ df[0] *= 2.0; }
                        if ( (df[2] == 0)  ){ df[3] *= 2.0; }
                        if ( (df[3] == 0)  ){ df[2] *= 2.0; }
                        if ( (df[4] == 0)  ){ df[5] *= 2.0; }
                        if ( (df[5] == 0)  ){ df[4] *= 2.0; }
                    }
                    
                    
                    // diffusion fluxes
                    double diffusionFluxIn  = ih2 * ( df[0] * lab(ix, iy, iz-1).phi +
                                                     df[1] * lab(ix, iy, iz+1).phi +
                                                     df[2] * lab(ix, iy-1, iz).phi +
                                                     df[3] * lab(ix, iy+1, iz).phi +
                                                     df[4] * lab(ix-1, iy, iz).phi +
                                                     df[5] * lab(ix+1, iy, iz).phi   );
                    
                    double diffusionFluxOut = -( ( df[0] + df[1] + df[2] + df[3] + df[4] + df[5] ) * lab(ix, iy, iz).phi * ih2 );
                    double reactionFlux		= rho * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi );
                    
                    o(ix, iy, iz).dphidt =   diffusionFluxOut + diffusionFluxIn + reactionFlux ;
                    
                }
        
        
    }
};

#endif
