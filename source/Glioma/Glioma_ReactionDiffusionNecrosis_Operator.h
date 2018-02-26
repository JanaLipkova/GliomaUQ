//
//  Glioma_ReactionDiffusionNecrosis_Operator.h
//  GliomaBrutusXcode
//
//  Created by Lipkova on 15/03/16.
//  Copyright (c) 2016 Lipkova. All rights reserved.
//

#ifndef _Glioma_ReactionDiffusionNecrosis_Operator_h
#define _Glioma_ReactionDiffusionNecrosis_Operator_h



template<int nDim = 3>
struct Glioma_ReactionDiffusionNecrosis_Operator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dw, Dg, rho, gamma;
    
    Glioma_ReactionDiffusionNecrosis_Operator(const Real Dw_, const Real Dg_, const Real rho_, const Real gamma_): Dw(Dw_), Dg(Dg_), rho(rho_),gamma(gamma_)    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    Glioma_ReactionDiffusionNecrosis_Operator(const Glioma_ReactionDiffusionNecrosis_Operator& copy): Dw(copy.Dw), Dg(copy.Dg), rho(copy.rho), gamma(copy.gamma)
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
        double df[6];
        double hypoxic_thr = 0.4;
        double latent_thr = 0.7;
        double gamma_rate;
        double Hypoxia;
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0;
                  
                     gamma_rate = 0.;
                     Hypoxia = 0;

                    // check if we are in the brain domain
                    if ( (lab(ix, iy).p_w > 0.0) || (lab(ix, iy).p_g > 0.0) )
                    {
                        
                        // Harmonic averages of piecewise constant diffusion coefficients
                        // Bernstein 2005 is wrong: factor of 2 is missing
                        // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                        // so we directly obtain diffusion zero for gp out of domain :)))
                        df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix-1, iy).p_w*Dw + lab(ix-1, iy).p_g*Dg) ) ) );
                        df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix+1, iy).p_w*Dw + lab(ix+1, iy).p_g*Dg) ) ) );
                        df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy-1).p_w*Dw + lab(ix, iy-1).p_g*Dg) ) ) );
                        df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy+1).p_w*Dw + lab(ix, iy+1).p_g*Dg) ) ) );
                        
                        
                        // Neumann no flux boundary condition, 2nd order using ghosts
                        // need correction by factor of 2
                        // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                        if ( (df[0] == 0) && ( lab(ix-1,iy).phi == 0) ){ df[1] *= 2.0; }
                        if ( (df[1] == 0) && ( lab(ix+1,iy).phi == 0) ){ df[0] *= 2.0; }
                        if ( (df[2] == 0) && ( lab(ix,iy-1).phi == 0) ){ df[3] *= 2.0; }
                        if ( (df[3] == 0) && ( lab(ix,iy+1).phi == 0) ){ df[2] *= 2.0; }
                        
                        // Hypoxia
                        if( (hypoxic_thr < lab(ix,iy).phi) && (lab(ix,iy).phi < latent_thr) )
                            Hypoxia = 1. - (latent_thr - lab(ix,iy).phi ) / ( latent_thr - hypoxic_thr);
                        else if ( lab(ix,iy).phi >= latent_thr )
                                Hypoxia = 1;
                        
                        double tmp = 30. * (latent_thr - lab(ix,iy).phi);
                        double tmp2rad = tmp * M_PI/180.;
                        
                        if(lab(ix,iy).phi >= latent_thr )
                            gamma_rate = gamma * 0.5 * (1. - tanh(tmp2rad));
                        
                        // diffusion fluxes
                        double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy).phi +
                                                         df[1]*lab(ix+1, iy).phi +
                                                         df[2]*lab(ix, iy-1).phi +
                                                         df[3]*lab(ix, iy+1).phi  );
                        
                        double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3]) * lab(ix, iy).phi * ih2 );
                        double reactionFlux		= rho * (1. - Hypoxia) * lab(ix,iy).phi * ( 1. - lab(ix,iy).phi );
                        double necroticFlux     = gamma_rate * lab(ix,iy).phi;
                        
                        o(ix, iy).dphidt   =   diffusionFluxOut + diffusionFluxIn + reactionFlux - necroticFlux;
                        o(ix, iy).dnecrodt =   necroticFlux;
                    }
                    
                    
                }
        }
        else
        {
            
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0; df[4] = 0.0; df[5] = 0.0;
                        gamma_rate = 0.;
                        Hypoxia = 0;
                        
                        // check if we are in the brain domain
                        if ( (lab(ix, iy, iz).p_w > 0.0) || (lab(ix, iy, iz).p_g > 0.0) )
                        {
                            // Harmonic averages of piecewise constant diffusion coefficients
                            // Bernstein 2005 is wrong: factor of 2 is missing
                            // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                            // so we directly obtain diffusion zero for gp out of domain :)))
                            df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix-1, iy, iz).p_w*Dw + lab(ix-1, iy, iz).p_g*Dg) ) ) );
                            df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix+1, iy, iz).p_w*Dw + lab(ix+1, iy, iz).p_g*Dg) ) ) );
                            df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy-1, iz).p_w*Dw + lab(ix, iy-1, iz).p_g*Dg) ) ) );
                            df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy+1, iz).p_w*Dw + lab(ix, iy+1, iz).p_g*Dg) ) ) );
                            df[4] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz-1).p_w*Dw + lab(ix, iy, iz-1).p_g*Dg) ) ) );
                            df[5] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz+1).p_w*Dw + lab(ix, iy, iz+1).p_g*Dg) ) ) );
                            
                            
                            // Neumann no flux boundary condition, 2nd order using ghosts
                            // need correction by factor of 2
                            // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                            if ( df[0] == 0 ){ df[1] *= 2.0; }
                            if ( df[1] == 0 ){ df[0] *= 2.0; }
                            if ( df[2] == 0 ){ df[3] *= 2.0; }
                            if ( df[3] == 0 ){ df[2] *= 2.0; }
                            if ( df[4] == 0 ){ df[5] *= 2.0; }
                            if ( df[5] == 0 ){ df[4] *= 2.0; }
                            
                            
                            // Hypoxia
//                            if( (hypoxic_thr < lab(ix,iy,iz).phi) && (lab(ix,iy,iz).phi < latent_thr) )
//                                Hypoxia = 1. - (latent_thr - lab(ix,iy,iz).phi ) / ( latent_thr - hypoxic_thr);
//                                else if ( lab(ix,iy,iz).phi >= latent_thr )
//                                    Hypoxia = 1;
//                            
                            
                            if ( lab(ix,iy,iz).phi > latent_thr)
                            {
                                double tmp = 30. * (latent_thr - lab(ix,iy,iz).phi);
                                double tmp2rad = tmp * M_PI/180.;
                                gamma_rate = gamma * 0.5 * (1. - tanh(tmp2rad));
                            }
                            
                            // diffusion fluxes
                            double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy, iz).phi +
                                                             df[1]*lab(ix+1, iy, iz).phi +
                                                             df[2]*lab(ix, iy-1, iz).phi +
                                                             df[3]*lab(ix, iy+1, iz).phi +
                                                             df[4]*lab(ix, iy, iz-1).phi +
                                                             df[5]*lab(ix, iy, iz+1).phi   );
                            
                            double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3] + df[4] + df[5]) * lab(ix, iy, iz).phi * ih2 );
                            double reactionFlux		= rho * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi - lab(ix,iy,iz).necro );
                            double necroticFlux     = gamma_rate * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi - lab(ix,iy,iz).necro  ) * ( 1. + 2. *  lab(ix,iy,iz).necro) ;

                            o(ix, iy, iz).dphidt   =   diffusionFluxOut + diffusionFluxIn + reactionFlux - necroticFlux;
                            o(ix, iy, iz).dnecrodt =   necroticFlux;

                        }
                    }
            
            
        }
    }
};




template<int nDim = 3>
struct TimeUpdate
{
    double dt;
    
    TimeUpdate(double dt_):dt(dt_)
    { }
    
    TimeUpdate(const TimeUpdate& copy):dt(copy.dt)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix, iy).phi   += dt * o(ix, iy).dphidt;
                    o(ix, iy).necro += dt * o(ix, iy).dnecrodt;
                    o(ix, iy).phi    = max((Real)0., o(ix,iy).phi);
                    o(ix, iy).phi    = min((Real)1., o(ix, iy).phi);
                }
            
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        o(ix, iy, iz).phi   += dt * o(ix, iy, iz).dphidt;
                        o(ix, iy, iz).necro += dt * o(ix, iy, iz).dnecrodt;
                        
                        o(ix, iy, iz).phi    = max((Real)0., o(ix,iy,iz).phi);
                        o(ix, iy, iz).phi    = min((Real)1., o(ix,iy,iz).phi);
                        o(ix, iy, iz).necro  = max((Real)0., o(ix,iy,iz).necro);
                        o(ix, iy, iz).necro  = min((Real)1., o(ix,iy,iz).necro);
                    }
            
        }
        
    }
};



#endif
