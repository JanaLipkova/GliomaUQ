
/*
 *  GliomaTypes.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#pragma once
#include "Glioma.h"
#include "Matrix.h"
#include "Glioma_ReactionDiffusionOperator.h"
#include "Glioma_MALA_SlopesOperator.h"
#include "Glioma_ReactionDiffusionNecrosis_Operator.h"
#include "Glioma_BMD_ReactionDiffusion_Operator.h"
#include "Glioma_PropagationStatisticsOperator.h"
#include "Turing_ReactionDiffusionOperator.h"

struct CCell_extended
{
    
    /* tumor */
    Real phi;
    Real dphidt;
    
    /* tissue percentage per voxel*/
	Real p_g, p_w;
	Real p_skull;
	Real p_csf;
    
    /* tissue concentration */
    Real wm, gm, csf;  //
    Real dwmdt, dgmdt, dcsfdt;
    
    // pressure eq. channels for velocity, pressure, rhs, phase field f., domain charac. f, exact solution
	Real ux,uy,uz;
	Real p;
    Real f;
    Real psi, psiT, dpsidt;   // pahse field function of whole anatomy, of tissue
    Real chi;      // domain characteristic function
    Real exact;
    Real mu;       // chemical potential for Cahn Hilliard computaion
    Real kappa;    // smooth relaxation factor

    
    // other functions
	Real omega;
	Real domegadt;
	Real eps;
	Real tmp;
    
    Real t1bc, t2bc;
    
    // UQ propagation
    Real mean, var, pdf;
    Real bmd;
    
    // UQ MALA - slope/gradient variables sRho = S^rho = @ phi \ @rho
    Real sRho, sD, dsRhodt, dsDdt;
    
    // Necorosis and Hypoxia
    Real necro, dnecrodt;
    
    // field for vp
    Real vp;
    
    
    
    
	CCell_extended()
	{
		phi		 = 0.0;
		dphidt	 = 0.0;
		p_g		 = 0.0;
		p_w		 = 0.0;
		p_skull	 = 0.0;
		p_csf    = 0.0;
        wm = gm = csf = 0.0;
        dwmdt = dgmdt = dcsfdt = 0.0;
		ux = uy = uz  = 0.0;
		p        = 0.0;
		omega	 = 0.0;
		domegadt = 0.0;
		eps      = 0.0;
		tmp      = 0.0;
		f        = 0.0;
		psi		 = 0.0;
        psiT	 = 0.0;
		dpsidt   = 0.0;
		mu       = 0.0;
        kappa    = 0.0;
        chi      = 0.0;
        exact    = 0.0;
        t1bc     = 0.0;
        t2bc     = 0.0;
        mean     = 0.0;
        var      = 0.0;
        pdf      = 0.0;
        bmd      = 0.0;
        necro    = 0.0;
        dnecrodt = 0.0;
        vp       = 0.0;
        sRho     = 0.0;
        sD       = 0.0;
        dsRhodt = dsDdt = 0.0;
    }
	
    CCell_extended(Real phi_, Real dphidt_, Real p_g_, Real p_w_, Real p_skull_, Real p_csf_, Real wm_, Real gm_, Real csf_, Real dwmdt_, Real dgmdt_, Real dcsfdt_, Real ux_, Real uy_, Real uz_, Real p_, Real omega_, Real domegadt_, Real eps_, Real tmp_ , Real f_, Real psi_, Real psiT_, Real dpsidt_, Real mu_, Real kappa_, Real chi_, Real exact_, Real t1bc_, Real t2bc_, Real mean_, Real var_, Real pdf_, Real bmd_, Real necro_, Real dnecrodt_, Real vp_, Real sRho_, Real sD_, Real dsRhodt_, Real dsDdt_)
	{
		phi		 = phi_	    ;
		dphidt	 = dphidt_  ;
		p_g		 = p_g_	    ;
		p_w		 = p_w_	    ;
		p_skull	 = p_skull_ ;
		p_csf	 = p_csf_   ;
        wm       = wm_      ;
        gm       = gm_      ;
        csf      = csf_     ;
        dwmdt    = dwmdt_   ;
        dgmdt    = dgmdt_   ;
        dcsfdt   = dcsfdt_  ;
		ux		 = ux_	    ;
		uy		 = uy_	    ;
		uz		 = uz_	    ;
		p		 = p_       ;
		omega    = omega_   ;
		domegadt = domegadt_;
		eps      = eps_     ;
		tmp      = tmp_     ;
		f        = f_       ;
		psi		 = psi_	    ;
        psiT	 = psiT_    ;
		dpsidt   = dpsidt_  ;
		mu       = mu_      ;
        kappa    = kappa_   ;
        chi      = chi_     ;
        exact    = exact_   ;
        t1bc     = t1bc_    ;
        t2bc     = t2bc_    ;
        mean     = mean_    ;
        var      = var_     ;
        pdf      = pdf_     ;
        bmd      = bmd_     ;
        necro    = necro_   ;
        dnecrodt = dnecrodt_;
        vp       = vp_      ;
        sRho     = sRho_    ;
        sD       = sD_      ;
        dsRhodt  = dsRhodt_ ;
        dsDdt    = dsDdt_   ;

    }
	
	void operator += (CCell_extended t) 
	{
		phi		 += t.phi	  ;
		dphidt	 += t.dphidt  ;
		p_g		 += t.p_g	  ;
		p_w		 += t.p_w	  ;
		p_skull	 += t.p_skull ;
		p_csf    += t.p_csf   ;
        wm       += t.wm      ;
        gm       += t.gm      ;
        csf      += t.csf     ;
        dwmdt    += t.dwmdt   ;
        dgmdt    += t.dgmdt   ;
        dcsfdt   += t.dcsfdt  ;
		ux		 += t.ux	  ;
		uy		 += t.uy	  ;
		uz		 += t.uz	  ;
		p		 += t.p 	  ;
		omega    += t.omega	  ;
		domegadt += t.domegadt;
		eps      += t.eps     ;
		tmp      += t.tmp	  ;
		f        += t.f		  ;
		psi		 += t.psi     ;
        psiT	 += t.psiT    ;
		dpsidt   += t.dpsidt  ;
		mu       += t.mu      ;
        kappa    += t.kappa   ;
        chi      += t.chi     ;
        exact    += t.exact   ;
        t1bc     += t.t1bc    ;
        t2bc     += t.t2bc    ;
        mean     += t.mean    ;
        var      += t.var     ;
        pdf      += t.pdf     ;
        bmd      += t.bmd     ;
        necro    += t.necro   ;
        dnecrodt += t.dnecrodt;
        vp       += t.vp      ;
        sRho     += t.sRho    ;
        sD       += t.sD      ;
        dsRhodt  += t.dsRhodt ;
        dsDdt    += t.dsDdt   ;
	}
	
	
	
	operator Real() 
	{
		return (Real)phi;
	}
	
	void integrate(float dt)
	{
		phi		+= dt*dphidt;
		dphidt	= 0.0;
	}
	
	template<int i>
	Real evaluate_concentration(double dt)
	{
		return  phi+dt*dphidt;
	}
	
	
	Real giveMe(int i, Real h=0)
	{
		switch(i)
		{
  
#ifdef HGG
            case 0: return phi;
            case 1: return phi + 0.1 * p_g + 0.2 * p_w + 2. * 0.3 * p_csf;
            case 2: return phi + 0.1 * p_g + 0.2 * p_w;
            case 3: return 0.1 * p_g + 0.2 * p_w;
            case 4: return t1bc;
            case 5: return t2bc;
            case 6: return p_g;
            case 7: return p_w;
            case 8: return p_csf;
            case 9: return omega;
#endif
        
                
#ifdef HGG_UQ_MALA
            case 0: return phi;
            case 1: return phi + 0.1 * p_g + 0.2 * p_w + 2. * 0.3 * p_csf;
            case 2: return phi + 0.1 * p_g + 0.2 * p_w;
            case 3: return 0.1 * p_g + 0.2 * p_w;
            case 4: return sRho;
            case 5: return sD;
#endif
 
#ifdef Propagation
            case 0: return 0.1 * p_g + 0.2 * p_w + 2. * 0.3 * p_csf;
            case 1: return 0.1 * p_g + 0.2 * p_w;
            case 2: return mean;
            case 3: return var;
#endif
                
#ifdef Visualisation
            case 0: return p_w;     // T1Gd
            case 1: return p_g;     // FLAIR
            case 2: return p_csf;   // PET
            case 3: return t1bc;    // T1 segm.
            case 4: return t2bc;    // T2 segm.
            case 5: return eps;     // distance map, RT
            case 6: return kappa;   // rec T1w
            case 7: return mu;      // rec FLAIR
            case 8: return psi;     // segm. rec. T1w
            case 9: return chi;     // segm. rec. FLAIR

#endif
   
                
#ifdef VP
            case 0: return vp;   // vp
            case 1: return psi;  // PET
            case 2: return phi;   // tumour
            case 3: return p_w;  // p_w
            case 4: return p_g;  // p_g
#endif
                
                
#ifdef Bone
            case 0: return phi;    // tumour
            case 1: return tmp;    // real bmd
            case 2: return 0.001 * tmp + phi;
            case 3: return bmd;     // 1 / (bmd * 0.01)^4
            case 4: return 0.1 * p_g + 0.2 * p_w + phi;
            case 5: return p_w;
            case 6: return p_g;
                
#endif
                
                
#ifdef Deformation
            case 0: return p*chi;
            case 1: return phi;
            case 2: return wm;
            case 3: return gm;
            case 4: return csf;
            case 5: return (1. - phi) * p_w;  // to remove tissue fraction from tumour region, that acted as a tissue memory
            case 6: return (1. - phi) * p_g;
            case 7: return p_csf;
            case 8: return ux;
            case 9: return uy;
            case 10: return uz;
            case 11: return psi;
            case 12: return psiT;
#endif
 
#ifdef CH
            case 0: return psi;
            case 1: return chi;
            case 2: return p_w;
            case 3: return p_g;
            case 4: return p_csf;
#endif
                
#ifdef Necrosis
            case 0: return phi;
            case 1: return phi + 0.1 * p_g + 0.2 * p_w;
            case 2: return necro;
            case 3: return max(phi - necro, (Real)0.);
#endif
        
                
#ifdef RT_margin
            case 0: return phi;
            case 1: return chi;
            case 2: return phi*chi;
            case 3: return t1bc;
            case 4: return t2bc;
#endif
                
#ifdef Turing
            case 0: return psi;
            case 1: return phi;
#endif



			default: abort(); return 0;
		}
	}
				
	
};

inline CCell_extended operator*(const CCell_extended& p, Real v)
{
	CCell_extended c;
	c.phi		= p.phi		 *v;
	c.dphidt	= p.dphidt	 *v;
	c.p_g		= p.p_g		 *v;
	c.p_w		= p.p_w		 *v;
	c.p_skull   = p.p_skull  *v;
	c.p_csf     = p.p_csf    *v;
    c.wm        = p.wm       *v;
    c.gm        = p.gm       *v;
    c.csf       = p.csf      *v;
    c.dwmdt     = p.dwmdt    *v;
    c.dgmdt     = p.dgmdt    *v;
    c.dcsfdt    = p.dcsfdt   *v;
	c.ux        = p.ux       *v;
	c.uy        = p.uy       *v;
	c.uz        = p.uz       *v;
	c.p         = p.p        *v;
	c.omega     = p.omega    *v;
	c.domegadt	= p.domegadt *v;
	c.eps       = p.eps      *v;
	c.tmp       = p.tmp      *v;
	c.f         = p.f        *v;
	c.psi		= p.psi		 *v;
    c.psiT		= p.psiT	 *v;
	c.dpsidt    = p.dpsidt   *v;
	c.mu        = p.mu       *v;
    c.kappa     = p.kappa    *v;
    c.chi       = p.chi      *v;
    c.exact     = p.exact    *v;
    c.t1bc      = p.t1bc     *v;
    c.t2bc      = p.t2bc     *v;
    c.mean      = p.mean     *v;
    c.var       = p.var      *v;
    c.pdf       = p.pdf      *v;
    c.bmd       = p.bmd      *v;
    c.necro     = p.necro    *v;
    c.dnecrodt  = p.dnecrodt *v;
    c.vp        = p.vp       *v;
    c.sRho      = p.sRho     *v;
    c.sD        = p.sD       *v;
    c.dsRhodt   = p.dsRhodt  *v;
    c.dsDdt     = p.dsDdt    *v;

	return c;
}


#pragma mark projectors

template <typename T, int i> 
inline Real RD_projector_impl_vtk(const T&t)
{
	//	return (Real)(0.2*t.p_w + 0.1*t.p_g );	
	return (Real)(0.1 * t.p_g + 0.2 * t.p_w + t.p_csf + t.phi );
}

template <typename T, int i> 
inline Real RD_projector_impl_wav(const T&t)
{
	//return i==0 ? (Real)(t.phi) : (Real)(t.p_w);  // for refinment w.r.t 2 channels
	return (Real)(t.phi) ;
}

make_projector(RD_Projector_VTK,      RD_projector_impl_vtk)
make_projector(RD_Projector_Wavelets, RD_projector_impl_wav)


#ifndef _DIM
#define _DIM 3
#endif

#ifndef _BLOCKSIZE_
#define _BLOCKSIZE_ 16
#endif

#ifndef _BLOCKSIZE_Z_
#define _BLOCKSIZE_Z_ _BLOCKSIZE_
#endif

#ifndef _BPD_
#define _BPD_ 4
#endif

#ifndef _MAXLEVEL_
#define _MAXLEVEL_ 3
#endif

static const int blockSize = _BLOCKSIZE_;
static const int blockSizeZ = _BLOCKSIZE_Z_;
static const int blocksPerDimension = _BPD_;


// Structural parameters
static const bool	bIsCellCentered = true;
static const bool	bVerbose		= true;

// Simulation parameters
// (static so they can be change during UQ process and also for different IC set up
static double L		;		// length of the brain in mm    21.7; 81 mm = lenght of the brain data,

static Real icp[3]	;	    // location of ICP point
static double ICP	;		// ICP in mmHg

static const double Omega = 1.0e7;		// concentration of cells per unif of volume
static double rho   ;


// Multiresolution parameters
//static const int maxLevel = 4;              // 16bpd has maxLevel 4 since 2^4
//static const int resJump			= 2;    // modulo(maxLevel,resJum) = 0, !!! and reJump < maxLevel, the bigger jump, the better adaptivity-> the faster

static const int maxLevel = _MAXLEVEL_;    // 4 for 16bpd, 3 for 8bpd has maxLevel 3 since 2^3
static const int resJump  = 1;    // modulo(maxLevel,resJum) = 0, !!! and reJump < maxLevel, the bigger jump, the better adaptivity-> the faster
const double refinement_tolerance	= 1e-4;
const double compression_tolerance	= 1e-5;

//typedef		Block< CCell_extended, blockSize, blockSize, blockSize>					B;
typedef		Block< CCell_extended, blockSize, blockSize, blockSizeZ>					B;
typedef		_WAVELET_TYPE															W;

typedef Matrix::D3D<double> MatrixD3D;
typedef Matrix::D2D<double> MatrixD2D;


#ifdef _MRAG_TBB
static const int nThreads = _MRAG_TBB_NTHREADS_HINT ;
typedef		Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, nThreads+1>	BlockProcessing;
#else
static const int nThreads = 1 ;
typedef		Multithreading::BlockProcessing_SingleCPU<B>							BlockProcessing;
#endif
