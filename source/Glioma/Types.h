/*
 *  Types.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/29/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


//static const int blockSize = 32;
//static const int blocksPerDimension = 4;//16;//8;//4;// 128/32;  //4
//static const int maxLevel = 3;
//
//
//// Structural parameters
//static const int	nDim			= 3;
//static const bool	bIsCellCentered = true;
//static const bool	bVerbose		= true;
//
//// Simulation parameters
//// (static so they can be change during UQ process and also for different IC set up, most of them are set to zero, so assign correct value w.r.t ic is required)
//static double L						= 81;								// length of the brain in mm    21.7; 81 mm = lenght of the brain data,
//static double D_g					= 0.;								// diffusion coefficient in grey matter, 1.3e-3 cm^2 / day, -> rescale for mm^2 scale, plus rescale by the length of the domain
//static double D_w					= 0.;								// diffusion coefficient in white matter
//static double rho					= 0.012;							// reaction rate, [rho] = 1 / day
//static const double Omega			= 1.0e7;							// concentration of cells per unif of volume
//
//// Multiresolution parameters
//static const int resJump			= 1;
//const double refinement_tolerance	= 1e-3;
//const double compression_tolerance	= 1e-4;
