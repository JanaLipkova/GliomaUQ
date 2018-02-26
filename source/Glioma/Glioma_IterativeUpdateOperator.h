/*
 *  Glioma_IterativeUpdateOperator.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 10/5/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

struct ItertativeUpdate
{

	template<typename BlockType>
	inline void operator()(const BlockInfo& info, BlockType& o) const
	{
		for(int iz=0; iz<BlockType::sizeZ; iz++)
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					o(ix,iy,iz).dphidt = abs(o(ix,iy,iz).psi - o(ix,iy,iz).tmp);

					o(ix, iy, iz).psi = o(ix,iy,iz).tmp;
					
				}
	}
};