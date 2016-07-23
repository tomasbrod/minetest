/*
Minetest
Copyright (C) 2010-2014 kwolekr, Ryan Kwolek <kwolekr@minetest.net>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Instead (inaddition) of this stupid license text that is already in
license.txt or similar here should be quick description of what this
file contains.

Ore spawning algorithms and OreManager.

*/

#include "mg_ore.h"
#include "mapgen.h"
#include "noise.h"
#include "util/numeric.h"
#include "map.h"
#include "log.h"

FlagDesc flagdesc_ore[] = {
	{"absheight",                 OREFLAG_ABSHEIGHT},
	{"puff_cliffs",               OREFLAG_PUFF_CLIFFS},
	{"puff_additive_composition", OREFLAG_PUFF_ADDITIVE},
	{NULL,                        0}
};


///////////////////////////////////////////////////////////////////////////////


OreManager::OreManager(IGameDef *gamedef) :
	ObjDefManager(gamedef, OBJDEF_ORE)
{
}


size_t OreManager::placeAllOres(Mapgen *mg, u32 blockseed, v3s16 nmin, v3s16 nmax)
{
	size_t nplaced = 0;

	for (size_t i = 0; i != m_objects.size(); i++) {
		Ore *ore = (Ore *)m_objects[i];
		if (!ore)
			continue;

		nplaced += ore->placeOre(mg, blockseed, nmin, nmax);
		blockseed++;
	}

	return nplaced;
}


void OreManager::clear()
{
	for (size_t i = 0; i < m_objects.size(); i++) {
		Ore *ore = (Ore *)m_objects[i];
		delete ore;
	}
	m_objects.clear();
}


///////////////////////////////////////////////////////////////////////////////


Ore::Ore()
{
	flags = 0;
	noise = NULL;
}


Ore::~Ore()
{
	delete noise;
}


void Ore::resolveNodeNames()
{
	getIdFromNrBacklog(&c_ore, "", CONTENT_AIR);
	getIdsFromNrBacklog(&c_wherein);
}


size_t Ore::placeOre(Mapgen *mg, u32 blockseed, v3s16 nmin, v3s16 nmax)
{
	int in_range = 0;

	in_range |= (nmin.Y <= y_max && nmax.Y >= y_min);
	if (flags & OREFLAG_ABSHEIGHT)
		in_range |= (nmin.Y >= -y_max && nmax.Y <= -y_min) << 1;
	if (!in_range)
		return 0;

	int actual_ymin, actual_ymax;
	if (in_range & ORE_RANGE_MIRROR) {
		actual_ymin = MYMAX(nmin.Y, -y_max);
		actual_ymax = MYMIN(nmax.Y, -y_min);
	} else {
		actual_ymin = MYMAX(nmin.Y, y_min);
		actual_ymax = MYMIN(nmax.Y, y_max);
	}
	if (clust_size >= actual_ymax - actual_ymin + 1)
		return 0;

	nmin.Y = actual_ymin;
	nmax.Y = actual_ymax;
	generate(mg->vm, mg->seed, blockseed, nmin, nmax, mg->biomemap);

	return 1;
}


///////////////////////////////////////////////////////////////////////////////


void OreScatter::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
{
	PcgRandom pr(blockseed);
	MapNode n_ore(c_ore, 0, ore_param2);

	u32 sizex  = (nmax.X - nmin.X + 1);
	u32 volume = (nmax.X - nmin.X + 1) *
				 (nmax.Y - nmin.Y + 1) *
				 (nmax.Z - nmin.Z + 1);
	u32 csize     = clust_size;
	u32 cvolume    = csize * csize * csize;
	u32 nclusters = volume / clust_scarcity;

	for (u32 i = 0; i != nclusters; i++) {
		int x0 = pr.range(nmin.X, nmax.X - csize + 1);
		int y0 = pr.range(nmin.Y, nmax.Y - csize + 1);
		int z0 = pr.range(nmin.Z, nmax.Z - csize + 1);

		if ((flags & OREFLAG_USE_NOISE) &&
			(NoisePerlin3D(&np, x0, y0, z0, mapseed) < nthresh))
			continue;

		if (biomemap && !biomes.empty()) {
			u32 index = sizex * (z0 - nmin.Z) + (x0 - nmin.X);
			std::set<u8>::iterator it = biomes.find(biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		for (u32 z1 = 0; z1 != csize; z1++)
		for (u32 y1 = 0; y1 != csize; y1++)
		for (u32 x1 = 0; x1 != csize; x1++) {
			if (pr.range(1, cvolume) > clust_num_ores)
				continue;

			u32 i = vm->m_area.index(x0 + x1, y0 + y1, z0 + z1);
			if (!CONTAINS(c_wherein, vm->m_data[i].getContent()))
				continue;

			vm->m_data[i] = n_ore;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////


void OreSheet::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
{
	PcgRandom pr(blockseed + 4234);
	MapNode n_ore(c_ore, 0, ore_param2);

	u16 max_height = column_height_max;
	int y_start_min = nmin.Y + max_height;
	int y_start_max = nmax.Y - max_height;

	int y_start = y_start_min < y_start_max ?
		pr.range(y_start_min, y_start_max) :
		(y_start_min + y_start_max) / 2;

	if (!noise) {
		int sx = nmax.X - nmin.X + 1;
		int sz = nmax.Z - nmin.Z + 1;
		noise = new Noise(&np, 0, sx, sz);
	}
	noise->seed = mapseed + y_start;
	noise->perlinMap2D(nmin.X, nmin.Z);

	size_t index = 0;
	for (int z = nmin.Z; z <= nmax.Z; z++)
	for (int x = nmin.X; x <= nmax.X; x++, index++) {
		float noiseval = noise->result[index];
		if (noiseval < nthresh)
			continue;

		if (biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		u16 height = pr.range(column_height_min, column_height_max);
		int ymidpoint = y_start + noiseval;
		int y0 = MYMAX(nmin.Y, ymidpoint - height * (1 - column_midpoint_factor));
		int y1 = MYMIN(nmax.Y, y0 + height - 1);

		for (int y = y0; y <= y1; y++) {
			u32 i = vm->m_area.index(x, y, z);
			if (!vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, vm->m_data[i].getContent()))
				continue;

			vm->m_data[i] = n_ore;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

OrePuff::OrePuff() :
	Ore()
{
	noise_puff_top    = NULL;
	noise_puff_bottom = NULL;
}


OrePuff::~OrePuff()
{
	delete noise_puff_top;
	delete noise_puff_bottom;
}


void OrePuff::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
{
	PcgRandom pr(blockseed + 4234);
	MapNode n_ore(c_ore, 0, ore_param2);

	int y_start = pr.range(nmin.Y, nmax.Y);

	if (!noise) {
		int sx = nmax.X - nmin.X + 1;
		int sz = nmax.Z - nmin.Z + 1;
		noise = new Noise(&np, 0, sx, sz);
		noise_puff_top = new Noise(&np_puff_top, 0, sx, sz);
		noise_puff_bottom = new Noise(&np_puff_bottom, 0, sx, sz);
	}

	noise->seed = mapseed + y_start;
	noise->perlinMap2D(nmin.X, nmin.Z);
	bool noise_generated = false;

	size_t index = 0;
	for (int z = nmin.Z; z <= nmax.Z; z++)
	for (int x = nmin.X; x <= nmax.X; x++, index++) {
		float noiseval = noise->result[index];
		if (noiseval < nthresh)
			continue;

		if (biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		if (!noise_generated) {
			noise_generated = true;
			noise_puff_top->perlinMap2D(nmin.X, nmin.Z);
			noise_puff_bottom->perlinMap2D(nmin.X, nmin.Z);
		}

		float ntop    = noise_puff_top->result[index];
		float nbottom = noise_puff_bottom->result[index];

		if (!(flags & OREFLAG_PUFF_CLIFFS)) {
			float ndiff = noiseval - nthresh;
			if (ndiff < 1.0f) {
				ntop *= ndiff;
				nbottom *= ndiff;
			}
		}

		int ymid = y_start;
		int y0 = ymid - nbottom;
		int y1 = ymid + ntop;

		if ((flags & OREFLAG_PUFF_ADDITIVE) && (y0 > y1))
			SWAP(int, y0, y1);

		for (int y = y0; y <= y1; y++) {
			u32 i = vm->m_area.index(x, y, z);
			if (!vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, vm->m_data[i].getContent()))
				continue;

			vm->m_data[i] = n_ore;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////


void OreBlob::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
{
	PcgRandom pr(blockseed + 2404);
	MapNode n_ore(c_ore, 0, ore_param2);

	u32 sizex  = (nmax.X - nmin.X + 1);
	u32 volume = (nmax.X - nmin.X + 1) *
				 (nmax.Y - nmin.Y + 1) *
				 (nmax.Z - nmin.Z + 1);
	u32 csize  = clust_size;
	u32 nblobs = volume / clust_scarcity;

	if (!noise)
		noise = new Noise(&np, mapseed, csize, csize, csize);

	for (u32 i = 0; i != nblobs; i++) {
		int x0 = pr.range(nmin.X, nmax.X - csize + 1);
		int y0 = pr.range(nmin.Y, nmax.Y - csize + 1);
		int z0 = pr.range(nmin.Z, nmax.Z - csize + 1);

		if (biomemap && !biomes.empty()) {
			u32 bmapidx = sizex * (z0 - nmin.Z) + (x0 - nmin.X);
			std::set<u8>::iterator it = biomes.find(biomemap[bmapidx]);
			if (it == biomes.end())
				continue;
		}

		bool noise_generated = false;
		noise->seed = blockseed + i;

		size_t index = 0;
		for (u32 z1 = 0; z1 != csize; z1++)
		for (u32 y1 = 0; y1 != csize; y1++)
		for (u32 x1 = 0; x1 != csize; x1++, index++) {
			u32 i = vm->m_area.index(x0 + x1, y0 + y1, z0 + z1);
			if (!CONTAINS(c_wherein, vm->m_data[i].getContent()))
				continue;

			// Lazily generate noise only if there's a chance of ore being placed
			// This simple optimization makes calls 6x faster on average
			if (!noise_generated) {
				noise_generated = true;
				noise->perlinMap3D(x0, y0, z0);
			}

			float noiseval = noise->result[index];

			float xdist = (s32)x1 - (s32)csize / 2;
			float ydist = (s32)y1 - (s32)csize / 2;
			float zdist = (s32)z1 - (s32)csize / 2;

			noiseval -= (sqrt(xdist * xdist + ydist * ydist + zdist * zdist) / csize);

			if (noiseval < nthresh)
				continue;

			vm->m_data[i] = n_ore;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

OreVein::OreVein() :
	Ore()
{
	noise2 = NULL;
}


OreVein::~OreVein()
{
	delete noise2;
}


void OreVein::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
{
	PcgRandom pr(blockseed + 520);
	MapNode n_ore(c_ore, 0, ore_param2);

	u32 sizex = (nmax.X - nmin.X + 1);

	if (!noise) {
		int sx = nmax.X - nmin.X + 1;
		int sy = nmax.Y - nmin.Y + 1;
		int sz = nmax.Z - nmin.Z + 1;
		noise  = new Noise(&np, mapseed, sx, sy, sz);
		noise2 = new Noise(&np, mapseed + 436, sx, sy, sz);
	}
	bool noise_generated = false;

	size_t index = 0;
	for (int z = nmin.Z; z <= nmax.Z; z++)
	for (int y = nmin.Y; y <= nmax.Y; y++)
	for (int x = nmin.X; x <= nmax.X; x++, index++) {
		u32 i = vm->m_area.index(x, y, z);
		if (!vm->m_area.contains(i))
			continue;
		if (!CONTAINS(c_wherein, vm->m_data[i].getContent()))
			continue;

		if (biomemap && !biomes.empty()) {
			u32 bmapidx = sizex * (z - nmin.Z) + (x - nmin.X);
			std::set<u8>::iterator it = biomes.find(biomemap[bmapidx]);
			if (it == biomes.end())
				continue;
		}

		// Same lazy generation optimization as in OreBlob
		if (!noise_generated) {
			noise_generated = true;
			noise->perlinMap3D(nmin.X, nmin.Y, nmin.Z);
			noise2->perlinMap3D(nmin.X, nmin.Y, nmin.Z);
		}

		// randval ranges from -1..1
		float randval   = (float)pr.next() / (pr.RANDOM_RANGE / 2) - 1.f;
		float noiseval  = contour(noise->result[index]);
		float noiseval2 = contour(noise2->result[index]);
		if (noiseval * noiseval2 + randval * random_factor < nthresh)
			continue;

		vm->m_data[i] = n_ore;
	}
}

///////////////////////////////////////////////////////////////////////////////

void OrePipe::placePipe(MMVManip *vm, v3s16 nmin, v3s16 nmax, PcgRandom pr, v3f pA, v3f pB, v3f pC)
/*
	Spawn a single "pipe" between points A and B with control point C
  $hint is nmin and nmax value in in vm object?
*/
{
  const float stepadj=1.7;
	MapNode n_ore(c_ore, 0, ore_param2);
  v3f sa_p=pA;
  v3f p;
  //*estimate a step size*
  float step=(stepadj*pipe_radius)/(pA.getDistanceFrom(pB)+pB.getDistanceFrom(pC));
  for(float t=0; t<=1; t=t+step) //drawing loop
  {
    char sa_c=0;
    while(sa_c++<6) //step adjust loop
    {
      //*compute curve point*
      p  = pA * ((1-t)*(1-t));
      p += pC * (2*t*(1-t));
      p += pB * (t*t);

      //*adjust step*
      if(t==0) break;
      float sa_e = p.getDistanceFrom(sa_p)-pipe_radius;
      if (sa_e>=1)
        { step=step*0.50; t=t-step; }
      else if (sa_e<0)
        { t=t+(step/2); step=step+(step/2); }
      else break;
      printf("step_adjust %d %f %f \n",sa_c,sa_e,step);
    }
    p.X=round(p.X);p.Y=round(p.Y);p.Z=round(p.Z);
    sa_p=p;

    //*convert float[3] to area index*
    u32 di = vm->m_area.index(p.X, p.Y, p.Z); //maybe round...

    //*check can place*
    if (!vm->m_area.contains(di))
      continue; //stupid style
    if (!CONTAINS(c_wherein, vm->m_data[di].getContent()))
      continue;

    //*place the actual pipe segment*
    //$note todo use brush
    vm->m_data[di] = n_ore;
  }
}

void OrePipe::generate(MMVManip *vm, int mapseed, u32 blockseed,
	v3s16 nmin, v3s16 nmax, u8 *biomemap)
/*
  Generate all "pipe" ores in current mg-block
  not using noise
*/
{
	PcgRandom pr(blockseed);

	u32 sizex  = (nmax.X - nmin.X + 1);
	u32 volume = (nmax.X - nmin.X + 1) *
				 (nmax.Y - nmin.Y + 1) *
				 (nmax.Z - nmin.Z + 1);

	float gencnt_raw = (float)volume / (float)clust_scarcity;
  u32 gencnt = ((float)(gencnt_raw) + (float)((((float)pr.next())/((float)pr.RANDOM_RANGE))));
  printf("gencnt: %d..%d %d\n", nmin.X,nmax.X, gencnt);
	for (u32 i = 0; i != gencnt; i++) {
    v3f A, B, C;
    float dir_theta_ar = dir_theta + ((float)pr.next()/pr.RANDOM_RANGE)*dir_theta_rnd;
    float dir_phi_ar = ((float)pr.next()/pr.RANDOM_RANGE)*3.1416;
    v3f dir_nv ( sin(dir_phi_ar)*cos(dir_theta_ar), sin(dir_phi_ar)*sin(dir_theta_ar), cos(dir_phi_ar) );
    float length = pipe_length + (pr.next()*pipe_length_rnd/pr.RANDOM_RANGE);
    v3f dirv ( round(dir_nv.X*length), round(dir_nv.Y*length), round(dir_nv.Z*length) );
    
    printf("dir theta=%f phi=%f length=%f nv=(%f,%f,%f) v=(%f.1,%f.1,%f.1) d=%f\n",dir_theta_ar,dir_phi_ar,length,dir_nv.X,dir_nv.Y,dir_nv.Z,dirv.X,dirv.Y,dirv.Z,dirv.getLength());

		A.X = pr.range(nmin.X, nmax.X +1 -dirv.X);
		A.Y = pr.range(nmin.Y, nmax.Y +1 -dirv.Y);
		if(dirv.Z>0)
			{A.Z = pr.range(nmin.Z, nmax.Z +1 -dirv.Z);}
			else
			{A.Z = pr.range(nmin.Z -dirv.Z, nmax.Z +1);}
		
		B = A + dirv;

    //$note todo calculate A,B based on dit_theta,pipe_length,etc
    /*A.X = pr.range(nmin.X, nmax.X + 1); //not sure about that +1
    A.Y = pr.range(nmin.Y, nmax.Y + 1);
    A.Z = pr.range(nmin.Z, nmax.Z + 1);
    B.X = pr.range(nmin.X, nmax.X + 1);
    B.Y = pr.range(nmin.Y, nmax.Y + 1);
    B.Z = pr.range(nmin.Z, nmax.Z + 1);*/
    C.X = pr.range(nmin.X, nmax.X + 1);
    C.Y = pr.range(nmin.Y, nmax.Y + 1);
    C.Z = pr.range(nmin.Z, nmax.Z + 1);

    C = A;

		if (biomemap && !biomes.empty()) {
			u32 index = sizex * (A.Z - nmin.Z) + (A.X - nmin.X);
			std::set<u8>::iterator it = biomes.find(biomemap[index]);
			if (it == biomes.end())
				continue;
		}
    printf("place (%f.0,%f.0,%f.0),(%f.0,%f.0,%f.0)\n",A.X,A.Y,A.Z,B.X,B.Y,B.Z);
    placePipe(vm, nmin, nmax, pr, A, B, C);

	}
}

