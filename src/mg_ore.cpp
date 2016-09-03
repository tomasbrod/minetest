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

OreManager::placeAllOres for each registered ore calls Ore::placeOre,
which does some checks and calls Ore::generate
which is responsible for spawing all ores of that type in mapblock
can call Ore::GenSingle
which spawns single cluster of the ore

Sub ores spawn only near or inside their parent ore.
Sub Ores generate in GenSingle method.

*/

#include "mg_ore.h"
#include "mapgen.h"
#include "noise.h"
#include "util/numeric.h"
#include "map.h"
#include "log.h"
#include "util/timetaker.h"

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
	parent = NULL;
}


Ore::~Ore()
{
	delete noise;
	//FIXME delete biomes?
}


void Ore::resolveNodeNames()
/*
	Called by NodeDefManager. Read resolved content ids.
	On sub ores replace wherein by parent's wherein and add parent's ore.
*/
{
	getIdFromNrBacklog(&c_ore, "", CONTENT_AIR);
	getIdsFromNrBacklog(&c_wherein);
	if(parent)
	{
		//delete c_wherein?
		c_wherein=parent->c_wherein;
		c_wherein.push_back(parent->c_ore);
		//FIXME delete biomes?
		biomes=parent->biomes;
	}
}

void Ore::generate(struct OreEnv env)
{
	return; //do nothing in case this is sub- only ore
}
void Ore::GenSingle(struct OreEnv env, v3s16 pA)
{
	return; //do nothing in case this has no sub- mode
}

size_t Ore::placeOre(Mapgen *mg, u32 blockseed, v3s16 nmin, v3s16 nmax)
/*
	Basically call generate method. But not on sub ores.
	Check y_min/max constraints, negate if using absheight
*/
{
	if(parent) return 0;
	int in_range = 0;
	struct OreEnv env;

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
	//FIXME move env initialization to placeAllOres
	env.vm=mg->vm;
  env.mapseed=mg->seed;
  env.blockseed=blockseed;
  env.nmin=nmin;
  env.nmax=nmax;
  env.biomemap=mg->biomemap;
  env.pr=NULL;
  env.size.X=(env.nmax.X - env.nmin.X + 1);
  env.size.Y=(env.nmax.Y - env.nmin.Y + 1);
  env.size.Z=(env.nmax.Z - env.nmin.Z + 1);
  env.volume=(nmax.X - nmin.X + 1) *
				 (nmax.Y - nmin.Y + 1) *
				 (nmax.Z - nmin.Z + 1);

	generate(env);

	return 1;
}

u32 Ore::GetOreGenCount(PcgRandom *pr,u32 volume,u32 scarcity)
/*
	Calculate ore clusters count in given volume
	Random contribution is <0;1)
*/
{
	float gencnt_raw = (float)volume / (float)scarcity;
  u32 gencnt = ((float)(gencnt_raw) + (float)((((float)pr->next())/((float)pr->RANDOM_RANGE))));
  return gencnt;
}

///////////////////////////////////////////////////////////////////////////////


void OreScatter::GenSingle(struct OreEnv env, v3s16 pA)
/*
	Generate cluster of ore. Cube of clust_size with only some
	nodes ore.
*/
{
	MapNode n_ore(c_ore, 0, ore_param2);
	u32 csize     = clust_size;
	u32 cvolume    = csize * csize * csize;

	if ((flags & OREFLAG_USE_NOISE) &&
		(NoisePerlin3D(&np, pA.X, pA.Y, pA.Z, env.mapseed) < nthresh))
		return;

	if (env.biomemap && !biomes.empty()) {
		u32 index = env.size.X * (pA.Z - env.nmin.Z) + (pA.X - env.nmin.X);
		std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
		if (it == biomes.end())
			return;
	}

	for (u32 z1 = 0; z1 != csize; z1++)
	for (u32 y1 = 0; y1 != csize; y1++)
	for (u32 x1 = 0; x1 != csize; x1++) {
		if (env.pr->range(1, cvolume) > clust_num_ores)
			continue;

		u32 i = env.vm->m_area.index(pA.X + x1, pA.Y + y1, pA.Z + z1);
		if (!env.vm->m_area.contains(i))
				continue;
		if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
			continue;

		env.vm->m_data[i] = n_ore;
	}
}

void OreScatter::generate(struct OreEnv env)
/*
	Place clusters randomly in mapblock
*/
{
	if(parent) return; //redundant
	PcgRandom pr(env.blockseed);
	env.pr=&pr;
	
	MapNode n_ore(c_ore, 0, ore_param2);
	u32 csize     = clust_size;
	u32 nclusters = GetOreGenCount(env.pr,env.volume,clust_scarcity);

	for (u32 i = 0; i != nclusters; i++) {
		int x0 = pr.range(env.nmin.X, env.nmax.X - csize + 1);
		int y0 = pr.range(env.nmin.Y, env.nmax.Y - csize + 1);
		int z0 = pr.range(env.nmin.Z, env.nmax.Z - csize + 1);

		if ((flags & OREFLAG_USE_NOISE) &&
			(NoisePerlin3D(&np, x0, y0, z0, env.mapseed) < nthresh))
			continue;

		GenSingle(env,v3s16(x0,y0,z0));
	}
}


///////////////////////////////////////////////////////////////////////////////


void OreSheet::generate(struct OreEnv env)
{
	PcgRandom pr(env.blockseed + 4234);
	MapNode n_ore(c_ore, 0, ore_param2);

	u16 max_height = column_height_max;
	int y_start_min = env.nmin.Y + max_height;
	int y_start_max = env.nmax.Y - max_height;

	int y_start = y_start_min < y_start_max ?
		pr.range(y_start_min, y_start_max) :
		(y_start_min + y_start_max) / 2;

	if (!noise) {
		noise = new Noise(&np, 0, env.size.X, env.size.Z);
	}
	noise->seed = env.mapseed + y_start;
	noise->perlinMap2D(env.nmin.X, env.nmin.Z);

	size_t index = 0;
	for (int z = env.nmin.Z; z <= env.nmax.Z; z++)
	for (int x = env.nmin.X; x <= env.nmax.X; x++, index++) {
		float noiseval = noise->result[index];
		if (noiseval < nthresh)
			continue;

		if (env.biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		u16 height = pr.range(column_height_min, column_height_max);
		int ymidpoint = y_start + noiseval;
		int y0 = MYMAX(env.nmin.Y, ymidpoint - height * (1 - column_midpoint_factor));
		int y1 = MYMIN(env.nmax.Y, y0 + height - 1);

		for (int y = y0; y <= y1; y++) {
			u32 i = env.vm->m_area.index(x, y, z);
			if (!env.vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
				continue;

			env.vm->m_data[i] = n_ore;
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


void OrePuff::generate(struct OreEnv env)
{
	PcgRandom pr(env.blockseed + 4234);
	MapNode n_ore(c_ore, 0, ore_param2);

	int y_start = pr.range(env.nmin.Y, env.nmax.Y);

	if (!noise) {
		noise = new Noise(&np, 0, env.size.X, env.size.X);
		noise_puff_top = new Noise(&np_puff_top, 0, env.size.X, env.size.Z);
		noise_puff_bottom = new Noise(&np_puff_bottom, 0, env.size.X, env.size.Z);
	}

	noise->seed = env.mapseed + y_start;
	noise->perlinMap2D(env.nmin.X, env.nmin.Z);
	bool noise_generated = false;

	size_t index = 0;
	for (int z = env.nmin.Z; z <= env.nmax.Z; z++)
	for (int x = env.nmin.X; x <= env.nmax.X; x++, index++) {
		float noiseval = noise->result[index];
		if (noiseval < nthresh)
			continue;

		if (env.biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		if (!noise_generated) {
			noise_generated = true;
			noise_puff_top->perlinMap2D(env.nmin.X, env.nmin.Z);
			noise_puff_bottom->perlinMap2D(env.nmin.X, env.nmin.Z);
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
			u32 i = env.vm->m_area.index(x, y, z);
			if (!env.vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
				continue;

			env.vm->m_data[i] = n_ore;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////


void OreBlob::generate(struct OreEnv env)
{
	PcgRandom pr(env.blockseed + 2404);
	env.pr=&pr;
	MapNode n_ore(c_ore, 0, ore_param2);

	u32 csize  = clust_size;
	u32 nblobs = env.volume / clust_scarcity;

	if (!noise)
		noise = new Noise(&np, env.mapseed, csize, csize, csize);

	for (u32 i = 0; i != nblobs; i++) {
		int x0 = pr.range(env.nmin.X, env.nmax.X - csize + 1);
		int y0 = pr.range(env.nmin.Y, env.nmax.Y - csize + 1);
		int z0 = pr.range(env.nmin.Z, env.nmax.Z - csize + 1);

		if (env.biomemap && !biomes.empty()) {
			u32 bmapidx = env.size.X * (z0 - env.nmin.Z) + (x0 - env.nmin.X);
			std::set<u8>::iterator it = biomes.find(env.biomemap[bmapidx]);
			if (it == biomes.end())
				continue;
		}

		bool noise_generated = false;
		noise->seed = env.blockseed + i;

		size_t index = 0;
		for (u32 z1 = 0; z1 != csize; z1++)
		for (u32 y1 = 0; y1 != csize; y1++)
		for (u32 x1 = 0; x1 != csize; x1++, index++) {
			u32 i = env.vm->m_area.index(x0 + x1, y0 + y1, z0 + z1);
			if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
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

			env.vm->m_data[i] = n_ore;
		}
		for (auto it = sub_ores.begin() ; it != sub_ores.end(); ++it)
		{
			Ore &sub = **it;
			u32 gencnt = GetOreGenCount(env.pr,csize*csize*csize,sub.clust_scarcity);
			for(;gencnt>0;gencnt--)
			{
				v3s16 p;
				p.X=pr.range(csize)+x0;
				p.Y=pr.range(csize)+y0;
				p.Z=pr.range(csize)+z0;
				sub.GenSingle(env,p);
			}
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


void OreVein::generate(struct OreEnv env)
{
	PcgRandom pr(env.blockseed + 520);
	MapNode n_ore(c_ore, 0, ore_param2);


	if (!noise) {
		int sx = env.nmax.X - env.nmin.X + 1;
		int sy = env.nmax.Y - env.nmin.Y + 1;
		int sz = env.nmax.Z - env.nmin.Z + 1;
		noise  = new Noise(&np, env.mapseed, sx, sy, sz);
		noise2 = new Noise(&np, env.mapseed + 436, sx, sy, sz);
	}
	bool noise_generated = false;

	size_t index = 0;
	for (int z = env.nmin.Z; z <= env.nmax.Z; z++)
	for (int y = env.nmin.Y; y <= env.nmax.Y; y++)
	for (int x = env.nmin.X; x <= env.nmax.X; x++, index++) {
		u32 i = env.vm->m_area.index(x, y, z);
		if (!env.vm->m_area.contains(i))
			continue;
		if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
			continue;

		if (env.biomemap && !biomes.empty()) {
			u32 bmapidx = env.size.X * (z - env.nmin.Z) + (x - env.nmin.X);
			std::set<u8>::iterator it = biomes.find(env.biomemap[bmapidx]);
			if (it == biomes.end())
				continue;
		}

		// Same lazy generation optimization as in OreBlob
		if (!noise_generated) {
			noise_generated = true;
			noise->perlinMap3D(env.nmin.X, env.nmin.Y, env.nmin.Z);
			noise2->perlinMap3D(env.nmin.X, env.nmin.Y, env.nmin.Z);
		}

		// randval ranges from -1..1
		float randval   = (float)pr.next() / (pr.RANDOM_RANGE / 2) - 1.f;
		float noiseval  = contour(noise->result[index]);
		float noiseval2 = contour(noise2->result[index]);
		if (noiseval * noiseval2 + randval * random_factor < nthresh)
			continue;

		env.vm->m_data[i] = n_ore;
	}
}

///////////////////////////////////////////////////////////////////////////////

static bool GetLine(v3s16 nmin, v3s16 nmax, PcgRandom *pr,
   short length, short length_rnd, float theta, float theta_rnd,
   v3f &start, v3f &dirv)
/*
	Calculate random line inside nmin..nmax
	resuli is start(point) and dirv(vector)
*/
{
	//*select r, theta, phi*
	float dtha = theta + ((float)pr->next()/pr->RANDOM_RANGE)*theta_rnd;
	float dpha = ((float)pr->next()/pr->RANDOM_RANGE)*3.1416;
	float lena = length + ((float)pr->next()*length_rnd/pr->RANDOM_RANGE);
	v3f dir_nv ( sin(dpha)*cos(dtha), sin(dpha)*sin(dtha), cos(dpha) );
	dirv = dir_nv * lena;
	//*select start point*
    //FIXME: ensure min<max !
	nmax.X -= dirv.X -1;
	nmax.Y -= dirv.Y -1;
	if(dirv.Z>0)
		{nmax.Z -= dirv.Z -1;}
		else
		{nmin.Z -= dirv.Z -1;}
	if((nmin.X>nmax.X)||(nmin.Y>nmax.Y)||(nmin.Z>nmax.Z)) return false;
	start.X = pr->range(nmin.X, nmax.X);
	start.Y = pr->range(nmin.Y, nmax.Y);
	start.Z = pr->range(nmin.Z, nmax.Z);
	return true;
}

void OrePipe::placePipe(struct OreEnv env,
   v3f pA, v3f pB, v3f pC)
/*
	Spawn a single "pipe" between points A and B with control point C
*/
{
	TimeTaker tt ("OrePipe place pipe with all sub-dists",NULL,PRECISION_MICRO);
	MapNode n_ore(c_ore, 0, ore_param2);
  v3f sa_p=pA;
  v3f p;
 	u32 volume = 0;
  //*estimate a step size*
  float step=1.7/(pA.getDistanceFrom(pB)+pB.getDistanceFrom(pC));
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
      if (sa_e>=1.7)
        { step=step*0.50; t=t-step; }
      else if (sa_e<0.0)
        { t=t+(step/2); step=step+(step/2); }
      else break;
      //printf("step_adjust %d %f %f \n",sa_c,sa_e,step);
    }
    //p.X=round(p.X);p.Y=round(p.Y);p.Z=round(p.Z);
    sa_p=p;

    for(int bz=-pipe_radius; bz<=pipe_radius; bz++)
    for(int by=-pipe_radius; by<=pipe_radius; by++)
    for(int bx=-pipe_radius; bx<=pipe_radius; bx++)
    {

			//*convert float[3] to area index*
			u32 di = env.vm->m_area.index(p.X+bx, p.Y+by, p.Z+bz);

			//*check can place*
			if (!env.vm->m_area.contains(di))
				continue;
			if(((bx*bx)+(by*by)+(bz*bz)-(pipe_radius*pipe_radius))>0)
				continue;
			if (!CONTAINS(c_wherein, env.vm->m_data[di].getContent()))
				continue;

			//*place the actual pipe segment*
			env.vm->m_data[di] = n_ore;
			volume++;
		}
  }
  //printf("no of sub-dists in this %s pipe dist %u volume: %lu\n",name.c_str(),volume,sub_ores.size());
  for (auto it = sub_ores.begin() ; it != sub_ores.end(); ++it)
  {
		Ore &sub = **it;
		u32 gencnt = GetOreGenCount(env.pr,volume,sub.clust_scarcity);
		//printf("this pipe dist: %u of %s\n",gencnt,sub.name.c_str());
		for(;gencnt>0;gencnt--)
		{
			float t = ((float)env.pr->next()/(env.pr->RANDOM_RANGE));
			p  = pA * ((1-t)*(1-t));
      p += pC * (2*t*(1-t));
      p += pB * (t*t);
			sub.GenSingle(env,
				v3s16(p.X+env.pr->range(-pipe_radius,+pipe_radius),
				   p.Y+env.pr->range(-pipe_radius,+pipe_radius),
				   p.Z+env.pr->range(-pipe_radius,+pipe_radius)
				)
		  );
		}
	}
  tt.stop(false);
}

void OrePipe::generate(struct OreEnv env)
/*
  Generate all "pipe" ores in current mg-block
  not using noise
*/
{
	PcgRandom pr(env.blockseed);
	env.pr=&pr;
	TimeTaker tt ("OrePipe generate",NULL,PRECISION_MICRO);

	//*calc count of pipes in this block*
  u32 gencnt = GetOreGenCount(env.pr,env.volume,clust_scarcity);
  //printf("gencnt: %d..%d %d\n", nmin.X,nmax.X, gencnt);
	for (u32 i = 0; i != gencnt; i++) {
    v3f A, B, C, dirv;

    //*select start and direction vectors*
    if(!
    GetLine(env.nmin, env.nmax, env.pr,
		   pipe_length, pipe_length_rnd, dir_theta, dir_theta_rnd,
       A, dirv)) break; //or continue?

		//*calc end of pipe*
		B = A + dirv;

		//*calc curve control point*
    int curv = ceil( curving * pipe_length * 0.5 );
    C = A + (dirv*0.5);
    C.X += pr.range(-curv, +curv);
    C.Y += pr.range(-curv, +curv);
    C.Z += pr.range(-curv, +curv);

		if (env.biomemap && !biomes.empty()) {
			u32 index = env.size.X * (A.Z - env.nmin.Z) + (A.X - env.nmin.X);
			std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
			if (it == biomes.end())
				continue;
		}

    //printf("place (%.0f,%.0f,%.0f),(%.0f,%.0f,%.0f),(%.0f,%.0f,%.0f)\n",A.X,A.Y,A.Z,C.X,C.Y,C.Z,B.X,B.Y,B.Z);
    placePipe(env, A, B, C);
	}
	tt.stop(false);
}

///////////////////////////////////////////////////////////////////////////////


void OreLayer::generate(struct OreEnv env)
/*
	Solid layer of ore between y_min and y_max,
	respects biomes and noise.
*/
{
	PcgRandom pr(env.blockseed + 873);
	MapNode n_ore(c_ore, 0, ore_param2);

	if (flags & OREFLAG_USE_NOISE) {
		if (!noise) {
			noise = new Noise(&np, 0, env.size.X, env.size.Z);
		}
		noise->seed = env.mapseed;
		noise->perlinMap2D(env.nmin.X, env.nmin.Z);
	}

	size_t index = 0;
	for (int z = env.nmin.Z; z <= env.nmax.Z; z++)
	for (int x = env.nmin.X; x <= env.nmax.X; x++, index++) {

		if ((flags & OREFLAG_USE_NOISE)&&(noise->result[index] < 0))
			continue;

		if (env.biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		for (int y = env.nmin.Y; y <= env.nmax.Y; y++) {
			u32 i = env.vm->m_area.index(x, y, z);
			if (!env.vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
				continue;

			env.vm->m_data[i] = n_ore;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

OreRegion::OreRegion() :
	Ore()
{
	noise_a = NULL;
	noise_b = NULL;
}


OreRegion::~OreRegion()
{
	delete noise_b;
	delete noise_a;
}

void OreRegion::resolveNodeNames()
{
	getIdFromNrBacklog(&c_ore2, "", CONTENT_IGNORE);
	getIdFromNrBacklog(&c_ore3, "", CONTENT_IGNORE);
	getIdFromNrBacklog(&c_ore4, "", c_ore2);
	Ore::resolveNodeNames();
	if(c_ore2==CONTENT_IGNORE) c_ore2=c_ore;
	if(c_ore3==CONTENT_IGNORE) c_ore3=c_ore2;
}

void OreRegion::generate(struct OreEnv env)
{
	PcgRandom pr(env.blockseed + 873);
	MapNode n_ore[4]={
		MapNode(c_ore,  0, ore_param2),
		MapNode(c_ore2, 0, ore_param2),
		MapNode(c_ore3, 0, ore_param2),
		MapNode(c_ore4, 0, ore_param2)
	};

	if (flags & OREFLAG_USE_NOISE) {
		if (!noise) {
			noise = new Noise(&np, 0, env.size.X, env.size.Z);
		}
		noise->seed = env.mapseed;
		noise->perlinMap2D(env.nmin.X, env.nmin.Z);
	}
	if(!noise_a) { //this is not a good solution
		noise_a = new Noise(&np_region, 0, env.size.X, env.size.Z);
		noise_b = new Noise(&np_region, 7782, env.size.X, env.size.Z);
	}
	noise_a->perlinMap2D(env.nmin.X, env.nmin.Z);
	noise_b->perlinMap2D(env.nmin.X, env.nmin.Z);

	size_t index = 0;
	for (int z = env.nmin.Z; z <= env.nmax.Z; z++)
	for (int x = env.nmin.X; x <= env.nmax.X; x++, index++) {

		if (env.biomemap && !biomes.empty()) {
			std::set<u8>::iterator it = biomes.find(env.biomemap[index]);
			if (it == biomes.end())
				continue;
		}

		unsigned short region=0;
		if(noise_a->result[index]>0) region|=2;
		if(noise_b->result[index]>0) region|=1;

		//TODO use noise for y (min=min+nv, max=max-nv)?

		for (int y = env.nmin.Y; y <= env.nmax.Y; y++) {
			u32 i = env.vm->m_area.index(x, y, z);

			if (!env.vm->m_area.contains(i))
				continue;
			if (!CONTAINS(c_wherein, env.vm->m_data[i].getContent()))
				continue;

			env.vm->m_data[i] = n_ore[region];

		}
	}
}

///////////////////////////////////////////////////////////////////////////////
