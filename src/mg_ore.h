/*
Minetest
Copyright (C) 2010-2013 kwolekr, Ryan Kwolek <kwolekr@minetest.net>

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
*/

#ifndef MG_ORE_HEADER
#define MG_ORE_HEADER

#include "objdef.h"
#include "noise.h"
#include "nodedef.h"

class Noise;
class Mapgen;
class MMVManip;

/////////////////// Ore generation flags

#define OREFLAG_ABSHEIGHT     0x01
#define OREFLAG_PUFF_CLIFFS   0x02
#define OREFLAG_PUFF_ADDITIVE 0x04
#define OREFLAG_USE_NOISE     0x08

#define ORE_RANGE_ACTUAL 1
#define ORE_RANGE_MIRROR 2

enum OreType {
	ORE_SCATTER,
	ORE_SHEET,
	ORE_PUFF,
	ORE_BLOB,
	ORE_VEIN,
  ORE_PIPE,
  ORE_SUB_SCATTER,
};

extern FlagDesc flagdesc_ore[];

class OreSub;
class Ore : public ObjDef, public NodeResolver {
public:
	static const bool NEEDS_NOISE = false;

	content_t c_ore;                  // the node to place
	std::vector<content_t> c_wherein; // the nodes to be placed in
	u32 clust_scarcity; // ore cluster has a 1-in-clust_scarcity chance of appearing at a node
	s16 clust_num_ores; // how many ore nodes are in a chunk
	s16 clust_size;     // how large (in nodes) a chunk of ore is
	s16 y_min;
	s16 y_max;
	u8 ore_param2;		// to set node-specific attributes
	u32 flags;          // attributes for this ore
	float nthresh;      // threshold for noise at which an ore is placed
	NoiseParams np;     // noise for distribution of clusters (NULL for uniform scattering)
	Noise *noise;
	std::set<u8> biomes;
	std::vector<OreSub*> sub_ores;

	Ore();
	virtual ~Ore();

	virtual void resolveNodeNames();

	size_t placeOre(Mapgen *mg, u32 blockseed, v3s16 nmin, v3s16 nmax);
	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap) = 0;
};

class OreSub : public Ore {
public:
	OreSub();
	Ore *parent;
	virtual void resolveNodeNames();
  virtual void place(MMVManip *vm, int mapseed, u8 *biomemap,
    v3s16 nmin, v3s16 nmax, PcgRandom &pr,
    v3s16 pA) = 0;
};

class OreScatter : public OreSub {
public:
	static const bool NEEDS_NOISE = false;

	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
  virtual void place(MMVManip *vm, int mapseed, u8 *biomemap,
    v3s16 nmin, v3s16 nmax, PcgRandom &pr,
    v3s16 pA);
};

class OreSheet : public Ore {
public:
	static const bool NEEDS_NOISE = true;

	u16 column_height_min;
	u16 column_height_max;
	float column_midpoint_factor;

	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
};

class OrePuff : public Ore {
public:
	static const bool NEEDS_NOISE = true;

	NoiseParams np_puff_top;
	NoiseParams np_puff_bottom;
	Noise *noise_puff_top;
	Noise *noise_puff_bottom;

	OrePuff();
	virtual ~OrePuff();

	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
};

class OreBlob : public Ore {
public:
	static const bool NEEDS_NOISE = true;

	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
};

class OreVein : public Ore {
public:
	static const bool NEEDS_NOISE = true;

	float random_factor;
	Noise *noise2;

	OreVein();
	virtual ~OreVein();

	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
};

class OrePipe : public Ore {
public:
	static const bool NEEDS_NOISE = false;

	s16 pipe_radius;   // radius of pipe in nodes
  float pipe_length;   // min distance A-B
  float pipe_length_rnd; // add 0..this to length
  float dir_theta;   // min inclination in radians
  float dir_theta_rnd; // add 0..this to theta
  float curving; // 0 is strait line

  //clusters todo

  void placePipe(MMVManip *vm, int mapseed, u8 *biomemap,
    v3s16 nmin, v3s16 nmax, PcgRandom &pr,
    v3f pA, v3f pB, v3f pC);
	virtual void generate(MMVManip *vm, int mapseed, u32 blockseed,
		v3s16 nmin, v3s16 nmax, u8 *biomemap);
};

class OreManager : public ObjDefManager {
public:
	OreManager(IGameDef *gamedef);
	virtual ~OreManager() {}

	const char *getObjectTitle() const
	{
		return "ore";
	}

	static Ore *create(OreType type)
	{
		switch (type) {
		case ORE_SCATTER:
			return new OreScatter;
		case ORE_SHEET:
			return new OreSheet;
		case ORE_PUFF:
			return new OrePuff;
		case ORE_BLOB:
			return new OreBlob;
		case ORE_VEIN:
			return new OreVein;
		case ORE_PIPE:
			return new OrePipe;
		case ORE_SUB_SCATTER:
			return new OreScatter;
		default:
			return NULL;
		}
	}

	void clear();

	size_t placeAllOres(Mapgen *mg, u32 blockseed, v3s16 nmin, v3s16 nmax);
};

#endif
