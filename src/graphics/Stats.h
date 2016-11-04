// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _STATS_H
#define _STATS_H

#include "SDL_stdinc.h"
#include "RefCounted.h"
#include <utility>
#include <vector>

namespace Graphics {

class Stats : public RefCounted
{
public:
	static const Uint32 MAX_FRAMES_STORE = 30U;
	enum StatType {
		// renderer entries
		STAT_DRAWCALL = 0,
		STAT_DRAWTRIS,
		STAT_DRAWPOINTSPRITES,

		// buffers
		STAT_CREATE_BUFFER,
		STAT_DESTROY_BUFFER,

		// objects
		STAT_BUILDINGS,
		STAT_CITIES,
		STAT_GROUNDSTATIONS,
		STAT_SPACESTATIONS,
		STAT_ATMOSPHERES,
		STAT_PATCHES,
		STAT_PLANETS,
		STAT_GASGIANTS,
		STAT_STARS,
		STAT_SHIPS,

		// scenegraph entries
		STAT_BILLBOARD,

		MAX_STAT
	};
	
	enum StatTotalType {
		// Memory
		STAT_MEM_TEXTURE_MB = 0,
		STAT_MEM_TEXTURES,

		MAX_TOTAL_STAT
	};
	
	struct TFrameData {
		Uint32 m_stats[MAX_STAT];
	};
	
	struct TTotalData {
		Uint32 m_stats[MAX_TOTAL_STAT];
	};

	Stats();
	~Stats() {}

	void AddToStatCount(const StatType type, const Uint32 count);
	void AddToStatTotal(const StatTotalType type, const Sint32 count);
	void NextFrame();

	const TFrameData& FrameStats() const { return m_frameStats[m_currentFrame]; }
	const TFrameData& FrameStatsPrevious() const;
	
	const TTotalData& TotalStats() const { return m_totalStats; }

private:

	TFrameData m_frameStats[MAX_FRAMES_STORE];
	TTotalData m_totalStats;
	Uint32 m_currentFrame;
};

}

#endif
