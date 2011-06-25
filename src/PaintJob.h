#ifndef __PAINTJOB_H__
#define __PAINTJOB_H__

#include "Lengyel.h"

class CPaintJob
{
public:
	CPaintJob() {}
	~CPaintJob() {}

	void AddDecal(const vector3d& center, const vector3d& normal, const vector3d& tangent, float width, float height, float depth);
private:
	std::vector<Decal> m_decals;
};

#endif // __PAINTJOB_H__