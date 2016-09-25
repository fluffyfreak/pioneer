// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GEOSPHEREANALYSE_H
#define _GEOSPHEREANALYSE_H

#include "JobQueue.h"

class GeoSphere;

void Analyse(GeoSphere *geo);

// ********************************************************************************
// Overloaded PureJob class to handle analysing each GeoSphere
// ********************************************************************************
class AnalyseJob : public Job
{
public:
	AnalyseJob() {};
	AnalyseJob(GeoSphere *geoSphere)
		: m_geoSphere(geoSphere) {}

	virtual void OnRun() override final { PROFILE_SCOPED_DESC("GeoSphere-AnalyseJob"); Analyse(m_geoSphere); }    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() override final {}
	virtual void OnCancel() override final {}

protected:
	GeoSphere	*m_geoSphere;
};

#endif /* _GEOSPHEREANALYSE_H */
