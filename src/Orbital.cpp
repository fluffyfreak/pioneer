#include "libs.h"
#include "Orbital.h"
#include "Frame.h"
#include "Pi.h"
#include "WorldView.h"
#include "Serializer.h"
#include "GeoRing.h"
#include "graphics/Renderer.h"
#include "perlin.h"
#include "GameSaveError.h"
#include "JsonUtils.h"

Orbital::Orbital() : Body(),
	m_mass(0),
	pos(0.0),
	sbody(nullptr),
	m_geoRing(nullptr)
{
}

Orbital::Orbital(SystemBody *sbody_): Body(),
	m_mass(0),
	pos(0.0),
	sbody(sbody_),
	m_geoRing(nullptr)
{
	Init();
}

void Orbital::Init()
{
	m_mass = sbody->GetMass();
	if (!m_geoRing) {
		m_geoRing = new GeoRing(sbody);
	}
	SetClipRadius(sbody->GetRadius());
}

bool Orbital::IsSuperType(SystemBody::BodySuperType t) const
{
	if (!sbody) return false;
	else return sbody->GetSuperType() == t;
}
	
void Orbital::SaveToJson(Json::Value &jsonObj, Space *space)
{
	Body::SaveToJson(jsonObj, space);
	
	Json::Value orbitalObj(Json::objectValue); // Create JSON object to contain projectile data.
	VectorToJson(orbitalObj, pos, "pos");
	orbitalObj["sbodyIdx"] = space->GetIndexForSystemBody(sbody);

	jsonObj["orbital"] = orbitalObj; // Add terrain body object to supplied object.
}

void Orbital::LoadFromJson(const Json::Value &jsonObj, Space *space)
{
	Body::LoadFromJson(jsonObj, space);

	if (!jsonObj.isMember("orbital")) throw SavedGameCorruptException();
	Json::Value orbitalObj = jsonObj["orbital"];

	if (!orbitalObj.isMember("sbodyIdx")) throw SavedGameCorruptException();

	JsonToVector(&pos, orbitalObj, "pos");
	sbody = space->GetSystemBodyByIndex(orbitalObj["sbodyIdx"].asInt());
	assert(sbody != nullptr);

	Init();
}

Orbital::~Orbital()
{
	if (m_geoRing) delete m_geoRing;
}

double Orbital::GetBoundingRadius() const
{
	// needs to include all terrain so culling works {in the future},
	// and size of rotating frame is correct
	return sbody->GetRadius() * (1.1+m_geoRing->GetMaxFeatureHeight());
}

vector3d Orbital::GetPosition() const
{
	return pos;
}

void Orbital::SetPosition(vector3d p)
{
	pos = p;
}

/*
 * dist = distance from centre
 * returns pressure in earth atmospheres
 */
void Orbital::GetAtmosphericState(double dist, double *outPressure, double *outDensity)
{
	Color c;
	double surfaceDensity;
	double atmosDist = dist/(sbody->GetRadius()*ATMOSPHERE_RADIUS);
	
	m_geoRing->GetAtmosphereFlavor(&c, &surfaceDensity);
	// kg / m^3
	// exp term should be the same as in AtmosLengthDensityProduct GLSL function
	*outDensity = 1.15*surfaceDensity * exp(-500.0 * (atmosDist - (2.0 - ATMOSPHERE_RADIUS)));
	// XXX using earth's molar mass of air...
	const double GAS_MOLAR_MASS = 28.97;
	const double GAS_CONSTANT = 8.314;
	const double KPA_2_ATMOS = 1.0 / 101.325;
	// atmospheres
	*outPressure = KPA_2_ATMOS*(*outDensity/GAS_MOLAR_MASS)*GAS_CONSTANT*double(sbody->GetAverageTemp());
}

double Orbital::GetTerrainHeight(const vector3d pos) const
{
	double radius = sbody->GetRadius();
	if (m_geoRing) {
		return radius * (1.0 - m_geoRing->GetHeight(pos));//GetDistFromSurface(pos));
	} else {
		assert(0);
		return radius;
	}
}

double Orbital::GetRingWidth() const 
{ 
	double ret=0.0; 
	(m_geoRing) ? ret=m_geoRing->GetRingWidth() : ret=0.0; 
	return ret; 
}

void Orbital::Render(Graphics::Renderer *renderer, const Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform)
{
	matrix4x4d ftran = viewTransform;
	vector3d fpos = viewCoords;
	double rad = sbody->GetRadius();

	float znear, zfar;
	Pi::renderer->GetNearFarRange(znear, zfar);

	vector3d campos = fpos;
	ftran.ClearToRotOnly();
	campos = ftran.Inverse() * campos;

	campos = campos * (1.0/rad);		// position of camera relative to planet "model"

	std::vector<Camera::Shadow> shadows;
	if( camera ) {
		camera->PrincipalShadows(this, 3, shadows);
		for (std::vector<Camera::Shadow>::iterator it = shadows.begin(), itEnd=shadows.end(); it!=itEnd; ++it) {
			it->centre = ftran * it->centre;
		}
	}

	ftran.Scale(rad, rad, rad);

	// translation not applied until patch render to fix jitter
	m_geoRing->Render(renderer, ftran, -campos, sbody->GetRadius(), shadows);

	ftran.Translate(campos.x, campos.y, campos.z);
	//SubRender(renderer, ftran, campos);
}

void Orbital::SetFrame(Frame *f)
{
	if (GetFrame()) {
		GetFrame()->SetPlanetGeom(0, 0);
	}
	Body::SetFrame(f);
	if (f) {
		GetFrame()->SetPlanetGeom(0, 0);
	}
}

