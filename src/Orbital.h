#ifndef _ORBITAL_H
#define _ORBITAL_H

#include "Body.h"
// only for SBody::BodySuperType enum...
#include "galaxy/StarSystem.h"

class Frame;
class SystemBody;
class GeoRing;

class Orbital: public Body {
public:
	OBJDEF(Orbital, Body, ORBITAL);

	Orbital(SystemBody*);
	Orbital();
	virtual ~Orbital();

	virtual void SetPosition(vector3d p);
	virtual vector3d GetPosition() const;
	virtual double GetBoundingRadius() const;
	virtual void Render(Graphics::Renderer *r, const Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform) override final;
	virtual void SetFrame(Frame *f) override final;
	virtual bool OnCollision(Object *b, Uint32 flags, double relVel)  override final { return true; }
	virtual double GetMass() const  override final { return m_mass; }
	virtual const SystemBody *GetSystemBody() const  override final { return sbody; }

	double GetTerrainHeight(const vector3d pos) const;
	void GetAtmosphericState(double dist, double *outPressure, double *outDensity);
	bool IsSuperType(SystemBody::BodySuperType t) const;

#if WITH_OBJECTVIEWER
	friend class ObjectViewerView;
#endif

	double GetRingWidth() const;

protected:
	virtual void SaveToJson(Json::Value &jsonObj, Space *space) override;
	virtual void LoadFromJson(const Json::Value &jsonObj, Space *space) override;

private:
	void Init();

	double m_mass;
	vector3d pos;
	SystemBody *sbody;
	GeoRing *m_geoRing;
};

#endif /* _ORBITAL_H */
