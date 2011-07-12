#include "libs.h"
#include "Orbital.h"
#include "Frame.h"
#include "Pi.h"
#include "WorldView.h"
#include "Serializer.h"
#include "StarSystem.h"
#include "GeoRing.h"
#include "Render.h"
#include "perlin.h"

Orbital::Orbital(): Body()
{
	pos = vector3d(0,0,0);
	this->sbody = 0;
	this->m_geoRing = 0;
}

Orbital::Orbital(SBody *sbody_): Body()
{
	pos = vector3d(0,0,0);
	this->sbody = sbody_;
	this->m_geoRing = 0;
	Init();
	m_hasDoubleFrame = true;
}

void Orbital::Init()
{
	m_mass = sbody->GetMass();
	if (!m_geoRing) {
		m_geoRing = new GeoRing(sbody);
	}
}

bool Orbital::IsSuperType(SBody::BodySuperType t) const
{
	if (!sbody) return false;
	else return sbody->GetSuperType() == t;
}
	
void Orbital::Save(Serializer::Writer &wr)
{
	Body::Save(wr);
	wr.Vector3d(pos);
	wr.Int32(Serializer::LookupSystemBody(sbody));
}

void Orbital::Load(Serializer::Reader &rd)
{
	Body::Load(rd);
	pos = rd.Vector3d();
	sbody = Serializer::LookupSystemBody(rd.Int32());
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

void Orbital::SetRadius(double radius)
{
	assert(0);
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
	*outPressure = KPA_2_ATMOS*(*outDensity/GAS_MOLAR_MASS)*GAS_CONSTANT*double(sbody->averageTemp);
}

double Orbital::GetTerrainHeight(const vector3d pos_) const
{
	double radius = sbody->GetRadius();
	if (m_geoRing) {
		double height = -(m_geoRing->GetDistFromSurface(pos_, radius));
		//*(vts++) = p * (height + 1.0);
		return radius * (1.0 + height);
	} else {
		assert(0);
		return radius;
	}
}


#define PLANET_AMBIENT	0.1f

static void SetMaterialColor(const float col[4])
{
	float mambient[4];
	mambient[0] = col[0]*PLANET_AMBIENT;
	mambient[1] = col[1]*PLANET_AMBIENT;
	mambient[2] = col[2]*PLANET_AMBIENT;
	mambient[3] = col[3];
	glMaterialfv (GL_FRONT, GL_AMBIENT, mambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, col);
}

static void _DrawAtmosphere(double rad1, double rad2, vector3d &pos, const float col[4])
{
	glPushMatrix();
	// face the camera dammit
	vector3d zaxis = (-pos).Normalized();
	vector3d xaxis = vector3d(0,1,0).Cross(zaxis).Normalized();
	vector3d yaxis = zaxis.Cross(xaxis);
	matrix4x4d rot = matrix4x4d::MakeInvRotMatrix(xaxis, yaxis, zaxis);
	glMultMatrixd(&rot[0]);

	matrix4x4f invViewRot;
	glGetFloatv(GL_MODELVIEW_MATRIX, &invViewRot[0]);
	invViewRot.ClearToRotOnly();
	invViewRot = invViewRot.InverseOf();

	const int numLights = Pi::worldView->GetNumLights();
	assert(numLights < 4);
	vector3f lightDir[4];
	float lightCol[4][4];
	// only 
	for (int i=0; i<numLights; i++) {
		float temp[4];
		glGetLightfv(GL_LIGHT0 + i, GL_DIFFUSE, lightCol[i]);
		glGetLightfv(GL_LIGHT0 + i, GL_POSITION, temp);
		lightDir[i] = (invViewRot * vector3f(temp[0], temp[1], temp[2])).Normalized();
	}

	const double angStep = M_PI/32;
	// find angle player -> centre -> tangent point
	// tangent is from player to surface of sphere
	float tanAng = float(acos(rad1 / pos.Length()));

	// then we can put the fucking atmosphere on the horizon
	vector3d r1(0.0, 0.0, rad1);
	vector3d r2(0.0, 0.0, rad2);
	rot = matrix4x4d::RotateYMatrix(tanAng);
	r1 = rot * r1;
	r2 = rot * r2;

	rot = matrix4x4d::RotateZMatrix(angStep);

	glDisable(GL_LIGHTING);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);	
	glEnable(GL_BLEND);
	glDisable(GL_CULL_FACE);
	glBegin(GL_TRIANGLE_STRIP);
	for (float ang=0; ang<2*M_PI; ang+=float(angStep)) {
		vector3d norm = r1.Normalized();
		glNormal3dv(&norm.x);
		float _col[4] = { 0,0,0,0 };
		for (int lnum=0; lnum<numLights; lnum++) {
			const float dot = norm.x*lightDir[lnum].x + norm.y*lightDir[lnum].y + norm.z*lightDir[lnum].z;
			_col[0] += dot*lightCol[lnum][0];
			_col[1] += dot*lightCol[lnum][1];
			_col[2] += dot*lightCol[lnum][2];
		}
		for (int i=0; i<3; i++) _col[i] = _col[i] * col[i];
		_col[3] = col[3];
		glColor4fv(_col);
		glVertex3dv(&r1.x);
		glColor4f(0,0,0,0);
		glVertex3dv(&r2.x);
		r1 = rot * r1;
		r2 = rot * r2;
	}
	
	glEnd();
	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glPopMatrix();
}

void Orbital::DrawAtmosphere(vector3d &apos)
{
	Color c;
	double density;
	m_geoRing->GetAtmosphereFlavor(&c, &density);
	
	_DrawAtmosphere(0.999, 1.05, apos, c);
}

void Orbital::Render(const vector3d &viewCoords, const matrix4x4d &viewTransform)
{

	matrix4x4d ftran = viewTransform;
	vector3d fpos = viewCoords;
	double rad = sbody->GetRadius();

	float znear, zfar;
	Pi::worldView->GetNearFarClipPlane(&znear, &zfar);

	double apparent_size = rad / fpos.Length();
	double len = fpos.Length();
	double origLen = len;
	int shrink = 0;
	double scale = 1.0f;

	double dist_to_horizon;
	for (;;) {
		if (len < rad) break;		// player inside radius case
		dist_to_horizon = sqrt(len*len - rad*rad);

		if (dist_to_horizon < zfar*0.5) break;

		rad *= 0.25;
		fpos = 0.25*fpos;
		len *= 0.25;
		scale *= 4.0f;
		shrink++;
	}
	//if (GetLabel() == "Earth") printf("Horizon %fkm, shrink %d\n", dist_to_horizon*0.001, shrink);

	glPushMatrix();
	glTranslatef(float(fpos.x), float(fpos.y), float(fpos.z));
	glColor3f(1,1,1);

	if (apparent_size < 0.001) {
		Render::State::UseProgram(0);
		/* XXX WRONG. need to pick light from appropriate turd. */
		GLfloat col[4];
		glGetLightfv(GL_LIGHT0, GL_DIFFUSE, col);
		// face the camera dammit
		vector3d zaxis = fpos.Normalized();
		vector3d xaxis = vector3d(0,1,0).Cross(zaxis).Normalized();
		vector3d yaxis = zaxis.Cross(xaxis);
		matrix4x4d rot = matrix4x4d::MakeInvRotMatrix(xaxis, yaxis, zaxis);
		glMultMatrixd(&rot[0]);

		glDisable(GL_LIGHTING);
		glDisable(GL_DEPTH_TEST);
		
		glEnable(GL_BLEND);
		glColor4f(col[0], col[1], col[2], 1);
		glBegin(GL_TRIANGLE_FAN);
		glVertex3f(0,0,0);
		glColor4f(col[0],col[1],col[2],0);
		
		const float spikerad = float(0.005*len +  1e1*(1.0*sbody->GetRadius()*len)/origLen);
		{
			/* bezier with (0,0,0) control points */
			vector3f p0(0,spikerad,0), p1(spikerad,0,0);
			float t=0.1f; for (int i=1; i<10; i++, t+= 0.1f) {
				vector3f p = (1-t)*(1-t)*p0 + t*t*p1;
				glVertex3fv(&p[0]);
			}
		}
		{
			vector3f p0(spikerad,0,0), p1(0,-spikerad,0);
			float t=0.1f; for (int i=1; i<10; i++, t+= 0.1f) {
				vector3f p = (1-t)*(1-t)*p0 + t*t*p1;
				glVertex3fv(&p[0]);
			}
		}
		{
			vector3f p0(0,-spikerad,0), p1(-spikerad,0,0);
			float t=0.1f; for (int i=1; i<10; i++, t+= 0.1f) {
				vector3f p = (1-t)*(1-t)*p0 + t*t*p1;
				glVertex3fv(&p[0]);
			}
		}
		{
			vector3f p0(-spikerad,0,0), p1(0,spikerad,0);
			float t=0.1f; for (int i=1; i<10; i++, t+= 0.1f) {
				vector3f p = (1-t)*(1-t)*p0 + t*t*p1;
				glVertex3fv(&p[0]);
			}
		}
		glEnd();
		glDisable(GL_BLEND);

		glEnable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
	} else {
		vector3d campos = -fpos;
		ftran.ClearToRotOnly();
		campos = ftran.InverseOf() * campos;
		glMultMatrixd(&ftran[0]);
		glEnable(GL_NORMALIZE);
		glPushMatrix();
		glScaled(rad, rad, rad);
		campos = campos * (1.0/rad);
		m_geoRing->Render(campos, sbody->GetRadius(), scale);
		
		fpos = ftran.InverseOf() * fpos;
		fpos *= (1.0/rad);
		if (!Render::AreShadersEnabled()) DrawAtmosphere(fpos);
		
		glPopMatrix();
		glDisable(GL_NORMALIZE);
		
		// if not using shader then z-buffer precision is hopeless and
		// we can't place objects on the terrain without awful z artifacts
		if (shrink || !Render::AreShadersEnabled()) {
			glClear(GL_DEPTH_BUFFER_BIT);
		}
	}
	glPopMatrix();
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

