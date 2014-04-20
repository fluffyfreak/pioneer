#ifndef _VOLUMETRIC_H
#define _VOLUMETRIC_H

#include "libs.h"

#define SORT_TOWARD 0
#define SORT_AWAY 1

struct CloudPuff
{
	float 		Size;
	int 		ID;
	float 		Angle;
	vector3f	Position;
	float		DistanceToCam;
	Color4f		Color;
	float		Life;	
};

struct VolumetricCloud
{
	std::vector<CloudPuff> Puffs;
	unsigned 	ImpostorTex;
	vector3f	BoundingBox1, BoundingBox2;
	vector3f 	Center;
	float		Radius;
	vector3f*	VertexBuffer;
	Color4f*	ColorBuffer;
	vector2f*	TexCoordBuffer;
	
	vector3f	LastCamera;
	vector3f	LastLight;
	float		DistanceFromCamera;
	int			ImpostorSize;
	vector3f	vx, vy;	//up and right vectors for impostor rendering
};

class VolumetricClouds
{
public:
	VolumetricClouds();
	
	int 	Create(int NumClouds, float PlaneSize, float PlaneHeight);	
	void	Update(vector3f Sun, vector3f Camera);
	void	Render(vector3f Sun, vector3f Camera);
	void	Destroy();
	void	GetInfo(int* Sprites, int* Impostors)
	{
		*Sprites = NumSprites;
		*Impostors = NumImpostors;
	}

private:

	void	UpdateCloud(VolumetricCloud* Cloud, vector3f Sun, vector3f Camera);
	void	RenderCloudImpostor(VolumetricCloud* Cloud, float alpha); 
	void	RenderCloud3D(VolumetricCloud* Cloud, vector3f Camera, vector3f Sun, float alpha);

	void	MakeCloudImpostor(VolumetricCloud* Cloud, vector3f Sun, vector3f Camera);
	
	void	LightCloud(VolumetricCloud* Cloud, vector3f Sun);
	void	GrowCloud(VolumetricCloud* Cloud, int level, float radius, vector3f Position);

	int		GetImpostorSize(float distance2);	
	void	GenerateTexture();
	void	SortParticles(VolumetricCloud* Cloud, int mode);

	class	SortAwayComparison
	{
	public:
		bool operator () (CloudPuff puff1, CloudPuff puff2)
		{
			return puff1.DistanceToCam < puff2.DistanceToCam;
		}
	} SortAway;

	class	SortTowardComparison
	{
	public:
		bool operator () (CloudPuff puff1, CloudPuff puff2)
		{
			return puff1.DistanceToCam > puff2.DistanceToCam;
		}
	} SortToward;

	class	SortCloudsTowardComparison
	{
	public:
		bool operator () (VolumetricCloud cloud1, VolumetricCloud cloud2)
		{
			return cloud1.DistanceFromCamera > cloud2.DistanceFromCamera;
		};
	} SortCloudToward;

	unsigned PuffTexture;
	unsigned PuffImage;

	std::vector<VolumetricCloud> Clouds;
	
	int	SplatBufferSize;
	int ImpostorSize;
	int NumSprites, NumImpostors;
	float Albedo, Extinction;
};

#endif
