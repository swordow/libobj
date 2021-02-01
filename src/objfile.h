/*
*	swordow, 2021/01/29 
*/

#ifndef __OBJECT_FILE_H__
#define __OBJECT_FILE_H__
#include <stdint.h>

struct ObjMtl
{
	char Name[256];
	char MapKd[256];
	char MapBump[256];
	ObjMtl* Next;
};

struct ObjMtlLib
{
	char Name[256];
	ObjMtl* Mtls;
};

struct ObjVector4
{
	union {
		float v[4];
		struct
		{
			float x, y, z, w;
		};
	};
};

struct ObjVector3
{
	union {
		float v[3];
		struct
		{
			float x, y, z;
		};
	};
};

struct ObjVector2
{
	float x, y;
};

struct ObjFace
{
	union
	{
		int32_t P[3];
		struct
		{
			int32_t Pa, Pb, Pc;
		};
	};
	ObjVector3 N;
};

struct ObjMesh
{
	char		Name[256];
	int32_t		FCount;
	int32_t		VCount;
	int32_t		UVLayerCount;
	ObjVector3* Points;
	ObjVector4* UVs[8];
	ObjFace*	Faces;
	ObjMtl*		Mtl;
	ObjMesh*	Next;
};

struct ObjGroup
{
	char	  Name[256];
	ObjMesh*  MeshList;
	ObjGroup* Next;
};

struct ObjectFile
{
	char	   File[256];
	ObjGroup*  GroupList;
	ObjMtlLib  MtlLib;
};

ObjectFile* LoadObjectFile(const char* file);

#endif
