/*
Copyright (c) 2009 David Bucciarelli (davibu@interfree.it)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include\\ply.h"

#define TINYOBJ_LOADER_C_IMPLEMENTATION

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#if defined(__linux__) || defined(__APPLE__) || defined(__ANDROID__)
#include <sys/time.h>
#elif defined (WIN32)
#include <windows.h>
#else
        Unsupported Platform !!!
#endif

#include "../assets/sdcard/include/camera.h"
#include "../assets/sdcard/include/geom.h"
#include "include/displayfunc.h"
#include "include/native-lib.h"
#include "../assets/sdcard/include/material.h"
#include "include/tinyobj_loader_c.h"

typedef struct {
	float gammaCorrection;
	Vec skyColor1, skyColor2;
	unsigned int materialCount;
} SceneTemp;

typedef enum { SPH, TRI, MOD } TypeTemp;
typedef struct { float radius; Vec center; } SphereTemp;
typedef struct { Vec p1, p2, p3; } TriangleTemp;
typedef struct {
	// for any type: type and material
	short type;
	short materialId;
	Vec emission;
	Vec color;

	//union
	//{
	SphereTemp s;
	TriangleTemp t;
	//};

	// for objects that are lights (ie, emit color):
	float area;
} ObjectTemp;

extern void ReInit(const int);
extern void ReInitScene();
extern void UpdateRendering();
extern void UpdateCamera();

extern Camera cameraLeft, cameraRight;
extern Shape *shapes;
extern unsigned int shapeCnt;
extern unsigned int lightCnt;

int amiSmallptCPU;

short width = 640;
short height = 480;
unsigned int *pixelsLeft, *pixelsRight, *pixelsTemp;
char captionBuffer[256];

Material *materials;
SceneTemp scene;

double WallClockTime() 
{
#if defined(__linux__) || defined(__APPLE__) || defined(__ANDROID__)
	struct timeval t;
	gettimeofday(&t, NULL);

	return t.tv_sec*1000 + t.tv_usec / 1000.0;
#elif defined (WIN32)
	return GetTickCount();
#else
	Unsupported Platform !!!
#endif
}

int getMaterialType(unsigned int id, int lineIndex) {
	switch (id) {
	case 0: return LAMBERTIAN; // DIFFUSE
	case 1: return CONDUCTOR; // SPECULAR
	case 2: return DIELECTRIC; // REFRACTION
	case 3: return MATTE;
	case 4: return PLASTIC;
	case 5: return METAL;
	default:
		LOGE("Failed to read material type for sphere #%d: %d\n", lineIndex, id);
		exit(-1);
	}
}

const char* ReadFile(size_t* fileLength, const char* fileName) {
	char *contents;

	FILE *file = fopen(fileName, "rb");
	if (!file) return NULL;

	fseek(file, 0, SEEK_END);
	*fileLength = ftell(file);
	rewind(file);
	contents = (char *)malloc(*fileLength * (sizeof(char)));
	fread(contents, sizeof(char), *fileLength, file);
	fclose(file);

	return contents;
}

bool ReadTxt(char *fileName, bool bvr) {
	int objectCount = 0;
	ObjectTemp *objects;

	LOGI("Reading scene: %s\n", fileName);

	FILE *f = fopen(fileName, "r");
	if (!f) {
		LOGE("Failed to open file: %s\n", fileName);
		exit(-1);
	}

	int c;

	// part 1: gamma correction, sky colors
	c = fscanf(f, "gamma %f  sky1 %f %f %f  sky2 %f %f %f\n", &scene.gammaCorrection,
		&scene.skyColor1.s[0], &scene.skyColor1.s[1], &scene.skyColor1.s[2],
		&scene.skyColor2.s[0], &scene.skyColor2.s[1], &scene.skyColor2.s[2]);
	if (c != 7) {
		scene.gammaCorrection = 2.2f;
		vinit(scene.skyColor1, 0.15f, 0.3f, 0.5f);
		vinit(scene.skyColor2, 1.f, 1.f, 1.f);
	}

	// part 2: camera configuration: camera eye.x eye.y eye.z  target.x target.y target.z
	c = fscanf(f, "camera %f %f %f  %f %f %f\n",
		&cameraLeft.orig.s[0], &cameraLeft.orig.s[1], &cameraLeft.orig.s[2],
		&cameraLeft.target.s[0], &cameraLeft.target.s[1], &cameraLeft.target.s[2]);

	cameraLeft.pitch = 0;
	cameraLeft.yaw = 0;
	cameraLeft.width = width;
	cameraLeft.height = height;

	if (bvr) {
		cameraRight = cameraLeft;

		cameraLeft.orig.s[0] -= DIFF_LEFTRIGHTEYE;
		cameraRight.orig.s[0] += 2 * DIFF_LEFTRIGHTEYE;
	}

	if (c != 6) {
		LOGE("Failed to read 6 camera parameters: %d\n", c);
		exit(-1);
	}

	// part 3: material counts: materials n
	c = fscanf(f, "materials %u\n", &scene.materialCount);
	if (c < 1) {
		LOGE("Failed to read the number of materials: %d\n", c);
		exit(-1);
	}
	// creates space for the 3 default materials
#define NUMBER_OF_DEFAULT_MATERIALS 3
	scene.materialCount += NUMBER_OF_DEFAULT_MATERIALS;

	// part 4: material descriptors
	materials = (Material *)malloc(sizeof(Material) * scene.materialCount);
	materials[0].type = LAMBERTIAN;
	materials[1].type = CONDUCTOR;
	materials[2].type = DIELECTRIC;

	unsigned int lineIndex;
	for (lineIndex = NUMBER_OF_DEFAULT_MATERIALS; lineIndex < scene.materialCount; lineIndex++) {
		Material *material = &materials[lineIndex];
		char materialTypeString[100];

		int readValuesCount = fscanf(f, "%s", materialTypeString);
		if (readValuesCount != 1) {
			LOGE("Failed to read material description type. Expected 1, found %d value(s)\n", readValuesCount);
			exit(-1);
		}

		if (strcmp(materialTypeString, "matte") == 0) {
			material->type = MATTE;
			readValuesCount = fscanf(f, "%f %f %f  %f\n",
				&material->kd.s[0], &material->kd.s[1], &material->kd.s[2],
				&material->sigma);
			if (readValuesCount != 4) {
				LOGE("Failed to read matte material descriptor. Expected 4, found %d value(s)\n", readValuesCount);
				exit(-1);
			}
			clamp(material->sigma, 0.f, 0.9f);
			//            material->sigma = ((float)FLOAT_PI/180.f) * (material->sigma);
			float sigma2 = material->sigma*material->sigma;
			material->A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
			material->B = 0.45f * sigma2 / (sigma2 + 0.09f);

		}
		else if (strcmp(materialTypeString, "plastic") == 0) {
			material->type = PLASTIC;
			readValuesCount = fscanf(f, "%f %f %f  %f %f %f  %f\n",
				&material->kd.s[0], &material->kd.s[1], &material->kd.s[2],
				&material->ks.s[0], &material->ks.s[1], &material->ks.s[2],
				&material->roughness);
			if (readValuesCount != 7) {
				LOGE("Failed to read plastic material descriptor. Expected 4, found %d value(s)\n", readValuesCount);
				exit(-1);
			}
		}
	}

	// part 5: object counts: objects n
	c = fscanf(f, "objects %u\n", &objectCount);
	if (c < 1) {
		LOGE("Failed to read object count: %d\n", c);
		exit(-1);
	}
	LOGE("Scene size: %d\n", objectCount);

	// part 6: object descriptors
	objects = (ObjectTemp *)malloc(sizeof(ObjectTemp) * objectCount);
	unsigned int addedObjects = 0;
	int lightCount = 0;

	for (lineIndex = 0; lineIndex < objectCount; lineIndex++) {
		ObjectTemp *obj = &objects[lineIndex + addedObjects];
		unsigned int mat;
		char objectTypeString[100];

		int readValuesCount = fscanf(f, "%s", objectTypeString);
		if (readValuesCount != 1) {
			LOGE("Failed to read object description type. Expected 1, found %d value(s)\n", readValuesCount);
			exit(-1);
		}
		//        printf("reading an object of type %s and saving on position %d.\n", objectTypeString, lineIndex + addedObjects);

		if (strcmp(objectTypeString, "sphere") == 0) {
			obj->type = SPH;
			// line x: sphere     radius  xyz  emission  color  reflection
			readValuesCount = fscanf(f, "%f  %f %f %f  %f %f %f  %f %f %f  %u\n",
				&obj->s.radius,
				&obj->s.center.s[0], &obj->s.center.s[1], &obj->s.center.s[2],
				&obj->emission.s[0], &obj->emission.s[1], &obj->emission.s[2],
				&obj->color.s[0], &obj->color.s[1], &obj->color.s[2],
				&mat);

			if (readValuesCount != 11) {
				LOGE("Failed to read sphere #%d: %d\n", lineIndex, readValuesCount);
				exit(-1);
			}

			obj->materialId = getMaterialType(mat, lineIndex);

			LOGI("Sphere material id: %d\n", obj->materialId);
			obj->area = 4 * FLOAT_PI * obj->s.radius*obj->s.radius;

			if (!viszero(obj->emission)) {
				lightCnt++;
			}
		}
		else if (strcmp(objectTypeString, "triangle") == 0) {
			obj->type = TRI;
			ObjectTemp *tri = obj;

			readValuesCount = fscanf(f, "%f %f %f  %f %f %f  %f %f %f  %f %f %f  %f %f %f  %d\n",
				&obj->t.p1.s[0], &obj->t.p1.s[1], &obj->t.p1.s[2],
				&obj->t.p2.s[0], &obj->t.p2.s[1], &obj->t.p2.s[2],
				&obj->t.p3.s[0], &obj->t.p3.s[1], &obj->t.p3.s[2],
				&obj->emission.s[0], &obj->emission.s[1], &obj->emission.s[2],
				&obj->color.s[0], &obj->color.s[1], &obj->color.s[2],
				&mat);

			if (readValuesCount != 16) {
				LOGE("Failed to read triangle #%d: %d\n", lineIndex, readValuesCount);
				exit(-1);
			}
			//            LOGI("triangle: %f %f %f   %f %f %f   %f %f %f   %f %f %f   %f %f %f    %d\n",
			//                   obj->p1.x, obj->p1.y, obj->p1.z, obj->p2.x, obj->p2.y, obj->p2.z, obj->p3.x, obj->p3.y, obj->p3.z, obj->emission.x, obj->emission.y, obj->emission.z, obj->color.x, obj->color.y, obj->color.z, obj->refl);
			obj->materialId = getMaterialType(mat, lineIndex);

			// the area of the circle outside the triangle is (abc)/4area
			float a, b, c;

			a = dist(tri->t.p2, tri->t.p1);
			b = dist(tri->t.p3, tri->t.p2);
			c = dist(tri->t.p1, tri->t.p3);

			Vec e1; vsub(e1, tri->t.p2, tri->t.p1);
			Vec e2; vsub(e2, tri->t.p3, tri->t.p1);
			Vec normal; vxcross(normal, e1, e2);
			tri->area = norm(normal) * 0.5f;
			vnorm(normal);

			LOGI("tri #%d normal: %.2f %.2f %.2f\n", lineIndex + addedObjects, normal.s[0], normal.s[1], normal.s[2]);

			if (!viszero(tri->emission)) {
				lightCnt++;
				//                LOGI("Tri emitindo luz: %f %f %f\t%f %f %f\t%f %f %f\n", tri->p1.x, tri->p1.y, tri->p1.z, tri->p2.x, tri->p2.y, tri->p2.z, tri->p3.x, tri->p3.y, tri->p3.z);
				//                LOGI("Color: %f %f %f\tEmission: %f %f %f\n", tri->color.x, tri->color.y, tri->color.z, tri->emission.x, tri->emission.y, tri->emission.z);
				//                LOGI("Area: %f\n", tri->area);
				//                LOGI("Raio: %f\n", tri->radius);
			}
		}
		else if (strcmp(objectTypeString, "model") == 0) {
			char modelFilePath[200];
			Vec modelPosition;
			Vec modelScale;

			obj->type = MOD;
			readValuesCount = fscanf(f, "%s  %f %f %f  %f %f %f  %f %f %f  %f %f %f  %d\n", modelFilePath,
				&modelPosition.s[0], &modelPosition.s[1], &modelPosition.s[2],
				&modelScale.s[0], &modelScale.s[1], &modelScale.s[2],
				&obj->emission.s[0], &obj->emission.s[1], &obj->emission.s[2],
				&obj->color.s[0], &obj->color.s[1], &obj->color.s[2],
				&mat);

			if (readValuesCount != 14) {
				LOGE("Failed to read model #%d: %d\n", lineIndex, readValuesCount);
				exit(-1);
			}

			obj->materialId = getMaterialType(mat, lineIndex);

			tinyobj_attrib_t attrib;
			tinyobj_shape_t* shapes = NULL;
			size_t num_shapes;
			tinyobj_material_t* materials = NULL;
			size_t num_materials;
			size_t data_len = 0;

			const char* data = ReadFile(&data_len, modelFilePath);
			if (data == NULL) {
				LOGE("Error loading obj. Exiting...");
				exit(-1);
			}

			unsigned int flags = (unsigned int)NULL;
			int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &materials,
				&num_materials, data, data_len, flags);
			if (ret != TINYOBJ_SUCCESS) {
				LOGE("Error loading obj. Exiting...");
				exit(1);
			}

			//            LOGI("# of shapes    = %d\n", (int)num_shapes);
			//            LOGI("# of materials = %d\n", (int)num_materials);

			unsigned int numTriangles = attrib.num_face_num_verts;
			LOGI("Started reading model with %d triangles.\n", numTriangles);

			objects = (ObjectTemp *)realloc(objects, sizeof(ObjectTemp) * (objectCount + addedObjects + numTriangles));
			LOGI("Reallocated objects to have %d positions.\n", objectCount + addedObjects + numTriangles);
			obj = &objects[lineIndex + addedObjects];

			unsigned int faceIndex;
			unsigned int faceOffset = 0;

			for (faceIndex = 0; faceIndex < numTriangles; faceIndex++) { //6, 7
				ObjectTemp *tri = &objects[lineIndex + addedObjects + faceIndex + 1];
				//                printf("adding new triangle and saving on position %d.\n", lineIndex + addedObjects + faceIndex  + 1);
				tri->type = TRI;
				tri->emission = obj->emission;
				tri->color = obj->color;
				tri->materialId = obj->materialId;

				int idxVertex1 = attrib.faces[faceOffset + 0].v_idx;
				int idxVertex2 = attrib.faces[faceOffset + 1].v_idx;
				int idxVertex3 = attrib.faces[faceOffset + 2].v_idx;

				Vec p1; vinit(p1, attrib.vertices[3 * idxVertex1 + 0], attrib.vertices[3 * idxVertex1 + 1], attrib.vertices[3 * idxVertex1 + 2]);
				Vec p2; vinit(p2, attrib.vertices[3 * idxVertex2 + 0], attrib.vertices[3 * idxVertex2 + 1], attrib.vertices[3 * idxVertex2 + 2]);
				Vec p3; vinit(p3, attrib.vertices[3 * idxVertex3 + 0], attrib.vertices[3 * idxVertex3 + 1], attrib.vertices[3 * idxVertex3 + 2]);

				// scales the model (according to .scn file)
				vmul(p1, modelScale, p1);
				vmul(p2, modelScale, p2);
				vmul(p3, modelScale, p3);

				// translates the vertices to the model position (according to .scn file)
				vadd(p1, p1, modelPosition);
				vadd(p2, p2, modelPosition);
				vadd(p3, p3, modelPosition);

				// saves the vertices coords into the triangle structure
				tri->t.p1 = p1;
				tri->t.p2 = p2;
				tri->t.p3 = p3;

				// the area of the circle outside the triangle is (abc)/4area
				float a, b, c;

				a = dist(tri->t.p2, tri->t.p1);
				b = dist(tri->t.p3, tri->t.p2);
				c = dist(tri->t.p1, tri->t.p3);

				Vec e1; vsub(e1, tri->t.p2, tri->t.p1);
				Vec e2; vsub(e2, tri->t.p3, tri->t.p1);
				Vec normal; vxcross(normal, e1, e2);

				tri->area = norm(normal) * 0.5f;

				if (!viszero(tri->emission)) {
					lightCnt++;
				}

				faceOffset += 3;
			}

			addedObjects += numTriangles;
		}
	}

	objectCount += addedObjects;

	LOGI("Finished parsing scene. Resulting objectCount = %d\n", objectCount);
	
	fclose(f);

	unsigned int poiCnt = 0, sphCnt = 0, triCnt = 0;

	for (int i = 0; i < objectCount; i++)
	{
		ObjectTemp *obj = &objects[i];

		if (obj->type == SPH) {
			shapeCnt++;
			sphCnt++;
		}
		else if (obj->type == TRI) 
		{
			shapeCnt++;
			triCnt++;
			poiCnt += 3;
		}
	}

	shapes = (Shape *)malloc(sizeof(Shape) * (shapeCnt));
	int curShape = 0;

	for (int i = 0; i < objectCount; i++)
	{
		ObjectTemp *obj = &objects[i];
		
		if (obj->type == SPH)
		{
			shapes[curShape].type = SPHERE;
			shapes[curShape].s.p = obj->s.center;
			shapes[curShape].s.rad = obj->s.radius;
			shapes[curShape].e = obj->emission;
			shapes[curShape].c = obj->color;
			
			if (obj->materialId == LAMBERTIAN || obj->materialId == MATTE) shapes[curShape].refl = DIFF;
			else if (obj->materialId == CONDUCTOR || obj->materialId == METAL || obj->materialId == PLASTIC) shapes[curShape].refl = SPEC;
			else if (obj->materialId == DIELECTRIC) shapes[curShape].refl = REFR;
			else shapes[curShape].refl = DIFF;

			shapes[curShape].area = obj->area;
			curShape++;
		}
		else if (obj->type == TRI)
		{
			shapes[curShape].type = TRIANGLE;
			shapes[curShape].t.p1 = obj->t.p1;// curPoi - 3;
			shapes[curShape].t.p2 = obj->t.p2;// curPoi - 2;
			shapes[curShape].t.p3 = obj->t.p3;// curPoi - 1;
			shapes[curShape].e = obj->emission;
			shapes[curShape].c = obj->color;

            if (obj->materialId == LAMBERTIAN || obj->materialId == MATTE) shapes[curShape].refl = DIFF;
            else if (obj->materialId == CONDUCTOR || obj->materialId == METAL || obj->materialId == PLASTIC) shapes[curShape].refl = SPEC;
            else if (obj->materialId == DIELECTRIC) shapes[curShape].refl = REFR;
            else shapes[curShape].refl = DIFF;

			shapes[curShape].area = obj->area;
			curShape++;
		}
	}

	LOGI("Point count: %u, Sphere count: %u, Triangle count: %u", poiCnt, sphCnt, triCnt);

	return true;
}

bool ReadScene(const char *fileName, bool bvr) {
	LOGI("Reading scene: %s\n", fileName);

	FILE *f = fopen(fileName, "r");
	if (!f) {
		LOGE("Failed to open file: %s, Errorno: %d\n", fileName, 1);//errno);
		return false;
	}

	/* Read the camera position */
	int c = fscanf(f,"camera %f %f %f %f %f %f\n", &cameraLeft.orig.s[0], &cameraLeft.orig.s[1], &cameraLeft.orig.s[2], &cameraLeft.target.s[0], &cameraLeft.target.s[1], &cameraLeft.target.s[2]);
	if (c != 6) {
        fclose(f);
		LOGE("Failed to read 6 camera parameters: %d\n", c);
		return false;
	}

    if (bvr) {
        cameraRight = cameraLeft;

        cameraLeft.orig.s[0] -= DIFF_LEFTRIGHTEYE;
        cameraRight.orig.s[0] += 2 * DIFF_LEFTRIGHTEYE;
    }

	/* Read the shape count */
	c = fscanf(f,"size %u\n", &shapeCnt);
	if (c != 1) {		
        fclose(f);
		LOGE("Failed to read sphere count: %d\n", c);
		return false;
	}

	LOGI("Scene size: %d\n", shapeCnt);

	/* Read all shapes */
	shapes = (Shape *)malloc(sizeof(Shape) * shapeCnt);
	unsigned int i;

	for (i = 0; i < shapeCnt; i++) {
		Shape *s = &shapes[i];
		s->type = SPHERE;

		int mat;
		int c = fscanf(f,"sphere %f  %f %f %f  %f %f %f  %f %f %f  %d\n",
				&s->s.rad,
				&s->s.p.s[0], &s->s.p.s[1], &s->s.p.s[2],
				&s->e.s[0], &s->e.s[1], &s->e.s[2],
				&s->c.s[0], &s->c.s[1], &s->c.s[2],
				&mat);

		switch (mat) {
			case 0:
				s->refl = DIFF;
				break;
			case 1:
				s->refl = SPEC;
				break;
			case 2:
				s->refl = REFR;
				break;
			default:
				LOGE("Failed to read material type for sphere #%d: %d\n", i, mat);
				//return false;
				break;
		}
		if (c != 11) {
			LOGE("Failed to read sphere #%d: %d\n", i, c);
            fclose(f);
			return false;
		}
	}

	fclose(f);
    return true;
}

bool ReadPly(char *fileName, bool bvr) {
	typedef struct Face {
		unsigned char intensity; /* this user attaches intensity to faces */
		unsigned char nverts;    /* number of vertex indices in list */
		int *verts;              /* vertex index list */
	} Face;

	/* information needed to describe the user's data to the PLY routines */
	char *elem_names[] = { /* list of the kinds of elements in the user's object */
		(char *)"vertex", (char *)"face"
	};

	PlyProperty vert_props[] = { /* list of property information for a vertex */
		{ (char *)"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vec,s[0]), 0, 0, 0, 0 },
		{ (char *)"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vec,s[1]), 0, 0, 0, 0 },
		{ (char *)"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vec,s[2]), 0, 0, 0, 0 },
	};

	PlyProperty face_props[] = { /* list of property information for a vertex */
		{ (char *)"intensity", PLY_UCHAR, PLY_UCHAR, offsetof(Face,intensity), 0, 0, 0, 0 },
		{ (char *)"vertex_indices", PLY_INT, PLY_INT, offsetof(Face,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(Face,nverts) },
	};

	int i, j;// , k;
	PlyFile *ply;
	int nelems;
	char **elist;
	int file_type;
	float version;
	int nprops;
	int num_elems;
	PlyProperty **plist;
	//Vec **vlist;
	Face **flist;
	char *elem_name;
	//int num_comments;
	//char **comments;
	//int num_obj_info;
	//char **obj_info;

	cameraLeft.orig.s[0] = 0.0f, cameraLeft.orig.s[1] = 0.0f, cameraLeft.orig.s[2] = 75.0f,
	cameraLeft.target.s[0] = 0.0f, cameraLeft.target.s[1] = 0.0f, cameraLeft.target.s[2] = 0.0f;

    if (bvr) {
        cameraRight = cameraLeft;

        cameraLeft.orig.s[0] -= DIFF_LEFTRIGHTEYE;
        cameraRight.orig.s[0] += 2 * DIFF_LEFTRIGHTEYE;
    }

	/* open a PLY file for reading */
	ply = ply_open_for_reading(fileName, &nelems, &elist, &file_type, &version);

	/* print what we found out about the file */
	LOGI("version %f\n", version);
	LOGI("type %d\n", file_type);

	unsigned int poiCnt = 0;
	Poi* pois;

	/* go through each kind of element that we learned is in the file */
	/* and read them */
	for (i = 0; i < nelems; i++) {
		/* get the description of the first element */
		elem_name = elist[i];
		plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

		/* if we're on vertex elements, read them in */
		if (equal_strings((char *)"vertex", elem_name)) {
			/* create a vertex list to hold all the vertices */
			//vlist = (Vec **)malloc(sizeof(Vec *) * num_elems);
			poiCnt = num_elems;
			pois = (Poi *)malloc(sizeof(Poi) * poiCnt);

			/* set up for getting vertex elements */
			ply_get_property(ply, elem_name, &vert_props[0]);
			ply_get_property(ply, elem_name, &vert_props[1]);
			ply_get_property(ply, elem_name, &vert_props[2]);

			/* grab all the vertex elements */
			for (j = 0; j < num_elems; j++) {
				/* grab and element from the file */
				//vlist[j] = (Vec *)malloc(sizeof(Vec));
				//ply_get_element(ply, (void *)vlist[j]);
				ply_get_element(ply, &pois[j].p);

				/* print out vertex x,y,z for debugging */
				//LOGI("vertex: %g %g %g\n", vlist[j]->x, vlist[j]->y, vlist[j]->z);
			}
		}

		/* if we're on face elements, read them in */
		if (equal_strings((char *)"face", elem_name)) {
			/* create a list to hold all the face elements */
			flist = (Face **)malloc(sizeof(Face *) * num_elems);

			/* set up for getting face elements */
			ply_get_property(ply, elem_name, &face_props[0]);
			ply_get_property(ply, elem_name, &face_props[1]);

			shapeCnt = num_elems;
			
			/* Read all triangles */
			shapes = (Shape *)malloc(sizeof(Shape) * (shapeCnt + MAX_WALLSUN));

			/* grab all the face elements */
			for (j = 0; j < num_elems; j++) {
				/* grab and element from the file */
				flist[j] = (Face *)malloc(sizeof(Face));
				ply_get_element(ply, (void *)flist[j]);

				shapes[j].type = TRIANGLE;
				shapes[j].refl = DIFF;
				shapes[j].t.p1 = pois[flist[j]->verts[0]].p;
				shapes[j].t.p2 = pois[flist[j]->verts[1]].p;
				shapes[j].t.p3 = pois[flist[j]->verts[2]].p;

				free(flist[j]);
			}
			free(flist);
		}
	}

	/* close the PLY file */
	ply_close(ply);

	return true;
}

char* GetFileExt(char * file_name)
{
	int file_name_len = strlen(file_name);
	file_name += file_name_len;

	char *file_ext = NULL;
	for (int i = 0; i <file_name_len; i++)
	{
		if (*file_name == '.')
		{
			file_ext = file_name + 1;
			break;
		}
		file_name--;
	}

	return file_ext;
}

bool Read(char *fileName, bool *walllight, bool bvr) {
	bool ret = true;
	char *fileExt = GetFileExt(fileName);

	if (!strcmp(fileExt, "scn")) ret = ReadScene(fileName, bvr);
	else if (!strcmp(fileExt, "ply")) ret = ReadPly(fileName, bvr);
	else if (!strcmp(fileExt, "txt"))
	{
		ret = ReadTxt(fileName, bvr);
		*walllight = false;
	}

	return ret;
}

void UpdateCamera(bool bvr) {
	float xDirection = sin(cameraLeft.yaw) * cos(cameraLeft.pitch);
	float yDirection = sin(cameraLeft.pitch);
	float zDirection = cos(cameraLeft.yaw) * cos(cameraLeft.pitch);

	Vec directionToCamera = { xDirection, yDirection, zDirection };
	vsmul(cameraLeft.dir, -1, directionToCamera);
	vnorm(cameraLeft.dir);

	const Vec up = { 0.f, 1.f, 0.f };
    const float fov = (FLOAT_PI / 180.f) * HARD_CODED_CAMERA_FOV;

	// axis x is orthogonal to the camera direction and up
	vxcross(cameraLeft.x, cameraLeft.dir, up);
	vnorm(cameraLeft.x);
	// multiplies x axis by aspectRatio * fov
	vsmul(cameraLeft.x, cameraLeft.width * fov / cameraLeft.height, cameraLeft.x);

	// axis y is orthogonal to x and camera direction
	vxcross(cameraLeft.y, cameraLeft.x, cameraLeft.dir);
	vnorm(cameraLeft.y);

	// multiplies y axis by the fov
	vsmul(cameraLeft.y, fov, cameraLeft.y);

    const float invWidth = 1.f / (float)width;
    const float invHeight = 1.f / (float)height;

    float r1 = - 0.5f;
    float r2 = - 0.5f;
    float kcx = r1 * invWidth - .5f;
    float kcy = r2 * invHeight - .5f;

    Vec rdir;
    vinit(rdir,
          cameraLeft.x.s[0] * kcx + cameraLeft.y.s[0] * kcy + cameraLeft.dir.s[0],
          cameraLeft.x.s[1] * kcx + cameraLeft.y.s[1] * kcy + cameraLeft.dir.s[1],
          cameraLeft.x.s[2] * kcx + cameraLeft.y.s[2] * kcy + cameraLeft.dir.s[2]);

    vadd(cameraLeft.start, cameraLeft.orig, rdir);

    r1 = (float) width - 0.5f;
    r2 = (float) height - 0.5f;
    kcx = r1 * invWidth - .5f;
    kcy = r2 * invHeight - .5f;

    vinit(rdir,
          cameraLeft.x.s[0] * kcx + cameraLeft.y.s[0] * kcy + cameraLeft.dir.s[0],
          cameraLeft.x.s[1] * kcx + cameraLeft.y.s[1] * kcy + cameraLeft.dir.s[1],
          cameraLeft.x.s[2] * kcx + cameraLeft.y.s[2] * kcy + cameraLeft.dir.s[2]);

    vadd(cameraLeft.end, cameraLeft.orig, rdir);

	if (bvr)
	{
		xDirection = sin(cameraRight.yaw) * cos(cameraRight.pitch);
		yDirection = sin(cameraRight.pitch);
		zDirection = cos(cameraRight.yaw) * cos(cameraRight.pitch);

		directionToCamera = { xDirection, yDirection, zDirection };
		vsmul(cameraRight.dir, -1, directionToCamera);
		vnorm(cameraRight.dir);

		//up = { 0.f, 1.f, 0.f };
		//fov = (FLOAT_PI / 180.f) * HARD_CODED_CAMERA_FOV;

		// axis x is orthogonal to the camera direction and up
		vxcross(cameraRight.x, cameraRight.dir, up);
		vnorm(cameraRight.x);
		// multiplies x axis by aspectRatio * fov
		vsmul(cameraRight.x, cameraRight.width * fov / cameraRight.height, cameraRight.x);

		// axis y is orthogonal to x and camera direction
		vxcross(cameraRight.y, cameraRight.x, cameraRight.dir);
		vnorm(cameraRight.y);

		// multiplies y axis by the fov
		vsmul(cameraRight.y, fov, cameraRight.y);

        r1 = - 0.5f;
        r2 = - 0.5f;
        kcx = r1 * invWidth - .5f;
        kcy = r2 * invHeight - .5f;

        vinit(rdir,
              cameraRight.x.s[0] * kcx + cameraRight.y.s[0] * kcy + cameraRight.dir.s[0],
              cameraRight.x.s[1] * kcx + cameraRight.y.s[1] * kcy + cameraRight.dir.s[1],
              cameraRight.x.s[2] * kcx + cameraRight.y.s[2] * kcy + cameraRight.dir.s[2]);

        vadd(cameraRight.start, cameraRight.orig, rdir);

        r1 = (float) width - 0.5f;
        r2 = (float) height - 0.5f;
        kcx = r1 * invWidth - .5f;
        kcy = r2 * invHeight - .5f;

        vinit(rdir,
              cameraRight.x.s[0] * kcx + cameraRight.y.s[0] * kcy + cameraRight.dir.s[0],
              cameraRight.x.s[1] * kcx + cameraRight.y.s[1] * kcy + cameraRight.dir.s[1],
              cameraRight.x.s[2] * kcx + cameraRight.y.s[2] * kcy + cameraRight.dir.s[2]);

        vadd(cameraRight.end, cameraRight.orig, rdir);
	}
}

#ifdef __ANDROID__
#define TWO_PI 6.28318530717958647693f
#define PI_OVER_TWO 1.57079632679489661923f

void touchFunc(int deltax, int deltay, bool bvr) {
	if (deltax != 0 || deltay != 0) {
		if (!bvr) {
			// rotate the camera using pitch (nodding movement) and yaw (nonono movement)
			cameraLeft.yaw += deltax * 0.0001;
			cameraLeft.yaw = cameraLeft.yaw - TWO_PI * floor(cameraLeft.yaw / TWO_PI);
			cameraLeft.pitch += -deltay * 0.0001;
			cameraLeft.pitch = clamp(cameraLeft.pitch, -PI_OVER_TWO, PI_OVER_TWO);

			ReInit(0);
		}
	}
}
#endif