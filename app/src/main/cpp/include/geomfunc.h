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

#ifndef _GEOMFUNC_H
#define	_GEOMFUNC_H

#include <pthread.h>
#include "../../assets/sdcard/include/geom.h"
#include "../../assets/sdcard/include/camera.h"
#include "../../assets/sdcard/include/KDNodeGPU.h"
#include "simplernd.h"

#ifndef SMALLPT_GPU
#define INTERSECT_STACK_SIZE (18)

pthread_mutex_t lock_diff;

inline void swap(float *a, float *b)
{
	float temp = *a;
	*a = *b;
	*b = temp;
}

float SphereIntersect(
#ifdef GPU_KERNEL
__constant
#endif
	const Sphere *s,
	const Ray *r) { /* returns distance, 0 if nohit */
	Vec op; /* Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 */
	vsub(op, s->p, r->o);

	float b = vdot(op, r->d);
	float det = b * b - vdot(op, op) + s->rad * s->rad;
	if (det < 0.f)
		return 0.f;
	else
		det = sqrt(det);

	float t = b - det;
	if (t >  EPSILON)
		return t;
	else {
		t = b + det;

		if (t >  EPSILON)
			return t;
		else
			return 0.f;
	}
}


//Written in the paper "Fast, Minimum Storage Ray/ Triangle Intersection".
float TriangleIntersect(
#ifdef GPU_KERNEL
	__constant
#endif
	const Triangle *tr,
	const Ray *r) { /* returns distance, 0 if nohit */
	Vec v0 = tr->p1, v1 = tr->p2, v2 = tr->p3, e1, e2, tvec, pvec, qvec;
	float t = 0.0f, u, v;
	float det, inv_det;

	vsub(e1, v1, v0);
	vsub(e2, v2, v0);

	vxcross(pvec, r->d, e2);
	det = vdot(e1, pvec);
#ifdef TEST_CULL
	if (det < MTEPSILON) return 0.0f;

	vsub(tvec, r->o, v0);
	u = vdot(tvec, pvec);
	if (u < 0.0 || u > det) return 0.0f;

	vxcross(qvec, tvec, e1);
	v = vdot(r->d, qvec);
	if (v < 0.0 || u + v > det) return 0.0f;

	t = vdot(e2, qvec);
	inv_det = 1.0 / det;

	t *= inv_det;
	u *= inv_det;
	v *= inv_det;
#else
	if (det == 0) return 0.0f;

	inv_det = 1.0 / det;

	vsub(tvec, r->o, v0);
	u = vdot(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0) return 0.0f;

	vxcross(qvec, tvec, e1);
	v = vdot(r->d, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0) return 0.0f;

	t = vdot(e2, qvec) * inv_det;
	if (t > EPSILON) return t;
#endif
	return 0.0f;
}

bool PlaneIntersect(const Vec n, const Vec p0, const Vec l0, const Vec l, float *t)
{
    // assuming vectors are all normalized
    float denom = vdot(n, l);
    if (denom > 0.0f) {
        //Vec p0l0 = p0 - l0;
        Vec p0l0;
        vsub(p0l0, p0, l0); //= p0 - l0;
        *t = vdot(p0l0, n) / denom;

        return (*t >= 0);
    }

    return false;
}

void UniformSampleSphere(const float u1, const float u2, Vec *v) {
	const float zz = 1.f - 2.f * u1;
	const float r = sqrt(max(0.f, 1.f - zz * zz));
	const float phi = 2.f * FLOAT_PI * u2;
	const float xx = r * cos(phi);
	const float yy = r * sin(phi);

	vinit(*v, xx, yy, zz);
}

bool intersection_bound_test(const Ray r, Bound bound) {
    float t_min, t_max, t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
    float x_a = 1.0/r.d.s[0], y_a = 1.0/r.d.s[1], z_a = 1.0/r.d.s[2];
    float  x_e = r.o.s[0], y_e = r.o.s[1], z_e = r.o.s[2];

	// calculate t interval in x-axis
	if (x_a >= 0) {
		t_xmin = (bound.min_x - x_e) * x_a;
		t_xmax = (bound.max_x - x_e) * x_a;
	}
	else {
		t_xmin = (bound.max_x - x_e) * x_a;
		t_xmax = (bound.min_x - x_e) * x_a;
	}

	// calculate t interval in y-axis
	if (y_a >= 0) {
		t_ymin = (bound.min_y - y_e) * y_a;
		t_ymax = (bound.max_y - y_e) * y_a;
	}
	else {
		t_ymin = (bound.max_y - y_e) * y_a;
		t_ymax = (bound.min_y - y_e) * y_a;
	}

	// calculate t interval in z-axis
	if (z_a >= 0) {
		t_zmin = (bound.min_z - z_e) * z_a;
		t_zmax = (bound.max_z - z_e) * z_a;
	}
	else {
		t_zmin = (bound.max_z - z_e) * z_a;
		t_zmax = (bound.min_z - z_e) * z_a;
	}

	// find if there an intersection among three t intervals
	t_min = max(t_xmin, max(t_ymin, t_zmin));
	t_max = min(t_xmax, min(t_ymax, t_zmax));

	return (t_min <= t_max);
}

int Intersect(
#ifdef GPU_KERNEL
__constant
#endif
	const Shape *shapes,
	const short shapeCnt,
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btn,
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
	__constant
#endif
	KDNodeGPU *kng, 
	short kngCnt,
#ifdef GPU_KERNEL
	__constant
#endif
	int *kn, 
	short knCnt,
#endif
#ifdef GPU_KERNEL
	__constant
#endif
	const Ray *r,
	float *t, unsigned int *id
	) {
#if (ACCELSTR == 0)
	int intersected = 0;
	float inf = (*t) = 1e20f;

	short i = shapeCnt - 1;
	for (; i--;) {
		float d = 0.0f;
		if (shapes[i].type == SPHERE ) d = SphereIntersect(&shapes[i].s, r);
		if (shapes[i].type == TRIANGLE ) d = TriangleIntersect(&shapes[i].t, r);
		if ((d != 0.f) && (d < *t)) {
			*t = d;
			*id = i;
			intersected = 1;
		}
	}

	return intersected;
#elif (ACCELSTR == 1)
	(*t) = 1e20f;
	
    // Use static allocation because malloc() can't be called in parallel
    // Use stack to traverse BVH to save space (cost is O(height))
    int stack[INTERSECT_STACK_SIZE];
    int topIndex = INTERSECT_STACK_SIZE;
    stack[--topIndex] = 0; //btn;
    int intersected = 0, status = 0;

    // Do while stack is not empty
    while (topIndex != INTERSECT_STACK_SIZE) {
        int n = stack[topIndex++];
		
        if (intersection_bound_test(*r, btn[n].bound)) {
			if (btn[n].leaf == 0) {				
                stack[--topIndex] = btn[n].nRight;
                stack[--topIndex] = btn[n].nLeft;

                if (topIndex < 0) {
                    //printf("Intersect stack not big enough. Increase INTERSECT_STACK_SIZE!\n");
                    return 0;
                }
			}
			else if (btn[n].leaf == 2) {
                float d = 0.0f;

			    if (shapes[btl[btn[n].nLeft].nShape].type == SPHERE ) d = SphereIntersect(&shapes[btl[btn[n].nLeft].nShape].s, r);
				if (shapes[btl[btn[n].nLeft].nShape].type == TRIANGLE ) d = TriangleIntersect(&shapes[btl[btn[n].nLeft].nShape].t, r);

                if (d != 0.0) {
                    if (d < *t) {
                        *t = d;
                        *id = btl[btn[n].nLeft].nShape;
                    }
                    intersected = 1;
                }
				
			    if (shapes[btl[btn[n].nRight].nShape].type == SPHERE ) d = SphereIntersect(&shapes[btl[btn[n].nRight].nShape].s, r);
				if (shapes[btl[btn[n].nRight].nShape].type == TRIANGLE ) d = TriangleIntersect(&shapes[btl[btn[n].nRight].nShape].t, r);

                if (d != 0.0) {
                    if (d < *t) {
                        *t = d;
                        *id = btl[btn[n].nRight].nShape;
                    }
                    intersected = 1;
                }				
			}
			else {
				//printf("Unknown node, %d\n", tn[n].leaf);
            }
        }
    }
	
    return intersected;
#elif (ACCELSTR == 2)
	(*t) = 1e20f;

	// Use static allocation because malloc() can't be called in parallel
	// Use stack to traverse BVH to save space (cost is O(height))
	int stack[INTERSECT_STACK_SIZE];
	int topIndex = INTERSECT_STACK_SIZE;
	stack[--topIndex] = 1; //tn;
	int intersected = 0, status = 0;

	// Do while stack is not empty
	while (topIndex != INTERSECT_STACK_SIZE) {
		int n = stack[topIndex++];

		if (intersection_bound_test(*r, kng[n].bound)) {
			if (kng[n].leaf == 0) {
				stack[--topIndex] = kng[n].nRight;
				stack[--topIndex] = kng[n].nLeft;

				if (topIndex < 0) {
					//printf("Intersect stack not big enough. Increase INTERSECT_STACK_SIZE!\n");
					return 0;
				}
			}
			else if (kng[n].leaf == 1) {
				float d = 0.0f;

				for (int i = kng[n].min; i < kng[n].max; i++) {
					if (shapes[kn[i]].type == SPHERE) d = SphereIntersect(&shapes[kn[i]].s, r);
					if (shapes[kn[i]].type == TRIANGLE) d = TriangleIntersect(&shapes[kn[i]].t, r);

					if (d != 0.0) {
						if (d < *t) {
							*t = d;
							*id = kn[i];
						}
						intersected = 1;
					}
				}				
			}
			else {
				//printf("Unknown node, %d\n", tn[n].leaf);
			}
		}
	}

	return intersected;
#endif
}

/// finds a random point on a triangle
///  u1, u2 are random [0,1]
///  (u,v) is the barycentric coordinates for the output point
///
static void UniformSampleTriangle(const float u1, const float u2, float *u, float *v) {
	float su1 = (float)sqrt(u1);
	*u = 1.0f - su1;
	*v = u2 * su1;
}

int IntersectP(
#ifdef GPU_KERNEL
__constant
#endif
	const Shape *shapes,
	const short shapeCnt,
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
__constant
#endif
 BVHNodeGPU *btn,
#ifdef GPU_KERNEL
		__constant
#endif
 BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
__constant
#endif
 KDNodeGPU *kng,
 short kngCnt,
#ifdef GPU_KERNEL
__constant
#endif
 int *kn,
 short knCnt,
#endif
	const Ray *r,
	const float maxt) {
#if (ACCELSTR == 0)
    short i = shapeCnt - 1;

    for (; i--;) {
        float d = 0.0f;
        if (shapes[i].type == SPHERE ) d = SphereIntersect(&shapes[i].s, r);
        if (shapes[i].type == TRIANGLE ) d = TriangleIntersect(&shapes[i].t, r);
        if ((d != 0.f) && (d < maxt)) 
            return 1;        
    }

    return 0;
#elif (ACCELSTR == 1)
    // Use static allocation because malloc() can't be called in parallel
    // Use stack to traverse BVH to save space (cost is O(height))
    int stack[INTERSECT_STACK_SIZE];
    int topIndex = INTERSECT_STACK_SIZE;
    stack[--topIndex] = 0; //btn;
    int intersected = 0, status = 0;

    // Do while stack is not empty
    while (topIndex != INTERSECT_STACK_SIZE) {
        int n = stack[topIndex++];

        if (intersection_bound_test(*r, btn[n].bound)) {
            if (btn[n].leaf == 0) {
                stack[--topIndex] = btn[n].nRight;
                stack[--topIndex] = btn[n].nLeft;

                if (topIndex < 0) {
                    //printf("Intersect stack not big enough. Increase INTERSECT_STACK_SIZE!\n");
                    return 0;
                }
			}
			else if (btn[n].leaf == 2) {
                float d = 0.0f;

			    if (shapes[btl[btn[n].nLeft].nShape].type == SPHERE ) d = SphereIntersect(&shapes[btl[btn[n].nLeft].nShape].s, r);
				if (shapes[btl[btn[n].nLeft].nShape].type == TRIANGLE ) d = TriangleIntersect(&shapes[btl[btn[n].nLeft].nShape].t, r);

                if (d != 0.0 && d < maxt) {
                    return 1;
                }

			    if (shapes[btl[btn[n].nRight].nShape].type == SPHERE ) d = SphereIntersect(&shapes[btl[btn[n].nRight].nShape].s, r);
				if (shapes[btl[btn[n].nRight].nShape].type == TRIANGLE ) d = TriangleIntersect(&shapes[btl[btn[n].nRight].nShape].t, r);

                if (d != 0.0 && d < maxt) {
                    return 1;
                }
			}
			else {
				//printf("Unknown node, %d\n", btn[n].leaf);
            }
        }
    }

    return 0;
#elif (ACCELSTR == 2)
	// Use static allocation because malloc() can't be called in parallel
	// Use stack to traverse BVH to save space (cost is O(height))
	int stack[INTERSECT_STACK_SIZE];
	int topIndex = INTERSECT_STACK_SIZE;
	stack[--topIndex] = 1; //tn;
	int intersected = 0, status = 0;

	// Do while stack is not empty
	while (topIndex != INTERSECT_STACK_SIZE) {
		int n = stack[topIndex++];

		if (intersection_bound_test(*r, kng[n].bound)) {
			if (kng[n].leaf == 0) {
				stack[--topIndex] = kng[n].nRight;
				stack[--topIndex] = kng[n].nLeft;

				if (topIndex < 0) {
					//printf("Intersect stack not big enough. Increase INTERSECT_STACK_SIZE!\n");
					return 0;
				}
			}
			else if (kng[n].leaf == 1) {
				float d = 0.0f;

				for (int i = kng[n].min; i < kng[n].max; i++) {
					if (shapes[kn[i]].type == SPHERE) d = SphereIntersect(&shapes[kn[i]].s, r);
					if (shapes[kn[i]].type == TRIANGLE) d = TriangleIntersect(&shapes[kn[i]].t, r);

					if (d != 0.0 && d < maxt) {
                        return 1;
                    }
				}
			}
			else {
				//printf("Unknown node, %d\n", btn[n].leaf);
			}
		}
	}

	return 0;
#endif
}

void SampleLights(
#ifdef GPU_KERNEL
__constant
#endif
	const Shape *shapes,
	const short shapeCnt,
	const short lightCnt, 
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
		__constant
#endif
 BVHNodeGPU *btn,
#ifdef GPU_KERNEL
		__constant
#endif
 BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
__constant
#endif
 KDNodeGPU *kng,
 short kngCnt,
#ifdef GPU_KERNEL
__constant
#endif
 int *kn,
 short knCnt,
#endif
	unsigned int *seed0, unsigned int *seed1,
	const Vec *hitPoint,
	const Vec *normal,
	Vec *result) {
	vclr(*result);

	/* For each light */
	short i;
	short lightsVisited;
	for (i = 0, lightsVisited = 0; i < shapeCnt && lightsVisited < lightCnt; i++) {
#ifdef GPU_KERNEL
__constant
#endif
		const Shape *light = &shapes[i];
		if (!viszero(light->e)) {
			// this is a light source (as it has an emission component)...
            lightsVisited++;

			// the shadow ray starts at the hitPoint and goes to a random point on the light
			Ray shadowRay;
			shadowRay.o = *hitPoint;

			// choose a random point over the area of the light source
			Vec lightPoint;
			Vec lightNormalWhereShadowRayHit;
			switch (light->type) {
            case TRIANGLE:
                {
                    float barycentricU, barycentricV;
					Vec e1, e2;
                    UniformSampleTriangle(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &barycentricU, &barycentricV);
                    vsub(e1, light->t.p2, light->t.p1)
                    vsub(e2, light->t.p3, light->t.p1)
                    vsmul(e1, barycentricU, e1);
                    vsmul(e2, barycentricV, e2);
                    vassign(lightPoint, light->t.p1);
                    vadd(lightPoint, lightPoint, e1);
                    vadd(lightPoint, lightPoint, e2);

                    // sets what is the normal on the triangle where the shadow ray hit it
                    // (it's the same regardless of the point because it's a triangle)
                    vxcross(lightNormalWhereShadowRayHit, e1, e2);
                    vnorm(lightNormalWhereShadowRayHit);
                }
                break;

            case SPHERE:
            default:
                {
                    Vec unitSpherePoint;
                    UniformSampleSphere(GetRandom(seed0, seed1), GetRandom(seed0, seed1), &unitSpherePoint);
                    vsmul(lightPoint, light->s.rad, unitSpherePoint);
                    vadd(lightPoint, lightPoint, light->s.p);

                    // sets what is the normal on the triangle where the shadow ray hit it
                    vassign(lightNormalWhereShadowRayHit, unitSpherePoint);
                }
			}

			/* Build the shadow ray direction */
			vsub(shadowRay.d, lightPoint, *hitPoint);
			const float len = sqrt(vdot(shadowRay.d, shadowRay.d));
			vsmul(shadowRay.d, 1.f / len, shadowRay.d);

			float wo = vdot(shadowRay.d, lightNormalWhereShadowRayHit);
			if (wo > 0.f) {
				/* It is on the other half of the sphere */
				continue;
			} else
				wo = -wo;

			wo = clamp(wo, 0.f, 1.f);
			
			/* Check if the light is visible */
			const float wi = vdot(shadowRay.d, *normal);
			if ((wi > 0.f) && (!IntersectP(shapes, shapeCnt, 
#if (ACCELSTR == 1)
        		btn, btl,
#elif (ACCELSTR == 2)
        		kng, kngCnt, kn, knCnt,
#endif
        		&shadowRay, len - EPSILON))) {
				Vec c; vassign(c, light->e);
				const float s = light->area * wi * wo / (len *len);
				vsmul(c, s, c);

				// last, we add the direct contribution of this light to the output Ld radiance
				vadd(*result, *result, c);
			}
		}
	}
}

bool findDiffPos(
#ifdef __ANDROID__
#ifdef GPU_KERNEL
    __global
#endif
#else
    __constant
#endif
const Shape *shapes, const short shapeCnt,
#if (ACCELSTR == 1)
    __constant

    BVHNodeGPU *btn,
    __constant

    BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
    __constant
#endif
    KDNodeGPU *kng,
    short kngCnt,
#ifdef GPU_KERNEL
    __constant
#endif
    int *kn,
    short knCnt,
#endif
#ifdef GPU_KERNEL
    __global
#endif
    FirstHitInfo *fhi,
#ifdef GPU_KERNEL
    __constant
#endif
    Camera *camera, Vec *ptCol)
{
    const Vec l0 = fhi->ptFirstHit;

    Vec p0;
    vassign(p0, camera->start);
    //vadd(p0, cameraDiff->start, cameraDiff->end);
    //vsmul(p0, 0.5f, p0);
    //vnorm(p0);

    Vec n;
    vsmul(n, -1, camera->dir);
    vnorm(n);

    Vec vl;
    vsub(vl, camera->orig, fhi->ptFirstHit);
    vnorm(vl);

    float t1;
    bool brint = PlaneIntersect(n, p0, l0, vl, &t1);

    if (!brint) return false;

    Vec p;
    vsmad(p, t1, vl, l0);

    if (p.x < camera->start.x || p.x > camera->end.x || p.y < camera->start.y || p.y > camera->end.y || p.z < camera->start.z || p.z > camera->end.z) return false;

    float t2;
    Ray ray;
    unsigned int id;

    Vec ve;
    vsub(ve, fhi->ptFirstHit, camera->orig);
    vnorm(ve);

    rinit(ray, camera->orig, ve);
    //rinit(ray, fhi[sgid].ptFirstHit, vl);

    int irint = Intersect(shapes, shapeCnt,
#if (ACCELSTR == 1)
        btn, btl,
#elif (ACCELSTR == 2)
        kng, kngCnt,
        kn, knCnt,
#endif
        &ray, &t2, &id);

    if (id != fhi->idxShape) return false;
    //if (irint && t1 >= t2) return;

    *ptCol = p;
    return true;
}

bool findSpecRefrPos(
#ifdef __ANDROID__
#ifdef GPU_KERNEL
        __global
#endif
#else
        __constant
#endif
const Shape *shapes, const short shapeCnt,
#if (ACCELSTR == 1)
    __constant

    BVHNodeGPU *btn,
    __constant

    BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
    __constant
#endif
    KDNodeGPU *kng,
    short kngCnt,
#ifdef GPU_KERNEL
    __constant
#endif
    int *kn,
    short knCnt,
#endif
#ifdef GPU_KERNEL
    __global
#endif
    FirstHitInfo *fhi,
#ifdef GPU_KERNEL
    __constant
#endif
    Camera *cameraOrg,
#ifdef GPU_KERNEL
    __constant
#endif
    Camera *cameraDiff, Vec *ptCol)
{
    const Vec l0 = fhi->ptFirstHit;

    Vec p0;
    vassign(p0, cameraDiff->start);
    //vadd(p0, cameraDiff->start, cameraDiff->end);
    //vsmul(p0, 0.5f, p0);
    //vnorm(p0);

    Vec l;
    vsub(l, l0, cameraOrg->orig);

    Vec n;
    vsmul(n, -1, cameraDiff->dir);
    vnorm(n);

    Vec r;
    vsmul(r, 2.0f * vdot(n, l), n);
    vsub(r, r, l);
    vnorm(r);

    float t1;
    bool brint = PlaneIntersect(n, p0, l0, r, &t1);

    if (!brint) return false;

    Vec p;
    vsmad(p, t1, r, l0);

    if (p.x < cameraDiff->start.x || p.x > cameraDiff->end.x || p.y < cameraDiff->start.y || p.y > cameraDiff->end.y || p.z < cameraDiff->start.z || p.z > cameraDiff->end.z) return false;

    float t2;
    Ray ray;
    unsigned int id;

    Vec ve;
    vsub(ve, fhi->ptFirstHit, cameraDiff->orig);
    vnorm(ve);

    rinit(ray, cameraDiff->orig, ve);
    //rinit(ray, fhi[sgid].ptFirstHit, vl);

    int irint = Intersect(shapes, shapeCnt,
#if (ACCELSTR == 1)
            btn, btl,
#elif (ACCELSTR == 2)
        kng, kngCnt,
        kn, knCnt,
#endif
        &ray, &t2, &id);

    if (id != fhi->idxShape) return false;
    //if (irint && t1 >= t2) return;

    *ptCol = p;
    return true;
}

void RadianceOnePathTracing(
#ifdef GPU_KERNEL
	__constant
#endif
	const Shape *shapes, const short shapeCnt, const short lightCnt,
	const short width, const short height, 
#ifdef PAPER_20190701
	const short midwidth, const short midheight, const float maxdist,
#endif
	const short depth,
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btn,
#ifdef GPU_KERNEL
	__constant
#endif	
	BVHNodeGPU *btl,
#elif (ACCELSTR ==2)
#ifdef GPU_KERNEL
	__constant
#endif
	KDNodeGPU *kng, 
	short kngCnt,
#ifdef GPU_KERNEL
	__constant
#endif
	int *kn, 
	short knCnt,
#endif
	Ray *currentRay,
	unsigned int *seed0, unsigned int *seed1, 
	Vec *throughput, char *specularBounce, char *terminated, 
	Result *result, FirstHitInfo *fhi, Camera *cameraOrg, Camera *cameraDiff, ToDiffInfo *tdi, int *currentSampleDiff
 ) {
	float t; /* distance to intersection */
	unsigned int id = 0; /* id of intersected object */

	if (!Intersect(shapes, shapeCnt, 
#if (ACCELSTR == 1)
		btn, btl,
#elif (ACCELSTR == 2)
		kng, kngCnt, kn, knCnt,
#endif
		currentRay, &t, &id)) {
		//*result = *rad;
		*terminated = 1;
		return;
	}

	Vec hitPoint;
	vsmul(hitPoint, t, currentRay->d);
	vadd(hitPoint, currentRay->o, hitPoint);

#ifdef PAPER_20190701
  // In the primary ray case, the current color is initialized to the last color...
  if (depth == 0 && (shapes[id].refl == DIFF || shapes[id].refl == SPEC)) {
    fhi->x = result->x;
    fhi->y = result->y;
    vassign(fhi->ptFirstHit, hitPoint);
    fhi->idxShape = id;

    Vec p;
    bool ret;

    if (shapes[id].refl == DIFF)
        ret = findDiffPos(shapes, shapeCnt,
#if (ACCELSTR == 1)
            btn, btl,
#elif (ACCELSTR == 2)
            kng, kngCnt, kn, knCnt,
#endif
            fhi, cameraDiff, &p);
    else if (shapes[id].refl == SPEC)
        ret = findSpecRefrPos(shapes, shapeCnt,
#if (ACCELSTR == 1)
        btn, btl,
#elif (ACCELSTR == 2)
        kng, kngCnt, kn, knCnt,
#endif
            fhi, cameraOrg, cameraDiff, &p);

    tdi->refl = shapes[id].refl;

    if (ret) // //&& ret
    {
        int xDiff = round(((p.x - cameraDiff->start.x) / (cameraDiff->end.x - cameraDiff->start.x)) * (float)(width - 1));
        int yDiff = round(((p.y - cameraDiff->start.y) / (cameraDiff->end.y - cameraDiff->start.y)) * (float)(height - 1));

        const int locPixelDiff = yDiff * width + xDiff;

        tdi->x = xDiff;
        tdi->y = yDiff;
        tdi->indexDiff = locPixelDiff;

        //atomic_inc(&currentSampleDiff[locPixelDiff]);
		//pthread_mutex_lock(&lock_diff);
        //currentSampleDiff[locPixelDiff]++;
		//pthread_mutex_unlock(&lock_diff);
        //__atomic_fetch_add(&currentSampleDiff[locPixelDiff], 1, __ATOMIC_SEQ_CST);
    }
  }
#endif

	Vec normal, col;
	enum Refl refl;
	Shape s = shapes[id];

	if (s.type == SPHERE)
	{
		vsub(normal, hitPoint, s.s.p);
	}
	else if (s.type == TRIANGLE)
	{
		Vec v0 = s.t.p1, v1 = s.t.p2, v2 = s.t.p3, e1, e2;

		vsub(e1, v1, v0);
		vsub(e2, v2, v0);

		vxcross(normal, e1, e2);
	}
	refl = s.refl;
	col = s.c;
	
	vnorm(normal);

	const float dp = vdot(normal, currentRay->d);

	Vec nl;
	// SIMT optimization
	const float invSignDP = -1.f * sign(dp);

	vsmul(nl, invSignDP, normal);

	/* Add emitted light */
	Vec eCol;

	vassign(eCol, s.e);
	if (!viszero(eCol)) {
		if (*specularBounce) {
			vsmul(eCol, fabs(dp), eCol);
			vmul(eCol, *throughput, eCol);
			vadd(result->p, result->p, eCol);
			//vadd(*rad, *rad, eCol);
		}

		//*result = *rad;
		*terminated = 1;

		return;
	}

	if (refl == DIFF) { /* Ideal DIFFUSE reflection */
		*specularBounce = 0;

		vmul(*throughput, *throughput, col);

		/* Direct lighting component */
		Vec Ld;

		SampleLights(shapes, shapeCnt, lightCnt, 
#if (ACCELSTR == 1)
		    btn, btl,
#elif (ACCELSTR == 2)
		    kng, kngCnt, kn, knCnt,
#endif
	    	seed0, seed1, &hitPoint, &nl, &Ld);
		vmul(Ld, *throughput, Ld);
		vadd(result->p, result->p, Ld);
		//vadd(*rad, *rad, Ld);

		/* Diffuse component */
		float r1 = 2.f * FLOAT_PI * GetRandom(seed0, seed1);
		float r2 = GetRandom(seed0, seed1);
		float r2s = sqrt(r2);

		Vec w; 
		vassign(w, nl);

		Vec u, a;
		if (fabs(w.s[0]) > .1f) {
			vinit(a, 0.f, 1.f, 0.f);
		} else {
			vinit(a, 1.f, 0.f, 0.f);
		}
		vxcross(u, a, w);
		vnorm(u);

		Vec v, newDir;
		vxcross(v, w, u);

		vsmul(u, cos(r1) * r2s, u);
		vsmul(v, sin(r1) * r2s, v);
		vadd(newDir, u, v);
		vsmul(w, sqrt(1 - r2), w);
		vadd(newDir, newDir, w);
		rinit(*currentRay, hitPoint, newDir);
		
		return;
	} else if (refl == SPEC) { /* Ideal SPECULAR reflection */
#ifdef PAPER_20190701
    	tdi->pureDiff = 0;
#endif
		*specularBounce = 1;

		Vec newDir;

		vsmul(newDir, 2.f * vdot(normal, currentRay->d), normal);
		vsub(newDir, currentRay->d, newDir);

		vmul(*throughput, *throughput, col);
		rinit(*currentRay, hitPoint, newDir);

		return;
	} else {
#ifdef PAPER_20190701
    	tdi->pureDiff = 0;
#endif
		*specularBounce = 1;

		Vec newDir;
		vsmul(newDir, 2.f * vdot(normal, currentRay->d), normal);
		vsub(newDir, currentRay->d, newDir);

		Ray reflRay; 
		rinit(reflRay, hitPoint, newDir); /* Ideal dielectric REFRACTION */
		int into = (vdot(normal, nl) > 0); /* Ray from outside going in? */

		float nc = 1.f;
		float nt = 1.5f;
		float nnt = into ? nc / nt : nt / nc;
		float ddn = vdot(currentRay->d, nl);
		float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

		if (cos2t < 0.f) { /* Total internal reflection */
			vmul(*throughput, *throughput, col);
			rassign(*currentRay, reflRay);

			return;
		}

		float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));

		Vec nkk, transDir;
		vsmul(nkk, kk, normal);

		vsmul(transDir, nnt, currentRay->d);
		vsub(transDir, transDir, nkk);
		vnorm(transDir);

		float a = nt - nc;
		float b = nt + nc;
		float R0 = a * a / (b * b);
		float c = 1 - (into ? -ddn : vdot(transDir, normal));

		float Re = R0 + (1 - R0) * c * c * c * c*c;
		float Tr = 1.f - Re;
		float P = .25f + .5f * Re;
		float RP = Re / P;
		float TP = Tr / (1.f - P);

		if (GetRandom(seed0, seed1) < P) { /* R.R. */
			vsmul(*throughput, RP, *throughput);
			vmul(*throughput, *throughput, col);
			rassign(*currentRay, reflRay);
#ifdef PAPER_20190701
    if (depth == 0 && refl == REFR) {
        Vec p;
        bool ret = findSpecRefrPos(shapes, shapeCnt,
        #if (ACCELSTR == 1)
                btn, btl,
        #elif (ACCELSTR == 2)
                kng, kngCnt, kn, knCnt,
        #endif
                fhi, cameraOrg, cameraDiff, &p);

        tdi->refl = refl;

        if (ret) // //&& ret
        {
            int xDiff = round(((p.x - cameraDiff->start.x) / (cameraDiff->end.x - cameraDiff->start.x)) * (float)(width - 1));
            int yDiff = round(((p.y - cameraDiff->start.y) / (cameraDiff->end.y - cameraDiff->start.y)) * (float)(height - 1));

            const int locPixelDiff = yDiff * width + xDiff;

            tdi->x = xDiff;
            tdi->y = yDiff;
            tdi->indexDiff = locPixelDiff;

            //atomic_inc(&currentSampleDiff[locPixelDiff]);
			//pthread_mutex_lock(&lock_diff);
			//currentSampleDiff[locPixelDiff]++;
			//pthread_mutex_unlock(&lock_diff);
        }
    }
#endif
			return;
		} else {
			vsmul(*throughput, TP, *throughput);
			vmul(*throughput, *throughput, col);
			rinit(*currentRay, hitPoint, transDir);
#ifdef PAPER_20190701
    if (depth == 0 && refl == REFR) {
        Vec p;
        bool ret = findDiffPos(shapes, shapeCnt,
        #if (ACCELSTR == 1)
                btn, btl,
        #elif (ACCELSTR == 2)
                kng, kngCnt, kn, knCnt,
        #endif
                fhi, cameraDiff, &p);

        tdi->refl = refl;

        if (ret) // //&& ret
        {
            int xDiff = round(((p.x - cameraDiff->start.x) / (cameraDiff->end.x - cameraDiff->start.x)) * (float)(width - 1));
            int yDiff = round(((p.y - cameraDiff->start.y) / (cameraDiff->end.y - cameraDiff->start.y)) * (float)(height - 1));

            const int locPixelDiff = yDiff * width + xDiff;

            tdi->x = xDiff;
            tdi->y = yDiff;
            tdi->indexDiff = locPixelDiff;

            //atomic_inc(&currentSampleDiff[locPixelDiff]);
			//pthread_mutex_lock(&lock_diff);
			//currentSampleDiff[locPixelDiff]++;
			//pthread_mutex_unlock(&lock_diff);
        }
    }
#endif
			return;
		}
	}
}

void RadiancePathTracing(
#ifdef GPU_KERNEL
	__constant
#endif
	const int gid,
	const Shape *shapes,
	const short shapeCnt,
	const short lightCnt, 
	const short width, const short height, 
#ifdef PAPER_20190701
	const short midwidth, const short midheight, const float mindist, const float maxdist,
#endif
	const short curdepth,
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btn,
#ifdef GPU_KERNEL
	__constant
#endif	
	BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
	__constant
#endif
	KDNodeGPU *kng, 
	short kngCnt,
#ifdef GPU_KERNEL
	__constant
#endif
	int *kn, 
	short knCnt,
#endif
	Ray *rays, unsigned int *seedsInput, Vec *throughput, char *specularBounce, char *terminated, Result *results, FirstHitInfo *fhi, ToDiffInfo *tdi, Camera *cameraOrg, Camera *cameraDiff, int *currentSampleDiff
 ) {
 const int x = results[gid].x;//gid % width; //
 const int y = results[gid].y;//gid / width; //
 const int sgid = y * width + x; //(height - y - 1) * width + x;
 const int sgid2 = sgid << 1;
 
 if (terminated[sgid] != 1)
 {
	Ray aray = rays[sgid];
	RadianceOnePathTracing(shapes, shapeCnt, lightCnt, width, height, 
#ifdef PAPER_20190701
            midwidth, midheight, maxdist,
#endif
			curdepth,
#if (ACCELSTR == 1)
			btn, btl, 
#elif (ACCELSTR == 2)
			kng, kngCnt, kn, knCnt, 
#endif
			&aray, &seedsInput[sgid2], &seedsInput[sgid2 + 1], &throughput[sgid], &specularBounce[sgid], &terminated[sgid], &results[sgid], &fhi[sgid], cameraOrg, cameraDiff, &tdi[sgid], currentSampleDiff
#ifdef DEBUG_INTERSECTIONS
		, debug1, debug2
#endif
	);
	rays[sgid] = aray;

    if (tdi[sgid].x > -1 && tdi[sgid].y > -1 && tdi[sgid].indexDiff > -1) { // && tdi[sgid].refl == DIFF
        vassign(tdi[sgid].colDiff, results[sgid].p);
        //vinit(tdi[sgid].colDiff, 1.0f, 0.0f, 0.0f);
    }

#ifdef PAPER_20190701
    if (curdepth < 1 || tdi[sgid].refl == REFR) return;
    if (tdi[sgid].pureDiff == 0) return;

    cl_float2 p0, p1;
    p0.x = x, p0.y = y;
    p1.x = midwidth, p1.y = midheight;
    const float d = dist2(p1, p0); //sqrt((float)(midwidth - x) * (midwidth - x) + (float)(midheight - y) * (midheight - y));

    if (d < mindist) return;
    const float prob = clamp(d / maxdist, 0.0f, 0.9f);

    float rand = GetRandom(&seedsInput[sgid2], &seedsInput[sgid2 + 1]);
    if (rand < prob) {
        results[gid].depth_stopped = 1; //curdepth;
        terminated[sgid] = 1;
        return;
    }
#endif
 }

}

void RadianceDirectLighting(
#ifdef GPU_KERNEL
	__constant
#endif
	const Shape *shapes,
	const short shapeCnt,
	const short lightCnt, 
#if (ACCELSTR == 1)
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btn,
#ifdef GPU_KERNEL
	__constant
#endif
	BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
#ifdef GPU_KERNEL
	__constant
#endif
	KDNodeGPU *kng, 
	short kngCnt,
#ifdef GPU_KERNEL
	__constant
#endif
	int *kn, 
	short knCnt,
#endif
#ifdef GPU_KERNEL
	__constant
#endif
	const Ray *startRay,
	unsigned int *seed0, unsigned int *seed1,
	Vec *result) {
	Ray currentRay; rassign(currentRay, *startRay);
	Vec rad; vinit(rad, 0.f, 0.f, 0.f);
	Vec throughput; vinit(throughput, 1.f, 1.f, 1.f);

	unsigned int depth = 0;
	char specularBounce = 1;
	for (;; ++depth) {
		// Removed Russian Roulette in order to improve execution on SIMT
		if (depth > MAX_DEPTH) {
			*result = rad;
			return;
		}

		float t; /* distance to intersection */
		unsigned int id = 0; /* id of intersected object */
		if (!Intersect(shapes, shapeCnt, 
#if (ACCELSTR == 1)
			btn, btl,
#elif (ACCELSTR == 2)
			kng, kngCnt, kn, knCnt,
#endif
			&currentRay, &t, &id)) {
			*result = rad; /* if miss, return */
			return;
		}

#ifdef GPU_KERNEL
__constant
#endif
		const Shape *obj = &shapes[id]; /* the hit object */

		Vec hitPoint;
		vsmul(hitPoint, t, currentRay.d);
		vadd(hitPoint, currentRay.o, hitPoint);

		// the normal at the intersection point
		Vec normal;
		switch (obj->type) {
		case SPHERE:
			vsub(normal, hitPoint, obj->s.p);			
			break;

		case TRIANGLE:
			{
			Vec e1; 
			vsub(e1, obj->t.p2, obj->t.p1);
			Vec e2; 
			vsub(e2, obj->t.p3, obj->t.p1);
			
			vxcross(normal, e1, e2);
			break;
			}
		}
		vnorm(normal);

		const float dp = vdot(normal, currentRay.d);

		Vec nl;
		// SIMT optimization
		const float invSignDP = -1.f * sign(dp);
		vsmul(nl, invSignDP, normal);

		/* Add emitted light */
		Vec eCol; 
		
		vassign(eCol, obj->e);
		if (!viszero(eCol)) {
			if (specularBounce) {
				vsmul(eCol, fabs(dp), eCol);
				vmul(eCol, throughput, eCol);
				vadd(rad, rad, eCol);
			}

			*result = rad;
			return;
		}

		if (obj->refl == DIFF) { /* Ideal DIFFUSE reflection */
			specularBounce = 0;
			vmul(throughput, throughput, obj->c);

			/* Direct lighting component */

			Vec Ld;
			SampleLights(shapes, shapeCnt, lightCnt, 
			#if (ACCELSTR == 1)
    			btn, btl,
#elif (ACCELSTR == 2)
    			kng, kngCnt, kn, knCnt,
#endif
    			seed0, seed1, &hitPoint, &nl, &Ld);
			vmul(Ld, throughput, Ld);
			vadd(rad, rad, Ld);

			*result = rad;
			return;
		} else if (obj->refl == SPEC) { /* Ideal SPECULAR reflection */
			specularBounce = 1;

			Vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			vmul(throughput, throughput, obj->c);

			rinit(currentRay, hitPoint, newDir);
			continue;
		} else {
			specularBounce = 1;

			Vec newDir;
			vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
			vsub(newDir, currentRay.d, newDir);

			Ray reflRay; rinit(reflRay, hitPoint, newDir); /* Ideal dielectric REFRACTION */
			int into = (vdot(normal, nl) > 0); /* Ray from outside going in? */

			float nc = 1.f;
			float nt = 1.5f;
			float nnt = into ? nc / nt : nt / nc;
			float ddn = vdot(currentRay.d, nl);
			float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

			if (cos2t < 0.f)  { /* Total internal reflection */
				vmul(throughput, throughput, obj->c);

				rassign(currentRay, reflRay);
				continue;
			}

			float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));
			Vec nkk;
			vsmul(nkk, kk, normal);
			Vec transDir;
			vsmul(transDir, nnt, currentRay.d);
			vsub(transDir, transDir, nkk);
			vnorm(transDir);

			float a = nt - nc;
			float b = nt + nc;
			float R0 = a * a / (b * b);
			float c = 1 - (into ? -ddn : vdot(transDir, normal));

			float Re = R0 + (1 - R0) * c * c * c * c*c;
			float Tr = 1.f - Re;
			float P = .25f + .5f * Re;
			float RP = Re / P;
			float TP = Tr / (1.f - P);

			if (GetRandom(seed0, seed1) < P) { /* R.R. */
				vsmul(throughput, RP, throughput);
				vmul(throughput, throughput, obj->c);

				rassign(currentRay, reflRay);
				continue;
			} else {
				vsmul(throughput, TP, throughput);
				vmul(throughput, throughput, obj->c);

				rinit(currentRay, hitPoint, transDir);
				continue;
			}
		}
	}
}
#endif

void FillPixel(
   const short gid, const short width, const short height, const short bleft, Result *results, 
#if 0
   ToDiffInfo *tdi,
#endif
   int* currentSample, int* currentSampleDiff, Vec *colors, Vec *colorsDiff, int *pixels) {
 const int x = results[gid].x; //gid % width;
 const int y = results[gid].y; //gid / width;

 if (x < 0 || x >= width || y < 0 || y >= height) return;

 const int sgid = y * width + x;
 const int sgid2 = sgid << 1;

 if (currentSample[sgid] == 0) {
  vassign(colors[sgid], results[sgid].p);
  //{ (colors[sgid]).x = (r).x; (colors[sgid]).y = (r).y; (colors[sgid]).z = (r).z; };
 } else {
  const float k1 = currentSample[sgid] ;
  const float k2 = 1.f / (k1 + 1.f);

  vmad(colors[sgid], k1, colors[sgid], results[sgid].p);
  vsmul(colors[sgid], k2, colors[sgid]);
  //colors[sgid].x = (colors[sgid].x * k1 + results[sgid].p.x) * k2;
  //colors[sgid].y = (colors[sgid].y * k1 + results[sgid].p.y) * k2;
  //colors[sgid].z = (colors[sgid].z * k1 + results[sgid].p.z) * k2;
 }

 currentSample[sgid]++;

#ifdef __ANDROID__
 const int locPixelInv = (height - y - 1) * width + x;

 pixels[locPixelInv] = (toInt(colors[sgid].x)  << 16) |
     (toInt(colors[sgid].y) << 8) |
     (toInt(colors[sgid].z)) | 0xff000000;
#else
 pixelsLeft[locPixelLeft] = (toInt(colors[sgid].x)) |
     (toInt(colors[sgid].y) << 8) |
     (toInt(colors[sgid].z) << 16) | 0xff000000;

 if (bleft) {
     pixelsRight[locPixelRight] = (toInt(colors[sgid].x)) |
         (toInt(colors[sgid].y) << 8) |
         (toInt(colors[sgid].z) << 16) | 0xff000000;
 }
#endif
}
#endif	/* _GEOMFUNC_H */
