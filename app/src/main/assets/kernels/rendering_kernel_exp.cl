
#include "clheader.h"
#include "/sdcard/gamemobile.kmu.ac.kr.vrapp3/include/geom.h"

inline float GetRandom(__global unsigned int *seed0, __global unsigned int *seed1) {
 *seed0 = 36969 * ((*seed0) & 65535) + ((*seed0) >> 16);
 *seed1 = 18000 * ((*seed1) & 65535) + ((*seed1) >> 16);

 unsigned int ires = ((*seed0) << 16) + (*seed1);

 union {
  float f;
  unsigned int ui;
 } res;
 res.ui = (ires & 0x007fffff) | 0x40000000;

 return (res.f - 2.f) / 2.f;
}

inline void swap(float *a, float *b)
{
	float temp = *a;
	*a = *b;
	*b = temp;
}

float SphereIntersect(
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Sphere *s,
 const Ray *r) {
 Vec op;
 vsub(op, s->p, r->o);
 //{ (op).x = (s->p).x - (r->o).x; (op).y = (s->p).y - (r->o).y; (op).z = (s->p).z - (r->o).z; };
 
 float b = vdot(op, r->d); //((op).x * (r->d).x + (op).y * (r->d).y + (op).z * (r->d).z);
 float det = b * b - vdot(op, op) + s->rad * s->rad; //b * b - ((op).x * (op).x + (op).y * (op).y + (op).z * (op).z) + s->rad * s->rad;
 if (det < 0.f)
  return 0.f;
 else
  det = sqrt(det);
 
 float t = b - det;
 if (t > EPSILON)
  return t;
 else { 
  t = b + det;
  
  if (t > EPSILON)
   return t;
  else
   return 0.f;
 } 
}

//Written in the paper "Fast, Minimum Storage Ray/ Triangle Intersection".
float TriangleIntersect(
#ifdef __ANDROID__
__global
#else
__constant
#endif
	const Triangle *tr,
	const Ray *r) { 
	/* returns distance, 0 if nohit */
	Vec v0 = tr->p1, v1 = tr->p2, v2 = tr->p3, e1, e2, tvec, pvec, qvec;
	float t = 0.0f, u, v;
	float det, inv_det;

	vsub(e1, v1, v0);
	//{ e1.x = (v1).x - (v0).x; e1.y = (v1).y - (v0).y; e1.z = (v1).z - (v0).z; }
	vsub(e2, v2, v0);
	//{ e2.x = (v2).x - (v0).x; e2.y = (v2).y - (v0).y; e2.z = (v2).z - (v0).z; }

	vxcross(pvec, r->d, e2);
	//{ pvec.x = (r->d).y * (e2).z - (r->d).z * (e2).y; pvec.y = (r->d).z * (e2).x - (r->d).x * (e2).z; pvec.z = (r->d).x * (e2).y - (r->d).y * (e2).x; }
	det = vdot(e1, pvec);
	//det = ((e1).x * (pvec).x + (e1).y * (pvec).y + (e1).z * (pvec).z);
#ifdef TEST_CULL
	if (det < MTEPSILON) return 0.0f;
	
    	vsub(tvec, r->o, v0);
	//{ tvec.x = (r->o).x - (v0).x; tvec.y = (r->o).y - (v0).y; tvec.z = (r->o).z - (v0).z; }
	u = vdot(tvec, pvec);
	//u = ((tvec).x * (pvec).x + (tvec).y * (pvec).y + (tvec).z * (pvec).z);
	if (u < 0.0 || u > det) return 0.0f;

	vxcross(qvec, tvec, e1);
	//{ qvec.x = (tvec).y * (e1).z - (tvec).z * (e1).y; qvec.y = (tvec).z * (e1).x - (tvec).x * (e1).z; qvec.z = (tvec).x * (e1).y - (tvec).y * (e1).x; }
	v = vdot(r->d, qvec);
	//v = ((r->d).x * (qvec).x + (r->d).y * (qvec).y + (r->d).z * (qvec).z);
	if (v < 0.0 || u + v > det) return 0.0f;

	t = vdot(e2, qvec);
	//t = ((e2).x * (qvec).x + (e2).y * (qvec).y + (e2).z * (qvec).z);
	inv_det = 1.0 / det;

	t *= inv_det;
	u *= inv_det;
	v *= inv_det;
#else
	if (det == 0) return 0.0f;

	inv_det = 1.0 / det;

	vsub(tvec, r->o, v0);
	//{ tvec.x = (r->o).x - (v0).x; tvec.y = (r->o).y - (v0).y; tvec.z = (r->o).z - (v0).z; }
	u = vdot(tvec, pvec) * inv_det;
	//u = ((tvec).x * (pvec).x + (tvec).y * (pvec).y + (tvec).z * (pvec).z) * inv_det;
	if (u < 0.0 || u > 1.0) return 0.0f;

	vxcross(qvec, tvec, e1);
	//{ qvec.x = (tvec).y * (e1).z - (tvec).z * (e1).y; qvec.y = (tvec).z * (e1).x - (tvec).x * (e1).z; qvec.z = (tvec).x * (e1).y - (tvec).y * (e1).x; }
	v = vdot(r->d, qvec) * inv_det;
	//v = ((r->d).x * (qvec).x + (r->d).y * (qvec).y + (r->d).z * (qvec).z) * inv_det;
	if (v < 0.0 || u + v > 1.0) return 0.0f;

	t = vdot(e2, qvec) * inv_det;
	//t = ((e2).x * (qvec).x + (e2).y * (qvec).y + (e2).z * (qvec).z) * inv_det;
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
 //{ (*v).x = xx; (*v).y = yy; (*v).z = zz; };
}
 
/**
 * Test if a ray intersect a bound
 */
bool intersection_bound_test(const Ray r, Bound bound
#ifdef DEBUG_INTERSECTIONS
	, __global float *debug2
#endif
	) {
    float t_min, t_max, t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
    float x_a = 1.0/r.d.x, y_a = 1.0/r.d.y, z_a = 1.0/r.d.z;
    float  x_e = r.o.x, y_e = r.o.y, z_e = r.o.z;

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
	
#ifdef DEBUG_INTERSECTIONS
    if (get_global_id(0) == 0) {
        debug2[0] = t_min, debug2[1] = t_max, debug2[2] = t_xmin, debug2[3] = t_xmax, debug2[4] = t_ymin, debug2[5] = t_ymax, debug2[6] = t_zmin, debug2[7] = t_zmax;
        debug2[8] = x_a, debug2[9] = y_a, debug2[10] = z_a;
        debug2[11] = bound.min_x, debug2[12] = bound.max_x, debug2[13] = bound.min_y, debug2[14] = bound.max_y, debug2[15] = bound.min_z, debug2[16] = bound.max_z;
        debug2[17] = r.o.x, debug2[18] = r.o.y, debug2[19] = r.o.z, debug2[20] = r.d.x, debug2[21] = r.d.y, debug2[22] = r.d.z;
    }
#endif
    return (t_min <= t_max);
}

int IntersectAll(
#ifdef __ANDROID__
    __global
#else
    __constant
#endif
const Shape *shapes,
const short shapeCnt,
const Ray *r,
float *t, unsigned int *id) {
    float inf = (*t) = 1e20f;
    short i = shapeCnt - 1;
    int intersected = 0;

    for (; i--;) {
        float d = 0.0f;
        if (shapes[i].type == SPHERE ) d = SphereIntersect(&shapes[i].s, r);
        if (shapes[i].type == TRIANGLE ) d = TriangleIntersect(&shapes[i].t, r);
        if ((d != 0.0f) && (d < *t)) {
            *t = d;
            *id = i;
            intersected = 1;
        }
    }

    return intersected;
}

int Intersect(
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes,
 const short shapeCnt,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn, 
__constant

 BVHNodeGPU *btl, 
#elif (ACCELSTR == 2)
 __constant
 
 KDNodeGPU *kng, 
 short kngCnt, 
__constant

 int *kn, 
 short knCnt, 
#endif
 const Ray *r,
 float *t, unsigned int *id
#ifdef DEBUG_INTERSECTIONS
 , __global int *debug1,
 __global float *debug2
#endif
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
#ifdef DEBUG_INTERSECTIONS
	debug1[get_global_id(0)] = *id;
	debug2[get_global_id(0)] = *t;
#endif
	return (*t < inf);
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
#ifdef DEBUG_INTERSECTIONS
        atomic_add(&debug1[n], 1);
#endif
 
        if (intersection_bound_test(*r, btn[n].bound
#ifdef DEBUG_INTERSECTIONS
			, debug2
#endif
			)) {
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
				//printf("Unknown node, %d\n", btn[n].leaf);
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
#ifdef DEBUG_INTERSECTIONS
        atomic_add(&debug1[n], 1);
#endif

		if (intersection_bound_test(*r, kng[n].bound
#ifdef DEBUG_INTERSECTIONS
			, debug2
#endif
		)) {
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
				//printf("Unknown node, %d\n", btn[n].leaf);
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
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes,
 const short shapeCnt,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn,
__constant

 BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
__constant

 KDNodeGPU *kng,
 short kngCnt,
__constant

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
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes,
 const short shapeCnt,
 const short lightCnt,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn,
__constant

 BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
__constant

 KDNodeGPU *kng,
 short kngCnt,
__constant

 int *kn,
 short knCnt,
#endif
 __global unsigned int *seed0, __global unsigned int *seed1,
 const Vec *hitPoint,
 const Vec *normal,
 Vec *result) {
 vclr(*result);
 //{ (*result).x = 0.f; (*result).y = 0.f; (*result).z = 0.f; };

 short i;
 short lightsVisited;
	for (i = 0, lightsVisited = 0; i < shapeCnt && lightsVisited < lightCnt; i++) {
#ifdef __ANDROID__
__global
#else
__constant
#endif
  const Shape *light = &shapes[i];

  if (!viszero(light->e)) {
    //if (!(((light->e).x == 0.f) && ((light->e).x == 0.f) && ((light->e).z == 0.f))) {
	// this is a light source (as it has an emission component)...
    lightsVisited++;

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
 
   vsub(shadowRay.d, lightPoint, *hitPoint);
   //{ (shadowRay.d).x = (spherePoint).x - (*hitPoint).x; (shadowRay.d).y = (spherePoint).y - (*hitPoint).y; (shadowRay.d).z = (spherePoint).z - (*hitPoint).z; };
   const float len = sqrt(vdot(shadowRay.d, shadowRay.d));
   //((shadowRay.d).x * (shadowRay.d).x + (shadowRay.d).y * (shadowRay.d).y + (shadowRay.d).z * (shadowRay.d).z));
   vsmul(shadowRay.d, 1.f / len, shadowRay.d);
   //{ float k = (1.f / len); { (shadowRay.d).x = k * (shadowRay.d).x; (shadowRay.d).y = k * (shadowRay.d).y; (shadowRay.d).z = k * (shadowRay.d).z; } };

   float wo = vdot(shadowRay.d, lightNormalWhereShadowRayHit);//((shadowRay.d).x * (unitSpherePoint).x + (shadowRay.d).y * (unitSpherePoint).y + (shadowRay.d).z * (unitSpherePoint).z);
   
   if (wo > 0.f) {
    continue;
   } else
    wo = -wo;
	wo = clamp(wo, 0.f, 1.f);
	
   const float wi = vdot(shadowRay.d, *normal);
   //((shadowRay.d).x * (*normal).x + (shadowRay.d).y * (*normal).y + (shadowRay.d).z * (*normal).z);
   if ((wi > 0.f) && (!IntersectP(shapes, shapeCnt,
#if (ACCELSTR == 1)
        btn, btl,
#elif (ACCELSTR == 2)
        kng, kngCnt, kn, knCnt,
#endif
        &shadowRay, len - EPSILON))) {
    Vec c; vassign(c, light->e);
    //{ (c).x = (light->e).x; (c).y = (light->e).y; (c).z = (light->e).z; };
    const float s = light->area * wi * wo / (len *len);

    vsmul(c, s, c);
	//{ float k = (s); { (c).x = k * (c).x; (c).y = k * (c).y; (c).z = k * (c).z; } };
    vadd(*result, *result, c);
	//{ (*result).x = (*result).x + (c).x; (*result).y = (*result).y + (c).y; (*result).z = (*result).z + (c).z; };
   }
  }
 }
}

void RadianceOnePathTracing(
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes, const short shapeCnt, const short lightCnt,
 const short width, const short height, const short depth,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn, 
__constant

 BVHNodeGPU *btl, 
#elif (ACCELSTR == 2)
__constant 
 
 KDNodeGPU *kng,
 short kngCnt, 
__constant 

 int *kn, 
 short knCnt, 
#endif
 Ray *currentRay,
 __global unsigned int *seed0, __global unsigned int *seed1, 
 __global Vec *throughput, __global char *specularBounce, __global char *terminated,
 __global Result *result, __global FirstHitInfo *fhi
#ifdef DEBUG_INTERSECTIONS
 , __global int *debug1,
 __global float *debug2
#endif
 ) {
  float t;
  unsigned int id = 0;

  if (!Intersect(shapes, shapeCnt,
#if (ACCELSTR == 1)
    btn, btl,
#elif (ACCELSTR == 2)
    kng, kngCnt, kn, knCnt,
#endif
    currentRay, &t, &id
#ifdef DEBUG_INTERSECTIONS
    , debug1, debug2
#endif
  )) {
   *terminated = 1;

   return;
  }

  Vec hitPoint;
  vsmul(hitPoint, t, currentRay->d);
  //{ float k = (t); hitPoint.x = t * (currentRay->d).x; hitPoint.y = t * (currentRay->d).y; hitPoint.z = t * (currentRay->d).z; }
  vadd(hitPoint, currentRay->o, hitPoint);
  //{ hitPoint.x = (currentRay->o).x + (hitPoint).x; hitPoint.y = (currentRay->o).y + (hitPoint).y; hitPoint.z = (currentRay->o).z + (hitPoint).z; } 

  // In the primary ray case, the current color is initialized to the last color...
  if (depth == 0) {
        fhi->x = result->x;
        fhi->y = result->y;
        vassign(fhi->ptFirstHit, hitPoint);
        fhi->idxShape = id;
  }

  Vec normal, col;
  enum Refl refl;
  
  Shape s = shapes[id];

  if (s.type == SPHERE)
  {
	vsub(normal, hitPoint, s.s.p);
	//{ normal.x = (hitPoint).x - (s.s.p).x; normal.y = (hitPoint).y - (s.s.p).y; normal.z = (hitPoint).z - (s.s.p).z; }
  }
  else if (s.type == TRIANGLE)
  {
	Vec v0 = s.t.p1, v1 = s.t.p2, v2 = s.t.p3, e1, e2;

	vsub(e1, v1, v0);
	//{ e1.x = (v1).x - (v0).x; e1.y = (v1).y - (v0).y; e1.z = (v1).z - (v0).z; }
	vsub(e2, v2, v0);
	//{ e2.x = (v2).x - (v0).x; e2.y = (v2).y - (v0).y; e2.z = (v2).z - (v0).z; }

	vxcross(normal, e1, e2);
	//{ normal.x = (e1).y * (e2).z - (e1).z * (e2).y; normal.y = (e1).z * (e2).x - (e1).x * (e2).z; normal.z = (e1).x * (e2).y - (e1).y * (e2).x; }
	//col.x = col.y = 0.0f; 	col.z = 0.6f;
  }

  refl = s.refl;
  col = s.c;

  vnorm(normal);
  //{ float l = 1.f / sqrt(((normal).x * (normal).x + (normal).y * (normal).y + (normal).z * (normal).z)); float k = (l); normal.x = k * (normal).x; normal.y = k * (normal).y; normal.z = k * (normal).z; }
  const float dp = vdot(normal, currentRay->d);
  //const float dp = ((normal).x * (currentRay->d).x + (normal).y * (currentRay->d).y + (normal).z * (currentRay->d).z);

  Vec nl;
  const float invSignDP = -1.f * sign(dp);

  vsmul(nl, invSignDP, normal); 
  //{ float k = (invSignDP); { (nl).x = k * (normal).x; (nl).y = k * (normal).y; (nl).z = k * (normal).z; } };

  Vec eCol;

  eCol = s.e; //vassign(eCol, s.e);
  //{ (eCol).x = (s.e).x; (eCol).y = (s.e).y; (eCol).z = (s.e).z; };

  if (!viszero(eCol)) {
  //if (!(((eCol).x == 0.f) && ((eCol).x == 0.f) && ((eCol).z == 0.f))) {
   if (*specularBounce) {
    vsmul(eCol, fabs(dp), eCol);
    //{ float k = (fabs(dp)); { (eCol).x = k * (eCol).x; (eCol).y = k * (eCol).y; (eCol).z = k * (eCol).z; } };
	vmul(eCol, *throughput, eCol);
    //{ (eCol).x = (*throughput).x * (eCol).x; (eCol).y = (*throughput).y * (eCol).y; (eCol).z = (*throughput).z * (eCol).z; };
	vadd(result->p, result->p, eCol);
    //{ (result)->x = (result)->x + (eCol).x; (result)->y = (result)->y + (eCol).y; (result)->z = (result)->z + (eCol).z; };
   }

   *terminated = 1;
   return;
  }

  if (refl == DIFF) {
   *specularBounce = 0;
   vmul(*throughput, *throughput, col);
   //{ (throughput)->x = (throughput)->x * (col).x; (throughput)->y = (throughput)->y * (col).y; (throughput)->z = (throughput)->z * (col).z; };

   Vec Ld;
   SampleLights(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
    btn, btl,
#elif (ACCELSTR == 2)
    kng, kngCnt, kn, knCnt,
#endif
    seed0, seed1, &hitPoint, &nl, &Ld);

   vmul(Ld, *throughput, Ld);
   //{ (Ld).x = (throughput)->x * (Ld).x; (Ld).y = (throughput)->y * (Ld).y; (Ld).z = (throughput)->z * (Ld).z; };
   vadd(result->p, result->p, Ld);
   //{ (result)->x = (result)->x + (Ld).x; (result)->y = (result)->y + (Ld).y; (result)->z = (result)->z + (Ld).z; };

   float r1 = 2.f * FLOAT_PI * GetRandom(seed0, seed1);
   float r2 = GetRandom(seed0, seed1);
   float r2s = sqrt(r2);

   Vec w;
   vassign(w, nl)
   // { (w).x = (nl).x; (w).y = (nl).y; (w).z = (nl).z; };

   Vec u, a;
   if (fabs(w.x) > .1f) {
    vinit(a, 0.f, 1.f, 0.f);
    //{ (a).x = 0.f; (a).y = 1.f; (a).z = 0.f; };
   } else {
    vinit(a, 1.f, 0.f, 0.f);
    //{ (a).x = 1.f; (a).y = 0.f; (a).z = 0.f; };
   }
   vxcross(u, a, w);
   //{ (u).x = (a).y * (w).z - (a).z * (w).y; (u).y = (a).z * (w).x - (a).x * (w).z; (u).z = (a).x * (w).y - (a).y * (w).x; };
   vnorm(u);
   //{ float l = 1.f / sqrt(((u).x * (u).x + (u).y * (u).y + (u).z * (u).z)); { float k = (l); { (u).x = k * (u).x; (u).y = k * (u).y; (u).z = k * (u).z; } }; };
 
   Vec v, newDir;
   vxcross(v, w, u);
   //{ (v).x = (w).y * (u).z - (w).z * (u).y; (v).y = (w).z * (u).x - (w).x * (u).z; (v).z = (w).x * (u).y - (w).y * (u).x; };
   vsmul(u, cos(r1) * r2s, u);
   //{ float k = (cos(r1) * r2s); { (u).x = k * (u).x; (u).y = k * (u).y; (u).z = k * (u).z; } };
   vsmul(v, sin(r1) * r2s, v);
   //{ float k = (sin(r1) * r2s); { (v).x = k * (v).x; (v).y = k * (v).y; (v).z = k * (v).z; } };
   vadd(newDir, u, v);
   //{ (newDir).x = (u).x + (v).x; (newDir).y = (u).y + (v).y; (newDir).z = (u).z + (v).z; };
   vsmul(w, sqrt(1 - r2), w); //sqrt(max(0.0f, 1 - r2))
   //{ float k = (sqrt(1 - r2)); { (w).x = k * (w).x; (w).y = k * (w).y; (w).z = k * (w).z; } };
   vadd(newDir, newDir, w);
   //{ (newDir).x = (newDir).x + (w).x; (newDir).y = (newDir).y + (w).y; (newDir).z = (newDir).z + (w).z; };

   rinit(*currentRay, hitPoint, newDir);
   
   return;
  } else if (refl == SPEC) {
   *specularBounce = 1;

   Vec newDir;
   vsmul(newDir, 2.f * vdot(normal, currentRay->d), normal);
   //{ float k = (2.f * ((normal).x * (currentRay->d).x + (normal).y * (currentRay->d).y + (normal).z * (currentRay->d).z)); { (newDir).x = k * (normal).x; (newDir).y = k * (normal).y; (newDir).z = k * (normal).z; } };
   vsub(newDir, currentRay->d, newDir);
   //{ (newDir).x = (currentRay->d).x - (newDir).x; (newDir).y = (currentRay->d).y - (newDir).y; (newDir).z = (currentRay->d).z - (newDir).z; };
   vmul(*throughput, *throughput, col);
   //{ (throughput)->x = (throughput)->x * (col).x; (throughput)->y = (throughput)->y * (col).y; (throughput)->z = (throughput)->z * (col).z; };
   rinit(*currentRay, hitPoint, newDir);
   //{ { ((currentRay)->o).x = (hitPoint).x; ((currentRay)->o).y = (hitPoint).y; ((currentRay)->o).z = (hitPoint).z; }; { ((currentRay)->d).x = (newDir).x; ((currentRay)->d).y = (newDir).y; ((currentRay)->d).z = (newDir).z; }; };
   
   return;
  } else {
   *specularBounce = 1;

   Vec newDir;
   vsmul(newDir, 2.f * vdot(normal, currentRay->d), normal);
   //{ float k = (2.f * ((normal).x * (currentRay->d).x + (normal).y * (currentRay->d).y + (normal).z * (currentRay->d).z)); { (newDir).x = k * (normal).x; (newDir).y = k * (normal).y; (newDir).z = k * (normal).z; } };
   vsub(newDir, currentRay->d, newDir);
   //{ (newDir).x = (currentRay->d).x - (newDir).x; (newDir).y = (currentRay->d).y - (newDir).y; (newDir).z = (currentRay->d).z - (newDir).z; };

   Ray reflRay;
   rinit(reflRay, hitPoint, newDir);
   //{ { ((reflRay).o).x = (hitPoint).x; ((reflRay).o).y = (hitPoint).y; ((reflRay).o).z = (hitPoint).z; }; { ((reflRay).d).x = (newDir).x; ((reflRay).d).y = (newDir).y; ((reflRay).d).z = (newDir).z; }; };
   int into = (vdot(normal, nl) > 0);
   //(((normal).x * (nl).x + (normal).y * (nl).y + (normal).z * (nl).z) > 0);

   float nc = 1.f;
   float nt = 1.5f;
   float nnt = into ? nc / nt : nt / nc;
   float ddn = vdot(currentRay->d, nl);
   //((currentRay->d).x * (nl).x + (currentRay->d).y * (nl).y + (currentRay->d).z * (nl).z);
   float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

   if (cos2t < 0.f) {
    vmul(*throughput, *throughput, col);
    //{ (throughput)->x = (throughput)->x * (col).x; (throughput)->y = (throughput)->y * (col).y; (throughput)->z = (throughput)->z * (col).z; };
    rassign(*currentRay, reflRay);
    //{ { ((currentRay)->o).x = ((reflRay).o).x; ((currentRay)->o).y = ((reflRay).o).y; ((currentRay)->o).z = ((reflRay).o).z; }; { ((currentRay)->d).x = ((reflRay).d).x; ((currentRay)->d).y = ((reflRay).d).y; ((currentRay)->d).z = ((reflRay).d).z; }; };
    
	return;
   }

   float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));
   
   Vec nkk, transDir;
   vsmul(nkk, kk, normal);
   //{ float k = (kk); { (nkk).x = k * (normal).x; (nkk).y = k * (normal).y; (nkk).z = k * (normal).z; } };

   vsmul(transDir, nnt, currentRay->d);
   //{ float k = (nnt); { (transDir).x = k * (currentRay->d).x; (transDir).y = k * (currentRay->d).y; (transDir).z = k * (currentRay->d).z; } };
   vsub(transDir, transDir, nkk);
   //{ (transDir).x = (transDir).x - (nkk).x; (transDir).y = (transDir).y - (nkk).y; (transDir).z = (transDir).z - (nkk).z; };
   vnorm(transDir);
   //{ float l = 1.f / sqrt(((transDir).x * (transDir).x + (transDir).y * (transDir).y + (transDir).z * (transDir).z)); { float k = (l); { (transDir).x = k * (transDir).x; (transDir).y = k * (transDir).y; (transDir).z = k * (transDir).z; } }; };

   float a = nt - nc;
   float b = nt + nc;
   float R0 = a * a / (b * b);
   float c = 1 - (into ? -ddn : vdot(transDir, normal));
   //((transDir).x * (normal).x + (transDir).y * (normal).y + (transDir).z * (normal).z));

   float Re = R0 + (1 - R0) * c * c * c * c * c;
   float Tr = 1.f - Re;
   float P = .25f + .5f * Re;
   float RP = Re / P;
   float TP = Tr / (1.f - P);

   if (GetRandom(seed0, seed1) < P) {
    vsmul(*throughput, RP, *throughput);
    //{ float k = (RP); { (throughput)->x = k * (throughput)->x; (throughput)->y = k * (throughput)->y; (throughput)->z = k * (throughput)->z; } };
    vmul(*throughput, *throughput, col);
	//{ (throughput)->x = (throughput)->x * (col).x; (throughput)->y = (throughput)->y * (col).y; (throughput)->z = (throughput)->z * (col).z; };

    rassign(*currentRay, reflRay);
    //{ { ((currentRay)->o).x = ((reflRay).o).x; ((currentRay)->o).y = ((reflRay).o).y; ((currentRay)->o).z = ((reflRay).o).z; }; { ((currentRay)->d).x = ((reflRay).d).x; ((currentRay)->d).y = ((reflRay).d).y; ((currentRay)->d).z = ((reflRay).d).z; }; };
    
	return;
   } else {
    vsmul(*throughput, TP, *throughput);
    //{ float k = (TP); { (throughput)->x = k * (throughput)->x; (throughput)->y = k * (throughput)->y; (throughput)->z = k * (throughput)->z; } };
    vmul(*throughput, *throughput, col);
	//{ (throughput)->x = (throughput)->x * (col).x; (throughput)->y = (throughput)->y * (col).y; (throughput)->z = (throughput)->z * (col).z; };

    rinit(*currentRay, hitPoint, transDir);
    //{ { ((currentRay)->o).x = (hitPoint).x; ((currentRay)->o).y = (hitPoint).y; ((currentRay)->o).z = (hitPoint).z; }; { ((currentRay)->d).x = (transDir).x; ((currentRay)->d).y = (transDir).y; ((currentRay)->d).z = (transDir).z; }; };

	return;
   }
 }
}

__kernel void RadiancePathTracing_exp(
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes,
 const short shapeCnt,
 const short lightCnt,
 const short width, const short height,
 const short curdepth,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn, 
__constant

 BVHNodeGPU *btl, 
#elif (ACCELSTR == 2)
__constant 
 
 KDNodeGPU *kng,
 short kngCnt, 
__constant 

 int *kn, 
 short knCnt, 
#endif
 __global Ray *rays, __global unsigned int *seedsInput, __global Vec *throughput, __global char *specularBounce, __global char *terminated, __global Result *results, __global FirstHitInfo *fhi
#ifdef DEBUG_INTERSECTIONS
 , __global int *debug1,
 __global float *debug2
#endif
 ) {
 const int gid = get_global_id(0);

 const int x = results[gid].x;//gid % width; //
 const int y = results[gid].y;//gid / width; //
 const int sgid = y * width + x;
 const int sgid2 = sgid << 1;

 if (terminated[sgid] != 1)
 {
#if 1
    int max_depth = MAX_DEPTH;
    const int midx = round((float)width / 2.0f), midy = round((float)height / 2.0f);
    const float dist = sqrt((float)(midx - x) * (midx - x) + (float)(midy - y) * (midy - y)), maxdist =  sqrt((float)midx * midx + (float)midy * midy), partdist = maxdist / 3.0f;
    const float partdepth = (float)MAX_DEPTH / 3.0f;

    if (dist >= 0 && dist < partdist) max_depth = MAX_DEPTH;
    else if (dist >= partdist && dist < 2.0f * partdist) max_depth = round(2.0f * partdepth);
    else if (dist >= 2.0f * partdist) max_depth = round(partdepth);

    if (curdepth >= max_depth ) {
        terminated[sgid] = 1;
        return;
    }
#endif
	Ray aray = rays[sgid];
	RadianceOnePathTracing(shapes, shapeCnt, lightCnt, width, height, curdepth,
#if (ACCELSTR == 1)
			btn, btl, 
#elif (ACCELSTR == 2)
			kng, kngCnt, kn, knCnt, 
#endif
			&aray, &seedsInput[sgid2], &seedsInput[sgid2 + 1], &throughput[sgid], &specularBounce[sgid], &terminated[sgid], &results[sgid], &fhi[sgid]
#ifdef DEBUG_INTERSECTIONS
		, debug1, debug2
#endif
	);
	rays[sgid] = aray;
 }
}

void RadianceDirectLighting(
#ifdef __ANDROID__
__global
#else
__constant
#endif
 const Shape *shapes,
 const short shapeCnt,
 const short lightCnt,
#if (ACCELSTR == 1)
__constant

 BVHNodeGPU *btn, 
__constant

 BVHNodeGPU *btl, 
#elif (ACCELSTR == 2)
__constant 
 
 KDNodeGPU *kng,
 short kngCnt,
__constant 

 int *kn, 
 short knCnt,
#endif
 const Ray *startRay,
 __global unsigned int *seed0, __global unsigned int *seed1, 
 Vec *result
#ifdef DEBUG_INTERSECTIONS
 , __global int *debug1,
 __global float *debug2
#endif
 ) {
 Ray currentRay; rassign(currentRay, *startRay);
 //{ { ((currentRay).o).x = ((*startRay).o).x; ((currentRay).o).y = ((*startRay).o).y; ((currentRay).o).z = ((*startRay).o).z; }; { ((currentRay).d).x = ((*startRay).d).x; ((currentRay).d).y = ((*startRay).d).y; ((currentRay).d).z = ((*startRay).d).z; }; };
 Vec rad; vclr(rad);
 //{ (rad).x = 0.f; (rad).y = 0.f; (rad).z = 0.f; };
 Vec throughput; vinit(throughput, 1.f, 1.f, 1.f);
 //{ (throughput).x = 1.f; (throughput).y = 1.f; (throughput).z = 1.f; };

 unsigned int depth = 0;
 char specularBounce = 1;
 
 for (;; ++depth) {
  if (depth > MAX_DEPTH) {
   *result = rad;
   return;
  }

  float t;
  unsigned int id = 0;
  if (!Intersect(shapes, shapeCnt, 
#if (ACCELSTR == 1)
  btn, btl, 
#elif (ACCELSTR == 2)
  kng, kngCnt, kn, knCnt, 
#endif
  &currentRay, &t, &id
#ifdef DEBUG_INTERSECTIONS
  , debug1, debug2
#endif
  )) {
   *result = rad;

   return;
  }
#ifdef __ANDROID__
__global
#else
__constant
#endif
  const Shape *obj = &shapes[id];
  Vec hitPoint;

  vsmul(hitPoint, t, currentRay.d);
  //{ float k = (t); { (hitPoint).x = k * (currentRay.d).x; (hitPoint).y = k * (currentRay.d).y; (hitPoint).z = k * (currentRay.d).z; } };
  vadd(hitPoint, currentRay.o, hitPoint);
  //{ (hitPoint).x = (currentRay.o).x + (hitPoint).x; (hitPoint).y = (currentRay.o).y + (hitPoint).y; (hitPoint).z = (currentRay.o).z + (hitPoint).z; };

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
		//{ normal.x = (e1).y * (e2).z - (e1).z * (e2).y; normal.y = (e1).z * (e2).x - (e1).x * (e2).z; normal.z = (e1).x * (e2).y - (e1).y * (e2).x; }
	
		break;
	}
  }

  vnorm(normal);
  //{ float l = 1.f / sqrt(((normal).x * (normal).x + (normal).y * (normal).y + (normal).z * (normal).z)); { float k = (l); { (normal).x = k * (normal).x; (normal).y = k * (normal).y; (normal).z = k * (normal).z; } }; };

  const float dp = vdot(normal, currentRay.d);//((normal).x * (currentRay.d).x + (normal).y * (currentRay.d).y + (normal).z * (currentRay.d).z);

  Vec nl;

  const float invSignDP = -1.f * sign(dp);
  vsmul(nl, invSignDP, normal);
  //{ float k = (invSignDP); { (nl).x = k * (normal).x; (nl).y = k * (normal).y; (nl).z = k * (normal).z; } };

  Vec eCol;

  vassign(eCol, obj->e);
  //{ (eCol).x = (obj->e).x; (eCol).y = (obj->e).y; (eCol).z = (obj->e).z; };
  if (!viszero(eCol)) {
  //if (!(((eCol).x == 0.f) && ((eCol).x == 0.f) && ((eCol).z == 0.f))) {
   if (specularBounce) {
    vsmul(eCol, fabs(dp), eCol);
    //{ float k = (fabs(dp)); { (eCol).x = k * (eCol).x; (eCol).y = k * (eCol).y; (eCol).z = k * (eCol).z; } };
    vmul(eCol, throughput, eCol);
	//{ (eCol).x = (throughput).x * (eCol).x; (eCol).y = (throughput).y * (eCol).y; (eCol).z = (throughput).z * (eCol).z; };
    vadd(rad, rad, eCol);
	//{ (rad).x = (rad).x + (eCol).x; (rad).y = (rad).y + (eCol).y; (rad).z = (rad).z + (eCol).z; };
   }

   *result = rad;
   return;
  }

  if (obj->refl == DIFF) {
   specularBounce = 0;
   vmul(throughput, throughput, obj->c);
   //{ (throughput).x = (throughput).x * (obj->c).x; (throughput).y = (throughput).y * (obj->c).y; (throughput).z = (throughput).z * (obj->c).z; };

   Vec Ld;
   SampleLights(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
    btn, btl,
#elif (ACCELSTR == 2)
    kng, kngCnt, kn, knCnt,
#endif
    seed0, seed1, &hitPoint, &nl, &Ld);

   vmul(Ld, throughput, Ld);
   //{ (Ld).x = (throughput).x * (Ld).x; (Ld).y = (throughput).y * (Ld).y; (Ld).z = (throughput).z * (Ld).z; };
   vadd(rad, rad, Ld);
   //{ (rad).x = (rad).x + (Ld).x; (rad).y = (rad).y + (Ld).y; (rad).z = (rad).z + (Ld).z; };

   *result = rad;
   return;
  } else if (obj->refl == SPEC) {
   specularBounce = 1;

   Vec newDir;

   vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
   //{ float k = (2.f * ((normal).x * (currentRay.d).x + (normal).y * (currentRay.d).y + (normal).z * (currentRay.d).z)); { (newDir).x = k * (normal).x; (newDir).y = k * (normal).y; (newDir).z = k * (normal).z; } };
   vsub(newDir, currentRay.d, newDir);
   //{ (newDir).x = (currentRay.d).x - (newDir).x; (newDir).y = (currentRay.d).y - (newDir).y; (newDir).z = (currentRay.d).z - (newDir).z; };

   vmul(throughput, throughput, obj->c);
   //{ (throughput).x = (throughput).x * (obj->c).x; (throughput).y = (throughput).y * (obj->c).y; (throughput).z = (throughput).z * (obj->c).z; };
   rinit(currentRay, hitPoint, newDir);
   //{ { ((currentRay).o).x = (hitPoint).x; ((currentRay).o).y = (hitPoint).y; ((currentRay).o).z = (hitPoint).z; }; { ((currentRay).d).x = (newDir).x; ((currentRay).d).y = (newDir).y; ((currentRay).d).z = (newDir).z; }; };
   continue;
  } else {
   specularBounce = 1;

   Vec newDir;
   vsmul(newDir,  2.f * vdot(normal, currentRay.d), normal);
   //{ float k = (2.f * ((normal).x * (currentRay.d).x + (normal).y * (currentRay.d).y + (normal).z * (currentRay.d).z)); { (newDir).x = k * (normal).x; (newDir).y = k * (normal).y; (newDir).z = k * (normal).z; } };
   vsub(newDir, currentRay.d, newDir);
   //{ (newDir).x = (currentRay.d).x - (newDir).x; (newDir).y = (currentRay.d).y - (newDir).y; (newDir).z = (currentRay.d).z - (newDir).z; };

   Ray reflRay; rinit(reflRay, hitPoint, newDir); //{ { ((reflRay).o).x = (hitPoint).x; ((reflRay).o).y = (hitPoint).y; ((reflRay).o).z = (hitPoint).z; }; { ((reflRay).d).x = (newDir).x; ((reflRay).d).y = (newDir).y; ((reflRay).d).z = (newDir).z; }; };
   int into = (vdot(normal, nl) > 0); //(((normal).x * (nl).x + (normal).y * (nl).y + (normal).z * (nl).z) > 0);

   float nc = 1.f;
   float nt = 1.5f;
   float nnt = into ? nc / nt : nt / nc;
   float ddn = vdot(currentRay.d, nl); //((currentRay.d).x * (nl).x + (currentRay.d).y * (nl).y + (currentRay.d).z * (nl).z);
   float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

   if (cos2t < 0.f) {
    vmul(throughput, throughput, obj->c);
    //{ (throughput).x = (throughput).x * (obj->c).x; (throughput).y = (throughput).y * (obj->c).y; (throughput).z = (throughput).z * (obj->c).z; };
    rassign(currentRay, reflRay);
    //{ { ((currentRay).o).x = ((reflRay).o).x; ((currentRay).o).y = ((reflRay).o).y; ((currentRay).o).z = ((reflRay).o).z; }; { ((currentRay).d).x = ((reflRay).d).x; ((currentRay).d).y = ((reflRay).d).y; ((currentRay).d).z = ((reflRay).d).z; }; };
    continue;
   }

   float kk = (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t));

   Vec nkk;
   vsmul(nkk, kk, normal);
   //{ float k = (kk); { (nkk).x = k * (normal).x; (nkk).y = k * (normal).y; (nkk).z = k * (normal).z; } };
   Vec transDir;
   vsmul(transDir, nnt, currentRay.d);
   //{ float k = (nnt); { (transDir).x = k * (currentRay.d).x; (transDir).y = k * (currentRay.d).y; (transDir).z = k * (currentRay.d).z; } };
   vsub(transDir, transDir, nkk);
   //{ (transDir).x = (transDir).x - (nkk).x; (transDir).y = (transDir).y - (nkk).y; (transDir).z = (transDir).z - (nkk).z; };
   vnorm(transDir);
   //{ float l = 1.f / sqrt(((transDir).x * (transDir).x + (transDir).y * (transDir).y + (transDir).z * (transDir).z)); { float k = (l); { (transDir).x = k * (transDir).x; (transDir).y = k * (transDir).y; (transDir).z = k * (transDir).z; } }; };

   float a = nt - nc;
   float b = nt + nc;
   float R0 = a * a / (b * b);
   float c = 1 - (into ? -ddn : vdot(transDir, normal));
   //((transDir).x * (normal).x + (transDir).y * (normal).y + (transDir).z * (normal).z));

   float Re = R0 + (1 - R0) * c * c * c * c * c;
   float Tr = 1.f - Re;
   float P = .25f + .5f * Re;
   float RP = Re / P;
   float TP = Tr / (1.f - P);

   if (GetRandom(seed0, seed1) < P) {
    vsmul(throughput, RP, throughput);
    //{ float k = (RP); { (throughput).x = k * (throughput).x; (throughput).y = k * (throughput).y; (throughput).z = k * (throughput).z; } };
    vmul(throughput, throughput, obj->c);
    //{ (throughput).x = (throughput).x * (obj->c).x; (throughput).y = (throughput).y * (obj->c).y; (throughput).z = (throughput).z * (obj->c).z; };
    rassign(currentRay, reflRay);
    //{ { ((currentRay).o).x = ((reflRay).o).x; ((currentRay).o).y = ((reflRay).o).y; ((currentRay).o).z = ((reflRay).o).z; }; { ((currentRay).d).x = ((reflRay).d).x; ((currentRay).d).y = ((reflRay).d).y; ((currentRay).d).z = ((reflRay).d).z; }; };
    continue;
   } else {
    vsmul(throughput, TP, throughput);
    //{ float k = (TP); { (throughput).x = k * (throughput).x; (throughput).y = k * (throughput).y; (throughput).z = k * (throughput).z; } };
    vmul(throughput, throughput, obj->c);
    //{ (throughput).x = (throughput).x * (obj->c).x; (throughput).y = (throughput).y * (obj->c).y; (throughput).z = (throughput).z * (obj->c).z; };
    rinit(currentRay, hitPoint, transDir);
    //{ { ((currentRay).o).x = (hitPoint).x; ((currentRay).o).y = (hitPoint).y; ((currentRay).o).z = (hitPoint).z; }; { ((currentRay).d).x = (transDir).x; ((currentRay).d).y = (transDir).y; ((currentRay).d).z = (transDir).z; }; };
    continue;
   }
  }
 }
}

__kernel void GenerateCameraRay_exp(
  __constant Camera *camera,
  __global unsigned int *seedsInput,
  const short width, const short height,
  __global Ray *rays, __global Vec *throughput, __global char *specularBounce, __global char *terminated, __global Result *results, __global FirstHitInfo *fhi) {
 const int gid = get_global_id(0);

 const int x = gid % width;
 const int y = gid / width;

 const int sgid = y * width + x;
 const int sgid2 = sgid << 1;

 const float invWidth = 1.f / (float)width;
 const float invHeight = 1.f / (float)height;
 
 const float r1 = GetRandom(&seedsInput[sgid2], &seedsInput[sgid2 + 1]) - .5f;
 const float r2 = GetRandom(&seedsInput[sgid2], &seedsInput[sgid2 + 1]) - .5f;
 const float kcx = (x + r1) * invWidth - .5f;
 const float kcy = (y + r2) * invHeight - .5f;

 vinit(throughput[gid], 1.f, 1.f, 1.f);
 specularBounce[gid] = 1;
 terminated[gid] = 0;
 results[gid].x = x, results[gid].y = y;
 vclr(results[gid].p);

 fhi[gid].x = x, fhi[gid].y = y;
 fhi[gid].idxShape = -1;
 vclr(fhi[gid].ptFirstHit);

 Vec rorig;
 rorig = camera->orig;
 //vsmul(rorig, 0.1f, rdir);
 //{ float k = (0.1f); { (rorig).x = k * (rdir).x; (rorig).y = k * (rdir).y; (rorig).z = k * (rdir).z; } };
 //vadd(rorig, rorig, camera->orig);
 //{ (rorig).x = (rorig).x + (camera->orig).x; (rorig).y = (rorig).y + (camera->orig).y; (rorig).z = (rorig).z + (camera->orig).z; }

 Vec rdir;
 vinit(rdir,
    camera->x.x * kcx + camera->y.x * kcy + camera->dir.x,
 	camera->x.y * kcx + camera->y.y * kcy + camera->dir.y,
 	camera->x.z * kcx + camera->y.z * kcy + camera->dir.z);
 //{ (rdir).x = camera->x.x * kcx + camera->y.x * kcy + camera->dir.x; (rdir).y = camera->x.y * kcx + camera->y.y * kcy + camera->dir.y; (rdir).z = camera->x.z * kcx + camera->y.z * kcy + camera->dir.z; };

 vnorm(rdir);
 //{ float l = 1.f / sqrt(((rdir).x * (rdir).x + (rdir).y * (rdir).y + (rdir).z * (rdir).z)); { float k = (l); { (rdir).x = k * (rdir).x; (rdir).y = k * (rdir).y; (rdir).z = k * (rdir).z; } }; };
 rinit(rays[gid], rorig, rdir);
 //{ { ((*ray).o).x = (rorig).x; ((*ray).o).y = (rorig).y; ((*ray).o).z = (rorig).z; }; { ((*ray).d).x = (rdir).x; ((*ray).d).y = (rdir).y; ((*ray).d).z = (rdir).z; }; };
}

__kernel void FillPixel_exp(
   const short width, const short height, const short bleft,
   __global int* currentSample, __global Vec *colors, __global int *pixels, __global Result *results) {
 const int gid = get_global_id(0);

 const int x = results[gid].x; //gid % width;
 const int y = results[gid].y; //gid / width;

 if (x < 0 || x >= width || y < 0 || y >= height) return;

 const int sgid = y * width + x;
 const int sgid2 = sgid << 1;

 if (currentSample[sgid] == 0) {
  vassign(colors[sgid], results[sgid].p);
  //{ (colors[sgid]).x = (r).x; (colors[sgid]).y = (r).y; (colors[sgid]).z = (r).z; };
 } else {
  const float k1 = currentSample[sgid];
  const float k2 = 1.f / (k1 + 1.f);

  vmad(colors[sgid], colors[sgid], k1, results[sgid].p);
  vsmul(colors[sgid], k2, colors[sgid]);
  //colors[sgid].x = (colors[sgid].x * k1 + results[sgid].p.x) * k2;
  //colors[sgid].y = (colors[sgid].y * k1 + results[sgid].p.y) * k2;
  //colors[sgid].z = (colors[sgid].z * k1 + results[sgid].p.z) * k2;
 }

 currentSample[sgid]++;
 //fhi[sgid].currentColor = colors[sgid];

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

__kernel void FillDiffPixel_exp(
        const short width, const short height, __global Result *results,
        __global FirstHitInfo *fhi, __global const Shape *shapes, const short shapeCnt, __global Camera *cameraDiff, __global ToDiffInfo *tdi,
#if (ACCELSTR == 1)
__constant

     BVHNodeGPU *btn,
    __constant

     BVHNodeGPU *btl
#elif (ACCELSTR == 2)
__constant

     KDNodeGPU *kng,
     short kngCnt,
    __constant

     int *kn,
     short knCnt
#endif
) {
    const int gid = get_global_id(0);

    const int x = results[gid].x; //gid % width;
    const int y = results[gid].y; //gid / width;

    if (x < 0 || x >= width || y < 0 || y >= height) return;

    const int sgid = y * width + x;

    tdi[sgid].x = -1;
    tdi[sgid].y = -1;
    vclr(tdi[sgid].colDiff);

    if (fhi[sgid].idxShape == -1) return;

    const Vec l0 = fhi[sgid].ptFirstHit;

    Vec p0;
    vassign(p0, cameraDiff->start);
    //vadd(p0, cameraDiff->start, cameraDiff->end);
    //vsmul(p0, 0.5f, p0);
    //vnorm(p0);

    Vec n;
    vsmul(n, -1, cameraDiff->dir);
    vnorm(n);

    Vec vl;
    vsub(vl, cameraDiff->orig, fhi[sgid].ptFirstHit);
    vnorm(vl);

    float t1;
    bool brint = PlaneIntersect(n, p0, l0, vl, &t1);

    if (!brint) return;

    float t2;
    Ray ray;
    unsigned int id;

    Vec ve;
    vsub(ve, fhi[sgid].ptFirstHit, cameraDiff->orig);
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

    if (id != fhi[sgid].idxShape) return;
    //if (irint && t1 >= t2) return;

    Vec p;
    vsmul(p, t1, vl);
    vadd(p, p, l0);

    if (p.x >= cameraDiff->start.x && p.x < cameraDiff->end.x && p.y >= cameraDiff->start.y && p.y < cameraDiff->end.y)
    {
        int xDiff = round(((p.x - cameraDiff->start.x) / (cameraDiff->end.x - cameraDiff->start.x)) * (float)(width - 1));
        int yDiff = round(((p.y - cameraDiff->start.y) / (cameraDiff->end.y - cameraDiff->start.y)) * (float)(height - 1));

        const int locPixelDiff = yDiff * width + xDiff;

        tdi[sgid].x = xDiff;
        tdi[sgid].y = yDiff;
        tdi[sgid].colDiff = results[sgid].p;
#if 0
        //lock(the_lock);
          if (currentSampleDiff[locPixelDiff] == 0) {
           vassign(colorsDiff[locPixelDiff], results[sgid].p); //colors[sgid]);
          }
          else {
           const float k1 = currentSampleDiff[locPixelDiff];
           const float k2= 1.f / (k1 + 1.f);

           //Vec redP;
           //vinit(redP, 1.0f, 0.0f, 0.0f);
           vmad(colorsDiff[locPixelDiff], k1, colorsDiff[locPixelDiff], results[sgid].p);
           vsmul(colorsDiff[locPixelDiff], k2, colorsDiff[locPixelDiff]);
          }

          currentSampleDiff[locPixelDiff]++;
          //unlock(the_lock);
#endif
    }
}

#define WINDOW_SIZE 3

void bubbleSort(unsigned int *v,int size)
{
    // bubble-sort
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (v[i] > v[j]) { /* swap? */
                unsigned int tmp = v[i];
                v[i] = v[j];
                v[j] = tmp;
            }
        }
    }
}

unsigned int findMedian(unsigned int *v,int size)
{
    if (size <= 0 ) return -1;
    if (size == 1) return v[0];

    int loc = round((float) size / 2.0f);
    unsigned int median = 0xff000000;

    for (int i = 0; i < loc; i++)
    {
        for (int j = i + 1; j < size; j++)
        {
            if ((v[i] & 0x00ff0000) > (v[j] & 0x00ff0000))
            {
                unsigned int temp = v[i];
                v[i] = v[j];
                v[j] = temp;
            }
        }
    }

    median |= (v[loc] & 0x00ff0000);

    for (int i = 0; i < loc; i++)
    {
        for (int j = i + 1; j < size; j++)
        {
            if ((v[i] & 0x0000ff00) > (v[j] & 0x0000ff00))
            {
                unsigned int temp = v[i];
                v[i] = v[j];
                v[j] = temp;
            }
        }
    }

    median |= (v[loc] & 0x0000ff00);

    for (int i = 0; i < loc; i++)
    {
        for (int j = i + 1; j < size; j++)
        {
            if ((v[i] & 0x000000ff) > (v[j] & 0x000000ff00))
            {
                unsigned int temp = v[i];
                v[i] = v[j];
                v[j] = temp;
            }
        }
    }

    median |= (v[loc] & 0x000000ff);

    return median;
}

// __constant int WINDOW_SIZE=(int)sqrt((float)ARRAY_SIZE_ARG);
__kernel void MedianFilter2D( __global unsigned int *input, __global unsigned int* output, short width, short height)
{
    int filter_offset = WINDOW_SIZE / 2;

    int gid = get_global_id(0);

    const int x = gid % width;
    const int y = gid / width;

    if(y > height || x > width) return;

    unsigned int window[WINDOW_SIZE * WINDOW_SIZE];
    for (int count = 0; count < WINDOW_SIZE * WINDOW_SIZE; count++) {
        window[count] = 0;
    }

    int count = 0;
    for(int k = y - filter_offset; k <= y + filter_offset; k++) {
        for (int l = x - filter_offset; l <= x + filter_offset; l++) {
            if(k >= 0 && l >= 0 && k < height && l < width) window[count++] = input[k * width + l];
        }
    }

    //bubbleSort(window, count + 1);
    unsigned int median = findMedian(window, count + 1);
    output[y * width + x] = median; //window[count / 2];
}

#ifdef CPU_PARTRENDERING
__kernel void RadiancePathTracing_expbox(
#ifdef __ANDROID__
        __global
#else
        __constant
#endif
        const Shape *shapes,
        const short shapeCnt,
        const short lightCnt,
        const short startx, const short starty,
        const short bwidth, const short bheight,
        const short twidth, const short theight,
#if (ACCELSTR == 1)
        __constant

         BVHNodeGPU *btn,
        __constant

         BVHNodeGPU *btl,
#elif (ACCELSTR == 2)
        __constant

         KDNodeGPU *kng,
         short kngCnt,
        __constant

         int *kn,
         short knCnt,
#endif
        __global Ray *rays, __global unsigned int *seedsInput, __global Vec *throughput, __global char *specularBounce, __global char *terminated, __global Result *results
#ifdef DEBUG_INTERSECTIONS
        , __global int *debug1,
         __global float *debug2
#endif
) {
    const int gid = get_global_id(0);
    if (gid >= bwidth * bheight) return ;

    const short x = results[gid].x - startx;//gid % width; //
    const short y = results[gid].y - starty;//gid / width; //

    const int sgid = y * bwidth + x;
    const int sgid2 = sgid << 1;

    if (terminated[sgid] != 1)
    {
        Ray aray = rays[sgid];

        RadianceOnePathTracing(shapes, shapeCnt, lightCnt,
#if (ACCELSTR == 1)
            btn, btl,
#elif (ACCELSTR == 2)
            kng, kngCnt, kn, knCnt,
#endif
            &aray, &seedsInput[sgid2], &seedsInput[sgid2+1], &throughput[sgid], &specularBounce[sgid], &terminated[sgid], &results[sgid].p
#ifdef DEBUG_INTERSECTIONS
        , debug1, debug2
#endif
        );

        rays[sgid] = aray;
    }
}

__kernel void GenerateCameraRay_expbox(
    __constant Camera *camera,
    __global unsigned int *seedsInput,
    const short startx, const short starty,
    const short bwidth, const short bheight,
    const short twidth, const short theight,
    __global Ray *rays, __global Vec *throughput, __global char *specularBounce, __global char *terminated, __global Result *results) {
    const int gid = get_global_id(0);
    if (gid >= bwidth * bheight) return ;

    const short x = (gid % bwidth) + startx;
    const short y = (gid / bwidth) + starty;

    const int sgid = (y - starty) * bwidth + (x - startx);
    const int sgid2 = sgid << 1;

    const float invWidth = 1.f / twidth;
    const float invHeight = 1.f / theight;

    const float r1 = GetRandom(&seedsInput[sgid2], &seedsInput[sgid2 + 1]) - .5f;
    const float r2 = GetRandom(&seedsInput[sgid2], &seedsInput[sgid2 + 1]) - .5f;
    const float kcx = (x + r1) * invWidth - .5f;
    const float kcy = (y + r2) * invHeight - .5f;

    throughput[gid].x = throughput[gid].y = throughput[gid].z = 1.f;
    specularBounce[gid] = 1;
    terminated[gid] = 0;
    results[gid].x = x, results[gid].y = y, results[gid].p.x = results[gid].p.y = results[gid].p.z = 0.f;

    Vec rdir;
    vinit(rdir,
            camera->x.x * kcx + camera->y.x * kcy + camera->dir.x,
            camera->x.y * kcx + camera->y.y * kcy + camera->dir.y,
            camera->x.z * kcx + camera->y.z * kcy + camera->dir.z);
    //{ (rdir).x = camera->x.x * kcx + camera->y.x * kcy + camera->dir.x; (rdir).y = camera->x.y * kcx + camera->y.y * kcy + camera->dir.y; (rdir).z = camera->x.z * kcx + camera->y.z * kcy + camera->dir.z; };

    Vec rorig;
    vsmul(rorig, 0.1f, rdir);
    //{ float k = (0.1f); { (rorig).x = k * (rdir).x; (rorig).y = k * (rdir).y; (rorig).z = k * (rdir).z; } };
    vadd(rorig, rorig, camera->orig);
    //{ (rorig).x = (rorig).x + (camera->orig).x; (rorig).y = (rorig).y + (camera->orig).y; (rorig).z = (rorig).z + (camera->orig).z; }

    vnorm(rdir);
    //{ float l = 1.f / sqrt(((rdir).x * (rdir).x + (rdir).y * (rdir).y + (rdir).z * (rdir).z)); { float k = (l); { (rdir).x = k * (rdir).x; (rdir).y = k * (rdir).y; (rdir).z = k * (rdir).z; } }; };
    rinit(rays[gid], rorig, rdir);
    //{ { ((*ray).o).x = (rorig).x; ((*ray).o).y = (rorig).y; ((*ray).o).z = (rorig).z; }; { ((*ray).d).x = (rdir).x; ((*ray).d).y = (rdir).y; ((*ray).d).z = (rdir).z; }; };
}

__kernel void FillPixel_expbox(
    const short startx, const short starty,
    const short bwidth, const short bheight, const short currentSample,
    __global Vec *colors, __global Result *results, __global int *pixels) {
    const int gid = get_global_id(0);
    if (gid >= bwidth * bheight) return ;
    
    const short x = results[gid].x - startx; //gid % width;
    const short y = results[gid].y - starty; //gid / width;

    const int sgid = y * bwidth + x;
    const int sgid2 = sgid << 1;
    const int locPixel = (bheight - y - 1) * width + x;
    
    if (y >= bheight)
        return;

    if (currentSample == 0) {
        vassign(colors[sgid], results[sgid].p);
        //{ (colors[sgid]).x = (r).x; (colors[sgid]).y = (r).y; (colors[sgid]).z = (r).z; };
    } else {
        const float k1 = currentSample;
        const float k2 = 1.f / (currentSample + 1.f);

        colors[sgid].x = (colors[sgid].x * k1 + results[sgid].p.x) * k2;
        colors[sgid].y = (colors[sgid].y * k1 + results[sgid].p.y) * k2;
        colors[sgid].z = (colors[sgid].z * k1 + results[sgid].p.z) * k2;
    }
#ifdef __ANDROID__
    pixels[locPixel] = (toInt(colors[sgid].x)  << 16) |
        (toInt(colors[sgid].y) << 8) |
        (toInt(colors[sgid].z)) | 0xff000000;
#else
    pixels[locPixel] = (toInt(colors[sgid].x)) |
        (toInt(colors[sgid].y) << 8) |
        (toInt(colors[sgid].z) << 16) | 0xff000000;
#endif
}
#endif