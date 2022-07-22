#ifndef _H_GEOMETRY_3D_
#define _H_GEOMETRY_3D_
#include <jo/jo.h>
#include "VecMath.h"
#include "MatriceMath.h"

typedef struct ColLine {
    vec3 start;
    vec3 end;
} ColLine;

typedef struct Rray {
    vec3 origin;
    vec3 direction;
} Rray;

typedef struct Sphere {
    vec3 position;
    FIXED radius;
} Sphere;

typedef struct AABB {
    vec3 position;
    vec3 size;
} AABB;

typedef struct OBB {
    vec3 position;
    vec3 size; 
    mat3 orientation;
} OBB;

typedef struct Triangle {
	union {
		struct {
			vec3 a;
			vec3 b;
			vec3 c;
		};
		struct {
			vec3 p1;
			vec3 p2;
			vec3 p3;
		};

		vec3 points[3];
		FIXED values[9];
	};
} Triangle;

typedef struct Plane {
vec3 normal;
FIXED distance;
} Plane;

typedef struct Interval {
FIXED min;
FIXED max;
} Interval;

ColLine constructline(vec3 *start, vec3 *end);
Rray constructray(vec3 *origin, vec3 *direction);
Rray FromPoints(vec3 *from, vec3 *to);
Sphere constructsphere(vec3 *position, FIXED radius);
AABB constructaabb(vec3 *position, vec3 *size);
OBB constructobb(vec3 *position, vec3 *size, mat3 *orientation);
Plane constructplane(vec3 *normal, FIXED distance);
Triangle constructtriangle(vec3 *p1, vec3 *p2, vec3 *p3);

vec3 GetMin(AABB *aabb);
vec3 GetMax(AABB *aabb);
AABB FromMinMax(vec3 *min, vec3 *max);

FIXED PlaneEquation(vec3 *pt, Plane *plane);

Interval GetIntervalAABB(AABB *rect, vec3 *axis);
Interval GetIntervalOBB(OBB *rect, vec3 *axis);

bool OverlapOnAxisAABBOBB(AABB *aabb, OBB *obb, vec3 *axis);

bool PointInSphere(vec3 *point, Sphere *sphere);
vec3 ClosestPointSphere(Sphere *sphere, vec3 *point);

bool PointInAABB(vec3 *point, AABB *aabb);
vec3 ClosestPointAABB(AABB *aabb, vec3 *point);

bool PointInOBB(vec3 *point,  OBB *obb);
vec3 ClosestPointOBB(OBB *obb, vec3 *point);

bool PointOnPlane(vec3 *point, Plane *plane);
vec3 ClosestPointPlane(Plane *plane, vec3 *point);

bool PointOnLine(vec3 *point, ColLine *line);
vec3 ClosestPointLine(ColLine *line, vec3 *point);

bool PointOnRay(vec3 *point, Rray *ray);
vec3 ClosestPointRay(Rray *ray, vec3 *point);

bool PointInTriangle(vec3 *p, Triangle *t);
Plane PlaneFromTriangle(Triangle *t);
vec3 ClosestPointTriangle(Triangle *t, vec3 *p, int *out);
bool OverlapOnAxisTriangleAABB(AABB *aabb, Triangle *triangle, vec3 *axis);
bool OverlapOnAxisTriangleOBB(OBB *obb, Triangle *triangle, vec3 *axis);
bool TriangleSphere(Triangle *t, Sphere *s);
bool TriangleAABB(Triangle *t, AABB *a);
bool TriangleOBB(Triangle *t, OBB *o);
bool TrianglePlane(Triangle *t, Plane *p);
bool OverlapOnAxisTriangleTriangle(Triangle *t1, Triangle *t2, vec3 *axis);
bool TriangleTriangle(Triangle *t1, Triangle *t2);
vec3 SatCrossEdge(vec3 *a, vec3 *b, vec3 *c, vec3 *d);
bool TriangleTriangleRobust(Triangle *t1, Triangle *t2);
vec3 Barycentric(vec3 *p, Triangle *t);




bool SphereSphere(Sphere *s1, Sphere *s2);
bool SphereAABB(Sphere *sphere, AABB *aabb);
bool SphereOBB(Sphere *sphere, OBB *obb);
bool SpherePlane(Sphere *sphere, Plane *plane);

bool AABBAABB(AABB *aabb1, AABB *aabb2);
bool AABBOBB(AABB *aabb, OBB *obb);

bool OBBPlane(OBB *obb, Plane *plane);

bool PlanePlane(Plane *plane1, Plane *plane2);

bool RaycastSphere(Sphere *sphere, Rray *ray, FIXED *output);
bool RaycastAABB(AABB *aabb, Rray *ray, FIXED *outter);
bool RaycastOBB(OBB *obb, Rray *ray, FIXED *output);
FIXED RaycastPlane(Plane *plane, Rray *ray);
bool RaycastTriangle(Triangle *triangle, Rray *ray);


bool LinetestSphere(Sphere *sphere, ColLine *line);
bool LinetestAABB(AABB *aabb, ColLine *line);
bool LinetestOBB(OBB *obb, ColLine *line);



#endif