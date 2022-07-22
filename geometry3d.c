#include <jo/jo.h>
#include "geometry3d.h"
#include "VecMath.h"
#include "MatriceMath.h"

ColLine constructline(vec3 *start, vec3 *end) {
    ColLine retline; 
    retline.start = *start;
    retline.end = *end;
    return retline;
}

FIXED lengthline3d(ColLine *line) {
    vec3 retmag = vec3sub(&line->start, &line->end);
    return magvec3(&retmag);
}

FIXED lengthline3dsq(ColLine *line) {
    vec3 retmag = vec3sub(&line->start, &line->end);
    return magsqvec3(&retmag);
}

Rray constructray(vec3 *origin, vec3 *direction) {
    Rray retray;
    retray.direction = *direction;
    retray.origin = *origin;
    return retray;
}

AABB constructaabb(vec3 *position, vec3 *size) {
    AABB retaabb;
    retaabb.position = *position;
    retaabb.size = *size;
    return retaabb;
}

OBB constructobb(vec3 *position, vec3 *size, mat3 *orientation) {
    OBB retobb;
    retobb.position = *position;
    retobb.size = *size;
    retobb.orientation = *orientation;
    return retobb;
}

Plane constructplane(vec3 *normal, FIXED distance) {
    Plane retplane;
    retplane.distance = distance;
    retplane.normal = *normal;
    return retplane;
}

Triangle constructtriangle(vec3 *p1, vec3 *p2, vec3 *p3) {
    Triangle rettri;
    rettri.a = *p1;
    rettri.b = *p2;
    rettri.c = *p3;
    return rettri;
}

Sphere constructsphere(vec3 *position, FIXED radius) {
    Sphere retsphere;
    retsphere.position = *position;
    retsphere.radius = radius;
    return retsphere;
}

Rray FromPoints(vec3 *from, vec3 *to) {
vec3 tbn = vec3sub(to, from);
vec3 dnb = normalizedvec3(&tbn);
return constructray(from, &dnb);
}

vec3 GetMin(AABB *aabb) { //returns the min point of an aabb
vec3 p1 = vec3add(&aabb->position, &aabb->size);
vec3 p2 = vec3sub(&aabb->position, &aabb->size);
vec3 ret = {MIN(p1.x, p2.x), MIN(p1.y, p2.y), MIN(p1.z, p2.z)};
return ret;
}

vec3 GetMax(AABB *aabb) { //returns the max point of an aabb
vec3 p1 = vec3add(&aabb->position, &aabb->size);
vec3 p2 = vec3sub(&aabb->position, &aabb->size);
vec3 ret = {MAX(p1.x, p2.x), MAX(p1.y, p2.y), MAX(p1.z, p2.z)};
return ret;
}

AABB FromMinMax(vec3 *min, vec3 *max) { //creates a bounding box from a minimum and maximum vec3 point
vec3 tbf = vec3add(min, max); 
tbf = vec3flt(&tbf, JO_FIXED_1_DIV_2);
vec3 tbf2 = vec3sub(max, min); 
tbf = vec3flt(&tbf2, JO_FIXED_1_DIV_2);
return constructaabb(&tbf,&tbf2);
}

FIXED PlaneEquation(vec3 *pt, Plane *plane) { //equation of the plane
    return dotvec3(pt, &plane->normal) - plane->distance;
}

bool PointInSphere(vec3 *point, Sphere *sphere) { //check if a 3d point is inside a sphere
    vec3 tbs = vec3sub(point, &sphere->position);
    FIXED magSq = magsqvec3(&tbs);
    FIXED radSq = sphere->radius;

    return magSq < radSq;
}

vec3 ClosestPointSphere(Sphere *sphere, vec3 *point) { //get the closest point on a sphere.
    vec3 sphereToPoint = vec3sub(point, &sphere->position);
    normalizevec3(&sphereToPoint);
    sphereToPoint = vec3flt(&sphereToPoint, sphere->radius);
    return vec3add(&sphereToPoint, &sphere->position);
}

bool PointInAABB(vec3 *point, AABB *aabb) { //check if a 3d point is inside an axis aligned bounding box
	vec3 min = GetMin(aabb);
	vec3 max = GetMax(aabb);

	if (point->x < min.x || point->y < min.y || point->z < min.z) {
		return false;
	}
	if (point->x > max.x || point->y > max.y || point->z > max.z) {
		return false;
	}

	return true;
}

vec3 ClosestPointAABB(AABB *aabb, vec3 *point) { //get the closest point on an axis aligned bounding box
    vec3 result = *point;
    vec3 min = GetMin(aabb);
    vec3 max = GetMax(aabb);

    result.x = (result.x < min.x) ? min.x : result.x;
    result.y = (result.y < min.x) ? min.y : result.y;
    result.z = (result.z < min.x) ? min.z : result.z;
    result.x = (result.x > max.x) ? max.x : result.x;
    result.y = (result.y > max.x) ? max.y : result.y;
    result.z = (result.z > max.x) ? max.z : result.z;

    return result;
}

bool PointInOBB(vec3 *point,  OBB *obb) {
    vec3 dir = vec3sub(point, &obb->position);

    for (int i = 0; i < 3; ++i) {
        const FIXED* orientation = &obb->orientation.asArray[i * 3];
        vec3 axis = {orientation[0],orientation[1],orientation[2]};
        FIXED distance = dotvec3(&dir, &axis);
        if (distance > obb->size.asArray[i]) {
        return false;
         }
        if (distance < -obb->size.asArray[i]) {
        return false;
        }
    }
    return true;
}

vec3 ClosestPointOBB(OBB *obb, vec3 *point) {
	vec3 result = obb->position;
	vec3 dir = vec3sub(point, &obb->position);

	for (int i = 0; i < 3; ++i) {
		const FIXED* orientation = &obb->orientation.asArray[i * 3];
		vec3 axis = {orientation[0], orientation[1], orientation[2]};

		FIXED distance = dotvec3ult(&axis, &dir);
        switch(i) {
            case 0 :
            jo_printf(0, 22, "distance : %d         ", jo_fixed2int(distance));
            jo_printf(0,14, "length : %d          ", jo_fixed2int(obb->size.asArray[i]));
            break;
            case 1 : 
            jo_printf(0, 24, "distance : %d      ", jo_fixed2int(distance));
            jo_printf(0,16, "length : %d       ", jo_fixed2int(obb->size.asArray[i]));
            break;
            case 2 : 
            jo_printf(0, 26, "distance : %d         ", jo_fixed2int(distance));
            jo_printf(0,18, "length : %d         ", jo_fixed2int(obb->size.asArray[i]));
            break;
        }
		if (distance > obb->size.asArray[i]) {
			distance = obb->size.asArray[i];
		}
		if (distance < -obb->size.asArray[i]) {
			distance = -obb->size.asArray[i];
		}

        vec3 tbm = vec3flt(&axis, distance);
		result = vec3add(&result, &tbm);
       // printf("distance : %f \n " ,distance);
	}

	return result;
}


bool PointOnPlane(vec3 *point, Plane *plane) {
    FIXED dot = dotvec3(point, &plane->normal);
    // To make this more robust, use an epsilon check
    // The CMP macro performs epsilon tests:
    // CMP(dot - plane.distance, 0.0f)
    return dot - plane->distance == 0;
}

vec3 ClosestPointPlane(Plane *plane, vec3 *point) {
    FIXED distance = dotvec3ult(&plane->normal, point) - plane->distance;
	// If the plane normal wasn't normalized, we'd need this:
	//distance = slDivFX(distance, dotvec3(&plane->normal, &plane->normal));
    vec3 tba = vec3flt(&plane->normal,distance);
    return vec3sub(point, &tba);
}


vec3 ClosestPointLine(ColLine *line, vec3 *point) {
    vec3 lVec = vec3sub(&line->end, &line->start); // ColLine Vector
    vec3 dbub = vec3sub(point, &line->start);
    FIXED t = slDivFX(dotvec3ult(&dbub, &lVec), dotvec3ult(&lVec, &lVec));
    t = MAX(t, 0); // Clamp to 0
    t = MIN(t, JO_FIXED_1); // Clamp to 1
    vec3 tbf = vec3flt(&lVec, t);
    return vec3add(&line->start, &tbf);
}

bool PointOnLine(vec3 *point, ColLine *line) {
    vec3 closest = ClosestPointLine(line, point);
    vec3 tbs = vec3sub(&closest, point);
    FIXED distanceSq = magsqvec3(&tbs);
    // Consider using an epsilon test here!
    // CMP(distanceSq, 0.0f);
    return distanceSq == 0;
}


bool PointOnRay(vec3 *point, Rray *ray) {
    if (point == &ray->origin) {
    return true;
    }

    vec3 norm = vec3sub(point, &ray->origin);
    normalizevec3(&norm);
    // We assume the ray direction is normalized
    FIXED diff = dotvec3(&norm, &ray->direction);
    // If BOTH vectors point in the same direction,
    // their dot product (diff) should be 1
    return diff == JO_FIXED_1; // Consider using epsilon!
}

vec3 ClosestPointRay(Rray *ray, vec3 *point) {
    vec3 tbf = vec3sub(point, &ray->origin);
    FIXED t = dotvec3(&tbf, &ray->direction);
    // We assume the direction of the ray is normalized
    t = MAX(t, 0);
    vec3 retvec3 = vec3add(&ray->origin, &ray->direction);
    return vec3flt(&retvec3, t);
}

bool SphereSphere(Sphere *s1, Sphere *s2) { //check a sphere against a sphere
    FIXED radiiSum = s1->radius + s2->radius;
    vec3 tbs = vec3sub(&s1->position, &s2->position);
    FIXED sqDistance = magsqvec3(&tbs);
    return sqDistance < radiiSum;
}

bool SphereAABB(Sphere *sphere, AABB *aabb) { //check a sphere against an axis aligned bounding box
    vec3 closestPoint = ClosestPointAABB(aabb, &sphere->position);
    vec3 tbs = vec3sub(&sphere->position, &closestPoint);
    FIXED distSq = magsqvec3(&tbs);
    FIXED radiusSq = toFIXED(jo_fixed2int(sphere->radius) * jo_fixed2int(sphere->radius));
    return distSq < radiusSq;
}

bool SphereOBB(Sphere *sphere, OBB *obb) { //check a sphere against an oriented bounding box
	vec3 closestPoint = ClosestPointOBB(obb, &sphere->position);
    vec3 tbs = vec3sub(&sphere->position, &closestPoint);
	FIXED distSq = (magsqvec3(&tbs));
	FIXED radiusSq = sphere->radius;
    jo_printf(0,8, "dist : %d         ", jo_fixed2int(distSq));

    return distSq<=radiusSq;
}

bool SpherePlane(Sphere *sphere, Plane *plane) { //check a sphere against an infinite plane
    vec3 closestPoint = ClosestPointPlane(plane, &sphere->position);
    vec3 tbs = vec3sub(&sphere->position, &closestPoint);
    FIXED distSq = magsqvec3ult(&tbs)*100*2;
    FIXED radiusSq = sphere->radius;
    jo_printf(0, 10, "distance : %d      \nradius : %d        ",jo_fixed2int(distSq),jo_fixed2int(radiusSq));
    return ABS(distSq) < radiusSq;
}

bool AABBAABB(AABB *aabb1, AABB *aabb2) {
    vec3 aMin = GetMin(aabb1);
    vec3 aMax = GetMax(aabb1);
    vec3 bMin = GetMin(aabb2);
    vec3 bMax = GetMax(aabb2);
    return (aMin.x <= bMax.x && aMax.x >= bMin.x) && (aMin.y <= bMax.y && aMax.y >= bMin.y) && (aMin.z <= bMax.z && aMax.z >= bMin.z);
}

Interval GetIntervalAABB(AABB *rect, vec3 *axis) {  
vec3 i = GetMin(rect);
vec3 a = GetMax(rect);

vec3 vertex[8] = {
    (vec3){i.x, a.y, a.z},
    (vec3){i.x, a.y, i.z},
    (vec3){i.x, i.y, a.z},
    (vec3){i.x, i.y, i.z},
    (vec3){a.x, a.y, a.z},
    (vec3){a.x, a.y, i.z},
    (vec3){a.x, i.y, a.z},
    (vec3){a.x, i.y, i.z}
};

Interval result;
result.min = result.max = dotvec3(axis, &vertex[0]);

for (int i = 1; i < 8; ++i) {
    FIXED projection = dotvec3(axis, &vertex[i]);
    result.min = (projection < result.min) ? projection : result.min;
    result.max = (projection > result.max) ? projection : result.max;
}
return result;
}

Interval GetIntervalOBB(OBB *obb, vec3 *axis) {
    vec3 vertex[8];
    vec3 C = obb->position; // OBB Center
    vec3 E = obb->size; // OBB Extents
    const FIXED* o = obb->orientation.asArray;
    vec3 A[] = { // OBB Axis
    (vec3){o[0], o[1], o[2]},
    (vec3){o[3], o[4], o[5]},
    (vec3){o[6], o[7], o[8]},
};
vec3 ml0 = vec3flt(&A[0], E.asArray[0]);
vec3 ml1 = vec3flt(&A[1], E.asArray[1]);
vec3 ml2 = vec3flt(&A[2], E.asArray[2]);
vec3 ml1pml2 = vec3add(&ml1,&ml2);
vec3 ml1sml2 = vec3sub(&ml1,&ml2);

vec3 mlalpha = vec3add(&ml0,&ml1pml2);
vec3 mlbeta = vec3sub(&ml0,&ml1pml2);
vec3 mldelta = vec3add(&ml0,&ml1sml2);
vec3 mlgamma = vec3sub(&ml0,&ml1sml2);


vertex[0] = vec3add(&C, &mlalpha);
vertex[1] = vec3sub(&C, &mlalpha);
vertex[2] = vec3add(&C, &mlbeta);
vertex[3] = vec3add(&C, &mldelta);
vertex[4] = vec3sub(&C, &mlgamma);
vertex[5] = vec3add(&C, &mlgamma);
vertex[6] = vec3sub(&C, &mldelta);
vertex[7] = vec3sub(&C, &mlbeta);
Interval result;
result.min = result.max = dotvec3(axis, &vertex[0]);

for (int i = 1; i < 8; ++i) {
    FIXED projection = dotvec3(axis, &vertex[i]);
    result.min = (projection < result.min) ? projection : result.min;
    result.max = (projection > result.max) ? projection : result.max;
}
return result;
}

bool OverlapOnAxisAABBOBB(AABB *aabb, OBB *obb, vec3 *axis) {
Interval a = GetIntervalAABB(aabb, axis);
Interval b = GetIntervalOBB(obb, axis);
return ((b.min <= a.max) && (a.min <= b.max));
}

bool AABBOBB(AABB *aabb, OBB *obb) {
    const FIXED* o = obb->orientation.asArray;
    vec3 test[15] = {
        (vec3){1, 0, 0}, // AABB axis 1
        (vec3){0, 1, 0}, // AABB axis 2
        (vec3){0, 0, 1}, // AABB axis 3
        (vec3){o[0], o[1], o[2]}, // OBB axis 1
        (vec3){o[3], o[4], o[5]}, // OBB axis 2
        (vec3){o[6], o[7], o[8]} // OBB axis 3
};
for (int i = 0; i < 3; ++i) { // Fill out rest of axis
    test[6 + i * 3 + 0] = cross(&test[i], &test[0]);
    test[6 + i * 3 + 1] = cross(&test[i], &test[1]);
    test[6 + i * 3 + 2] = cross(&test[i], &test[2]);
    }
for (int i = 0; i < 15; ++i) {
    if (!OverlapOnAxisAABBOBB(aabb, obb, &test[i])) {
    return false; // Seperating axis found
}
}
return true; // Seperating axis not found
}

bool OBBPlane(OBB *obb, Plane *plane) {
    const FIXED* o = obb->orientation.asArray;
vec3 rot[] = { // rotation / orientation
    (vec3){o[0], o[1], o[2]},
    (vec3){o[3], o[4], o[5]},
    (vec3){o[6], o[7], o[8]},
};
vec3 normal = plane->normal;
    FIXED pLen = toFIXED(jo_fixed2int(obb->size.x) * jo_fixed2int(ABS(dotvec3(&normal, &rot[0]))) + jo_fixed2int(obb->size.y) * jo_fixed2int(ABS(dotvec3(&normal, &rot[1]))) + jo_fixed2int(obb->size.z) * jo_fixed2int(ABS(dotvec3(&normal, &rot[2]))));
    FIXED dot = dotvec3(&plane->normal, &obb->position);
    FIXED dist = dot - plane->distance;
    return ABS(dist) <= pLen;
}

bool PlanePlane(Plane *plane1, Plane *plane2) {
    vec3 d = cross(&plane1->normal, &plane2->normal);
    return dotvec3(&d, &d) != 0; // Consider using an epsilon!
}

bool RaycastSphere(Sphere *sphere, Rray *ray, FIXED *output) {
	vec3 e = vec3sub(&sphere->position, &ray->origin);
	FIXED rSq = toFIXED(jo_fixed2int(sphere->radius) * jo_fixed2int(sphere->radius));

	FIXED eSq = magvec3(&e);
	FIXED a = dotvec3(&e, &ray->direction); // ray.direction is assumed to be normalized
	FIXED bSq = eSq - toFIXED(jo_fixed2int(a) * jo_fixed2int(a));
	FIXED f = toFIXED(slSquart(jo_fixed2int(ABS((rSq)- bSq))));

	// Assume normal intersection!
	FIXED t = a - f;

	// No collision has happened
	if (rSq - (eSq - toFIXED(jo_fixed2int(a) * jo_fixed2int(a))) < 0) {
		return false;
	}
	// Ray starts inside the sphere
	else if (eSq < rSq) {
		// Just reverse direction
		t = a + f;
	}
    *output = t;
	return true;
}


bool RaycastAABB(AABB *aabb, Rray *ray, FIXED *outter) {
	vec3 min = GetMin(aabb);
	vec3 max = GetMax(aabb);

	// Any component of direction could be 0!
	// Address this by using a small number, close to
	// 0 in case any of directions components are 0
	FIXED t1 = slDivFX((min.x - ray->origin.x),  ray->direction.x);
	FIXED t2 = slDivFX((max.x - ray->origin.x),  ray->direction.x);
	FIXED t3 = slDivFX((min.y - ray->origin.y),  ray->direction.y);
	FIXED t4 = slDivFX((max.y - ray->origin.y),  ray->direction.y);
	FIXED t5 = slDivFX((min.z - ray->origin.z),  ray->direction.z);
	FIXED t6 = slDivFX((max.z - ray->origin.z),  ray->direction.z);

	FIXED tmin = MAX(MAX(MIN(t1, t2), MIN(t3, t4)), MIN(t5, t6));
	FIXED tmax = MIN(MIN(MAX(t1, t2), MAX(t3, t4)), MAX(t5, t6));

	// if tmax < 0, ray is intersecting AABB
	// but entire AABB is behing it's origin
if (tmax< 0) {
return false;
}
if (tmin>tmax) {
return false;
}
if (tmin< 0) {
    *outter = tmin;
    return true;
}
    *outter = tmax;
	return true;
}

bool RaycastOBB(OBB *obb, Rray *ray, FIXED *output) {
    const FIXED* o = obb->orientation.asArray;
    const FIXED* size = obb->size.asArray;
    // X, Y and Z axis of OBB
    vec3 Xx = {o[0], o[1], o[2]};
    vec3 Yy = {o[3], o[4], o[5]};
    vec3 Zz = {o[6], o[7], o[8]};
    
    vec3 p = vec3sub(&obb->position, &ray->origin);

    vec3 f = {
    dotvec3(&Xx, &ray->direction),
    dotvec3(&Yy, &ray->direction),
    dotvec3(&Zz, &ray->direction)
    };
    
    vec3 e = {
    dotvec3(&Xx, &p),
    dotvec3(&Yy, &p),
    dotvec3(&Zz, &p)
    };

    FIXED t[6] = { 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < 3; i++) {
            if (cmp(f.asArray[i], 0)) {
            if (-e.asArray[i] - size[i]>0 || -e.asArray[i] + size[i]<0) {
            return false;
        }
        f.asArray[i] = JO_FIXED_EPSILON; // Avoid div by 0!
    }
        t[i * 2 + 0] = slDivFX((e.asArray[i] + size[i]), f.asArray[i]); // min
        t[i * 2 + 1] = slDivFX((e.asArray[i] - size[i]), f.asArray[i]); // max
    }
    FIXED tmin = MAX(MAX(MIN(t[0], t[1]),MIN(t[2], t[3])),MIN(t[4], t[5]));
    FIXED tmax = MIN(MIN(MAX(t[0], t[1]),MAX(t[2], t[3])),MAX(t[4], t[5]));
    if (tmax< 0) {
        return false;
    }

    if (tmin>tmax) {
        return false;
    }

    if (tmin< 0) {
        *output = tmin;
        return true;
        }
*output = tmax;
return true;
}

FIXED RaycastPlane(Plane *plane, Rray *ray) {
	FIXED nd = dotvec3(&ray->direction, &plane->normal);
	FIXED pn = dotvec3(&ray->origin, &plane->normal);
	// nd must be negative, and not 0
	// if nd is positive, the ray and plane normals
	// point in the same direction. No intersection.
	if (nd >= 0) {
		return -JO_FIXED_1;
	}

	FIXED t = slDivFX((plane->distance - pn), nd);
    //printf("tval %f \n", t);
	// t must be positive
	if (t >= 0) {
		return t;
	}

	return -JO_FIXED_1;
}

bool LinetestSphere(Sphere *sphere, ColLine *line) {
    vec3 closest = ClosestPointLine(line, &sphere->position);
    vec3 tbs = vec3sub(&sphere->position, &closest);
    FIXED distSq = magsqvec3(&tbs);
    return distSq <= (sphere->radius);
}   

bool LinetestAABB(AABB *aabb, ColLine *line) {
	Rray ray;
	ray.origin = line->start; 
    vec3 tbs = vec3sub(&line->end, &line->start);
	ray.direction = normalizedvec3(&tbs);
	FIXED output;
	if (!RaycastAABB(aabb, &ray, &output)) {
		return false;
	}
   
	return (output >= 0 && output <= lengthline3d(line));
}

bool LinetestOBB(OBB *obb, ColLine *line) {
if(PointInOBB(&line->start,obb) || PointInOBB(&line->end,obb)) {
    return true;
}
Rray ray;
ray.origin = line->start;
vec3 tbs = vec3sub(&line->end, &line->start);
ray.direction = normalizedvec3(&tbs);
FIXED output;
if(!RaycastOBB(obb, &ray, &output)) {
    return false;
}

return output > 0 && output <= lengthline3d(line);
}

bool LinetestPlane(Plane *plane, ColLine *line) {
    vec3 ab = vec3sub(&line->end, &line->start);
    FIXED nA = dotvec3(&plane->normal, &line->start);
    FIXED nAB = dotvec3(&plane->normal, &ab);
    // If the line and plane are parallel, nAB will be 0
    // This will cause a divide by 0 exception below
    // If you plan on testing parallel lines and planes
    // it is sage to early out when nAB is 0.
    FIXED t = slDivFX((plane->distance - nA), nAB);
return t >= 0 && t <= JO_FIXED_1;
}

bool PointInTriangle(vec3 *p, Triangle *t) {
	// Move the triangle so that the point is  
	// now at the origin of the triangle
	vec3 a = vec3sub(&t->a, p);
	vec3 b = vec3sub(&t->b, p);
	vec3 c = vec3sub(&t->c, p);

	// The point should be moved too, so they are both
	// relative, but because we don't use p in the
	// equation anymore, we don't need it!
	// p -= p; // This would just equal the zero vector!

	vec3 normPBC = crossult(&b, &c); // Normal of PBC (u)
	vec3 normPCA = crossult(&c, &a); // Normal of PCA (v)
	vec3 normPAB = crossult(&a, &b); // Normal of PAB (w)

	// Test to see if the normals are facing 
	// the same direction, return false if not
	if (dotvec3ult(&normPBC, &normPCA) < 0) {
		return false;
	}
	else if (dotvec3ult(&normPBC, &normPAB) < 0) {
		return false;
	}

	// All normals facing the same way, return true
	return true;
}
Plane PlaneFromTriangle(Triangle *t) {
	Plane result;
    vec3 c1 = vec3sub(&t->b, &t->a);
    vec3 c2 = vec3sub(&t->c, &t->a);
    vec3 c4 = vec3sub(&c1, &c2);
    vec3 c3 = crossult(&c1, &c2);
	result.normal = normalizedvec3(&c3);
	result.distance = dotvec3ult(&t->a, &result.normal);
    //printf("distance %f \n", result.distance);
	return result;
}

vec3 ClosestPointTriangle(Triangle *t, vec3 *p, int *out) {
	Plane plane = PlaneFromTriangle(t);
	vec3 closest = ClosestPointPlane(&plane, p);

	// Closest point was inside triangle
	if (PointInTriangle(&closest, t)) {
        *out = 1;
	//	return closest;
	}

    ColLine lina = constructline(&t->a, &t->b);
    ColLine linb = constructline(&t->b, &t->c);
    ColLine linc = constructline(&t->c, &t->a);

	vec3 c1 = ClosestPointLine(&lina, &closest); // ColLine AB
	vec3 c2 = ClosestPointLine(&linb, &closest); // ColLine BC
	vec3 c3 = ClosestPointLine(&linc, &closest); // ColLine CA

    vec3 clos1 = vec3sub(&closest, &c1);
    vec3 clos2 = vec3sub(&closest, &c2);
    vec3 clos3 = vec3sub(&closest, &c3);

	FIXED magSq1 = magsqvec3(&clos1);
	FIXED magSq2 = magsqvec3(&clos2);
	FIXED magSq3 = magsqvec3(&clos3);

    *out = 1;

	if (magSq1 < magSq2 && magSq1 < magSq3) {
		return c1;
	}
	else if (magSq2 < magSq1 && magSq2 < magSq3) {
		return c2;
	}

	return c3;
}

bool TriangleSphere(Triangle *t, Sphere *s) {
    int isit = 0;
	vec3 closest = ClosestPointTriangle(t, &s->position,&isit);
    //if(LinetestSphere)
    //if(PointOnPlane()) = 
    vec3 tbs = vec3sub(&s->position, &closest);
	FIXED magSq = magsqvec3ult(&tbs); //else if(isit == 0) magsqvec3ult(&tbs);
    jo_printf(0,10,"distance : %d    \nradius : %d   ",jo_fixed2int(magSq),jo_fixed2int(s->radius));
    return magSq <= s->radius;
}

Interval GetIntervalTriangle(Triangle *triangle, vec3 *axis) {
    Interval result;
    result.min = dotvec3(axis, &triangle->points[0]);
    result.max = result.min;
        for (int i = 1; i < 3; ++i) {
            FIXED value = dotvec3(axis, &triangle->points[i]);
            result.min = MIN(result.min, value);
            result.max = MAX(result.max, value);
        }
return result;
}

bool OverlapOnAxisTriangleAABB(AABB *aabb, Triangle *triangle, vec3 *axis) {
    Interval a = GetIntervalAABB(aabb, axis);
    Interval b = GetIntervalTriangle(triangle, axis);
    return ((b.min <= a.max) && (a.min <= b.max));
}

bool TriangleAABB(Triangle *t, AABB *a) {
    vec3 f0 = vec3sub(&t->b, &t->a);
    vec3 f1 = vec3sub(&t->c, &t->b);
    vec3 f2 = vec3sub(&t->a, &t->c);

    vec3 u0 = {JO_FIXED_1, 0, 0};
    vec3 u1 = {0, JO_FIXED_1, 0};
    vec3 u2 = {0, 0, JO_FIXED_1};

    vec3 test[13] = {
        u0, // AABB Axis 1
        u1, // AABB Axis 2
        u2, // AABB Axis 3
        cross(&f0, &f1), cross(&u0, &f0), cross(&u0, &f1), cross(&u0, &f2), cross(&u1, &f0), cross(&u1, &f1), cross(&u1, &f2), cross(&u2, &f0), cross(&u2, &f1), cross(&u2, &f2)
};
    for (int i = 0; i < 13; ++i) {
        if (!OverlapOnAxisTriangleAABB(a, t, &test[i])) {
        return false; // Separating axis found
        }
    }
    return true; // Separating axis not found
}

bool OverlapOnAxisTriangleOBB(OBB *obb, Triangle *triangle, vec3 *axis) {
    Interval a = GetIntervalOBB(obb, axis);
    Interval b = GetIntervalTriangle(triangle, axis);
    return ((b.min <= a.max) && (a.min <= b.max));
}

bool TriangleOBB(Triangle *t, OBB *o) {
    vec3 f0 = vec3sub(&t->b, &t->a);
    vec3 f1 = vec3sub(&t->c, &t->b);
    vec3 f2 = vec3sub(&t->a, &t->c);

    const FIXED* orientation = o->orientation.asArray;
    vec3 u0 = {orientation[0], orientation[1], orientation[2]}; 
    vec3 u1 = {orientation[3], orientation[4], orientation[5]}; 
    vec3 u2 = {orientation[6], orientation[7], orientation[8]};

    vec3 test[13] = {
        u0, // OBB Axis 1
        u1, // OBB Axis 2
        u2, // OBB Axis 3
        cross(&f0, &f1), // Normal of the Triangle
        cross(&u0, &f0), cross(&u0, &f1), cross(&u0, &f2),
        cross(&u1, &f0), cross(&u1, &f1), cross(&u1, &f2),
        cross(&u2, &f0), cross(&u2, &f1), cross(&u2, &f2)
    };

    for (int i = 0; i < 13; i++) {
        if (!OverlapOnAxisTriangleOBB(o, t, &test[i])) {
        return false; // Separating axis found
        }
    }
    return true; // Separating axis not found
}

bool TrianglePlane(Triangle *t, Plane *p) {
    FIXED side1 = PlaneEquation(&t->a, p);
    FIXED side2 = PlaneEquation(&t->b, p);
    FIXED side3 = PlaneEquation(&t->c, p);

    if (cmp(side1, 0) && cmp(side2, 0) && cmp(side3, 0)) {
        return true;
    }

    if (side1 > 0 && side2 > 0 && side3 > 0) {
        return false;
    }

    if (side1 < 0 && side2 < 0 && side3 < 0) {
        return false;
    }
return true; // Intersection
}

bool OverlapOnAxisTriangleTriangle(Triangle *t1, Triangle *t2, vec3 *axis) {
    Interval a = GetIntervalTriangle(t1, axis);
    Interval b = GetIntervalTriangle(t2, axis);
    return ((b.min <= a.max) && (a.min <= b.max));
}

bool TriangleTriangle(Triangle *t1, Triangle *t2) {
    vec3 t1_f0 = vec3sub(&t1->b, &t1->a); // Triangle 1, Edge 0
    vec3 t1_f1 = vec3sub(&t1->c, &t1->b); // Triangle 1, Edge 1
    vec3 t1_f2 = vec3sub(&t1->a, &t1->c); // Triangle 1, Edge 2
    vec3 t2_f0 = vec3sub(&t2->b, &t2->a); // Triangle 2, Edge 0
    vec3 t2_f1 = vec3sub(&t2->c, &t2->b); // Triangle 2, Edge 1
    vec3 t2_f2 = vec3sub(&t2->a, &t2->c); // Triangle 2, Edge 2

    vec3 axisToTest[] = {
        cross(&t1_f0, &t1_f1),
        cross(&t2_f0, &t2_f1),
        cross(&t2_f0, &t1_f0), cross(&t2_f0, &t1_f1),
        cross(&t2_f0, &t1_f2), cross(&t2_f1, &t1_f0),
        cross(&t2_f1, &t1_f1), cross(&t2_f1, &t1_f2),
        cross(&t2_f2, &t1_f0), cross(&t2_f2, &t1_f1),
        cross(&t2_f2, &t1_f2),
        };
    for (int i = 0; i < 11; i++) {
        if (!OverlapOnAxisTriangleTriangle(t1, t2, &axisToTest[i])) {
        return false; // Seperating axis found
        }
    }
    return true; // Seperating axis not found   
}

vec3 SatCrossEdge(vec3 *a, vec3 *b, vec3 *c, vec3 *d) {
    vec3 ab = vec3sub(a, b);
    vec3 cd = vec3sub(c, d);
    vec3 result = cross(&ab, &cd);

    if (!cmp(magsqvec3(&result), 0)) {
            return result; // Not parallel!
        } else { // ab and cd are paralle
    vec3 tbs = vec3sub(c, a);
    vec3 axis = cross(&ab, &tbs);
    result = cross(&ab, &axis);

    if (!cmp(magsqvec3(&result), 0)) {
    return result; // Not parallel
        }
    }

    return (vec3){0,0,0};
}

bool TriangleTriangleRobust(Triangle *t1, Triangle *t2) {
vec3 axisToTest[] = {
    // Triangle 1, Normal
    SatCrossEdge(&t1->a, &t1->b, &t1->b, &t1->c),
    // Triangle 2, Normal
    SatCrossEdge(&t2->a, &t2->b, &t2->b, &t2->c),
    SatCrossEdge(&t2->a, &t2->b, &t1->a, &t1->b),
    SatCrossEdge(&t2->a, &t2->b, &t1->b, &t1->c),
    SatCrossEdge(&t2->a, &t2->b, &t1->c, &t1->a),
    SatCrossEdge(&t2->b, &t2->c, &t1->a, &t1->b),
    SatCrossEdge(&t2->b, &t2->c, &t1->b, &t1->c),
    SatCrossEdge(&t2->b, &t2->c, &t1->c, &t1->a),
    SatCrossEdge(&t2->c, &t2->a, &t1->a, &t1->b),
    SatCrossEdge(&t2->c, &t2->a, &t1->b, &t1->c),
    SatCrossEdge(&t2->c, &t2->a, &t1->c, &t1->a)
    };
    for (int i = 0; i < 11; ++i) {
        if (!OverlapOnAxisTriangleTriangle(t1, t2, &axisToTest[i])) {
            if (!cmp(magsqvec3(&axisToTest[i]), 0)) {
                return false; // Seperating axis found
                    }
                }
            }
return true; // Seperating axis not found
}

vec3 Barycentric(vec3 *p, Triangle *t) {
	vec3 ap = vec3sub(p, &t->a);
	vec3 bp = vec3sub(p, &t->b);
	vec3 cp = vec3sub(p, &t->c);

	vec3 ab = vec3sub(&t->b, &t->a);
	vec3 ac = vec3sub(&t->c, &t->a);
	vec3 bc = vec3sub(&t->c, &t->b);
	vec3 cb = vec3sub(&t->b, &t->c);
	vec3 ca = vec3sub(&t->a, &t->c);

    vec3 proj1 = projvec3(&ab, &cb);
	vec3 v = vec3sub(&ab, &proj1);
	FIXED a = JO_FIXED_1 - slDivFX(dotvec3(&v, &ap), dotvec3(&v, &ab));

    vec3 proj2 = projvec3(&bc, &ac);
	v = vec3sub(&bc, &proj2);
	FIXED b = JO_FIXED_1 - slDivFX(dotvec3(&v, &bp), dotvec3(&v, &bc));

    vec3 proj3 = projvec3(&ca, &ab);
	v = vec3sub(&ca, &proj3);
	FIXED c = JO_FIXED_1 - slDivFX(dotvec3(&v, &cp), dotvec3(&v, &ca));

    vec3 resultoid = {a,b,c};
	return resultoid;
}

bool RaycastTriangle(Triangle *triangle, Rray *ray) {
	Plane plane = PlaneFromTriangle(triangle);
    FIXED t;
	if (RaycastPlane(&plane, ray) == -JO_FIXED_1) {
        return false;
	}
    t = RaycastPlane(&plane, ray);

    vec3 tbf = vec3add(&ray->direction,&ray->origin);
	vec3 result = vec3flt(&tbf,t);
	vec3 barycentric = Barycentric(&result, triangle);
	if (barycentric.x >= 0 && barycentric.x <= JO_FIXED_1 &&
		barycentric.y >= 0 && barycentric.y <= JO_FIXED_1 &&
		barycentric.z >= 0 && barycentric.z <= JO_FIXED_1) {
        

		return true;
	}

	return false;
}

