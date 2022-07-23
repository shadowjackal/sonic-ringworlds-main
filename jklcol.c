#include <jo/jo.h>
#include "jklmath.h"
#include "jklcol.h"

tri tridef(const vec3 pnt1, const vec3 pnt2, const vec3 pnt3) {
    tri retri;
    retri.a = pnt1;
    retri.b = pnt2;
    retri.c = pnt3;
    return retri;
}

void tridefptr(tri* returntri,vec3 *pnt1, const vec3 *pnt2, const vec3 *pnt3) {
    returntri->a = *pnt1;
    returntri->b = *pnt2;
    returntri->c = *pnt3;
}

Plane PlaneFromTri(tri t) {
    Plane result;
    tri retri = t;
    VECTOR UVEC;
    VECTOR VVEC;
    UVEC[0] = retri.b.asArray[0] - retri.a.asArray[0];
	UVEC[1] = retri.b.asArray[1] - retri.a.asArray[1];
	UVEC[2] = retri.b.asArray[2] - retri.a.asArray[2];
	VVEC[0] = retri.c.asArray[0] - retri.a.asArray[0];
	VVEC[1] = retri.c.asArray[1] - retri.a.asArray[1];
	VVEC[2] = retri.c.asArray[2] - retri.a.asArray[2];
	result.normal.asArray[0] = (slMulFX(UVEC[1], VVEC[2]) - (slMulFX(UVEC[2], VVEC[1])));
	result.normal.asArray[1] = (slMulFX(UVEC[2], VVEC[0]) - (slMulFX(UVEC[0], VVEC[2])));
	result.normal.asArray[2] = (slMulFX(UVEC[0], VVEC[1]) - (slMulFX(UVEC[1], VVEC[0])));
    result.normal = normalize(result.normal);
    result.distance = slInnerProduct((FIXED *)result.normal.asArray, (FIXED *)retri.a.asArray);
    return result;
}


vec3 ClosestPoint(const Plane *plane, const vec3 *point){
FIXED dot = dotvec3(&plane->normal, point);
FIXED distance = dot - plane->distance;
return vec3subn(*point, vec3flt(&plane->normal,distance));
}

bool SpherePlane(const Sphere *s, const Plane *p) {
    vec3 closestPoint = ClosestPoint(p,&s->position);
    int dif = (JO_ABS(s->position.asArray[X] - closestPoint.asArray[X]) + JO_ABS(s->position.asArray[Y] - closestPoint.asArray[Y]) + JO_ABS(s->position.asArray[Z] - closestPoint.asArray[Z]));
    return dif < s->radius;
}


bool PointInTri(tri *t, vec3 *point) {    
    vec3 a, b, c;
    a = vec3sub(&t->a, point);
    b = vec3sub(&t->b, point);
    c = vec3sub(&t->c, point);
    
    // Prevent fixed-point overflow
    a = normalize(a);
    b = normalize(b);
    c = normalize(c);
    
    vec3 u, v, w;
    u = cross_fixed(b, c);
    v = cross_fixed(c, a);
    w = cross_fixed(a, b);

    slPrintFX(slInnerProduct(&u.asArray, &v.asArray),slLocate(0,5));
    slPrintFX(slInnerProduct(&u.asArray, &w.asArray),slLocate(0,6));

    
    if (slInnerProduct((FIXED *)&u.asArray, (FIXED *)&v.asArray) < 0)
        return false;
    if (slInnerProduct((FIXED *)&u.asArray, (FIXED *)&w.asArray) < 0)
        return false;
    slPrint("in tri", slLocate(0,24));
    return true;
}

    

vec3 ClosestPointLine(const vec3 *start,const vec3 *end, const vec3 point) {
    vec3 AP, AB;
    AP = vec3sub(&point, start);
    AB = vec3sub(end, start);
    
    FIXED ap_dot_ab = slInnerProduct(AP.asArray, AB.asArray);
    
    // Point is behind start of the segment, so perpendicular distance is not viable
    if (ap_dot_ab <= toFIXED(0.0)) {
        return *start;
    }
    
    FIXED ab_dot_ab = slInnerProduct(AB.asArray, AB.asArray);
    
    // Point is beyond end of the segment, so perpendicular distance is not viable
    if (ap_dot_ab >= ab_dot_ab) {
        return *end;
    }
    
    FIXED t = slDivFX(ab_dot_ab, ap_dot_ab);
    
    vec3 out = vec3flt(&AB, t);
    return vec3add(&out, start);
}

vec3 ClosestPointTri(tri *t, vec3 *point) {
	Plane plane = PlaneFromTri(*t);
	vec3 closest = ClosestPoint(&plane, point);

	// Closest point was inside triangle
	if (PointInTri(t, &closest)) {
		return closest;
	}
    slPrint("out of tri", slLocate(0,24));

	vec3 c1 = ClosestPointLine(&t->a, &t->b, closest); // Line AB
	vec3 c2 = ClosestPointLine(&t->b, &t->c, closest); // Line BC
	vec3 c3 = ClosestPointLine(&t->c, &t->a, closest); // Line CA

	FIXED magSq1 = magsqvec3n(vec3sub(&closest,&c1));//(JO_ABS(closest.asArray[X] - c1.asArray[X]) + JO_ABS(closest.asArray[Y] - c1.asArray[Y]) + JO_ABS(closest.asArray[Z] - c1.asArray[Z]));
	FIXED magSq2 = magsqvec3n(vec3sub(&closest,&c2));//(JO_ABS(closest.asArray[X] - c2.asArray[X]) + JO_ABS(closest.asArray[Y] - c2.asArray[Y]) + JO_ABS(closest.asArray[Z] - c2.asArray[Z]));
	FIXED magSq3 = magsqvec3n(vec3sub(&closest,&c3));//(JO_ABS(closest.asArray[X] - c3.asArray[X]) + JO_ABS(closest.asArray[Y] - c3.asArray[Y]) + JO_ABS(closest.asArray[Z] - c3.asArray[Z]));

	if (magSq1 < magSq2 && magSq1 < magSq3) {
		return c1;
	}
	else if (magSq2 < magSq1 && magSq2 < magSq3) {
		return c2;
	}
	return c3;
}

bool TriangleSphere(tri *t, Sphere *s) {
vec3 closest = ClosestPointTri(t, &s->position);
FIXED magSq = (JO_ABS(closest.asArray[X] - s->position.asArray[X]) + JO_ABS(closest.asArray[Y] - s->position.asArray[Y]) + JO_ABS(closest.asArray[Z] - s->position.asArray[Z]));

return magSq < s->radius;
}

vec3 linetodir(vec3 *from, vec3 *to) {
    return normalize(vec3sub(to, from));
}


PDATA   coltri2mesh(tri *t) {
    PDATA retdata;
    retdata.nbPoint = 4;
    retdata.nbPolygon = 1;
    retdata.pntbl = jo_malloc(sizeof(POINT)*retdata.nbPoint);
    retdata.pltbl = jo_malloc(sizeof(POLYGON)*retdata.nbPolygon);
    retdata.attbl = jo_malloc(sizeof(ATTR) * retdata.nbPolygon);

    retdata.pntbl[0][0] = t->a.x;
    retdata.pntbl[0][1] = t->a.y;
    retdata.pntbl[0][2] = t->a.z;
    retdata.pntbl[1][0] = t->b.x;
    retdata.pntbl[1][1] = t->b.y;
    retdata.pntbl[1][2] = t->b.z;
    retdata.pntbl[2][0] = t->c.x;
    retdata.pntbl[2][1] = t->c.y;
    retdata.pntbl[2][2] = t->c.z;
    retdata.pntbl[3][0] = t->c.x;
    retdata.pntbl[3][1] = t->c.y;
    retdata.pntbl[3][2] = t->c.z;
    retdata.pltbl[0].Vertices[0] = 0;
    retdata.pltbl[0].Vertices[1] = 1;
    retdata.pltbl[0].Vertices[2] = 2;
    retdata.pltbl[0].Vertices[3] = 3;

    ATTR pro_attribute = ATTRIBUTE(Dual_Plane,SORT_CEN|No_Gouraud,No_Texture,C_RGB(31,31,0),0,MESHoff,sprPolygon,No_Option);
    retdata.attbl[0] = pro_attribute;
    return retdata;
}

