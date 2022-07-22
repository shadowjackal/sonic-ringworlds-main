#ifndef JKLCOL_H
#define JKLCOL_H
//requires jklmath

typedef struct Plane {
vec3 normal;
FIXED distance;
} Plane;

typedef struct Sphere {
    vec3 position;
    FIXED radius;
} Sphere;

typedef struct tri {
    union {
    struct {
        vec3 a;
        vec3 b;
        vec3 c;
    };
    vec3 points[3];
    FIXED values[9];
    }
} tri;

extern tri tridef(const vec3 pnt1, const vec3 pnt2, const vec3 pnt3);
extern void tridefptr(tri* returntri,vec3 *pnt1, const vec3 *pnt2, const vec3 *pnt3);

extern Plane PlaneFromTri(tri t);
extern vec3 ClosestPoint(const Plane *plane, const vec3 *point);
extern bool SpherePlane(const Sphere *s, const Plane *p);

extern bool PointInTri(tri *t, vec3 *point);
extern vec3 ClosestPointLine(const vec3 *start,const vec3 *end, const vec3 point);

extern vec3 ClosestPointTri(tri *t, vec3 *point);
extern bool TriangleSphere(tri *t, Sphere *s);
extern vec3 linetodir(vec3 *from, vec3 *to);
extern PDATA   coltri2mesh(tri *t);
#endif