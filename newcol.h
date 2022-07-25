#ifndef __NEW_COL_H__
#define __NEW_COL_H__

typedef struct {
    Uint32 num_tris;
    POINT *points;
} Collision;

typedef struct {
    VECTOR normal;
    FIXED distance;
} Plane;

typedef struct {
    POINT pos;
    FIXED radius;
} Sphere;

PDATA   coltri2mesh(POINT *t);
void   mesh2coltri(jklmesh *m, Collision *out);
void    tridef(const POINT p1, const POINT p2, const POINT p3, POINT* out);
bool Collision_SphereColResolve(Sphere *sphere, Collision *col);
bool Collision_SpherePlaneResolve(Sphere *sphere, Plane *plane);
bool Collision_SphereCol_bool(Sphere *sphere, Collision *col);
bool Collision_SpherePlane_bool(Sphere *sphere, Plane *plane);

#endif
