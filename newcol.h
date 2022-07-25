#ifndef __NEW_COL_H__
#define __NEW_COL_H__

typedef struct {
    Uint32 num_tris;
    POINT *points;
} Collision;

typedef struct {
    POINT pos;
    FIXED radius;
} Sphere;

PDATA   coltri2mesh(POINT *t);
void    tridef(const POINT p1, const POINT p2, const POINT p3, POINT* out);
bool Collision_SphereColResolve(Sphere *sphere, Collision *col);

#endif
