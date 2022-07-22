/*
** Jo Sega Saturn Engine
** Copyright (c) 2012-2021, Johannes Fetz (johannesfetz@gmail.com)
** All rights reserved.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are met:
**     * Redistributions of source code must retain the above copyright
**       notice, this list of conditions and the following disclaimer.
**     * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in the
**       documentation and/or other materials provided with the distribution.
**     * Neither the name of the Johannes Fetz nor the
**       names of its contributors may be used to endorse or promote products
**       derived from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
** ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
** WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
** DISCLAIMED. IN NO EVENT SHALL Johannes Fetz BE LIABLE FOR ANY
** DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
** (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
** LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
** ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <jo/jo.h>
#include <string.h>
#include "ssv.h"

#define LWRAM 0x00200000

#define LWRAM_HEAP_SIZE 1310720


jo_camera               cam;
char                    *file_name;
int                     qicount;
int                     big_chungometer;
unsigned int            megainfo;

FIXED time;
jklmesh player;
jklmesh player0;

unsigned int rx;
unsigned int ry;
unsigned int rz;
int frame;
float steel = 2.0f;
int ja = 0;
jklmesh ssvlist0;
char *ssvarray[25];
FIXED timz;
FIXED zpitch;
FIXED zroll;
FIXED zyaw;
FIXED livec[XYZ];
struct ssvHeader {
    int totalvertices;
    int totalpolys;
} ssvHeader;

struct quad {
    int x[4];
    int y[4];
    int z[4];
};

typedef struct vec3
{
  union 
  {
    struct 
    {
      FIXED x;
      FIXED y;
      FIXED z;
    };
  FIXED asArray[3];
  };
} vec3;

typedef struct vec3ang
{
  union 
  {
    struct 
    {
      ANGLE x;
      ANGLE y;
      ANGLE z;
    };
  ANGLE asArray[3];
  };
} vec3ang;

typedef struct {
    VECTOR	 norm ;			/* Normal vector */
    int Vertices[4] ;		/* The vertex number that comprises the polygon */
} POLYGON32 ;



struct vertice {
	int x;
	int y;
	int z;
} vertice;

typedef struct playerobject {
    FIXED x,y,z;
    int state;
    ANGLE yrot;
    vec3ang orientation;
} playerobject;

playerobject sonic;

typedef float NORML[3];
typedef unsigned int face[4]; 

static Uint8 GourWork[4096];
static GOURAUDTBL Gour[4096];

Uint16 GRreal_ptr = 0xe000;



#define GRTBL(r,g,b) (((b&0x1f)<<10) | ((g&0x1f)<<5) | (r&0x1f))
static Uint16 GourTbl[32] = {
    GRTBL( 0,  0,  0),
    GRTBL( 0,  0,  0),
    GRTBL( 1,  1,  1),
    GRTBL( 1,  1,  1),
    GRTBL( 2,  2,  2),
    GRTBL( 2,  2,  2),
    GRTBL( 3,  3,  3),
    GRTBL( 3,  3,  3),
    GRTBL( 4,  4,  4),
    GRTBL( 4,  4,  4),
    GRTBL( 5,  5,  5),
    GRTBL( 5,  5,  5),
    GRTBL( 6,  6,  6),
    GRTBL( 6,  6,  6),
    GRTBL( 7,  7,  7),
    GRTBL( 7,  7,  7),
    GRTBL(11, 11, 11),
    GRTBL(13, 13, 13),
    GRTBL(14, 14, 14),
    GRTBL(16, 16, 16),
    GRTBL(17, 17, 17),
    GRTBL(18, 18, 18),
    GRTBL(19, 19, 19),
    GRTBL(20, 20, 20),
    GRTBL(21, 21, 21),
    GRTBL(22, 22, 22),
    GRTBL(23, 23, 23),
    GRTBL(25, 25, 25),
    GRTBL(26, 26, 26),
    GRTBL(26, 26, 26),
    GRTBL(27, 27, 27),
    GRTBL(27, 27, 27),
    GRTBL(29, 29, 29),
};

#define    GRrealBase    0xe000

jo_palette                  image_pal;
jo_palette                  image_pal2;


vec3 vec3add(vec3 *left, vec3 *right) {
    vec3 newvec3 = { left->x + right->x, left->y + right->y, left->z + right->z};
    return newvec3;
}

vec3 vec3addn(vec3 left, vec3 right) {
    vec3 newvec3 = { left.x + right.x, left.y + right.y, left.z + right.z};
    return newvec3;
}

vec3 vec3sub(vec3 *left, vec3 *right) {
    vec3 newvec3 = { left->x - right->x, left->y - right->y, left->z - right->z};
    return newvec3;
}

vec3 vec3subn(vec3 left, vec3 right) {
    vec3 newvec3 = { left.x - right.x, left.y - right.y, left.z - right.z};
    return newvec3;
}

vec3 vec3mul(vec3 *left, vec3 *right) {
    vec3 newvec3 = {toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->x)), toFIXED(jo_fixed2int(left->y) * jo_fixed2int(right->y)), toFIXED(jo_fixed2int(left->z) * jo_fixed2int(right->z))};
    return newvec3;
}

vec3 vec3flt(vec3 *left, FIXED flt) {
    vec3 newvec3 = {slMulFX(left->x, flt), slMulFX(left->y, flt), slMulFX(left->z, flt)};
    return newvec3;
}

FIXED dotvec3(vec3 *left, vec3 *right) {
return slMulFX(left->x, right->x) + slMulFX(left->y, right->y) + slMulFX(left->z, right->z);
}

FIXED dotvec3n(vec3 left, vec3 right) {
return slMulFX(left.x, right.x) + slMulFX(left.y, right.y) + slMulFX(left.z, right.z);
}

int dotvec3int(vec3 *left, vec3 *right) {
return ((left->x/JO_FIXED_1) * (right->x/JO_FIXED_1)) + ((left->y/JO_FIXED_1) * (right->y/JO_FIXED_1)) + ((left->z/JO_FIXED_1) * (right->z/JO_FIXED_1));
}

float dotvec3flt(vec3 *left, vec3 *right) {
return ((jo_fixed2float(left->x)) * (jo_fixed2float(right->x))) + ((jo_fixed2float(left->y)) * (jo_fixed2float(right->y))) + ((jo_fixed2float(left->z)) * (jo_fixed2float(right->z)));
}



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


FIXED magsqvec3(vec3 *v) {
return (slSquartFX(dotvec3(v, v)));
}

FIXED magsqvec3n(vec3 v) {
return (slSquartFX(dotvec3(&v, &v)));
}

int magsqint(vec3 *v) {
return (slSquart(dotvec3int(v, v)));
}

vec3 normalize(vec3 p)
{
    FIXED w = slSquartFX(slMulFX(p.asArray[0], p.asArray[0]) + slMulFX(p.asArray[1], p.asArray[1]) + slMulFX(p.asArray[2], p.asArray[2]));
    p.asArray[0] = slDivFX(p.asArray[0], w);
    p.asArray[1] = slDivFX(p.asArray[1], w);
    p.asArray[2] = slDivFX(p.asArray[2], w);
    return p;
}

vec3	cross_fixed(vec3 vector1, vec3 vector2)
{
    vec3 output;
	output.asArray[X] = slMulFX(vector1.asArray[Y], vector2.asArray[Z]) - slMulFX(vector1.asArray[Z], vector2.asArray[Y]);
	output.asArray[Y] = slMulFX(vector1.asArray[Z], vector2.asArray[X]) - slMulFX(vector1.asArray[X], vector2.asArray[Z]);
	output.asArray[Z] = slMulFX(vector1.asArray[X], vector2.asArray[Y]) - slMulFX(vector1.asArray[Y], vector2.asArray[X]);
    return output;
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
    slPrintFX(result.normal.x,slLocate(0,12));
    slPrintFX(result.normal.y,slLocate(0,13));
    slPrintFX(result.normal.z,slLocate(0,14));
    return result;
}


vec3 ClosestPoint(const Plane *plane, const vec3 *point){
FIXED dot = dotvec3(&plane->normal, point);
FIXED distance = dot - plane->distance;

//slPrintFX(vec3subn(*point, vec3flt(&plane->normal, distance)).x,slLocate(0,4));
//slPrintFX(vec3subn(*point, vec3flt(&plane->normal, distance)).y,slLocate(0,5));
//slPrintFX(vec3subn(*point, vec3flt(&plane->normal, distance)).z,slLocate(0,6));
return vec3subn(*point, vec3flt(&plane->normal,distance));
}

bool SpherePlane(const Sphere *s, const Plane *p) {
    vec3 closestPoint = ClosestPoint(p,&s->position);
    int dif = (JO_ABS(s->position.asArray[X] - closestPoint.asArray[X]) + JO_ABS(s->position.asArray[Y] - closestPoint.asArray[Y]) + JO_ABS(s->position.asArray[Z] - closestPoint.asArray[Z]));
    slPrintFX(dif,slLocate(0,17));
    return dif < s->radius;
}


tri testtri;

PDATA trimesh;

bool PointInTri(tri *t, vec3 *point) {
    vec3 a = vec3sub(&t->a,point); 
    vec3 b = vec3sub(&t->b,point); 
    vec3 c = vec3sub(&t->c,point);  

    vec3 u = cross_fixed(b, c);
    u.x = u.x>>8;
    u.y = u.y>>8;
    u.z = u.z>>8;

    vec3 v = cross_fixed(c, a);
    v.x = v.x>>8;
    v.y = v.y>>8;
    v.z = v.z>>8;

    vec3 w = cross_fixed(a, b);   
    w.x = w.x>>8;
    w.y = w.y>>8;
    w.z = w.z>>8;


    slPrintFX(dotvec3(&u, &v),slLocate(0,5));
    slPrintFX(dotvec3(&u, &w),slLocate(0,6));
    
    if (dotvec3(&u, &v) < toFIXED(0)) {
      return false;
    }

    if (dotvec3(&u, &w) < toFIXED(0)) {
      return false;
    }

    return true;
}
 
bool pntintri(tri *t, const vec3 *p) {
    vec3 a = vec3sub(&t->a, p);
    vec3 c = vec3sub(&t->c, p);
    vec3 b = vec3sub(&t->b, p);

    FIXED chk;
    //vec3 magsum = vec3sub()
    if(magsqint(&a) > 100) return false;
    
    vec3 normPBC = cross_fixed(b, c); // Normal of PBC (u)
    normPBC = normalize(normPBC);
    vec3 normPCA = cross_fixed(c, a); // Normal of PCA (v)
    normPCA = normalize(normPCA);
    vec3 normPAB = cross_fixed(a, b); // Normal of PAB (w)
    normPAB = normalize(normPAB);

    if (slInnerProduct((FIXED *)normPBC.asArray, (FIXED *)normPCA.asArray) < 1) {
        return false;
    } else if (slInnerProduct((FIXED *)normPBC.asArray, (FIXED *)normPAB.asArray) < 1) {
        return false;
    }
    return true;
}

vec3 ClosestPointLine(const vec3 *start,const vec3 *end, const vec3 point) {
	vec3 lVec = vec3sub(end, start); // Line Vector
	// Project "point" onto the "Line Vector", computing:
	// closest(t) = start + t * (end - start)
	// T is how far along the line the projected point is
	FIXED t = slDivFX(dotvec3n(vec3sub(&point, start), lVec), dotvec3(&lVec, &lVec));
	// Clamp t to the 0 to 1 range
	t = JO_MAX(t, toFIXED(0.0f));
	t = JO_MIN(t, toFIXED(1.0f));
	// Return projected position of t
    slPrintFX(vec3addn(*start,vec3flt(&lVec,t)).x,slLocate(0,22));
    slPrintFX(vec3addn(*start,vec3flt(&lVec,t)).y,slLocate(0,23));
    slPrintFX(vec3addn(*start,vec3flt(&lVec,t)).z,slLocate(0,24));

	return vec3addn(*start, vec3flt(&lVec, t));
}

vec3 ClosestPointTri(tri *t, vec3 *point) {
	Plane plane = PlaneFromTri(*t);
	vec3 closest = ClosestPoint(&plane, point);

	// Closest point was inside triangle
	if (PointInTri(t, &closest)) {
        slPrintFX(closest.x,slLocate(0,18));
        slPrintFX(closest.y,slLocate(0,19));
        slPrintFX(closest.z,slLocate(0,20));
		return closest;
	}

	vec3 c1 = ClosestPointLine(&t->a, &t->b, closest); // Line AB
	vec3 c2 = ClosestPointLine(&t->b, &t->c, closest); // Line BC
	vec3 c3 = ClosestPointLine(&t->c, &t->a, closest); // Line CA

	FIXED magSq1 = magsqvec3n(vec3sub(&closest,&c1));//(JO_ABS(closest.asArray[X] - c1.asArray[X]) + JO_ABS(closest.asArray[Y] - c1.asArray[Y]) + JO_ABS(closest.asArray[Z] - c1.asArray[Z]));
	FIXED magSq2 = magsqvec3n(vec3sub(&closest,&c2));//(JO_ABS(closest.asArray[X] - c2.asArray[X]) + JO_ABS(closest.asArray[Y] - c2.asArray[Y]) + JO_ABS(closest.asArray[Z] - c2.asArray[Z]));
	FIXED magSq3 = magsqvec3n(vec3sub(&closest,&c3));//(JO_ABS(closest.asArray[X] - c3.asArray[X]) + JO_ABS(closest.asArray[Y] - c3.asArray[Y]) + JO_ABS(closest.asArray[Z] - c3.asArray[Z]));

	if (magSq1 < magSq2 && magSq1 < magSq3) {
        slPrintFX(c1.x,slLocate(0,18));
        slPrintFX(c1.y,slLocate(0,19));
        slPrintFX(c1.z,slLocate(0,20));
		return c1;
	}
	else if (magSq2 < magSq1 && magSq2 < magSq3) {
        slPrintFX(c2.x,slLocate(0,18));
        slPrintFX(c2.y,slLocate(0,19));
        slPrintFX(c2.z,slLocate(0,20));
		return c2;
	}
        slPrintFX(c3.x,slLocate(0,18));
        slPrintFX(c3.y,slLocate(0,19));
        slPrintFX(c3.z,slLocate(0,20));
	return c3;
}

bool TriangleSphere(tri *t, Sphere *s) {
vec3 closest = ClosestPointTri(t, &s->position);
FIXED magSq = (JO_ABS(closest.asArray[X] - s->position.asArray[X]) + JO_ABS(closest.asArray[Y] - s->position.asArray[Y]) + JO_ABS(closest.asArray[Z] - s->position.asArray[Z]));

return magSq <= s->radius;
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


int                ssv_update(jklmesh *mesh, int frame) {
    if(mesh->framecount < frame) {
        return -1;
    }

    mesh->data.pntbl = mesh->rdata[frame].pntbl;
    
    for(int i = 0; i < mesh->data.nbPolygon; i++) {
        mesh->data.pltbl[i].norm[0] = -ANORMS[mesh->rdata[frame].fnorm[i]][0];
        mesh->data.pltbl[i].norm[1] = -ANORMS[mesh->rdata[frame].fnorm[i]][1];
        mesh->data.pltbl[i].norm[2] = -ANORMS[mesh->rdata[frame].fnorm[i]][2];

    }

    for(int i = 0; i < mesh->data.nbPoint; i++) {
        mesh->data.vntbl[i][0] = -ANORMS[mesh->rdata[frame].vnorm[i]][0];
        mesh->data.vntbl[i][1] = -ANORMS[mesh->rdata[frame].vnorm[i]][1];
        mesh->data.vntbl[i][2] = -ANORMS[mesh->rdata[frame].vnorm[i]][2];
    }
}

void                vec3orbit(vec3 *position, vec3 *target, FIXED distance, ANGLE angle) {
        position->x = target->x + slMulFX(slCos(angle),distance);
        position->z = target->z + slMulFX(slSin(angle),distance);
}


jklmesh                ssv_load(char *file_input, int textureoffset, bool gouraud, bool staticd)
{
    jklmesh              mesh;
   // mesh = (jklmesh *)jo_malloc_with_behaviour(sizeof(*mesh), JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    //loads in file

    char *stream = jo_fs_read_file(file_input,NULL);

    int offset = 0;
    int qi;
    int texids[2000];
    int colorrgb[2000][3];

    memcpy(&mesh.data.nbPoint, &stream[offset], sizeof(unsigned int));
    jo_printf(0,2,"vertex %d",mesh.data.nbPoint);
    offset += sizeof(unsigned int);

    memcpy(&mesh.data.nbPolygon, &stream[offset], sizeof(unsigned int));
    jo_printf(0,3,"polygon %d",mesh.data.nbPolygon);
    offset += sizeof(unsigned int);

    memcpy(&mesh.framecount, &stream[offset], sizeof(unsigned int));
    jo_printf(0,4,"frame %d",mesh.framecount);
    offset += sizeof(unsigned int);

    mesh.rdata = (JKLR *)jo_malloc_with_behaviour(mesh.framecount * sizeof(*mesh.rdata),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    mesh.data.pntbl = (POINT *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.pntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    mesh.data.vntbl = (VECTOR *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.vntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    for(int i = 0; i < mesh.framecount; i++) {
        mesh.rdata[i].pntbl = (POINT *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.pntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }
    
    
    mesh.data.pltbl = (POLYGON *)jo_malloc_with_behaviour(mesh.data.nbPolygon * sizeof(*mesh.data.pltbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    for(int i = 0; i < mesh.data.nbPolygon; i++) {
    for(int ii = 0; ii < 4; ii++){
        memcpy(&mesh.data.pltbl[i].Vertices[ii], &stream[offset], sizeof(short));
        offset += sizeof(short);
    }
        memcpy(&texids[i],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][0],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][1],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][2],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

    }


    ////memcpy(mesh.data.pntbl, &stream[offset], sizeof(POINT) * mesh.data.nbPoint);
    for(int i = 0; i < mesh.framecount; i++) {
    memcpy(mesh.rdata[i].pntbl, &stream[offset], sizeof(POINT) * mesh.data.nbPoint);
    offset += sizeof(POINT)*mesh.data.nbPoint;
    }

    mesh.data.pntbl = mesh.rdata[0].pntbl;


    for(int i = 0; i < mesh.framecount; i++) {
    mesh.rdata[i].fnorm = jo_malloc_with_behaviour(sizeof(Uint32)*mesh.data.nbPolygon,JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }

    for(int i = 0; i < mesh.framecount; i++) {
    for(int j = 0; j < mesh.data.nbPolygon; j++) {
        memcpy(&mesh.rdata[i].fnorm[j],&stream[offset],sizeof(Uint32));
        offset += sizeof(Uint32);
    }
    }
//
    for(int i = 0; i < mesh.framecount; i++) {
    mesh.rdata[i].vnorm = jo_malloc_with_behaviour(sizeof(Uint32)*mesh.data.nbPoint,JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }
//
    for(int i = 0; i < mesh.framecount; i++) {
    for(int j = 0; j < mesh.data.nbPoint; j++) {
        memcpy(&mesh.rdata[i].vnorm[j],&stream[offset],sizeof(Uint32));
        offset += sizeof(Uint32);
    }
    }
//
//
    for(int i = 0; i < mesh.data.nbPolygon; i++) {
        mesh.data.pltbl[i].norm[0] = -ANORMS[mesh.rdata[0].fnorm[i]][0];
        mesh.data.pltbl[i].norm[1] = -ANORMS[mesh.rdata[0].fnorm[i]][1];
        mesh.data.pltbl[i].norm[2] = -ANORMS[mesh.rdata[0].fnorm[i]][2];

    }

    for(int i = 0; i < mesh.data.nbPoint; i++) {
        mesh.data.vntbl[i][0] = -ANORMS[mesh.rdata[0].vnorm[i]][0];
        mesh.data.vntbl[i][1] = -ANORMS[mesh.rdata[0].vnorm[i]][1];
        mesh.data.vntbl[i][2] = -ANORMS[mesh.rdata[0].vnorm[i]][2];
    }

    //sets attributes of faces (transparency, mesh effect, etc)
    mesh.data.attbl = (ATTR *)jo_malloc_with_behaviour(mesh.data.nbPolygon * sizeof(*mesh.data.attbl), JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    for(int i = 0; i < mesh.data.nbPolygon; i++)
    {
        //jo_color coloroid = JO_COLOR_SATURN_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]);
        //if(i != 0) {
        ATTR proxy_attribute = ATTRIBUTE(Single_Plane,SORT_CEN, texids[i]+textureoffset, No_Palet, 0, CL32KRGB|MESHoff, sprNoflip, No_Option);
        ATTR pro_attribute = ATTRIBUTE(Single_Plane,SORT_CEN,No_Texture,C_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]),0,MESHoff,sprPolygon,No_Option);
        
        if(gouraud == true)proxy_attribute.atrb |= CL_Gouraud; 
        if(gouraud == true)proxy_attribute.sort |= UseGouraud; else proxy_attribute.sort |= No_Gouraud;
        if(gouraud == true)proxy_attribute.gstb = GRreal_ptr++; else proxy_attribute.gstb = NULL;
        if(gouraud == true)pro_attribute.sort |= UseGouraud; else pro_attribute.sort |= No_Gouraud;
        if(gouraud == true)pro_attribute.gstb = GRreal_ptr++; else pro_attribute.gstb = NULL;
        if(gouraud == true)pro_attribute.atrb |= CL_Gouraud;
        if(texids[i] != -1) {
        mesh.data.attbl[i] = proxy_attribute;
        } else {
        mesh.data.attbl[i] = pro_attribute;
        }
        //} else {
        //ATTR proxy_attribute = ATTRIBUTE(Dual_Plane,SORT_MAX, texids[i], No_Palet, 0, CL32KRGB|MESHoff, sprNoflip, No_Option);
        //ATTR pro_attribute = ATTRIBUTE(Dual_Plane,SORT_MAX,No_Texture,C_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]),0,MESHoff,sprPolygon,No_Option);
        

        if(texids[i] != -1) {
        mesh.data.attbl[i] = proxy_attribute;
        } else {
        mesh.data.attbl[i] = pro_attribute;
        }
        
    }    
    if(staticd == true) jo_free(mesh.rdata->pntbl);
    return mesh;
}

float fframe =0;
int b = 0;
ANGLE rotation;
bool jump;

void                my_gamepad(void)
{
    switch (jo_get_input_direction_pressed(0))
    {
    case LEFT: 
    sonic.x -= toFIXED(5); 
    break;
    case RIGHT: 
    sonic.x += toFIXED(5); 
    break;
    case UP: 
    sonic.z += toFIXED(5); 
    break;
    case DOWN: 
    sonic.z -= toFIXED(5); 
    break;
    case UP_LEFT:  
    sonic.x -= toFIXED(5); 
    sonic.z += toFIXED(5);  
    break;
    case UP_RIGHT:  
    sonic.x += toFIXED(5); 
    sonic.z += toFIXED(5); 
    break;
    case DOWN_LEFT: 
    sonic.x -= toFIXED(5); 
    sonic.z -= toFIXED(5); 
    break;
    case DOWN_RIGHT:  
    sonic.x += toFIXED(5); 
    sonic.z -= toFIXED(5); 
    break;
    case NONE: break;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_L)) {
        rotation -= 400;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_R)) {
        rotation += 400;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_X))
        zpitch += toFIXED(3);
    if (jo_is_input_key_pressed(0, JO_KEY_Y))
        zroll += toFIXED(3);
    if (jo_is_input_key_pressed(0, JO_KEY_Z))
        zyaw += toFIXED(3);

    if (jo_is_input_key_down(0, JO_KEY_A))
        sonic.y -= toFIXED(1);
    if (jo_is_input_key_down(0, JO_KEY_B))
        sonic.y += toFIXED(1);

    if(jump) {
        sonic.y += 1;
    }
}

vec3ang    dirtoeuler(vec3 *direction) {
    vec3ang rot;
    rot.x = -slAtan(-direction->y, direction->z);
    if (-direction->y >= 0) {
       rot.z = slAtan(slMulFX(direction->x, slCos(rot.x)), direction->z);
    }else{
       rot.z = slAtan(slMulFX(direction->x, slCos(rot.x)), -direction->z );
    }
    return rot;
}

void			    my_draw(void)
{
    my_gamepad();
 
    Plane plna = PlaneFromTri(testtri);
    cam.viewpoint_pos.y = toFIXED(-35);
    vec3 campos;
    vec3 target;
    target.x = sonic.x;
    target.y = sonic.y;
    target.z = sonic.z;
    vec3orbit(&campos,&target,toFIXED(60),rotation);
    


    cam.viewpoint_pos.x = campos.x;
    cam.viewpoint_pos.z = campos.z;

    cam.target_pos.x = sonic.x;
    cam.target_pos.y = sonic.y;
    cam.target_pos.z = sonic.z;

    fframe += 0.25;
    if(fframe >= 80) fframe = 0;
   // jo_printf(12, 1, "*ANIM TEST*");
   ssv_update(&player,(int)fframe);

    Plane floor;
    floor.normal.asArray[0] = toFIXED(0);
    floor.normal.asArray[1] = toFIXED(-1);
    floor.normal.asArray[2] = toFIXED(0);
    floor.distance = 0;
    Sphere camsphere;
    camsphere.position.asArray[0] = sonic.x;
    camsphere.position.asArray[1] = sonic.y+toFIXED(-6);
    camsphere.position.asArray[2] = sonic.z;
    camsphere.radius = toFIXED(6);




    vec3 globalup;
    globalup.x = toFIXED(0);
    globalup.y = toFIXED(0);
    globalup.z = toFIXED(0);
    //separateAngles(floor.normal,globalup, &sonic.orientation);

    sonic.orientation = dirtoeuler(&plna.normal);


    slPrintFX((sonic.x),slLocate(0,1));
    slPrintFX((sonic.y),slLocate(0,2));
    slPrintFX((sonic.z),slLocate(0,3));
    if(TriangleSphere(&testtri, &camsphere)) jo_printf(0,4,"BOOL : True"); else jo_printf(0,4,"BOOL : False");
    b += (1);
    jo_3d_camera_look_at(&cam);
    slPushMatrix();
    {
        slTranslate(sonic.x, sonic.y+toFIXED(-12), sonic.z);
        slScale(toFIXED(1),toFIXED(1),toFIXED(1));
        slRotX(sonic.orientation.asArray[0]);
        slRotY(sonic.orientation.asArray[1]);
        slRotZ(sonic.orientation.asArray[2]);
        slPutPolygonX(&player.data, livec);
    }
    slPopMatrix();

    slPushMatrix();
    {
        slPutPolygonX(&trimesh,(VECTOR){0,0,0});
    }
    slPopMatrix();

    jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(0), 0);
        slScale(toFIXED(0.75),toFIXED(0.75),toFIXED(0.75));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
        jo_background_3d_plane_a_draw(true);
        //slPutPolygon(&player.data);
    }
    jo_3d_pop_matrix();

         jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(-200), 0);
        slScale(toFIXED(1),toFIXED(1),toFIXED(1));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
        jo_background_3d_plane_b_draw(true);
        //slPutPolygon(&player.data);
    }
    jo_3d_pop_matrix();
}

jo_palette          *my_tga_palette_handling(void)
{
    // We create a new palette for each image. It's not optimal but OK for a demo.
    jo_create_palette(&image_pal);
    return (&image_pal);
}

jo_palette          *my_tga_palette_handling2(void)
{
    // We create a new palette for each image. It's not optimal but OK for a demo.
    jo_create_palette(&image_pal2);
    return (&image_pal2);
}

void			jo_main(void)
{      

    slInitGouraud(Gour, 4096, GRrealBase, GourWork);
    slIntFunction(slGouraudTblCopy); // Copy Gour to VRAM each frame
    slSetGouraudTbl(GourTbl);

    jo_add_memory_zone((unsigned char *)LWRAM, LWRAM_HEAP_SIZE);

	jo_core_init(JO_COLOR_Blue);
    jo_memory_fragmentation();
    
    jo_3d_camera_init(&cam);

    cam.viewpoint_pos.x = 0;
    cam.viewpoint_pos.y = toFIXED(-20);
    cam.viewpoint_pos.z = 0;
    cam.target_pos.x = sonic.x;
    cam.target_pos.y = sonic.y+toFIXED(-12);
    cam.target_pos.z = sonic.z;
    timz = toFIXED(0.08f);
    //slDynamicFrame(1);
    //ssvlist0 = jo_3d_create_mesh(ssv_quad_count("CBE.SSV"));
    jo_printf(0,0,"Loading : 0 %%");
    jo_fixed_point_time();
    player = ssv_load("OBJ.SSV",0,true,false);


    //player0 = ssv_load("BAL.SSV",14,true,false);
    jo_printf(0,0,"Loading : 25 %%");
    jo_printf(0,0,"Loading : 50 %%");
    jo_printf(0,0,"Loading : 75 %%");
    jo_printf(0,0,"Loading : 100 %%");

    //est_load(ssvlist0);

   // jo_printf(1, 4, "Thing: %d", ssv_quad_count(file_name));

    livec[0] = toFIXED(-0.5f);
    livec[1] = toFIXED(0.5f);
    livec[2] = toFIXED(0);
    jo_sprite_add_tga("00", "S01.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("01", "S02.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("02", "S03.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("03", "S04.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("04", "S05.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("05", "S06.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("06", "S07.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("07", "S08.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("08", "S09.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("09", "S10.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("10", "S11.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("11", "S12.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("12", "S14.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("13", "S15.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("14", "BL0.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("15", "BL1.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("16", "BL2.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("17", "BL3.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("18","CHP.TGA",JO_COLOR_Transparent);
    testtri = tridef((vec3){toFIXED(0),toFIXED(-50),toFIXED(0)},(vec3){toFIXED(50),toFIXED(0),toFIXED(0)},(vec3){toFIXED(0),toFIXED(0),toFIXED(50)});
    trimesh = coltri2mesh(&testtri);

    jo_set_tga_palette_handling(my_tga_palette_handling);    //jo_set_tga_palette_handling(&image_pal);
    jo_img_8bits    img;

    jo_enable_background_3d_plane(JO_COLOR_Black);

    // FLOOR
    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "GRS.TGA", 0);
    jo_background_3d_plane_a_img(&img, image_pal.id, true, true);
    jo_free_img(&img);

    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "SKY.TGA", 0);
    jo_background_3d_plane_b_img(&img, image_pal.id, true, true);
    jo_free_img(&img);

    //jo_core_add_callback(my_gamepad);
    jo_core_add_vblank_callback(slGouraudTblCopy);
	jo_core_add_callback(my_draw);
	jo_core_run();
}
