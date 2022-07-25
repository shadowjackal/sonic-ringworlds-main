#include <jo/jo.h>
#include "newmath.h"
#include "newcol.h"

PDATA   coltri2mesh(POINT *t) {
    PDATA retdata;
    retdata.nbPoint = 4;
    retdata.nbPolygon = 1;
    retdata.pntbl = jo_malloc(sizeof(POINT)*retdata.nbPoint);
    retdata.pltbl = jo_malloc(sizeof(POLYGON)*retdata.nbPolygon);
    retdata.attbl = jo_malloc(sizeof(ATTR) * retdata.nbPolygon);

    retdata.pntbl[0][0] = t[0][0];
    retdata.pntbl[0][1] = t[0][1];
    retdata.pntbl[0][2] = t[0][2];
    retdata.pntbl[1][0] = t[1][0];
    retdata.pntbl[1][1] = t[1][1];
    retdata.pntbl[1][2] = t[1][2];
    retdata.pntbl[2][0] = t[2][0];
    retdata.pntbl[2][1] = t[2][1];
    retdata.pntbl[2][2] = t[2][2];
    retdata.pntbl[3][0] = t[2][0];
    retdata.pntbl[3][1] = t[2][1];
    retdata.pntbl[3][2] = t[2][2];
    retdata.pltbl[0].Vertices[0] = 0;
    retdata.pltbl[0].Vertices[1] = 1;
    retdata.pltbl[0].Vertices[2] = 2;
    retdata.pltbl[0].Vertices[3] = 3;

    ATTR pro_attribute = ATTRIBUTE(Dual_Plane,SORT_CEN|No_Gouraud,No_Texture,C_RGB(31,31,0),0,MESHoff,sprPolygon,No_Option);
    retdata.attbl[0] = pro_attribute;
    return retdata;
}

void    tridef(const POINT p1, const POINT p2, const POINT p3, POINT* out) {
   // if(out[0] == NULL) &out[0] = jo_malloc(sizeof(POINT)*3);
    out[0][0] = p1[0];
    out[0][1] = p1[1];
    out[0][2] = p1[2];
    out[1][0] = p2[0];
    out[1][1] = p2[1];
    out[1][2] = p2[2];
    out[2][0] = p3[0];
    out[2][1] = p3[1];
    out[2][2] = p3[2];
}


static inline void plane_point_closest_pt(POINT *triangle, VECTOR P, VECTOR out) {
    VECTOR plane_normal;
    slNormalVector(triangle[0], triangle[1], triangle[2], plane_normal);
    
    FIXED plane_distance = slInnerProduct(plane_normal, triangle[0]);
    FIXED point_plane_dist = slInnerProduct(plane_normal, P) - plane_distance;
    
    vec_mul_scalar(plane_normal, point_plane_dist, out);
    vec_sub(P, out, out);
}

static inline bool tri_point_intersects(POINT *triangle, VECTOR P) {
    VECTOR a, b, c;
    vec_sub(triangle[0], P, a);
    vec_sub(triangle[1], P, b);
    vec_sub(triangle[2], P, c);
    
    // Prevent fixed-point overflow
    vec_normalize(a);
    vec_normalize(b);
    vec_normalize(c);
    
    VECTOR u, v, w;
    vec_cross(b, c, u);
    vec_cross(c, a, v);
    vec_cross(a, b, w);
    
    if (slInnerProduct(u, v) < 0)
        return false;
    if (slInnerProduct(u, w) < 0)
        return false;
    
    return true;
}

static inline FIXED dotvec3n(VECTOR vec1, VECTOR vec2) {
    VECTOR tmp1 = { vec1[X] >> 4, vec1[Y] >> 4, vec1[Z] >> 4 };
    VECTOR tmp2 = { vec2[X] >> 4, vec2[Y] >> 4, vec2[Z] >> 4 };
    
    return slInnerProduct(tmp1, tmp2);
}

static inline void line_point_closest_pt(VECTOR A, VECTOR B, VECTOR P, VECTOR out) {
    VECTOR AP, AB;
    vec_sub(P, A, AP);
    vec_sub(B, A, AB);
    
    FIXED ap_dot_ab = dotvec3n(AP, AB);
    
    // Point is behind start of the segment, so perpendicular distance is not viable
    if (ap_dot_ab <= toFIXED(0.0)) {
        vec_copy(A, out);
        return;
    }
    
    FIXED ab_dot_ab = dotvec3n(AB, AB);
    
    // Point is beyond end of the segment, so perpendicular distance is not viable
    if (ap_dot_ab >= ab_dot_ab) {
        vec_copy(B, out);
        return;
    }
    
    FIXED t = slDivFX(ab_dot_ab, ap_dot_ab);
    
    vec_mul_scalar(AB, t, out);
    vec_add(out, A, out);
}

static inline void tri_point_closest_pt(POINT P, POINT *triangle, VECTOR out) {
    // Closest tri plane point test
    VECTOR plane_closest_pt;
    plane_point_closest_pt(triangle, P, plane_closest_pt);
    
    // Is tri plane the closest point?
    if(tri_point_intersects(triangle, plane_closest_pt)) {
        vec_copy(plane_closest_pt, out);
        return;
    }
    
    // Closest line point tests
    VECTOR A, B, C;
    line_point_closest_pt(triangle[0], triangle[1], P, A);
    line_point_closest_pt(triangle[0], triangle[2], P, B);
    line_point_closest_pt(triangle[1], triangle[2], P, C);
    
    // Length vectors
    VECTOR AP, BP, CP;
    vec_sub(P, A, AP);
    vec_sub(P, B, BP);
    vec_sub(P, C, CP);
    
    FIXED AP_len = dotvec3n(AP, AP);
    FIXED BP_len = dotvec3n(BP, BP);
    FIXED CP_len = dotvec3n(CP, CP);
    
    // Return the closest of the 3 line points
    if(AP_len < BP_len && AP_len < CP_len)
        vec_copy(A, out);
    else if(BP_len < CP_len)
        vec_copy(B, out);
    else
        vec_copy(C, out);
}

bool Collision_SphereColResolve(Sphere *sphere, Collision *col) {
    bool is_hit = false;
    
    // Resolve collision against each triangle
    for(int i = 0; i < col->num_tris; i++) {
        // Get closest point on tri to the sphere centre
        VECTOR closest_pt;
        tri_point_closest_pt(sphere->pos, &col->points[i*3], closest_pt);
        
        VECTOR delta;
        vec_sub(sphere->pos, closest_pt, delta);
        
        FIXED delta_dist = slSquartFX(dotvec3n(delta, delta));
        
        if(delta_dist < (sphere->radius >> 4)) {
            vec_normalize(delta);
            
            VECTOR push_back;
            vec_mul_scalar(delta, sphere->radius - (delta_dist << 4), push_back);
            
            vec_add(sphere->pos, push_back, sphere->pos);
            is_hit = true;
        }
    }
    
    return is_hit;
}