#include <jo/jo.h>
#include "ssv.h"
#include "newmath.h"
#include "newcol.h"

int lasthitid = -1;
int tobecheck[20];
int hitcnt = 0;

Uint8	solve_domain(FIXED normal[XYZ]){
	if(normal[X] >= 0 && normal[Z] >= 0){
		//PP
		return 1;
	} else if(normal[X] >= 0 && normal[Z] < 0){
		//PN
		return 2;
	} else if(normal[X] < 0 && normal[Z] >= 0){
		//NP
		return 3;
	} else if(normal[X] < 0 && normal[Z] < 0){
		//NN
		return 4;
	};
	/*
	3	-	1
	4	-	2
	*/
	return 0;
}

__jo_force_inline FIXED		fxm(FIXED d1, FIXED d2) //Fixed Point Multiplication
{
	register volatile FIXED rtval;
	asm(
	"dmuls.l %[d1],%[d2];"
	"sts MACH,r1;"		// Store system register [sts] , high of 64-bit register MAC to r1
	"sts MACL,%[out];"	// Low of 64-bit register MAC to the register of output param "out"
	"xtrct r1,%[out];" 	//This whole procress gets the middle 32-bits of 32 * 32 -> (2x32 bit registers)
    :    [out] "=r" (rtval)       		 //OUT
    :    [d1] "r" (d1), [d2] "r" (d2)    //IN
	:		"r1"						//CLOBBERS
	);
	return rtval;
}

//Note: this function loves to crash the system, dunno why
__jo_force_inline FIXED	fxdot(FIXED * ptA, FIXED * ptB) //Fixed-point dot product
{
	register volatile FIXED rtval;
	asm(
		"clrmac;"
		"mac.l @%[ptr1]+,@%[ptr2]+;"
		"mac.l @%[ptr1]+,@%[ptr2]+;"
		"mac.l @%[ptr1]+,@%[ptr2]+;"
		"sts MACH,r1;"
		"sts MACL,%[ox];"
		"xtrct r1,%[ox];"
		: 	[ox] "=r" (rtval)											//OUT
		:	[ptr1] "r" (ptA) , [ptr2] "r" (ptB)							//IN
		:	"r1"														//CLOBBERS
	);
	return rtval;
}



static inline FIXED dotvec3n(VECTOR vec1, VECTOR vec2) {
    VECTOR tmp1 = { vec1[X] >> 4, vec1[Y] >> 4, vec1[Z] >> 4 };
    VECTOR tmp2 = { vec2[X] >> 4, vec2[Y] >> 4, vec2[Z] >> 4 };
    
    return slInnerProduct(tmp1, tmp2);
}

void	separateAngles(FIXED unitA[XYZ], FIXED plUN[XYZ], int degreeOut[XYZ])
{	
	FIXED crossXZ[XYZ];
	FIXED crossYZ[XYZ];
	FIXED crossYX[XYZ];
	FIXED mcXZ;
	FIXED mcYZ;
	FIXED mcYX;
	FIXED dotXZ;
	FIXED dotYZ;
	FIXED dotYX;
	FIXED smuXZ[XYZ]; 
	FIXED smuYZ[XYZ]; 
	FIXED smuYX[XYZ]; 
	FIXED spnXZ[XYZ]; 
	FIXED spnYZ[XYZ]; 
	FIXED spnYX[XYZ];
	smuXZ[X] = unitA[X]; smuXZ[Y] = 0; smuXZ[Z] = unitA[Z];

	smuYZ[X] = 0; smuYZ[Y] = unitA[Y]; smuYZ[Z] = unitA[Z];

	smuYX[X] = unitA[X]; smuYX[Y] = unitA[Y]; smuYX[Z] = 0;

	spnXZ[X] = plUN[X]; spnXZ[Y] = 0; spnXZ[Z] = plUN[Z];
	
	spnYZ[X] = 0; spnYZ[Y] = plUN[Y]; spnYZ[Z] = plUN[Z];
	
	spnYX[X] = plUN[X]; spnYX[Y] = plUN[Y]; spnYX[Z] = 0;

	vec_cross(smuXZ, spnXZ, crossXZ);
	vec_cross(smuYZ, spnYZ, crossYZ);
	vec_cross(smuYX, spnYX, crossYX);
	/////////////////////////////////////////
	// Efficient (original) version
	// These appear to be functionally interchangeable, but I will make a note:
	// the vector cross-products (XZ, YZ, and YX) won't be normalized, because cross-products of normals are less than normals.
	/////////////////////////////////////////
	mcXZ = slSquartFX(fxdot(crossXZ, crossXZ));
	mcYZ = slSquartFX(fxdot(crossYZ, crossYZ));
	mcYX = slSquartFX(fxdot(crossYX, crossYX));
	/////////////////////////////////////////
	// Mathematically correct (new) version
	/////////////////////////////////////////
	// int tDot = fxdot(crossXZ, crossXZ);
	// mcXZ = slMulFX(fxdiv(1<<16, slSquartFX(tDot)), tDot);
	// tDot = fxdot(crossYZ, crossYZ);
	// mcYZ = slMulFX(fxdiv(1<<16, slSquartFX(tDot)), tDot);
	// tDot = fxdot(crossYX, crossYX);
	// mcYX = slMulFX(fxdiv(1<<16, slSquartFX(tDot)), tDot);
	//////////////////////////////////////////
	dotXZ = fxm(unitA[X],plUN[X]) + fxm(unitA[Z],plUN[Z]);
	dotYZ = fxm(unitA[Y],plUN[Y]) + fxm(unitA[Z],plUN[Z]);
	dotYX = fxm(unitA[Y],plUN[Y]) + fxm(unitA[X],plUN[X]);
	if(mcXZ != 0 || dotXZ != 0) degreeOut[X] = (slAtan(mcXZ, dotXZ));
	if(mcYZ != 0 || dotYZ != 0) degreeOut[Y] = (slAtan(mcYZ, dotYZ));
	if(mcYX != 0 || dotYX != 0) degreeOut[Z] = (slAtan(mcYX, dotYX));
	//This is currently producing angle numbers that appear like this:
	//Degrees are relative to the plane normal as if a plane. 0 degrees is parallel to the plane, but perpendicular to the planar normal direction.
	//Degrees 360 to 270 appear to be for angles facing away from the plane.
	//Degrees 0 to 90 appear to be for angles facing into the plane.
	//The sign (+/-) of the angle does not appear to be spoken of.
	
	//Three components: Y rot [X-Z]. Xrot [Y-Z]. Zrot [Y-X]. The exact definition of these rotations depends on the axis.
}

//What is this doing?
//It is re-processing the X and Z values of the output as if it were rotated X and Z _after_ it is rotated by output's Y.
//This is very strange, when I think about it.
void	sort_angle_to_domain(FIXED unitNormal[XYZ], FIXED unitOrient[XYZ], int output[XYZ])
{
static int angleComponents[XYZ];
separateAngles(unitOrient, unitNormal, angleComponents);
Uint8 domain = solve_domain(unitNormal);
//jo_printf(0, 20, "(%i)", domain);
// deg * 182 = angle
if(domain == 1){ //++
output[X] = (fxm(fxm(slSin(angleComponents[Z]), (angleComponents[Z] - 49140) ), slSin((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), (angleComponents[Y] - 49140) ), slSin((output[Y] + (90 * 182))) )); 
output[Z] = (fxm(fxm(slSin(angleComponents[Z]), angleComponents[Z] - 49140), slCos((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), angleComponents[Y] - 49140), slCos((output[Y] + (90 * 182))) )); 
						//return;
} else if(domain == 2){ //-+
output[X] = (fxm(fxm(slSin(angleComponents[Z]), (angleComponents[Z] - 49140) ), slSin((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), (angleComponents[Y] - 49140) ), slSin((output[Y] - (90 * 182))) )); 
output[Z] = (fxm(fxm(slSin(angleComponents[Z]), angleComponents[Z] - 49140), slCos((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), angleComponents[Y] - 49140), slCos((output[Y] - (90 * 182))) )); 
						//return;
} else if(domain == 3){ //+-
output[X] = -(fxm(fxm(slSin(angleComponents[Z]), (angleComponents[Z] - 49140) ), slSin((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), (angleComponents[Y] - 49140) ), slSin((output[Y] - (90 * 182))) )); 
output[Z] = -(fxm(fxm(slSin(angleComponents[Z]), angleComponents[Z] - 49140), slCos((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), angleComponents[Y] - 49140), slCos((output[Y] - (90 * 182))) )); 
						//return;
} else if(domain == 4){ //--
output[X] = -(fxm(fxm(slSin(angleComponents[Z]), (angleComponents[Z] - 49140) ), slSin((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), (angleComponents[Y] - 49140) ), slSin((output[Y] + (90 * 182))) )); 
output[Z] = -(fxm(fxm(slSin(angleComponents[Z]), angleComponents[Z] - 49140), slCos((output[Y]) - angleComponents[X])) +
						fxm(fxm(slSin(angleComponents[Y]), angleComponents[Y] - 49140), slCos((output[Y] + (90 * 182))) )); 
						//return;
}	
}

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


void   mesh2coltri(jklmesh *m, Collision *out) {
    out->points = jo_malloc((sizeof(POINT) * 3) * m->data.nbPolygon*2);
    out->num_tris = m->data.nbPolygon*2;
    int offset = 0;
    for(int i = 0; i < m->data.nbPolygon; i++) {
    tridef(m->data.pntbl[m->data.pltbl[i].Vertices[0]],m->data.pntbl[m->data.pltbl[i].Vertices[1]],m->data.pntbl[m->data.pltbl[i].Vertices[2]],&out->points[offset]);
    tridef(m->data.pntbl[m->data.pltbl[i].Vertices[2]],m->data.pntbl[m->data.pltbl[i].Vertices[3]],m->data.pntbl[m->data.pltbl[i].Vertices[0]],&out->points[offset+3]);
    offset += 6;
    }
}

static inline void plane_point_closest_pt(POINT *triangle, VECTOR P, VECTOR out) {
    VECTOR plane_normal;
    slNormalVector(triangle[0], triangle[1], triangle[2], plane_normal);
    
    FIXED plane_distance = slInnerProduct(plane_normal, triangle[0]);
    FIXED point_plane_dist = slInnerProduct(plane_normal, P) - plane_distance;
    
    vec_mul_scalar(plane_normal, point_plane_dist, out);
    vec_sub(P, out, out);
}

static inline void plane_point_closest_pt_plane(Plane *plane, VECTOR P, VECTOR out) {
    

    FIXED point_plane_dist = slInnerProduct(plane->normal, P) - plane->distance;
    
    vec_mul_scalar(plane->normal, point_plane_dist, out);
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
        if(manhattandistance(sphere->pos,col->points[i*3]) < toFIXED(150)) {
        VECTOR closest_pt;
        tri_point_closest_pt(sphere->pos, &col->points[i*3], closest_pt);
        
        VECTOR delta;
        vec_sub(sphere->pos, closest_pt, delta);
        
        FIXED delta_dist = slSquartFX(dotvec3n(delta, delta));
        
        if(delta_dist < (sphere->radius >> 4)) {
            if(is_hit != true) hitcnt = 0;
            vec_normalize(delta);
            
            VECTOR push_back;
            vec_mul_scalar(delta, sphere->radius - (delta_dist << 4), push_back);
            
            vec_add(sphere->pos, push_back, sphere->pos);
            tobecheck[hitcnt] = i;
            is_hit = true;
            hitcnt += 1;
        }}

    }
    
    if(is_hit == false) lasthitid = -1;
    return is_hit;
}

bool Collision_SpherePlaneResolve(Sphere *sphere, Plane *plane) {
    bool is_hit = false;
    VECTOR closest_pt;
    plane_point_closest_pt_plane(plane,sphere->pos,closest_pt);

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
        
    return false;
}

FIXED orient[3] = {0, -1<<16, 0};

bool Collision_SphereCol_bool(Sphere *sphere, Collision *col) {
    bool is_hit = false;
    
    // Resolve collision against each triangle
    for(int i = 0; i < col->num_tris; i++) {
        if(manhattandistance(sphere->pos,col->points[i]) < toFIXED(300)) {
        // Get closest point on tri to the sphere centre
        VECTOR closest_pt;
        tri_point_closest_pt(sphere->pos, &col->points[i], closest_pt);
        
        VECTOR delta;
        vec_sub(sphere->pos, closest_pt, delta);
        
        FIXED delta_dist = slSquartFX(dotvec3n(delta, delta));
        
        if(delta_dist < (sphere->radius >> 4)) {
            is_hit = true;
        }
        if(is_hit == true) return is_hit;
    }}
    
    return is_hit;
}


bool Collision_SphereCol_bool_special(Sphere *sphere, Collision *col) {
    bool is_hit = false;
    // Resolve collision against each triangle
    for(int i = 0; i < hitcnt; i++) {
        if(manhattandistance(sphere->pos,col->points[tobecheck[i]*3]) < toFIXED(300)) {
        // Get closest point on tri to the sphere centre
        VECTOR closest_pt;
        tri_point_closest_pt(sphere->pos, &col->points[tobecheck[i]*3], closest_pt);
        
        VECTOR delta;
        vec_sub(sphere->pos, closest_pt, delta);
        
        FIXED delta_dist = slSquartFX(dotvec3n(delta, delta));
        
        if(delta_dist < (sphere->radius >> 4)) {
            is_hit = true;
            FIXED normies[3];

            VECTOR U;
            vec_sub(col->points[(tobecheck[i]+1)*3],col->points[tobecheck[i]*3],U);
            VECTOR V;
            vec_sub(col->points[(tobecheck[i]+2)*3],col->points[tobecheck[i]*3],V);

            normies[X] = (slMulFX(U[Y],V[Z])) - (slMulFX(U[Z],V[Y]));
            normies[Y] = (slMulFX(U[Z],V[X])) - (slMulFX(U[X],V[Z]));
            normies[Z] = (slMulFX(U[X],V[Y])) - (slMulFX(U[Y],V[X]));

            vec_normalize(normies);

            slPrintFX(toFIXED(normies[X]), slLocate(20,4));
            slPrintFX(toFIXED(normies[Y]), slLocate(20,5));
            slPrintFX(toFIXED(normies[Z]), slLocate(20,6));


            sonic.orientation[X] = LerpAng(sonic.orientation[X], slAtan(JO_ABS(normies[Y]), normies[Z]),toFIXED(0.25));
            sonic.orientation[Z] = LerpAng(sonic.orientation[Z], slAtan(JO_ABS(normies[Y]), normies[X]),toFIXED(0.25));

            return is_hit;
        }
    }}
    
    return is_hit;
}

bool Collision_SpherePlane_bool(Sphere *sphere, Plane *plane) {
    bool is_hit = false;
    VECTOR closest_pt;
    plane_point_closest_pt_plane(plane,sphere->pos,closest_pt);

    VECTOR delta;
    vec_sub(sphere->pos, closest_pt, delta);

    FIXED delta_dist = slSquartFX(dotvec3n(delta, delta));

            if(delta_dist < (sphere->radius >> 4)) {
            
            is_hit = true;
        }
        
    return false;
}