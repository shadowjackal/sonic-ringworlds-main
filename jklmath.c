#include <jo/jo.h>
#include "jklmath.h"

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

void                vec3orbit(vec3 *position, vec3 *target, FIXED distance, ANGLE angle) {
        position->x = target->x + slMulFX(slCos(angle),distance);
        position->z = target->z + slMulFX(slSin(angle),distance);
}

vec3 shift8(vec3 *tbs) {
    vec3 ret;
    ret.x = tbs->x>>8;
    ret.y = tbs->y>>8;
    ret.z = tbs->z>>8;
    return ret;
}

vec3ang    dirtoeuler(vec3 *direction) {
    vec3ang rot;
    //rot.x = -slAtan(-direction->y, direction->z);
    //if (-direction->y >= 0) {
    //   rot.z = slAtan(slMulFX(direction->x, slCos(rot.x)), direction->z);
    //}else{
    //   rot.z = slAtan(slMulFX(direction->x, slCos(rot.x)), -direction->z );
    //}
    rot.z=slAtan(-direction->y, direction->x);
    rot.x=-slAtan(-direction->y, direction->z);
    return rot;
}
