#include <jo/jo.h>
#include "VecMath.h"


bool cmp(FIXED x, FIXED y) {
    return (ABS(x - y) <= slMulFX(JO_FIXED_EPSILON, MAX(JO_FIXED_1,MAX(ABS(x), ABS(y)))));
}

vec2 vec2add(vec2 *left, vec2 *right) {
    vec2 newvec2 = { left->x + right->x, left->y + right->y };
    return newvec2;
}

vec3 vec3add(vec3 *left, vec3 *right) {
    vec3 newvec3 = { left->x + right->x, left->y + right->y, left->z + right->z};
    return newvec3;
}

vec2 vec2sub(vec2 *left, vec2 *right) {
    vec2 newvec2 = { left->x - right->x, left->y - right->y };
    return newvec2;
}

vec3 vec3sub(vec3 *left, vec3 *right) {
    vec3 newvec3 = { left->x - right->x, left->y - right->y, left->z - right->z};
    return newvec3;
}

vec2 vec2mul(vec2 *left, vec2 *right) {
    vec2 newvec2 = {toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->x)), toFIXED(jo_fixed2int(left->y) * jo_fixed2int(right->y))};
    return newvec2;
}

vec3 vec3mul(vec3 *left, vec3 *right) {
    vec3 newvec3 = {toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->x)), toFIXED(jo_fixed2int(left->y) * jo_fixed2int(right->y)), toFIXED(jo_fixed2int(left->z) * jo_fixed2int(right->z))};
    return newvec3;
}

vec2 vec2flt(vec2 *left, FIXED flt) {
    vec2 newvec2 = {toFIXED(jo_fixed2int(left->x) * jo_fixed2int(flt)), toFIXED(jo_fixed2int(left->y) * jo_fixed2int(flt))};
    return newvec2;
}

vec3 vec3flt(vec3 *left, FIXED flt) {
    vec3 newvec3 = {slMulFX(left->x, flt), slMulFX(left->y, flt), slMulFX(left->z, flt)};
    return newvec3;
}

vec3 vec3div(vec3 *left, FIXED flt) {
    vec3 newvec3 = {toFIXED(jo_fixed2int(left->x) / jo_fixed2int(flt)), toFIXED(jo_fixed2int(left->y) / jo_fixed2int(flt)), toFIXED(jo_fixed2int(left->z) / jo_fixed2int(flt))};
    return newvec3;
}

bool vec2cmp(vec2 *left, vec2 *right) {
    return cmp(left->x, right->x) && cmp(left->y,right->y);
}

bool vec3cmp(vec3 *left, vec3 *right) {
    return cmp(left->x, right->x) && cmp(left->y,right->y) && cmp(left->z,right->z);
}

FIXED dotvec2(vec2 *left, vec2 *right) {
return  toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->x) + jo_fixed2int(left->y) * jo_fixed2int(right->y));
}

FIXED dotvec3(vec3 *left, vec3 *right) {
return toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->x) + jo_fixed2int(left->y) * jo_fixed2int(right->y) + jo_fixed2int(left->z) * jo_fixed2int(right->z));
}

FIXED dotvec3ult(vec3 *left, vec3 *right) {
return slMulFX(left->x, right->x) + slMulFX(left->y, right->y) + slMulFX(left->z, right->z);
}

FIXED magsqvec2(vec2 *v) {
return jo_int2fixed(slSquart(jo_fixed2int(dotvec2(v, v))));
}

FIXED magsqvec3(vec3 *v) {
return jo_int2fixed(slSquart(jo_fixed2int(dotvec3(v, v))));
}

FIXED magsqvec3ult(vec3 *v) {
return (slSquartFX(dotvec3ult(v, v)));
}

FIXED magvec2(vec2 *v) {
return dotvec2(v, v);
}

FIXED magvec3(vec3 *v) {
return dotvec3(v, v);
}

FIXED distvec2(vec2 *p1, vec2 *p2) {
vec2 t = vec2sub(p1, p2);
return magvec2(&t);
}

FIXED distvec3(vec3 *p1, vec3 *p2) {
vec3 t = vec3sub(p1, p2);
return magvec3(&t);
}

void normalizevec2(vec2 *v) {
*v = vec2flt(v, slDivFX(JO_FIXED_1, magvec2(v)));
}
void normalizevec3(vec3 *v) {
*v = vec3flt(v, slDivFX(JO_FIXED_1, magvec3(v)));
}
vec2 normalizedvec2(vec2 *v) {
return vec2flt(v, slDivFX(JO_FIXED_1, magvec2(v)));
}
vec3 normalizedvec3(vec3 *v) {
return vec3flt(v, slDivFX(JO_FIXED_1, magvec3(v)));
}

vec3 cross(vec3 *left, vec3 *right) {
vec3 result;
result.x = toFIXED(jo_fixed2int(left->y) * jo_fixed2int(right->z) - jo_fixed2int(left->z) * jo_fixed2int(right->y));
result.y = toFIXED(jo_fixed2int(left->z) * jo_fixed2int(right->x) - jo_fixed2int(left->x) * jo_fixed2int(right->z));
result.z = toFIXED(jo_fixed2int(left->x) * jo_fixed2int(right->y) - jo_fixed2int(left->y) * jo_fixed2int(right->x));
return result; // Done
}

vec3 crossult(vec3 *left, vec3 *right) {
vec3 result;
result.x = slMulFX(left->y, right->z) - slMulFX(left->z, right->y);
result.y = slMulFX(left->z, right->x) - slMulFX(left->x, right->z);
result.z = slMulFX(left->x, right->y) - slMulFX(left->y, right->x);
return result; // Done
}

FIXED angvec2(vec2 *left, vec2 *right) {
int m = slSquart(jo_fixed2int(magsqvec2(left)) * jo_fixed2int(magsqvec2(right)));
return toFIXED(jo_acos_radf(jo_fixed2int(dotvec2(left, right)) / m));
}

FIXED angvec3(vec3 *left, vec3 *right) {
int m = slSquart(jo_fixed2int(magsqvec3(left)) * jo_fixed2int(magsqvec3(right)));
return toFIXED(jo_acos_radf(jo_fixed2int(dotvec3(left, right)) / m));
}

vec2 projvec2(vec2 *length, vec2 *direction) {
FIXED dot = dotvec2(length, direction);
FIXED magSq = magsqvec2(direction);
return vec2flt(direction, slDivFX(dot, magSq));
}

vec3 projvec3(vec3 *length, vec3 *direction) {
FIXED dot = dotvec3(length, direction);
FIXED magSq = magsqvec3(direction);
return vec3flt(direction, slDivFX(dot, magSq));
}

vec2 perpenvec2(vec2 *len, vec2 *dir) {
    vec2 proj = projvec2(len, dir);
    return vec2sub(len, &proj);
}

vec3 perpenvec3(vec3 *len, vec3 *dir) {
    vec3 proj = projvec3(len, dir);
return vec3sub(len, &proj);
}

vec2 reflectvec2(vec2 *vec, vec2 *normal) {
FIXED d = dotvec2(vec, normal);
vec2 othervec = vec2flt(normal, slMulFX(d, JO_FIXED_2 ));
return vec2sub(vec, &othervec);
}

vec3 reflectvec3(vec3 *vec,  vec3 *normal) {
FIXED d = dotvec3(vec, normal);
vec3 othervec = vec3flt(normal, slMulFX(d, JO_FIXED_2));
return vec3sub(vec, &othervec);
}