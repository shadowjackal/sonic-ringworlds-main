#ifndef _H_MATH_VECTORS_
#define _H_MATH_VECTORS_
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

bool cmp(FIXED x, FIXED y);

typedef struct vec2 
{
  union 
  {
    struct 
    {
      FIXED x;
      FIXED y;
    };
  FIXED asArray[2];
  };
} vec2;

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

vec2 vec2add(vec2 *left, vec2 *right);
vec3 vec3add(vec3 *left, vec3 *right);
vec2 vec2sub(vec2 *left, vec2 *right);
vec3 vec3sub(vec3 *left, vec3 *right);
vec2 vec2mul(vec2 *left, vec2 *right);
vec3 vec3mul(vec3 *left, vec3 *right);
vec2 vec2flt(vec2 *left, FIXED right);
vec3 vec3flt(vec3 *left, FIXED right);
vec3 vec3div(vec3 *left, FIXED flt);

bool vec2cmp(vec2 *left, vec2 *right);
bool vec3cmp(vec3 *left, vec3 *right);
bool vec2ncmp(vec2 *left, vec2 *right);
bool vec3ncmp(vec3 *left, vec3 *right);

FIXED dotvec2(vec2 *left, vec2 *right);
FIXED dotvec3(vec3 *left, vec3 *right);
FIXED dotvec3ult(vec3 *left, vec3 *right);


FIXED magvec2(vec2 *v);
FIXED magvec3(vec3 *v);

FIXED magsqvec2(vec2 *v);
FIXED magsqvec3(vec3 *v);

FIXED magsqvec3ult(vec3 *v);

FIXED distvec2(vec2 *p1, vec2 *p2);
FIXED distvec3(vec3 *p1, vec3 *p2);

void normalizevec2(vec2 *v);
void normalizevec3(vec3 *v);
vec2 normalizedvec2(vec2 *v);
vec3 normalizedvec3(vec3 *v);

vec3 cross(vec3 *left, vec3 *right);
vec3 crossult(vec3 *left, vec3 *right);


FIXED angvec2(vec2 *left, vec2 *right);
FIXED angvec3(vec3 *left, vec3 *right);


vec2 projvec2(vec2 *length, vec2 *direction);
vec3 projvec3(vec3 *length, vec3 *direction);

vec2 perpenvec2(vec2 *len, vec2 *dir);
vec3 perpenvec3(vec3 *len, vec3 *dir);

vec2 reflectvec2(vec2 *vec, vec2 *normal);
vec3 reflectvec3(vec3 *vec, vec3 *normal);
// Method declarations
#endif