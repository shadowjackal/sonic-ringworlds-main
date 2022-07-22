#ifndef JKLMATH_H
#define JKLMATH_H

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

extern vec3 vec3add(vec3 *left, vec3 *right);
extern vec3 vec3addn(vec3 left, vec3 right);
extern vec3 vec3sub(vec3 *left, vec3 *right);
extern vec3 vec3subn(vec3 left, vec3 right);
extern vec3 vec3mul(vec3 *left, vec3 *right);
extern vec3 vec3flt(vec3 *left, FIXED flt);
extern FIXED dotvec3(vec3 *left, vec3 *right);
extern FIXED dotvec3n(vec3 left, vec3 right);
extern int dotvec3int(vec3 *left, vec3 *right);
extern float dotvec3flt(vec3 *left, vec3 *right);

extern FIXED magsqvec3(vec3 *v);
extern FIXED magsqvec3n(vec3 v);
extern int magsqint(vec3 *v);
extern vec3 normalize(vec3 p);
extern vec3	cross_fixed(vec3 vector1, vec3 vector2);

extern void vec3orbit(vec3 *position, vec3 *target, FIXED distance, ANGLE angle);
extern vec3ang dirtoeuler(vec3 *direction);


#endif 