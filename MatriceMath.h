#ifndef _H_MATH_MATRICES_
#define _H_MATH_MATRICES_
#include "VecMath.h"
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Structure definitions
typedef struct mat2 {
union {
struct {
FIXED _11, _12,
_21, _22;
};
FIXED asArray[4];
};
} mat2;

typedef struct mat3 {
union {
struct {
FIXED _11, _12, _13,
_21, _22, _23,
_31, _32, _33;
};
FIXED asArray[9];
};
} mat3;

typedef struct mat4 {
union {
struct {
FIXED _11, _12, _13, _14,
_21, _22, _23, _24,
_31, _32, _33, _34,
_41, _42, _43, _44;
};
FIXED asArray[16];
};
} mat4;

void Transpose(FIXED *srcMat, FIXED *dstMat, int srcRows, int srcCols);
mat2 TransposeMat2(mat2 *matrix);
mat3 TransposeMat3(mat3 *matrix);
mat4 TransposeMat4(mat4 *matrix);

mat2 mat2flt(mat2 *matrix, FIXED scalar);
mat3 mat3flt(mat3 *matrix, FIXED scalar);
mat4 mat4flt(mat4 *matrix, FIXED scalar);

bool Multiply(FIXED* out, FIXED* matA, int aRows, int aCols, FIXED* matB, int bRows, int bCols);
mat2 mat2mul(mat2 *matA, mat2 *matB);
mat3 mat3mul(mat3 *matA, mat3 *matB);
mat4 mat4mul(mat4 *matA, mat4 *matB);

FIXED Determinant(mat2 *matrix);

void mat2default(mat2 *mat);
void mat3default(mat3 *mat);
void mat4default(mat4 *mat);

void mat2construct(mat2 *matdes, FIXED f11, FIXED f12, FIXED f21, FIXED f22);
void mat3construct(mat3 *matdes, FIXED f11, FIXED f12, FIXED f13, FIXED f21, FIXED f22, FIXED f23, FIXED f31, FIXED f32, FIXED f33);
void mat4construct(mat4 *matdes, FIXED f11, FIXED f12, FIXED f13, FIXED f14, FIXED f21, FIXED f22, FIXED f23, FIXED f24, FIXED f31, FIXED f32, FIXED f33, FIXED f34, FIXED f41, FIXED f42, FIXED f43, FIXED f44);

FIXED mat2determinant(mat2 *matrix);
FIXED mat3determinant(mat3 *mat);
FIXED mat4determinant(mat4 *mat);

mat3  RotateUnit(FIXED pitch, FIXED roll, FIXED yaw);

mat2 mat3cut(mat3 *mat, int row, int col);
mat3 mat4cut(mat4 *mat, int row, int col);

mat2 mat2minor(mat2 *mat);
mat3 mat3minor(mat3 *mat);
mat4 mat4minor(mat4 *mat);

void gencofactor(FIXED* out, FIXED* minor, int rows, int cols);
mat3 mat3cofactor(mat3 *mat);
mat2 mat2cofactor(mat2 *mat);
mat4 mat4cofactor(mat4 *mat);

mat2 mat2adjugate(mat2 *mat);
mat3 mat3adjugate(mat3 *mat);
mat4 mat4adjugate(mat4 *mat);

mat2 mat2inverse(mat2 *mat);
mat3 mat3inverse(mat3 *mat);
mat4 mat4inverse(mat4 *mat);

//stuff that does translation rotation etc
mat4 postranslation(FIXED x, FIXED y, FIXED z);
mat4 vectranslation(vec3 *pos);
vec3 gettranslation(mat4 *mat);

mat4 posscale(FIXED x, FIXED y, FIXED z);
mat4 vecscale(vec3 *vec);
vec3 getscale(mat4 *mat);

mat4 mat4zrotation(FIXED angle);
mat3 mat3zrotation(FIXED angle);
mat4 mat4xrotation(FIXED angle);
mat3 mat3xrotation(FIXED angle);
mat4 mat4yrotation(FIXED angle);
mat3 mat3yrotation(FIXED angle);

mat4 mat4rotation(FIXED pitch, FIXED yaw, FIXED roll);
mat3 mat3rotation(FIXED pitch, FIXED yaw, FIXED roll);

mat4 mat4axisangle(vec3 *axis, FIXED angle);
mat3 mat3axisangle(vec3 *axis, FIXED angle);

vec3 mat4multpoint(vec3 *vec, mat4 *mat);
vec3 mat4multvector(vec3 *vec, mat4 *mat);
vec3 mat3multvector(vec3 *vec, mat3 *mat);

mat4 transformeuler(vec3 *scale, vec3 *eulerRotation, vec3 *translate);
mat4 transformaxis(vec3 *scale, vec3 *rotationAxis,FIXED rotationAngle, vec3 *translate);

mat4 lookatpoint(vec3 *position, vec3 *target, vec3 *up);

mat4 Projection(FIXED fov, FIXED aspect, FIXED zNear, FIXED zFar);
mat4 Ortho(FIXED left, FIXED right, FIXED bottom, FIXED top, FIXED zNear, FIXED zFar);
#endif