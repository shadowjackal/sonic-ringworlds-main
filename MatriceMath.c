#include <jo/jo.h>
#include "MatriceMath.h"

void Transpose(FIXED *srcMat, FIXED *dstMat, int srcRows, int srcCols) {
for (int i = 0; i < srcRows * srcCols; i++) {
int row = i / srcRows;
int col = i % srcRows;
dstMat[i] = srcMat[srcCols * col + row];
}
}

mat2 TransposeMat2(mat2 *matrix) {
mat2 result;
Transpose(matrix->asArray, result.asArray, 2, 2);
return result;
}

mat3 TransposeMat3(mat3 *matrix) {
mat3 result;
Transpose(matrix->asArray, result.asArray, 3, 3);
return result;
}

mat4 TransposeMat4(mat4 *matrix) {
mat4 result;
Transpose(matrix->asArray, result.asArray, 4, 4);
return result;
}

mat2 mat2flt(mat2 *matrix, FIXED scalar) {
mat2 result;
for (int i = 0; i < 4; ++i) {
result.asArray[i] = slMulFX(matrix->asArray[i], scalar);
}
return result;
}

mat3 mat3flt(mat3 *matrix, FIXED scalar) {
mat3 result;
for (int i = 0; i < 9; ++i) {
result.asArray[i] = slMulFX(matrix->asArray[i], scalar);
}
return result;
}

mat4 mat4flt(mat4 *matrix, FIXED scalar) {
mat4 result;
for (int i = 0; i < 16; ++i) {
result.asArray[i] = slMulFX(matrix->asArray[i], scalar);
}
return result;
}

bool Multiply(FIXED* out, FIXED* matA, int aRows, int aCols, FIXED* matB, int bRows, int bCols) {
	if (aCols != bRows) {
		return false;
	}

	for (int i = 0; i < aRows; ++i) {
		for (int j = 0; j < bCols; ++j) {
			out[bCols * i + j] = 0;
		for (int k = 0; k < bRows; ++k) {
			int a = aCols * i + k;
			int b = bCols * k + j;
			out[bCols * i + j] += slMulFX(matA[a], matB[b]);
			}
		}
	}
return true;
}

mat2 mat2mul(mat2 *matA, mat2 *matB) {
mat2 res;
Multiply(res.asArray, matA->asArray, 2, 2, matB->asArray, 2, 2);
return res;
}

mat3 mat3mul(mat3 *matA, mat3 *matB) {
mat3 res;
Multiply(res.asArray, matA->asArray, 3, 3, matB->asArray, 3, 3);
return res;
}

mat4 mat4mul(mat4 *matA, mat4 *matB) {
mat4 res;
Multiply(res.asArray, matA->asArray, 4, 4, matB->asArray, 4, 4);
return res;
}

FIXED Determinant(mat2 *matrix) {
return slMulFX(matrix->_11, matrix->_22) - slMulFX(matrix->_12, matrix->_21);
}

void mat2default(mat2 *mat) {
mat->_11 = mat->_22 = JO_FIXED_1;
mat->_12 = mat->_21 = 0;
}

void mat3default(mat3 *mat) {
mat->_11 = mat->_22 = mat->_33 = JO_FIXED_1;
mat->_12 = mat->_13 = mat->_21 = 0;
mat->_23 = mat->_31 = mat->_32 = 0;
}

void mat4default(mat4 *mat) {
mat->_11 = mat->_22 = mat->_33 = mat->_44 = JO_FIXED_1;
mat->_12 = mat->_13 = mat->_14 = mat->_21 = 0;
mat->_23 = mat->_24 = mat->_31 = mat->_32 = 0;
mat->_34 = mat->_41 = mat->_42 = mat->_43 = 0;
}

void mat2construct(mat2 *matdes, FIXED f11, FIXED f12, FIXED f21, FIXED f22) {
matdes->_11 = f11; matdes->_12 = f12;
matdes->_21 = f21; matdes->_22 = f22;
}

void mat3construct(mat3 *matdes, FIXED f11, FIXED f12, FIXED f13, FIXED f21, FIXED f22, FIXED f23, FIXED f31, FIXED f32, FIXED f33) {
matdes->_11 = f11; matdes->_12 = f12; matdes->_13 = f13;
matdes->_21 = f21; matdes->_22 = f22; matdes->_23 = f23;
matdes->_31 = f31; matdes->_32 = f32; matdes->_33 = f33;
}

void mat4construct(mat4 *matdes, FIXED f11, FIXED f12, FIXED f13, FIXED f14, FIXED f21, FIXED f22, FIXED f23, FIXED f24, FIXED f31, FIXED f32, FIXED f33, FIXED f34, FIXED f41, FIXED f42, FIXED f43, FIXED f44) {
matdes->_11 = f11; matdes->_12 = f12; matdes->_13 = f13; matdes->_14 = f14;
matdes->_21 = f21; matdes->_22 = f22; matdes->_23 = f23; matdes->_24 = f24;
matdes->_31 = f31; matdes->_32 = f32; matdes->_33 = f33; matdes->_34 = f34;
matdes->_41 = f41; matdes->_42 = f42; matdes->_43 = f43; matdes->_44 = f44;
}

FIXED mat2determinant(mat2 *matrix) {
return slMulFX(matrix->_11, matrix->_22) - slMulFX(matrix->_12, matrix->_21);
}

mat2 matcut(mat3 *mat, int row, int col) {
mat2 result;
int index = 0;
for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
if (i == row || j == col) {
continue;
}
int target = index++;
int source = 3 * i + j;
result.asArray[target] = mat->asArray[source];
}
}
return result;
}

mat2 mat2minor(mat2 *mat) {
mat2 matnew;
mat2construct(&matnew, mat->_22, mat->_21, mat->_12, mat->_11);
return matnew;
}

mat3 mat3minor(mat3 *mat) {
mat3 result;
for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
mat2 curmat = matcut(mat, i, j);
result.asArray[3 * i + j] = mat2determinant(&curmat);
}
}
return result;
}

void gencofactor(FIXED* out, FIXED* minor, int rows, int cols) {
for (int i = 0; i < rows; ++i) {
for (int j = 0; j < cols; ++j) {
int t = cols * j + i; // Target index
int s = cols * j + i; // Source index
FIXED sign = jo_fixed_pow(-JO_FIXED_1, toFIXED(i + j)); // + or –
out[t] = slMulFX(minor[s], sign);
}}
}

mat2 mat2cofactor(mat2 *mat) {
mat2 result;
gencofactor(result.asArray, mat2minor(mat).asArray, 2, 2);
return result;
}

mat3 mat3cofactor(mat3 *mat) {
mat3 result;
gencofactor(result.asArray, mat3minor(mat).asArray, 3, 3);
return result;
}

FIXED mat3determinant(mat3 *mat) {
FIXED result = 0;
mat3 gofactor = mat3cofactor(mat);
for (int j = 0; j < 3; ++j) {
int index = 3 * 0 + j;
result += mat->asArray[index] * gofactor.asArray[j*3];
}
return result;
}

mat3 mat4cut(mat4 *mat, int row, int col) {
mat3 result;
int index = 0;
for (int i = 0; i < 4; ++i) {
for (int j = 0; j < 4; ++j) {
if (i == row || j == col) {
continue;
}
int target = index++;
int source = 4 * i + j;
result.asArray[target] = mat->asArray[source];
}
}
return result;
}

mat4 mat4minor(mat4 *mat) {
mat4 result;
for (int i = 0; i <4; ++i) {
for (int j = 0; j <4; ++j) {
mat3 curmat = mat4cut(mat, i, j);
result.asArray[4 * i + j] = mat3determinant(&curmat);
}
}
return result;
}

mat4 mat4cofactor(mat4 *mat) {
mat4 result;
gencofactor(result.asArray, mat4minor(mat).asArray, 4, 4);
return result;
}

FIXED mat4determinant(mat4 *mat) {
FIXED result = 0;
mat4 cofactor = mat4cofactor(mat);
for (int j = 0; j < 4; ++j) {
result += slMulFX(mat->asArray[4 * 0 + j], cofactor.asArray[j * 4]);
}
return result;
}

mat2 mat2adjugate(mat2 *mat) {
mat2 cofact = mat2cofactor(mat);
return TransposeMat2(&cofact);
}

mat3 mat3adjugate(mat3 *mat) {
mat3 cofact = mat3cofactor(mat);
return TransposeMat3(&cofact);
}

mat4 mat4adjugate(mat4 *mat) {
mat4 cofact = mat4cofactor(mat);
return TransposeMat4(&cofact);
}

mat2 mat2inverse(mat2 *mat) {
FIXED det = mat2determinant(mat);
if (cmp(det, 0)) { mat2 retmatdef; mat2default(&retmatdef); return retmatdef; }
mat2 matadj = mat2adjugate(mat);
return mat2flt(&matadj, slDivFX(JO_FIXED_1, det));
}

mat3 mat3inverse(mat3 *mat) {
FIXED det = mat3determinant(mat);
if (cmp(det, 0)) { mat3 retmatdef; mat3default(&retmatdef); return retmatdef; }
mat3 matadj = mat3adjugate(mat);
return mat3flt(&matadj, slDivFX(JO_FIXED_1, det));
}

mat4 mat4inverse(mat4 *mat) {
FIXED det = mat4determinant(mat);
if (cmp(det, 0)) { mat4 retmatdef; mat4default(&retmatdef); return retmatdef; }
mat4 matadj = mat4adjugate(mat);
return mat4flt(&matadj, slDivFX(JO_FIXED_1, det));
}

mat4 Translation(FIXED x, FIXED y, FIXED z) {
mat4 retmat;
 mat4construct(&retmat, JO_FIXED_1, 0, 0, 0, 0, JO_FIXED_1, 0, 0, 0, 0, JO_FIXED_1, 0, x, y, z, JO_FIXED_1
);
return retmat;
}

mat4 vectranslation(vec3 *pos) {
mat4 retmat;
 mat4construct(&retmat, JO_FIXED_1, 0, 0, 0, 0, JO_FIXED_1, 0, 0, 0, 0, JO_FIXED_1, 0, pos->x, pos->y, pos->z,JO_FIXED_1);
return retmat;
}

vec3 gettranslation(mat4 *mat) {
 vec3 retvec;
retvec.x = mat->_41;
retvec.y = mat->_42;
retvec.z = mat->_43;
return retvec;  
}

mat4 posscale(FIXED x, FIXED y, FIXED z) {
mat4 retmat;
mat4construct(&retmat, x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, JO_FIXED_1
);
return retmat;
}

mat4 vecscale(vec3 *vec) {
mat4 retmat;
mat4construct(&retmat, vec->x, 0, 0, 0, 0, vec->y,0, 0, 0, 0, vec->z,0, 0, 0, 0, JO_FIXED_1
);
return retmat;
}

vec3 getscale(mat4 *mat) {
vec3 retvec;
retvec.x = mat->_11;
retvec.y = mat->_22;
retvec.z = mat->_33;
return retvec;
}

mat4 mat4zrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat4 retmat;
mat4construct(&retmat,slCos(angle), slSin(angle), 0, 0, -slSin(angle), slCos(angle), 0, 0, 0, 0, JO_FIXED_1, 0, 0, 0, 0, JO_FIXED_1);
return retmat;
}

mat3 mat3zrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat3 retmat;
mat3construct(&retmat, slCos(angle), slSin(angle), 0, -slSin(angle), slCos(angle), 0, 0, 0, JO_FIXED_1);
return retmat;
}

mat4 mat4xrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat4 retmat;
mat4construct(&retmat, JO_FIXED_1, 0, 0, 0, 0, slCos(angle), slSin(angle), 0, 0, -slSin(angle), slCos(angle), 0, 0, 0, 0, JO_FIXED_1);
return retmat;
}

mat3 mat3xrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat3 retmat;
mat3construct(&retmat,
JO_FIXED_1, 0, 0,
0, slCos(angle), slSin(angle),
0, -slSin(angle), slCos(angle)
);
return retmat;
}

mat4 mat4yrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat4 retmat;
mat4construct(&retmat,
slCos(angle), 0, -slSin(angle), 0,
0, JO_FIXED_1, 0, 0,
slSin(angle), 0, slCos(angle), 0,
0, 0, 0, JO_FIXED_1
);
return retmat;
}

mat3 mat3yrotation(FIXED angle) {
angle = DEGtoANG(angle);
mat3 retmat;
mat3construct(&retmat,
slCos(angle), 0, -slSin(angle),
0, JO_FIXED_1, 0,
slSin(angle), 0, slCos(angle)
);
return retmat;
}


mat4 mat4rotation(FIXED pitch, FIXED yaw, FIXED roll) {
mat4 ZRot = mat4zrotation(roll);
mat4 XRot = mat4xrotation(pitch);
mat4 YRot = mat4yrotation(yaw);
mat4 ZXRot = mat4mul(&ZRot,&XRot);
return mat4mul(&ZXRot,&YRot);
}


mat3  RotateUnit(FIXED pitch, FIXED roll, FIXED yaw) {
    MATRIX tmatrix;

    jo_3d_push_matrix();
	{
        slGetMatrix(tmatrix);
        slRotX(DEGtoANG(pitch));
        slRotY(DEGtoANG(roll));
        slRotZ(DEGtoANG(yaw));
        slGetMatrix(tmatrix);
    }
    jo_3d_pop_matrix();
	mat3 retmat;
    retmat.asArray[0] = tmatrix[X][X];    
    retmat.asArray[1] = tmatrix[X][Y];
    retmat.asArray[2] = tmatrix[X][Z];    
    retmat.asArray[3] = tmatrix[Y][X];    
    retmat.asArray[4] = tmatrix[Y][Y];
    retmat.asArray[5] = tmatrix[Y][Z];    
    retmat.asArray[6] = tmatrix[Z][X];    
    retmat.asArray[7] = tmatrix[Z][Y];
    retmat.asArray[8] = tmatrix[Z][Z];    
	return retmat;
}

mat3 mat3rotation(FIXED pitch, FIXED yaw, FIXED roll) { 
    //mat3 rx = mat3xrotation(pitch); 
	//mat3 ry = mat3yrotation(yaw); 
	//mat3 rz = mat3zrotation(roll);
//
    //mat3 rot1 = mat3mul(&rx,&ry);
    //mat3 rot2 = mat3mul(&rz,&rot1);

	return RotateUnit(pitch,yaw,roll);
}


mat4 mat4axisangle(vec3 *axis, FIXED angle) {
angle = DEGtoANG(angle);
FIXED c = slCos(angle);
FIXED s = slSin(angle);
FIXED t = JO_FIXED_1 - slCos(angle);
FIXED x = axis->x;
FIXED y = axis->y;
FIXED z = axis->z;
if (!cmp(magsqvec3(axis), JO_FIXED_1)) {
FIXED inv_len = slDivFX(JO_FIXED_1, magvec3(axis));
x = slMulFX(x, inv_len); // Normalize x
y = slMulFX(y, inv_len); // Normalize y
z = slMulFX(z, inv_len); // Normalize z
} // x, y, and z are a normalized vector
mat4 retmat;
mat4construct(&retmat,
toFIXED(jo_fixed2int(t)*(jo_fixed2int(x)*jo_fixed2int(x)) + jo_fixed2int(c)),
toFIXED(jo_fixed2int(t)*jo_fixed2int(x)*jo_fixed2int(y) + jo_fixed2int(s)*jo_fixed2int(z)),
toFIXED(jo_fixed2int(t)*jo_fixed2int(x)*jo_fixed2int(z) - jo_fixed2int(s)*jo_fixed2int(y)), 
0,
toFIXED(jo_fixed2int(t)*jo_fixed2int(x)*jo_fixed2int(y) - jo_fixed2int(s)*jo_fixed2int(z)), 
toFIXED(jo_fixed2int(t)*(jo_fixed2int(y)*jo_fixed2int(y)) + jo_fixed2int(c)), 
toFIXED(jo_fixed2int(t)*jo_fixed2int(y)*jo_fixed2int(z) + jo_fixed2int(s)*jo_fixed2int(x)), 
0,
toFIXED(jo_fixed2int(t)*jo_fixed2int(x)*jo_fixed2int(z) + jo_fixed2int(s)*jo_fixed2int(y)), 
toFIXED(jo_fixed2int(t)*jo_fixed2int(y)*jo_fixed2int(z) - jo_fixed2int(s)*jo_fixed2int(x)), 
toFIXED(jo_fixed2int(t)*(jo_fixed2int(z)*jo_fixed2int(z)) + jo_fixed2int(c)), 
0,
0, 0, 0, JO_FIXED_1
);
return retmat;
}

mat3 mat3axisangle(vec3 *axis, FIXED angle) {
	angle = DEGtoANG(angle);
	FIXED c = slCos(angle);
	FIXED s = slSin(angle);
	FIXED t = JO_FIXED_1 - slCos(angle);

	FIXED x = axis->x;
	FIXED y = axis->y;
	FIXED z = axis->z;
	if (!cmp(magsqvec3(axis), JO_FIXED_1)) {
		FIXED inv_len = slDivFX(JO_FIXED_1, magsqvec3(axis));
		x = slMulFX(x ,inv_len);
		y = slMulFX(y ,inv_len);
		z = slMulFX(z ,inv_len);
	}

	return (mat3){
		toFIXED(jo_fixed2int(t) * (jo_fixed2int(x) * jo_fixed2int(x)) + jo_fixed2int(c)),
		toFIXED(jo_fixed2int(t) * jo_fixed2int(x) * jo_fixed2int(y) + jo_fixed2int(s) * jo_fixed2int(z)),
		toFIXED(jo_fixed2int(t) * jo_fixed2int(x) * jo_fixed2int(z) - jo_fixed2int(s) * jo_fixed2int(y)), 
		toFIXED(jo_fixed2int(t) * jo_fixed2int(x) * jo_fixed2int(y) - jo_fixed2int(s) * jo_fixed2int(z)),
		toFIXED(jo_fixed2int(t) * (jo_fixed2int(y) * jo_fixed2int(y)) + jo_fixed2int(c)),
		toFIXED(jo_fixed2int(t) * jo_fixed2int(y) * jo_fixed2int(z) + jo_fixed2int(s) * jo_fixed2int(x)), 
		toFIXED(jo_fixed2int(t) * jo_fixed2int(x) * jo_fixed2int(z) + jo_fixed2int(s) * jo_fixed2int(y)),
		toFIXED(jo_fixed2int(t) * jo_fixed2int(y) * jo_fixed2int(z) - jo_fixed2int(s) * jo_fixed2int(x)),
		toFIXED(jo_fixed2int(t) * (jo_fixed2int(z) * jo_fixed2int(z)) + jo_fixed2int(c))
    };
}


vec3 mat4multpoint(vec3 *vec, mat4 *mat) {
    vec3 result;
result.x = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_11) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_21) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_31) + 1 * jo_fixed2int(mat->_41));
result.y = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_12) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_22) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_32) + 1 * jo_fixed2int(mat->_42));
result.z = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_13) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_23) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_33) + 1 * jo_fixed2int(mat->_43));
return result;
}


vec3 mat4multvector(vec3 *vec, mat4 *mat) {
vec3 result;
result.x = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_11) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_21) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_31) + 0 * jo_fixed2int(mat->_41));
result.y = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_12) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_22) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_32) + 0 * jo_fixed2int(mat->_42));
result.z = toFIXED(jo_fixed2int(vec->x) * jo_fixed2int(mat->_13) + jo_fixed2int(vec->y) * jo_fixed2int(mat->_23) + jo_fixed2int(vec->z) * jo_fixed2int(mat->_33) + 0 * jo_fixed2int(mat->_43));
return result;
}

vec3 mat3multvector(vec3 *vec, mat3 *mat) {
vec3 result;
vec3 vect1 = {mat->_11, mat->_21, mat->_31};
vec3 vect2 = {mat->_12, mat->_22, mat->_32};
vec3 vect3 = {mat->_13, mat->_23, mat->_33};
result.x = dotvec3(vec, &vect1);
result.y = dotvec3(vec, &vect2);
result.z = dotvec3(vec, &vect3);
return result;
}
    
mat4 transformeuler(vec3 *scale, vec3 *eulerRotation, vec3 *translate) {
mat4 matscal = vecscale(scale);
mat4 matrot = mat4rotation(eulerRotation->x, eulerRotation->y, eulerRotation->z);
mat4 matpos = vectranslation(translate);
mat4 firstmul = mat4mul(&matscal, &matrot);
return mat4mul(&firstmul, &matpos);
}

mat4 transformaxis(vec3 *scale, vec3 *rotationAxis, FIXED rotationAngle, vec3 *translate) {
mat4 matscal = vecscale(scale);
mat4 matrot = mat4axisangle(rotationAxis,rotationAngle);
mat4 matpos = vectranslation(translate);
mat4 firstmul = mat4mul(&matscal,&matrot);
return mat4mul(&firstmul, &matpos);
}

mat4 lookatpoint(vec3 *position, vec3 *target, vec3 *up) {
vec3 subvec = vec3sub(target, position);
vec3 forward = normalizedvec3(&subvec);
vec3 crossrighter = cross(up, &forward);
vec3 right = normalizedvec3(&crossrighter);
vec3 newUp = cross(&forward, &right);
mat4 retmat;
mat4construct(&retmat,right.x, newUp.x, forward.x, 0, right.y, newUp.y, forward.y, 0, right.z, newUp.z, forward.z, 0, -dotvec3(&right, position), -dotvec3(&newUp, position), -dotvec3(&forward, position), JO_FIXED_1);
return retmat;
}

mat4 Projection(FIXED fov, FIXED aspect, FIXED zNear, FIXED zFar) {
FIXED tanHalfFov = slTan(DEGtoANG(slMulFX(fov, JO_FIXED_1_DIV_2)));
FIXED fovY = slDivFX(JO_FIXED_1, tanHalfFov); 
FIXED fovX = slDivFX(fovY, aspect); 
mat4 result;
result._11 = fovX;
result._22 = fovY;
result._33 = slDivFX(zFar, (zFar - zNear));
result._34 = JO_FIXED_1;
result._43 = slMulFX(-zNear, result._33);
result._44 = 0;
return result;
}

mat4 Ortho(FIXED left, FIXED right, FIXED bottom, FIXED top, FIXED zNear, FIXED zFar) { 
FIXED _11 = slDivFX(JO_FIXED_2, (right - left));
FIXED _22 = slDivFX(JO_FIXED_2, (top - bottom));
FIXED _33 = slDivFX(JO_FIXED_1, (zFar - zNear));
FIXED _41 = slDivFX((left + right), (left - right));
FIXED _42 = slDivFX((top + bottom), (bottom - top));
FIXED _43 = slDivFX((zNear), (zNear - zFar));
mat4 retmat;
mat4construct(&retmat, _11, 0, 0, 0, 0, _22, 0, 0, 0, 0, _33, 0, _41, _42, _43, JO_FIXED_1);
return retmat;
}

