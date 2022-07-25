#ifndef __NEW_MATH_H__
#define __NEW_MATH_H__

#include <jo/jo.h>

typedef ANGLE ROTATE[XYZ];
#define ROTtoANG(x, y, z)   { DEGtoANG(x), DEGtoANG(y), DEGtoANG(z) }

////////// FIXED-POINT FUNCTIONS //////////

#define FLOAT2FIXED(x)  ((FIXED)((x) * (float)toFIXED(1.0)))

static inline FIXED Lerp(const FIXED A, const FIXED B, const FIXED F) {
    return A + slMulFX(B - A, F);
}

static inline ANGLE LerpAng(const ANGLE A, const ANGLE B, const FIXED F) {
    FIXED A_to_FX = 360 * (Uint16)A;
    FIXED B_to_FX = 360 * (Uint16)B;
    
    FIXED diff = (B_to_FX - A_to_FX);
    
    if(diff > toFIXED(180.0))
        A_to_FX += toFIXED(360.0);
    else if(diff < toFIXED(-180.0))
        B_to_FX += toFIXED(360.0);
    
    FIXED res = Lerp(A_to_FX, B_to_FX, F);
    
    return (ANGLE)(res / 360);
}

static inline ANGLE Arccos(const FIXED x) {
	if(x > toFIXED(1.0))
        return DEGtoANG(0.0);
    else if(x < toFIXED(-1.0))
        return DEGtoANG(180.0);
    
	return slAtan(slDivFX(slSquartFX(toFIXED(1.0) - slMulFX(x, x)), x), toFIXED(1.0));
}

////////// VECTOR FUNCTIONS //////////

static const VECTOR vec_zero    = POStoFIXED(0, 0, 0);

static const VECTOR vec_right   = POStoFIXED(1, 0, 0);
static const VECTOR vec_up      = POStoFIXED(0, 1, 0);
static const VECTOR vec_forward = POStoFIXED(0, 0, 1);

static const VECTOR vec_left = POStoFIXED(-1, 0, 0);
static const VECTOR vec_down = POStoFIXED(0, -1, 0);
static const VECTOR vec_back = POStoFIXED(0, 0, -1);

static inline FIXED vec_length(VECTOR vec) {
    return slSquartFX(slInnerProduct(vec, vec));
}

static inline ANGLE vec_angle(VECTOR A, VECTOR B) {
    return Arccos(slDivFX(slMulFX(vec_length(A), vec_length(B)), slInnerProduct(A, B)));
}

static inline void vec_set(const FIXED x, const FIXED y, const FIXED z, VECTOR out) {
    out[X] = x;
    out[Y] = y;
    out[Z] = z;
}

static inline void vec_copy(const VECTOR vec, VECTOR out) {
    out[X] = vec[X];
    out[Y] = vec[Y];
    out[Z] = vec[Z];
}

static inline void vec_mul_scalar(const VECTOR vec, const FIXED scalar, VECTOR out) {
    out[X] = slMulFX(vec[X], scalar);
    out[Y] = slMulFX(vec[Y], scalar);
    out[Z] = slMulFX(vec[Z], scalar);
}

static inline void vec_add(const VECTOR vec1, const VECTOR vec2, VECTOR out) {
    out[X] = vec1[X] + vec2[X];
    out[Y] = vec1[Y] + vec2[Y];
    out[Z] = vec1[Z] + vec2[Z];
}

static inline void vec_sub(const VECTOR vec1, const VECTOR vec2, VECTOR out) {
    out[X] = vec1[X] - vec2[X];
    out[Y] = vec1[Y] - vec2[Y];
    out[Z] = vec1[Z] - vec2[Z];
}

static inline void vec_mul(const VECTOR vec1, const VECTOR vec2, VECTOR out) {
    out[X] = slMulFX(vec1[X], vec2[X]);
    out[Y] = slMulFX(vec1[Y], vec2[Y]);
    out[Z] = slMulFX(vec1[Z], vec2[Z]);
}

static inline void vec_lerp(const VECTOR vec1, const VECTOR vec2, const FIXED F, VECTOR out) {
    out[X] = Lerp(vec1[X], vec2[X], F);
    out[Y] = Lerp(vec1[Y], vec2[Y], F);
    out[Z] = Lerp(vec1[Z], vec2[Z], F);
}

static inline void vec_cross(const VECTOR vec1, const VECTOR vec2, VECTOR out) {
    out[X] = slMulFX(vec1[Y], vec2[Z]) - slMulFX(vec1[Z], vec2[Y]);
    out[Y] = slMulFX(vec1[Z], vec2[X]) - slMulFX(vec1[X], vec2[Z]);
    out[Z] = slMulFX(vec1[X], vec2[Y]) - slMulFX(vec1[Y], vec2[X]);
}

static inline void vec_reflect(VECTOR inDirection, VECTOR inNormal, VECTOR out) {
    FIXED factor = slMulFX(slInnerProduct(inNormal, inDirection), toFIXED(-2.0));
    vec_mul_scalar(inNormal, factor, out);
    vec_add(out, inDirection, out);
}

static inline void vec_sqrt(const VECTOR vec, VECTOR out) {
    out[X] = slSquartFX(ABS(vec[X]));
    out[Y] = slSquartFX(ABS(vec[Y]));
    out[Z] = slSquartFX(ABS(vec[Z]));
}

static inline void vec_normalize(VECTOR vec) {
    FIXED k = slDivFX(vec_length(vec), toFIXED(1.0));
    vec_mul_scalar(vec, k, vec);
}

static inline void vec_rotate_about(VECTOR p, VECTOR v, ANGLE a, VECTOR out) {
    FIXED ca = slCos(a), sa = slSin(a);
    FIXED t = toFIXED(1.0) - ca;
    FIXED x = v[X], y = v[Y], z = v[Z];
    
    VECTOR A = { ca + slMulFX(slMulFX(x, x), t),             slMulFX(slMulFX(x, y), t) - slMulFX(z, sa), slMulFX(slMulFX(x, z), t) + slMulFX(y, sa) };
    VECTOR B = { slMulFX(slMulFX(x, y), t) + slMulFX(z, sa), ca + slMulFX(slMulFX(y, y), t),             slMulFX(slMulFX(y, z), t) - slMulFX(x, sa) };
    VECTOR C = { slMulFX(slMulFX(z, x), t) - slMulFX(y, sa), slMulFX(slMulFX(z, y), t) + slMulFX(x, sa), ca + slMulFX(slMulFX(z, z), t) };
    
    out[X] = slInnerProduct(A, p);
    out[Y] = slInnerProduct(B, p);
    out[Z] = slInnerProduct(C, p);
}

////////// MATRIX FUNCTIONS //////////

static inline void set_view_mtx(const VECTOR pos, const ROTATE ang) {
    slRotX(-ang[X]);
    slRotY(-ang[Y]);
    slRotZ(-ang[Z]);
    slTranslate(-pos[X], -pos[Y], -pos[Z]);
}

static inline void set_model_mtx(const VECTOR pos, const ROTATE ang) {
    slTranslate(pos[X], pos[Y], pos[Z]);
    slRotZ(ang[Z]);
    slRotY(ang[Y]);
    slRotX(ang[X]);
}

#endif