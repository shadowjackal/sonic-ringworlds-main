#ifndef SSV_H
#define SSV_H

extern VECTOR ANORMS[162];

typedef struct {
    POINT *pntbl;
    Uint32 *fnorm;
    Uint32 *vnorm;
} JKLR;

typedef struct
{
    XPDATA       data;
    JKLR        *rdata;
    int framecount;
}               jklmesh;

extern Uint16 GRreal_ptr;

extern jklmesh                ssv_load(char *file_input, int textureoffset, bool gouraud, bool staticd);
extern int                ssv_update(jklmesh *mesh, int frame);

#define GRTBL(r,g,b) (((b&0x1f)<<10) | ((g&0x1f)<<5) | (r&0x1f))
extern Uint16 GourTbl[32];
#endif
