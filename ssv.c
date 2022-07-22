#include <jo/jo.h>
#include "ssv.h"

Uint16 GRreal_ptr = 0xe000;


VECTOR ANORMS[]=
{
    POStoFIXED(0, -1, 0),
    POStoFIXED(0.723608, -0.447221, 0.525724),
    POStoFIXED(-0.276387, -0.447221, 0.850649),
    POStoFIXED(-0.894427, -0.447214, 0),
    POStoFIXED(-0.276387, -0.447221, -0.850649),
    POStoFIXED(0.723608, -0.447221, -0.525724),
    POStoFIXED(0.276387, 0.447221, 0.850649),
    POStoFIXED(-0.723608, 0.447221, 0.525724),
    POStoFIXED(-0.723608, 0.447221, -0.525724),
    POStoFIXED(0.276387, 0.447221, -0.850649),
    POStoFIXED(0.894427, 0.447214, 0),
    POStoFIXED(0, 1, 0),
    POStoFIXED(-0.23282, -0.657519, 0.716564),
    POStoFIXED(-0.162454, -0.850655, 0.499995),
    POStoFIXED(-0.0776063, -0.96795, 0.23885),
    POStoFIXED(0.203179, -0.96795, 0.147618),
    POStoFIXED(0.425322, -0.850655, 0.309011),
    POStoFIXED(0.609548, -0.657518, 0.442857),
    POStoFIXED(0.531942, -0.502302, 0.681711),
    POStoFIXED(0.262868, -0.525739, 0.809011),
    POStoFIXED(-0.0296365, -0.502302, 0.864184),
    POStoFIXED(0.81273, -0.502301, -0.295235),
    POStoFIXED(0.850649, -0.525735, 0),
    POStoFIXED(0.81273, -0.502301, 0.295235),
    POStoFIXED(0.203179, -0.96795, -0.147618),
    POStoFIXED(0.425322, -0.850655, -0.309011),
    POStoFIXED(0.609548, -0.657518, -0.442857),
    POStoFIXED(-0.753442, -0.657514, 0),
    POStoFIXED(-0.525728, -0.850653, 0),
    POStoFIXED(-0.251145, -0.96795, 0),
    POStoFIXED(-0.483971, -0.502301, 0.716565),
    POStoFIXED(-0.688189, -0.525736, 0.499998),
    POStoFIXED(-0.831052, -0.502297, 0.23885),
    POStoFIXED(-0.23282, -0.657519, -0.716564),
    POStoFIXED(-0.162454, -0.850655, -0.499995),
    POStoFIXED(-0.0776063, -0.96795, -0.23885),
    POStoFIXED(-0.831052, -0.502297, -0.23885),
    POStoFIXED(-0.688189, -0.525736, -0.499998),
    POStoFIXED(-0.483971, -0.502301, -0.716565),
    POStoFIXED(-0.0296365, -0.502302, -0.864184),
    POStoFIXED(0.262868, -0.525739, -0.809011),
    POStoFIXED(0.531942, -0.502302, -0.681711),
    POStoFIXED(0.956626, 0.251149, 0.147618),
    POStoFIXED(0.951059, 0, 0.30901),
    POStoFIXED(0.860698, -0.251149, 0.442858),
    POStoFIXED(0.860698, -0.251149, -0.442858),
    POStoFIXED(0.951059, 0, -0.30901),
    POStoFIXED(0.956626, 0.251149, -0.147618),
    POStoFIXED(0.155213, 0.251149, 0.955423),
    POStoFIXED(0, 0, 1),
    POStoFIXED(-0.155213, -0.251149, 0.955423),
    POStoFIXED(0.687159, -0.251149, 0.681716),
    POStoFIXED(0.587786, 0, 0.809016),
    POStoFIXED(0.436006, 0.251149, 0.864189),
    POStoFIXED(-0.860698, 0.251149, 0.442858),
    POStoFIXED(-0.951059, 0, 0.30901),
    POStoFIXED(-0.956626, -0.251149, 0.147618),
    POStoFIXED(-0.436006, -0.251149, 0.864189),
    POStoFIXED(-0.587786, 0, 0.809016),
    POStoFIXED(-0.687159, 0.251149, 0.681716),
    POStoFIXED(-0.687159, 0.251149, -0.681716),
    POStoFIXED(-0.587786, 0, -0.809016),
    POStoFIXED(-0.436006, -0.251149, -0.864189),
    POStoFIXED(-0.956626, -0.251149, -0.147618),
    POStoFIXED(-0.951059, 0, -0.30901),
    POStoFIXED(-0.860698, 0.251149, -0.442858),
    POStoFIXED(0.436006, 0.251149, -0.864189),
    POStoFIXED(0.587786, 0, -0.809016),
    POStoFIXED(0.687159, -0.251149, -0.681716),
    POStoFIXED(-0.155213, -0.251149, -0.955423),
    POStoFIXED(0, 0, -1),
    POStoFIXED(0.155213, 0.251149, -0.955423),
    POStoFIXED(0.831052, 0.502297, 0.23885),
    POStoFIXED(0.688189, 0.525736, 0.499998),
    POStoFIXED(0.483971, 0.502301, 0.716565),
    POStoFIXED(0.0296365, 0.502302, 0.864184),
    POStoFIXED(-0.262868, 0.525739, 0.809011),
    POStoFIXED(-0.531942, 0.502302, 0.681711),
    POStoFIXED(-0.81273, 0.502301, 0.295235),
    POStoFIXED(-0.850649, 0.525735, 0),
    POStoFIXED(-0.81273, 0.502301, -0.295235),
    POStoFIXED(-0.531942, 0.502302, -0.681711),
    POStoFIXED(-0.262868, 0.525739, -0.809011),
    POStoFIXED(0.0296365, 0.502302, -0.864184),
    POStoFIXED(0.483971, 0.502301, -0.716565),
    POStoFIXED(0.688189, 0.525736, -0.499998),
    POStoFIXED(0.831052, 0.502297, -0.23885),
    POStoFIXED(0.0776063, 0.96795, 0.23885),
    POStoFIXED(0.162454, 0.850655, 0.499995),
    POStoFIXED(0.23282, 0.657519, 0.716564),
    POStoFIXED(0.753442, 0.657514, 0),
    POStoFIXED(0.525728, 0.850653, 0),
    POStoFIXED(0.251145, 0.96795, 0),
    POStoFIXED(-0.203179, 0.96795, 0.147618),
    POStoFIXED(-0.425322, 0.850655, 0.309011),
    POStoFIXED(-0.609548, 0.657518, 0.442857),
    POStoFIXED(-0.203179, 0.96795, -0.147618),
    POStoFIXED(-0.425322, 0.850655, -0.309011),
    POStoFIXED(-0.609548, 0.657518, -0.442857),
    POStoFIXED(0.0776063, 0.96795, -0.23885),
    POStoFIXED(0.162454, 0.850655, -0.499995),
    POStoFIXED(0.23282, 0.657519, -0.716564),
    POStoFIXED(0.361798, 0.894431, -0.26286),
    POStoFIXED(0.638192, 0.723611, -0.262864),
    POStoFIXED(0.44721, 0.723612, -0.525728),
    POStoFIXED(-0.138195, 0.894431, -0.425317),
    POStoFIXED(-0.052788, 0.723612, -0.688185),
    POStoFIXED(-0.361803, 0.723613, -0.587778),
    POStoFIXED(-0.447209, 0.89443, 0),
    POStoFIXED(-0.670816, 0.723612, -0.162457),
    POStoFIXED(-0.670816, 0.723612, 0.162457),
    POStoFIXED(-0.138195, 0.894431, 0.425317),
    POStoFIXED(-0.361803, 0.723613, 0.587778),
    POStoFIXED(-0.052788, 0.723612, 0.688185),
    POStoFIXED(0.361798, 0.894431, 0.26286),
    POStoFIXED(0.44721, 0.723612, 0.525728),
    POStoFIXED(0.638192, 0.723611, 0.262864),
    POStoFIXED(0.861805, 0.276395, -0.425321),
    POStoFIXED(0.80902, 0, -0.587782),
    POStoFIXED(0.670821, 0.276395, -0.68819),
    POStoFIXED(-0.138199, 0.276395, -0.951056),
    POStoFIXED(-0.309015, 0, -0.951057),
    POStoFIXED(-0.447214, 0.276395, -0.85065),
    POStoFIXED(-0.947214, 0.276394, -0.162457),
    POStoFIXED(-1, 0, 0),
    POStoFIXED(-0.947214, 0.276394, 0.162457),
    POStoFIXED(-0.447214, 0.276395, 0.85065),
    POStoFIXED(-0.309015, 0, 0.951057),
    POStoFIXED(-0.138199, 0.276395, 0.951056),
    POStoFIXED(0.670821, 0.276395, 0.68819),
    POStoFIXED(0.80902, 0, 0.587782),
    POStoFIXED(0.861805, 0.276395, 0.425321),
    POStoFIXED(0.309015, 0, -0.951057),
    POStoFIXED(0.447213, -0.276398, -0.850649),
    POStoFIXED(0.138199, -0.276398, -0.951055),
    POStoFIXED(-0.80902, 0, -0.587782),
    POStoFIXED(-0.670819, -0.276394, -0.688192),
    POStoFIXED(-0.861803, -0.276395, -0.425325),
    POStoFIXED(-0.80902, 0, 0.587782),
    POStoFIXED(-0.861803, -0.276395, 0.425325),
    POStoFIXED(-0.670819, -0.276394, 0.688192),
    POStoFIXED(0.309015, 0, 0.951057),
    POStoFIXED(0.138199, -0.276398, 0.951055),
    POStoFIXED(0.447213, -0.276398, 0.850649),
    POStoFIXED(1, 0, 0),
    POStoFIXED(0.947214, -0.276394, 0.162457),
    POStoFIXED(0.947214, -0.276394, -0.162457),
    POStoFIXED(0.361803, -0.723613, -0.587778),
    POStoFIXED(0.138195, -0.89443, -0.42532),
    POStoFIXED(0.052788, -0.723612, -0.688185),
    POStoFIXED(-0.44721, -0.723612, -0.525728),
    POStoFIXED(-0.361798, -0.894431, -0.26286),
    POStoFIXED(-0.638195, -0.72361, -0.262861),
    POStoFIXED(-0.638194, -0.723609, 0.262864),
    POStoFIXED(-0.361799, -0.894429, 0.262865),
    POStoFIXED(-0.44721, -0.723612, 0.525728),
    POStoFIXED(0.670816, -0.723612, -0.162457),
    POStoFIXED(0.670816, -0.723612, 0.162457),
    POStoFIXED(0.447211, -0.894429, 0),
    POStoFIXED(0.052788, -0.723612, 0.688185),
    POStoFIXED(0.138199, -0.894429, 0.42532),
    POStoFIXED(0.361806, -0.723612, 0.587778),
};

Uint16 GourTbl[32] = {
    GRTBL( 0,  0,  0),
    GRTBL( 0,  0,  0),
    GRTBL( 1,  1,  1),
    GRTBL( 1,  1,  1),
    GRTBL( 2,  2,  2),
    GRTBL( 2,  2,  2),
    GRTBL( 3,  3,  3),
    GRTBL( 3,  3,  3),
    GRTBL( 4,  4,  4),
    GRTBL( 4,  4,  4),
    GRTBL( 5,  5,  5),
    GRTBL( 5,  5,  5),
    GRTBL( 6,  6,  6),
    GRTBL( 6,  6,  6),
    GRTBL( 7,  7,  7),
    GRTBL( 7,  7,  7),
    GRTBL(11, 11, 11),
    GRTBL(13, 13, 13),
    GRTBL(14, 14, 14),
    GRTBL(16, 16, 16),
    GRTBL(17, 17, 17),
    GRTBL(18, 18, 18),
    GRTBL(19, 19, 19),
    GRTBL(20, 20, 20),
    GRTBL(21, 21, 21),
    GRTBL(22, 22, 22),
    GRTBL(23, 23, 23),
    GRTBL(25, 25, 25),
    GRTBL(26, 26, 26),
    GRTBL(26, 26, 26),
    GRTBL(27, 27, 27),
    GRTBL(27, 27, 27),
    GRTBL(29, 29, 29),
};

jklmesh                ssv_load(char *file_input, int textureoffset, bool gouraud, bool staticd)
{
    jklmesh              mesh;
   // mesh = (jklmesh *)jo_malloc_with_behaviour(sizeof(*mesh), JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    //loads in file

    char *stream = jo_fs_read_file(file_input,NULL);

    int offset = 0;
    int qi;
    int texids[2000];
    int colorrgb[2000][3];

    memcpy(&mesh.data.nbPoint, &stream[offset], sizeof(unsigned int));
    jo_printf(0,2,"vertex %d",mesh.data.nbPoint);
    offset += sizeof(unsigned int);

    memcpy(&mesh.data.nbPolygon, &stream[offset], sizeof(unsigned int));
    jo_printf(0,3,"polygon %d",mesh.data.nbPolygon);
    offset += sizeof(unsigned int);

    memcpy(&mesh.framecount, &stream[offset], sizeof(unsigned int));
    jo_printf(0,4,"frame %d",mesh.framecount);
    offset += sizeof(unsigned int);

    mesh.rdata = (JKLR *)jo_malloc_with_behaviour(mesh.framecount * sizeof(*mesh.rdata),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    mesh.data.pntbl = (POINT *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.pntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    mesh.data.vntbl = (VECTOR *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.vntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);

    for(int i = 0; i < mesh.framecount; i++) {
        mesh.rdata[i].pntbl = (POINT *)jo_malloc_with_behaviour(mesh.data.nbPoint * sizeof(*mesh.data.pntbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }
    
    
    mesh.data.pltbl = (POLYGON *)jo_malloc_with_behaviour(mesh.data.nbPolygon * sizeof(*mesh.data.pltbl),JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    for(int i = 0; i < mesh.data.nbPolygon; i++) {
    for(int ii = 0; ii < 4; ii++){
        memcpy(&mesh.data.pltbl[i].Vertices[ii], &stream[offset], sizeof(short));
        offset += sizeof(short);
    }
        memcpy(&texids[i],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][0],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][1],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

        memcpy(&colorrgb[i][2],&stream[offset],sizeof(unsigned int));
        offset += sizeof(unsigned int);

    }


    ////memcpy(mesh.data.pntbl, &stream[offset], sizeof(POINT) * mesh.data.nbPoint);
    for(int i = 0; i < mesh.framecount; i++) {
    memcpy(mesh.rdata[i].pntbl, &stream[offset], sizeof(POINT) * mesh.data.nbPoint);
    offset += sizeof(POINT)*mesh.data.nbPoint;
    }

    mesh.data.pntbl = mesh.rdata[0].pntbl;


    for(int i = 0; i < mesh.framecount; i++) {
    mesh.rdata[i].fnorm = jo_malloc_with_behaviour(sizeof(Uint32)*mesh.data.nbPolygon,JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }

    for(int i = 0; i < mesh.framecount; i++) {
    for(int j = 0; j < mesh.data.nbPolygon; j++) {
        memcpy(&mesh.rdata[i].fnorm[j],&stream[offset],sizeof(Uint32));
        offset += sizeof(Uint32);
    }
    }
//
    for(int i = 0; i < mesh.framecount; i++) {
    mesh.rdata[i].vnorm = jo_malloc_with_behaviour(sizeof(Uint32)*mesh.data.nbPoint,JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    }
//
    for(int i = 0; i < mesh.framecount; i++) {
    for(int j = 0; j < mesh.data.nbPoint; j++) {
        memcpy(&mesh.rdata[i].vnorm[j],&stream[offset],sizeof(Uint32));
        offset += sizeof(Uint32);
    }
    }
//
//
    for(int i = 0; i < mesh.data.nbPolygon; i++) {
        mesh.data.pltbl[i].norm[0] = -ANORMS[mesh.rdata[0].fnorm[i]][0];
        mesh.data.pltbl[i].norm[1] = -ANORMS[mesh.rdata[0].fnorm[i]][1];
        mesh.data.pltbl[i].norm[2] = -ANORMS[mesh.rdata[0].fnorm[i]][2];

    }

    for(int i = 0; i < mesh.data.nbPoint; i++) {
        mesh.data.vntbl[i][0] = -ANORMS[mesh.rdata[0].vnorm[i]][0];
        mesh.data.vntbl[i][1] = -ANORMS[mesh.rdata[0].vnorm[i]][1];
        mesh.data.vntbl[i][2] = -ANORMS[mesh.rdata[0].vnorm[i]][2];
    }

    //sets attributes of faces (transparency, mesh effect, etc)
    mesh.data.attbl = (ATTR *)jo_malloc_with_behaviour(mesh.data.nbPolygon * sizeof(*mesh.data.attbl), JO_MALLOC_TRY_REUSE_SAME_BLOCK_SIZE);
    for(int i = 0; i < mesh.data.nbPolygon; i++)
    {
        //jo_color coloroid = JO_COLOR_SATURN_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]);
        //if(i != 0) {
        ATTR proxy_attribute = ATTRIBUTE(Single_Plane,SORT_CEN, texids[i]+textureoffset, No_Palet, 0, CL32KRGB|MESHoff, sprNoflip, No_Option);
        ATTR pro_attribute = ATTRIBUTE(Single_Plane,SORT_CEN,No_Texture,C_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]),0,MESHoff,sprPolygon,No_Option);
        
        if(gouraud == true)proxy_attribute.atrb |= CL_Gouraud; 
        if(gouraud == true)proxy_attribute.sort |= UseGouraud; else proxy_attribute.sort |= No_Gouraud;
        if(gouraud == true)proxy_attribute.gstb = GRreal_ptr++; else proxy_attribute.gstb = NULL;
        if(gouraud == true)pro_attribute.sort |= UseGouraud; else pro_attribute.sort |= No_Gouraud;
        if(gouraud == true)pro_attribute.gstb = GRreal_ptr++; else pro_attribute.gstb = NULL;
        if(gouraud == true)pro_attribute.atrb |= CL_Gouraud;
        if(texids[i] != -1) {
        mesh.data.attbl[i] = proxy_attribute;
        } else {
        mesh.data.attbl[i] = pro_attribute;
        }
        //} else {
        //ATTR proxy_attribute = ATTRIBUTE(Dual_Plane,SORT_MAX, texids[i], No_Palet, 0, CL32KRGB|MESHoff, sprNoflip, No_Option);
        //ATTR pro_attribute = ATTRIBUTE(Dual_Plane,SORT_MAX,No_Texture,C_RGB(colorrgb[i][0],colorrgb[i][1],colorrgb[i][2]),0,MESHoff,sprPolygon,No_Option);
        

        if(texids[i] != -1) {
        mesh.data.attbl[i] = proxy_attribute;
        } else {
        mesh.data.attbl[i] = pro_attribute;
        }
        
    }    
    if(staticd == true) jo_free(mesh.rdata->pntbl);
    return mesh;
}

int                ssv_update(jklmesh *mesh, int frame) {
    if(mesh->framecount < frame) {
        return -1;
    }

    mesh->data.pntbl = mesh->rdata[frame].pntbl;
    
    for(int i = 0; i < mesh->data.nbPolygon; i++) {
        mesh->data.pltbl[i].norm[0] = -ANORMS[mesh->rdata[frame].fnorm[i]][0];
        mesh->data.pltbl[i].norm[1] = -ANORMS[mesh->rdata[frame].fnorm[i]][1];
        mesh->data.pltbl[i].norm[2] = -ANORMS[mesh->rdata[frame].fnorm[i]][2];

    }

    for(int i = 0; i < mesh->data.nbPoint; i++) {
        mesh->data.vntbl[i][0] = -ANORMS[mesh->rdata[frame].vnorm[i]][0];
        mesh->data.vntbl[i][1] = -ANORMS[mesh->rdata[frame].vnorm[i]][1];
        mesh->data.vntbl[i][2] = -ANORMS[mesh->rdata[frame].vnorm[i]][2];
    }
}
