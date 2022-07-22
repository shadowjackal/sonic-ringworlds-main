/*
** Jo Sega Saturn Engine
** Copyright (c) 2012-2021, Johannes Fetz (johannesfetz@gmail.com)
** All rights reserved.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are met:
**     * Redistributions of source code must retain the above copyright
**       notice, this list of conditions and the following disclaimer.
**     * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in the
**       documentation and/or other materials provided with the distribution.
**     * Neither the name of the Johannes Fetz nor the
**       names of its contributors may be used to endorse or promote products
**       derived from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
** ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
** WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
** DISCLAIMED. IN NO EVENT SHALL Johannes Fetz BE LIABLE FOR ANY
** DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
** (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
** LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
** ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <jo/jo.h>
#include <string.h>
#include "ssv.h"
#include "jklmath.h"
#include "jklcol.h"

#define LWRAM 0x00200000
#define LWRAM_HEAP_SIZE 1310720


jo_camera               cam;
char                    *file_name;
int                     qicount;
int                     big_chungometer;
unsigned int            megainfo;

FIXED time;
jklmesh player;
jklmesh player0;

int frame;
int ja = 0;
FIXED zpitch;
FIXED zroll;
FIXED zyaw;
FIXED livec[XYZ];


typedef struct playerobject {
    FIXED x,y,z;
    int state;
    ANGLE yrot;
    vec3ang orientation;
} playerobject;

playerobject sonic;

typedef float NORML[3];

static Uint8 GourWork[4096];
static GOURAUDTBL Gour[4096];

#define    GRrealBase    0xe000

jo_palette                  image_pal;
jo_palette                  image_pal2;


tri testtri;

PDATA trimesh;


float fframe =0;
int b = 0;
ANGLE rotation;
bool jump;

void                my_gamepad(void)
{
    switch (jo_get_input_direction_pressed(0))
    {
    case LEFT: 
    sonic.x -= toFIXED(5); 
    break;
    case RIGHT: 
    sonic.x += toFIXED(5); 
    break;
    case UP: 
    sonic.z += toFIXED(5); 
    break;
    case DOWN: 
    sonic.z -= toFIXED(5); 
    break;
    case UP_LEFT:  
    sonic.x -= toFIXED(5); 
    sonic.z += toFIXED(5);  
    break;
    case UP_RIGHT:  
    sonic.x += toFIXED(5); 
    sonic.z += toFIXED(5); 
    break;
    case DOWN_LEFT: 
    sonic.x -= toFIXED(5); 
    sonic.z -= toFIXED(5); 
    break;
    case DOWN_RIGHT:  
    sonic.x += toFIXED(5); 
    sonic.z -= toFIXED(5); 
    break;
    case NONE: break;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_L)) {
        rotation -= 400;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_R)) {
        rotation += 400;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_X))
        zpitch += toFIXED(3);
    if (jo_is_input_key_pressed(0, JO_KEY_Y))
        zroll += toFIXED(3);
    if (jo_is_input_key_pressed(0, JO_KEY_Z))
        zyaw += toFIXED(3);

    if (jo_is_input_key_down(0, JO_KEY_A))
        sonic.y -= toFIXED(1);
    if (jo_is_input_key_down(0, JO_KEY_B))
        sonic.y += toFIXED(1);

    if(jump) {
        sonic.y += 1;
    }
}


void			    my_draw(void)
{
    my_gamepad();
 
    Plane plna = PlaneFromTri(testtri);
    cam.viewpoint_pos.y = toFIXED(-35);
    vec3 campos;
    vec3 target;
    target.x = sonic.x;
    target.y = sonic.y;
    target.z = sonic.z;
    vec3orbit(&campos,&target,toFIXED(60),rotation);

    cam.viewpoint_pos.x = campos.x;
    cam.viewpoint_pos.z = campos.z;

    cam.target_pos.x = sonic.x;
    cam.target_pos.y = sonic.y;
    cam.target_pos.z = sonic.z;

    fframe += 0.25;
    if(fframe >= 80) fframe = 0;
   // jo_printf(12, 1, "*ANIM TEST*");
   ssv_update(&player,(int)fframe);

    Plane floor;
    floor.normal.asArray[0] = toFIXED(0);
    floor.normal.asArray[1] = toFIXED(-1);
    floor.normal.asArray[2] = toFIXED(0);
    floor.distance = 0;
    Sphere camsphere;
    camsphere.position.asArray[0] = sonic.x;
    camsphere.position.asArray[1] = sonic.y+toFIXED(-6);
    camsphere.position.asArray[2] = sonic.z;
    camsphere.radius = toFIXED(6);

    sonic.orientation = dirtoeuler(&plna.normal);

    slPrintFX((sonic.x),slLocate(0,1));
    slPrintFX((sonic.y),slLocate(0,2));
    slPrintFX((sonic.z),slLocate(0,3));
    if(TriangleSphere(&testtri, &camsphere)) jo_printf(0,4,"BOOL : True"); else jo_printf(0,4,"BOOL : False");
    b += (1);
    jo_3d_camera_look_at(&cam);
    slPushMatrix();
    {
        slTranslate(sonic.x, sonic.y+toFIXED(-12), sonic.z);
        slScale(toFIXED(1),toFIXED(1),toFIXED(1));
        slRotX(sonic.orientation.asArray[0]);
        slRotY(sonic.orientation.asArray[1]);
        slRotZ(sonic.orientation.asArray[2]);
        slPutPolygonX(&player.data, livec);
    }
    slPopMatrix();

    slPushMatrix();
    {
        slPutPolygonX(&trimesh,(VECTOR){0,0,0});
    }
    slPopMatrix();

    jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(0), 0);
        slScale(toFIXED(0.75),toFIXED(0.75),toFIXED(0.75));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
        jo_background_3d_plane_a_draw(true);
    }
    jo_3d_pop_matrix();

         jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(-200), 0);
        slScale(toFIXED(1),toFIXED(1),toFIXED(1));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
        jo_background_3d_plane_b_draw(true);
    }
    jo_3d_pop_matrix();
}

jo_palette          *my_tga_palette_handling(void)
{
    // We create a new palette for each image. It's not optimal but OK for a demo.
    jo_create_palette(&image_pal);
    return (&image_pal);
}

void			jo_main(void)
{      

    slInitGouraud(Gour, 4096, GRrealBase, GourWork);
    slIntFunction(slGouraudTblCopy); // Copy Gour to VRAM each frame
    slSetGouraudTbl(GourTbl);

    jo_add_memory_zone((unsigned char *)LWRAM, LWRAM_HEAP_SIZE);

	jo_core_init(JO_COLOR_Blue);
    jo_memory_fragmentation();
    
    jo_3d_camera_init(&cam);

    cam.viewpoint_pos.x = 0;
    cam.viewpoint_pos.y = toFIXED(-20);
    cam.viewpoint_pos.z = 0;
    cam.target_pos.x = sonic.x;
    cam.target_pos.y = sonic.y+toFIXED(-12);
    cam.target_pos.z = sonic.z;
    //slDynamicFrame(1);
    //ssvlist0 = jo_3d_create_mesh(ssv_quad_count("CBE.SSV"));
    jo_printf(0,0,"Loading : 0 %%");
    jo_fixed_point_time();
    player = ssv_load("OBJ.SSV",0,true,false);


    //player0 = ssv_load("BAL.SSV",14,true,false);
    jo_printf(0,0,"Loading : 25 %%");
    jo_printf(0,0,"Loading : 50 %%");
    jo_printf(0,0,"Loading : 75 %%");
    jo_printf(0,0,"Loading : 100 %%");

    //est_load(ssvlist0);

   // jo_printf(1, 4, "Thing: %d", ssv_quad_count(file_name));

    livec[0] = toFIXED(-0.5f);
    livec[1] = toFIXED(0.5f);
    livec[2] = toFIXED(0);
    jo_sprite_add_tga("00", "S01.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("01", "S02.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("02", "S03.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("03", "S04.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("04", "S05.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("05", "S06.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("06", "S07.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("07", "S08.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("08", "S09.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("09", "S10.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("10", "S11.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("11", "S12.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("12", "S14.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("13", "S15.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("14", "BL0.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("15", "BL1.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("16", "BL2.TGA", JO_COLOR_Transparent);
    //jo_sprite_add_tga("17", "BL3.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("18","CHP.TGA",JO_COLOR_Transparent);
    testtri = tridef((vec3){toFIXED(0),toFIXED(-50),toFIXED(0)},(vec3){toFIXED(50),toFIXED(0),toFIXED(0)},(vec3){toFIXED(0),toFIXED(0),toFIXED(50)});
    trimesh = coltri2mesh(&testtri);

    jo_set_tga_palette_handling(my_tga_palette_handling);    //jo_set_tga_palette_handling(&image_pal);
    jo_img_8bits    img;

    jo_enable_background_3d_plane(JO_COLOR_Black);

    // FLOOR
    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "GRS.TGA", 0);
    jo_background_3d_plane_a_img(&img, image_pal.id, true, true);
    jo_free_img(&img);

    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "SKY.TGA", 0);
    jo_background_3d_plane_b_img(&img, image_pal.id, true, true);
    jo_free_img(&img);

    //jo_core_add_callback(my_gamepad);
    jo_core_add_vblank_callback(slGouraudTblCopy);
	jo_core_add_callback(my_draw);
	jo_core_run();
}
