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
#include "newmath.h"
#include "newcol.h"

#define LWRAM 0x00200000
#define LWRAM_HEAP_SIZE 1048576
#define HWRAM 0x06040000
#define HWRAM_SIZE 262144


jo_camera               cam;
char                    *file_name;
int                     qicount;
int                     big_chungometer;
unsigned int            megainfo;

FIXED time;
jklmesh player;
jklmesh player0;
jklmesh orbinautmodel;
jklmesh lvlmodel;

int frame;
int ja = 0;
FIXED zpitch;
FIXED zroll;
FIXED zyaw;
FIXED livec[XYZ];
bool CollisionBool = true;
bool viewbool = true;


typedef struct orbinaut {
    POINT pos;
    ROTATE orientation;
    FIXED rot;
} orbinaut;

orbinaut orb1;

playerobject sonic;

typedef float NORML[3];

static Uint8 GourWork[4096];
static GOURAUDTBL Gour[4096];

#define    GRrealBase    0xe000

jo_palette                  image_pal;
jo_palette                  image_pal2;


PDATA trimesh;
Collision colmesh;

float fframe =0;
int b = 0;
bool jump;

FIXED camera_speed_x;
FIXED camera_speed_z;

void                camera_rotate(ANGLE angle_change) {
        //we can't just change camera angle because it will affect camera's chasing parameters
        //calculate old angle and dist
        FIXED dist_sq = (slMulFX((cam.target_pos.x - cam.viewpoint_pos.x),(cam.target_pos.x - cam.viewpoint_pos.x))+
                     slMulFX((cam.target_pos.z - cam.viewpoint_pos.z),(cam.target_pos.z - cam.viewpoint_pos.z))); 
        FIXED dist = slSquartFX(dist_sq);
        ANGLE camera_rotation = slAtan(-(cam.target_pos.x - cam.viewpoint_pos.x),-(cam.target_pos.z - cam.viewpoint_pos.z));
        camera_rotation += angle_change;
        //update camera coordinates acconding to new angle
        cam.viewpoint_pos.x = cam.target_pos.x + slMulFX(slCos(camera_rotation),dist);
        cam.viewpoint_pos.z = cam.target_pos.z + slMulFX(slSin(camera_rotation),dist);
}

void                my_gamepad(void)
{
    FIXED tmp_sin;
    FIXED tmp_cos;
   // bool bMovementActive = false;
    ANGLE camera_rotation = slAtan(-(cam.target_pos.x - cam.viewpoint_pos.x),-(cam.target_pos.z - cam.viewpoint_pos.z));
    switch (jo_get_input_direction_pressed(0))
    {
    case LEFT: 
        //bMovementActive = true;
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
        sonic.spd[X] += tmp_sin; 
        sonic.spd[Z] -= tmp_cos;
    break;
    case RIGHT:
        //bMovementActive = true;
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
        sonic.spd[X] -= tmp_sin;
        sonic.spd[Z] += tmp_cos;
    break;
    case UP: 
       // bMovementActive = true;
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
        sonic.spd[X] -= tmp_cos;
        sonic.spd[Z] -= tmp_sin;
    break;
    case DOWN: 
       // bMovementActive = true;
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
        sonic.spd[X] += tmp_cos;
        sonic.spd[Z] += tmp_sin; 
    break;
    case UP_LEFT:  
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
       // bMovementActive = true;
        sonic.spd[X] += tmp_sin;
        sonic.spd[X] -= tmp_cos;
        sonic.spd[Z] -= tmp_cos; 
        sonic.spd[Z] -= tmp_sin;
    break;
    case UP_RIGHT:  
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
        //bMovementActive = true;
        sonic.spd[X] -= tmp_sin;
        sonic.spd[X] -= tmp_cos;
        sonic.spd[Z] += tmp_cos;
        sonic.spd[Z] -= tmp_sin;
    break;
    case DOWN_LEFT: 
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
       // bMovementActive = true;
        sonic.spd[X] += tmp_sin; 
        sonic.spd[X] += tmp_cos; 
        sonic.spd[Z] -= tmp_cos;
        sonic.spd[Z] += tmp_sin;
    break;
    case DOWN_RIGHT:  
        tmp_sin = slMulFX(slSin(camera_rotation),sonic.acceleration);
        tmp_cos = slMulFX(slCos(camera_rotation),sonic.acceleration);
       // bMovementActive = true;
        sonic.spd[X] -= tmp_sin;
        sonic.spd[X] += tmp_cos; 
        sonic.spd[Z] += tmp_cos;
        sonic.spd[Z] += tmp_sin;
    break;
    case NONE:
        break;
    }
    if (jo_is_input_key_pressed(0, JO_KEY_L)) {
        //rotate camera around player, but keep player in place 
        camera_rotate(400);
    }
    if (jo_is_input_key_pressed(0, JO_KEY_R)) {
        //rotate camera around player, but keep player in place 
        camera_rotate(-400);
    }


    if (jo_is_input_key_pressed(0, JO_KEY_X)) {
        frame -= 1;
    };

    if (jo_is_input_key_down(0, JO_KEY_Y))
        CollisionBool = !CollisionBool;
    if (jo_is_input_key_down(0, JO_KEY_Z))
        viewbool = !viewbool;


    if (jo_is_input_key_down(0, JO_KEY_A) && sonic.gnd == true) {
        sonic.spd[Y] = -toFIXED(4);
        sonic.spd[X] += slMulFX(toFIXED(5), slSin(sonic.orientation[Z]));
        sonic.spd[Z] += slMulFX(toFIXED(5), slSin(sonic.orientation[X]));
        sonic.gnd = false;
        sonic.jmp = true;
    }


    if (jo_is_input_key_down(0, JO_KEY_B))
        sonic.pos[Y] -= toFIXED(25);

    if (jo_is_input_key_down(0, JO_KEY_C))
        frame += 1;

  

    //applying traction among all axises with a normal rate (no brakes)
    for (int i=0;i<3;i++)
    {
        if(sonic.spd[i] < 0) {   
            sonic.spd[i] += sonic.traction; 
            if (sonic.spd[i] > 0) sonic.spd[i] = 0;//speed zeroing check
        }
        else {
            sonic.spd[i] -= sonic.traction;
            if (sonic.spd[i] < 0) sonic.spd[i] = 0;//speed zeroing check
        }
    }

    //speed saturation check
    for (int i=0;i<3;i++)
    {
        if(JO_ABS(sonic.spd[i]) > sonic.max_speed)
        {
            if (sonic.spd[i] > 0)
                sonic.spd[i] = sonic.max_speed; 
            else
                sonic.spd[i] = -sonic.max_speed; 
        } 
    }

    //calculating player's rotation based on his speed along x and z axises (he does'n rotate along y yet)
    //only updating rotation when actually moving
    if ( (sonic.spd[X]!=0) && (sonic.spd[Z]!=0) )
        sonic.orientation[Y] = slAtan(-sonic.spd[X],sonic.spd[Z]);

    if(sonic.gnd == false) {
        sonic.orientation[X] = LerpAng(sonic.orientation[X],toFIXED(0),toFIXED(0.10));
        sonic.orientation[Z] = LerpAng(sonic.orientation[Z],toFIXED(0),toFIXED(0.10));
    };

    //coordinate update
    for (int i=0;i<3;i+=2){
        sonic.pos[i] += sonic.spd[i];
    }


    sonic.pos[Y] += slMulFX(sonic.spd[X], slSin(sonic.orientation[Z])) + slMulFX(sonic.spd[Z], slSin(sonic.orientation[X])) ;
    sonic.pos[Y]+= sonic.spd[Y];

}

Plane gplane;
    Sphere undersonc;   

rungles ringlist[20];


void othercollision(void) {

}
    FIXED bnb = 0;

    float ringframe = 0;

void			    my_draw(void)
{
    my_gamepad();

    
 
    cam.viewpoint_pos.y = toFIXED(-35);
    //FIXED campos[XYZ];
   // FIXED target[XYZ];
    //target[0] = sonic.pos[0];
    //target[1] = sonic.pos[1];
    //target[2] = sonic.pos[2];
    sonic.col.pos[0] = sonic.pos[0];
    sonic.col.pos[1] = sonic.pos[1];
    sonic.col.pos[2] = sonic.pos[2];

    cam.target_pos.x = sonic.pos[0];
    cam.target_pos.y = sonic.pos[1]+toFIXED(-12);
    cam.target_pos.z = sonic.pos[2];

    //using chasing camera that always focuses on player
    if ( (camera_speed_x != 0) && (camera_speed_z!=0) )
    {
        cam.viewpoint_pos.x += camera_speed_x;
        cam.viewpoint_pos.z += camera_speed_z;
        //update camera rotation change due to chasing
        //camera_rotation = slAtan(-(cam.target_pos.x - cam.viewpoint_pos.x),-(cam.target_pos.z - cam.viewpoint_pos.z));
    }
    //updating speed
    FIXED dist_sq = (slMulFX((cam.target_pos.x - cam.viewpoint_pos.x),(cam.target_pos.x - cam.viewpoint_pos.x))+
                     slMulFX((cam.target_pos.z - cam.viewpoint_pos.z),(cam.target_pos.z - cam.viewpoint_pos.z))); 
    if (dist_sq > toFIXED(2500))
    {
        //move camera closer1
        camera_speed_x = slMulFX((cam.target_pos.x - cam.viewpoint_pos.x),toFIXED(0.08));
        camera_speed_z = slMulFX((cam.target_pos.z - cam.viewpoint_pos.z),toFIXED(0.08));
    }
    else
    {
        camera_speed_x = 0;
        camera_speed_z = 0;
    }

if(CollisionBool == true) {
    if(Collision_SphereColResolve(&sonic.col,&colmesh)) {
   // sonic.gnd = true;
    } //else sonic.gnd = false;
} 
    //Collision_SpherePlaneResolve(&sonic.col,&gplane);

    bnb -= toFIXED(0.65);
    if(bnb < toFIXED(-5)) bnb = 0;

    undersonc.pos[0] = sonic.col.pos[0];
    undersonc.pos[1] = sonic.col.pos[1] + toFIXED(11);
    undersonc.pos[2] = sonic.col.pos[2];
    undersonc.radius = toFIXED(2);

    //slPrintFX(undersonc.pos[0],slLocate(20,7));
    //slPrintFX(undersonc.pos[1],slLocate(20,8));
    //slPrintFX(undersonc.pos[2],slLocate(20,9));
        jo_3d_camera_look_at(&cam);

    if(CollisionBool == true) {
        if(Collision_SphereCol_bool_special(&undersonc, &colmesh) || Collision_SpherePlane_bool(&undersonc,&gplane)) {sonic.gnd = true; 
        sonic.jmp = false;}else sonic.gnd = true;
    }
    if(sonic.gnd == false)  {sonic.spd[1] += toFIXED(0.15);} else {
        sonic.spd[1] = 0;
        //sonic.jmp == false;
    }
    sonic.pos[0] = sonic.col.pos[0];
    sonic.pos[1] = sonic.col.pos[1];
    sonic.pos[2] = sonic.col.pos[2];



    //jo_core_exec_on_slave(othercollision);

    int animstart = 0;
    int animend = 0;

    if((JO_ABS(sonic.spd[X]) >= toFIXED(0) || JO_ABS(sonic.spd[Z]) >= toFIXED(0)) && ((JO_ABS(sonic.spd[X]) < toFIXED(0.049) || JO_ABS(sonic.spd[Z]) < toFIXED(0.049)))) {
            animstart = 0; 
            animend = 8;
    };

    if((JO_ABS(sonic.spd[X]) >= toFIXED(0.05) || JO_ABS(sonic.spd[Z]) >= toFIXED(0.05)) && ((JO_ABS(sonic.spd[X]) < toFIXED(3) || JO_ABS(sonic.spd[Z]) < toFIXED(3)))) {
            animstart = 21; 
            animend = 37;
    };

    if((JO_ABS(sonic.spd[X]) > toFIXED(3) || JO_ABS(sonic.spd[Z]) > toFIXED(3)) && ((JO_ABS(sonic.spd[X]) < toFIXED(5) || JO_ABS(sonic.spd[Z]) < toFIXED(5)))) {
            animstart = 38; 
            animend = 58;
    };

    if((JO_ABS(sonic.spd[X]) > toFIXED(5) || JO_ABS(sonic.spd[Z]) > toFIXED(5)) && ((JO_ABS(sonic.spd[X]) < toFIXED(10) || JO_ABS(sonic.spd[Z]) < toFIXED(10)))) {
            animstart = 59; 
            animend = 63;
    };
    //if()
    //if()
    jo_vdp2_clear_bitmap_nbg1(JO_COLOR_Black);


    fframe += 0.25;
    if(fframe >= animend) fframe = animstart;
    if(fframe < animstart) fframe = animstart;


    ssv_update(&player,fframe);

    slPrintFX((sonic.pos[0]),slLocate(0,0));
    slPrintFX((sonic.pos[1]),slLocate(0,1));
    slPrintFX((sonic.pos[2]),slLocate(0,2));
    slPrintFX(slAng2FX(sonic.orientation[0]), slLocate(20,0));
    slPrintFX(slAng2FX(sonic.orientation[Y]), slLocate(20,1));
    slPrintFX(slAng2FX(sonic.orientation[2]), slLocate(20,2));

    b += (1);

    FIXED pos[][XYZS] =
    {    
        {0,0,0,toFIXED(0.15f)},   
    };
    
    SPR_ATTR attr[] =
    {
        SPR_ATTRIBUTE(18, 0, No_Gouraud, CL32KRGB | ECdis | SPenb, sprNoflip | FUNC_Sprite | _ZmCC),
    };

    ringframe += 0.15f;
        
    SPR_ATTR ringattr[] =
    {
        SPR_ATTRIBUTE((int)ringframe+21, 0, No_Gouraud, CL32KRGB | ECdis | SPenb, sprNoflip | FUNC_Sprite | _ZmCC),
    };

    //ANGLE camera_rotation = slAtan(-(cam.target_pos.x - cam.viewpoint_pos.x),-(cam.target_pos.z - cam.viewpoint_pos.z));
    

    slPushMatrix();
    {
        slTranslate(sonic.pos[0], sonic.pos[1], sonic.pos[2]);
        if(sonic.jmp != true)slScale(toFIXED(1),toFIXED(1),toFIXED(1)); else slScale(toFIXED(2),toFIXED(2),toFIXED(2));
        if(sonic.jmp != true)slRotZ(sonic.orientation[Z]);
        if(sonic.jmp != true)slRotX(sonic.orientation[X]); else slRotX(DEGtoANG(bnb));
        if(sonic.jmp != true)slRotY(sonic.orientation[Y]); else slRotY(sonic.orientation[Y]+DEGtoANG(90));

        if(sonic.jmp == false)slPutPolygonX(&player.data, livec); else slPutPolygonX(&player0.data, livec);
    }
    slPopMatrix();

    if(ringframe > 4) ringframe = 0;

    slPushMatrix();
    {
        slTranslate(orb1.pos[0],orb1.pos[1],orb1.pos[2]);
    //    slPutPolygonX(&orbinautmodel.data,livec);
    //    slPutSprite((FIXED *)pos[0], (SPR_ATTR *)(&(attr[0].texno)), 0);
    }
    slPopMatrix();

    for(int i = 0; i < 20; i++) {
        if(ringlist[i].collectedornot == false)ringlist[i].collectedornot = Collision_SphereSphere(&sonic.col, &ringlist[i].col);
    }


    for(int i = 0; i < 20; i++) {
    if(ringlist[i].collectedornot == false) {
    slPushMatrix();
    {
    FIXED POSSY[][XYZS] =
    {    
        {ringlist[i].col.pos[0], ringlist[i].col.pos[1], ringlist[i].col.pos[2],toFIXED(0.5)},   
    };        slPutSprite((FIXED *)POSSY[0], (SPR_ATTR *)(&(ringattr[0].texno)), 0);
    }
    slPopMatrix();
    }}


 
        if(viewbool)slPutPolygonX(&lvlmodel.data,livec);
        //slPutPolygonX(&trimesh,(VECTOR){0,0,0});
        //slPutPolygonX(&orbinautmodel,livec);


    jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(0), 0);
       // slScale(toFIXED(0.75),toFIXED(0.75),toFIXED(0.75));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
        // /jo_background_3d_plane_a_draw(true);
    }
    jo_3d_pop_matrix();

         jo_3d_push_matrix();
    {
        slTranslate(0, toFIXED(-200), 0);
        slScale(toFIXED(1),toFIXED(1),toFIXED(1));
        slRotX(DEGtoANG(90));
        slRotY(DEGtoANG(0));
        slRotZ(DEGtoANG(0));
       // jo_background_3d_plane_b_draw(true);
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

    jo_3d_camera_init(&cam);

    jo_set_tga_palette_handling(my_tga_palette_handling);    //jo_set_tga_palette_handling(&image_pal);

    jo_img_8bits    img;

   // img.data = JO_NULL;
   // jo_tga_8bits_loader(&img, JO_ROOT_DIR, "BGI.TGA", 0);
   // jo_vdp2_set_nbg1_8bits_image(&img, image_pal.id, true);
   // jo_free_img(&img);
jo_vdp2_clear_bitmap_nbg1(JO_COLOR_Black);

    //colmesh.points = HWRAM;
    sonic.acceleration = toFIXED(0.1f);
    sonic.traction = toFIXED(0.05f);
    sonic.max_speed = toFIXED(8);

    cam.viewpoint_pos.x = toFIXED(-20);
    cam.viewpoint_pos.y = toFIXED(-20);
    cam.viewpoint_pos.z = toFIXED(-20);
    sonic.pos[0] = toFIXED(5);
    sonic.pos[1] = toFIXED(-11);
    sonic.pos[2] = toFIXED(0);
    cam.target_pos.x = sonic.pos[0];
    cam.target_pos.y = sonic.pos[1]+toFIXED(-12);
    cam.target_pos.z = sonic.pos[2];
    slDynamicFrame(1);
    //ssvlist0 = jo_3d_create_mesh(ssv_quad_count("CBE.SSV"));
    jo_printf(0,0,"Loading : 0 %%");
    jo_fixed_point_time();
    slPrint("bfr",slLocate(0,6));

    player = ssv_load("OBJ.SSV",0,true,false);
    slPrint("SNC",slLocate(0,6));
    player0 = ssv_load("BAL.SSV",14,true,false);
        slPrint("BAL",slLocate(0,6));

   // orbinautmodel = ssv_load("ORB.SSV",19,false,false);
        slPrint("ORB",slLocate(0,6));

    //lvlmodel = ssv_load("TNL.SSV",20,true,false);
    //    slPrint("LVL",slLocate(0,6));
//
    //mesh2coltri(&lvlmodel,&colmesh);
//
    //jo_free(lvlmodel.rdata->fnorm);
    //jo_free(lvlmodel.rdata->vnorm);
    //jo_free(lvlmodel.rdata->pntbl);
    //jo_free(lvlmodel.data.pntbl);
    //jo_free(lvlmodel.data.vntbl);
    //jo_free(lvlmodel.data.pltbl);
  //  jo_free(lvlmodel.data.attbl);
    lvlmodel = ssv_load("TNL.SSV",25,false,false);


    orb1.pos[0] = toFIXED(-10);
    orb1.pos[1] = toFIXED(-10);
    orb1.pos[2] = toFIXED(-10);
    livec[0] = toFIXED(-0.5f);
    livec[1] = toFIXED(0.5f);
    livec[2] = toFIXED(0);
    gplane.normal[0] = 0;
    gplane.normal[1] = -toFIXED(1);
    gplane.normal[2] = 0;
    gplane.distance = 0;

    sonic.pos[1] = toFIXED(-11);
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
    jo_sprite_add_tga("14", "BL0.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("15", "BL1.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("16", "BL2.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("17", "BL3.TGA", JO_COLOR_Transparent);
    jo_sprite_add_tga("18","CHP.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("19","EYE.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("20","TIL.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("21","RG0.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("22","RG1.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("23","RG2.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("24","RG3.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("25","18.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("26","20.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("27","26.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("28","13.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("29","25.TGA",JO_COLOR_Transparent);
    //jo_sprite_add_tga("30","16.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("25","TN1.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("26","TN2.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("27","TN3.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("28","TN4.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("29","TN5.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("30","TN6.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("31","TN7.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("32","TN8.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("33","TN9.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("34","TN10.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("35","TN11.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("36","TN12.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("37","TN13.TGA",JO_COLOR_Transparent);
    jo_sprite_add_tga("38","TN14.TGA",JO_COLOR_Transparent);



    for(int i = 0; i < (int)lvlmodel.data.nbPolygon; i++) {
    lvlmodel.data.attbl[i].sort |= SORT_MAX; 
    }
//
    sonic.col.radius = toFIXED(11);

    jo_enable_background_3d_plane(JO_COLOR_Black);


    	//jo_printf(0, 0, jo_get_last_error());

    // FLOOR
    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "GRS.TGA", 0);
    jo_background_3d_plane_a_img(&img, image_pal.id, true, true);
    jo_free_img(&img);

    img.data = JO_NULL;
    jo_tga_8bits_loader(&img, JO_ROOT_DIR, "SKY.TGA", 1);
    jo_background_3d_plane_b_img(&img, image_pal.id, true, true);
    //jo_core_add_callback(my_gamepad);
    jo_free_img(&img);

    for(int i = 0; i < 20; i++) {
    ringlist[i].collectedornot = false;
    ringlist[i].col.radius = toFIXED(10);
    }

    ringlist[0].col.pos[0] = toFIXED(-30);
    ringlist[0].col.pos[1] = toFIXED(-10);
    ringlist[0].col.pos[2] = toFIXED(0);

    ringlist[1].col.pos[0] = toFIXED(-30);
    ringlist[1].col.pos[1] = toFIXED(-10);
    ringlist[1].col.pos[2] = toFIXED(-50);

    ringlist[2].col.pos[0] = toFIXED(-30);
    ringlist[2].col.pos[1] = toFIXED(-10);
    ringlist[2].col.pos[2] = toFIXED(-100);

    ringlist[3].col.pos[0] = toFIXED(-30);
    ringlist[3].col.pos[1] = toFIXED(-10);
    ringlist[3].col.pos[2] = toFIXED(-150);

    ringlist[4].col.pos[0] = toFIXED(-30);
    ringlist[4].col.pos[1] = toFIXED(-10);
    ringlist[4].col.pos[2] = toFIXED(-200);

    ringlist[5].col.pos[0] = toFIXED(-45);
    ringlist[5].col.pos[1] = toFIXED(-15);
    ringlist[5].col.pos[2] = toFIXED(-400);

    ringlist[6].col.pos[0] = toFIXED(-45);
    ringlist[6].col.pos[1] = toFIXED(-15);
    ringlist[6].col.pos[2] = toFIXED(-450);

    ringlist[7].col.pos[0] = toFIXED(-45);
    ringlist[7].col.pos[1] = toFIXED(-15);
    ringlist[7].col.pos[2] = toFIXED(-500);

    ringlist[8].col.pos[0] = toFIXED(-45);
    ringlist[8].col.pos[1] = toFIXED(-15);
    ringlist[8].col.pos[2] = toFIXED(-550);

    ringlist[9].col.pos[0] = toFIXED(-45);
    ringlist[9].col.pos[1] = toFIXED(-15);
    ringlist[9].col.pos[2] = toFIXED(-600);

    ringlist[10].col.pos[0] = toFIXED(-50);
    ringlist[10].col.pos[1] = toFIXED(-10);
    ringlist[10].col.pos[2] = toFIXED(200);

    ringlist[11].col.pos[0] = toFIXED(-50);
    ringlist[11].col.pos[1] = toFIXED(-10);
    ringlist[11].col.pos[2] = toFIXED(150);

    ringlist[12].col.pos[0] = toFIXED(-50);
    ringlist[12].col.pos[1] = toFIXED(-10);
    ringlist[12].col.pos[2] = toFIXED(200);

    ringlist[13].col.pos[0] = toFIXED(-50);
    ringlist[13].col.pos[1] = toFIXED(-15);
    ringlist[13].col.pos[2] = toFIXED(250);

    ringlist[14].col.pos[0] = toFIXED(-50);
    ringlist[14].col.pos[1] = toFIXED(-15);
    ringlist[14].col.pos[2] = toFIXED(300);

    ringlist[15].col.pos[0] = toFIXED(-50);
    ringlist[15].col.pos[1] = toFIXED(-15);
    ringlist[15].col.pos[2] = toFIXED(350);

    ringlist[16].col.pos[0] = toFIXED(-50);
    ringlist[16].col.pos[1] = toFIXED(-15);
    ringlist[16].col.pos[2] = toFIXED(400);

    ringlist[17].col.pos[0] = toFIXED(-50);
    ringlist[17].col.pos[1] = toFIXED(-10);
    ringlist[17].col.pos[2] = toFIXED(450);

    ringlist[18].col.pos[0] = toFIXED(-50);
    ringlist[18].col.pos[1] = toFIXED(-10);
    ringlist[18].col.pos[2] = toFIXED(200);

    ringlist[19].col.pos[0] = toFIXED(-50);
    ringlist[19].col.pos[1] = toFIXED(-10);
    ringlist[19].col.pos[2] = toFIXED(200);



//jo_vdp2_zoom_nbg1(0.3f);

    jo_core_add_vblank_callback(slGouraudTblCopy);
	jo_core_add_callback(my_draw);
	jo_core_run();
}

//5 is width/height of the orb so 10 for entire radius.


