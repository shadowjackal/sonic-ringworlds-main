#include <jo/jo.h>
#include "meshprim.h"
#include "VecMath.h"

void        meshcreatecube(jo_3d_mesh *obium, int xwidth, int ywidth, int height) {

        jo_3d_set_mesh_vertice(obium, -xwidth, height, ywidth, 0);
        jo_3d_set_mesh_vertice(obium, xwidth, height, ywidth, 1);
        jo_3d_set_mesh_vertice(obium, xwidth, height, -ywidth, 2);
        jo_3d_set_mesh_vertice(obium, -xwidth, height, -ywidth, 3);   

        jo_3d_set_mesh_vertice(obium, -xwidth, -height, ywidth, 4);
        jo_3d_set_mesh_vertice(obium,  xwidth, -height, ywidth, 5);
        jo_3d_set_mesh_vertice(obium,  xwidth, -height, -ywidth, 6);
        jo_3d_set_mesh_vertice(obium, -xwidth, -height, -ywidth, 7);

        jo_3d_set_mesh_vertice(obium, -xwidth, -height,   ywidth, 8);
        jo_3d_set_mesh_vertice(obium, -xwidth,  height,   ywidth, 9);
        jo_3d_set_mesh_vertice(obium, -xwidth,  height, -ywidth, 10);
        jo_3d_set_mesh_vertice(obium, -xwidth, -height, -ywidth, 11);

        jo_3d_set_mesh_vertice(obium, xwidth, -height, ywidth, 12);
        jo_3d_set_mesh_vertice(obium, xwidth, height, ywidth, 13);
        jo_3d_set_mesh_vertice(obium, xwidth, height, -ywidth, 14);
        jo_3d_set_mesh_vertice(obium, xwidth, -height, -ywidth, 15);

        jo_3d_set_mesh_vertice(obium, -xwidth, height, ywidth, 16);
        jo_3d_set_mesh_vertice(obium, xwidth, height, ywidth, 17);
        jo_3d_set_mesh_vertice(obium, xwidth, -height, ywidth, 18);
        jo_3d_set_mesh_vertice(obium, -xwidth, -height, ywidth, 19);

        jo_3d_set_mesh_vertice(obium, -xwidth, height, -ywidth, 20);
        jo_3d_set_mesh_vertice(obium, xwidth, height, -ywidth, 21);
        jo_3d_set_mesh_vertice(obium, xwidth, -height, -ywidth, 22);
        jo_3d_set_mesh_vertice(obium, -xwidth, -height, -ywidth, 23);

        jo_3d_set_mesh_color(obium, JO_COLOR_Green);
}