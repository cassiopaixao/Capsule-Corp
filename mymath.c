#include "mymath.h"

#include <stdlib.h>
#include <stdio.h>

void v3d_set(v3d *v, double x, double y, double z) {
	v->x = x;
	v->y = y;
	v->z = z;
}


void v3d_print(v3d *v) {
	printf("(%lf, %lf, %lf)", v->x, v->y, v->z);
}
