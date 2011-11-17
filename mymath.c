#include "mymath.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

void v3d_set(v3d *v, double x, double y, double z) {

	assert(v != NULL);

	v->x = x;
	v->y = y;
	v->z = z;
}

void v3d_minus(v3d *a, v3d *b, v3d *res) {

	assert(a != NULL);
	assert(b != NULL);
	assert(res != NULL);

	res->x = a->x - b->x;
	res->y = a->y - b->y;
	res->z = a->z - b->z;
}


double v3d_length(v3d *v) {

	assert(v != NULL);

	return sqrt(
	 (v->x * v->x) +
	 (v->y * v->y) +
	 (v->z * v->z));
}

void v3d_copy(v3d *res, v3d *a) {
	
	assert(res != NULL);
	assert(a != NULL);

	v3d_set(res, a->x, a->y, a->z);
}

void v3d_cross(v3d *a, v3d *b, v3d *res) {

	assert(a != NULL);
	assert(b != NULL);
	assert(res != NULL);

	res->x = a->y * b->z - a->z * b->y;
	res->y = a->z * b->x - a->x * b->z;
	res->z = a->x * b->y - a->y * b->x;
}

void v3d_print(v3d *v) {

	assert(v != NULL);

	printf("(%lf, %lf, %lf)", v->x, v->y, v->z);
}
