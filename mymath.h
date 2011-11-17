#ifndef MYMATH_H
#define MYMATH_H

typedef struct {
	double x, y, z;
} v3d;

void v3d_set(v3d *v, double x, double y, double z);
void v3d_print(v3d *v);

#endif
