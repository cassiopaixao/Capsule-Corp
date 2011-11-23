#ifndef MYMATH_H
#define MYMATH_H

typedef struct {
	double x, y, z;
} v3d;

#define MM_PI 3.1415926535897932

// v = (x, y, z)
void v3d_set(v3d *v, double x, double y, double z);
// res = a - b
void v3d_minus(v3d *a, v3d *b, v3d *res);
// |v|
double v3d_length(v3d *v);
// res = a
void v3d_copy(v3d *res, v3d *a);
void v3d_print(v3d *v);
// res = a x b
void v3d_cross(v3d *a, v3d *b, v3d *res);
// res = a . b
double v3d_dot(v3d *a, v3d *b);
#endif
