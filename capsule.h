
#ifndef CAPSULE_H
#define CAPSULE_H

#include "mymath.h"

typedef struct {

	double h,         // altura
	       a,
	       d,         // lado da pastilha
	       alpha,     // parâmetro da função de atrito
	       delta,     // parâmetro da função de atrito
	       t_0,       // tempo inicial
	       t_inicial,
	       theta_crit,
	       theta_0;

	v3d    pos, // posição inicial
	       vel;

	unsigned int steps;

} Capsule;

void capsule_print_params(Capsule *capsule);

#endif
