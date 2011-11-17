#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"
#include "mymath.h"

static Capsule capsule;

static inline void parse_input(void) {

	capsule.h = 122.;
	capsule.a = 1.;
	capsule.d = 2.;
	capsule.alpha = 2.;
	capsule.delta = 1.;
	capsule.t_0 = 0.;
	capsule.t_inicial = 0.;
	capsule.theta_crit = 1050.;
	capsule.theta_0 = 1000.;

	v3d_set(&capsule.pos, 1., 1., 1.);
	v3d_set(&capsule.vel, 1., 1., 1.);

	capsule.steps = 20;
}

int main(int argc, char **argv) {

	parse_input();

	capsule_print_params(&capsule);

	return EXIT_SUCCESS;
}
