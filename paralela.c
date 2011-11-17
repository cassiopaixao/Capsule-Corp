#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"
#include "mymath.h"

static Capsule capsule;

static inline void parse_input() {
/*
    i["h"] = 122.0
    i["a"] = 1.0
    i["d"] = 2.0
    i["alpha"] = 2.0
    i["delta"] = 1.0
    i["t_0"] = 0.0
    i["t_inicial"] = 0
    i["theta_crit"] = 1050 #?
    i["theta_0"] = 1000.0
    i["pos"] = (1.0, 1.0, 1.0)
    i["vel"] = (1.0, 1.0, 1.0)
    i["steps"] = 20
*/
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
