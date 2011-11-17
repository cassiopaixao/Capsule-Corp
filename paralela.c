#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"
#include "mymath.h"

#define DEBUG

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

	// configura parametros
	parse_input();

	#ifdef DEBUG
	// imprime parametros
	capsule_print_params(&capsule);
	#endif /* DEBUG */

	// inicializa a capsula (t = 0)
	capsule_init(&capsule);

	//executa as iterações
	capsule_iterate(&capsule);

	// imprime saida
	capsule_output(&capsule);
	
	return EXIT_SUCCESS;
}
