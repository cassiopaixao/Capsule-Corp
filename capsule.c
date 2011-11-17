
#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"

void capsule_print_params(Capsule *c) {

	printf("Parâmetros da cápsula:\n");

	printf("h = %lf \n" \
	 "a = %lf \n" \
	 "d = %lf \n" \
	 "alpha = %lf \n" \
	 "delta = %lf \n" \
	 "t_0 = %lf \n" \
	 "t_inicial = %lf \n" \
	 "theta_crit = %lf \n" \
	 "theta_0 = %lf \n", 
	 c->a, c->d, c->alpha, c->delta,
	 c->t_0, c->t_inicial, c->theta_crit,
	 c->theta_0);

	printf("pos = ");
	v3d_print(&c->pos);
	printf("\n");

	printf("vel = ");
	v3d_print(&c->vel);
	printf("\n");
}
