#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"
#include "mymath.h"

static Capsule capsule;

double prox_linha(FILE *arq){

	double val;

	//char Linha[100];
	//fgets(Linha, 100, arq);  // o 'fgets' lê até 99 caracteres ou até o '\n'
	//return atof (Linha);
	
	fscanf(arq, "%lf", &val);
	return val;
}

static inline void parse_input(const char* filename) {
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
	
	FILE *arq;

	arq = fopen(filename,"rt");   /* Arquivo ASCII, para leitura */


	capsule.h = prox_linha(arq);
	capsule.a = prox_linha(arq);
	capsule.d = prox_linha(arq);
	capsule.alpha = prox_linha(arq);
	capsule.t_0 = prox_linha(arq);
	capsule.delta = prox_linha(arq);
	capsule.theta_crit = prox_linha(arq);
	capsule.theta_0 = prox_linha(arq);

	double x, y, z;
	fscanf(arq, "%lf %lf %lf\n", &x, &y, &z);
	v3d_set(&capsule.pos, x, y, z);

	fscanf(arq, "%lf %lf %lf\n", &x, &y, &z);
	v3d_set(&capsule.vel, x, y, z);

	//capsule.steps = prox_linha(arq);
	fscanf(arq, "%lu", &capsule.steps);

	fclose(arq);
}

int main(int argc, char **argv) {

    // configura parametros
    parse_input("entrada.txt");

    #ifdef DEBUG
    // imprime parametros
    capsule_print_params(&capsule);
    #endif /* DEBUG */

    // inicializa a capsula (t = 0)
    capsule_init(&capsule);

    //executa as iterações
    capsule_iterate(&capsule);

    // imprime saida
    capsule_output(&capsule, "saida.txt");
   
    return EXIT_SUCCESS;
}
