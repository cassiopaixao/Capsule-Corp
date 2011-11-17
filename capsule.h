
#ifndef CAPSULE_H
#define CAPSULE_H

#include "mymath.h"

typedef struct STile {
	double t_0, // tempo inicial
	       t,   // tempo corrente
	       alpha, // equação de atrito
	       delta, // equação de atrito
	       last_temp, // temperatura inicial,
	       max_temp; // temperatura em que a pastilha explode
	       v3d vel, // velocidade
	           pos; // posição
	double d, // base da pastilha
	       dl; // altura da pastilha
	int    bursted; // 1 se a pastilha estourou
	v3d    normal;
	struct SRing *ring;
	v3d    edges[4];
	struct STile *left, *right;
	double new_temp;
} Tile;

typedef struct SRing {

} Ring;

typedef struct {

	double h,         // altura
	       a,         // fator de forma do parabolóide
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

//-------------*-----------
// Pastilha
// ------------*-----------
// Inicializa pastilha
void tile_init(Tile *t, Capsule *cap, v3d *a, v3d *b, v3d *c, v3d *d, Ring *ring);
// Define pastilhas vizinhas
void tile_link(Tile *t, Tile *left, Tile *right);
// Calcula temperatura do próximo timestep
void tile_calc_temp(Tile *t);

//-------------*-----------
// Anel
// ------------*-----------

void ring_neighborhood_temp(Ring *ring, double *t1, double *t2);

// --------------*-------------
// Cápsula
// --------------*-------------
void capsule_print_params(Capsule *capsule);
void capsule_init(Capsule *capsule);
void capsule_iterate(Capsule *capsule);
void capsule_output(Capsule *capsule);

#endif
