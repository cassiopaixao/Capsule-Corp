
#ifndef CAPSULE_H
#define CAPSULE_H

#include "mymath.h"

typedef struct STile {
	double t_0,             // tempo inicial
	       t,               // tempo corrente
	       alpha,           // equação de atrito
	       delta,           // equação de atrito
	       last_temp,       // temperatura inicial,
	       max_temp;        // temperatura em que a pastilha explode
	       v3d vel,         // velocidade
	           pos;         // posição

	double d,       // base da pastilha
	       dl;      // altura da pastilha

	int    bursted; // 1 se a pastilha estourou

	v3d    normal;  // vetor normal

	struct SRing *ring;     // anel em que a pastilha está

	v3d    edges[4];

	struct STile *left, *right;     // pastilhas ao lado

	double new_temp;        // nova temperatura
} Tile;

typedef struct SRing {
	Tile *tiles; // array de tiles
	unsigned int n_tiles;
	struct SCapsule *capsule;

        double  r0,             // raio menor (inferior)
                r1,             // raio maior (superior)
                temp;        // temperatura média das pastilhas

	struct SRing *next_ring, *previous_ring;
} Ring;

typedef struct SCapsule {

	double h,               // altura
	       a,               // fator de forma do parabolóide
	       d,               // lado da pastilha
	       alpha,           // parâmetro da função de atrito
	       delta,           // parâmetro da função de atrito
	       t_0,             // tempo inicial
	       t_inicial,       
	       theta_crit,      // temperatura em que as pastilhas explodem
	       theta_0;         // temperatura inicial das pastilhas

	v3d    pos,     // posição inicial
	       vel;     // velocidade

	unsigned int steps;

        Ring* anel;     // aneis
        
        int qtd_aneis;  // quantidade de anéis

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
// Atualiza temperatura da pastilha
void tile_update_temp(Tile *t);

//-------------*-----------
// Anel
// ------------*-----------

void ring_neighborhood_temp(Ring *ring, double *t1, double *t2);

// --------------*-------------
// Cápsula
// --------------*-------------
// Imprime parâmetros da capsula
void capsule_print_params(Capsule *capsule);
// Inicializa a capsula (situação t=0)
void capsule_init(Capsule *capsule);
// Realiza as iterações
void capsule_iterate(Capsule *capsule);
// Imprime a saída
void capsule_output(Capsule *capsule);

#endif
