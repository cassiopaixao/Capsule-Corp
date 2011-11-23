
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

	struct SRing *next_ring, *prev_ring;
} Ring;

typedef struct SCover {
	Ring ring;

	unsigned int bursted;
	v3d normal;
	double t; // tempo
	double new_temp;
	double last_temp;
} Cover;

typedef struct SMesh {
	Ring *rings;
	Cover cover;
	unsigned int n_rings;
	struct SCapsule *cap;
} Mesh;

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

	Mesh mesh;
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
double tile_update_temp(Tile *t);

//-------------*-----------
// Anel
// ------------*-----------

void ring_init(Ring *ring, Capsule *cap, double l, double L);
void ring_calc_temp(Ring *ring);
void ring_neighborhood_temp(Ring *ring, double *t1, double *t2);
void ring_print(Ring *ring, FILE* file);
void ring_update_temp(Ring *ring);

// --------------*-------------
// Calota
// --------------*-------------

void cover_init(Cover *c, Capsule *capsule);
void cover_calc_temp(Cover *c);
double cover_update_temp(Cover *c);
void cover_print(Cover *c, FILE* file);

// --------------*-------------
// Mesh
// --------------*-------------

void mesh_init(Mesh *m, Capsule *cap);
void mesh_print(Mesh *m, FILE* file);

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
void capsule_output(Capsule *capsule, const char* filename);

#endif
