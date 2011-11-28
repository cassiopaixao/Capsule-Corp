
#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"

#include <math.h>
#include <assert.h>

void tile_init(Tile *t, Capsule *cap, v3d *a, v3d *b, v3d *c, v3d *d, Ring *ring) {
	v3d vdl, p, q;

	#ifdef DEBUG
	assert(t != NULL);
	assert(cap != NULL);
	assert(a != NULL);
	assert(b != NULL);
	assert(c != NULL);
	assert(d != NULL);
	#endif

	t->t_0 = cap->t_0;
	t->t = t->t_0;
	t->alpha = cap->alpha;
	t->delta = cap->delta;
	t->temp[0] = t->temp[1] = cap->theta_0;
	t->max_temp = cap->theta_crit;

	v3d_set(&t->vel, cap->vel.x, cap->vel.y, cap->vel.z);
	v3d_set(&t->pos, cap->pos.x, cap->pos.y, cap->pos.z);

	v3d_minus(c, a, &vdl);
	t->dl = v3d_length(&vdl);
	t->d = cap->d;
	t->bursted = 0;

	/*
	v3d_copy(&t->edges[0], a);
	v3d_copy(&t->edges[1], b);
	v3d_copy(&t->edges[2], c);
	v3d_copy(&t->edges[3], d);
	*/

	t->ring = ring;

	v3d_minus(b, a, &p);
	v3d_minus(d, a, &q);

	v3d_cross(&p, &q, &t->normal);

	v3d_normalize(&t->normal);
}


void tile_link(Tile *t, Tile *left, Tile *right) {

	#ifdef DEBUG
	assert(t != NULL);
	assert(left != NULL);
	assert(right != NULL);
	#endif

	t->left = left;
	t->right = right;
}

// média das temperaturas na vizinhança
static double _tile_perimenter_temp(Tile *t, const int not_step_mod_2) {

	#ifdef DEBUG
	assert(t != NULL);
	assert(t->left != NULL);
	assert(t->right != NULL);
	#endif

	return ((t->left->temp[not_step_mod_2] + t->right->temp[not_step_mod_2]) * t->dl +
		    (ring_neighborhood_temp_plus(t->ring)) * t->d
		    ) / (2.*(t->dl + t->d));
}

void tile_calc_temp(Tile *t, const int step_mod_2) {

	#ifdef DEBUG
	assert(t != NULL);
	assert(t->left != NULL);
	assert(t->right != NULL);
	#endif

	t->t += 1;
	if (t->bursted) { // estourou?
		// só rejunte
		t->temp[step_mod_2] = _tile_perimenter_temp(t, (step_mod_2+1)%2);
	} else {
		double vn,
			delta_atrito = 0.,
			delta_dissip;

		vn = v3d_dot(&t->normal, &t->vel);
		if (vn > 0.) {
			double val;
		
			val = t->t - t->t_0;
			delta_atrito = t->alpha * vn * atan(val*val);
		}
	
		delta_dissip = t->delta * abs(vn);
		t->temp[step_mod_2] = _tile_perimenter_temp(t, (step_mod_2+1)%2) + delta_atrito - delta_dissip;

		if (t->temp[step_mod_2] > t->ring->capsule->theta_crit) {
			t->bursted = 1; // estourou
			t->temp[step_mod_2] = (t->left->temp[(step_mod_2+1)%2] + t->right->temp[(step_mod_2+1)%2])/2.;
		}
	}
}

// BEGIN [RING]

inline unsigned int _ring_n_tiles(Ring *ring, double l) {
	#ifdef DEBUG
	assert(ring != NULL);
	assert(ring->capsule != NULL);
	#endif

	return ((unsigned int)
	 floor((2. * MM_PI * sqrt(l / ring->capsule->a)) / ring->capsule->d));
}

void ring_init(Ring *ring, Capsule *cap, double l, double L) {
	unsigned int n_tiles;
	double z0, z1, r0, r1;
	double alpha0, alpha1, divd;
	double beta0, beta1;
	static double two_pi = 2. * MM_PI;
	unsigned int i;
	double x0, y0, X0, Y0;
	double x1, y1, X1, Y1;
	v3d a, b, c, d;
	unsigned int next, prev;

	#ifdef DEBUG
	assert(ring != NULL);
	assert(cap != NULL);
	#endif

	// número de pastilhas

	ring->capsule = cap;
	n_tiles = _ring_n_tiles(ring, l);

	ring->n_tiles = n_tiles;
	ring->tiles = (Tile*) malloc(sizeof(Tile) * n_tiles);

	z0 = l;
	z1 = L + l;

	// raios
	r0 = sqrt(z0 / cap->a);
	r1 = sqrt(z1 / cap->a);

	// ângulo de cada pastilha
	divd = cap->d / 2.;
	alpha0 = 2 * asin(divd / r0);
	alpha1 = 2 * asin(divd / r1);

	// ângulo da diferença entre cada pastilha
	beta0 = (two_pi - (n_tiles * alpha0)) / n_tiles;
	beta1 = (two_pi - (n_tiles * alpha1)) / n_tiles;

	// gera cada uma das pastilhas
	for (i = 0; i < n_tiles; i++) {
		x0 = r0 * cos( i * (alpha0 + beta0) );
		y0 = r0 * sin( i * (alpha0 + beta0) );
	
		X0 = r0 * cos( (i+1) * (alpha0 + beta0) );
		Y0 = r0 * sin( (i+1) * (alpha0 + beta0) );

		x1 = r1 * cos( i * (alpha1 + beta1) );
		y1 = r1 * sin( i * (alpha1 + beta1) );

		X1 = r1 * cos( (i+1) * (alpha1 + beta1) );
		Y1 = r1 * sin( (i+1) * (alpha1 + beta1) );

		v3d_set(&a, x0, y0, z0);
		v3d_set(&b, X0, Y0, z0);
		v3d_set(&c, x1, y1, z1);
		v3d_set(&d, X1, Y1, z1);

		tile_init(&ring->tiles[i], cap, &a, &b, &c, &d, ring);
	}

	for (i = 0; i < n_tiles; i++) {
		next = (i + 1) % n_tiles;
		prev = (i - 1 + n_tiles) % n_tiles;

		tile_link(&ring->tiles[i], &ring->tiles[prev], &ring->tiles[next]);
	}

	ring->temp = cap->theta_0;
	ring->next_ring = NULL;
	ring->prev_ring = NULL;
}


void ring_free(Ring *ring) {
	free(ring->tiles);
}


void ring_calc_temp(Ring *ring, const int step_mod_2) {
	unsigned int i, n_tiles;

	#ifdef DEBUG
	assert(ring != NULL);
	#endif
	
	n_tiles = ring->n_tiles;
	
	#pragma omp parallel for \
		private (i) \
		shared (n_tiles, ring) \
		num_threads(NUM_THREADS) \
		if (n_tiles > MIN_TILES_TO_PARALLEL)
	for (i = 0; i < n_tiles; i++) {
		tile_calc_temp(&ring->tiles[i], step_mod_2);
	}
}

// 
void ring_update_temp(Ring *ring, const int step_mod_2) {
	double s;
	unsigned int i, n_tiles;

	n_tiles = ring->n_tiles;

	s = 0.;
	for (i = 0; i < n_tiles; i++) {
		s += ring->tiles[i].temp[step_mod_2];
	}

	ring->temp = s / ((double) n_tiles);
}

void ring_neighborhood_temp(Ring *ring, double *t1, double *t2) {
	#ifdef DEBUG
	assert(ring != NULL);
	assert(t1 != NULL && t2 != NULL);
	assert(ring->next_ring != NULL);
	assert(ring->prev_ring != NULL);
	#endif

	*t1 = ring->next_ring->temp;
	*t2 = ring->prev_ring->temp;
}

double inline ring_neighborhood_temp_plus(Ring *ring) {
	#ifdef DEBUG
	assert(ring != NULL);
	assert(ring->next_ring != NULL);
	assert(ring->prev_ring != NULL);
	#endif

	return ring->next_ring->temp + ring->prev_ring->temp;
}


void ring_print(Ring *ring, FILE *file, const int step_mod_2) {
	/* ===== arquivo de saída =====
	[OK] Primeira linha: a, h e d.
	[OK] Segunda linha: temperatura da calota.
	[OK] Terceira linha: número de anéis.
	Linhas subsequentes: um anel por linha, com o número de pastilhas na primeira
		posição, seguido pelas temperaturas de cada uma. Suponha as pastilhas
		dispostas no sentido horário, quando se olha de frente para o vértice
		do paraboloide. Comece pelo topo do anel.
	
	Temperatura positiva indica presença, um valor negativo indica a temperatura
	do rejunte na posição (com o valor trocado). Valores em Celsius.
	*/

	unsigned int i, n_tiles;

	assert(ring != NULL);

	n_tiles = ring->n_tiles;

	fprintf(file, "%d", ring->n_tiles);

	for (i = 0; i < n_tiles; i++) {
		fprintf(file, " %.2lf",
			((ring->tiles[i].bursted ? -1 : 1) * ring->tiles[i].temp[step_mod_2]));
	}

	fprintf(file, "\n");
}

// END [RING]

// BEGIN [COVER]

void cover_init(Cover *c, Capsule *capsule) {

	assert(c != NULL);
	assert(capsule != NULL);

	c->ring.capsule = capsule;
	v3d_set(&c->normal,
	 -capsule->pos.x, -capsule->pos.y, -capsule->pos.z);

	c->ring.next_ring = NULL;
	c->ring.temp = capsule->theta_0;
	c->t = capsule->t_0;
	c->bursted = 0;
	c->last_temp = c->new_temp = c->ring.temp;
}

void cover_calc_temp(Cover *c) {

	#ifdef DEBUG
	assert(c != NULL);
	assert(c->ring.next_ring != NULL);
	#endif

	c->t += 1.;
	if (c->bursted) { // estourou?
		c->new_temp = c->ring.next_ring->temp;
	} else {
		double vn;
		double delta_atrito, delta_dissip;
		double val;

		vn = v3d_dot(&c->normal, &c->ring.capsule->vel);

		if (vn > 0.) {
			val = c->t - c->ring.capsule->t_0;
			delta_atrito = c->ring.capsule->alpha * vn * atan(val * val);
		} else {
			delta_atrito = 0.;
		}

		delta_dissip = c->ring.capsule->delta * abs(vn);
		c->new_temp = c->ring.next_ring->temp + delta_atrito - delta_dissip;

		if (c->new_temp > c->ring.capsule->theta_crit) {
			c->bursted = 1;
			c->new_temp = c->ring.next_ring->temp;
		}
	}
}

double cover_update_temp(Cover *c) {
	#ifdef DEBUG
	assert(c != NULL);
	#endif

	return (c->ring.temp = c->last_temp = c->new_temp);
}

void cover_print(Cover *c, FILE *file) {
	assert(c != NULL);

	fprintf(file, "%.2f\n",
		((c->bursted) ? -1 : 1 ) * c->last_temp);
}

// END [COVER]

// BEGIN [MESH]

void mesh_init(Mesh *m, Capsule *cap) {
	double L, l, tmp;
	Ring *prev_ring;
	unsigned int i;

	assert(m != NULL);
	assert(cap != NULL);

	tmp = 3 * (cap->d / MM_PI);
	L = cap->a * ( tmp * tmp );
	l = L;

	m->cap = cap;

	m->n_rings = 0;
	cover_init(&m->cover, cap);
	while ((l + L) < cap->h) {
		m->n_rings++;
		l += L;
	}

	m->rings = (Ring*) malloc(sizeof(Ring) * m->n_rings);

	i = 0;
	l = L;

	prev_ring = (Ring*) (&m->cover);

	while ((l + L) < cap->h) {
		ring_init(&m->rings[i], cap, l, L);
		m->rings[i].prev_ring = prev_ring;
		if (i < (m->n_rings-1)) {
			m->rings[i].next_ring = &m->rings[i+1];
		}
		prev_ring = &m->rings[i];
		l += L;
		i++;
	}

	m->rings[m->n_rings-1].next_ring = &m->rings[m->n_rings-1];
	m->cover.ring.next_ring = &m->rings[0];

}


void mesh_free(Mesh *m) {
	unsigned int i, n_rings;
	n_rings = m->n_rings;
	for (i=0; i<n_rings; i++) {
		ring_free(&m->rings[i]);
	}
	free(m->rings);
}


void mesh_print(Mesh *m, FILE* file, const int step_mod_2) {
	/* ===== arquivo de saída =====
	[OK] Primeira linha: a, h e d.
	Segunda linha: temperatura da calota.
	Terceira linha: número de anéis.
	Linhas subsequentes: um anel por linha, com o número de pastilhas na primeira
		posição, seguido pelas temperaturas de cada uma. Suponha as pastilhas
		dispostas no sentido horário, quando se olha de frente para o vértice
		do paraboloide. Comece pelo topo do anel.
	
	Temperatura positiva indica presença, um valor negativo indica a temperatura
	do rejunte na posição (com o valor trocado). Valores em Celsius.
	*/

	unsigned int i;

	cover_print(&m->cover, file);

	fprintf(file, "%u\n", m->n_rings);

	for (i = 0; i < m->n_rings; i++) {
		ring_print(&m->rings[i], file, step_mod_2);
	}
}

void mesh_step(Mesh *m, const int step_mod_2) {
	unsigned int i, n_rings;

	n_rings = m->n_rings;	
	
	#pragma omp parallel for \
		private (i) \
		shared (n_rings, m) \
		num_threads(NUM_THREADS)
	for (i = 0; i < n_rings; i++) {
		ring_calc_temp(&m->rings[i], step_mod_2);
	}
	cover_calc_temp(&m->cover);

	#pragma omp parallel for \
		private (i) \
		shared (n_rings, m) \
		num_threads(NUM_THREADS)
	for (i = 0; i < n_rings; i++) {
		ring_update_temp(&m->rings[i], step_mod_2);
	}
	cover_update_temp(&m->cover);
}

double mesh_temp_media_pastilhas(Mesh *m, const int step_mod_2) {
	unsigned int i, j, qtd_pastilhas = 0, n_rings;
	double soma = 0;

	n_rings = m->n_rings;

	if (!m->cover.bursted) {
		qtd_pastilhas++;
		soma += m->cover.last_temp;
	}

	for (i=0; i < n_rings; i++) {
		for (j=0; j < m->rings[i].n_tiles; j++) {
			if (!m->rings[i].tiles[j].bursted) {
				soma += m->rings[i].tiles[j].temp[step_mod_2];
				qtd_pastilhas++;
			}
		}
	}

	if (qtd_pastilhas == 0) {
		return 0;
	}
	return soma / qtd_pastilhas;
}

double mesh_temp_media_rejunte(Mesh *m) {
	unsigned int i, n_rings;
	double soma = 0;

	n_rings = m->n_rings;

	for (i=0; i < n_rings; i++) {
		soma += m->rings[i].temp;
	}
	return soma / n_rings;
}

// END [MESH]

void capsule_print_params(Capsule *c) {

	assert(c != NULL);

	printf("Parâmetros da cápsula:\n");

	printf("h = %lf \n" \
	 "a = %lf \n" \
	 "d = %lf \n" \
	 "alpha = %lf \n" \
	 "delta = %lf \n" \
	 "t_0 = %lf \n" \
	 "theta_crit = %lf \n" \
	 "theta_0 = %lf \n", 
	 c->h, c->a, c->d, c->alpha, c->delta,
	 c->t_0, c->theta_crit,
	 c->theta_0);

	printf("pos = ");
	v3d_print(&c->pos);
	printf("\n");

	printf("vel = ");
	v3d_print(&c->vel);
	printf("\n");

	printf("steps = %lu\n", c->steps);
}

void capsule_init(Capsule *capsule) {

/* início da rotação de vetores */
	double alfa, beta, x_linha, y_linha, z_linha;

	// rotacao em relacao ao eixo z (a fim de zerar x)
	alfa = (capsule->pos.x == 0)
		? M_PI / 2.0
		: (-1) * atan( capsule->pos.y / capsule->pos.x );

	// rotação do vetor velocidade
	x_linha = capsule->vel.x * cos(alfa) - capsule->vel.y * sin(alfa);
	y_linha = capsule->vel.x * sin(alfa) + capsule->vel.y * cos(alfa);

	capsule->vel.x = x_linha;
	capsule->vel.y = y_linha;

	// rotação do vetor posicao
	x_linha = capsule->pos.x * cos(alfa) - capsule->pos.y * sin(alfa);
	y_linha = capsule->pos.x * sin(alfa) + capsule->pos.y * cos(alfa);

	capsule->pos.x = x_linha;
	capsule->pos.y = y_linha;

	// rotacao em relacao ao eixo x (a fim de zerar y)
	beta = (capsule->pos.z == 0)
		? M_PI / 2.0
	 	:(-1) * atan( capsule->pos.x / capsule->pos.z );

	// rotação do vetor velocidade
	x_linha = capsule->vel.z * sin(beta) + capsule->vel.x * cos(beta);
	z_linha = capsule->vel.z * cos(beta) - capsule->vel.x * sin(beta);	

	capsule->vel.x = x_linha;
	capsule->vel.z = z_linha;	

	// rotação do vetor posicao
	x_linha = capsule->pos.z * sin(beta) + capsule->pos.x * cos(beta);
	z_linha = capsule->pos.z * cos(beta) - capsule->pos.x * sin(beta);	

	capsule->pos.x = x_linha;
	capsule->pos.z = z_linha;
	
	// se o vetor posição aponta para cima, inverte os dois vetores.
	if (capsule->pos.z > 0) {
		capsule->vel.x = -capsule->vel.x;
		capsule->vel.y = -capsule->vel.y;
		capsule->vel.z = -capsule->vel.z;

		capsule->pos.z = -capsule->pos.z;
	}
/* fim da rotação de vetores */

	mesh_init(&capsule->mesh, capsule);
}

void capsule_free(Capsule *capsule) {
	mesh_free(&capsule->mesh);
}

void capsule_iterate(Capsule *capsule) {
	unsigned long int i;

	for (i = 0; i < capsule->steps; i++) {
		mesh_step(&capsule->mesh, i%2);
	}
}

void capsule_output(Capsule *capsule, const char* filename) {
	FILE *file;
	file = fopen(filename, "w+");

	/* ===== arquivo de saída =====
	Primeira linha: a, h e d.
	Segunda linha: temperatura da calota.
	Terceira linha: número de anéis.
	Linhas subsequentes: um anel por linha, com o número de pastilhas na primeira
		posição, seguido pelas temperaturas de cada uma. Suponha as pastilhas
		dispostas no sentido horário, quando se olha de frente para o vértice
		do paraboloide. Comece pelo topo do anel.
	
	Temperatura positiva indica presença, um valor negativo indica a temperatura
	do rejunte na posição (com o valor trocado). Valores em Celsius.
	*/

	fprintf(file, "%f %f %f\n", capsule->a, capsule->h, capsule->d);

	mesh_print(&capsule->mesh, file, (capsule->steps-1)%2);

	fclose(file);

	printf("Temperatura média das pastilhas: %f\nTemperatura média do rejunte: %f\n",
		mesh_temp_media_pastilhas(&capsule->mesh, (capsule->steps-1)%2),
		mesh_temp_media_rejunte(&capsule->mesh));
}

