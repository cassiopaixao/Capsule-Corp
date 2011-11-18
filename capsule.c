
#include <stdio.h>
#include <stdlib.h>

#include "capsule.h"

#include <math.h>
#include <assert.h>

void tile_init(Tile *t, Capsule *cap, v3d *a, v3d *b, v3d *c, v3d *d, Ring *ring) {
	v3d vdl, p, q;

	assert(t != NULL);
	assert(cap != NULL);
	assert(a != NULL);
	assert(b != NULL);
	assert(c != NULL);
	assert(d != NULL);

	t->t_0 = cap->t_0;
	t->t = t->t_0;
	t->alpha = cap->alpha;
	t->delta = cap->delta;
	t->last_temp = cap->theta_0;
	t->max_temp = cap->theta_crit;
	
	v3d_set(&t->vel, cap->vel.x, cap->vel.y, cap->vel.z);
	v3d_set(&t->pos, cap->pos.x, cap->pos.y, cap->pos.z);
	
	v3d_minus(c, a, &vdl);
	t->dl = v3d_length(&vdl);
	t->d = cap->d;
	t->bursted = 0;

	v3d_copy(&t->edges[0], a);
	v3d_copy(&t->edges[1], b);
	v3d_copy(&t->edges[2], c);
	v3d_copy(&t->edges[3], d);

	t->ring = ring;

	v3d_minus(b, a, &p);
	v3d_minus(d, a, &q);

	v3d_cross(&p, &q, &t->normal);
}


void tile_link(Tile *t, Tile *left, Tile *right) {

	assert(t != NULL);
	assert(left != NULL);
	assert(right != NULL);

	t->left = left;
	t->right = right;
}

static double _tile_perimenter_temp(Tile *t) {
	double t1, t2;

	assert(t != NULL);
	assert(t->left != NULL);
	assert(t->right != NULL);

	ring_neighborhood_temp(t->ring, &t1, &t2);

	return ((t->left->last_temp * t->dl) +
	        (t->right->last_temp * t->dl) +
	        (t1 * t->d) +
	        (t2 * t->d)) / (t->dl * 2. + t->d * 2.);
}

void tile_calc_temp(Tile *t) {
	assert(t != NULL);
	assert(t->left != NULL);
	assert(t->right != NULL);

	t->t += 1;
	if (t->bursted) { // estourou?
		// só rejunte
		t->new_temp = _tile_perimenter_temp(t);
	} else {
		double vn, delta_atrito, delta_dissip;

		vn = v3d_dot(&t->normal, &t->vel);
		if (vn > 0.) {
			double val;
			
			val = t->t - t->t_0;
			delta_atrito = t->alpha * vn * atan(val*val);
		} else {
			delta_atrito = 0.;
		}
		
		delta_dissip = t->delta * abs(vn);
		t->new_temp = _tile_perimenter_temp(t) + delta_atrito - delta_dissip;

		if (t->new_temp > t->ring->capsule->theta_crit) {
			t->bursted = 1; // estourou
			t->new_temp = (t->left->last_temp + t->right->last_temp)/2.;
		}
	}
}

double tile_update_temp(Tile *t) {
	assert(t != NULL);

	t->last_temp = t->new_temp;
	return t->last_temp;
}

// BEGIN [RING]

unsigned int _ring_n_tiles(Ring *ring, double l) {
	static double pi = 0.;
	unsigned int pi_ok = 0;
	//FIXME: Calcula PI no código mais de uma vez
	if (!pi_ok) {
		pi = 4. * atan(1.);
		pi_ok = 1;
	}

	assert(ring != NULL);
	assert(ring->capsule != NULL);
	
	return ((unsigned int)
	 floor((2. * pi * sqrt(l / ring->capsule->a)) / ring->capsule->d));
}

void ring_init(Ring *ring, Capsule *cap, double l, double L) {
	unsigned int n_tiles;
	double z0, z1, r0, r1;
	double alpha0, alpha1, divd;
	double beta0, beta1;
	//FIXME: Calcula PI no código mais de uma vez
	static double pi = 0.;
	static double two_pi = 0.;
	unsigned int pi_ok = 0;
	unsigned int i;
	double x0, y0, X0, Y0;
	double x1, y1, X1, Y1;
	v3d a, b, c, d;
	unsigned int next, prev;

	if (!pi_ok) {
		pi = 4. * atan(1.);
		pi_ok = 1;
		two_pi = 2. * pi;
	}
	assert(ring != NULL);
	assert(cap != NULL);

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


void ring_calc_temp(Ring *ring) {
	unsigned int i;

	assert(ring != NULL);

	for (i = 0; i < ring->n_tiles; i++) {
		tile_calc_temp(&ring->tiles[i]);
	}
}

// 
void ring_update_temp(Ring *ring) {
	double s;
	unsigned int i;	

	s = 0.;
	for (i = 0; i < ring->n_tiles; i++) {
		s += tile_update_temp(&ring->tiles[i]);
	}

	ring->temp = s / ((double) ring->n_tiles);
}

void ring_neighborhood_temp(Ring *ring, double *t1, double *t2) {
	assert(ring != NULL);
	assert(t1 != NULL && t2 != NULL);
	assert(ring->next_ring != NULL);
	assert(ring->prev_ring != NULL);

	*t1 = ring->next_ring->temp;
	*t2 = ring->prev_ring->temp;
}

void ring_print(Ring *ring) {
	unsigned int i;

	assert(ring != NULL);

	printf("[nº of tiles = %d] ", ring->n_tiles);	

	for (i = 0; i < ring->n_tiles; i++) {
		printf("%.2lf[%d] ", ring->tiles[i].last_temp, ring->tiles[i].bursted);
	}

	printf("\n");
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
	c->last_temp = c->ring.temp;
}

void cover_calc_temp(Cover *c) {
	
	assert(c != NULL);
	assert(c->ring.next_ring != NULL);

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

// FIXME: temp não é atualizado
double cover_update_temp(Cover *c) {
	assert(c != NULL);

	return (c->ring.temp = c->last_temp = c->new_temp);
}

void cover_print(Cover *c) {
	assert(c != NULL);

	if (c->bursted) {
		printf("[bursted]");
	}

	printf("cover: %.2f\n", c->last_temp);
}

// END [COVER]

// BEGIN [MESH]

void mesh_init(Mesh *m, Capsule *cap) {
	double L, l, tmp;
	double pi = 4. * atan(1.);
	Ring *prev_ring;
	unsigned int i;

	assert(m != NULL);
	assert(cap != NULL);
	
	pi = 4. * atan(1.);

	tmp = 3 * (cap->d / pi);
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


void mesh_print(Mesh *m) {
	unsigned int i;

	cover_print(&m->cover);
	printf("number of rings: %u\n", m->n_rings);
	
	for (i = 0; i < m->n_rings; i++) {
		ring_print(&m->rings[i]);
		printf("\n");
	}
}

void mesh_step(Mesh *m) {
	unsigned int i;
	
	for (i = 0; i < m->n_rings; i++) {
		ring_calc_temp(&m->rings[i]);
	}
	cover_calc_temp(&m->cover);

	for (i = 0; i < m->n_rings; i++) {
		ring_update_temp(&m->rings[i]);
	}
	cover_update_temp(&m->cover);
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
	 "t_inicial = %lf \n" \
	 "theta_crit = %lf \n" \
	 "theta_0 = %lf \n", 
	 c->h, c->a, c->d, c->alpha, c->delta,
	 c->t_0, c->t_inicial, c->theta_crit,
	 c->theta_0);

	printf("pos = ");
	v3d_print(&c->pos);
	printf("\n");

	printf("vel = ");
	v3d_print(&c->vel);
	printf("\n");

	printf("steps = %u\n", c->steps);
}

void capsule_init(Capsule *capsule) {
	mesh_init(&capsule->mesh, capsule);
}

void capsule_iterate(Capsule *capsule) {
	unsigned int i;

	for (i = 0; i < capsule->steps; i++) {
		mesh_step(&capsule->mesh);
	}
}

void capsule_output(Capsule *capsule) {
	mesh_print(&capsule->mesh);
}

