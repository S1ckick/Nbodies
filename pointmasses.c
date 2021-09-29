#include "pointmasses.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int earthNum = 3;
const int moonNum = 10;
const int barrier = 16;
const int n_bodies_with_partials = 10;

#define SQR(x) ((x) * (x))

//#define PEPDEBUG
//#define NO_RELATIVITY

#ifdef __MINGW32__
#include <fenv.h>
#define printf __mingw_printf
#endif

typedef struct {
  int offset;
  int weighted_sum_offset;
} MassPartialDesc;

typedef struct {
  DOUBLE *masses;
  int n_objects;

  DOUBLE *dx, *dy, *dz;
  DOUBLE *dist, *dist2, *dist3;
  DOUBLE *temp_acc;
  MassPartialDesc *mass_partials;

  double mass0_initial;
  double mass0_trend;

  int ab_start, ab_end;

} PointMasses;

void pointmassesDestroy(void *pm) {
  PointMasses *self = (PointMasses *)pm;
  if (pm == NULL) return;

  free(self->masses);
  free(self->dx);
  free(self->dy);
  free(self->dz);
  free(self->dist);
  free(self->dist2);
  free(self->dist3);
  free(self->temp_acc);
  free(self->mass_partials);

  free(pm);
}

int _allocatePM(PointMasses *pm, int n_objects) {
#define ALLOC(name, size) \
  if ((pm->name = (DOUBLE *)malloc(sizeof(DOUBLE) * size)) == NULL) return 0
  ALLOC(masses, n_objects);
  ALLOC(dx, barrier * barrier);
  ALLOC(dy, barrier * barrier);
  ALLOC(dz, barrier * barrier);
  ALLOC(dist, barrier * barrier);
  ALLOC(dist2, barrier * barrier);
  ALLOC(dist3, barrier * barrier);
  ALLOC(temp_acc, n_objects * 3);
#undef ALLOC

  if ((pm->mass_partials = malloc(sizeof(MassPartialDesc) * n_objects *
                                  n_bodies_with_partials)) == NULL)
    return 0;

  return 1;
}

void *pointmassesCreate(int n_objects, double *masses) {
  int i;
  PointMasses *pm = (PointMasses *)malloc(sizeof(PointMasses));
  if (pm == NULL) return NULL;

  memset(pm, 0, sizeof(PointMasses));

  if (!_allocatePM(pm, n_objects)) {
    pointmassesDestroy(pm);
    return NULL;
  }

  for (i = 0; i < n_objects; i++) pm->masses[i] = masses[i];

  pm->n_objects = n_objects;
  memset(pm->mass_partials, 0,
         n_objects * n_bodies_with_partials * sizeof(MassPartialDesc));

#ifdef PEPDEBUG
  pm->masses[0] = 0.00029591220828559110253e0L;  // sun
  pm->masses[1] = 4.9124804503647600162e-11L;
  pm->masses[2] = 7.2434523326441197223e-10L;
  pm->masses[3] = 8.8876924488373404697e-10L;  // earth
  pm->masses[4] = 9.549548695550770507e-11L;
  pm->masses[5] = 2.8253458259715502583e-07L;
  pm->masses[6] = 8.4597060732450306266e-08L;
  pm->masses[7] = 1.2920265796290799292e-08L;
  pm->masses[8] = 1.5243573478851101461e-08L;
  pm->masses[9] = 2.1750990867748499446e-12L;
  pm->masses[10] = 1.0931894537290947847e-11L;  // moon
  for (int j = 0; j < 11; j++) printf("GM[%2d] = % .20Le\n", j, pm->masses[j]);
#endif

  pm->mass0_initial = pm->masses[0];

  return pm;
}

int pointmassesSetMassTrend(void *pm, int body, double dmass_dt) {
  PointMasses *self = (PointMasses *)pm;

  if (body != 0) return -1;

  self->mass0_trend = dmass_dt;

  return 0;
}

void pointmassesSetAB(void *pm, int ab_start, int ab_end) {
  PointMasses *self = (PointMasses *)pm;

  self->ab_start = ab_start;
  self->ab_end = ab_end;
}

void pointmassesCalculateXdot(DOUBLE *x, double t, DOUBLE *f, void *userdata) {
  PointMasses *self = (PointMasses *)userdata;
  DOUBLE *masses = self->masses;
  int n_objects = self->n_objects;
  DOUBLE *dist = self->dist, *dist2 = self->dist2, *dist3 = self->dist3;
  DOUBLE *dx = self->dx, *dy = self->dy, *dz = self->dz;
  DOUBLE *temp_acc = self->temp_acc;
  int i, j, k;

#ifdef __MINGW32__
  fenv_t fenv;
  fegetenv(&fenv);
  fesetenv(FE_PC64_ENV);
#endif

  // Copy velocities of bodies
  for (i = 0; i < self->n_objects; i++) {
    f[6 * i] = x[6 * i + 3];
    f[6 * i + 1] = x[6 * i + 4];
    f[6 * i + 2] = x[6 * i + 5];
  }

  // Account for Sun's mass trend
  masses[0] = self->mass0_initial + self->mass0_trend * t;

  // Convert positions and velocities from EMB to barycentric in case of Earth
  // and Moon
  DOUBLE earth_x = x[6 * earthNum], earth_y = x[6 * earthNum + 1],
         earth_z = x[6 * earthNum + 2];
  DOUBLE earth_vx = x[6 * earthNum + 3], earth_vy = x[6 * earthNum + 4],
         earth_vz = x[6 * earthNum + 5];
  DOUBLE gcmoon_x = x[6 * moonNum], gcmoon_y = x[6 * moonNum + 1],
         gcmoon_z = x[6 * moonNum + 2];
  DOUBLE gcmoon_vx = x[6 * moonNum + 3], gcmoon_vy = x[6 * moonNum + 4],
         gcmoon_vz = x[6 * moonNum + 5];

  x[6 * moonNum] = earth_x + gcmoon_x;
  x[6 * moonNum + 1] = earth_y + gcmoon_y;
  x[6 * moonNum + 2] = earth_z + gcmoon_z;
  x[6 * moonNum + 3] = earth_vx + gcmoon_vx;
  x[6 * moonNum + 4] = earth_vy + gcmoon_vy;
  x[6 * moonNum + 5] = earth_vz + gcmoon_vz;

  // Precalculate dx, dy, dz, distances
  for (i = 0; i < barrier; i++) {
    DOUBLE *i_dx_ptr = dx + i * barrier;
    DOUBLE *i_dy_ptr = dy + i * barrier;
    DOUBLE *i_dz_ptr = dz + i * barrier;
    DOUBLE *i_dist_ptr = dist + i * barrier;
    DOUBLE *i_dist2_ptr = dist2 + i * barrier;
    DOUBLE *i_dist3_ptr = dist3 + i * barrier;

    for (j = i + 1; j < barrier; j++) {
      DOUBLE _dx = x[6 * j] - x[6 * i], _dy = x[6 * j + 1] - x[6 * i + 1],
             _dz = x[6 * j + 2] - x[6 * i + 2];
      DOUBLE _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
      DOUBLE _dist = SQRT(_dist2);
#ifdef __MINGW32__
      fesetenv(FE_PC64_ENV);
#endif

      DOUBLE _dist3 = _dist * _dist2;

      i_dx_ptr[j] = _dx;
      i_dy_ptr[j] = _dy;
      i_dz_ptr[j] = _dz;
      dx[j * barrier + i] = -_dx;
      dy[j * barrier + i] = -_dy;
      dz[j * barrier + i] = -_dz;

      i_dist_ptr[j] = _dist;
      i_dist2_ptr[j] = _dist2;
      i_dist3_ptr[j] = _dist3;

      dist[j * barrier + i] = _dist;
      dist2[j * barrier + i] = _dist2;
      dist3[j * barrier + i] = _dist3;
    }
  }
  // Go backwards so that the acceleration from the Sun is added last
  for (i = barrier - 1; i >= 0; i--) {
    DOUBLE *i_dx_ptr = dx + i * barrier;
    DOUBLE *i_dy_ptr = dy + i * barrier;
    DOUBLE *i_dz_ptr = dz + i * barrier;
    DOUBLE *i_dist3_ptr = dist3 + i * barrier;

    for (j = i + 1; j < barrier; j++) {
      DOUBLE _dx = i_dx_ptr[j], _dy = i_dy_ptr[j], _dz = i_dz_ptr[j];
      DOUBLE _dist3 = i_dist3_ptr[j];

      DOUBLE k_dx = _dx * _dist3;
      DOUBLE k_dy = _dy * _dist3;
      DOUBLE k_dz = _dz * _dist3;

      f[6 * i + 3] += masses[j] * k_dx;
      f[6 * i + 4] += masses[j] * k_dy;
      f[6 * i + 5] += masses[j] * k_dz;
      f[6 * j + 3] -= masses[i] * k_dx;
      f[6 * j + 4] -= masses[i] * k_dy;
      f[6 * j + 5] -= masses[i] * k_dz;
    }

    // simplified continuation of the above loop
    for (j = barrier; j < n_objects; j++) {
      // Five individual asteroids must not perturb the discrete asteroid ring
      if ((j >= self->ab_start) && (j < self->ab_end) && (i > moonNum))
        continue;

      DOUBLE _dx = x[6 * j] - x[6 * i], _dy = x[6 * j + 1] - x[6 * i + 1],
             _dz = x[6 * j + 2] - x[6 * i + 2];
      DOUBLE _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
      DOUBLE _dist = SQRT(_dist2);
#ifdef __MINGW32__
      fesetenv(FE_PC64_ENV);
#endif

      DOUBLE _dist3 = _dist * _dist2;

      DOUBLE k_dx = _dx * _dist3;
      DOUBLE k_dy = _dy * _dist3;
      DOUBLE k_dz = _dz * _dist3;

      f[6 * i + 3] += masses[j] * k_dx;
      f[6 * i + 4] += masses[j] * k_dy;
      f[6 * i + 5] += masses[j] * k_dz;

      f[6 * j + 3] -= masses[i] * k_dx;
      f[6 * j + 4] -= masses[i] * k_dy;
      f[6 * j + 5] -= masses[i] * k_dz;
    }
  }

  // Relativistic forces computation
#ifndef NO_RELATIVITY
  for (i = 0; i < barrier; ++i) {
    temp_acc[3 * i] = f[6 * i + 3];
    temp_acc[3 * i + 1] = f[6 * i + 4];
    temp_acc[3 * i + 2] = f[6 * i + 5];
  }

  for (i = 0; i < barrier; ++i) {
    DOUBLE vi_sqr = x[6 * i + 3] * x[6 * i + 3] + x[6 * i + 4] * x[6 * i + 4] +
                    x[6 * i + 5] * x[6 * i + 5];

    DOUBLE sum_i = 0.0;
    for (k = 0; k < barrier; ++k) {
      if (i == k) continue;
      sum_i += masses[k] * dist[barrier * i + k];
    }

    for (j = 0; j < barrier; ++j) {
      if (i == j) continue;

      int ij = barrier * i + j;

      DOUBLE vj_sqr = x[6 * j + 3] * x[6 * j + 3] +
                      x[6 * j + 4] * x[6 * j + 4] + x[6 * j + 5] * x[6 * j + 5];

      DOUBLE vi_dot_vj = x[6 * i + 3] * x[6 * j + 3] +
                         x[6 * i + 4] * x[6 * j + 4] +
                         x[6 * i + 5] * x[6 * j + 5];

      DOUBLE sum_j = 0.0;
      for (k = 0; k < barrier; ++k) {
        if (j == k) continue;
        sum_j += masses[k] * dist[barrier * j + k];
      }

      DOUBLE t1 =
          -4.0 * sum_i - sum_j + vi_sqr + 2.0 * vj_sqr - 4.0 * vi_dot_vj;
      DOUBLE t2 = dist[ij] * (-dx[ij] * x[6 * j + 3] + -dy[ij] * x[6 * j + 4] +
                              -dz[ij] * x[6 * j + 5]);
      t2 = -1.5 * t2 * t2;
      DOUBLE t3 =
          0.5 * (dx[ij] * temp_acc[3 * j] + dy[ij] * temp_acc[3 * j + 1] +
                 dz[ij] * temp_acc[3 * j + 2]);

      DOUBLE c1 = revLS2 * masses[j] * dist3[ij] * (t1 + t2 + t3);
      f[6 * i + 3] += dx[ij] * c1;
      f[6 * i + 4] += dy[ij] * c1;
      f[6 * i + 5] += dz[ij] * c1;

      DOUBLE p1 = dist[ij] * (-dx[ij] * (4 * x[6 * i + 3] - 3 * x[6 * j + 3]) +
                              -dy[ij] * (4 * x[6 * i + 4] - 3 * x[6 * j + 4]) +
                              -dz[ij] * (4 * x[6 * i + 5] - 3 * x[6 * j + 5]));

      DOUBLE c2 = revLS2 * masses[j] * dist2[ij] * p1;

      DOUBLE dv_x = x[6 * i + 3] - x[6 * j + 3];
      DOUBLE dv_y = x[6 * i + 4] - x[6 * j + 4];
      DOUBLE dv_z = x[6 * i + 5] - x[6 * j + 5];
      f[6 * i + 3] += dv_x * c2;
      f[6 * i + 4] += dv_y * c2;
      f[6 * i + 5] += dv_z * c2;

      DOUBLE c3 = 3.5 * revLS2 * masses[j] * dist[ij];
      f[6 * i + 3] += temp_acc[3 * j] * c3;
      f[6 * i + 4] += temp_acc[3 * j + 1] * c3;
      f[6 * i + 5] += temp_acc[3 * j + 2] * c3;
    }
  }
#endif

  // Recover state of geocentric Moon and convert acceleration
  x[6 * moonNum] = gcmoon_x;
  x[6 * moonNum + 1] = gcmoon_y;
  x[6 * moonNum + 2] = gcmoon_z;
  x[6 * moonNum + 3] = gcmoon_vx;
  x[6 * moonNum + 4] = gcmoon_vy;
  x[6 * moonNum + 5] = gcmoon_vz;

  f[6 * moonNum + 3] -= f[6 * earthNum + 3];
  f[6 * moonNum + 4] -= f[6 * earthNum + 4];
  f[6 * moonNum + 5] -= f[6 * earthNum + 5];

#ifdef __MINGW32__
  fesetenv(&fenv);
#endif
}
