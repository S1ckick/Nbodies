#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>

#include "../def.h"

long double to_double(const long double &x)
{
  return x;
}

using helper_type = double;

template <typename Type>
class ObjectsData
{
public:
  ObjectsData(std::vector<double> &data)
  {
    this->n_objects = data.size();
    masses.resize(this->n_objects);
    for (int i = 0; i < data.size(); i++)
    {
      this->masses[i] = data[i];
    }
    this->dx.resize(this->n_objects * (this->n_objects + 1));
    this->dy.resize(this->n_objects * (this->n_objects + 1));
    this->dz.resize(this->n_objects * (this->n_objects + 1));
    this->dist.resize(this->n_objects * (this->n_objects + 1));
    this->dist2.resize(this->n_objects * (this->n_objects + 1));
    this->dist3.resize(this->n_objects * (this->n_objects + 1));

    this->fx.resize(this->n_objects, 0);
    this->fy.resize(this->n_objects, 0);
    this->fz.resize(this->n_objects, 0);
    this->temp_acc.resize(this->n_objects * 3);
  }

  std::vector<double> masses;
  int n_objects;
  std::vector<helper_type> dx, dy, dz;
  std::vector<helper_type> dist, dist2, dist3;
  std::vector<Type> temp_acc;
  std::vector<helper_type> fx, fy, fz;
};
#define NO_RELATIVITY

template <typename ABMD_DOUBLE>
struct ContextData
{
  ObjectsData<ABMD_DOUBLE> *objects;
  ABMD_DOUBLE *sol;
  ABMD_DOUBLE *sol_back;
  double *callback_t;
  ABMD_DOUBLE *energy;
  ABMD_DOUBLE *impulse;
  ABMD_DOUBLE *center;
  ABMD_DOUBLE init_energy;
  ABMD_DOUBLE *init_impulse;
  ABMD_DOUBLE *init_center;
  int i;
  int dim;
  std::ofstream f;
  std::ofstream fb;
};

template <typename ABMD_DOUBLE>
void pointmassesCalculateXdot_tmp(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *f, void *context)
{

  int i, j, k;
  //userdata()
  // Copy velocities of bodies
  ObjectsData<ABMD_DOUBLE> *userdata = static_cast<ContextData<ABMD_DOUBLE> *>(context)->objects;
  memset(&userdata->fx[0], 0, sizeof(helper_type) * userdata->n_objects);
  memset(&userdata->fy[0], 0, sizeof(helper_type) * userdata->n_objects);
  memset(&userdata->fz[0], 0, sizeof(helper_type) * userdata->n_objects);
  for (i = 0; i < userdata->n_objects; i++)
  {
    f[6 * i] = x[6 * i + 3];
    f[6 * i + 1] = x[6 * i + 4];
    f[6 * i + 2] = x[6 * i + 5];
  }

  // Convert positions and velocities from EMB to barycentric in case of Earth
  // and Moon
  ABMD_DOUBLE earth_x = x[6 * earthNum], earth_y = x[6 * earthNum + 1],
              earth_z = x[6 * earthNum + 2];
  ABMD_DOUBLE earth_vx = x[6 * earthNum + 3], earth_vy = x[6 * earthNum + 4],
              earth_vz = x[6 * earthNum + 5];
  ABMD_DOUBLE gcmoon_x = x[6 * moonNum], gcmoon_y = x[6 * moonNum + 1],
              gcmoon_z = x[6 * moonNum + 2];
  ABMD_DOUBLE gcmoon_vx = x[6 * moonNum + 3], gcmoon_vy = x[6 * moonNum + 4],
              gcmoon_vz = x[6 * moonNum + 5];

  x[6 * moonNum] = earth_x + gcmoon_x;
  x[6 * moonNum + 1] = earth_y + gcmoon_y;
  x[6 * moonNum + 2] = earth_z + gcmoon_z;
  x[6 * moonNum + 3] = earth_vx + gcmoon_vx;
  x[6 * moonNum + 4] = earth_vy + gcmoon_vy;
  x[6 * moonNum + 5] = earth_vz + gcmoon_vz;

  ABMD_DOUBLE _dx_me, _dy_me, _dz_me, _dist2_me, _dist_me, _dist3_me;

  // Precalculate dx, dy, dz, distances
  for (i = 0; i < barrier; i++)
  {
    for (j = i + 1; j < barrier; j++)
    {
      helper_type _dx, _dy, _dz, _dist2, _dist, _dist3;

      if (i == 3 && j == 10)
      {
        _dx_me = gcmoon_x;
        _dy_me = gcmoon_y;
        _dz_me = gcmoon_z;
        _dist2_me = 1.0 / (_dx_me * _dx_me + _dy_me * _dy_me + _dz_me * _dz_me);
        _dist_me = sqrt(_dist2_me);
        _dist3_me = _dist_me * _dist2_me;

        userdata->dx[i * barrier + j] = to_double(_dx_me);
        userdata->dy[i * barrier + j] = to_double(_dy_me);
        userdata->dz[i * barrier + j] = to_double(_dz_me);

        userdata->dx[j * barrier + i] = -to_double(_dx_me);
        userdata->dy[j * barrier + i] = -to_double(_dy_me);
        userdata->dz[j * barrier + i] = -to_double(_dz_me);
      }
      else
      {
        _dx = to_double(x[6 * j]) - to_double(x[6 * i]);
        _dy = to_double(x[6 * j + 1]) - to_double(x[6 * i + 1]);
        _dz = to_double(x[6 * j + 2]) - to_double(x[6 * i + 2]);
        _dist2 = 1.0 / ((_dx * _dx) + (_dy * _dy) + (_dz * _dz));
        _dist = sqrt(_dist2);
        _dist3 = _dist * _dist2;

        userdata->dx[i * barrier + j] = _dx;
        userdata->dy[i * barrier + j] = _dy;
        userdata->dz[i * barrier + j] = _dz;

        userdata->dx[j * barrier + i] = -_dx;
        userdata->dy[j * barrier + i] = -_dy;
        userdata->dz[j * barrier + i] = -_dz;

        userdata->dist[i * barrier + j] = _dist;
        userdata->dist2[i * barrier + j] = _dist2;
        userdata->dist3[i * barrier + j] = _dist3;

        userdata->dist[j * barrier + i] = _dist;
        userdata->dist2[j * barrier + i] = _dist2;
        userdata->dist3[j * barrier + i] = _dist3;
      }
    }
  }

  // Go backwards so that the acceleration from the Sun is added last
  for (i = barrier - 1; i >= 0; i--)
  {
    for (j = i + 1; j < barrier; j++)
    {
      if (i == earthNum && j == moonNum)
      {
        // ABMD_DOUBLE k_dx = _dx_me * _dist3_me;
        // ABMD_DOUBLE k_dy = _dy_me * _dist3_me;
        // ABMD_DOUBLE k_dz = _dz_me * _dist3_me;

        // userdata->fx[i] += to_double(userdata->masses[j] * k_dx);
        // userdata->fy[i] += to_double(userdata->masses[j] * k_dy);
        // userdata->fz[i] += to_double(userdata->masses[j] * k_dz);
        // userdata->fx[j] -= to_double(userdata->masses[i] * k_dx);
        // userdata->fy[j] -= to_double(userdata->masses[i] * k_dy);
        // userdata->fz[j] -= to_double(userdata->masses[i] * k_dz);

        // f[6 * i + 3] += userdata->masses[j] * k_dx;
        // f[6 * i + 4] += userdata->masses[j] * k_dy;
        // f[6 * i + 5] += userdata->masses[j] * k_dz;
        // f[6 * j + 3] -= userdata->masses[i] * k_dx;
        // f[6 * j + 4] -= userdata->masses[i] * k_dy;
        // f[6 * j + 5] -= userdata->masses[i] * k_dz;
      }
      else
      {
        helper_type _dx = userdata->dx[i * barrier + j],
                   _dy = userdata->dy[i * barrier + j],
                   _dz = userdata->dz[i * barrier + j];
        helper_type _dist3 = userdata->dist3[i * barrier + j];

        helper_type k_dx = _dx * _dist3;
        helper_type k_dy = _dy * _dist3;
        helper_type k_dz = _dz * _dist3;
        
        userdata->fx[i] += userdata->masses[j] * k_dx;
        userdata->fy[i] += userdata->masses[j] * k_dy;
        userdata->fz[i] += userdata->masses[j] * k_dz;
        userdata->fx[j] -= userdata->masses[i] * k_dx;
        userdata->fy[j] -= userdata->masses[i] * k_dy;
        userdata->fz[j] -= userdata->masses[i] * k_dz;
      }
    }
  }

  for (int i = 0; i < barrier; i++)
  {
    f[6 * i + 3] = userdata->fx[i];
    f[6 * i + 4] = userdata->fy[i];
    f[6 * i + 5] = userdata->fz[i];
  }

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

  ABMD_DOUBLE k_dx = _dx_me * _dist3_me;
  ABMD_DOUBLE k_dy = _dy_me * _dist3_me;
  ABMD_DOUBLE k_dz = _dz_me * _dist3_me;

  f[6 * earthNum + 3] += to_double(userdata->masses[moonNum] * k_dx);
  f[6 * earthNum + 4] += to_double(userdata->masses[moonNum] * k_dy);
  f[6 * earthNum + 5] += to_double(userdata->masses[moonNum] * k_dz);
  f[6 * moonNum + 3] -= to_double((userdata->masses[earthNum] + userdata->masses[moonNum]) * k_dx);
  f[6 * moonNum + 4] -= to_double((userdata->masses[earthNum] + userdata->masses[moonNum]) * k_dy);
  f[6 * moonNum + 5] -= to_double((userdata->masses[earthNum] + userdata->masses[moonNum]) * k_dz);

  for (i = 0; i < userdata->n_objects; i++)
  {
    f[6 * i] = to_double(f[6 * i]);
    f[6 * i + 1] = to_double(f[6 * i + 1]);
    f[6 * i + 2] = to_double(f[6 * i + 2]);

    f[6 * i + 3] = to_double(f[6 * i + 3]);
    f[6 * i + 4] = to_double(f[6 * i + 4]);
    f[6 * i + 5] = to_double(f[6 * i + 5]);
  }

  double revLS2 = pow(1731456.84 / 299792458, 2);
  // RELATIVISTIC
  for (i = 0; i < barrier; ++i)
  {
    // луна должна быть в барицентре прибавить ускорение земли если луна
    userdata->temp_acc[3 * i] = to_double(f[6 * i + 3]);
    userdata->temp_acc[3 * i + 1] = to_double(f[6 * i + 4]);
    userdata->temp_acc[3 * i + 2] = to_double(f[6 * i + 5]);

    if (i == 10)
    {
      userdata->temp_acc[3 * i] += to_double(f[6 * 3 + 3]);
      userdata->temp_acc[3 * i + 1] += to_double(f[6 * 3 + 4]);
      userdata->temp_acc[3 * i + 2] += to_double(f[6 * 3 + 5]);
    }
  }

  for (i = 0; i < barrier; ++i)
  {
    double vi_sqr = to_double(x[6 * i + 3]) * to_double(x[6 * i + 3]) +
                    to_double(x[6 * i + 4]) * to_double(x[6 * i + 4]) +
                    to_double(x[6 * i + 5]) * to_double(x[6 * i + 5]);

    double sum_i = 0.0;
    for (k = 0; k < barrier; ++k)
    {
      if (i == k)
        continue;
      sum_i += userdata->masses[k] * userdata->dist[barrier * i + k];
    }

    for (j = 0; j < barrier; ++j)
    {
      if (i == j)
        continue;

      int ij = barrier * i + j;

      double vj_sqr = to_double(x[6 * j + 3]) * to_double(x[6 * j + 3]) +
                    to_double(x[6 * j + 4]) * to_double(x[6 * j + 4]) +
                    to_double(x[6 * j + 5]) * to_double(x[6 * j + 5]);

      double vi_dot_vj = to_double(x[6 * i + 3]) * to_double(x[6 * j + 3]) +
                    to_double(x[6 * i + 4]) * to_double(x[6 * j + 4]) +
                    to_double(x[6 * i + 5]) * to_double(x[6 * j + 5]);

      double sum_j = 0.0;
      for (k = 0; k < barrier; ++k)
      {
        if (j == k)
          continue;
        sum_j += userdata->masses[k] * userdata->dist[barrier * j + k];
      }

      double t1 = -4.0 * sum_i - sum_j + vi_sqr + 2.0 * vj_sqr - 4.0 * vi_dot_vj;
      double t2 = userdata->dist[ij] * (-userdata->dx[ij] * to_double(x[6 * j + 3]) +
                                        -userdata->dy[ij] * to_double(x[6 * j + 4]) +
                                        -userdata->dz[ij] * to_double(x[6 * j + 5]));
      t2 = -1.5 * t2 * t2;
      double t3 = 0.5 * to_double(userdata->dx[ij] * userdata->temp_acc[3 * j] +
                         userdata->dy[ij] * userdata->temp_acc[3 * j + 1] +
                         userdata->dz[ij] * userdata->temp_acc[3 * j + 2]);

      double c1 = revLS2 * userdata->masses[j] * userdata->dist3[ij] * (t1 + t2 + t3);
      f[6 * i + 3] += userdata->dx[ij] * c1;
      f[6 * i + 4] += userdata->dy[ij] * c1;
      f[6 * i + 5] += userdata->dz[ij] * c1;

      double p1 = userdata->dist[ij] * (-userdata->dx[ij] * (4 * to_double(x[6 * i + 3] - 3 * x[6 * j + 3])) +
                                        -userdata->dy[ij] * (4 * to_double(x[6 * i + 4] - 3 * x[6 * j + 4])) +
                                        -userdata->dz[ij] * (4 * to_double(x[6 * i + 5] - 3 * x[6 * j + 5])));

      double c2 = revLS2 * userdata->masses[j] * userdata->dist2[ij] * p1;

      double dv_x = to_double(x[6 * i + 3] - x[6 * j + 3]);
      double dv_y = to_double(x[6 * i + 4] - x[6 * j + 4]);
      double dv_z = to_double(x[6 * i + 5] - x[6 * j + 5]);
      f[6 * i + 3] += dv_x * c2;
      f[6 * i + 4] += dv_y * c2;
      f[6 * i + 5] += dv_z * c2;

      double c3 = 3.5 * revLS2 * userdata->masses[j] * userdata->dist[ij];
      f[6 * i + 3] += userdata->temp_acc[3 * j] * c3;
      f[6 * i + 4] += userdata->temp_acc[3 * j + 1] * c3;
      f[6 * i + 5] += userdata->temp_acc[3 * j + 2] * c3;
    }
  }
}