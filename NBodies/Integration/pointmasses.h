#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>
const int earthNum = 3;
const int moonNum = 10;
const int barrier = 16;

long double to_double(const long double &x)
{
  return x;
}

using debug_type = double;

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
  std::vector<debug_type> dx, dy, dz;
  std::vector<debug_type> dist, dist2, dist3;
  std::vector<Type> temp_acc;
  std::vector<debug_type> fx, fy, fz;
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
  memset(&userdata->fx[0], 0, sizeof(debug_type) * userdata->n_objects);
  memset(&userdata->fy[0], 0, sizeof(debug_type) * userdata->n_objects);
  memset(&userdata->fz[0], 0, sizeof(debug_type) * userdata->n_objects);
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
      debug_type _dx, _dy, _dz, _dist2, _dist, _dist3;

      if (i == 3 && j == 10)
      {
        _dx_me = gcmoon_x;
        _dy_me = gcmoon_y;
        _dz_me = gcmoon_z;
        _dist2_me = 1.0 / (_dx_me * _dx_me + _dy_me * _dy_me + _dz_me * _dz_me);
        _dist_me = sqrt(_dist2_me);
        _dist3_me = _dist_me * _dist2_me;
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

        userdata->dist3[i * barrier + j] = _dist3;

        userdata->dist3[j * barrier + i] = _dist3;
      }
    }
  }

  // Go backwards so that the acceleration from the Sun is added last
  for (i = barrier - 1; i >= 0; i--)
  {
    for (j = i + 1; j < barrier; j++)
    {
      if (i == 3 && j == 10)
      {
        ABMD_DOUBLE k_dx = _dx_me * _dist3_me;
        ABMD_DOUBLE k_dy = _dy_me * _dist3_me;
        ABMD_DOUBLE k_dz = _dz_me * _dist3_me;

        userdata->fx[i] += to_double(userdata->masses[j] * k_dx);
        userdata->fy[i] += to_double(userdata->masses[j] * k_dy);
        userdata->fz[i] += to_double(userdata->masses[j] * k_dz);
        userdata->fx[j] -= to_double(userdata->masses[i] * k_dx);
        userdata->fy[j] -= to_double(userdata->masses[i] * k_dy);
        userdata->fz[j] -= to_double(userdata->masses[i] * k_dz);

        // f[6 * i + 3] += userdata->masses[j] * k_dx;
        // f[6 * i + 4] += userdata->masses[j] * k_dy;
        // f[6 * i + 5] += userdata->masses[j] * k_dz;
        // f[6 * j + 3] -= userdata->masses[i] * k_dx;
        // f[6 * j + 4] -= userdata->masses[i] * k_dy;
        // f[6 * j + 5] -= userdata->masses[i] * k_dz;
      }
      else
      {
        debug_type _dx = userdata->dx[i * barrier + j],
                   _dy = userdata->dy[i * barrier + j],
                   _dz = userdata->dz[i * barrier + j];
        debug_type _dist3 = userdata->dist3[i * barrier + j];

        debug_type k_dx = _dx * _dist3;
        debug_type k_dy = _dy * _dist3;
        debug_type k_dz = _dz * _dist3;
        
        userdata->fx[i] += userdata->masses[j] * k_dx;
        userdata->fy[i] += userdata->masses[j] * k_dy;
        userdata->fz[i] += userdata->masses[j] * k_dz;
        userdata->fx[j] -= userdata->masses[i] * k_dx;
        userdata->fy[j] -= userdata->masses[i] * k_dy;
        userdata->fz[j] -= userdata->masses[i] * k_dz;

        // f[6 * i + 3] += ABMD_DOUBLE(userdata->masses[j] * k_dx);
        // f[6 * i + 4] += ABMD_DOUBLE(userdata->masses[j] * k_dy);
        // f[6 * i + 5] += ABMD_DOUBLE(userdata->masses[j] * k_dz);
        // f[6 * j + 3] -= ABMD_DOUBLE(userdata->masses[i] * k_dx);
        // f[6 * j + 4] -= ABMD_DOUBLE(userdata->masses[i] * k_dy);
        // f[6 * j + 5] -= ABMD_DOUBLE(userdata->masses[i] * k_dz);
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

  for (i = 0; i < userdata->n_objects; i++)
  {

    f[6 * i] = to_double(f[6 * i]);
    f[6 * i + 1] = to_double(f[6 * i + 1]);
    f[6 * i + 2] = to_double(f[6 * i + 2]);

    f[6 * i + 3] = to_double(f[6 * i + 3]);
    f[6 * i + 4] = to_double(f[6 * i + 4]);
    f[6 * i + 5] = to_double(f[6 * i + 5]);
  }
}