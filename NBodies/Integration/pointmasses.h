#include <cmath>
#include <iostream>
#include <vector>
const int earthNum = 3;
const int moonNum = 10;
const int barrier = 16;

template <typename Type>
class ObjectsData
{
public:
  ObjectsData(std::vector<Type> &data)
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
    this->temp_acc.resize(this->n_objects * 3);
  }

  std::vector<Type> masses;
  int n_objects;
  std::vector<Type> dx, dy, dz;
  std::vector<Type> dist, dist2, dist3;
  std::vector<Type> temp_acc;
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

  // Precalculate dx, dy, dz, distances
  for (i = 0; i < barrier; i++)
  {
    for (j = i + 1; j < barrier; j++)
    {
      ABMD_DOUBLE _dx, _dy, _dz, _dist2, _dist, _dist3;
      if (i == 3 && j == 10)
      {
        _dx = x[6 * j] - x[6 * i];
        _dy = x[6 * j + 1] - x[6 * i + 1];
        _dz = x[6 * j + 2] - x[6 * i + 2];
        _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
      }
      else
      {
        _dx = ABMD_DOUBLE(to_double(x[6 * j]) - to_double(x[6 * i]));
        _dy = ABMD_DOUBLE(to_double(x[6 * j + 1]) - to_double(x[6 * i + 1]));
        _dz = ABMD_DOUBLE(to_double(x[6 * j + 2]) - to_double(x[6 * i + 2]));
        _dist2 = 1.0 / ABMD_DOUBLE(to_double(_dx * _dx) + to_double(_dy * _dy) + to_double(_dz * _dz));
      }

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

  // Go backwards so that the acceleration from the Sun is added last
  for (i = barrier - 1; i >= 0; i--)
  {
    for (j = i + 1; j < barrier; j++)
    {
      ABMD_DOUBLE _dx = userdata->dx[i * barrier + j],
                  _dy = userdata->dy[i * barrier + j],
                  _dz = userdata->dz[i * barrier + j];
      ABMD_DOUBLE _dist3 = userdata->dist3[i * barrier + j];

      ABMD_DOUBLE k_dx = _dx * _dist3;
      ABMD_DOUBLE k_dy = _dy * _dist3;
      ABMD_DOUBLE k_dz = _dz * _dist3;

      f[6 * i + 3] += ABMD_DOUBLE(to_double(userdata->masses[j]) * to_double(k_dx));
      f[6 * i + 4] += ABMD_DOUBLE(to_double(userdata->masses[j]) * to_double(k_dy));
      f[6 * i + 5] += ABMD_DOUBLE(to_double(userdata->masses[j]) * to_double(k_dz));
      f[6 * j + 3] -= ABMD_DOUBLE(to_double(userdata->masses[i]) * to_double(k_dx));
      f[6 * j + 4] -= ABMD_DOUBLE(to_double(userdata->masses[i]) * to_double(k_dy));
      f[6 * j + 5] -= ABMD_DOUBLE(to_double(userdata->masses[i]) * to_double(k_dz));
    }
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