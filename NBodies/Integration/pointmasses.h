#include <vector>
const int earthNum = 3;
const int moonNum = 10;
const int barrier = 16;

template <typename Type>
class ObjectsData {
 public:
  ObjectsData(std::vector<Type> &data) {
    this->n_objects = data.size();
    masses.resize(this->n_objects);
    for (int i = 0; i < data.size(); i++) {
      this->masses[i] = data[i];
    }
    this->dx.resize(this->n_objects);
    this->dy.resize(this->n_objects);
    this->dz.resize(this->n_objects);
    this->dist.resize(this->n_objects);
    this->dist2.resize(this->n_objects);
    this->dist3.resize(this->n_objects);
    this->temp_acc.resize(this->n_objects * 3);
  }

  std::vector<Type> masses;
  int n_objects;
  std::vector<Type> dx, dy, dz;
  std::vector<Type> dist, dist2, dist3;
  std::vector<Type> temp_acc;
};

template <typename Type>
void pointmassesCalculateXdot(std::vector<Type> &x, std::vector<Type> &f,
                              ObjectsData<Type> &userdata) {
  int i, j, k;

  // Copy velocities of bodies
  for (i = 0; i < userdata.n_objects; i++) {
    f[6 * i] = x[6 * i + 3];
    f[6 * i + 1] = x[6 * i + 4];
    f[6 * i + 2] = x[6 * i + 5];
  }

  // Convert positions and velocities from EMB to barycentric in case of Earth
  // and Moon
  Type earth_x = x[6 * earthNum], earth_y = x[6 * earthNum + 1],
       earth_z = x[6 * earthNum + 2];
  Type earth_vx = x[6 * earthNum + 3], earth_vy = x[6 * earthNum + 4],
       earth_vz = x[6 * earthNum + 5];
  Type gcmoon_x = x[6 * moonNum], gcmoon_y = x[6 * moonNum + 1],
       gcmoon_z = x[6 * moonNum + 2];
  Type gcmoon_vx = x[6 * moonNum + 3], gcmoon_vy = x[6 * moonNum + 4],
       gcmoon_vz = x[6 * moonNum + 5];

  x[6 * moonNum] = earth_x + gcmoon_x;
  x[6 * moonNum + 1] = earth_y + gcmoon_y;
  x[6 * moonNum + 2] = earth_z + gcmoon_z;
  x[6 * moonNum + 3] = earth_vx + gcmoon_vx;
  x[6 * moonNum + 4] = earth_vy + gcmoon_vy;
  x[6 * moonNum + 5] = earth_vz + gcmoon_vz;

  // Precalculate dx, dy, dz, distances
  for (i = 0; i < barrier; i++) {
    for (j = i + 1; j < barrier; j++) {
      Type _dx = x[6 * j] - x[6 * i], _dy = x[6 * j + 1] - x[6 * i + 1],
           _dz = x[6 * j + 2] - x[6 * i + 2];
      Type _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
      Type _dist = std::sqrt(_dist2);
      Type _dist3 = _dist * _dist2;

      userdata.dx[i * barrier + j] = _dx;
      userdata.dy[i * barrier + j] = _dy;
      userdata.dz[i * barrier + j] = _dz;

      userdata.dx[j * barrier + i] = -_dx;
      userdata.dy[j * barrier + i] = -_dy;
      userdata.dz[j * barrier + i] = -_dz;

      userdata.dist[i * barrier + j] = _dist;
      userdata.dist2[i * barrier + j] = _dist2;
      userdata.dist3[i * barrier + j] = _dist3;

      userdata.dist[j * barrier + i] = _dist;
      userdata.dist2[j * barrier + i] = _dist2;
      userdata.dist3[j * barrier + i] = _dist3;
    }
  }

  // Go backwards so that the acceleration from the Sun is added last
  for (i = barrier - 1; i >= 0; i--) {
    for (j = i + 1; j < barrier; j++) {
      Type _dx = userdata.dx[i * barrier + j],
           _dy = userdata.dy[i * barrier + j],
           _dz = userdata.dz[i * barrier + j];
      Type _dist3 = userdata.dist3[i * barrier + j];

      Type k_dx = _dx * _dist3;
      Type k_dy = _dy * _dist3;
      Type k_dz = _dz * _dist3;

      f[6 * i + 3] += userdata.masses[j] * k_dx;
      f[6 * i + 4] += userdata.masses[j] * k_dy;
      f[6 * i + 5] += userdata.masses[j] * k_dz;
      f[6 * j + 3] -= userdata.masses[i] * k_dx;
      f[6 * j + 4] -= userdata.masses[i] * k_dy;
      f[6 * j + 5] -= userdata.masses[i] * k_dz;
    }

    // simplified continuation of the above loop
    for (j = barrier; j < userdata.n_objects; j++) {
      Type _dx = x[6 * j] - x[6 * i], _dy = x[6 * j + 1] - x[6 * i + 1],
           _dz = x[6 * j + 2] - x[6 * i + 2];
      Type _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
      Type _dist = std::sqrt(_dist2);

      Type _dist3 = _dist * _dist2;

      Type k_dx = _dx * _dist3;
      Type k_dy = _dy * _dist3;
      Type k_dz = _dz * _dist3;

      f[6 * i + 3] += userdata.masses[j] * k_dx;
      f[6 * i + 4] += userdata.masses[j] * k_dy;
      f[6 * i + 5] += userdata.masses[j] * k_dz;

      f[6 * j + 3] -= userdata.masses[i] * k_dx;
      f[6 * j + 4] -= userdata.masses[i] * k_dy;
      f[6 * j + 5] -= userdata.masses[i] * k_dz;
    }
  }

  // Relativistic forces computation
#ifndef NO_RELATIVITY
  for (i = 0; i < barrier; ++i) {
    userdata.temp_acc[3 * i] = f[6 * i + 3];
    userdata.temp_acc[3 * i + 1] = f[6 * i + 4];
    userdata.temp_acc[3 * i + 2] = f[6 * i + 5];
  }

  for (i = 0; i < barrier; ++i) {
    Type vi_sqr = x[6 * i + 3] * x[6 * i + 3] + x[6 * i + 4] * x[6 * i + 4] +
                  x[6 * i + 5] * x[6 * i + 5];

    Type sum_i = 0.0;
    for (k = 0; k < barrier; ++k) {
      if (i == k) continue;
      sum_i += userdata.masses[k] * userdata.dist[barrier * i + k];
    }

    for (j = 0; j < barrier; ++j) {
      if (i == j) continue;

      int ij = barrier * i + j;

      Type vj_sqr = x[6 * j + 3] * x[6 * j + 3] + x[6 * j + 4] * x[6 * j + 4] +
                    x[6 * j + 5] * x[6 * j + 5];

      Type vi_dot_vj = x[6 * i + 3] * x[6 * j + 3] +
                       x[6 * i + 4] * x[6 * j + 4] +
                       x[6 * i + 5] * x[6 * j + 5];

      Type sum_j = 0.0;
      for (k = 0; k < barrier; ++k) {
        if (j == k) continue;
        sum_j += userdata.masses[k] * userdata.dist[barrier * j + k];
      }

      Type t1 = -4.0 * sum_i - sum_j + vi_sqr + 2.0 * vj_sqr - 4.0 * vi_dot_vj;
      Type t2 = userdata.dist[ij] * (-userdata.dx[ij] * x[6 * j + 3] +
                                     -userdata.dy[ij] * x[6 * j + 4] +
                                     -userdata.dz[ij] * x[6 * j + 5]);
      t2 = -1.5 * t2 * t2;
      Type t3 = 0.5 * (userdata.dx[ij] * userdata.temp_acc[3 * j] +
                       userdata.dy[ij] * userdata.temp_acc[3 * j + 1] +
                       userdata.dz[ij] * userdata.temp_acc[3 * j + 2]);

      Type c1 = userdata.masses[j] * userdata.dist3[ij] * (t1 + t2 + t3);
      f[6 * i + 3] += userdata.dx[ij] * c1;
      f[6 * i + 4] += userdata.dy[ij] * c1;
      f[6 * i + 5] += userdata.dz[ij] * c1;

      Type p1 = userdata.dist[ij] *
                (-userdata.dx[ij] * (4 * x[6 * i + 3] - 3 * x[6 * j + 3]) +
                 -userdata.dy[ij] * (4 * x[6 * i + 4] - 3 * x[6 * j + 4]) +
                 -userdata.dz[ij] * (4 * x[6 * i + 5] - 3 * x[6 * j + 5]));

      Type c2 = userdata.masses[j] * userdata.dist2[ij] * p1;

      Type dv_x = x[6 * i + 3] - x[6 * j + 3];
      Type dv_y = x[6 * i + 4] - x[6 * j + 4];
      Type dv_z = x[6 * i + 5] - x[6 * j + 5];
      f[6 * i + 3] += dv_x * c2;
      f[6 * i + 4] += dv_y * c2;
      f[6 * i + 5] += dv_z * c2;

      Type c3 = 3.5 * userdata.masses[j] * userdata.dist[ij];
      f[6 * i + 3] += userdata.temp_acc[3 * j] * c3;
      f[6 * i + 4] += userdata.temp_acc[3 * j + 1] * c3;
      f[6 * i + 5] += userdata.temp_acc[3 * j + 2] * c3;
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
}