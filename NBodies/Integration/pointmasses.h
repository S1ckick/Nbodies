#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>

#include "../def.h"

#define dotProduct(x1, y1, z1, x2, y2, z2) ((x1) * (x2) + (y1) * (y2) + (z1) * (z2))
#define vecLen2(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))

long double to_double(const long double &x)
{
  return x;
}


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
  std::vector<helper_type> temp_acc;
  std::vector<helper_type> fx, fy, fz;
};
#define NO_RELATIVITY

template <typename ABMD_DOUBLE>
struct ContextData
{
  ObjectsData<ABMD_DOUBLE> *objects;
  double *callback_t;
  double *center;
  int i;
  int dim;
  std::ofstream f;
  std::ofstream fb;
  std::ofstream fi;
  std::ofstream fen;
};

helper_type tay_arr[] = {-3.0 / 2.0, 3.0 / 8.0, 1. / 16., 3. / 128., 3. / 256., 7. / 1024., 9. / 2048., 99. / 32768., 143./65536., 429./262144.};

helper_type tay(helper_type x, int n)
{
  helper_type res = 0;
  helper_type temp_x = x;
  for (int i = 0; i < n; i++)
  {
    res += tay_arr[i] * x;
    x *= temp_x;
  }
  return res;
}

template <typename ABMD_DOUBLE>
void pointmassesCalculateXdot_tmp(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *f, void *context)
{

  int i, j, k;
  // userdata()
  //  Copy velocities of bodies
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

  // Земля в барицентре
  ABMD_DOUBLE earth_x = x[6 * earthNum], earth_y = x[6 * earthNum + 1],
              earth_z = x[6 * earthNum + 2];
  ABMD_DOUBLE earth_vx = x[6 * earthNum + 3], earth_vy = x[6 * earthNum + 4],
              earth_vz = x[6 * earthNum + 5];

  // Луна в геоцентре, числа на 3 порядка меньше, чем координаты Земли в барицентре
  ABMD_DOUBLE gcmoon_x = x[6 * moonNum], gcmoon_y = x[6 * moonNum + 1],
              gcmoon_z = x[6 * moonNum + 2];
  ABMD_DOUBLE gcmoon_vx = x[6 * moonNum + 3], gcmoon_vy = x[6 * moonNum + 4],
              gcmoon_vz = x[6 * moonNum + 5];

  // Переводим координаты Луны в барицентрическую систему
  // Теперь здесь числа такого же порядка как и у Земли, Луна близка к Земле
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
      if (i == earthNum && j == moonNum)
      {
        // Так мы берем расстояние между Землей и Луной (Луна в геоцентре должна быть с предыдущих шагов)
        // !!! маленькие числа
        _dx_me = gcmoon_x;
        _dy_me = gcmoon_y;
        _dz_me = gcmoon_z;
        _dist2_me = 1.0 / (_dx_me * _dx_me + _dy_me * _dy_me + _dz_me * _dz_me);
        _dist_me = sqrt(_dist2_me);

        // сохраняется 1/||r_earth - r_moon||^3
        _dist3_me = _dist_me * _dist2_me;

        // сохраняем (r_moon - r_earth)
        userdata->dx[i * barrier + j] = to_double(_dx_me);
        userdata->dy[i * barrier + j] = to_double(_dy_me);
        userdata->dz[i * barrier + j] = to_double(_dz_me);

        // сохраняем (r_earth - r_moon)
        userdata->dx[j * barrier + i] = -to_double(_dx_me);
        userdata->dy[j * barrier + i] = -to_double(_dy_me);
        userdata->dz[j * barrier + i] = -to_double(_dz_me);

        userdata->dist[i * barrier + j] = to_double(_dist_me);
        userdata->dist2[i * barrier + j] = to_double(_dist2_me);
        userdata->dist3[i * barrier + j] = to_double(_dist3_me);

        userdata->dist[j * barrier + i] = to_double(_dist_me);
        userdata->dist2[j * barrier + i] = to_double(_dist2_me);
        userdata->dist3[j * barrier + i] = to_double(_dist3_me);
      }
      else
      {
        helper_type _dx, _dy, _dz, _dist2, _dist, _dist3;
        _dx = to_double(x[6 * j] - x[6 * i]);
        _dy = to_double(x[6 * j + 1] - x[6 * i + 1]);
        _dz = to_double(x[6 * j + 2] - x[6 * i + 2]);
        _dist2 = 1.0 / ((_dx * _dx) + (_dy * _dy) + (_dz * _dz));
        _dist = sqrt(_dist2);
        // честно посчитали 1/||r1-r2||^3
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
#ifdef TAYLOR
    if(i == moonNum)
      continue;
#endif
    for (j = i + 1; j < barrier; j++)
    {
#ifdef TAYLOR
      if(j == moonNum)
        continue;
#endif
      // Если не пара Луна-Земля
      if (i != earthNum || j != moonNum)
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

    // На этом этапе к Луне и к Земле взаимодействие Луна-Земля еще не добавлены

    // Астероиды

    for (j = barrier; j < userdata->n_objects; j++)
    {
      // Five individual asteroids must not perturb the discrete asteroid ring
      if ((j >= 390) && (j < 570) && (i > moonNum))
        continue;

      helper_type k_dx, k_dy, k_dz;
#ifdef TAYLOR
      if (i != moonNum)
      {
#endif
        helper_type _dx = to_double(x[6 * j]) - to_double(x[6 * i]),
                    _dy = to_double(x[6 * j + 1]) - to_double(x[6 * i + 1]),
                    _dz = to_double(x[6 * j + 2]) - to_double(x[6 * i + 2]);
        helper_type _dist2 = 1.0 / (_dx * _dx + _dy * _dy + _dz * _dz);
        helper_type _dist = sqrt(_dist2);
        helper_type _dist3 = _dist * _dist2;

        k_dx = _dx * _dist3;
        k_dy = _dy * _dist3;
        k_dz = _dz * _dist3;

        userdata->fx[i] += userdata->masses[j] * k_dx;
        userdata->fy[i] += userdata->masses[j] * k_dy;
        userdata->fz[i] += userdata->masses[j] * k_dz;
        userdata->fx[j] -= userdata->masses[i] * k_dx;
        userdata->fy[j] -= userdata->masses[i] * k_dy;
        userdata->fz[j] -= userdata->masses[i] * k_dz;
#ifdef TAYLOR
      }
#endif
    }
  }

  for (int i = 0; i < userdata->n_objects; i++)
  {
    f[6 * i + 3] = userdata->fx[i];
    f[6 * i + 4] = userdata->fy[i];
    f[6 * i + 5] = userdata->fz[i];
  }

#ifdef TAYLOR
  // i == moonNum
  for (j = userdata->n_objects - 1; j >= 0; j--)
  {
    if(j == moonNum || j == earthNum)
      continue;
    // _dx, _dy, _dz -- r_am
    // _dist3 -- ||r_am||^3

    // r_em:  gcmoon_x
    // r_ae:  x[6 * j] - earth_x

    helper_type rem_x = to_double(gcmoon_x);
    helper_type rem_y = to_double(gcmoon_y);
    helper_type rem_z = to_double(gcmoon_z);

    helper_type rae_x = (to_double(x[6 * j]) - to_double(earth_x));
    helper_type rae_y = (to_double(x[6 * j + 1]) - to_double(earth_y));
    helper_type rae_z = (to_double(x[6 * j + 2]) - to_double(earth_z));

    helper_type bmoon_x = to_double(gcmoon_x) + to_double(earth_x);
    helper_type bmoon_y = to_double(gcmoon_y) + to_double(earth_y);
    helper_type bmoon_z = to_double(gcmoon_z) + to_double(earth_z);

    helper_type ram_x = (to_double(x[6 * j]) - to_double(bmoon_x));
    helper_type ram_y = (to_double(x[6 * j + 1]) - to_double(bmoon_y));
    helper_type ram_z = (to_double(x[6 * j + 2]) - to_double(bmoon_z));

    helper_type dot_emae = dotProduct(rem_x, rem_y, rem_z, rae_x, rae_y, rae_z);
    helper_type dot_em = vecLen2(rem_x, rem_y, rem_z);
    helper_type dot_ae = vecLen2(rae_x, rae_y, rae_z);

    helper_type r_x = (2 * dot_emae - dot_em) / dot_ae;

    helper_type tay_res = tay(r_x, 4);
    

    // В строчках ниже rae_x - должно быть большое число, в то время как rem_x малое
    // Но мы от малого rem_x отнимаем rae_x * tay_res, где tay_res малое
    // Так делаем ошибку округления меньше
    helper_type _dx = -rem_x - rae_x * tay_res;
    helper_type _dy = -rem_y - rae_y * tay_res;
    helper_type _dz = -rem_z - rae_z * tay_res;

    helper_type _dist2 = 1.0 / vecLen2(ram_x, ram_y, ram_z);
    helper_type _dist = sqrt(_dist2);
    helper_type _dist3 = _dist * _dist2;

    helper_type dist_ae2 = 1.0 / vecLen2(rae_x, rae_y, rae_z);
    helper_type dist_ae = sqrt(dist_ae2);
    helper_type dist_ae3 = dist_ae * dist_ae2;

    helper_type to_moon_x = userdata->masses[j] * _dx * _dist3;
    helper_type to_moon_y = userdata->masses[j] * _dy * _dist3;
    helper_type to_moon_z = userdata->masses[j] * _dz * _dist3;

    f[moonNum * 6 + 3] += to_moon_x;
    f[moonNum * 6 + 4] += to_moon_y;
    f[moonNum * 6 + 5] += to_moon_z;
    f[j * 6 + 3] -= userdata->masses[moonNum] * ram_x * _dist3;
    f[j * 6 + 4] -= userdata->masses[moonNum] * ram_y * _dist3;
    f[j * 6 + 5] -= userdata->masses[moonNum] * ram_z * _dist3;

  }
#endif

  // Добавили взаимодействие Луна-Астероиды

  // В самом конце добавляем взаимодействие Луна-Земля
  ABMD_DOUBLE k_dx = _dx_me * _dist3_me;
  ABMD_DOUBLE k_dy = _dy_me * _dist3_me;
  ABMD_DOUBLE k_dz = _dz_me * _dist3_me;

#ifndef TAYLOR
  f[6 * moonNum + 3] -= f[6 * earthNum + 3];
  f[6 * moonNum + 4] -= f[6 * earthNum + 4];
  f[6 * moonNum + 5] -= f[6 * earthNum + 5];
#endif

  f[6 * earthNum + 3] += to_double(userdata->masses[moonNum] * k_dx);
  f[6 * earthNum + 4] += to_double(userdata->masses[moonNum] * k_dy);
  f[6 * earthNum + 5] += to_double(userdata->masses[moonNum] * k_dz);
  // Сумма масс, потому что в f ускорение Луны в геоцентр мы перевели
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
#ifdef RELATIVISTIC
  double revLS2 = 3.335661199676477670e-005;
  // // ------------------------------------RELATIVISTIC------------------------------------
  for (i = 0; i < barrier; ++i)
  {
    userdata->temp_acc[3 * i] = to_double(f[6 * i + 3]);
    userdata->temp_acc[3 * i + 1] = to_double(f[6 * i + 4]);
    userdata->temp_acc[3 * i + 2] = to_double(f[6 * i + 5]);
    // в f ускорение Луны в геоцентре
    // надо перевести ее в барицентр
    if (i == moonNum)
    {
      userdata->temp_acc[3 * i] += to_double(f[6 * earthNum + 3]);
      userdata->temp_acc[3 * i + 1] += to_double(f[6 * earthNum + 4]);
      userdata->temp_acc[3 * i + 2] += to_double(f[6 * earthNum + 5]);
    }
  }

  helper_type acc_e_x = 0.0;
  helper_type acc_e_y = 0.0;
  helper_type acc_e_z = 0.0;

  for (i = 0; i < barrier; ++i)
  {
    // (v_A)^2
    helper_type vi_sqr = to_double(x[6 * i + 3]) * to_double(x[6 * i + 3]) +
                    to_double(x[6 * i + 4]) * to_double(x[6 * i + 4]) +
                    to_double(x[6 * i + 5]) * to_double(x[6 * i + 5]);
    // (sum_{C != A} Gm_C/r_AC)
    helper_type sum_i = 0.0;
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

      // (v_B)^2
      helper_type vj_sqr = to_double(x[6 * j + 3]) * to_double(x[6 * j + 3]) +
                      to_double(x[6 * j + 4]) * to_double(x[6 * j + 4]) +
                      to_double(x[6 * j + 5]) * to_double(x[6 * j + 5]);
      // v_A dot v_B
      helper_type vi_dot_vj = to_double(x[6 * i + 3]) * to_double(x[6 * j + 3]) +
                         to_double(x[6 * i + 4]) * to_double(x[6 * j + 4]) +
                         to_double(x[6 * i + 5]) * to_double(x[6 * j + 5]);

      // (sum_{C != B} Gm_C/r_BC)
      helper_type sum_j = 0.0;
      for (k = 0; k < barrier; ++k)
      {
        if (j == k)
          continue;
        sum_j += userdata->masses[k] * userdata->dist[barrier * j + k];
      }

      helper_type t1 = -4.0 * sum_i - sum_j + vi_sqr + 2.0 * vj_sqr - 4.0 * vi_dot_vj;
      helper_type t2 = userdata->dist[ij] * (-userdata->dx[ij] * to_double(x[6 * j + 3]) +
                                        -userdata->dy[ij] * to_double(x[6 * j + 4]) +
                                        -userdata->dz[ij] * to_double(x[6 * j + 5]));
      t2 = -1.5 * t2 * t2;
      helper_type t3 = 0.5 * to_double(userdata->dx[ij] * userdata->temp_acc[3 * j] +
                                  userdata->dy[ij] * userdata->temp_acc[3 * j + 1] +
                                  userdata->dz[ij] * userdata->temp_acc[3 * j + 2]);

      helper_type c1 = revLS2 * userdata->masses[j] * userdata->dist3[ij] * (t1 + t2 + t3);
      // ok
      f[6 * i + 3] += userdata->dx[ij] * c1;
      f[6 * i + 4] += userdata->dy[ij] * c1;
      f[6 * i + 5] += userdata->dz[ij] * c1;

      helper_type p1 = userdata->dist[ij] * (-userdata->dx[ij] * (4. * to_double(x[6 * i + 3]) - 3. * to_double(x[6 * j + 3])) +
                                        -userdata->dy[ij] * (4. * to_double(x[6 * i + 4]) - 3. * to_double(x[6 * j + 4])) +
                                        -userdata->dz[ij] * (4. * to_double(x[6 * i + 5]) - 3. * to_double(x[6 * j + 5])));

      helper_type c2 = revLS2 * userdata->masses[j] * userdata->dist2[ij] * p1;

      helper_type dv_x = to_double(x[6 * i + 3]) - to_double(x[6 * j + 3]);
      helper_type dv_y = to_double(x[6 * i + 4]) - to_double(x[6 * j + 4]);
      helper_type dv_z = to_double(x[6 * i + 5]) - to_double(x[6 * j + 5]);
      //ok
      f[6 * i + 3] += dv_x * c2;
      f[6 * i + 4] += dv_y * c2;
      f[6 * i + 5] += dv_z * c2;

      helper_type c3 = 3.5 * revLS2 * userdata->masses[j] * userdata->dist[ij];
      //ok
      f[6 * i + 3] += userdata->temp_acc[3 * j] * c3;
      f[6 * i + 4] += userdata->temp_acc[3 * j + 1] * c3;
      f[6 * i + 5] += userdata->temp_acc[3 * j + 2] * c3;

      if(i == earthNum){
        acc_e_x += userdata->dx[ij] * c1 + dv_x * c2 + userdata->temp_acc[3 * j] * c3;
        acc_e_y += userdata->dy[ij] * c1 + dv_y * c2 + userdata->temp_acc[3 * j + 1] * c3;
        acc_e_z += userdata->dz[ij] * c1 + dv_z * c2 + userdata->temp_acc[3 * j + 2] * c3;
      }
    }
  }

  f[6 * moonNum + 3] -= acc_e_x;
  f[6 * moonNum + 4] -= acc_e_y;
  f[6 * moonNum + 5] -= acc_e_z;
#endif

        // Переводим Луну в геоцентр
  // !!!Теперь ее координаты малы
  x[6 * moonNum] = gcmoon_x;
  x[6 * moonNum + 1] = gcmoon_y;
  x[6 * moonNum + 2] = gcmoon_z;
  x[6 * moonNum + 3] = gcmoon_vx;
  x[6 * moonNum + 4] = gcmoon_vy;
  x[6 * moonNum + 5] = gcmoon_vz;


}