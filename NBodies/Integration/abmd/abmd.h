#ifndef ABMD_H
#define ABMD_H

#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>

std::vector<std::string> split(std::string s, std::string delimiter) {
  std::size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }

  res.push_back(s.substr(pos_start));
  return res;
}

#include "../pointmasses.h"
#include "abmd_rk.h"

template <typename ABMD_DOUBLE>
class coefss {
 public:
  static ABMD_DOUBLE PREDICTOR_COEFFS[ABMD_MAX_ORDER];
};
#ifdef NUMBER_DOUBLE_DOUBLE
template <>
dd_real coefss<dd_real>::PREDICTOR_COEFFS[ABMD_MAX_ORDER] = {
    dd_real(1.0),
    dd_real(1.0) / dd_real(2.0),
    dd_real(5.0) / dd_real(12.0),
    dd_real(3.0) / dd_real(8.0),
    dd_real(251.0) / dd_real(720.0),
    dd_real(95.0) / dd_real(288.0),
    dd_real(19087.0) / dd_real(60480.0),
    dd_real(5257.0) / dd_real(17280.0),
    dd_real(1070017.0) / dd_real(3628800.0),
    dd_real(25713.0) / dd_real(89600.0),
    dd_real(26842253.0) / dd_real(95800320.0),
    dd_real(4777223.0) / dd_real(17418240.0),
    dd_real(703604254357.0) / dd_real(2615348736000.0),
    dd_real(106364763817.0) / dd_real(402361344000.0),
    dd_real(1166309819657.0) / dd_real(4483454976000.0),
    dd_real(2.5221445e7) / dd_real(9.8402304e7),
    dd_real(8.092989203533249e15) / dd_real(3.201186852864e16),
    dd_real(8.5455477715379e13) / dd_real(3.4237292544e14),
    dd_real(1.2600467236042756559e19) / dd_real(5.109094217170944e19)};
#elif NUMBER_DOUBLE == 1
template <>
long double coefss<long double>::PREDICTOR_COEFFS[ABMD_MAX_ORDER] = {
    (long double)(1.0),
    (long double)(1.0) / (long double)(2.0),
    (long double)(5.0) / (long double)(12.0),
    (long double)(3.0) / (long double)(8.0),
    (long double)(251.0) / (long double)(720.0),
    (long double)(95.0) / (long double)(288.0),
    (long double)(19087.0) / (long double)(60480.0),
    (long double)(5257.0) / (long double)(17280.0),
    (long double)(1070017.0) / (long double)(3628800.0),
    (long double)(25713.0) / (long double)(89600.0),
    (long double)(26842253.0) / (long double)(95800320.0),
    (long double)(4777223.0) / (long double)(17418240.0),
    (long double)(703604254357.0) / (long double)(2615348736000.0),
    (long double)(106364763817.0) / (long double)(402361344000.0),
    (long double)(1166309819657.0) / (long double)(4483454976000.0),
    (long double)(2.5221445e7) / (long double)(9.8402304e7),
    (long double)(8.092989203533249e15) / (long double)(3.201186852864e16),
    (long double)(8.5455477715379e13) / (long double)(3.4237292544e14),
    (long double)(1.2600467236042756559e19) /
        (long double)(5.109094217170944e19)};
#else
template <>
double coefss<double>::PREDICTOR_COEFFS[ABMD_MAX_ORDER] = {
    (double)(1.0),
    (double)(1.0) / (double)(2.0),
    (double)(5.0) / (double)(12.0),
    (double)(3.0) / (double)(8.0),
    (double)(251.0) / (double)(720.0),
    (double)(95.0) / (double)(288.0),
    (double)(19087.0) / (double)(60480.0),
    (double)(5257.0) / (double)(17280.0),
    (double)(1070017.0) / (double)(3628800.0),
    (double)(25713.0) / (double)(89600.0),
    (double)(26842253.0) / (double)(95800320.0),
    (double)(4777223.0) / (double)(17418240.0),
    (double)(703604254357.0) / (double)(2615348736000.0),
    (double)(106364763817.0) / (double)(402361344000.0),
    (double)(1166309819657.0) / (double)(4483454976000.0),
    (double)(2.5221445e7) / (double)(9.8402304e7),
    (double)(8.092989203533249e15) / (double)(3.201186852864e16),
    (double)(8.5455477715379e13) / (double)(3.4237292544e14),
    (double)(1.2600467236042756559e19) / (double)(5.109094217170944e19)};
#endif
/********** ABMD method **********/
template <typename ABMD_DOUBLE>
struct ABMData {
  ABMD<ABMD_DOUBLE> input;
  double rk4_h;
  ABMD_DOUBLE *temp;
  Queue<ABMD_DOUBLE> *queue;
  ABMD_DOUBLE *xs_delayed;
  ABMD_DOUBLE *xs_delayed_inner;
  ABMD_DOUBLE *xs_delayed_tmp;
  ABMD_DOUBLE *dxs_delayed;
  ABMD_DOUBLE *rk_memory;
  ABMD_DOUBLE *inner_rk_memory;
  ABMD_DOUBLE *hoho;

  ABMData(){};
  ~ABMData() {
    free(temp);
    free(xs_delayed);
    free(xs_delayed_inner);
    free(xs_delayed_tmp);
    free(dxs_delayed);
    free(rk_memory);
    free(inner_rk_memory);
  }

  void predict() {
    int dim = this->input.dim;
    int abm_order = this->input.abm_order;
    double h = this->input.h;
    Queue<ABMD_DOUBLE> *queue = this->queue;

    queue->pop();
    ABMD_DOUBLE *prev = queue->peek_right_x();
    ABMD_DOUBLE *out = queue->push();

    memset(out, 0, sizeof(ABMD_DOUBLE) * dim);

    helper_type *diffs = queue->get_diffs_r();
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < abm_order; j++) {
        out[i] += diffs[i * abm_order + j] * this->hoho[j];
      }
      out[i] += prev[i];
    }
  }

  void correct(ABMD_DOUBLE *x_predicted) {
    int dim = this->input.dim;
    int abm_order = this->input.abm_order;
    double h = this->input.h;
    Queue<ABMD_DOUBLE> *queue = this->queue;
    ABMD_DOUBLE *out = queue->peek_right_x();

    if (x_predicted == NULL) x_predicted = out;

    ABMD_DOUBLE *diff = queue->get_last_diff();
    ABMD_DOUBLE ch = h * coefss<ABMD_DOUBLE>::PREDICTOR_COEFFS[abm_order];
    for (int k = 0; k < dim; k++) {
      out[k] = x_predicted[k] + ch * diff[k];
    }
  }
};

template <typename ABMD_DOUBLE>
void copy_delayed_states(ABMD_DOUBLE x[], ABMD_DOUBLE xs_delayed[], int ndelays,
                         int **idxs, int *idxs_lens) {
  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < idxs_lens[i]; j++) {
      int idx = (idxs[i] == NULL) ? j : idxs[i][j];
      *xs_delayed++ = x[idx];
    }
  }
}

template <typename ABMD_DOUBLE>
void rhs(ABMD_DOUBLE x[], ABMD_DOUBLE xs_delayed[], ABMD_DOUBLE dxs_delayed[],
         double t, ABMD_DOUBLE *out, void *abm_data) {
  auto st1 = std::chrono::high_resolution_clock::now();
  ABMData<ABMD_DOUBLE> *data = (ABMData<ABMD_DOUBLE> *)abm_data;
  int dim = data->input.dim;
  ABMD_DOUBLE *temp = data->temp;

  memset(out, 0, sizeof(ABMD_DOUBLE) * dim);
  auto st2 = std::chrono::high_resolution_clock::now();

  if (data->input.f1 == NULL) {
    data->input.f2(x, xs_delayed, dxs_delayed, t, out, data->input.context);
    return;
  }

  data->input.f1(x, t, out, data->input.context);
  auto st3 = std::chrono::high_resolution_clock::now();

  if (data->input.f2 == NULL) return;

  ABMD_DOUBLE *out2 = &temp[0];
  memset(out2, 0, sizeof(ABMD_DOUBLE) * dim);
  data->input.f2(x, xs_delayed, dxs_delayed, t, out2, data->input.context);

  for (int i = 0; i < dim; i++) {
    out[i] += out2[i];
  }
}

template <typename ABMD_DOUBLE>
static void rhs_rk_inner(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *out,
                         void *abm_data) {
  ABMData<ABMD_DOUBLE> *data = (ABMData<ABMD_DOUBLE> *)abm_data;
  int ndelays = data->input.ndelays;
  int **idxs = data->input.delayed_idxs;
  int *idxs_lens = data->input.delayed_idxs_lens;
  ABMD_DOUBLE *xs_delayed = data->xs_delayed_inner;
  copy_delayed_states(x, xs_delayed, ndelays, idxs, idxs_lens);
  rhs<ABMD_DOUBLE>(x, xs_delayed, NULL, t, out, data);
}

template <typename ABMD_DOUBLE>
static void rhs_rk4(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *out,
                    void *abm_data) {
  ABMData<ABMD_DOUBLE> *data = (ABMData<ABMD_DOUBLE> *)abm_data;
  int dim = data->input.dim;
  int ndelays = data->input.ndelays;
  int **idxs = data->input.delayed_idxs;
  int *idxs_lens = data->input.delayed_idxs_lens;

  ABMD_DOUBLE *xs_delayed = data->xs_delayed;
  ABMD_DOUBLE *xs_delayed_tmp = data->xs_delayed_tmp;

  for (int i = 0; i < ndelays; i++) {
    double delay = data->input.delays[i];
    rk_step<ABMD_DOUBLE>(rhs_rk_inner, -delay, t, x, dim, data,
                         &xs_delayed_tmp[i * dim], NULL, &data->inner_rk_memory,
                         METHOD_RK4);
  }

  int k = 0;
  for (int i = 0; i < ndelays; i++) {
    for (int j = 0; j < idxs_lens[i]; j++) {
      int idx = (idxs[i] == NULL) ? j : idxs[i][j];
      xs_delayed[k] = xs_delayed_tmp[i * dim + idx];
      k++;
    }
  }

  rhs<ABMD_DOUBLE>(x, xs_delayed, NULL, t, out, data);
}

template <typename ABMD_DOUBLE>
void get_delayed_states(ABMData<ABMD_DOUBLE> *abm_data, double ti,
                        int last_dx_known, ABMD_DOUBLE *x_out,
                        ABMD_DOUBLE *dx_out) {
  Queue<ABMD_DOUBLE> *q = abm_data->queue;
  int ndelays = abm_data->input.ndelays;
  double *delays = abm_data->input.delays;
  int **idxs = abm_data->input.delayed_idxs;
  int *idxs_lens = abm_data->input.delayed_idxs_lens;
  int x_start_idx = 0;

  for (int i = 0; i < ndelays; i++) {
    double delay = delays[i];
    double t = ti - delay;

    q->evaluate_x_idxs(t, idxs[i], idxs_lens[i], &x_out[x_start_idx]);
    x_start_idx += idxs_lens[i];
  }

  if (dx_out == NULL) {
    return;
  }

  int dx_start_idx = 0;
  int *dx_delay_idxs = abm_data->input.dx_delays_idxs;

  for (int i = 0; i < abm_data->input.dx_delays_len; i++) {
    int delay_idx = dx_delay_idxs == NULL ? i : dx_delay_idxs[i];
    double t = ti - delays[delay_idx];
    q->evaluate_dx(t, idxs[i], idxs_lens[i], last_dx_known,
                   &dx_out[dx_start_idx]);
    dx_start_idx += idxs_lens[i];
  }
}

template <typename ABMD_DOUBLE>
void ABMD_run(ABMD<ABMD_DOUBLE> *abm) {
  int abm_order = abm->abm_order;

  // if (!(1 <= abm_order && abm_order <= 19)) {
  //   abm->error = "ABM order must be not less than 1 and not greater than 19";
  //   return 1;
  // }

  int ndelays = abm->ndelays;
  int delays_degree = abm->delays_poly_degree;
  int pointsave_degree = abm->pointsave_poly_degree;

  // if (ndelays > 0 && !(1 <= delays_degree && delays_degree <= abm_order)) {
  //   abm->error = "Delay interpolation degree must be not less than 1 and "
  //                "not greater than ABM order";
  //   return 1;
  // }

  // if (!(1 <= pointsave_degree && pointsave_degree <= abm_order)) {
  //   abm->error = "Pointsave interpolation degree must be not less than 1 and
  //   "
  //                "not greater than ABM order";
  //   return 1;
  // }

  // if (abm->f1 == NULL && abm->f2 == NULL) {
  //   abm->error = "Both RHSs are NULL";
  //   return 1;
  // }
  auto start = std::chrono::high_resolution_clock::now();
  int dim = abm->dim;
  double t0 = abm->t0;
  double t1 = abm->t1;

  int hsgn = (t1 > t0) - (t1 < t0);
  abm->h *= hsgn;
  double h = abm->h;

  ABMD_DOUBLE *init = (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * dim);
  for (int i = 0; i < dim; i++) {
    init[i] = abm->init[i];
  }

  if (abm->init_call != NULL) {
    abm->init_call(init, abm->context);
  }

  int RK_STEPS_IN_ABM = 8;

  int n = (int)(1 + (t1 - t0) / h);
  int rk4_i1 = abm_order - 1;
  rk4_i1 += 1;  // One more step to fill the whole queue with RK data
  double rk4_h = h / (double)RK_STEPS_IN_ABM;
  int rk4_n = rk4_i1 * RK_STEPS_IN_ABM;
  double rk4_t1 = t0 + rk4_i1 * h;

  int rk_size = rk4_n + 1;  // One more block to store the initial condition
  ABMD_DOUBLE *rk4_sol =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * rk_size * dim);
  ABMD_DOUBLE *rk4_rhss =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * rk_size * dim);
  ABMD_DOUBLE *rhs_temp = (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * 2 * dim);
  ABMD_DOUBLE *xs_delayed_tmp =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * ndelays * dim);

  int x_total_delays_len = 0;
  int dx_total_delays_len = 0;
  for (int i = 0; i < ndelays; i++) {
    x_total_delays_len += abm->delayed_idxs_lens[i];
  }
  for (int i = 0; i < abm->dx_delays_len; i++) {
    int *idxs = abm->dx_delays_idxs;
    int delay_idx = idxs == NULL ? i : idxs[i];
    dx_total_delays_len += abm->delayed_idxs_lens[delay_idx];
  }
  size_t x_delayed_size = x_total_delays_len * sizeof(ABMD_DOUBLE);

  ABMD_DOUBLE *xs_delayed = (ABMD_DOUBLE *)malloc(x_delayed_size);
  ABMD_DOUBLE *xs_delayed_inner = (ABMD_DOUBLE *)malloc(x_delayed_size);
  ABMD_DOUBLE *dxs_delayed =
      (ABMD_DOUBLE *)malloc(dx_total_delays_len * sizeof(ABMD_DOUBLE));

  int queue_size = abm_order + 1;
  Queue<ABMD_DOUBLE> *queue = new Queue<ABMD_DOUBLE>(queue_size, dim);
  //*queue = // = (queue_size, dim);

  queue->delays_poly_degree = delays_degree;
  queue->pointsave_poly_degree = pointsave_degree;
  queue->t0 = t0;
  queue->set_step(h);

  ABMData<ABMD_DOUBLE> abm_data;

  ABMD_DOUBLE *hoho =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * ABMD_MAX_ORDER);
  for (int i = 0; i < ABMD_MAX_ORDER; i++) {
    hoho[i] = coefss<ABMD_DOUBLE>::PREDICTOR_COEFFS[i] * h;
  }

  abm_data.input = *abm;
  abm_data.rk4_h = rk4_h;
  abm_data.temp = rhs_temp;
  abm_data.queue = queue;
  abm_data.xs_delayed = xs_delayed;
  abm_data.xs_delayed_inner = xs_delayed_inner;
  abm_data.xs_delayed_tmp = xs_delayed_tmp;
  abm_data.dxs_delayed = dxs_delayed;
  abm_data.rk_memory = NULL;
  abm_data.inner_rk_memory = NULL;
  abm_data.hoho = &hoho[0];

  // Setting initial conditions for RK4 solution
  for (int i = 0; i < dim; i++) {
    rk4_sol[i] = init[i];
  }

  // Doing rk4_n RK4 steps
  for (int i = 1; i < rk4_n + 1; i++) {
    double t = t0 + rk4_h * (i - 1);
    rk_step(rhs_rk4, rk4_h, t, &rk4_sol[(i - 1) * dim], dim, &abm_data,
            &rk4_sol[i * dim], &rk4_rhss[(i - 1) * dim], &abm_data.rk_memory,
            METHOD_DOPRI8);
  }

  // Computing last ABMD_RHS
  int last = rk4_n * dim;
  rhs_rk4<ABMD_DOUBLE>(&rk4_sol[last], t0 + rk4_h * rk4_n, &rk4_rhss[last],
                       &abm_data);

  // Writing data from RK4 to the queue
  int k = 0;
  for (int i = 0; i < rk4_n + 1; i += RK_STEPS_IN_ABM) {
    ABMD_DOUBLE *sol_address = queue->push();
    memcpy(sol_address, &rk4_sol[i * dim], dim * sizeof(ABMD_DOUBLE));
    memcpy(&sol_address[dim], &rk4_rhss[i * dim], dim * sizeof(ABMD_DOUBLE));
    queue->update_diffs();
    queue->swap_diffs();
    k++;
  }
  free(rk4_rhss);

  int run_callback = abm->callback && abm->callback_t;
  // double *callback_state =
  //     (double *)malloc(sizeof(double) * dim);
  ABMD_DOUBLE *callback_state_l =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  double write_step = 1;
#ifdef ALL_STEPS
  write_step = h * hsgn;
#else
  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) 
#endif
  {
    queue->evaluate_x_all(*callback_t, callback_state_l);
    // for (int i = 0; i < dim; i++) {
    //   #ifdef NUMBER_DOUBLE_DOUBLE
    //   callback_state[i] = to_double(callback_state_l[i]);
    //   #else
    //   callback_state[i] = callback_state_l[i];
    //   #endif
    // }
    run_callback = abm->callback(callback_t, callback_state_l, abm->context, write_step);
  }

  free(rk4_sol);
  free(init);

  // Main ABM loop
  int start_index = rk4_i1 + 1;
  double step = (n - start_index) / 100.0;
  double last_step = 0;
  int steps = 1;
  for (int i = start_index; i < n; i++) {
    int cur_step = i - start_index;
    if (last_step <= cur_step) {
      last_step += step;
      std::cout << "Progress: " << steps << "/100\n";
      steps++;
    }
    double t = t0 + i * h;

    abm_data.predict();
    ABMD_DOUBLE *rhs_out = queue->peek_right_dx();
    get_delayed_states(&abm_data, t, 0, xs_delayed, dxs_delayed);
    rhs(queue->peek_right_x(), xs_delayed, dxs_delayed, t, rhs_out, &abm_data);

    queue->update_diffs();
    ABMD_DOUBLE *predicted = queue->peek_right_x();
    queue->backup_last_x();
    abm_data.correct(predicted);

    get_delayed_states(&abm_data, t, 1, xs_delayed, dxs_delayed);
    rhs(queue->peek_right_x(), xs_delayed, dxs_delayed, t, rhs_out, &abm_data);

    queue->update_diffs();
    queue->restore_last_x();
    abm_data.correct(nullptr);
    queue->swap_diffs();
#ifndef ALL_STEPS
    while (run_callback && (t - h) * hsgn <= *callback_t * hsgn &&
           *callback_t * hsgn <= t * hsgn) 
#endif
    {
      queue->evaluate_x_all(*callback_t, callback_state_l);
      // for (int j = 0; j < dim; j++) {
      //   callback_state[j] = (ABMD_DOUBLE)callback_state_l[j];
      // }
      run_callback = abm->callback(callback_t, callback_state_l, abm->context, write_step);
    }
  }

  #ifndef ALL_STEPS
  ABMD_DOUBLE *final_arr = queue->peek_right_x();
  #else
  ABMD_DOUBLE *final_arr = callback_state_l;
  #endif

  for (int i = 0; i < dim; i++) {
    abm->final_state[i] = (ABMD_DOUBLE)final_arr[i];
  }

  // destroy_abm_data(abm_data);
  // free(callback_state);
  free(callback_state_l);
}

template <typename ABMD_DOUBLE>
void moonToBarycenter(ABMD_DOUBLE *rr, ABMD_DOUBLE *res) {
  ABMD_DOUBLE earth_x = rr[earthNum * 6];
  ABMD_DOUBLE earth_y = rr[earthNum * 6 + 1];
  ABMD_DOUBLE earth_z = rr[earthNum * 6 + 2];
  ABMD_DOUBLE earth_vx = rr[earthNum * 6 + 3];
  ABMD_DOUBLE earth_vy = rr[earthNum * 6 + 4];
  ABMD_DOUBLE earth_vz = rr[earthNum * 6 + 5];

  ABMD_DOUBLE gcmoon_x = rr[moonNum * 6];
  ABMD_DOUBLE gcmoon_y = rr[moonNum * 6 + 1];
  ABMD_DOUBLE gcmoon_z = rr[moonNum * 6 + 2];
  ABMD_DOUBLE gcmoon_vx = rr[moonNum * 6 + 3];
  ABMD_DOUBLE gcmoon_vy = rr[moonNum * 6 + 4];
  ABMD_DOUBLE gcmoon_vz = rr[moonNum * 6 + 5];

  rr[moonNum * 6] = earth_x + gcmoon_x;
  rr[moonNum * 6 + 1] = earth_y + gcmoon_y;
  rr[moonNum * 6 + 2] = earth_z + gcmoon_z;
  rr[moonNum * 6 + 3] = earth_vx + gcmoon_vx;
  rr[moonNum * 6 + 4] = earth_vy + gcmoon_vy;
  rr[moonNum * 6 + 5] = earth_vz + gcmoon_vz;

  res[0] = gcmoon_x;
  res[1] = gcmoon_y;
  res[2] = gcmoon_z;
  res[3] = gcmoon_vx;
  res[4] = gcmoon_vy;
  res[5] = gcmoon_vz;
}

#define dist(x1, y1, z1, x2, y2, z2)                                \
  (sqrt(((x1) - (x2)) * ((x1)-x2) + ((y1) - (y2)) * ((y1) - (y2)) + \
        ((z1) - (z2)) * ((z1) - (z2))))
#define vecLen2(x, y, z) (((x) * (x)) + (y) * (y) + (z) * (z))

template <typename ABMD_DOUBLE>
void calc_invariants(ABMD_DOUBLE *rr, double *masses, int num,
                     ABMD_DOUBLE *new_center, ABMD_DOUBLE *inv_energy,
                     ABMD_DOUBLE *inv_momentum) {
  /// for energy fix need change moon center
  ABMD_DOUBLE *gcmoon = (ABMD_DOUBLE *)malloc(6 * sizeof(ABMD_DOUBLE));
  moonToBarycenter(rr, gcmoon);
#ifdef NUMBER_DOUBLE_DOUBLE
  ABMD_DOUBLE revLS2 =
      ABMD_DOUBLE("3.335661199676477669820576657401368279e-05");
#else
#if NUMBER_DOUBLE == 1
  ABMD_DOUBLE revLS2 = 3.335661199676477669820577e-05L;
#else
  ABMD_DOUBLE revLS2 = 3.335661199676477670e-005;
#endif
#endif
#ifdef RELATIVISTIC
  //------------relativistic barycentre-------------
  std::vector<ABMD_DOUBLE> b(3, 0);
  ABMD_DOUBLE sum_mu_star = 0;
  for (int i = 0; i < num; i++) {
    ABMD_DOUBLE mu_star = 0;
    ABMD_DOUBLE sum = 0;
    for (int j = 0; j < num; j++) {
      if (j != i) {
        ABMD_DOUBLE dist = dist(rr[j * 6], rr[j * 6 + 1], rr[j * 6 + 2],
                                rr[i * 6], rr[i * 6 + 1], rr[i * 6 + 2]);
        sum += ((ABMD_DOUBLE)masses[j]) / dist;
      }
    }
    mu_star = ((ABMD_DOUBLE)masses[i]) *
              (1.0 + 0.5 * revLS2 *
                         (vecLen2(rr[i * 6 + 3], rr[i * 6 + 4], rr[i * 6 + 5]) -
                          sum));
    sum_mu_star += mu_star;
    b[0] += mu_star * rr[i * 6];
    b[1] += mu_star * rr[i * 6 + 1];
    b[2] += mu_star * rr[i * 6 + 2];
  }

  b[0] /= sum_mu_star;
  b[1] /= sum_mu_star;
  b[2] /= sum_mu_star;
  b[0] *= 149597870700;
  b[1] *= 149597870700;
  b[2] *= 149597870700;

  *new_center = sqrt(vecLen2(b[0], b[1], b[2]));

//----------END of relativistic barycentre-------------
#else
  ABMD_DOUBLE pos_x = 0;
  ABMD_DOUBLE pos_y = 0;
  ABMD_DOUBLE pos_z = 0;
  ABMD_DOUBLE total_mass = 0;
  for (int i = 0; i < num; i++) {
    total_mass += masses[i];
    pos_x += rr[i * 6] * masses[i];
    pos_y += rr[i * 6 + 1] * masses[i];
    pos_z += rr[i * 6 + 2] * masses[i];
  }
  *new_center =
      sqrt(vecLen2(pos_x / total_mass, pos_y / total_mass, pos_z / total_mass));
#endif

  //-----------------------energy------------------------

  ABMD_DOUBLE total_energy = 0;
  for (int i = 0; i < num; i++) {
    total_energy +=
        masses[i] * (1 / sqrt(1 - revLS2 * vecLen2(rr[i * 6 + 3], rr[i * 6 + 4],
                                                   rr[i * 6 + 5])) -
                     1);
  }
  total_energy *= revLS2;
  *inv_energy = total_energy;

  //--------------------END of energy--------------------

  //-----------------momentum---------------------------
  ABMD_DOUBLE inv_v_x = 0;
  ABMD_DOUBLE inv_v_y = 0;
  ABMD_DOUBLE inv_v_z = 0;
  for (int i = 0; i < num; i++) {
    ABMD_DOUBLE nu_m =
        masses[i] /
        sqrt(1 - revLS2 * vecLen2(rr[i * 6 + 3], rr[i * 6 + 4], rr[i * 6 + 5]));
    inv_v_x += nu_m * rr[i * 6 + 3];
    inv_v_y += nu_m * rr[i * 6 + 4];
    inv_v_z += nu_m * rr[i * 6 + 5];
  }
  inv_momentum[0] = inv_v_x;
  inv_momentum[1] = inv_v_y;
  inv_momentum[2] = inv_v_z;

  //------------------END of momentum-------------------

  rr[moonNum * 6] = gcmoon[0];
  rr[moonNum * 6 + 1] = gcmoon[1];
  rr[moonNum * 6 + 2] = gcmoon[2];
  rr[moonNum * 6 + 3] = gcmoon[3];
  rr[moonNum * 6 + 4] = gcmoon[4];
  rr[moonNum * 6 + 5] = gcmoon[5];

  free(gcmoon);
}

template <typename ABMD_DOUBLE>
int callback_there(double *t, ABMD_DOUBLE *state, void *context, double h) {
  ContextData<ABMD_DOUBLE> *abm_test = (ContextData<ABMD_DOUBLE> *)context;
#ifdef SAVE_STEPS
  for (int planet : DIFF_PLANETS) {
    abm_test->f << std::setprecision(20) << state[6 * planet] << " "
                << state[6 * planet + 1] << " " << state[6 * planet + 2] << " "
                << state[6 * planet + 3] << " " << state[6 * planet + 4] << " "
                << state[6 * planet + 5] << " ";
  }
  abm_test->f << std::endl;

#endif

#ifdef SAVE_INV
  ABMD_DOUBLE inv_barycentre = 0;
  ABMD_DOUBLE inv_energy = 0;
  ABMD_DOUBLE inv_momentum[3];
  calc_invariants(state, abm_test->objects->masses.data(),
                  abm_test->objects->n_objects, &inv_barycentre, &inv_energy,
                  inv_momentum);
  abm_test->fi << " " << inv_barycentre << "\n";
  abm_test->fen << " " << std::setprecision(40) << inv_momentum[0] << " "
                << inv_momentum[1] << " " << inv_momentum[2] << "\n";
#endif

  abm_test->i++;
  t[0] += h;

  return 1;
}

template <typename ABMD_DOUBLE>
int callback_back(double *t, ABMD_DOUBLE *state, void *context, double h) {
  ContextData<ABMD_DOUBLE> *abm_test = (ContextData<ABMD_DOUBLE> *)context;
  int dim = abm_test->dim;

#ifdef SAVE_STEPS
  for (int planet : DIFF_PLANETS) {
    abm_test->fb << std::setprecision(20) << state[6 * planet] << " "
                 << state[6 * planet + 1] << " " << state[6 * planet + 2] << " "
                 << state[6 * planet + 3] << " " << state[6 * planet + 4] << " "
                 << state[6 * planet + 5] << " ";
  }
  abm_test->fb << std::endl;
#endif

  abm_test->i++;
  t[0] -= h;

  return 1;
}

template <typename ABMD_DOUBLE>
void ABMD_calc_diff(std::vector<ABMD_DOUBLE> &x, std::vector<double> masses,
                    double h, double t1) {
  int order = 13;
  double t0 = 0;
  int dim = 6 * masses.size();

  int n = (int)(1 + (t1 - t0) / h);
  int sol_size = n;
  std::cout << "sol_size : " << sol_size << std::endl;
  double callback_t = 0;

  ObjectsData<ABMD_DOUBLE> *objects = new ObjectsData<ABMD_DOUBLE>(masses);
  ContextData<ABMD_DOUBLE> abm_test{
      .objects = objects,
      .callback_t = &callback_t,
      .center = nullptr,
      .i = 0,
      .dim = dim,
      .f = std::ofstream("res_out.txt"),     // fopen("res_out.txt", "wt"),
      .fb = std::ofstream("res_b_out.txt"),  // fopen("res_b_out.txt", "wt")};
      .fi = std::ofstream("res_inv.txt"),
      .fen = std::ofstream("energy_inv.txt")};
  ABMD<ABMD_DOUBLE> *abm =
      abmd_create(pointmassesCalculateXdot_tmp, dim, t0, t1, h, x.data());

  abm->callback = callback_there;
  abm->callback_t = &callback_t;
  abm->context = &abm_test;

  ABMD_DOUBLE *prev_final_state =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * dim);

  auto start = std::chrono::high_resolution_clock::now();
  ABMD_run(abm);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";
  for (int i = 0; i < dim; i++) {
    prev_final_state[i] = abm->final_state[i];
  }

  abmd_destroy(abm);

  abm = abmd_create(pointmassesCalculateXdot_tmp, dim, t1, t0, h,
                    prev_final_state);
  abm->callback = callback_back;
  abm->callback_t = &callback_t;
  abm->context = &abm_test;
  callback_t = t1;
  abm_test.i = 0;

  ABMD_run(abm);
  abmd_destroy(abm);

#ifdef SAVE_DIFF
  std::ofstream diff_file = std::ofstream("Mars_diff.txt");
  std::ifstream b_f("res_b_out.txt");
  std::ifstream o_f("res_out.txt");
  std::string str;
  std::vector<std::string> vec_str;

  while (std::getline(b_f, str)) {
    if (str.size() > 3) {
      vec_str.push_back(str);
    }
  }

  int j = vec_str.size() - 1;

  while (std::getline(o_f, str)) {
    if (str.size() > 3) {
      std::vector<std::string> res_split = split(str, " ");
      std::vector<std::string> res_b_split = split(vec_str[j], " ");
      for (int i = 0; i < DIFF_PLANETS.size(); i++) {
        double x1 = std::stod(res_split[6 * i]);
        double x2 = std::stod(res_b_split[6 * i]);

        double y1 = std::stod(res_split[6 * i + 2]);
        double y2 = std::stod(res_b_split[6 * i + 2]);

        double z1 = std::stod(res_split[6 * i + 4]);
        double z2 = std::stod(res_b_split[6 * i + 4]);

        diff_file << " " << std::setprecision(30)
                  << sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
      }
      diff_file << std::endl;
      j--;
    }
  }
#endif
}

#endif  // ABMD_H