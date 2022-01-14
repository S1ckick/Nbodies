#ifndef ABMD_H
#define ABMD_H

#include <cstdlib>
#include <string>

#include "../pointmasses.h"
#include "abmd_rk.h"

template <typename ABMD_DOUBLE>
ABMD_DOUBLE PREDICTOR_COEFFS[ABMD_MAX_ORDER] = {
    1.0L,
    1.0 / 2.0L,
    5.0 / 12.0L,
    3.0 / 8.0L,
    251.0 / 720.0L,
    95.0 / 288.0L,
    19087.0 / 60480.0L,
    5257.0 / 17280.0L,
    1070017.0 / 3628800.0L,
    25713.0 / 89600.0L,
    26842253.0 / 95800320.0L,
    4777223.0 / 17418240.0L,
    703604254357.0 / 2615348736000.0L,
    106364763817.0 / 402361344000.0L,
    1166309819657.0 / 4483454976000.0L,
    2.5221445e7 / 9.8402304e7L,
    8.092989203533249e15 / 3.201186852864e16L,
    8.5455477715379e13 / 3.4237292544e14L,
    1.2600467236042756559e19 / 5.109094217170944e19L};

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

    ABMD_DOUBLE *diffs = queue->get_diffs_r();
    for (int i = 0; i < dim; i++) {
      ABMD_DOUBLE *c = PREDICTOR_COEFFS<ABMD_DOUBLE>;
      for (int j = 0; j < abm_order; j++) {
        ABMD_DOUBLE ch = *c++ * h;
        out[i] += *diffs++ * ch;
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
    ABMD_DOUBLE ch = h * PREDICTOR_COEFFS<ABMD_DOUBLE>[abm_order];
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
  ABMData<ABMD_DOUBLE> *data = (ABMData<ABMD_DOUBLE> *)abm_data;
  int dim = data->input.dim;
  ABMD_DOUBLE *temp = data->temp;

  memset(out, 0, sizeof(ABMD_DOUBLE) * dim);

  if (data->input.f1 == NULL) {
    data->input.f2(x, xs_delayed, dxs_delayed, t, out, data->input.context);
    return;
  }

  data->input.f1(x, t, out, data->input.context);

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
  ABMD_DOUBLE *callback_state =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * dim);
  ABMD_DOUBLE *callback_state_l =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * dim);
  double *callback_t = abm->callback_t;

  while (run_callback && *callback_t * hsgn <= rk4_t1 * hsgn) {
    queue->evaluate_x_all(*callback_t, callback_state_l);
    for (int i = 0; i < dim; i++) {
      callback_state[i] = (double)callback_state_l[i];
    }
    run_callback = abm->callback(callback_t, callback_state, abm->context);
  }

  free(rk4_sol);
  free(init);

  // Main ABM loop
  int start_index = rk4_i1 + 1;
  for (int i = start_index; i < n; i++) {
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

    while (run_callback && (t - h) * hsgn <= *callback_t * hsgn &&
           *callback_t * hsgn <= t * hsgn) {
      queue->evaluate_x_all(*callback_t, callback_state_l);
      for (int j = 0; j < dim; j++) {
        callback_state[j] = (ABMD_DOUBLE)callback_state_l[j];
      }
      run_callback = abm->callback(callback_t, callback_state, abm->context);
    }
  }

  for (int i = 0; i < dim; i++) {
    abm->final_state[i] = (ABMD_DOUBLE)queue->peek_right_x()[i];
  }

  // destroy_abm_data(abm_data);
  free(callback_state);
  free(callback_state_l);
}

template <typename ABMD_DOUBLE>
int callback_there(double *t, ABMD_DOUBLE *state, void *context) {
  ContextData<ABMD_DOUBLE> *abm_test = (ContextData<ABMD_DOUBLE> *)context;
  int dim = abm_test->dim;
  fprintf(abm_test->f, "%lf ", *t);
  for (int i = 0; i < dim; i++) {
    fprintf(abm_test->f, " %.16le", state[i]);
  }
  fprintf(abm_test->f, "\n");
  t[0] += 1 / 32.0;
  return 1;
}

template <typename ABMD_DOUBLE>
int callback_back(double *t, ABMD_DOUBLE *state, void *context) {
  ContextData<ABMD_DOUBLE> *abm_test = (ContextData<ABMD_DOUBLE> *)context;
  int dim = abm_test->dim;
  memcpy(&abm_test->sol_back[abm_test->i * dim], state,
         dim * sizeof(ABMD_DOUBLE));
  abm_test->i++;
  t[0] -= 1 / 32.0;
  //  t[0] -= 1 / 16.0;
  return 1;
}

template <typename ABMD_DOUBLE>
void ABMD_calc_diff(std::vector<ABMD_DOUBLE> &x,
                    std::vector<ABMD_DOUBLE> masses, ABMD_DOUBLE h) {
  int order = 11;
  double t0 = 0;
  double t1 = 365;
  // double h = 1 / 16.0;
  int dim = 6 * masses.size();

  int n = (int)(1 + (t1 - t0) / h);
  int sol_size = 2 * n - 1;
  //  int sol_size = n;
  ABMD_DOUBLE *sol =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * sol_size * dim);
  ABMD_DOUBLE *sol_back =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * sol_size * dim);
  double callback_t = 0;

  ObjectsData<ABMD_DOUBLE> *objects = new ObjectsData<ABMD_DOUBLE>(masses);
  ContextData<ABMD_DOUBLE> abm_test{.objects = objects,
                                    .sol = sol,
                                    .sol_back = sol_back,
                                    .callback_t = &callback_t,
                                    .i = 0,
                                    .dim = dim,
                                    .f = fopen("res_out.txt", "wt")};
  ABMD<ABMD_DOUBLE> *abm =
      abmd_create(pointmassesCalculateXdot_tmp, dim, t0, t1, h, x.data());

  abm->callback = callback_there;
  abm->callback_t = &callback_t;
  abm->context = &abm_test;

  ABMD_run(abm);
  // printf("Final: %e %e\n", abmd_get_final_state(abm)[0],
  // abmd_get_final_state(abm)[2]);
  abmd_destroy(abm);

  ABMD_DOUBLE *sol_reversed =
      (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * sol_size * dim);

  for (int i = 0; i < sol_size; i++) {
    memcpy(&sol_reversed[i * dim], &sol[(sol_size - i - 1) * dim],
           sizeof(ABMD_DOUBLE) * dim);
    sol_reversed[i * dim + 1] *= -1;
    sol_reversed[i * dim + 3] *= -1;
  }

  /*
    abm = abmd_create(pointmassesCalculateXdot_tmp, dim, t1, t0, h,
                      &sol[(sol_size - 1) * dim]);
    abm->callback = callback_back;
    abm->callback_t = &callback_t;
    abm->context = &abm_test;
    callback_t = t1;
    abm_test.i = 0;

    ABMD_run(abm);
    abmd_destroy(abm);
    ABMD_DOUBLE *diff = (ABMD_DOUBLE *)malloc(sizeof(ABMD_DOUBLE) * sol_size *
    2);

    for (int i = 0; i < sol_size; i++) {
      ABMD_DOUBLE x1 = sol_reversed[i * dim];
      ABMD_DOUBLE x2 = sol_back[i * dim];

      ABMD_DOUBLE y1 = sol_reversed[i * dim + 2];
      ABMD_DOUBLE y2 = sol_back[i * dim + 2];
      diff[i * 2] = t0 + i * h;
      diff[i * 2 + 1] = (ABMD_DOUBLE)sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
    }
  */
  free(sol);
  free(sol_reversed);
  // free(sol_back);
  // free(diff);
}

#endif  // ABMD_H