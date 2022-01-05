#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ABMD_MAX_ORDER 19
#define ABMD_DEFAULT_ORDER 11

enum {
  METHOD_DOPRI8 = 1,
  METHOD_RK4 = 2
};

template <typename ABMD_DOUBLE>
using ABMD_RHS = void (*)(ABMD_DOUBLE x[], double t, ABMD_DOUBLE *out, void *context);

template <typename ABMD_DOUBLE>
using ABMD_RHSD = void (*)(ABMD_DOUBLE x[], ABMD_DOUBLE xs_delayed[], ABMD_DOUBLE dxs_delayed[],
                    double t, ABMD_DOUBLE *out, void *context);

template <typename ABMD_DOUBLE>
struct ABMD {
  ABMD_RHS<ABMD_DOUBLE> f1;
  ABMD_RHSD<ABMD_DOUBLE> f2;
  int dim;
  double t0;
  double t1;
  double h;
  double *init;
  double *delays;
  int ndelays;
  int abm_order;
  int delays_poly_degree;
  int pointsave_poly_degree;
  double *final_state;
  void *context;
  void(*init_call)(ABMD_DOUBLE[], void*);
  double *callback_t;
  int(*callback)(double *t, double state[], void *context);
  int **delayed_idxs;
  int *delayed_idxs_lens;
  int *dx_delays_idxs;
  int dx_delays_len;
  char *error;
};

template <typename ABMD_DOUBLE>
ABMD<ABMD_DOUBLE> *abmd_create(ABMD_RHS<ABMD_DOUBLE> f, int dim, double t0, double t1, double h, double *init) {
  ABMD<ABMD_DOUBLE> *abm = (ABMD<ABMD_DOUBLE> *) malloc(sizeof(ABMD<ABMD_DOUBLE>));
  double *final_state = (double *) malloc(sizeof(double) * dim);
  char *error = (char *) malloc(256 * sizeof(char));
  *abm = (ABMD<ABMD_DOUBLE>) {
          .f1=f,
          .f2=NULL,
          .dim=dim,
          .t0=t0,
          .t1=t1,
          .h=h,
          .init=init,
          .delays=NULL,
          .ndelays=0,
          .abm_order=ABMD_DEFAULT_ORDER,
          .delays_poly_degree=ABMD_DEFAULT_ORDER,
          .pointsave_poly_degree=ABMD_DEFAULT_ORDER,
          .final_state=final_state,
          .context=NULL,
          .init_call=NULL,
          .callback_t=NULL,
          .callback=NULL,
          .delayed_idxs=NULL,
          .delayed_idxs_lens=NULL,
          .dx_delays_idxs=NULL,
          .dx_delays_len=0,
          .error=error
  };
  return abm;
}

template <typename ABMD_DOUBLE>
void abmd_destroy(ABMD<ABMD_DOUBLE> *abm) {
  free(abm->delays);
  free(abm->final_state);
  for (int i = 0; i < abm->ndelays; i++) {
    free(abm->delayed_idxs[i]);
  }
  free(abm->delayed_idxs);
  free(abm->delayed_idxs_lens);
  free(abm->dx_delays_idxs);
  free(abm->error);
  free(abm);
}

template <typename ABMD_DOUBLE>
struct Queue {
  double t0, h;
  int head, tail, size, dim;
  int capacity, block_size;
  ABMD_DOUBLE* _array;
  ABMD_DOUBLE **xarray, **dxarray;
  ABMD_DOUBLE *x_backup;
  ABMD_DOUBLE* diffs_r;
  ABMD_DOUBLE* diffs_w;
  ABMD_DOUBLE* last_diff;
  int delays_poly_degree;
  int pointsave_poly_degree;
  double *lgr_delay_ws;
  double *lgr_delay_ws_nolast;
  double *lgr_pointsave_ws;
  ABMD_DOUBLE *lgr_nom;

  Queue(int capacity, int dim) {
    int block_size = 2 * dim;
    //Queue<ABMD_DOUBLE>* queue = (Queue<ABMD_DOUBLE> *) malloc(sizeof(Queue<ABMD_DOUBLE>));
    t0 = 0;
    h = 1;
    capacity = capacity;
    block_size = block_size;
    dim = dim;
    head = size = 0;
    tail = capacity - 1;
    _array = (ABMD_DOUBLE *) malloc(capacity * block_size * sizeof(ABMD_DOUBLE));
    xarray = (ABMD_DOUBLE **) malloc(capacity * sizeof(ABMD_DOUBLE *));
    dxarray = (ABMD_DOUBLE **) malloc(capacity * sizeof(ABMD_DOUBLE *));
    x_backup = (ABMD_DOUBLE *) malloc(dim * sizeof(ABMD_DOUBLE));
    for (int i = 0; i < capacity; i++) {
      xarray[i] = &this->_array[i * block_size];
      dxarray[i] = &this->_array[i * block_size + dim];
    }

    diffs_r = (ABMD_DOUBLE *) malloc((capacity - 1) * dim * sizeof(ABMD_DOUBLE));
    diffs_w = (ABMD_DOUBLE *) malloc((capacity - 1) * dim * sizeof(ABMD_DOUBLE));
    last_diff = (ABMD_DOUBLE *) malloc(dim * sizeof(ABMD_DOUBLE));
    delays_poly_degree = capacity - 1;
    pointsave_poly_degree = capacity - 1;
    lgr_delay_ws = (double *) malloc(capacity * sizeof(double));
    lgr_delay_ws_nolast = (double *) malloc((capacity - 1) * sizeof(double));
    lgr_pointsave_ws = (double *) malloc(capacity * sizeof(double));
    lgr_nom = (ABMD_DOUBLE *) malloc(dim * sizeof(ABMD_DOUBLE));
  }

  ~Queue() {
    free(_array);
    free(xarray);
    free(dxarray);
    free(x_backup);
    free(diffs_r);
    free(diffs_w);
    free(last_diff);
    free(lgr_delay_ws);
    free(lgr_delay_ws_nolast);
    free(lgr_pointsave_ws);
    free(lgr_nom);
  }

  int is_full() {
    return size == capacity;
  }

  int is_empty() {
    return size == 0;
  }

  // int get_capacity() {
  //   return capacity;
  // }

  // void set_t0(double new_t0) {
  //   this->t0 = new_t0;
  // }

  // void qset_delays_poly_degree(Queue *q, int deg) {
  //   q->delays_poly_degree = deg;
  // }

  // void qset_pointsave_poly_degree(Queue *q, int deg) {
  //   q->pointsave_poly_degree = deg;
  // }

  double _compute_wj(int j, int len) {
    /*
    * Computes j-th barycentric Lagrange weight (eq. 3.2)
    */
    double w = 1;
    double h = this->h;
    for (int k = 0; k < len; k++) {
      if (k == j) continue;
      w /= (j - k) * h;  // ts[j] - ts[k] simplified
    }
    return w;
  }

  void _initialize_weights() {
    int n = delays_poly_degree + 1;
    for (int i = 0; i < n - 1; i++) {
      lgr_delay_ws[i] = _compute_wj(i, n);
      lgr_delay_ws_nolast[i] = _compute_wj(i, n - 1);
    }
    lgr_delay_ws[n - 1] = _compute_wj(n - 1, n);

    int m = pointsave_poly_degree + 1;
    for (int i = 0; i < m; i++) {
      lgr_pointsave_ws[i] = _compute_wj(i, m);
    }
  }

  void set_step(double h) {
    this->h = h;
    _initialize_weights();
  }

  double _get_t_index(double t, int last_known) {
    double t0 = this->t0;
    double h = this->h;
    int n = this->capacity;
    if (!last_known) {
      n -= 1;
    }
    double t1 = t0 + (n - 1) * h;
    int hsgn = (t1 > t0) - (t1 < t0);
    t0 *= hsgn;
    t1 *= hsgn;
    t *= hsgn;
    h *= hsgn;

    if (t0 <= t && t <= t1) {
      return (t - t0) / h;
    }

    return -1;
  }

  ABMD_DOUBLE* get_x(int block_idx) {
    return xarray[(head + block_idx) % capacity];
  }

  ABMD_DOUBLE* get_dx(int block_idx) {
    return dxarray[(head + block_idx) % capacity];
  }

  void backup_last_x() {
    xarray[tail] = x_backup;
  }

  void restore_last_x() {
    int i = (tail * block_size) % (capacity * block_size);
    xarray[tail] = &_array[i];
  }

  void _evaluate(double t, int *idxs, int idxs_len, int n_points,
                double *ws, int last_known, ABMD_DOUBLE *(*get)(int),
                ABMD_DOUBLE *out) {

    double t_idx = _get_t_index(t, last_known);
    if (t_idx != -1 && fmod(t_idx, 1) < 1e-13) {
      ABMD_DOUBLE *x = get((int) round(t_idx));
      if (idxs == NULL) {
        memcpy(out, x, idxs_len * sizeof(ABMD_DOUBLE));
        return;
      }
      for (int j = 0; j < idxs_len; j++) {
        out[j] = x[idxs[j]];
      }
      return;
    }

    int left = capacity - n_points;
    if (!last_known) {
      left -= 1;
    }
    if (t_idx != -1 && t_idx < left) {
      left = (int) t_idx;
    }

    ABMD_DOUBLE denom = 0;
    ABMD_DOUBLE coefs[ABMD_MAX_ORDER + 1];
    ABMD_DOUBLE *xs[ABMD_MAX_ORDER + 1];
    
    for (int i = 0; i < n_points; i++) {
      double tti = t - (t0 + (i + left) * h);
      coefs[i] = ws[i] / tti;
      denom += coefs[i];
    }
    for (int i = 0; i < n_points; i++)
    {
      coefs[i] /= denom;
      xs[i] = get(i + left);
    }

    if (idxs == NULL) {
      for (int j = 0; j < idxs_len; j++) {
        ABMD_DOUBLE res = 0;
        for (int i = 0; i < n_points; i++) {
          res += coefs[i] * xs[i][j];
        }
        out[j] = res;
      }
    } else {
      for (int j = 0; j < idxs_len; j++) {
        ABMD_DOUBLE res = 0;
        for (int i = 0; i < n_points; i++) {
          res += coefs[i] * xs[i][idxs[j]];
        }
        out[j] = res;
      }
    }
  }

  void evaluate_x_all(double t, ABMD_DOUBLE *out) {
    _evaluate(t, NULL, dim, pointsave_poly_degree + 1,
              lgr_pointsave_ws, 1, get_x, out);
  }

  void evaluate_x_idxs(double t, int *idxs, int idxs_len, ABMD_DOUBLE *out) {
    _evaluate(t, idxs, idxs_len, delays_poly_degree + 1,
              lgr_delay_ws, 1, get_x, out);
  }

  void evaluate_dx(double t, int *idxs, int idxs_len,
                  int last_known, ABMD_DOUBLE *out) {
    int deg = delays_poly_degree;
    double *ws = lgr_delay_ws;
    int n_points = deg + 1;
    if (!last_known && deg == capacity - 1) {
      ws = lgr_delay_ws_nolast;
      n_points -= 1;
    }
    _evaluate(t, idxs, idxs_len, n_points, ws, last_known, get_dx, out);
  }

  ABMD_DOUBLE* push() {
    if (is_full())
      return NULL;
    tail = (tail + 1) % capacity;
    size += 1;
    return xarray[tail];
  }

  ABMD_DOUBLE* pop() {
    if (is_empty())
      return NULL;
    ABMD_DOUBLE *address = xarray[head];
    head = (head + 1) % capacity;
    size -= 1;
    t0 += h;
    return address;
  }

  ABMD_DOUBLE* peek_left() {
    if (is_empty())
      return NULL;
    return xarray[head];
  }

  ABMD_DOUBLE* peek_right_x() {
    if (is_empty())
      return NULL;
    return xarray[tail];
  }

  ABMD_DOUBLE* peek_right_dx() {
    if (is_empty())
      return NULL;
    return dxarray[tail];
  }

  ABMD_DOUBLE* get_diffs_r() {
    return diffs_r;
  }

  ABMD_DOUBLE* get_diffs_w() {
    return diffs_w;
  }

  ABMD_DOUBLE* get_last_diff() {
    return last_diff;
  }

  void swap_diffs() {
    ABMD_DOUBLE *diffs_r = this->diffs_r;
    this->diffs_r = this->diffs_w;
    this->diffs_w = diffs_r;
  }

  void update_diffs() {
    int size = this->size;
    int diffs_len = this->capacity - 1;
    int dim = this->dim;

    ABMD_DOUBLE *right = peek_right_dx();
    ABMD_DOUBLE *diffs_r = this->diffs_r;
    ABMD_DOUBLE *diffs_w = this->diffs_w;

    if (!is_full()) {
      for (int i = 0; i < dim; i++) {
        ABMD_DOUBLE new_d = *right++;
        for (int j = 0; j < size - 1; j++) {
          *diffs_w++ = new_d;
          new_d -= *diffs_r++;
        }
        *diffs_w = new_d;
        diffs_w += diffs_len - size + 1;
        diffs_r += diffs_len - size + 1;
      }
      return;
    }

    ABMD_DOUBLE *last_diff = last_diff;

    for (int i = 0; i < dim; i++) {
      ABMD_DOUBLE new_d = *right++;
      for (int j = 0; j < size - 1; j++) {
        *diffs_w++ = new_d;
        new_d -= *diffs_r++;
      }
      *last_diff++ = new_d;
    }
  }
};
