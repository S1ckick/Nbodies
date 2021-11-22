#include <cstdlib>
#include "types.h"

double c2  = 0.526001519587677318785587544488E-01,
       c3  = 0.789002279381515978178381316732E-01,
       c4  = 0.118350341907227396726757197510E+00,
       c5  = 0.281649658092772603273242802490E+00,
       c6  = 0.333333333333333333333333333333E+00,
       c7  = 0.25E+00,
       c8  = 0.307692307692307692307692307692E+00,
       c9  = 0.651282051282051282051282051282E+00,
       c10 = 0.6E+00,
       c11 = 0.857142857142857142857142857142E+00;

/// fixfixfix
template <typename ABMD_DOUBLE>
struct coefs {
ABMD_DOUBLE b1 =   5.42937341165687622380535766363E-2L,
       b6 =   4.45031289275240888144113950566E0L,
       b7 =   1.89151789931450038304281599044E0L,
       b8 =  -5.8012039600105847814672114227E0L,
       b9 =   3.1116436695781989440891606237E-1L,
       b10 = -1.52160949662516078556178806805E-1L,
       b11 =  2.01365400804030348374776537501E-1L,
       b12 =  4.47106157277725905176885569043E-2L,


       a21 =    5.26001519587677318785587544488E-2L,
       a31 =    1.97250569845378994544595329183E-2L,
       a32 =    5.91751709536136983633785987549E-2L,
       a41 =    2.95875854768068491816892993775E-2L,
       a43 =    8.87627564304205475450678981324E-2L,
       a51 =    2.41365134159266685502369798665E-1L,
       a53 =   -8.84549479328286085344864962717E-1L,
       a54 =    9.24834003261792003115737966543E-1L,
       a61 =    3.7037037037037037037037037037E-2L,
       a64 =    1.70828608729473871279604482173E-1L,
       a65 =    1.25467687566822425016691814123E-1L,
       a71 =    3.7109375E-2L,
       a74 =    1.70252211019544039314978060272E-1L,
       a75 =    6.02165389804559606850219397283E-2L,
       a76 =   -1.7578125E-2L,

       a81 =    3.70920001185047927108779319836E-2L,
       a84 =    1.70383925712239993810214054705E-1L,
       a85 =    1.07262030446373284651809199168E-1L,
       a86 =   -1.53194377486244017527936158236E-2L,
       a87 =    8.27378916381402288758473766002E-3L,
       a91 =    6.24110958716075717114429577812E-1L,
       a94 =   -3.36089262944694129406857109825E0L,
       a95 =   -8.68219346841726006818189891453E-1L,
       a96 =    2.75920996994467083049415600797E1L,
       a97 =    2.01540675504778934086186788979E1L,
       a98 =   -4.34898841810699588477366255144E1L,
       a101 =   4.77662536438264365890433908527E-1L,
       a104 =  -2.48811461997166764192642586468E0L,
       a105 =  -5.90290826836842996371446475743E-1L,
       a106 =   2.12300514481811942347288949897E1L,
       a107 =   1.52792336328824235832596922938E1L,
       a108 =  -3.32882109689848629194453265587E1L,
       a109 =  -2.03312017085086261358222928593E-2L,

       a111 =  -9.3714243008598732571704021658E-1L,
       a114 =   5.18637242884406370830023853209E0L,
       a115 =   1.09143734899672957818500254654E0L,
       a116 =  -8.14978701074692612513997267357E0L,
       a117 =  -1.85200656599969598641566180701E1L,
       a118 =   2.27394870993505042818970056734E1L,
       a119 =   2.49360555267965238987089396762E0L,
       a1110 = -3.0467644718982195003823669022E0L,
       a121 =   2.27331014751653820792359768449E0L,
       a124 =  -1.05344954667372501984066689879E1L,
       a125 =  -2.00087205822486249909675718444E0L,
       a126 =  -1.79589318631187989172765950534E1L,
       a127 =   2.79488845294199600508499808837E1L,
       a128 =  -2.85899827713502369474065508674E0L,
       a129 =  -8.87285693353062954433549289258E0L,
       a1210 =  1.23605671757943030647266201528E1L,
       a1211 =  6.43392746015763530355970484046E-1L;
};

template <typename ABMD_DOUBLE>
void rk4_step(ABMD_RHS<ABMD_DOUBLE> f, double h, double t, ABMD_DOUBLE *x, int dim, void *context,
              ABMD_DOUBLE *out, ABMD_DOUBLE *rhs_out, ABMD_DOUBLE **memory) {

  if (*memory == NULL) {
    *memory = (ABMD_DOUBLE *) malloc(sizeof(ABMD_DOUBLE) * (5 * dim));
  }
  ABMD_DOUBLE *data = *memory;

  ABMD_DOUBLE *k1 = data;
  ABMD_DOUBLE *k2 = &data[dim];
  ABMD_DOUBLE *k3 = &data[dim * 2];
  ABMD_DOUBLE *k4 = &data[dim * 3];

  ABMD_DOUBLE *input = &data[dim * 4];

  memcpy(input, x, dim * sizeof(ABMD_DOUBLE));
  f(input, t, k1, context);

  if (rhs_out != NULL) {
    memcpy(rhs_out, k1, dim * sizeof(ABMD_DOUBLE));
  }

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * k1[i] / 2;
  }
  f(input, t + h / 2, k2, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * k2[i] / 2;
  }
  f(input, t + h / 2, k3, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * k3[i];
  }
  f(input, t + h, k4, context);

  for (int i = 0; i < dim; i++) {
    out[i] = x[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
  }
}

template <typename ABMD_DOUBLE>
void dopri8_step(ABMD_RHS<ABMD_DOUBLE> f, double h, double t, ABMD_DOUBLE *x, int dim, void *context,
                 ABMD_DOUBLE *out, ABMD_DOUBLE *rhs_out, ABMD_DOUBLE **memory) {
  coefs<ABMD_DOUBLE> cfs();
                
  if (*memory == NULL) {
    *memory = (ABMD_DOUBLE *) malloc(sizeof(ABMD_DOUBLE) * (11 * dim));
  }
  ABMD_DOUBLE *data = *memory;

  ABMD_DOUBLE *k1 = data;
  ABMD_DOUBLE *k2 = &data[dim];
  ABMD_DOUBLE *k3 = &data[dim * 2];
  ABMD_DOUBLE *k4 = &data[dim * 3];
  ABMD_DOUBLE *k5 = &data[dim * 4];
  ABMD_DOUBLE *k6 = &data[dim * 5];
  ABMD_DOUBLE *k7 = &data[dim * 6];
  ABMD_DOUBLE *k8 = &data[dim * 7];
  ABMD_DOUBLE *k9 = &data[dim * 8];
  ABMD_DOUBLE *k10 = &data[dim * 9];

  ABMD_DOUBLE *input = &data[dim * 10];

  memcpy(input, x, dim * sizeof(ABMD_DOUBLE));
  f(input, t, k1, context);

  if (rhs_out != NULL) {
    memcpy(rhs_out, k1, dim * sizeof(ABMD_DOUBLE));
  }

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * cfs.a21 * k1[i];
  }
  f(input, t + c2 * h, k2, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a31 * k1[i] + cfs.a32 * k2[i]);
  }
  f(input, t + c3 * h, k3, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a41 * k1[i] + cfs.a43 * k3[i]);
  }
  f(input, t + c4 * h, k4, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a51 * k1[i] + cfs.a53 * k3[i] + cfs.a54 * k4[i]);
  }
  f(input, t + c5 * h, k5, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a61 * k1[i] + cfs.a64 * k4[i] + cfs.a65 * k5[i]);
  }
  f(input, t + c6 * h, k6, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a71 * k1[i] + cfs.a74 * k4[i] +
                           cfs.a75 * k5[i] + cfs.a76 * k6[i]);
  }
  f(input, t + c7 * h, k7, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a81 * k1[i] + cfs.a84 * k4[i] + cfs.a85 * k5[i] +
                           cfs.a86 * k6[i] + cfs.a87 * k7[i]);
  }
  f(input, t + c8 * h, k8, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a91 * k1[i] + cfs.a94 * k4[i] + cfs.a95 * k5[i] +
                           cfs.a96 * k6[i] + cfs.a97 * k7[i] + cfs.a98 * k8[i]);
  }
  f(input, t + c9 * h, k9, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a101 * k1[i] + cfs.a104 * k4[i] + cfs.a105 * k5[i] +
                           cfs.a106 * k6[i] + cfs.a107 * k7[i] +
                           cfs.a108 * k8[i] + cfs.a109 * k9[i]);
  }
  f(input, t + c10 * h, k10, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a111 * k1[i] + cfs.a114 * k4[i] + cfs.a115 * k5[i] +
                           cfs.a116 * k6[i] + cfs.a117 * k7[i] + cfs.a118 * k8[i] +
                           cfs.a119 * k9[i] + cfs.a1110 * k10[i]);
  }
  f(input, t + c11 * h, k2, context);

  for (int i = 0; i < dim; i++) {
    input[i] = x[i] + h * (cfs.a121 * k1[i] + cfs.a124 * k4[i] + cfs.a125 * k5[i] +
                             cfs.a126 * k6[i] + cfs.a127 * k7[i] + cfs.a128 * k8[i] +
                             cfs.a129 * k9[i] + cfs.a1210 * k10[i] + cfs.a1211 * k2[i]);
  }
  f(input, t + h, k3, context);

  for (int i = 0; i < dim; i++) {
    k4[i] = cfs.b1 * k1[i] + cfs.b6 * k6[i] + cfs.b7 * k7[i] + cfs.b8 * k8[i] + cfs.b9 * k9[i] +
            cfs.b10 * k10[i] + cfs.b11 * k2[i] + cfs.b12 * k3[i];
    k5[i] = x[i] + h * k4[i];
  }

  memcpy(out, k5, dim * sizeof(ABMD_DOUBLE));
}

template <typename ABMD_DOUBLE>
void rk_step(ABMD_RHS<ABMD_DOUBLE> f, double h, double t, ABMD_DOUBLE *state, int dim, void *context,
             ABMD_DOUBLE *out, ABMD_DOUBLE *rhs_out, ABMD_DOUBLE **memory, int method) {
  if (method == METHOD_DOPRI8) {
    dopri8_step(f, h, t, state, dim, context, out, rhs_out, memory);
  } else if (method == METHOD_RK4) {
    rk4_step(f, h, t, state, dim, context, out, rhs_out, memory);
  }
}
