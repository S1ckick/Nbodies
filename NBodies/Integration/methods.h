#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>

#include "../Utils/helper.h"
#endif

#include <vector>

//#include "abmd/abmd.h"
#include "pointmasses.h"

template <typename Type>
void Euler(std::vector<Type> &x, std::vector<Type> &masses, Type h) {
  ObjectsData<Type> *objects = new ObjectsData<Type>(masses);
  std::vector<Type> k_1;
  k_1.resize(x.size());

  pointmassesCalculateXdot(x, k_1, objects);

  for (int i = 0; i < x.size(); i++) {
    x[i] = x[i] + k_1[i] * h;
  }
}

template <typename Type>
void RungeKutta4(std::vector<Type> &x, std::vector<Type> &masses, Type h) {
  ObjectsData<Type> *objects = new ObjectsData<Type>(masses);
  std::vector<Type> prev_x, new_x, k_1, k_2, k_3, k_4;
  prev_x.resize(x.size());
  new_x.resize(x.size());
  k_1.resize(x.size());
  k_2.resize(x.size());
  k_3.resize(x.size());
  k_4.resize(x.size());

  pointmassesCalculateXdot(x, k_1, objects);

  for (int i = 0; i < x.size(); i++) {
    prev_x[i] = x[i] + k_1[i] * Type(0.5) * h;
  }

  pointmassesCalculateXdot(prev_x, k_2, objects);

  for (int i = 0; i < x.size(); i++) {
    prev_x[i] = x[i] + k_2[i] * Type(0.5) * h;
  }

  pointmassesCalculateXdot(prev_x, k_3, objects);

  for (int i = 0; i < x.size(); i++) {
    prev_x[i] = x[i] + k_3[i] * h;
  }

  pointmassesCalculateXdot(prev_x, k_4, objects);

  for (int i = 0; i < x.size(); i++) {
    x[i] =
        x[i] + (k_1[i] * (Type(1) / Type(6)) + k_2[i] * (Type(1) / Type(3)) +
                k_3[i] * (Type(1) / Type(3)) + k_4[i] * (Type(1) / Type(6))) *
                   h;
  }
}

template <typename Type>
std::vector<Type> initDDCoef() {
  Type b1, b6, b7, b8, b9, b10, b11, b12, a21, a31, a32, a41, a43, a51, a53,
      a54, a61, a64, a65, a71, a74, a75, a76, a81, a84, a85, a86, a87, a91, a94,
      a95, a96, a97, a98, a101;
  Type a104, a105, a106, a107, a108, a109, a111, a114, a115, a116, a117, a118,
      a119, a1110, a121, a124, a125, a126, a127, a128, a129, a1210, a1211;

#ifdef NUMBER_DOUBLE_DOUBLE
  read("5.42937341165687622380535766363E-2", b1);
  read("4.45031289275240888144113950566E0", b6);
  read("1.89151789931450038304281599044E0", b7);
  read("-5.8012039600105847814672114227E0", b8);
  read("3.1116436695781989440891606237E-1", b9);
  read("-1.52160949662516078556178806805E-1", b10);
  read("2.01365400804030348374776537501E-1", b11);
  read("4.47106157277725905176885569043E-2", b12);

  read("5.26001519587677318785587544488E-2", a21);
  read("1.97250569845378994544595329183E-2", a31);
  read("5.91751709536136983633785987549E-2", a32);
  read("2.95875854768068491816892993775E-2", a41);
  read("8.87627564304205475450678981324E-2", a43);
  read("2.41365134159266685502369798665E-1", a51);
  read("-8.84549479328286085344864962717E-1", a53);
  read("9.24834003261792003115737966543E-1", a54);
  read("3.7037037037037037037037037037E-2", a61);
  read("1.70828608729473871279604482173E-1", a64);
  read("1.25467687566822425016691814123E-1", a65);
  read("3.7109375E-2", a71);
  read("1.70252211019544039314978060272E-1", a74);
  read("6.02165389804559606850219397283E-2", a75);
  read("-1.7578125E-2", a76);

  read("3.70920001185047927108779319836E-2", a81);
  read("1.70383925712239993810214054705E-1", a84);
  read("1.07262030446373284651809199168E-1", a85);
  read("-1.53194377486244017527936158236E-2", a86);
  read("8.27378916381402288758473766002E-3", a87);
  read("6.24110958716075717114429577812E-1", a91);
  read("-3.36089262944694129406857109825E0", a94);
  read("-8.68219346841726006818189891453E-1", a95);
  read("2.75920996994467083049415600797E1", a96);
  read("2.01540675504778934086186788979E1", a97);
  read("-4.34898841810699588477366255144E1", a98);
  read("4.77662536438264365890433908527E-1", a101);

  read("-2.48811461997166764192642586468E0", a104);
  read("-5.90290826836842996371446475743E-1", a105);
  read("2.12300514481811942347288949897E1", a106);
  read("1.52792336328824235832596922938E1", a107);
  read("-3.32882109689848629194453265587E1", a108);
  read("-2.03312017085086261358222928593E-2", a109);

  read("-9.3714243008598732571704021658E-1", a111);
  read("5.18637242884406370830023853209E0", a114);
  read("1.09143734899672957818500254654E0", a115);
  read("-8.14978701074692612513997267357E0", a116);
  read("-1.85200656599969598641566180701E1", a117);
  read("2.27394870993505042818970056734E1", a118);
  read("2.49360555267965238987089396762E0", a119);
  read("-3.0467644718982195003823669022E0", a1110);
  read("2.27331014751653820792359768449E0", a121);
  read("-1.05344954667372501984066689879E1", a124);
  read("-2.00087205822486249909675718444E0", a125);
  read("-1.79589318631187989172765950534E1", a126);
  read("2.79488845294199600508499808837E1", a127);
  read("-2.85899827713502369474065508674E0", a128);
  read("-8.87285693353062954433549289258E0", a129);
  read("1.23605671757943030647266201528E1", a1210);
  read("6.43392746015763530355970484046E-1", a1211);
#endif
#ifdef NUMBER_DOUBLE

  b1 = 5.42937341165687622380535766363E-2L;
  b6 = 4.45031289275240888144113950566E0L;
  b7 = 1.89151789931450038304281599044E0L;
  b8 = -5.8012039600105847814672114227E0L;
  b9 = 3.1116436695781989440891606237E-1L;
  b10 = -1.52160949662516078556178806805E-1L;
  b11 = 2.01365400804030348374776537501E-1L;
  b12 = 4.47106157277725905176885569043E-2L;

  a21 = 5.26001519587677318785587544488E-2L;
  a31 = 1.97250569845378994544595329183E-2L;
  a32 = 5.91751709536136983633785987549E-2L;
  a41 = 2.95875854768068491816892993775E-2L;
  a43 = 8.87627564304205475450678981324E-2L;
  a51 = 2.41365134159266685502369798665E-1L;
  a53 = -8.84549479328286085344864962717E-1L;
  a54 = 9.24834003261792003115737966543E-1L;
  a61 = 3.7037037037037037037037037037E-2L;
  a64 = 1.70828608729473871279604482173E-1L;
  a65 = 1.25467687566822425016691814123E-1L;
  a71 = 3.7109375E-2L;
  a74 = 1.70252211019544039314978060272E-1L;
  a75 = 6.02165389804559606850219397283E-2L;
  a76 = -1.7578125E-2L;

  a81 = 3.70920001185047927108779319836E-2L;
  a84 = 1.70383925712239993810214054705E-1L;
  a85 = 1.07262030446373284651809199168E-1L;
  a86 = -1.53194377486244017527936158236E-2L;
  a87 = 8.27378916381402288758473766002E-3L;
  a91 = 6.24110958716075717114429577812E-1L;
  a94 = -3.36089262944694129406857109825E0L;
  a95 = -8.68219346841726006818189891453E-1L;
  a96 = 2.75920996994467083049415600797E1L;
  a97 = 2.01540675504778934086186788979E1L;
  a98 = -4.34898841810699588477366255144E1L;
  a101 = 4.77662536438264365890433908527E-1L;
  a104 = -2.48811461997166764192642586468E0L;
  a105 = -5.90290826836842996371446475743E-1L;
  a106 = 2.12300514481811942347288949897E1L;
  a107 = 1.52792336328824235832596922938E1L;
  a108 = -3.32882109689848629194453265587E1L;
  a109 = -2.03312017085086261358222928593E-2L;

  a111 = -9.3714243008598732571704021658E-1L;
  a114 = 5.18637242884406370830023853209E0L;
  a115 = 1.09143734899672957818500254654E0L;
  a116 = -8.14978701074692612513997267357E0L;
  a117 = -1.85200656599969598641566180701E1L;
  a118 = 2.27394870993505042818970056734E1L;
  a119 = 2.49360555267965238987089396762E0L;
  a1110 = -3.0467644718982195003823669022E0L;
  a121 = 2.27331014751653820792359768449E0L;
  a124 = -1.05344954667372501984066689879E1L;
  a125 = -2.00087205822486249909675718444E0L;
  a126 = -1.79589318631187989172765950534E1L;
  a127 = 2.79488845294199600508499808837E1L;
  a128 = -2.85899827713502369474065508674E0L;
  a129 = -8.87285693353062954433549289258E0L;
  a1210 = 1.23605671757943030647266201528E1L;
  a1211 = 6.43392746015763530355970484046E-1L;

#endif

  std::vector<Type> items;
  items.push_back(a21);
  items.push_back(a31);
  items.push_back(a32);
  items.push_back(a41);
  items.push_back(a43);
  items.push_back(a51);
  items.push_back(a53);
  items.push_back(a54);
  items.push_back(a61);
  items.push_back(a64);
  items.push_back(a65);
  items.push_back(a71);
  items.push_back(a74);
  items.push_back(a75);
  items.push_back(a76);
  items.push_back(a81);
  items.push_back(a84);
  items.push_back(a85);
  items.push_back(a86);
  items.push_back(a87);
  items.push_back(a91);
  items.push_back(a94);
  items.push_back(a95);
  items.push_back(a96);
  items.push_back(a97);
  items.push_back(a98);
  items.push_back(a101);
  items.push_back(a104);
  items.push_back(a105);
  items.push_back(a106);
  items.push_back(a107);
  items.push_back(a108);
  items.push_back(a109);
  items.push_back(a111);
  items.push_back(a114);
  items.push_back(a115);
  items.push_back(a116);
  items.push_back(a117);
  items.push_back(a118);
  items.push_back(a119);
  items.push_back(a1110);
  items.push_back(a121);
  items.push_back(a124);
  items.push_back(a125);
  items.push_back(a126);
  items.push_back(a127);
  items.push_back(a128);
  items.push_back(a129);
  items.push_back(a1210);
  items.push_back(a1211);

  items.push_back(b1);
  items.push_back(b6);
  items.push_back(b7);
  items.push_back(b8);
  items.push_back(b9);
  items.push_back(b10);
  items.push_back(b11);
  items.push_back(b12);

  return items;
}