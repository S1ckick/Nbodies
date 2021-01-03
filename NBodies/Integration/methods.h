//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_METHODS_H
#define NBODIES_METHODS_H

#include "../Nbodies/nbodies.h"

template <typename Type>
void RungeKutta4(std::vector<Body<Type>> &bodies, double h) {

    std::vector<Body<Type>> k_1;
    copyBodies(bodies, k_1);

    std::vector<Body<Type>> temp;
    copyBodies(bodies, temp);


    f(temp,k_1);

    //temp = bodies + 0.5*k1*h
    copyBodies(bodies + k_1*(Type(0.5)*h), temp);

    std::vector<Body<Type>> k_2;
    copyBodies(bodies, k_2);
    f(temp, k_2);

    //temp = bodies + 0.5*k2*h
    copyBodies(bodies + k_2*(h*Type(0.5)), temp);

    std::vector<Body<Type>> k_3;
    copyBodies(bodies, k_3);
    f(temp, k_3);

    //temp = bodies + 1.0*k3*h
    copyBodies(bodies + k_3*h, temp);

    std::vector<Body<Type>> k_4;
    copyBodies(bodies, k_4);
    f(temp, k_4);

    copyBodies(k_1*(Type(1)/Type(6)), k_1);
    copyBodies(k_2*(Type(1)/Type(3)), k_2);
    copyBodies(k_3*(Type(1)/Type(3)), k_3);
    copyBodies(k_4*(Type(1)/Type(6)), k_4);

    //y += 	dt( k_1/6 + k_2/3 + k_3/3 + k_4/6 )
    copyBodies( bodies + (k_1 + k_2 + k_3 + k_4)*h ,bodies);
}




template <typename Type>
void dormanPrince8(std::vector<Body<Type>> &bodies, double h) {

    Type b1 =   5.42937341165687622380535766363E-2L,
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

    std::vector<Body<Type>> k_1;
    copyBodies(bodies, k_1);

    std::vector<Body<Type>> temp;
    copyBodies(bodies, temp);


    f(temp,k_1);

    //temp = bodies + a21*k1*h
    copyBodies(bodies + k_1*(h*a21), temp);

    std::vector<Body<Type>> k_2;
    copyBodies(bodies, k_2);
    f(temp, k_2);

    //temp = bodies + h*(a31*k1 + a32*k2)
    copyBodies(bodies + (k_1*a31 + k_2*a32)*h, temp);

    std::vector<Body<Type>> k_3;
    copyBodies(bodies, k_3);
    f(temp, k_3);

    //temp = bodies + h*( a41*k1 + a43*k3 )
    copyBodies(bodies + (k_1*a41 + k_3*a43)*h, temp);

    std::vector<Body<Type>> k_4;
    copyBodies(bodies, k_4);
    f(temp, k_4);

    //temp = bodies + h*( a51*k1 + a53*k3 + a54*k4)
    copyBodies(bodies + (k_1*a51 + k_3*a53 + k_4*a54)*h,temp);

    std::vector<Body<Type>> k_5;
    copyBodies(bodies, k_5);
    f(temp, k_5);

    //temp = bodies + h*( a61*k1 + a64*k4 + a65*k5)
    copyBodies(bodies + (k_1*a61 + k_4*a64 + k_5*a65)*h,temp);

    std::vector<Body<Type>> k_6;
    copyBodies(bodies, k_6);
    f(temp, k_6);

    //temp = bodies + h*( a71*k1 + a74*k4 + a75*k5 + a76*k6)
    copyBodies(bodies + (k_1*a71 + k_4*a74 + k_5*a75 + k_6*a76)*h,temp);

    std::vector<Body<Type>> k_7;
    copyBodies(bodies, k_7);
    f(temp, k_7);

    //temp = bodies + h*( a81*k1 + a84*k4 + a85*k5 + a86*k6 + a87*k7)
    copyBodies(bodies + (k_1*a81 + k_4*a84 + k_5*a85 + k_6*a86 + k_7*a87)*h,temp);

    std::vector<Body<Type>> k_8;
    copyBodies(bodies, k_8);
    f(temp, k_8);

    //temp = bodies + h*( a91*k1 + a94*k4 + a95*k5 + a96*k6 + a97*k7 + a98*k8)
    copyBodies(bodies + (k_1*a91 + k_4*a94 + k_5*a95 + k_6*a96 + k_7*a97 + k_8*a98)*h,temp);

    std::vector<Body<Type>> k_9;
    copyBodies(bodies, k_9);
    f(temp, k_9);

    //temp = bodies + h*( a101*k1 + a104*k4 + a105*k5 + a106*k6 + a107*k7 + a108*k8 + a109*k9)
    copyBodies(bodies + (k_1*a101 + k_4*a104 + k_5*a105 + k_6*a106 + k_7*a107 + k_8*a108 + k_9*a109)*h,temp);

    std::vector<Body<Type>> k_10;
    copyBodies(bodies, k_10);
    f(temp, k_10);

    //temp = bodies + h*( a111*k1 + a114*k4 + a115*k5 + a116*k6 + a117*k7 + a118*k8 + a119*k9 + a1110*k_10)
    copyBodies(bodies + (k_1*a111 + k_4*a114 + k_5*a115 + k_6*a116 + k_7*a117 + k_8*a118 + k_9*a119 + k_10*a1110)*h,temp);

    copyBodies(bodies, k_2);
    f(temp, k_2);

    //temp = bodies + h*( a121*k1 + a124*k4 + a125*k5 + a126*k6 + a127*k7 + a128*k8 + a129*k9 + a1210*k_10 + a1211*k_2)
    copyBodies(bodies + (k_1*a121 + k_4*a124 + k_5*a125 + k_6*a126 + k_7*a127 + k_8*a128 + k_9*a129 + k_10*a1210 + k_2*a1211)*h,temp);

    copyBodies(bodies, k_3);
    f(temp, k_3);

    copyBodies(k_1*b1 + k_6*b6 + k_7*b7 + k_8*b8 + k_9*b9 + k_10*b10 + k_2*b11 + k_3*b12,k_4);
    copyBodies(bodies + k_4*h,bodies);

}




#endif //NBODIES_METHODS_H
