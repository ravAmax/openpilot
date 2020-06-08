/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1536666512115734907);
void inv_err_fun(double *nom_x, double *true_x, double *out_8250147961189473733);
void H_mod_fun(double *state, double *out_9153956784655755710);
void f_fun(double *state, double dt, double *out_6592157343640020629);
void F_fun(double *state, double dt, double *out_5617018183117879922);
void h_3(double *state, double *unused, double *out_1645097081930375361);
void H_3(double *state, double *unused, double *out_7087457282334079905);
void h_4(double *state, double *unused, double *out_1912949227850046409);
void H_4(double *state, double *unused, double *out_1262420423677199833);
void h_9(double *state, double *unused, double *out_6377980340510303867);
void H_9(double *state, double *unused, double *out_1966925501450689608);
void h_10(double *state, double *unused, double *out_8013204234311813899);
void H_10(double *state, double *unused, double *out_6776713216419999281);
void h_12(double *state, double *unused, double *out_7160788152741406783);
void H_12(double *state, double *unused, double *out_6954421532273645384);
void h_31(double *state, double *unused, double *out_7668855218514234924);
void H_31(double *state, double *unused, double *out_6849722377292652688);
void h_32(double *state, double *unused, double *out_807070843776782916);
void H_32(double *state, double *unused, double *out_3840592621874159454);
void h_13(double *state, double *unused, double *out_5327367460280228198);
void H_13(double *state, double *unused, double *out_8437814592451181722);
void h_14(double *state, double *unused, double *out_6377980340510303867);
void H_14(double *state, double *unused, double *out_1966925501450689608);
void h_19(double *state, double *unused, double *out_131762747836752563);
void H_19(double *state, double *unused, double *out_2162886027585538044);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);