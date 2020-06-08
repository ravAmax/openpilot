/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2048977256418502276);
void inv_err_fun(double *nom_x, double *true_x, double *out_6705258261903374480);
void H_mod_fun(double *state, double *out_3692715011050714555);
void f_fun(double *state, double dt, double *out_2001184518181836126);
void F_fun(double *state, double dt, double *out_1137646979211532694);
void h_25(double *state, double *unused, double *out_3182774736060427257);
void H_25(double *state, double *unused, double *out_1031646437293131827);
void h_24(double *state, double *unused, double *out_2461486246405790661);
void H_24(double *state, double *unused, double *out_8386198841853398353);
void h_30(double *state, double *unused, double *out_6777121964931226171);
void H_30(double *state, double *unused, double *out_176468745696846407);
void h_26(double *state, double *unused, double *out_7180283939846045230);
void H_26(double *state, double *unused, double *out_82437915401707589);
void h_27(double *state, double *unused, double *out_9040611449494797098);
void H_27(double *state, double *unused, double *out_4913572462411010281);
void h_29(double *state, double *unused, double *out_4309499505199190808);
void H_29(double *state, double *unused, double *out_8414405033380959203);
void h_28(double *state, double *unused, double *out_3744830019118191349);
void H_28(double *state, double *unused, double *out_1058322534344106699);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
