
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2048977256418502276) {
   out_2048977256418502276[0] = delta_x[0] + nom_x[0];
   out_2048977256418502276[1] = delta_x[1] + nom_x[1];
   out_2048977256418502276[2] = delta_x[2] + nom_x[2];
   out_2048977256418502276[3] = delta_x[3] + nom_x[3];
   out_2048977256418502276[4] = delta_x[4] + nom_x[4];
   out_2048977256418502276[5] = delta_x[5] + nom_x[5];
   out_2048977256418502276[6] = delta_x[6] + nom_x[6];
   out_2048977256418502276[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6705258261903374480) {
   out_6705258261903374480[0] = -nom_x[0] + true_x[0];
   out_6705258261903374480[1] = -nom_x[1] + true_x[1];
   out_6705258261903374480[2] = -nom_x[2] + true_x[2];
   out_6705258261903374480[3] = -nom_x[3] + true_x[3];
   out_6705258261903374480[4] = -nom_x[4] + true_x[4];
   out_6705258261903374480[5] = -nom_x[5] + true_x[5];
   out_6705258261903374480[6] = -nom_x[6] + true_x[6];
   out_6705258261903374480[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3692715011050714555) {
   out_3692715011050714555[0] = 1.0;
   out_3692715011050714555[1] = 0.0;
   out_3692715011050714555[2] = 0.0;
   out_3692715011050714555[3] = 0.0;
   out_3692715011050714555[4] = 0.0;
   out_3692715011050714555[5] = 0.0;
   out_3692715011050714555[6] = 0.0;
   out_3692715011050714555[7] = 0.0;
   out_3692715011050714555[8] = 0.0;
   out_3692715011050714555[9] = 1.0;
   out_3692715011050714555[10] = 0.0;
   out_3692715011050714555[11] = 0.0;
   out_3692715011050714555[12] = 0.0;
   out_3692715011050714555[13] = 0.0;
   out_3692715011050714555[14] = 0.0;
   out_3692715011050714555[15] = 0.0;
   out_3692715011050714555[16] = 0.0;
   out_3692715011050714555[17] = 0.0;
   out_3692715011050714555[18] = 1.0;
   out_3692715011050714555[19] = 0.0;
   out_3692715011050714555[20] = 0.0;
   out_3692715011050714555[21] = 0.0;
   out_3692715011050714555[22] = 0.0;
   out_3692715011050714555[23] = 0.0;
   out_3692715011050714555[24] = 0.0;
   out_3692715011050714555[25] = 0.0;
   out_3692715011050714555[26] = 0.0;
   out_3692715011050714555[27] = 1.0;
   out_3692715011050714555[28] = 0.0;
   out_3692715011050714555[29] = 0.0;
   out_3692715011050714555[30] = 0.0;
   out_3692715011050714555[31] = 0.0;
   out_3692715011050714555[32] = 0.0;
   out_3692715011050714555[33] = 0.0;
   out_3692715011050714555[34] = 0.0;
   out_3692715011050714555[35] = 0.0;
   out_3692715011050714555[36] = 1.0;
   out_3692715011050714555[37] = 0.0;
   out_3692715011050714555[38] = 0.0;
   out_3692715011050714555[39] = 0.0;
   out_3692715011050714555[40] = 0.0;
   out_3692715011050714555[41] = 0.0;
   out_3692715011050714555[42] = 0.0;
   out_3692715011050714555[43] = 0.0;
   out_3692715011050714555[44] = 0.0;
   out_3692715011050714555[45] = 1.0;
   out_3692715011050714555[46] = 0.0;
   out_3692715011050714555[47] = 0.0;
   out_3692715011050714555[48] = 0.0;
   out_3692715011050714555[49] = 0.0;
   out_3692715011050714555[50] = 0.0;
   out_3692715011050714555[51] = 0.0;
   out_3692715011050714555[52] = 0.0;
   out_3692715011050714555[53] = 0.0;
   out_3692715011050714555[54] = 1.0;
   out_3692715011050714555[55] = 0.0;
   out_3692715011050714555[56] = 0.0;
   out_3692715011050714555[57] = 0.0;
   out_3692715011050714555[58] = 0.0;
   out_3692715011050714555[59] = 0.0;
   out_3692715011050714555[60] = 0.0;
   out_3692715011050714555[61] = 0.0;
   out_3692715011050714555[62] = 0.0;
   out_3692715011050714555[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2001184518181836126) {
   out_2001184518181836126[0] = state[0];
   out_2001184518181836126[1] = state[1];
   out_2001184518181836126[2] = state[2];
   out_2001184518181836126[3] = state[3];
   out_2001184518181836126[4] = state[4];
   out_2001184518181836126[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2001184518181836126[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2001184518181836126[7] = state[7];
}
void F_fun(double *state, double dt, double *out_1137646979211532694) {
   out_1137646979211532694[0] = 1;
   out_1137646979211532694[1] = 0;
   out_1137646979211532694[2] = 0;
   out_1137646979211532694[3] = 0;
   out_1137646979211532694[4] = 0;
   out_1137646979211532694[5] = 0;
   out_1137646979211532694[6] = 0;
   out_1137646979211532694[7] = 0;
   out_1137646979211532694[8] = 0;
   out_1137646979211532694[9] = 1;
   out_1137646979211532694[10] = 0;
   out_1137646979211532694[11] = 0;
   out_1137646979211532694[12] = 0;
   out_1137646979211532694[13] = 0;
   out_1137646979211532694[14] = 0;
   out_1137646979211532694[15] = 0;
   out_1137646979211532694[16] = 0;
   out_1137646979211532694[17] = 0;
   out_1137646979211532694[18] = 1;
   out_1137646979211532694[19] = 0;
   out_1137646979211532694[20] = 0;
   out_1137646979211532694[21] = 0;
   out_1137646979211532694[22] = 0;
   out_1137646979211532694[23] = 0;
   out_1137646979211532694[24] = 0;
   out_1137646979211532694[25] = 0;
   out_1137646979211532694[26] = 0;
   out_1137646979211532694[27] = 1;
   out_1137646979211532694[28] = 0;
   out_1137646979211532694[29] = 0;
   out_1137646979211532694[30] = 0;
   out_1137646979211532694[31] = 0;
   out_1137646979211532694[32] = 0;
   out_1137646979211532694[33] = 0;
   out_1137646979211532694[34] = 0;
   out_1137646979211532694[35] = 0;
   out_1137646979211532694[36] = 1;
   out_1137646979211532694[37] = 0;
   out_1137646979211532694[38] = 0;
   out_1137646979211532694[39] = 0;
   out_1137646979211532694[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1137646979211532694[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1137646979211532694[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1137646979211532694[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1137646979211532694[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1137646979211532694[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1137646979211532694[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1137646979211532694[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1137646979211532694[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1137646979211532694[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1137646979211532694[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1137646979211532694[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1137646979211532694[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1137646979211532694[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1137646979211532694[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1137646979211532694[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1137646979211532694[56] = 0;
   out_1137646979211532694[57] = 0;
   out_1137646979211532694[58] = 0;
   out_1137646979211532694[59] = 0;
   out_1137646979211532694[60] = 0;
   out_1137646979211532694[61] = 0;
   out_1137646979211532694[62] = 0;
   out_1137646979211532694[63] = 1;
}
void h_25(double *state, double *unused, double *out_3182774736060427257) {
   out_3182774736060427257[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1031646437293131827) {
   out_1031646437293131827[0] = 0;
   out_1031646437293131827[1] = 0;
   out_1031646437293131827[2] = 0;
   out_1031646437293131827[3] = 0;
   out_1031646437293131827[4] = 0;
   out_1031646437293131827[5] = 0;
   out_1031646437293131827[6] = 1;
   out_1031646437293131827[7] = 0;
}
void h_24(double *state, double *unused, double *out_2461486246405790661) {
   out_2461486246405790661[0] = state[4];
   out_2461486246405790661[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8386198841853398353) {
   out_8386198841853398353[0] = 0;
   out_8386198841853398353[1] = 0;
   out_8386198841853398353[2] = 0;
   out_8386198841853398353[3] = 0;
   out_8386198841853398353[4] = 1;
   out_8386198841853398353[5] = 0;
   out_8386198841853398353[6] = 0;
   out_8386198841853398353[7] = 0;
   out_8386198841853398353[8] = 0;
   out_8386198841853398353[9] = 0;
   out_8386198841853398353[10] = 0;
   out_8386198841853398353[11] = 0;
   out_8386198841853398353[12] = 0;
   out_8386198841853398353[13] = 1;
   out_8386198841853398353[14] = 0;
   out_8386198841853398353[15] = 0;
}
void h_30(double *state, double *unused, double *out_6777121964931226171) {
   out_6777121964931226171[0] = state[4];
}
void H_30(double *state, double *unused, double *out_176468745696846407) {
   out_176468745696846407[0] = 0;
   out_176468745696846407[1] = 0;
   out_176468745696846407[2] = 0;
   out_176468745696846407[3] = 0;
   out_176468745696846407[4] = 1;
   out_176468745696846407[5] = 0;
   out_176468745696846407[6] = 0;
   out_176468745696846407[7] = 0;
}
void h_26(double *state, double *unused, double *out_7180283939846045230) {
   out_7180283939846045230[0] = state[7];
}
void H_26(double *state, double *unused, double *out_82437915401707589) {
   out_82437915401707589[0] = 0;
   out_82437915401707589[1] = 0;
   out_82437915401707589[2] = 0;
   out_82437915401707589[3] = 0;
   out_82437915401707589[4] = 0;
   out_82437915401707589[5] = 0;
   out_82437915401707589[6] = 0;
   out_82437915401707589[7] = 1;
}
void h_27(double *state, double *unused, double *out_9040611449494797098) {
   out_9040611449494797098[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4913572462411010281) {
   out_4913572462411010281[0] = 0;
   out_4913572462411010281[1] = 0;
   out_4913572462411010281[2] = 0;
   out_4913572462411010281[3] = 1;
   out_4913572462411010281[4] = 0;
   out_4913572462411010281[5] = 0;
   out_4913572462411010281[6] = 0;
   out_4913572462411010281[7] = 0;
}
void h_29(double *state, double *unused, double *out_4309499505199190808) {
   out_4309499505199190808[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8414405033380959203) {
   out_8414405033380959203[0] = 0;
   out_8414405033380959203[1] = 1;
   out_8414405033380959203[2] = 0;
   out_8414405033380959203[3] = 0;
   out_8414405033380959203[4] = 0;
   out_8414405033380959203[5] = 0;
   out_8414405033380959203[6] = 0;
   out_8414405033380959203[7] = 0;
}
void h_28(double *state, double *unused, double *out_3744830019118191349) {
   out_3744830019118191349[0] = state[5];
   out_3744830019118191349[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1058322534344106699) {
   out_1058322534344106699[0] = 0;
   out_1058322534344106699[1] = 0;
   out_1058322534344106699[2] = 0;
   out_1058322534344106699[3] = 0;
   out_1058322534344106699[4] = 0;
   out_1058322534344106699[5] = 1;
   out_1058322534344106699[6] = 0;
   out_1058322534344106699[7] = 0;
   out_1058322534344106699[8] = 0;
   out_1058322534344106699[9] = 0;
   out_1058322534344106699[10] = 0;
   out_1058322534344106699[11] = 0;
   out_1058322534344106699[12] = 0;
   out_1058322534344106699[13] = 0;
   out_1058322534344106699[14] = 1;
   out_1058322534344106699[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
