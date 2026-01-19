//
/*
Parameter guide for running simulations under different physical conditions:

The simulation conditions can be modified by changing the following parameters:
- Lattice constant
- Optical field strength
- Pump frequency
- Pulse length
- Hubbard U
- Hopping t
- Electron–environment scattering rate

All parameters listed above must be modified consistently in both
`double _Complex *hhg(double alpha)` and
`double _Complex *test_f(...)`,
with identical values used in both functions.

The electron-environment scattering rate is an exception and should be modified only in `double _Complex *test_f(...)`.
This parameter is defined around line 3106 of the source code.
All physical quantities are expressed in atomic units unless otherwise noted.

*/

//
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include <lapacke.h>

# include <gsl/gsl_blas.h>
# include <gsl/gsl_complex_math.h>

# include "rk4_kxky.h"

double _Complex *test_f (double t, int m, double _Complex u[], double tau, double kx, double ky, double alpha );

double _Complex *hhg(double alpha){

    int m = 24010000, j = 0, q = 0, rx, ry;// m = hamiltonian size
    int totn, dian, NCn = 70, sn, hopi;
    int rhi, rhj;
    double t0 = -2067.1*2.0, t1 , dt = 1.4, tau0 = 1.0;// solve ODE, pulse start t0
    double kx=0.0, ky=0.0;

    double latt_const = 7.56;// lattice constant (a.u.)

    double _Complex *sol;// rk4 solution
    double _Complex *ini;// initial condition
    double _Complex *Jave;// main data
    double pi = 3.14159265358979323846;
    double Ax;
    double ellip = 0.0;
    double Ex = ( 3.5/6.5 )*0.00126410*alpha/(sqrt(1.0+ellip*ellip));//optical field strength
    double wpump = 0.26/27.2114; //pump frequency
    double Tfull = 4134.2*2.0; //pulse length

    double u_rep = 0.0*0.52/27.2114;// Hubbard U
    double t_hub = 0.52/27.2114;// hopping t

    lapack_int n = 4900;
    lapack_int lda = n;
    lapack_int info;
    double E_eigval[n];
    lapack_complex_double H_0 [n*n];

    int size = 4900;
    gsl_matrix_complex *U = NULL;
    gsl_matrix_complex *RHO_TILDE = NULL;
    gsl_matrix_complex *MIDDLE = NULL;
    gsl_matrix_complex *RHO = NULL;

    gsl_matrix_complex *J_X = NULL;
    gsl_matrix_complex *J_Y = NULL;
    gsl_matrix_complex *CURRENT_X = NULL;
    gsl_matrix_complex *CURRENT_Y = NULL;

    U = gsl_matrix_complex_alloc(size, size);
    RHO_TILDE = gsl_matrix_complex_alloc(size, size);
    MIDDLE = gsl_matrix_complex_alloc(size,size);
    RHO = gsl_matrix_complex_alloc(size, size);

    J_X = gsl_matrix_complex_alloc(size, size);
    J_Y = gsl_matrix_complex_alloc(size, size);
    CURRENT_X = gsl_matrix_complex_alloc(size,size);
    CURRENT_Y = gsl_matrix_complex_alloc(size, size);

    gsl_matrix_complex_set_all(CURRENT_X, GSL_COMPLEX_ZERO);
    gsl_matrix_complex_set_all(CURRENT_Y, GSL_COMPLEX_ZERO);

// 

    double sing_s1[NCn], sing_s2[NCn], sing_s3[NCn], sing_s4[NCn], sing_s5[NCn], sing_s6[NCn], sing_s7[NCn], sing_s8[NCn];

    sing_s1[0] = 1.0,  sing_s2[0] = 1.0,  sing_s3[0] = 1.0,  sing_s4[0] = 1.0,  sing_s5[0] = 0.0,  sing_s6[0] = 0.0,  sing_s7[0] = 0.0,  sing_s8[0] = 0.0;
    sing_s1[1] = 1.0,  sing_s2[1] = 1.0,  sing_s3[1] = 1.0,  sing_s4[1] = 0.0,  sing_s5[1] = 1.0,  sing_s6[1] = 0.0,  sing_s7[1] = 0.0,  sing_s8[1] = 0.0;
    sing_s1[2] = 1.0,  sing_s2[2] = 1.0,  sing_s3[2] = 0.0,  sing_s4[2] = 1.0,  sing_s5[2] = 1.0,  sing_s6[2] = 0.0,  sing_s7[2] = 0.0,  sing_s8[2] = 0.0;
    sing_s1[3] = 1.0,  sing_s2[3] = 0.0,  sing_s3[3] = 1.0,  sing_s4[3] = 1.0,  sing_s5[3] = 1.0,  sing_s6[3] = 0.0,  sing_s7[3] = 0.0,  sing_s8[3] = 0.0;
    sing_s1[4] = 0.0,  sing_s2[4] = 1.0,  sing_s3[4] = 1.0,  sing_s4[4] = 1.0,  sing_s5[4] = 1.0,  sing_s6[4] = 0.0,  sing_s7[4] = 0.0,  sing_s8[4] = 0.0;
    sing_s1[5] = 1.0,  sing_s2[5] = 1.0,  sing_s3[5] = 1.0,  sing_s4[5] = 0.0,  sing_s5[5] = 0.0,  sing_s6[5] = 1.0,  sing_s7[5] = 0.0,  sing_s8[5] = 0.0;
    sing_s1[6] = 1.0,  sing_s2[6] = 1.0,  sing_s3[6] = 0.0,  sing_s4[6] = 1.0,  sing_s5[6] = 0.0,  sing_s6[6] = 1.0,  sing_s7[6] = 0.0,  sing_s8[6] = 0.0;
    sing_s1[7] = 1.0,  sing_s2[7] = 0.0,  sing_s3[7] = 1.0,  sing_s4[7] = 1.0,  sing_s5[7] = 0.0,  sing_s6[7] = 1.0,  sing_s7[7] = 0.0,  sing_s8[7] = 0.0;
    sing_s1[8] = 0.0,  sing_s2[8] = 1.0,  sing_s3[8] = 1.0,  sing_s4[8] = 1.0,  sing_s5[8] = 0.0,  sing_s6[8] = 1.0,  sing_s7[8] = 0.0,  sing_s8[8] = 0.0;
    sing_s1[9] = 1.0,  sing_s2[9] = 1.0,  sing_s3[9] = 0.0,  sing_s4[9] = 0.0,  sing_s5[9] = 1.0,  sing_s6[9] = 1.0,  sing_s7[9] = 0.0,  sing_s8[9] = 0.0;
    sing_s1[10] = 1.0, sing_s2[10] = 0.0, sing_s3[10] = 1.0, sing_s4[10] = 0.0, sing_s5[10] = 1.0, sing_s6[10] = 1.0, sing_s7[10] = 0.0, sing_s8[10] = 0.0;
    sing_s1[11] = 0.0, sing_s2[11] = 1.0, sing_s3[11] = 1.0, sing_s4[11] = 0.0, sing_s5[11] = 1.0, sing_s6[11] = 1.0, sing_s7[11] = 0.0, sing_s8[11] = 0.0;
    sing_s1[12] = 1.0, sing_s2[12] = 0.0, sing_s3[12] = 0.0, sing_s4[12] = 1.0, sing_s5[12] = 1.0, sing_s6[12] = 1.0, sing_s7[12] = 0.0, sing_s8[12] = 0.0;
    sing_s1[13] = 0.0, sing_s2[13] = 1.0, sing_s3[13] = 0.0, sing_s4[13] = 1.0, sing_s5[13] = 1.0, sing_s6[13] = 1.0, sing_s7[13] = 0.0, sing_s8[13] = 0.0;
    sing_s1[14] = 0.0, sing_s2[14] = 0.0, sing_s3[14] = 1.0, sing_s4[14] = 1.0, sing_s5[14] = 1.0, sing_s6[14] = 1.0, sing_s7[14] = 0.0, sing_s8[14] = 0.0;
    sing_s1[15] = 1.0, sing_s2[15] = 1.0, sing_s3[15] = 1.0, sing_s4[15] = 0.0, sing_s5[15] = 0.0, sing_s6[15] = 0.0, sing_s7[15] = 1.0, sing_s8[15] = 0.0;
    sing_s1[16] = 1.0, sing_s2[16] = 1.0, sing_s3[16] = 0.0, sing_s4[16] = 1.0, sing_s5[16] = 0.0, sing_s6[16] = 0.0, sing_s7[16] = 1.0, sing_s8[16] = 0.0;
    sing_s1[17] = 1.0, sing_s2[17] = 0.0, sing_s3[17] = 1.0, sing_s4[17] = 1.0, sing_s5[17] = 0.0, sing_s6[17] = 0.0, sing_s7[17] = 1.0, sing_s8[17] = 0.0;
    sing_s1[18] = 0.0, sing_s2[18] = 1.0, sing_s3[18] = 1.0, sing_s4[18] = 1.0, sing_s5[18] = 0.0, sing_s6[18] = 0.0, sing_s7[18] = 1.0, sing_s8[18] = 0.0;
    sing_s1[19] = 1.0, sing_s2[19] = 1.0, sing_s3[19] = 0.0, sing_s4[19] = 0.0, sing_s5[19] = 1.0, sing_s6[19] = 0.0, sing_s7[19] = 1.0, sing_s8[19] = 0.0;
    sing_s1[20] = 1.0, sing_s2[20] = 0.0, sing_s3[20] = 1.0, sing_s4[20] = 0.0, sing_s5[20] = 1.0, sing_s6[20] = 0.0, sing_s7[20] = 1.0, sing_s8[20] = 0.0;
    sing_s1[21] = 0.0, sing_s2[21] = 1.0, sing_s3[21] = 1.0, sing_s4[21] = 0.0, sing_s5[21] = 1.0, sing_s6[21] = 0.0, sing_s7[21] = 1.0, sing_s8[21] = 0.0;
    sing_s1[22] = 1.0, sing_s2[22] = 0.0, sing_s3[22] = 0.0, sing_s4[22] = 1.0, sing_s5[22] = 1.0, sing_s6[22] = 0.0, sing_s7[22] = 1.0, sing_s8[22] = 0.0;
    sing_s1[23] = 0.0, sing_s2[23] = 1.0, sing_s3[23] = 0.0, sing_s4[23] = 1.0, sing_s5[23] = 1.0, sing_s6[23] = 0.0, sing_s7[23] = 1.0, sing_s8[23] = 0.0;
    sing_s1[24] = 0.0, sing_s2[24] = 0.0, sing_s3[24] = 1.0, sing_s4[24] = 1.0, sing_s5[24] = 1.0, sing_s6[24] = 0.0, sing_s7[24] = 1.0, sing_s8[24] = 0.0;
    sing_s1[25] = 1.0, sing_s2[25] = 1.0, sing_s3[25] = 0.0, sing_s4[25] = 0.0, sing_s5[25] = 0.0, sing_s6[25] = 1.0, sing_s7[25] = 1.0, sing_s8[25] = 0.0;
    sing_s1[26] = 1.0, sing_s2[26] = 0.0, sing_s3[26] = 1.0, sing_s4[26] = 0.0, sing_s5[26] = 0.0, sing_s6[26] = 1.0, sing_s7[26] = 1.0, sing_s8[26] = 0.0;
    sing_s1[27] = 0.0, sing_s2[27] = 1.0, sing_s3[27] = 1.0, sing_s4[27] = 0.0, sing_s5[27] = 0.0, sing_s6[27] = 1.0, sing_s7[27] = 1.0, sing_s8[27] = 0.0;
    sing_s1[28] = 1.0, sing_s2[28] = 0.0, sing_s3[28] = 0.0, sing_s4[28] = 1.0, sing_s5[28] = 0.0, sing_s6[28] = 1.0, sing_s7[28] = 1.0, sing_s8[28] = 0.0;
    sing_s1[29] = 0.0, sing_s2[29] = 1.0, sing_s3[29] = 0.0, sing_s4[29] = 1.0, sing_s5[29] = 0.0, sing_s6[29] = 1.0, sing_s7[29] = 1.0, sing_s8[29] = 0.0;
    sing_s1[30] = 0.0, sing_s2[30] = 0.0, sing_s3[30] = 1.0, sing_s4[30] = 1.0, sing_s5[30] = 0.0, sing_s6[30] = 1.0, sing_s7[30] = 1.0, sing_s8[30] = 0.0;
    sing_s1[31] = 1.0, sing_s2[31] = 0.0, sing_s3[31] = 0.0, sing_s4[31] = 0.0, sing_s5[31] = 1.0, sing_s6[31] = 1.0, sing_s7[31] = 1.0, sing_s8[31] = 0.0;
    sing_s1[32] = 0.0, sing_s2[32] = 1.0, sing_s3[32] = 0.0, sing_s4[32] = 0.0, sing_s5[32] = 1.0, sing_s6[32] = 1.0, sing_s7[32] = 1.0, sing_s8[32] = 0.0;
    sing_s1[33] = 0.0, sing_s2[33] = 0.0, sing_s3[33] = 1.0, sing_s4[33] = 0.0, sing_s5[33] = 1.0, sing_s6[33] = 1.0, sing_s7[33] = 1.0, sing_s8[33] = 0.0;
    sing_s1[34] = 0.0, sing_s2[34] = 0.0, sing_s3[34] = 0.0, sing_s4[34] = 1.0, sing_s5[34] = 1.0, sing_s6[34] = 1.0, sing_s7[34] = 1.0, sing_s8[34] = 0.0;

    sing_s1[35] = 1.0, sing_s2[35] = 1.0, sing_s3[35] = 1.0, sing_s4[35] = 0.0, sing_s5[35] = 0.0, sing_s6[35] = 0.0, sing_s7[35] = 0.0, sing_s8[35] = 1.0;
    sing_s1[36] = 1.0, sing_s2[36] = 1.0, sing_s3[36] = 0.0, sing_s4[36] = 1.0, sing_s5[36] = 0.0, sing_s6[36] = 0.0, sing_s7[36] = 0.0, sing_s8[36] = 1.0;
    sing_s1[37] = 1.0, sing_s2[37] = 0.0, sing_s3[37] = 1.0, sing_s4[37] = 1.0, sing_s5[37] = 0.0, sing_s6[37] = 0.0, sing_s7[37] = 0.0, sing_s8[37] = 1.0;
    sing_s1[38] = 0.0, sing_s2[38] = 1.0, sing_s3[38] = 1.0, sing_s4[38] = 1.0, sing_s5[38] = 0.0, sing_s6[38] = 0.0, sing_s7[38] = 0.0, sing_s8[38] = 1.0;
    sing_s1[39] = 1.0, sing_s2[39] = 1.0, sing_s3[39] = 0.0, sing_s4[39] = 0.0, sing_s5[39] = 1.0, sing_s6[39] = 0.0, sing_s7[39] = 0.0, sing_s8[39] = 1.0;
    sing_s1[40] = 1.0, sing_s2[40] = 0.0, sing_s3[40] = 1.0, sing_s4[40] = 0.0, sing_s5[40] = 1.0, sing_s6[40] = 0.0, sing_s7[40] = 0.0, sing_s8[40] = 1.0;
    sing_s1[41] = 0.0, sing_s2[41] = 1.0, sing_s3[41] = 1.0, sing_s4[41] = 0.0, sing_s5[41] = 1.0, sing_s6[41] = 0.0, sing_s7[41] = 0.0, sing_s8[41] = 1.0;
    sing_s1[42] = 1.0, sing_s2[42] = 0.0, sing_s3[42] = 0.0, sing_s4[42] = 1.0, sing_s5[42] = 1.0, sing_s6[42] = 0.0, sing_s7[42] = 0.0, sing_s8[42] = 1.0;
    sing_s1[43] = 0.0, sing_s2[43] = 1.0, sing_s3[43] = 0.0, sing_s4[43] = 1.0, sing_s5[43] = 1.0, sing_s6[43] = 0.0, sing_s7[43] = 0.0, sing_s8[43] = 1.0;
    sing_s1[44] = 0.0, sing_s2[44] = 0.0, sing_s3[44] = 1.0, sing_s4[44] = 1.0, sing_s5[44] = 1.0, sing_s6[44] = 0.0, sing_s7[44] = 0.0, sing_s8[44] = 1.0;
    sing_s1[45] = 1.0, sing_s2[45] = 1.0, sing_s3[45] = 0.0, sing_s4[45] = 0.0, sing_s5[45] = 0.0, sing_s6[45] = 1.0, sing_s7[45] = 0.0, sing_s8[45] = 1.0;
    sing_s1[46] = 1.0, sing_s2[46] = 0.0, sing_s3[46] = 1.0, sing_s4[46] = 0.0, sing_s5[46] = 0.0, sing_s6[46] = 1.0, sing_s7[46] = 0.0, sing_s8[46] = 1.0;
    sing_s1[47] = 0.0, sing_s2[47] = 1.0, sing_s3[47] = 1.0, sing_s4[47] = 0.0, sing_s5[47] = 0.0, sing_s6[47] = 1.0, sing_s7[47] = 0.0, sing_s8[47] = 1.0;
    sing_s1[48] = 1.0, sing_s2[48] = 0.0, sing_s3[48] = 0.0, sing_s4[48] = 1.0, sing_s5[48] = 0.0, sing_s6[48] = 1.0, sing_s7[48] = 0.0, sing_s8[48] = 1.0;
    sing_s1[49] = 0.0, sing_s2[49] = 1.0, sing_s3[49] = 0.0, sing_s4[49] = 1.0, sing_s5[49] = 0.0, sing_s6[49] = 1.0, sing_s7[49] = 0.0, sing_s8[49] = 1.0;
    sing_s1[50] = 0.0, sing_s2[50] = 0.0, sing_s3[50] = 1.0, sing_s4[50] = 1.0, sing_s5[50] = 0.0, sing_s6[50] = 1.0, sing_s7[50] = 0.0, sing_s8[50] = 1.0;
    sing_s1[51] = 1.0, sing_s2[51] = 0.0, sing_s3[51] = 0.0, sing_s4[51] = 0.0, sing_s5[51] = 1.0, sing_s6[51] = 1.0, sing_s7[51] = 0.0, sing_s8[51] = 1.0;
    sing_s1[52] = 0.0, sing_s2[52] = 1.0, sing_s3[52] = 0.0, sing_s4[52] = 0.0, sing_s5[52] = 1.0, sing_s6[52] = 1.0, sing_s7[52] = 0.0, sing_s8[52] = 1.0;
    sing_s1[53] = 0.0, sing_s2[53] = 0.0, sing_s3[53] = 1.0, sing_s4[53] = 0.0, sing_s5[53] = 1.0, sing_s6[53] = 1.0, sing_s7[53] = 0.0, sing_s8[53] = 1.0;
    sing_s1[54] = 0.0, sing_s2[54] = 0.0, sing_s3[54] = 0.0, sing_s4[54] = 1.0, sing_s5[54] = 1.0, sing_s6[54] = 1.0, sing_s7[54] = 0.0, sing_s8[54] = 1.0;
    sing_s1[55] = 1.0, sing_s2[55] = 1.0, sing_s3[55] = 0.0, sing_s4[55] = 0.0, sing_s5[55] = 0.0, sing_s6[55] = 0.0, sing_s7[55] = 1.0, sing_s8[55] = 1.0;
    sing_s1[56] = 1.0, sing_s2[56] = 0.0, sing_s3[56] = 1.0, sing_s4[56] = 0.0, sing_s5[56] = 0.0, sing_s6[56] = 0.0, sing_s7[56] = 1.0, sing_s8[56] = 1.0;
    sing_s1[57] = 0.0, sing_s2[57] = 1.0, sing_s3[57] = 1.0, sing_s4[57] = 0.0, sing_s5[57] = 0.0, sing_s6[57] = 0.0, sing_s7[57] = 1.0, sing_s8[57] = 1.0;
    sing_s1[58] = 1.0, sing_s2[58] = 0.0, sing_s3[58] = 0.0, sing_s4[58] = 1.0, sing_s5[58] = 0.0, sing_s6[58] = 0.0, sing_s7[58] = 1.0, sing_s8[58] = 1.0;
    sing_s1[59] = 0.0, sing_s2[59] = 1.0, sing_s3[59] = 0.0, sing_s4[59] = 1.0, sing_s5[59] = 0.0, sing_s6[59] = 0.0, sing_s7[59] = 1.0, sing_s8[59] = 1.0;
    sing_s1[60] = 0.0, sing_s2[60] = 0.0, sing_s3[60] = 1.0, sing_s4[60] = 1.0, sing_s5[60] = 0.0, sing_s6[60] = 0.0, sing_s7[60] = 1.0, sing_s8[60] = 1.0;
    sing_s1[61] = 1.0, sing_s2[61] = 0.0, sing_s3[61] = 0.0, sing_s4[61] = 0.0, sing_s5[61] = 1.0, sing_s6[61] = 0.0, sing_s7[61] = 1.0, sing_s8[61] = 1.0;
    sing_s1[62] = 0.0, sing_s2[62] = 1.0, sing_s3[62] = 0.0, sing_s4[62] = 0.0, sing_s5[62] = 1.0, sing_s6[62] = 0.0, sing_s7[62] = 1.0, sing_s8[62] = 1.0;
    sing_s1[63] = 0.0, sing_s2[63] = 0.0, sing_s3[63] = 1.0, sing_s4[63] = 0.0, sing_s5[63] = 1.0, sing_s6[63] = 0.0, sing_s7[63] = 1.0, sing_s8[63] = 1.0;
    sing_s1[64] = 0.0, sing_s2[64] = 0.0, sing_s3[64] = 0.0, sing_s4[64] = 1.0, sing_s5[64] = 1.0, sing_s6[64] = 0.0, sing_s7[64] = 1.0, sing_s8[64] = 1.0;
    sing_s1[65] = 1.0, sing_s2[65] = 0.0, sing_s3[65] = 0.0, sing_s4[65] = 0.0, sing_s5[65] = 0.0, sing_s6[65] = 1.0, sing_s7[65] = 1.0, sing_s8[65] = 1.0;
    sing_s1[66] = 0.0, sing_s2[66] = 1.0, sing_s3[66] = 0.0, sing_s4[66] = 0.0, sing_s5[66] = 0.0, sing_s6[66] = 1.0, sing_s7[66] = 1.0, sing_s8[66] = 1.0;
    sing_s1[67] = 0.0, sing_s2[67] = 0.0, sing_s3[67] = 1.0, sing_s4[67] = 0.0, sing_s5[67] = 0.0, sing_s6[67] = 1.0, sing_s7[67] = 1.0, sing_s8[67] = 1.0;
    sing_s1[68] = 0.0, sing_s2[68] = 0.0, sing_s3[68] = 0.0, sing_s4[68] = 1.0, sing_s5[68] = 0.0, sing_s6[68] = 1.0, sing_s7[68] = 1.0, sing_s8[68] = 1.0;
    sing_s1[69] = 0.0, sing_s2[69] = 0.0, sing_s3[69] = 0.0, sing_s4[69] = 0.0, sing_s5[69] = 1.0, sing_s6[69] = 1.0, sing_s7[69] = 1.0, sing_s8[69] = 1.0;

// 

    ini = ( double _Complex * ) calloc ( m , sizeof ( double _Complex ));
    Jave = ( double _Complex * ) calloc ( (7337)  , sizeof ( double _Complex ));

    for (ry = 0; ry < 1; ry++) {
        for (rx = 0; rx < 1; rx++) {

// 

            for ( totn = 0; totn < m; totn++ ) {
                H_0[totn] = 0.0;
            }


            for ( dian = 0; dian < size; dian++ ) {
                H_0[ ( size + 1 )*( dian ) ] = 2.0 * u_rep;
            }


            for ( sn = 0; sn < NCn; sn++ ) {
                H_0[ ( size + 1 )*( NCn + 1 )*( sn ) ] = 4.0 * u_rep;
            }


            for ( sn = 0; sn < NCn; sn++ ) {
                H_0[ ( size + 1 )*( NCn - 1 )*( sn + 1 ) ] = 0.0;
            }

// 

            for ( sn = 0; sn < 1; sn++ ) {

                H_0[ ( size + 1 )*( 31 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 31 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 32 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 32 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 33 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 33 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 34 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 34 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*1 + 28 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 29 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 30 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 34 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*1 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*2 + 26 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 27 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 30 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 33 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*2 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*3 + 25 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 27 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 29 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 32 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*3 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*4 + 25 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 26 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 28 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 31 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*4 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*5 + 22 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 23 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 24 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 34 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*5 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*6 + 20 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 21 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 24 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 33 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*6 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*7 + 19 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 21 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 23 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 32 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*7 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*8 + 19 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 20 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 22 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 31 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*8 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*9 + 17 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 18 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 24 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 30 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*9 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*10 + 16 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 18 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 23 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 29 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*10 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*11 + 16 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 17 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 22 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 28 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*11 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*12 + 15 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 18 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 21 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 27 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*12 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*13 + 15 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 17 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 20 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 26 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*13 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*14 + 15 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 16 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 19 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 25 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*14 + NCn - 1 - 66 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*15 + 12 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 13 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 14 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 34 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*15 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*16 + 10 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 11 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 14 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 33 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*16 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*17 + 9 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 11 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 13 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 32 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*17 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*18 + 9 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 10 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 12 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 31 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 69 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*18 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*19 + 7 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 8 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 14 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 30 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*19 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*20 + 6 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 8 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 13 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 29 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*20 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*20 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*21 + 6 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 7 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 12 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 28 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*21 + 68 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*21 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*22 + 5 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 8 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 11 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 27 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*22 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*22 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*23 + 5 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 7 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 10 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 26 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*23 + 67 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*23 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*24 + 5 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 6 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 9 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 25 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 65 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*24 + 66 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*24 + NCn - 1 - 66 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*25 + 3 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 4 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 14 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 24 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*25 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*25 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*26 + 2 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 4 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 13 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 23 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*26 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*26 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*27 + 2 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 3 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 12 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 22 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 54 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*27 + 64 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*27 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*28 + 1 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 4 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 11 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 21 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*28 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*28 + NCn - 1 - 63 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*29 + 1 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 3 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 10 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 20 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 53 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*29 + 63 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*29 + NCn - 1 - 63 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*30 + 1 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 2 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 9 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 19 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 51 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 52 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 61 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*30 + 62 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*30 + NCn - 1 - 62 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*31 + 0 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 4 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 8 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 18 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*31 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*31 + NCn - 1 - 60 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*32 + 0 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 3 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 7 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 17 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 44 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 50 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*32 + 60 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*32 + NCn - 1 - 60 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*33 + 0 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 2 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 6 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 16 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 35 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 42 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 43 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 48 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 49 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 58 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*33 + 59 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*33 + NCn - 1 - 59 ) ] = 3.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*34 + 0 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 1 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 5 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 15 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 36 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 37 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 38 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 39 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 40 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 41 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 45 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 46 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 47 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 55 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 56 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*34 + 57 ) ] = u_rep, H_0[ ( size + 1 )*( NCn*34 + NCn - 1 - 57 ) ] = 3.0 * u_rep;

// 

                H_0[ ( size + 1 )*( NCn*35 + 0 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 0 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 1 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 1 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 5 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 5 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 15 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 15 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 36 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 37 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 38 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 39 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 40 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 41 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 45 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 46 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 47 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 55 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 56 ) ] = u_rep;
                H_0[ ( size + 1 )*( NCn*35 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*35 + NCn - 1 - 57 ) ] = u_rep;

                H_0[ ( size + 1 )*( NCn*36 + 0 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 2 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 6 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 16 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*36 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*36 + NCn - 1 - 59 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*37 + 0 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 3 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 7 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 17 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*37 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*37 + NCn - 1 - 60 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*38 + 0 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 4 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 8 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 18 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*38 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*38 + NCn - 1 - 60 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*39 + 1 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 2 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 9 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 19 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*39 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*39 + NCn - 1 - 62 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*40 + 1 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 3 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 10 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 20 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*40 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*40 + NCn - 1 - 63 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*41 + 1 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 4 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 11 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 21 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*41 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*41 + NCn - 1 - 63 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*42 + 2 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 3 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 12 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 22 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*42 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*42 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*43 + 2 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 4 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 13 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 23 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*43 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*43 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*44 + 3 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 4 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 14 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 24 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*44 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*44 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*45 + 5 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 6 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 9 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 25 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*45 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*45 + NCn - 1 - 66 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*46 + 5 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 7 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 10 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 26 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*46 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*46 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*47 + 5 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 8 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 11 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 27 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*47 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*47 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*48 + 6 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 7 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 12 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 28 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*48 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*48 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*49 + 6 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 8 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 13 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 29 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*49 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*49 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*50 + 7 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 8 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 14 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 30 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*50 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*50 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*51 + 9 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 10 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 12 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 31 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*51 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*51 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*52 + 9 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 11 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 13 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 32 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*52 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*52 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*53 + 10 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 11 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 14 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 33 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*53 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*53 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*54 + 12 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 13 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 14 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 34 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*54 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*54 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*55 + 15 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 16 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 19 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 25 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*55 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*55 + NCn - 1 - 66 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*56 + 15 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 17 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 20 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 26 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*56 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*56 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*57 + 15 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 18 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 21 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 27 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 35 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*57 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*57 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*58 + 16 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 17 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 22 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 28 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*58 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*58 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*59 + 16 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 18 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 23 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 29 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 36 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*59 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*59 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*60 + 17 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 18 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 24 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 30 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 37 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 38 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*60 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*60 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*61 + 19 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 20 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 22 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 31 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*61 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*61 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*62 + 19 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 21 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 23 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 32 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 39 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*62 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*62 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*63 + 20 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 21 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 24 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 33 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 40 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 41 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*63 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*63 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*64 + 22 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 23 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 24 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 34 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 42 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 43 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 44 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*64 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*64 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*65 + 25 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 26 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 28 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 31 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*65 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*65 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*66 + 25 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 27 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 29 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 32 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 45 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 55 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*66 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*66 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*67 + 26 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 27 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 30 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 33 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 46 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 47 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 56 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 57 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*67 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*67 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*68 + 28 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 29 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 30 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 34 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 48 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 49 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 50 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 58 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 59 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 60 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*68 + 69 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*68 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

                H_0[ ( size + 1 )*( NCn*69 + 31 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 32 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 33 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 34 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 51 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 52 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 53 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 54 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 61 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 62 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 63 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 64 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 65 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 66 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 67 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*69 + 68 ) ] = 3.0 * u_rep, H_0[ ( size + 1 )*( NCn*69 + NCn - 1 - 68 ) ] = 1.0 * u_rep;


            }

// 

            for ( hopi = 0; hopi < NCn; hopi++ ) {

                H_0[ size*(NCn*hopi) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi) + 38 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+1) + 0 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+1) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+1) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+1) + 41 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+2) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+2) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+2) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+2) + 43 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+3) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+3) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+3) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+3) + 44 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+4) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+4) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*(NCn*hopi+5) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+5) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+5) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+5) + 47 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+6) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+6) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+6) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+6) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+6) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+6) + 49 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+7) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+7) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+7) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+7) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+7) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+7) + 50 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+8) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+8) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+8) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+8) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+9) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+9) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+9) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+9) + 52 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+10) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+10) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+10) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+10) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+10) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+10) + 53 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+11) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+11) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+11) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+11) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+12) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+12) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+12) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+12) + 54 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+13) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+13) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+13) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+13) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+14) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+14) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+15) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+15) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+15) + 35 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+15) + 57 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+16) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+16) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+16) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+16) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+16) + 36 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+16) + 59 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+17) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+17) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+17) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+17) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+17) + 37 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+17) + 60 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+18) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+18) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+18) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+18) + 38 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+19) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+19) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+19) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+19) + 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+19) + 39 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+19) + 62 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+20) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+20) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+20) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+20) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+20) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+20) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+20) + 40 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+20) + 63 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+21) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+21) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+21) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+21) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+21) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+21) + 41 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+22) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+22) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+22) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+22) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+22) + 42 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+22) + 64 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+23) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+23) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+23) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+23) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+23) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+23) + 43 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+24) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+24) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+24) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+24) + 44 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+25) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+25) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+25) + 45 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+25) + 66 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+26) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+26) + 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+26) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+26) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+26) + 46 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+26) + 67 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+27) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+27) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+27) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+27) + 47 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+28) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+28) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+28) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+28) + 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+28) + 48 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+28) + 68 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+29) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+29) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+29) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+29) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+29) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+29) + 49 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+30) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+30) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+30) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+30) + 50 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+31) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+31) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+31) + 51 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+31) + 69 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+32) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+32) + 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+32) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+32) + 52 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+33) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+33) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+33) + 34 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+33) + 53 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+34) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+34) + 54 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

//

                H_0[ size*(NCn*hopi+35) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+35) + NCn - 1 - 54 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+36) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+36) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+36) + NCn - 1 - 34 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+36) + NCn - 1 - 53 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+37) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+37) + NCn - 1 - 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+37) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+37) + NCn - 1 - 52 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+38) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+38) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+38) + NCn - 1 - 51 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+38) + NCn - 1 - 69 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+39) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+39) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+39) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+39) + NCn - 1 - 50 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+40) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+40) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+40) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+40) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+40) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+40) + NCn - 1 - 49 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+41) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+41) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+41) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+41) + NCn - 1 - 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+41) + NCn - 1 - 48 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+41) + NCn - 1 - 68 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+42) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+42) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+42) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+42) + NCn - 1 - 47 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+43) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+43) + NCn - 1 - 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+43) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+43) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+43) + NCn - 1 - 46 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+43) + NCn - 1 - 67 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+44) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+44) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+44) + NCn - 1 - 45 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+44) + NCn - 1 - 66 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+45) + NCn - 1 - 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+45) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+45) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+45) + NCn - 1 - 44 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+46) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+46) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+46) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+46) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+46) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+46) + NCn - 1 - 43 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+47) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+47) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+47) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+47) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+47) + NCn - 1 - 42 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+47) + NCn - 1 - 64 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+48) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+48) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+48) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+48) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+48) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+48) + NCn - 1 - 41 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+49) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 40 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+49) + NCn - 1 - 63 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+50) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+50) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+50) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+50) + NCn - 1 - 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+50) + NCn - 1 - 39 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+50) + NCn - 1 - 62 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+51) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+51) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+51) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+51) + NCn - 1 - 38 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );

                H_0[ size*(NCn*hopi+52) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+52) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+52) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+52) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+52) + NCn - 1 - 37 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+52) + NCn - 1 - 60 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+53) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+53) + NCn - 1 - 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+53) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+53) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+53) + NCn - 1 - 36 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+53) + NCn - 1 - 59 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+54) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+54) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+54) + NCn - 1 - 35 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] );
                H_0[ size*(NCn*hopi+54) + NCn - 1 - 57 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+55) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+55) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+56) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+56) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+56) + NCn - 1 - 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+56) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+57) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+57) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+57) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+57) + NCn - 1 - 54 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+58) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+58) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+58) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+58) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+59) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+59) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+59) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+59) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+59) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+59) + NCn - 1 - 53 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+60) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+60) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+60) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+60) + NCn - 1 - 52 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+61) + NCn - 1 - 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+61) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+61) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+61) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*(NCn*hopi+62) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+62) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+62) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+62) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+62) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+62) + NCn - 1 - 50 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+63) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+63) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+63) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+63) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+63) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+63) + NCn - 1 - 49 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+64) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+64) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+64) + NCn - 1 - 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*(NCn*hopi+64) + NCn - 1 - 47 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+65) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+65) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*(NCn*hopi+66) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+66) + NCn - 1 - 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+66) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+66) + NCn - 1 - 44 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+67) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+67) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+67) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+67) + NCn - 1 - 43 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+68) + NCn - 1 - 0 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+68) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+68) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+68) + NCn - 1 - 41 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );

                H_0[ size*(NCn*hopi+69) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+69) + NCn - 1 - 38 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) );


            }

// 

            for ( hopi = 0; hopi < NCn; hopi++ ) {

                H_0[ size*hopi + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*hopi + 38*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*1 ) + 0*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 41*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*2 ) + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 43*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*3 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*3 ) + 4*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*3 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*3 ) + 44*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*4 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*4 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*( hopi + NCn*5 ) + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 15*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 47*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*6 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 49*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*7 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 50*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*8 ) + 4*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*9 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*9 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*9 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*9 ) + 52*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*10 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 53*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*11 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*12 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*12 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*12 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*12 ) + 54*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*13 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 14*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*14 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*14 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*15 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*15 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*15 ) + 35*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*15 ) + 57*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*16 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 15*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 36*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 59*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*17 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 37*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 60*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*18 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*18 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*18 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*18 ) + 38*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*19 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*19 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*19 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*19 ) + 25*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*19 ) + 39*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*19 ) + 62*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*20 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 40*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*20 ) + 63*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*21 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*21 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*21 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*21 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*21 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*21 ) + 41*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*22 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*22 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*22 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*22 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*22 ) + 42*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*22 ) + 64*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*23 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*23 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*23 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*23 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*23 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*23 ) + 43*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*24 ) + 14*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*24 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*24 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*24 ) + 44*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*25 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*25 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*25 ) + 45*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*25 ) + 66*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*26 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*26 ) + 25*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*26 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*26 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*26 ) + 46*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*26 ) + 67*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*27 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*27 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*27 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*27 ) + 47*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*28 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*28 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*28 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*28 ) + 31*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*28 ) + 48*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*28 ) + 68*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*29 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*29 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*29 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*29 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*29 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*29 ) + 49*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*30 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*30 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*30 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*30 ) + 50*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*31 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*31 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*31 ) + 51*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*31 ) + 69*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*32 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*32 ) + 31*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*32 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*32 ) + 52*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*33 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*33 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*33 ) + 34*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*33 ) + 53*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*34 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*34 ) + 54*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

//

                H_0[ size*( hopi + NCn*35 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*35 ) + ( NCn - 1 - 54 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*36 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*36 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*36 ) + ( NCn - 1 - 34 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*36 ) + ( NCn - 1 - 53 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*37 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*37 ) + ( NCn - 1 - 31 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*37 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*37 ) + ( NCn - 1 - 52 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*38 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*38 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*38 ) + ( NCn - 1 - 51 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*38 ) + ( NCn - 1 - 69 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*39 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*39 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*39 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*39 ) + ( NCn - 1 - 50 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*40 ) + ( NCn - 1 - 49 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 31 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 48 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*41 ) + ( NCn - 1 - 68 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*42 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*42 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*42 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*42 ) + ( NCn - 1 - 47 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 25 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 46 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*43 ) + ( NCn - 1 - 67 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*44 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*44 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*44 ) + ( NCn - 1 - 45 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*44 ) + ( NCn - 1 - 66 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*45 ) + ( NCn - 1 - 14 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*45 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*45 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*45 ) + ( NCn - 1 - 44 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*46 ) + ( NCn - 1 - 43 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 42 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*47 ) + ( NCn - 1 - 64 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*48 ) + ( NCn - 1 - 41 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 40 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*49 ) + ( NCn - 1 - 63 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 25 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 39 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*50 ) + ( NCn - 1 - 62 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*51 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*51 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*51 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*51 ) + ( NCn - 1 - 38 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );

                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 37 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*52 ) + ( NCn - 1 - 60 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 15 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 36 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*53 ) + ( NCn - 1 - 59 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*54 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*54 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*54 ) + ( NCn - 1 - 35 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] );
                H_0[ size*( hopi + NCn*54 ) + ( NCn - 1 - 57 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*55 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*55 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*56 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*56 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*56 ) + ( NCn - 1 - 14 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*56 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*57 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*57 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*57 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*57 ) + ( NCn - 1 - 54 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*58 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*58 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*58 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*58 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*59 ) + ( NCn - 1 - 53 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*60 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*60 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*60 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*60 ) + ( NCn - 1 - 52 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*61 ) + ( NCn - 1 - 4 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*61 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*61 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*61 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*62 ) + ( NCn - 1 - 50 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*63 ) + ( NCn - 1 - 49 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*64 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*64 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*64 ) + ( NCn - 1 - 15 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*( hopi + NCn*64 ) + ( NCn - 1 - 47 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*65 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*65 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*( hopi + NCn*66 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*66 ) + ( NCn - 1 - 4 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*66 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*66 ) + ( NCn - 1 - 44 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*67 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*67 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*67 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*67 ) + ( NCn - 1 - 43 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*68 ) + ( NCn - 1 - 0 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*68 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*68 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*68 ) + ( NCn - 1 - 41 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );

                H_0[ size*( hopi + NCn*69 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*69 ) + ( NCn - 1 - 38 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) );


            }

// end Hamiltonian


            info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', n, H_0, lda, E_eigval);
            if(info !=0){printf("LAPACKE_zheev error %d\n",info); }

            for ( rhi = 0 ; rhi < size ; rhi++ ) {
                for ( rhj = 0 ; rhj < size ; rhj++ ) {
                    gsl_matrix_complex_set(U, rhi, rhj, gsl_complex_rect(creal(H_0[rhi * size + rhj]),cimag(H_0[rhi * size + rhj])));
                }
            }

            gsl_matrix_complex_set_all(RHO_TILDE, GSL_COMPLEX_ZERO);

// ground state

            gsl_matrix_complex_set(RHO_TILDE, 0, 0, gsl_complex_rect(1.0,0.0));

            gsl_matrix_complex_set_all(MIDDLE, GSL_COMPLEX_ZERO);
            gsl_matrix_complex_set_all(RHO, GSL_COMPLEX_ZERO);

            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, RHO_TILDE,
                    GSL_COMPLEX_ZERO, MIDDLE);

            gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, MIDDLE, U,
                    GSL_COMPLEX_ZERO, RHO);

            for ( rhi = 0 ; rhi < size ; rhi++ ) {
                for ( rhj = 0 ; rhj < size ; rhj++ ) {
                    ini[rhi * size + rhj] = GSL_REAL(gsl_matrix_complex_get(RHO,rhi,rhj)) + GSL_IMAG(gsl_matrix_complex_get(RHO,rhi,rhj))*I;
                }
            }

            for (q = 0; q < 7337; q++) {// time evolve

                if (t0<Tfull/2.0) {
                    Ax = -(Ex/wpump) * (cos(wpump*t0+pi*0.25)) * (pow(cos((pi/Tfull)*t0),4.0));
                }
                else{
                    Ax = 0.0;
                }

// current operator


                gsl_matrix_complex_set_all(J_X, GSL_COMPLEX_ZERO);


                for ( hopi = 0; hopi < NCn; hopi++ ) {

                    gsl_matrix_complex_set(J_X, 0 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 0 + NCn*hopi , 38 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 0 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 41 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 43 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 44 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 47 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 49 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 50 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 52 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 53 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 54 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 35 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 57 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 36 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 59 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 37 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 60 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 38 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 25 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 39 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 62 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 40 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 20 + NCn*hopi , 63 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 21 + NCn*hopi , 41 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 42 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 22 + NCn*hopi , 64 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 23 + NCn*hopi , 43 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 24 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 24 + NCn*hopi , 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 24 + NCn*hopi , 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 24 + NCn*hopi , 44 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 25 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 25 + NCn*hopi , 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 25 + NCn*hopi , 45 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 25 + NCn*hopi , 66 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 25 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 46 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 26 + NCn*hopi , 67 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 27 + NCn*hopi , 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 27 + NCn*hopi , 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 27 + NCn*hopi , 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 27 + NCn*hopi , 47 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 31 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 48 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 28 + NCn*hopi , 68 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 29 + NCn*hopi , 49 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 30 + NCn*hopi , 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 30 + NCn*hopi , 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 30 + NCn*hopi , 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 30 + NCn*hopi , 50 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 31 + NCn*hopi , 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 31 + NCn*hopi , 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 31 + NCn*hopi , 51 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 31 + NCn*hopi , 69 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 32 + NCn*hopi , 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 32 + NCn*hopi , 31 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 32 + NCn*hopi , 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 32 + NCn*hopi , 52 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 33 + NCn*hopi , 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 33 + NCn*hopi , 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 33 + NCn*hopi , 34 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 33 + NCn*hopi , 53 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 34 + NCn*hopi , 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 34 + NCn*hopi , 54 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));


                    gsl_matrix_complex_set(J_X, 35 + NCn*hopi , NCn - 1 - 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 35 + NCn*hopi , NCn - 1 - 54 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 36 + NCn*hopi , NCn - 1 - 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 36 + NCn*hopi , NCn - 1 - 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 36 + NCn*hopi , NCn - 1 - 34 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 36 + NCn*hopi , NCn - 1 - 53 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 37 + NCn*hopi , NCn - 1 - 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 37 + NCn*hopi , NCn - 1 - 31 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 37 + NCn*hopi , NCn - 1 - 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 37 + NCn*hopi , NCn - 1 - 52 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 38 + NCn*hopi , NCn - 1 - 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 38 + NCn*hopi , NCn - 1 - 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 38 + NCn*hopi , NCn - 1 - 51 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 38 + NCn*hopi , NCn - 1 - 69 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 39 + NCn*hopi , NCn - 1 - 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 39 + NCn*hopi , NCn - 1 - 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 39 + NCn*hopi , NCn - 1 - 33 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 39 + NCn*hopi , NCn - 1 - 50 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 32 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 40 + NCn*hopi , NCn - 1 - 49 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 31 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 48 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 41 + NCn*hopi , NCn - 1 - 68 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 42 + NCn*hopi , NCn - 1 - 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 42 + NCn*hopi , NCn - 1 - 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 42 + NCn*hopi , NCn - 1 - 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 42 + NCn*hopi , NCn - 1 - 47 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 25 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 46 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 43 + NCn*hopi , NCn - 1 - 67 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 44 + NCn*hopi , NCn - 1 - 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 44 + NCn*hopi , NCn - 1 - 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 44 + NCn*hopi , NCn - 1 - 45 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 44 + NCn*hopi , NCn - 1 - 66 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 45 + NCn*hopi , NCn - 1 - 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 45 + NCn*hopi , NCn - 1 - 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 45 + NCn*hopi , NCn - 1 - 30 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 45 + NCn*hopi , NCn - 1 - 44 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 29 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 46 + NCn*hopi , NCn - 1 - 43 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 28 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 42 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 47 + NCn*hopi , NCn - 1 - 64 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 27 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 48 + NCn*hopi , NCn - 1 - 41 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 26 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 40 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 49 + NCn*hopi , NCn - 1 - 63 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 25 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 39 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 50 + NCn*hopi , NCn - 1 - 62 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 51 + NCn*hopi , NCn - 1 - 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 51 + NCn*hopi , NCn - 1 - 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 51 + NCn*hopi , NCn - 1 - 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 51 + NCn*hopi , NCn - 1 - 38 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 37 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 52 + NCn*hopi , NCn - 1 - 60 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 36 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 53 + NCn*hopi , NCn - 1 - 59 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 54 + NCn*hopi , NCn - 1 - 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 54 + NCn*hopi , NCn - 1 - 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 54 + NCn*hopi , NCn - 1 - 35 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s8[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 54 + NCn*hopi , NCn - 1 - 57 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 55 + NCn*hopi , NCn - 1 - 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 55 + NCn*hopi , NCn - 1 - 24 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 56 + NCn*hopi , NCn - 1 - 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 56 + NCn*hopi , NCn - 1 - 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 56 + NCn*hopi , NCn - 1 - 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 56 + NCn*hopi , NCn - 1 - 23 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 57 + NCn*hopi , NCn - 1 - 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 57 + NCn*hopi , NCn - 1 - 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 57 + NCn*hopi , NCn - 1 - 22 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 57 + NCn*hopi , NCn - 1 - 54 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 58 + NCn*hopi , NCn - 1 - 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 58 + NCn*hopi , NCn - 1 - 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 58 + NCn*hopi , NCn - 1 - 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 58 + NCn*hopi , NCn - 1 - 21 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 20 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 59 + NCn*hopi , NCn - 1 - 53 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 60 + NCn*hopi , NCn - 1 - 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 60 + NCn*hopi , NCn - 1 - 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 60 + NCn*hopi , NCn - 1 - 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 60 + NCn*hopi , NCn - 1 - 52 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 61 + NCn*hopi , NCn - 1 - 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 61 + NCn*hopi , NCn - 1 - 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 61 + NCn*hopi , NCn - 1 - 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 61 + NCn*hopi , NCn - 1 - 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 62 + NCn*hopi , NCn - 1 - 50 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 63 + NCn*hopi , NCn - 1 - 49 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 64 + NCn*hopi , NCn - 1 - 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 64 + NCn*hopi , NCn - 1 - 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 64 + NCn*hopi , NCn - 1 - 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 64 + NCn*hopi , NCn - 1 - 47 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 65 + NCn*hopi , NCn - 1 - 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 65 + NCn*hopi , NCn - 1 - 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 66 + NCn*hopi , NCn - 1 - 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 66 + NCn*hopi , NCn - 1 - 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 66 + NCn*hopi , NCn - 1 - 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 66 + NCn*hopi , NCn - 1 - 44 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 67 + NCn*hopi , NCn - 1 - 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 67 + NCn*hopi , NCn - 1 - 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 67 + NCn*hopi , NCn - 1 - 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 67 + NCn*hopi , NCn - 1 - 43 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 68 + NCn*hopi , NCn - 1 - 0 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 68 + NCn*hopi , NCn - 1 - 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 68 + NCn*hopi , NCn - 1 - 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 68 + NCn*hopi , NCn - 1 - 41 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 69 + NCn*hopi , NCn - 1 - 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 69 + NCn*hopi , NCn - 1 - 38 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) ));



                }

//

                for ( hopi = 0; hopi < NCn; hopi++ ) {

                    gsl_matrix_complex_set(J_X, NCn*0 + hopi , 1*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*0 + hopi , 38*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 0*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 2*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 5*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 41*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 1*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 3*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 6*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 43*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 2*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 4*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 7*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 44*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 3*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 8*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 1*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 6*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 15*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 47*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] )  ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 2*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 5*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 7*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 9*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 16*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 49*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 3*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 6*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 8*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 10*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 17*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 50*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 4*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 7*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 11*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 18*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 6*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 10*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 19*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 52*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 7*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 9*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 11*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 12*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 53*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 8*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 10*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 13*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 21*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 10*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 13*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 22*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 54*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 11*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 12*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 14*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 23*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 13*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 24*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 5*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 16*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 35*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 57*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 6*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 15*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 17*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 19*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 36*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 59*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 7*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 16*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 18*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 37*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 60*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 8*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 17*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 21*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 38*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 9*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 16*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 25*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 39*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 62*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 10*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 17*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 19*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 21*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 22*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 26*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 40*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*20 + hopi , 63*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 11*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 18*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 23*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 27*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*21 + hopi , 41*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 12*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 23*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 28*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 42*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*22 + hopi , 64*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 13*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 21*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 22*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 24*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 29*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*23 + hopi , 43*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*24 + hopi , 14*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*24 + hopi , 23*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*24 + hopi , 30*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*24 + hopi , 44*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*25 + hopi , 19*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*25 + hopi , 26*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*25 + hopi , 45*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*25 + hopi , 66*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 20*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 25*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 27*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 28*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 46*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*26 + hopi , 67*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*27 + hopi , 21*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*27 + hopi , 26*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*27 + hopi , 29*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*27 + hopi , 47*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 22*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 26*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 29*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 31*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 48*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*28 + hopi , 68*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 23*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 27*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 28*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 30*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 32*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*29 + hopi , 49*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*30 + hopi , 24*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*30 + hopi , 29*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*30 + hopi , 33*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*30 + hopi , 50*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*31 + hopi , 28*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*31 + hopi , 32*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*31 + hopi , 51*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*31 + hopi , 69*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*32 + hopi , 29*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*32 + hopi , 31*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*32 + hopi , 33*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*32 + hopi , 52*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*33 + hopi , 30*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*33 + hopi , 32*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*33 + hopi , 34*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*33 + hopi , 53*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*34 + hopi , 33*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*34 + hopi , 54*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));


                    gsl_matrix_complex_set(J_X, NCn*35 + hopi , ( NCn - 1 - 33 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*35 + hopi , ( NCn - 1 - 54 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*36 + hopi , ( NCn - 1 - 30 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*36 + hopi , ( NCn - 1 - 32 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*36 + hopi , ( NCn - 1 - 34 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*36 + hopi , ( NCn - 1 - 53 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*37 + hopi , ( NCn - 1 - 29 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*37 + hopi , ( NCn - 1 - 31 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*37 + hopi , ( NCn - 1 - 33 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*37 + hopi , ( NCn - 1 - 52 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*38 + hopi , ( NCn - 1 - 28 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*38 + hopi , ( NCn - 1 - 32 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*38 + hopi , ( NCn - 1 - 51 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*38 + hopi , ( NCn - 1 - 69 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*39 + hopi , ( NCn - 1 - 24 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*39 + hopi , ( NCn - 1 - 29 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*39 + hopi , ( NCn - 1 - 33 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*39 + hopi , ( NCn - 1 - 50 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 23 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 27 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 28 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 30 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 32 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*40 + hopi , ( NCn - 1 - 49 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 22 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 26 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 29 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 31 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 48 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*41 + hopi , ( NCn - 1 - 68 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*42 + hopi , ( NCn - 1 - 21 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*42 + hopi , ( NCn - 1 - 26 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*42 + hopi , ( NCn - 1 - 29 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*42 + hopi , ( NCn - 1 - 47 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 25 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 27 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 28 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 46 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*43 + hopi , ( NCn - 1 - 67 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*44 + hopi , ( NCn - 1 - 19 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*44 + hopi , ( NCn - 1 - 26 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*44 + hopi , ( NCn - 1 - 45 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*44 + hopi , ( NCn - 1 - 66 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*45 + hopi , ( NCn - 1 - 14 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*45 + hopi , ( NCn - 1 - 23 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*45 + hopi , ( NCn - 1 - 30 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*45 + hopi , ( NCn - 1 - 44 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 13 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 21 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 22 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 24 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 29 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*46 + hopi , ( NCn - 1 - 43 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 12 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 23 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 28 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 42 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*47 + hopi , ( NCn - 1 - 64 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 11 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 18 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 23 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 27 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*48 + hopi , ( NCn - 1 - 41 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 10 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 17 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 19 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 21 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 22 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 26 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 40 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*49 + hopi , ( NCn - 1 - 63 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 9 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 16 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 25 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 39 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*50 + hopi , ( NCn - 1 - 62 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*51 + hopi , ( NCn - 1 - 8 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*51 + hopi , ( NCn - 1 - 17 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*51 + hopi , ( NCn - 1 - 21 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*51 + hopi , ( NCn - 1 - 38 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 7 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 16 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 18 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 37 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*52 + hopi , ( NCn - 1 - 60 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 ,  ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 6 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 15 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 17 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 19 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 36 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*53 + hopi , ( NCn - 1 - 59 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*54 + hopi , ( NCn - 1 - 5 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*54 + hopi , ( NCn - 1 - 16 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*54 + hopi , ( NCn - 1 - 35 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s7[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*54 + hopi , ( NCn - 1 - 57 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*55 + hopi , ( NCn - 1 - 13 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*55 + hopi , ( NCn - 1 - 24 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*56 + hopi , ( NCn - 1 - 11 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*56 + hopi , ( NCn - 1 - 12 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*56 + hopi , ( NCn - 1 - 14 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*56 + hopi , ( NCn - 1 - 23 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*57 + hopi , ( NCn - 1 - 10 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*57 + hopi , ( NCn - 1 - 13 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*57 + hopi , ( NCn - 1 - 22 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*57 + hopi , ( NCn - 1 - 54 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*58 + hopi , ( NCn - 1 - 8 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*58 + hopi , ( NCn - 1 - 10 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*58 + hopi , ( NCn - 1 - 13 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*58 + hopi , ( NCn - 1 - 21 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 7 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 9 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 11 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 12 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 20 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*59 + hopi , ( NCn - 1 - 53 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*60 + hopi , ( NCn - 1 - 6 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*60 + hopi , ( NCn - 1 - 10 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*60 + hopi , ( NCn - 1 - 19 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*60 + hopi , ( NCn - 1 - 52 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*61 + hopi , ( NCn - 1 - 4 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*61 + hopi , ( NCn - 1 - 7 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*61 + hopi , ( NCn - 1 - 11 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*61 + hopi , ( NCn - 1 - 18 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 3 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 6 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 8 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 10 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 17 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*62 + hopi , ( NCn - 1 - 50 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 2 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 5 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 7 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 9 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 16 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*63 + hopi , ( NCn - 1 - 49 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*64 + hopi , ( NCn - 1 - 1 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*64 + hopi , ( NCn - 1 - 6 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*64 + hopi , ( NCn - 1 - 15 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*64 + hopi , ( NCn - 1 - 47 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*65 + hopi , ( NCn - 1 - 3 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*65 + hopi , ( NCn - 1 - 8 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*66 + hopi , ( NCn - 1 - 2 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*66 + hopi , ( NCn - 1 - 4 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*66 + hopi , ( NCn - 1 - 7 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*66 + hopi , ( NCn - 1 - 44 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*67 + hopi , ( NCn - 1 - 1 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*67 + hopi , ( NCn - 1 - 3 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*67 + hopi , ( NCn - 1 - 6 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*67 + hopi , ( NCn - 1 - 43 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*68 + hopi , ( NCn - 1 - 0 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*68 + hopi , ( NCn - 1 - 2 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*68 + hopi , ( NCn - 1 - 5 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*68 + hopi , ( NCn - 1 - 41 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*69 + hopi , ( NCn - 1 - 1 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*69 + hopi , ( NCn - 1 - 38 )*NCn + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) ));



                }

// end of current

                gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, J_X, RHO,
                        GSL_COMPLEX_ZERO, CURRENT_X);


                for ( dian = 0; dian < size; dian++ ) {

                    Jave[q] = Jave[q] + GSL_REAL(gsl_matrix_complex_get(CURRENT_X, dian , dian ));

                }



                t1 = t0 + dt;
                sol = rk4veckxky(t0, m, ini, dt, tau0, kx, ky, alpha, test_f);
                t0 = t1;
                for (j = 0; j < m; j++) {
                    ini[j] = sol[j];
                }
                free(sol);

                for ( rhi = 0 ; rhi < size ; rhi++ ) {
                    for ( rhj = 0 ; rhj < size ; rhj++ ) {
                        gsl_matrix_complex_set(RHO, rhi, rhj, gsl_complex_rect(creal(ini[rhi * size + rhj]),cimag(ini[rhi * size + rhj])));
                    }
                }

            }
            t0 = -2067.1*2.0;

        }
    }

    free(ini);
    gsl_matrix_complex_free(U);
    gsl_matrix_complex_free(RHO_TILDE);
    gsl_matrix_complex_free(MIDDLE);
    gsl_matrix_complex_free(RHO);

    gsl_matrix_complex_free(J_X);
    gsl_matrix_complex_free(J_Y);
    gsl_matrix_complex_free(CURRENT_X);
    gsl_matrix_complex_free(CURRENT_Y);

    return Jave;

}
double _Complex *test_f ( double t, int m, double _Complex u[], double tau, double kx, double ky, double alpha ){

    double _Complex *arr;

    double *sub_D;

    double pi = 3.14159265358979323846;

    double latt_const = 7.56;// lattice constant (a.u.)
    int rhi, rhj, oparai, oparaj;

    double ellip = 0.0;
    double Ex = ( 3.5/6.5 )*0.00126410*alpha/(sqrt(1.0+ellip*ellip));//optical field strength
    double wpump = 0.26/27.2114; //pump frequency
    double Tfull = 4134.2*2.0; //pulse length
    double Ax;

    double u_rep = 0.0*0.52/27.2114;// Hubbard U
    double t_hub = 0.52/27.2114;// hopping t

    double gamma = 0.0 * 0.52/27.2114;// electron-environment scattering rate


    arr = ( double _Complex * ) calloc ( m , sizeof ( double _Complex ) );

    sub_D = ( double * ) calloc ( m , sizeof ( double ) );

    if (t<Tfull/2.0) {
        Ax = -(Ex/wpump) * (cos(wpump*t+pi*0.25)) * (pow(cos((pi/Tfull)*t),4.0));
    }
    else{
        Ax = 0.0;
    }

    lapack_int n2 = 4900;
    lapack_complex_double H_KA [n2*n2];//time-dep-hamil

// for hamiltonian

    int totn, dian, sn, NCn = 70, hopi;
    int size2 = 4900;

// 8 site

    double site1[n2], site2[n2], site3[n2], site4[n2], site5[n2], site6[n2], site7[n2], site8[n2];
    double sing_s1[NCn], sing_s2[NCn], sing_s3[NCn], sing_s4[NCn], sing_s5[NCn], sing_s6[NCn], sing_s7[NCn], sing_s8[NCn];

    sing_s1[0] = 1.0,  sing_s2[0] = 1.0,  sing_s3[0] = 1.0,  sing_s4[0] = 1.0,  sing_s5[0] = 0.0,  sing_s6[0] = 0.0,  sing_s7[0] = 0.0,  sing_s8[0] = 0.0;
    sing_s1[1] = 1.0,  sing_s2[1] = 1.0,  sing_s3[1] = 1.0,  sing_s4[1] = 0.0,  sing_s5[1] = 1.0,  sing_s6[1] = 0.0,  sing_s7[1] = 0.0,  sing_s8[1] = 0.0;
    sing_s1[2] = 1.0,  sing_s2[2] = 1.0,  sing_s3[2] = 0.0,  sing_s4[2] = 1.0,  sing_s5[2] = 1.0,  sing_s6[2] = 0.0,  sing_s7[2] = 0.0,  sing_s8[2] = 0.0;
    sing_s1[3] = 1.0,  sing_s2[3] = 0.0,  sing_s3[3] = 1.0,  sing_s4[3] = 1.0,  sing_s5[3] = 1.0,  sing_s6[3] = 0.0,  sing_s7[3] = 0.0,  sing_s8[3] = 0.0;
    sing_s1[4] = 0.0,  sing_s2[4] = 1.0,  sing_s3[4] = 1.0,  sing_s4[4] = 1.0,  sing_s5[4] = 1.0,  sing_s6[4] = 0.0,  sing_s7[4] = 0.0,  sing_s8[4] = 0.0;
    sing_s1[5] = 1.0,  sing_s2[5] = 1.0,  sing_s3[5] = 1.0,  sing_s4[5] = 0.0,  sing_s5[5] = 0.0,  sing_s6[5] = 1.0,  sing_s7[5] = 0.0,  sing_s8[5] = 0.0;
    sing_s1[6] = 1.0,  sing_s2[6] = 1.0,  sing_s3[6] = 0.0,  sing_s4[6] = 1.0,  sing_s5[6] = 0.0,  sing_s6[6] = 1.0,  sing_s7[6] = 0.0,  sing_s8[6] = 0.0;
    sing_s1[7] = 1.0,  sing_s2[7] = 0.0,  sing_s3[7] = 1.0,  sing_s4[7] = 1.0,  sing_s5[7] = 0.0,  sing_s6[7] = 1.0,  sing_s7[7] = 0.0,  sing_s8[7] = 0.0;
    sing_s1[8] = 0.0,  sing_s2[8] = 1.0,  sing_s3[8] = 1.0,  sing_s4[8] = 1.0,  sing_s5[8] = 0.0,  sing_s6[8] = 1.0,  sing_s7[8] = 0.0,  sing_s8[8] = 0.0;
    sing_s1[9] = 1.0,  sing_s2[9] = 1.0,  sing_s3[9] = 0.0,  sing_s4[9] = 0.0,  sing_s5[9] = 1.0,  sing_s6[9] = 1.0,  sing_s7[9] = 0.0,  sing_s8[9] = 0.0;
    sing_s1[10] = 1.0, sing_s2[10] = 0.0, sing_s3[10] = 1.0, sing_s4[10] = 0.0, sing_s5[10] = 1.0, sing_s6[10] = 1.0, sing_s7[10] = 0.0, sing_s8[10] = 0.0;
    sing_s1[11] = 0.0, sing_s2[11] = 1.0, sing_s3[11] = 1.0, sing_s4[11] = 0.0, sing_s5[11] = 1.0, sing_s6[11] = 1.0, sing_s7[11] = 0.0, sing_s8[11] = 0.0;
    sing_s1[12] = 1.0, sing_s2[12] = 0.0, sing_s3[12] = 0.0, sing_s4[12] = 1.0, sing_s5[12] = 1.0, sing_s6[12] = 1.0, sing_s7[12] = 0.0, sing_s8[12] = 0.0;
    sing_s1[13] = 0.0, sing_s2[13] = 1.0, sing_s3[13] = 0.0, sing_s4[13] = 1.0, sing_s5[13] = 1.0, sing_s6[13] = 1.0, sing_s7[13] = 0.0, sing_s8[13] = 0.0;
    sing_s1[14] = 0.0, sing_s2[14] = 0.0, sing_s3[14] = 1.0, sing_s4[14] = 1.0, sing_s5[14] = 1.0, sing_s6[14] = 1.0, sing_s7[14] = 0.0, sing_s8[14] = 0.0;
    sing_s1[15] = 1.0, sing_s2[15] = 1.0, sing_s3[15] = 1.0, sing_s4[15] = 0.0, sing_s5[15] = 0.0, sing_s6[15] = 0.0, sing_s7[15] = 1.0, sing_s8[15] = 0.0;
    sing_s1[16] = 1.0, sing_s2[16] = 1.0, sing_s3[16] = 0.0, sing_s4[16] = 1.0, sing_s5[16] = 0.0, sing_s6[16] = 0.0, sing_s7[16] = 1.0, sing_s8[16] = 0.0;
    sing_s1[17] = 1.0, sing_s2[17] = 0.0, sing_s3[17] = 1.0, sing_s4[17] = 1.0, sing_s5[17] = 0.0, sing_s6[17] = 0.0, sing_s7[17] = 1.0, sing_s8[17] = 0.0;
    sing_s1[18] = 0.0, sing_s2[18] = 1.0, sing_s3[18] = 1.0, sing_s4[18] = 1.0, sing_s5[18] = 0.0, sing_s6[18] = 0.0, sing_s7[18] = 1.0, sing_s8[18] = 0.0;
    sing_s1[19] = 1.0, sing_s2[19] = 1.0, sing_s3[19] = 0.0, sing_s4[19] = 0.0, sing_s5[19] = 1.0, sing_s6[19] = 0.0, sing_s7[19] = 1.0, sing_s8[19] = 0.0;
    sing_s1[20] = 1.0, sing_s2[20] = 0.0, sing_s3[20] = 1.0, sing_s4[20] = 0.0, sing_s5[20] = 1.0, sing_s6[20] = 0.0, sing_s7[20] = 1.0, sing_s8[20] = 0.0;
    sing_s1[21] = 0.0, sing_s2[21] = 1.0, sing_s3[21] = 1.0, sing_s4[21] = 0.0, sing_s5[21] = 1.0, sing_s6[21] = 0.0, sing_s7[21] = 1.0, sing_s8[21] = 0.0;
    sing_s1[22] = 1.0, sing_s2[22] = 0.0, sing_s3[22] = 0.0, sing_s4[22] = 1.0, sing_s5[22] = 1.0, sing_s6[22] = 0.0, sing_s7[22] = 1.0, sing_s8[22] = 0.0;
    sing_s1[23] = 0.0, sing_s2[23] = 1.0, sing_s3[23] = 0.0, sing_s4[23] = 1.0, sing_s5[23] = 1.0, sing_s6[23] = 0.0, sing_s7[23] = 1.0, sing_s8[23] = 0.0;
    sing_s1[24] = 0.0, sing_s2[24] = 0.0, sing_s3[24] = 1.0, sing_s4[24] = 1.0, sing_s5[24] = 1.0, sing_s6[24] = 0.0, sing_s7[24] = 1.0, sing_s8[24] = 0.0;
    sing_s1[25] = 1.0, sing_s2[25] = 1.0, sing_s3[25] = 0.0, sing_s4[25] = 0.0, sing_s5[25] = 0.0, sing_s6[25] = 1.0, sing_s7[25] = 1.0, sing_s8[25] = 0.0;
    sing_s1[26] = 1.0, sing_s2[26] = 0.0, sing_s3[26] = 1.0, sing_s4[26] = 0.0, sing_s5[26] = 0.0, sing_s6[26] = 1.0, sing_s7[26] = 1.0, sing_s8[26] = 0.0;
    sing_s1[27] = 0.0, sing_s2[27] = 1.0, sing_s3[27] = 1.0, sing_s4[27] = 0.0, sing_s5[27] = 0.0, sing_s6[27] = 1.0, sing_s7[27] = 1.0, sing_s8[27] = 0.0;
    sing_s1[28] = 1.0, sing_s2[28] = 0.0, sing_s3[28] = 0.0, sing_s4[28] = 1.0, sing_s5[28] = 0.0, sing_s6[28] = 1.0, sing_s7[28] = 1.0, sing_s8[28] = 0.0;
    sing_s1[29] = 0.0, sing_s2[29] = 1.0, sing_s3[29] = 0.0, sing_s4[29] = 1.0, sing_s5[29] = 0.0, sing_s6[29] = 1.0, sing_s7[29] = 1.0, sing_s8[29] = 0.0;
    sing_s1[30] = 0.0, sing_s2[30] = 0.0, sing_s3[30] = 1.0, sing_s4[30] = 1.0, sing_s5[30] = 0.0, sing_s6[30] = 1.0, sing_s7[30] = 1.0, sing_s8[30] = 0.0;
    sing_s1[31] = 1.0, sing_s2[31] = 0.0, sing_s3[31] = 0.0, sing_s4[31] = 0.0, sing_s5[31] = 1.0, sing_s6[31] = 1.0, sing_s7[31] = 1.0, sing_s8[31] = 0.0;
    sing_s1[32] = 0.0, sing_s2[32] = 1.0, sing_s3[32] = 0.0, sing_s4[32] = 0.0, sing_s5[32] = 1.0, sing_s6[32] = 1.0, sing_s7[32] = 1.0, sing_s8[32] = 0.0;
    sing_s1[33] = 0.0, sing_s2[33] = 0.0, sing_s3[33] = 1.0, sing_s4[33] = 0.0, sing_s5[33] = 1.0, sing_s6[33] = 1.0, sing_s7[33] = 1.0, sing_s8[33] = 0.0;
    sing_s1[34] = 0.0, sing_s2[34] = 0.0, sing_s3[34] = 0.0, sing_s4[34] = 1.0, sing_s5[34] = 1.0, sing_s6[34] = 1.0, sing_s7[34] = 1.0, sing_s8[34] = 0.0;

    sing_s1[35] = 1.0, sing_s2[35] = 1.0, sing_s3[35] = 1.0, sing_s4[35] = 0.0, sing_s5[35] = 0.0, sing_s6[35] = 0.0, sing_s7[35] = 0.0, sing_s8[35] = 1.0;
    sing_s1[36] = 1.0, sing_s2[36] = 1.0, sing_s3[36] = 0.0, sing_s4[36] = 1.0, sing_s5[36] = 0.0, sing_s6[36] = 0.0, sing_s7[36] = 0.0, sing_s8[36] = 1.0;
    sing_s1[37] = 1.0, sing_s2[37] = 0.0, sing_s3[37] = 1.0, sing_s4[37] = 1.0, sing_s5[37] = 0.0, sing_s6[37] = 0.0, sing_s7[37] = 0.0, sing_s8[37] = 1.0;
    sing_s1[38] = 0.0, sing_s2[38] = 1.0, sing_s3[38] = 1.0, sing_s4[38] = 1.0, sing_s5[38] = 0.0, sing_s6[38] = 0.0, sing_s7[38] = 0.0, sing_s8[38] = 1.0;
    sing_s1[39] = 1.0, sing_s2[39] = 1.0, sing_s3[39] = 0.0, sing_s4[39] = 0.0, sing_s5[39] = 1.0, sing_s6[39] = 0.0, sing_s7[39] = 0.0, sing_s8[39] = 1.0;
    sing_s1[40] = 1.0, sing_s2[40] = 0.0, sing_s3[40] = 1.0, sing_s4[40] = 0.0, sing_s5[40] = 1.0, sing_s6[40] = 0.0, sing_s7[40] = 0.0, sing_s8[40] = 1.0;
    sing_s1[41] = 0.0, sing_s2[41] = 1.0, sing_s3[41] = 1.0, sing_s4[41] = 0.0, sing_s5[41] = 1.0, sing_s6[41] = 0.0, sing_s7[41] = 0.0, sing_s8[41] = 1.0;
    sing_s1[42] = 1.0, sing_s2[42] = 0.0, sing_s3[42] = 0.0, sing_s4[42] = 1.0, sing_s5[42] = 1.0, sing_s6[42] = 0.0, sing_s7[42] = 0.0, sing_s8[42] = 1.0;
    sing_s1[43] = 0.0, sing_s2[43] = 1.0, sing_s3[43] = 0.0, sing_s4[43] = 1.0, sing_s5[43] = 1.0, sing_s6[43] = 0.0, sing_s7[43] = 0.0, sing_s8[43] = 1.0;
    sing_s1[44] = 0.0, sing_s2[44] = 0.0, sing_s3[44] = 1.0, sing_s4[44] = 1.0, sing_s5[44] = 1.0, sing_s6[44] = 0.0, sing_s7[44] = 0.0, sing_s8[44] = 1.0;
    sing_s1[45] = 1.0, sing_s2[45] = 1.0, sing_s3[45] = 0.0, sing_s4[45] = 0.0, sing_s5[45] = 0.0, sing_s6[45] = 1.0, sing_s7[45] = 0.0, sing_s8[45] = 1.0;
    sing_s1[46] = 1.0, sing_s2[46] = 0.0, sing_s3[46] = 1.0, sing_s4[46] = 0.0, sing_s5[46] = 0.0, sing_s6[46] = 1.0, sing_s7[46] = 0.0, sing_s8[46] = 1.0;
    sing_s1[47] = 0.0, sing_s2[47] = 1.0, sing_s3[47] = 1.0, sing_s4[47] = 0.0, sing_s5[47] = 0.0, sing_s6[47] = 1.0, sing_s7[47] = 0.0, sing_s8[47] = 1.0;
    sing_s1[48] = 1.0, sing_s2[48] = 0.0, sing_s3[48] = 0.0, sing_s4[48] = 1.0, sing_s5[48] = 0.0, sing_s6[48] = 1.0, sing_s7[48] = 0.0, sing_s8[48] = 1.0;
    sing_s1[49] = 0.0, sing_s2[49] = 1.0, sing_s3[49] = 0.0, sing_s4[49] = 1.0, sing_s5[49] = 0.0, sing_s6[49] = 1.0, sing_s7[49] = 0.0, sing_s8[49] = 1.0;
    sing_s1[50] = 0.0, sing_s2[50] = 0.0, sing_s3[50] = 1.0, sing_s4[50] = 1.0, sing_s5[50] = 0.0, sing_s6[50] = 1.0, sing_s7[50] = 0.0, sing_s8[50] = 1.0;
    sing_s1[51] = 1.0, sing_s2[51] = 0.0, sing_s3[51] = 0.0, sing_s4[51] = 0.0, sing_s5[51] = 1.0, sing_s6[51] = 1.0, sing_s7[51] = 0.0, sing_s8[51] = 1.0;
    sing_s1[52] = 0.0, sing_s2[52] = 1.0, sing_s3[52] = 0.0, sing_s4[52] = 0.0, sing_s5[52] = 1.0, sing_s6[52] = 1.0, sing_s7[52] = 0.0, sing_s8[52] = 1.0;
    sing_s1[53] = 0.0, sing_s2[53] = 0.0, sing_s3[53] = 1.0, sing_s4[53] = 0.0, sing_s5[53] = 1.0, sing_s6[53] = 1.0, sing_s7[53] = 0.0, sing_s8[53] = 1.0;
    sing_s1[54] = 0.0, sing_s2[54] = 0.0, sing_s3[54] = 0.0, sing_s4[54] = 1.0, sing_s5[54] = 1.0, sing_s6[54] = 1.0, sing_s7[54] = 0.0, sing_s8[54] = 1.0;
    sing_s1[55] = 1.0, sing_s2[55] = 1.0, sing_s3[55] = 0.0, sing_s4[55] = 0.0, sing_s5[55] = 0.0, sing_s6[55] = 0.0, sing_s7[55] = 1.0, sing_s8[55] = 1.0;
    sing_s1[56] = 1.0, sing_s2[56] = 0.0, sing_s3[56] = 1.0, sing_s4[56] = 0.0, sing_s5[56] = 0.0, sing_s6[56] = 0.0, sing_s7[56] = 1.0, sing_s8[56] = 1.0;
    sing_s1[57] = 0.0, sing_s2[57] = 1.0, sing_s3[57] = 1.0, sing_s4[57] = 0.0, sing_s5[57] = 0.0, sing_s6[57] = 0.0, sing_s7[57] = 1.0, sing_s8[57] = 1.0;
    sing_s1[58] = 1.0, sing_s2[58] = 0.0, sing_s3[58] = 0.0, sing_s4[58] = 1.0, sing_s5[58] = 0.0, sing_s6[58] = 0.0, sing_s7[58] = 1.0, sing_s8[58] = 1.0;
    sing_s1[59] = 0.0, sing_s2[59] = 1.0, sing_s3[59] = 0.0, sing_s4[59] = 1.0, sing_s5[59] = 0.0, sing_s6[59] = 0.0, sing_s7[59] = 1.0, sing_s8[59] = 1.0;
    sing_s1[60] = 0.0, sing_s2[60] = 0.0, sing_s3[60] = 1.0, sing_s4[60] = 1.0, sing_s5[60] = 0.0, sing_s6[60] = 0.0, sing_s7[60] = 1.0, sing_s8[60] = 1.0;
    sing_s1[61] = 1.0, sing_s2[61] = 0.0, sing_s3[61] = 0.0, sing_s4[61] = 0.0, sing_s5[61] = 1.0, sing_s6[61] = 0.0, sing_s7[61] = 1.0, sing_s8[61] = 1.0;
    sing_s1[62] = 0.0, sing_s2[62] = 1.0, sing_s3[62] = 0.0, sing_s4[62] = 0.0, sing_s5[62] = 1.0, sing_s6[62] = 0.0, sing_s7[62] = 1.0, sing_s8[62] = 1.0;
    sing_s1[63] = 0.0, sing_s2[63] = 0.0, sing_s3[63] = 1.0, sing_s4[63] = 0.0, sing_s5[63] = 1.0, sing_s6[63] = 0.0, sing_s7[63] = 1.0, sing_s8[63] = 1.0;
    sing_s1[64] = 0.0, sing_s2[64] = 0.0, sing_s3[64] = 0.0, sing_s4[64] = 1.0, sing_s5[64] = 1.0, sing_s6[64] = 0.0, sing_s7[64] = 1.0, sing_s8[64] = 1.0;
    sing_s1[65] = 1.0, sing_s2[65] = 0.0, sing_s3[65] = 0.0, sing_s4[65] = 0.0, sing_s5[65] = 0.0, sing_s6[65] = 1.0, sing_s7[65] = 1.0, sing_s8[65] = 1.0;
    sing_s1[66] = 0.0, sing_s2[66] = 1.0, sing_s3[66] = 0.0, sing_s4[66] = 0.0, sing_s5[66] = 0.0, sing_s6[66] = 1.0, sing_s7[66] = 1.0, sing_s8[66] = 1.0;
    sing_s1[67] = 0.0, sing_s2[67] = 0.0, sing_s3[67] = 1.0, sing_s4[67] = 0.0, sing_s5[67] = 0.0, sing_s6[67] = 1.0, sing_s7[67] = 1.0, sing_s8[67] = 1.0;
    sing_s1[68] = 0.0, sing_s2[68] = 0.0, sing_s3[68] = 0.0, sing_s4[68] = 1.0, sing_s5[68] = 0.0, sing_s6[68] = 1.0, sing_s7[68] = 1.0, sing_s8[68] = 1.0;
    sing_s1[69] = 0.0, sing_s2[69] = 0.0, sing_s3[69] = 0.0, sing_s4[69] = 0.0, sing_s5[69] = 1.0, sing_s6[69] = 1.0, sing_s7[69] = 1.0, sing_s8[69] = 1.0;


// 


    for ( totn = 0; totn < m; totn++ ) {
        H_KA[totn] = 0.0;
    }


    for ( dian = 0; dian < size2; dian++ ) {
        H_KA[ ( size2 + 1 )*( dian ) ] = 2.0 * u_rep;
    }


    for ( sn = 0; sn < NCn; sn++ ) {
        H_KA[ ( size2 + 1 )*( NCn + 1 )*( sn ) ] = 4.0 * u_rep;
    }


    for ( sn = 0; sn < NCn; sn++ ) {
        H_KA[ ( size2 + 1 )*( NCn - 1 )*( sn + 1 ) ] = 0.0;
    }


    for ( sn = 0; sn < 1; sn++ ) {

        H_KA[ ( size2 + 1 )*( 31 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 31 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 32 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 32 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 33 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 33 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 34 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 34 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*1 + 28 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 29 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 30 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 34 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*1 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*2 + 26 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 27 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 30 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 33 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*2 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*3 + 25 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 27 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 29 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 32 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*3 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*4 + 25 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 26 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 28 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 31 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*4 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*5 + 22 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 23 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 24 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 34 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*5 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*6 + 20 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 21 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 24 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 33 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*6 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*7 + 19 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 21 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 23 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 32 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*7 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*8 + 19 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 20 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 22 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 31 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*8 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*9 + 17 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 18 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 24 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 30 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*9 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*10 + 16 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 18 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 23 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 29 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*10 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*11 + 16 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 17 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 22 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 28 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*11 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*12 + 15 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 18 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 21 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 27 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*12 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*13 + 15 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 17 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 20 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 26 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*13 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*14 + 15 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 16 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 19 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 25 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*14 + NCn - 1 - 66 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*15 + 12 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 13 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 14 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 34 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 34 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 64 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 68 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*15 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*16 + 10 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 11 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 14 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 33 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 33 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*16 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*17 + 9 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 11 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 13 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 32 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 32 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*17 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*18 + 9 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 10 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 12 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 31 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 31 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 69 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*18 + NCn - 1 - 69 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*19 + 7 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 8 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 14 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 30 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 30 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 67 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*19 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*20 + 6 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 8 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 13 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 29 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 29 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*20 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*20 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*21 + 6 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 7 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 12 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 28 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 28 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*21 + 68 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*21 + NCn - 1 - 68 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*22 + 5 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 8 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 11 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 27 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 27 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 66 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*22 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*22 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*23 + 5 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 7 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 10 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 26 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 26 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*23 + 67 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*23 + NCn - 1 - 67 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*24 + 5 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 6 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 9 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 25 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 25 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 65 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 65 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*24 + 66 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*24 + NCn - 1 - 66 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*25 + 3 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 4 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 14 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 14 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 24 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 24 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 60 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 63 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*25 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*25 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*26 + 2 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 4 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 13 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 13 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 23 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 23 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*26 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*26 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*27 + 2 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 3 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 12 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 12 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 22 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 22 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 54 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 54 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*27 + 64 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*27 + NCn - 1 - 64 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*28 + 1 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 4 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 11 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 11 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 21 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 21 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 62 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*28 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*28 + NCn - 1 - 63 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*29 + 1 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 3 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 10 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 10 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 20 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 20 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 53 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 53 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*29 + 63 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*29 + NCn - 1 - 63 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*30 + 1 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 2 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 9 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 9 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 19 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 19 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 51 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 51 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 52 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 52 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 61 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 61 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*30 + 62 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*30 + NCn - 1 - 62 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*31 + 0 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 4 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 4 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 8 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 8 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 18 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 18 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 57 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 59 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*31 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*31 + NCn - 1 - 60 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*32 + 0 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 3 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 3 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 7 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 7 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 17 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 17 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 44 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 44 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 50 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 50 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*32 + 60 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*32 + NCn - 1 - 60 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*33 + 0 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 2 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 2 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 6 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 6 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 16 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 16 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 35 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 35 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 42 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 42 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 43 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 43 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 48 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 48 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 49 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 49 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 58 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 58 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*33 + 59 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*33 + NCn - 1 - 59 ) ] = 3.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*34 + 0 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 0 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 1 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 1 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 5 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 5 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 15 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 15 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 36 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 36 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 37 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 37 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 38 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 38 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 39 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 39 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 40 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 40 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 41 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 41 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 45 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 45 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 46 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 46 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 47 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 47 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 55 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 55 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 56 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 56 ) ] = 3.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*34 + 57 ) ] = u_rep, H_KA[ ( size2 + 1 )*( NCn*34 + NCn - 1 - 57 ) ] = 3.0 * u_rep;

// 

        H_KA[ ( size2 + 1 )*( NCn*35 + 0 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 0 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 1 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 1 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 5 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 5 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 15 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 15 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 36 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 37 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 38 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 39 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 40 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 41 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 45 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 46 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 47 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 55 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 56 ) ] = u_rep;
        H_KA[ ( size2 + 1 )*( NCn*35 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*35 + NCn - 1 - 57 ) ] = u_rep;

        H_KA[ ( size2 + 1 )*( NCn*36 + 0 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 2 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 6 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 16 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*36 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*36 + NCn - 1 - 59 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*37 + 0 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 3 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 7 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 17 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*37 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*37 + NCn - 1 - 60 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*38 + 0 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 0 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 4 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 8 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 18 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*38 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*38 + NCn - 1 - 60 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*39 + 1 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 2 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 9 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 19 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*39 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*39 + NCn - 1 - 62 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*40 + 1 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 3 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 10 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 20 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*40 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*40 + NCn - 1 - 63 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*41 + 1 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 1 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 4 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 11 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 21 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*41 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*41 + NCn - 1 - 63 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*42 + 2 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 3 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 12 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 22 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*42 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*42 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*43 + 2 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 2 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 4 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 13 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 23 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*43 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*43 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*44 + 3 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 3 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 4 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 4 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 14 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 24 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*44 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*44 + NCn - 1 - 64 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*45 + 5 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 6 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 9 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 25 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*45 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*45 + NCn - 1 - 66 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*46 + 5 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 7 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 10 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 26 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*46 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*46 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*47 + 5 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 5 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 8 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 11 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 27 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*47 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*47 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*48 + 6 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 7 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 12 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 28 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*48 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*48 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*49 + 6 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 6 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 8 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 13 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 29 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*49 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*49 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*50 + 7 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 7 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 8 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 8 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 14 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 30 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*50 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*50 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*51 + 9 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 10 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 12 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 31 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*51 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*51 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*52 + 9 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 9 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 11 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 13 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 32 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*52 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*52 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*53 + 10 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 10 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 11 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 11 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 14 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 33 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*53 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*53 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*54 + 12 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 12 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 13 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 13 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 14 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 14 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 34 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*54 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*54 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*55 + 15 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 16 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 19 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 25 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*55 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*55 + NCn - 1 - 66 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*56 + 15 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 17 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 20 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 26 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*56 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*56 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*57 + 15 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 15 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 18 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 21 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 27 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 35 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 35 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*57 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*57 + NCn - 1 - 67 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*58 + 16 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 17 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 22 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 28 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*58 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*58 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*59 + 16 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 16 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 18 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 23 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 29 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 36 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 36 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*59 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*59 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*60 + 17 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 17 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 18 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 18 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 24 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 30 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 37 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 37 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 38 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 38 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*60 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*60 + NCn - 1 - 68 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*61 + 19 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 20 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 22 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 31 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*61 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*61 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*62 + 19 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 19 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 21 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 23 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 32 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 39 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 39 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*62 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*62 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*63 + 20 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 20 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 21 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 21 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 24 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 33 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 40 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 40 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 41 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 41 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*63 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*63 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*64 + 22 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 22 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 23 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 23 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 24 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 24 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 34 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 42 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 42 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 43 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 43 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 44 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 44 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*64 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*64 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*65 + 25 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 26 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 28 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 31 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*65 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*65 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*66 + 25 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 25 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 27 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 29 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 32 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 45 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 45 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 55 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 55 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*66 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*66 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*67 + 26 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 26 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 27 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 27 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 30 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 33 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 46 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 46 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 47 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 47 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 56 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 56 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 57 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 57 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 68 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*67 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*67 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*68 + 28 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 28 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 29 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 29 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 30 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 30 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 34 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 48 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 48 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 49 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 49 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 50 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 50 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 58 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 58 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 59 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 59 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 60 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 60 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*68 + 69 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*68 + NCn - 1 - 69 ) ] = 1.0 * u_rep;

        H_KA[ ( size2 + 1 )*( NCn*69 + 31 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 31 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 32 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 32 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 33 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 33 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 34 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 34 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 51 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 51 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 52 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 52 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 53 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 53 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 54 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 54 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 61 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 61 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 62 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 62 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 63 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 63 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 64 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 64 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 65 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 65 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 66 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 66 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 67 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 67 ) ] = 1.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*69 + 68 ) ] = 3.0 * u_rep, H_KA[ ( size2 + 1 )*( NCn*69 + NCn - 1 - 68 ) ] = 1.0 * u_rep;


    }

// 




    for ( hopi = 0; hopi < NCn; hopi++ ) {


        H_KA[ size2*(NCn*hopi) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi) + 38 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+1) + 0 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+1) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+1) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+1) + 41 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+2) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+2) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+2) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+2) + 43 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+3) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+3) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+3) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+3) + 44 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+4) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+4) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+5) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+5) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+5) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+5) + 47 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+6) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+6) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+6) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+6) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+6) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+6) + 49 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+7) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+7) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+7) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+7) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+7) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+7) + 50 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+8) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+8) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+8) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+8) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+9) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+9) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+9) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+9) + 52 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+10) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+10) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+10) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+10) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+10) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+10) + 53 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+11) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+11) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+11) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+11) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+12) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+12) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+12) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+12) + 54 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+13) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+13) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+13) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+13) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+14) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+14) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+15) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+15) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+15) + 35 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+15) + 57 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+16) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+16) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+16) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+16) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+16) + 36 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+16) + 59 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+17) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+17) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+17) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+17) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+17) + 37 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+17) + 60 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+18) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+18) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+18) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+18) + 38 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+19) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+19) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+19) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+19) + 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+19) + 39 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+19) + 62 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+20) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 40 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+20) + 63 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+21) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+21) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+21) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+21) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+21) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+21) + 41 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+22) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+22) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+22) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+22) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+22) + 42 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+22) + 64 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+23) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+23) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+23) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+23) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+23) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+23) + 43 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+24) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+24) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+24) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+24) + 44 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+25) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+25) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+25) + 45 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+25) + 66 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+26) + 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+26) + 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+26) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+26) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+26) + 46 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+26) + 67 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+27) + 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+27) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+27) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+27) + 47 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+28) + 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+28) + 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+28) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+28) + 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+28) + 48 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+28) + 68 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+29) + 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+29) + 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+29) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+29) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+29) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+29) + 49 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+30) + 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+30) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+30) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+30) + 50 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+31) + 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+31) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+31) + 51 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+31) + 69 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+32) + 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+32) + 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+32) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+32) + 52 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+33) + 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+33) + 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+33) + 34 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+33) + 53 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+34) + 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+34) + 54 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

// 

        H_KA[ size2*(NCn*hopi+35) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+35) + NCn - 1 - 54 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+36) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+36) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+36) + NCn - 1 - 34 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+36) + NCn - 1 - 53 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+37) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+37) + NCn - 1 - 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+37) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+37) + NCn - 1 - 52 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+38) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+38) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+38) + NCn - 1 - 51 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+38) + NCn - 1 - 69 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+39) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+39) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+39) + NCn - 1 - 33 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+39) + NCn - 1 - 50 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 32 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+40) + NCn - 1 - 49 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 31 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 48 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+41) + NCn - 1 - 68 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+42) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+42) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+42) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+42) + NCn - 1 - 47 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 46 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+43) + NCn - 1 - 67 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+44) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+44) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+44) + NCn - 1 - 45 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+44) + NCn - 1 - 66 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+45) + NCn - 1 - 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+45) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+45) + NCn - 1 - 30 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+45) + NCn - 1 - 44 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 29 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+46) + NCn - 1 - 43 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 28 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 42 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+47) + NCn - 1 - 64 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 27 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+48) + NCn - 1 - 41 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 26 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 40 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+49) + NCn - 1 - 63 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 25 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 39 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+50) + NCn - 1 - 62 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+51) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+51) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+51) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+51) + NCn - 1 - 38 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 37 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+52) + NCn - 1 - 60 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 36 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+53) + NCn - 1 - 59 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+54) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+54) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+54) + NCn - 1 - 35 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s8[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+54) + NCn - 1 - 57 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+55) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+55) + NCn - 1 - 24 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+56) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+56) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+56) + NCn - 1 - 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+56) + NCn - 1 - 23 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+57) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+57) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+57) + NCn - 1 - 22 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+57) + NCn - 1 - 54 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+58) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+58) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+58) + NCn - 1 - 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+58) + NCn - 1 - 21 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 20 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+59) + NCn - 1 - 53 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+60) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+60) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+60) + NCn - 1 - 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+60) + NCn - 1 - 52 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+61) + NCn - 1 - 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+61) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+61) + NCn - 1 - 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+61) + NCn - 1 - 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+62) + NCn - 1 - 50 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+63) + NCn - 1 - 49 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+64) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+64) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+64) + NCn - 1 - 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+64) + NCn - 1 - 47 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+65) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+65) + NCn - 1 - 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+66) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+66) + NCn - 1 - 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+66) + NCn - 1 - 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+66) + NCn - 1 - 44 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+67) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+67) + NCn - 1 - 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+67) + NCn - 1 - 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+67) + NCn - 1 - 43 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+68) + NCn - 1 - 0 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+68) + NCn - 1 - 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+68) + NCn - 1 - 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+68) + NCn - 1 - 41 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*(NCn*hopi+69) + NCn - 1 - 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*(NCn*hopi+69) + NCn - 1 - 38 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] + sing_s8[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );


    }

// 

    for ( hopi = 0; hopi < NCn; hopi++ ) {

        H_KA[ size2*hopi + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*hopi + 38*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s6[0] + sing_s7[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*1 ) + 0*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*1 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*1 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*1 ) + 41*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s6[1] + sing_s7[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*2 ) + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*2 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*2 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*2 ) + 43*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s6[2] + sing_s7[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*3 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*3 ) + 4*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*3 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*3 ) + 44*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[3] + sing_s3[3] + sing_s4[3] + sing_s5[3] + sing_s6[3] + sing_s7[3] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*4 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*4 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*5 ) + 1*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*5 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*5 ) + 15*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*5 ) + 47*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s6[5] + sing_s7[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*6 ) + 2*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*6 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*6 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*6 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*6 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*6 ) + 49*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[6] + sing_s3[6] + sing_s4[6] + sing_s5[6] + sing_s6[6] + sing_s7[6] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*7 ) + 3*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*7 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*7 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*7 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*7 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*7 ) + 50*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s6[7] + sing_s7[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*8 ) + 4*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*8 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*8 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*8 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*9 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*9 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*9 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*9 ) + 52*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[9] + sing_s3[9] + sing_s4[9] + sing_s5[9] + sing_s6[9] + sing_s7[9] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*10 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*10 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*10 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*10 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*10 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*10 ) + 53*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[10] + sing_s3[10] + sing_s4[10] + sing_s5[10] + sing_s6[10] + sing_s7[10] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*11 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*11 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*11 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*11 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*12 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*12 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*12 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*12 ) + 54*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s6[12] + sing_s7[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*13 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*13 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*13 ) + 14*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*13 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*14 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*14 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*15 ) + 5*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*15 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*15 ) + 35*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*15 ) + 57*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s6[15] + sing_s7[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*16 ) + 6*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*16 ) + 15*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*16 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*16 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*16 ) + 36*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*16 ) + 59*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[16] + sing_s3[16] + sing_s4[16] + sing_s5[16] + sing_s6[16] + sing_s7[16] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*17 ) + 7*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*17 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*17 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*17 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*17 ) + 37*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*17 ) + 60*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s6[17] + sing_s7[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*18 ) + 8*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*18 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*18 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*18 ) + 38*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*19 ) + 9*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*19 ) + 16*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*19 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*19 ) + 25*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*19 ) + 39*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*19 ) + 62*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s6[19] + sing_s7[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*20 ) + 10*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 17*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 40*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*20 ) + 63*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[20] + sing_s3[20] + sing_s4[20] + sing_s5[20] + sing_s6[20] + sing_s7[20] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*21 ) + 11*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*21 ) + 18*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*21 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*21 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*21 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*21 ) + 41*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*22 ) + 12*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*22 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*22 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*22 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*22 ) + 42*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*22 ) + 64*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[22] + sing_s3[22] + sing_s4[22] + sing_s5[22] + sing_s6[22] + sing_s7[22] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*23 ) + 13*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*23 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*23 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*23 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*23 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*23 ) + 43*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*24 ) + 14*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*24 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*24 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*24 ) + 44*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*25 ) + 19*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*25 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*25 ) + 45*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*25 ) + 66*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[25] + sing_s3[25] + sing_s4[25] + sing_s5[25] + sing_s6[25] + sing_s7[25] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*26 ) + 20*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*26 ) + 25*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*26 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*26 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*26 ) + 46*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*26 ) + 67*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[26] + sing_s3[26] + sing_s4[26] + sing_s5[26] + sing_s6[26] + sing_s7[26] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*27 ) + 21*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*27 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*27 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*27 ) + 47*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*28 ) + 22*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*28 ) + 26*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*28 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*28 ) + 31*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*28 ) + 48*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*28 ) + 68*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[28] + sing_s3[28] + sing_s4[28] + sing_s5[28] + sing_s6[28] + sing_s7[28] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*29 ) + 23*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*29 ) + 27*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*29 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*29 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*29 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*29 ) + 49*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*30 ) + 24*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*30 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*30 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*30 ) + 50*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*31 ) + 28*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*31 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*31 ) + 51*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*31 ) + 69*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[31] + sing_s3[31] + sing_s4[31] + sing_s5[31] + sing_s6[31] + sing_s7[31] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*32 ) + 29*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*32 ) + 31*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*32 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*32 ) + 52*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*33 ) + 30*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*33 ) + 32*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*33 ) + 34*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*33 ) + 53*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*34 ) + 33*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*34 ) + 54*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

// 

        H_KA[ size2*( hopi + NCn*35 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*35 ) + ( NCn - 1 - 54 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*36 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*36 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*36 ) + ( NCn - 1 - 34 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*36 ) + ( NCn - 1 - 53 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*37 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*37 ) + ( NCn - 1 - 31 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*37 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*37 ) + ( NCn - 1 - 52 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*38 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*38 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*38 ) + ( NCn - 1 - 51 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*38 ) + ( NCn - 1 - 69 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[38] + sing_s3[38] + sing_s4[38] + sing_s5[38] + sing_s6[38] + sing_s7[38] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*39 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*39 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*39 ) + ( NCn - 1 - 33 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*39 ) + ( NCn - 1 - 50 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 32 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*40 ) + ( NCn - 1 - 49 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 31 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 48 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*41 ) + ( NCn - 1 - 68 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[41] + sing_s3[41] + sing_s4[41] + sing_s5[41] + sing_s6[41] + sing_s7[41] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*42 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*42 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*42 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*42 ) + ( NCn - 1 - 47 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 25 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 46 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*43 ) + ( NCn - 1 - 67 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[43] + sing_s3[43] + sing_s4[43] + sing_s5[43] + sing_s6[43] + sing_s7[43] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*44 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*44 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*44 ) + ( NCn - 1 - 45 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*44 ) + ( NCn - 1 - 66 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[44] + sing_s3[44] + sing_s4[44] + sing_s5[44] + sing_s6[44] + sing_s7[44] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*45 ) + ( NCn - 1 - 14 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*45 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*45 ) + ( NCn - 1 - 30 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*45 ) + ( NCn - 1 - 44 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 29 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*46 ) + ( NCn - 1 - 43 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 28 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 42 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*47 ) + ( NCn - 1 - 64 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[47] + sing_s3[47] + sing_s4[47] + sing_s5[47] + sing_s6[47] + sing_s7[47] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 27 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*48 ) + ( NCn - 1 - 41 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 26 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 40 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*49 ) + ( NCn - 1 - 63 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[49] + sing_s3[49] + sing_s4[49] + sing_s5[49] + sing_s6[49] + sing_s7[49] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 25 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 39 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*50 ) + ( NCn - 1 - 62 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[50] + sing_s3[50] + sing_s4[50] + sing_s5[50] + sing_s6[50] + sing_s7[50] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*51 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*51 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*51 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*51 ) + ( NCn - 1 - 38 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 37 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*52 ) + ( NCn - 1 - 60 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[52] + sing_s3[52] + sing_s4[52] + sing_s5[52] + sing_s6[52] + sing_s7[52] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 15 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 36 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*53 ) + ( NCn - 1 - 59 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[53] + sing_s3[53] + sing_s4[53] + sing_s5[53] + sing_s6[53] + sing_s7[53] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*54 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*54 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*54 ) + ( NCn - 1 - 35 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s7[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*54 ) + ( NCn - 1 - 57 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[54] + sing_s3[54] + sing_s4[54] + sing_s5[54] + sing_s6[54] + sing_s7[54] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*55 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*55 ) + ( NCn - 1 - 24 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*56 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*56 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*56 ) + ( NCn - 1 - 14 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*56 ) + ( NCn - 1 - 23 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*57 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*57 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*57 ) + ( NCn - 1 - 22 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*57 ) + ( NCn - 1 - 54 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[57] + sing_s3[57] + sing_s4[57] + sing_s5[57] + sing_s6[57] + sing_s7[57] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*58 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*58 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*58 ) + ( NCn - 1 - 13 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*58 ) + ( NCn - 1 - 21 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 12 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 20 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*59 ) + ( NCn - 1 - 53 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[59] + sing_s3[59] + sing_s4[59] + sing_s5[59] + sing_s6[59] + sing_s7[59] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*60 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*60 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*60 ) + ( NCn - 1 - 19 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*60 ) + ( NCn - 1 - 52 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[60] + sing_s3[60] + sing_s4[60] + sing_s5[60] + sing_s6[60] + sing_s7[60] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*61 ) + ( NCn - 1 - 4 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*61 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*61 ) + ( NCn - 1 - 11 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*61 ) + ( NCn - 1 - 18 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 10 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 17 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*62 ) + ( NCn - 1 - 50 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[62] + sing_s3[62] + sing_s4[62] + sing_s5[62] + sing_s6[62] + sing_s7[62] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 9 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 16 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*63 ) + ( NCn - 1 - 49 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[63] + sing_s3[63] + sing_s4[63] + sing_s5[63] + sing_s6[63] + sing_s7[63] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*64 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*64 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*64 ) + ( NCn - 1 - 15 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*64 ) + ( NCn - 1 - 47 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[64] + sing_s3[64] + sing_s4[64] + sing_s5[64] + sing_s6[64] + sing_s7[64] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*65 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*65 ) + ( NCn - 1 - 8 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*66 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*66 ) + ( NCn - 1 - 4 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*66 ) + ( NCn - 1 - 7 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*66 ) + ( NCn - 1 - 44 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[66] + sing_s3[66] + sing_s4[66] + sing_s5[66] + sing_s6[66] + sing_s7[66] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*67 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*67 ) + ( NCn - 1 - 3 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*67 ) + ( NCn - 1 - 6 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*67 ) + ( NCn - 1 - 43 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[67] + sing_s3[67] + sing_s4[67] + sing_s5[67] + sing_s6[67] + sing_s7[67] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*68 ) + ( NCn - 1 - 0 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*68 ) + ( NCn - 1 - 2 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*68 ) + ( NCn - 1 - 5 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*68 ) + ( NCn - 1 - 41 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[68] + sing_s3[68] + sing_s4[68] + sing_s5[68] + sing_s6[68] + sing_s7[68] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );

        H_KA[ size2*( hopi + NCn*69 ) + ( NCn - 1 - 1 )*NCn + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] ) * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) );
        H_KA[ size2*( hopi + NCn*69 ) + ( NCn - 1 - 38 )*NCn + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[69] + sing_s3[69] + sing_s4[69] + sing_s5[69] + sing_s6[69] + sing_s7[69] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] + sing_s7[hopi] ) ) * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) );


    }

// 


    gsl_matrix_complex *H_T = NULL;
    gsl_matrix_complex *RHO_T = NULL;
    gsl_matrix_complex *COMM1 = NULL;
    gsl_matrix_complex *COMM2 = NULL;

    gsl_matrix_complex *D = NULL;

    H_T = gsl_matrix_complex_alloc(size2, size2);
    RHO_T = gsl_matrix_complex_alloc(size2, size2);
    COMM1 = gsl_matrix_complex_alloc(size2, size2);
    COMM2 = gsl_matrix_complex_alloc(size2, size2);

    D = gsl_matrix_complex_alloc(size2, size2);

    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {
            gsl_matrix_complex_set(H_T, rhi, rhj, gsl_complex_rect(creal(H_KA[rhi * size2 + rhj]),cimag(H_KA[rhi * size2 + rhj])));
        }
    }

    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {
            gsl_matrix_complex_set(RHO_T, rhi, rhj, gsl_complex_rect(creal(u[rhi * size2 + rhj]),cimag(u[rhi * size2 + rhj])));
        }
    }

    gsl_matrix_complex_set_all(COMM1, GSL_COMPLEX_ZERO);
    gsl_matrix_complex_set_all(COMM2, GSL_COMPLEX_ZERO);

    gsl_matrix_complex_set_all(D, GSL_COMPLEX_ZERO);



    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, H_T, RHO_T,
                    GSL_COMPLEX_ZERO, COMM1);


    for ( rhi = 0 ; rhi < NCn ; rhi++ ) {

        for ( rhj = 0 ; rhj < NCn ; rhj++ ) {

            site1[ rhi * NCn + rhj ] = sing_s1[ rhi ] + sing_s1[ rhj ];
            site2[ rhi * NCn + rhj ] = sing_s2[ rhi ] + sing_s2[ rhj ];
            site3[ rhi * NCn + rhj ] = sing_s3[ rhi ] + sing_s3[ rhj ];
            site4[ rhi * NCn + rhj ] = sing_s4[ rhi ] + sing_s4[ rhj ];
            site5[ rhi * NCn + rhj ] = sing_s5[ rhi ] + sing_s5[ rhj ];
            site6[ rhi * NCn + rhj ] = sing_s6[ rhi ] + sing_s6[ rhj ];
            site7[ rhi * NCn + rhj ] = sing_s7[ rhi ] + sing_s7[ rhj ];
            site8[ rhi * NCn + rhj ] = sing_s8[ rhi ] + sing_s8[ rhj ];

        }

    }

//

    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {

            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site1[ rhi ] * site1[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site2[ rhi ] * site2[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site3[ rhi ] * site3[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site4[ rhi ] * site4[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site5[ rhi ] * site5[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site6[ rhi ] * site6[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site7[ rhi ] * site7[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] + site8[ rhi ] * site8[ rhj ];

        }
    }

// 

    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {

            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site1[ rhi ] * site1[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site2[ rhi ] * site2[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site3[ rhi ] * site3[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site4[ rhi ] * site4[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site5[ rhi ] * site5[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site6[ rhi ] * site6[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site7[ rhi ] * site7[ rhi ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site8[ rhi ] * site8[ rhi ];

        }
    }

// 

    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {

            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site1[ rhj ] * site1[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site2[ rhj ] * site2[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site3[ rhj ] * site3[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site4[ rhj ] * site4[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site5[ rhj ] * site5[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site6[ rhj ] * site6[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site7[ rhj ] * site7[ rhj ];
            sub_D[ rhi * size2 + rhj ] = sub_D[ rhi * size2 + rhj ] - 0.5 * site8[ rhj ] * site8[ rhj ];

        }
    }

// 


    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {

            gsl_matrix_complex_set(D, rhi, rhj, gsl_complex_rect( sub_D[ rhi * size2 + rhj ] * GSL_REAL(gsl_matrix_complex_get(RHO_T,rhi,rhj)) , sub_D[ rhi * size2 + rhj ] * GSL_IMAG(gsl_matrix_complex_get(RHO_T,rhi,rhj)) ));

        }
    }


// main eq. 

    for ( oparai = 0 ; oparai < size2 ; oparai++ ) {
        for ( oparaj = 0 ; oparaj < size2 ; oparaj++ ) {
            arr[oparai * size2 + oparaj] = -1.0*I * ( GSL_REAL(gsl_matrix_complex_get(COMM1,oparai,oparaj)) + GSL_IMAG(gsl_matrix_complex_get(COMM1,oparai,oparaj))*I - GSL_REAL(gsl_matrix_complex_get(COMM1,oparaj,oparai)) + GSL_IMAG(gsl_matrix_complex_get(COMM1,oparaj,oparai))*I )
                                                + gamma * ( GSL_REAL(gsl_matrix_complex_get(D,oparai,oparaj)) + GSL_IMAG(gsl_matrix_complex_get(D,oparai,oparaj))*I );
        }
    }


    free(sub_D);

    gsl_matrix_complex_free(H_T);
    gsl_matrix_complex_free(RHO_T);
    gsl_matrix_complex_free(COMM1);
    gsl_matrix_complex_free(COMM2);

    gsl_matrix_complex_free(D);

    return arr;
}
