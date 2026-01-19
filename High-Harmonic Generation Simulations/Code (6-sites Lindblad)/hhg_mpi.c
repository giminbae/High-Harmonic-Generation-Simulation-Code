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
This parameter is defined around line 880 of the source code.
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

    int m = 160000, j = 0, q = 0, rx, ry;// m = hamiltonian size
    int totn, dian, NCn = 20, sn, hopi;
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
    double Ex = ( 3.0/6.5 )*0.00126410*alpha/(sqrt(1.0+ellip*ellip));//optical field strength
    double wpump = 0.26/27.2114; //pump frequency
    double Tfull = 4134.2*2.0; //pulse length

    double u_rep = 2.0*0.52/27.2114;// Hubbard U
    double t_hub = 0.52/27.2114;// hopping t

    lapack_int n = 400;
    lapack_int lda = n;
    lapack_int info;
    double E_eigval[n];
    lapack_complex_double H_0 [n*n];

    int size = 400;
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

    double sing_s1[NCn], sing_s2[NCn], sing_s3[NCn], sing_s4[NCn], sing_s5[NCn], sing_s6[NCn];

    sing_s1[0] = 1.0,  sing_s2[0] = 1.0,  sing_s3[0] = 1.0,  sing_s4[0] = 0.0,  sing_s5[0] = 0.0,  sing_s6[0] = 0.0;
    sing_s1[1] = 1.0,  sing_s2[1] = 1.0,  sing_s3[1] = 0.0,  sing_s4[1] = 1.0,  sing_s5[1] = 0.0,  sing_s6[1] = 0.0;
    sing_s1[2] = 1.0,  sing_s2[2] = 0.0,  sing_s3[2] = 1.0,  sing_s4[2] = 1.0,  sing_s5[2] = 0.0,  sing_s6[2] = 0.0;
    sing_s1[3] = 0.0,  sing_s2[3] = 1.0,  sing_s3[3] = 1.0,  sing_s4[3] = 1.0,  sing_s5[3] = 0.0,  sing_s6[3] = 0.0;
    sing_s1[4] = 1.0,  sing_s2[4] = 1.0,  sing_s3[4] = 0.0,  sing_s4[4] = 0.0,  sing_s5[4] = 1.0,  sing_s6[4] = 0.0;
    sing_s1[5] = 1.0,  sing_s2[5] = 0.0,  sing_s3[5] = 1.0,  sing_s4[5] = 0.0,  sing_s5[5] = 1.0,  sing_s6[5] = 0.0;
    sing_s1[6] = 0.0,  sing_s2[6] = 1.0,  sing_s3[6] = 1.0,  sing_s4[6] = 0.0,  sing_s5[6] = 1.0,  sing_s6[6] = 0.0;
    sing_s1[7] = 1.0,  sing_s2[7] = 0.0,  sing_s3[7] = 0.0,  sing_s4[7] = 1.0,  sing_s5[7] = 1.0,  sing_s6[7] = 0.0;
    sing_s1[8] = 0.0,  sing_s2[8] = 1.0,  sing_s3[8] = 0.0,  sing_s4[8] = 1.0,  sing_s5[8] = 1.0,  sing_s6[8] = 0.0;
    sing_s1[9] = 0.0,  sing_s2[9] = 0.0,  sing_s3[9] = 1.0,  sing_s4[9] = 1.0,  sing_s5[9] = 1.0,  sing_s6[9] = 0.0;
    sing_s1[10] = 1.0, sing_s2[10] = 1.0, sing_s3[10] = 0.0, sing_s4[10] = 0.0, sing_s5[10] = 0.0, sing_s6[10] = 1.0;
    sing_s1[11] = 1.0, sing_s2[11] = 0.0, sing_s3[11] = 1.0, sing_s4[11] = 0.0, sing_s5[11] = 0.0, sing_s6[11] = 1.0;
    sing_s1[12] = 0.0, sing_s2[12] = 1.0, sing_s3[12] = 1.0, sing_s4[12] = 0.0, sing_s5[12] = 0.0, sing_s6[12] = 1.0;
    sing_s1[13] = 1.0, sing_s2[13] = 0.0, sing_s3[13] = 0.0, sing_s4[13] = 1.0, sing_s5[13] = 0.0, sing_s6[13] = 1.0;
    sing_s1[14] = 0.0, sing_s2[14] = 1.0, sing_s3[14] = 0.0, sing_s4[14] = 1.0, sing_s5[14] = 0.0, sing_s6[14] = 1.0;
    sing_s1[15] = 0.0, sing_s2[15] = 0.0, sing_s3[15] = 1.0, sing_s4[15] = 1.0, sing_s5[15] = 0.0, sing_s6[15] = 1.0;
    sing_s1[16] = 1.0, sing_s2[16] = 0.0, sing_s3[16] = 0.0, sing_s4[16] = 0.0, sing_s5[16] = 1.0, sing_s6[16] = 1.0;
    sing_s1[17] = 0.0, sing_s2[17] = 1.0, sing_s3[17] = 0.0, sing_s4[17] = 0.0, sing_s5[17] = 1.0, sing_s6[17] = 1.0;
    sing_s1[18] = 0.0, sing_s2[18] = 0.0, sing_s3[18] = 1.0, sing_s4[18] = 0.0, sing_s5[18] = 1.0, sing_s6[18] = 1.0;
    sing_s1[19] = 0.0, sing_s2[19] = 0.0, sing_s3[19] = 0.0, sing_s4[19] = 1.0, sing_s5[19] = 1.0, sing_s6[19] = 1.0;

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
                H_0[ ( size + 1 )*( dian ) ] = u_rep;
            }


            for ( sn = 0; sn < NCn; sn++ ) {
                H_0[ ( size + 1 )*( NCn + 1 )*( sn ) ] = 3.0 * u_rep;
            }


            for ( sn = 0; sn < NCn; sn++ ) {
                H_0[ ( size + 1 )*( NCn - 1 )*( sn + 1 ) ] = 0.0;
            }

// 

            for ( sn = 0; sn < 1; sn++ ) {

                H_0[ ( size + 1 )*( 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( 12 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*1 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*1 + 14 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*2 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*2 + 15 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*3 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*3 + 15 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*4 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*4 + 17 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*5 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*5 + 18 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*6 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*6 + 18 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*7 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*7 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*8 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*8 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*9 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 18 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*9 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*10 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*10 + 17 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*11 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*11 + 18 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*12 + 0 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*12 + 18 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*13 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*13 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*14 + 1 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*14 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*15 + 2 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 3 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 18 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*15 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*16 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 18 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*16 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*17 + 4 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 10 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 18 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*17 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*18 + 5 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 6 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 11 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 12 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*18 + 19 ) ] = 2.0 * u_rep;


                H_0[ ( size + 1 )*( NCn*19 + 7 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 8 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 9 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 13 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 14 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 15 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 16 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 17 ) ] = 2.0 * u_rep;
                H_0[ ( size + 1 )*( NCn*19 + 18 ) ] = 2.0 * u_rep;
            }

//

            for ( hopi = 0; hopi < NCn; hopi++ ) {

                H_0[ size*(NCn*hopi) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+1) + 0 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+1) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+1) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+1) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+2) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+2) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+2) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+2) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+3) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+3) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*(NCn*hopi+4) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+4) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+4) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+4) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+5) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+5) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+5) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+5) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+5) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+5) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+6) + 3 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+6) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+6) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+6) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*(NCn*hopi+7) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+7) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+7) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+7) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

                H_0[ size*(NCn*hopi+8) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+8) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+8) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+8) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*(NCn*hopi+9) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+9) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );

                H_0[ size*(NCn*hopi+10) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+10) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

                H_0[ size*(NCn*hopi+11) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+11) + 10 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+11) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+11) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*(NCn*hopi+12) + 0 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+12) + 6 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+12) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+12) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*(NCn*hopi+13) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+13) + 11 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+13) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+13) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*(NCn*hopi+14) + 1 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+14) + 8 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+14) + 12 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*(NCn*hopi+14) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+14) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+14) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*(NCn*hopi+15) + 2 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+15) + 9 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s6[hopi] );
                H_0[ size*(NCn*hopi+15) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+15) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*(NCn*hopi+16) + 13 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+16) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );

                H_0[ size*(NCn*hopi+17) + 4 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+17) + 14 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+17) + 16 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*(NCn*hopi+17) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

                H_0[ size*(NCn*hopi+18) + 5 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+18) + 15 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*(NCn*hopi+18) + 17 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*(NCn*hopi+18) + 19 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*(NCn*hopi+19) + 7 + NCn*hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
                H_0[ size*(NCn*hopi+19) + 18 + NCn*hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

            }

// 

            for ( hopi = 0; hopi < NCn; hopi++ ) {

                H_0[ size*hopi + 20 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*hopi + 240 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*1 ) + 0 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 40 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 80 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*1 ) + 280 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*2 ) + 20 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 60 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 100 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*2 ) + 300 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*3 ) + 40 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*3 ) + 120 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*( hopi + NCn*4 ) + 20 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*4 ) + 100 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*4 ) + 200 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*4 ) + 340 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*5 ) + 40 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 80 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 120 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 140 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 220 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*5 ) + 360 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*6 ) + 60 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 100 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 160 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*6 ) + 240 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*( hopi + NCn*7 ) + 100 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 160 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 260 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*7 ) + 380 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

                H_0[ size*( hopi + NCn*8 ) + 120 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 140 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 180 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*8 ) + 280 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*( hopi + NCn*9 ) + 160 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*9 ) + 300 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );

                H_0[ size*( hopi + NCn*10 ) + 80 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*10 ) + 220 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );

                H_0[ size*( hopi + NCn*11 ) + 100 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 200 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 240 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*11 ) + 260 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

                H_0[ size*( hopi + NCn*12 ) + 0 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*12 ) + 120 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*12 ) + 220 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*12 ) + 280 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

                H_0[ size*( hopi + NCn*13 ) + 140 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 220 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 280 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*13 ) + 320 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*( hopi + NCn*14 ) + 20 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*14 ) + 160 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*14 ) + 240 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );
                H_0[ size*( hopi + NCn*14 ) + 260 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*14 ) + 300 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*14 ) + 340 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*( hopi + NCn*15 ) + 40 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*15 ) + 180 + hopi ] = - t_hub * pow( -1.0 , sing_s5[hopi] );
                H_0[ size*( hopi + NCn*15 ) + 280 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*15 ) + 360 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );

                H_0[ size*( hopi + NCn*16 ) + 260 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*16 ) + 340 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );

                H_0[ size*( hopi + NCn*17 ) + 80 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*17 ) + 280 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 320 + hopi ] = - t_hub * pow( -1.0 , sing_s1[hopi] );
                H_0[ size*( hopi + NCn*17 ) + 360 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );

                H_0[ size*( hopi + NCn*18 ) + 100 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*18 ) + 300 + hopi ] = - t_hub * pow( -1.0 , sing_s4[hopi] );
                H_0[ size*( hopi + NCn*18 ) + 340 + hopi ] = - t_hub * pow( -1.0 , sing_s2[hopi] );
                H_0[ size*( hopi + NCn*18 ) + 380 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

                H_0[ size*( hopi + NCn*19 ) + 140 + hopi ] = - t_hub * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
                H_0[ size*( hopi + NCn*19 ) + 360 + hopi ] = - t_hub * pow( -1.0 , sing_s3[hopi] );

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



            for (q = 0; q < 7337; q++) {

                if (t0<Tfull/2.0) {
                    Ax = -(Ex/wpump) * (cos(wpump*t0+pi*0.25)) * (pow(cos((pi/Tfull)*t0),4.0));
                }
                else{
                    Ax = 0.0;
                }

// current operator 


                gsl_matrix_complex_set_all(J_X, GSL_COMPLEX_ZERO);


                for ( hopi = 0; hopi < NCn; hopi++ ) {

                    gsl_matrix_complex_set(J_X, 0 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 0 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 0 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 1 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 2 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 3 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 4 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 5 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 3 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 6 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 7 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 8 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 9 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 10 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 10 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 11 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 0 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 6 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 12 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 11 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 13 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 1 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 8 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 12 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 14 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 2 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 9 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s6[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 15 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 13 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 16 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 4 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 14 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 16 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 17 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 5 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 15 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 17 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, 18 + NCn*hopi , 19 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 7 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, 19 + NCn*hopi , 18 + NCn*hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                }

// 

                for ( hopi = 0; hopi < NCn; hopi++ ) {

                    gsl_matrix_complex_set(J_X, NCn*0 + hopi , 20 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*0 + hopi , 240 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 0 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 40 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 80 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*1 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 20 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 60 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*2 + hopi , 300 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 40 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*3 + hopi , 120 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 20 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 200 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*4 + hopi , 340 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 40 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 80 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 120 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 140 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 220 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*5 + hopi , 360 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 60 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 160 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*6 + hopi , 240 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 160 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 260 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*7 + hopi , 380 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));

                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 120 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 140 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 180 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*8 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 160 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*9 + hopi , 300 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 80 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*10 + hopi , 220 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 200 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 240 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*11 + hopi , 260 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 0 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 120 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 220 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*12 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 140 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 220 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*13 + hopi , 320 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 20 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 160 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 240 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 260 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 300 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*14 + hopi , 340 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 40 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 180 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s5[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*15 + hopi , 360 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 260 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*16 + hopi , 340 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 80 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 280 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 320 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s1[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*17 + hopi , 360 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 100 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 300 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s4[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 340 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s2[hopi] ) ));
                    gsl_matrix_complex_set(J_X, NCn*18 + hopi , 380 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 140 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ),+latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) ) ));
                    gsl_matrix_complex_set(J_X, NCn*19 + hopi , 360 + hopi, gsl_complex_rect( -latt_const*t_hub*sin(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ),-latt_const*t_hub*cos(latt_const*Ax) * pow( -1.0 , sing_s3[hopi] ) ));

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

    double latt_const = 7.56;// lattice constant
    int rhi, rhj;

    double ellip = 0.0;
    double Ex = ( 3.0/6.5 )*0.00126410*alpha/(sqrt(1.0+ellip*ellip));// optical field strength
    double wpump = 0.26/27.2114; //pump frequency
    double Tfull = 4134.2*2.0; //pulse length
    double Ax;

    double u_rep = 2.0*0.52/27.2114;// Hubbard U
    double t_hub = 0.52/27.2114;// hopping t

    double gamma = 0.0 * 0.0316 * 0.52/27.2114;// electron-environment scattering rate


    arr = ( double _Complex * ) calloc ( m , sizeof ( double _Complex ) );

    sub_D = ( double * ) calloc ( m , sizeof ( double ) );

    if (t<Tfull/2.0) {
        Ax = -(Ex/wpump) * (cos(wpump*t+pi*0.25)) * (pow(cos((pi/Tfull)*t),4.0));
    }
    else{
        Ax = 0.0;
    }

    lapack_int n2 = 400;
    lapack_complex_double H_KA [n2*n2];//time-dep-hamil

// for hamiltonian

    int totn, dian, sn, NCn = 20, hopi;
    int size2 = 400;

// 

    double site1[n2], site2[n2], site3[n2], site4[n2], site5[n2], site6[n2];
    double sing_s1[NCn], sing_s2[NCn], sing_s3[NCn], sing_s4[NCn], sing_s5[NCn], sing_s6[NCn];


    sing_s1[0] = 1.0,  sing_s2[0] = 1.0,  sing_s3[0] = 1.0,  sing_s4[0] = 0.0,  sing_s5[0] = 0.0,  sing_s6[0] = 0.0;
    sing_s1[1] = 1.0,  sing_s2[1] = 1.0,  sing_s3[1] = 0.0,  sing_s4[1] = 1.0,  sing_s5[1] = 0.0,  sing_s6[1] = 0.0;
    sing_s1[2] = 1.0,  sing_s2[2] = 0.0,  sing_s3[2] = 1.0,  sing_s4[2] = 1.0,  sing_s5[2] = 0.0,  sing_s6[2] = 0.0;
    sing_s1[3] = 0.0,  sing_s2[3] = 1.0,  sing_s3[3] = 1.0,  sing_s4[3] = 1.0,  sing_s5[3] = 0.0,  sing_s6[3] = 0.0;
    sing_s1[4] = 1.0,  sing_s2[4] = 1.0,  sing_s3[4] = 0.0,  sing_s4[4] = 0.0,  sing_s5[4] = 1.0,  sing_s6[4] = 0.0;
    sing_s1[5] = 1.0,  sing_s2[5] = 0.0,  sing_s3[5] = 1.0,  sing_s4[5] = 0.0,  sing_s5[5] = 1.0,  sing_s6[5] = 0.0;
    sing_s1[6] = 0.0,  sing_s2[6] = 1.0,  sing_s3[6] = 1.0,  sing_s4[6] = 0.0,  sing_s5[6] = 1.0,  sing_s6[6] = 0.0;
    sing_s1[7] = 1.0,  sing_s2[7] = 0.0,  sing_s3[7] = 0.0,  sing_s4[7] = 1.0,  sing_s5[7] = 1.0,  sing_s6[7] = 0.0;
    sing_s1[8] = 0.0,  sing_s2[8] = 1.0,  sing_s3[8] = 0.0,  sing_s4[8] = 1.0,  sing_s5[8] = 1.0,  sing_s6[8] = 0.0;
    sing_s1[9] = 0.0,  sing_s2[9] = 0.0,  sing_s3[9] = 1.0,  sing_s4[9] = 1.0,  sing_s5[9] = 1.0,  sing_s6[9] = 0.0;
    sing_s1[10] = 1.0, sing_s2[10] = 1.0, sing_s3[10] = 0.0, sing_s4[10] = 0.0, sing_s5[10] = 0.0, sing_s6[10] = 1.0;
    sing_s1[11] = 1.0, sing_s2[11] = 0.0, sing_s3[11] = 1.0, sing_s4[11] = 0.0, sing_s5[11] = 0.0, sing_s6[11] = 1.0;
    sing_s1[12] = 0.0, sing_s2[12] = 1.0, sing_s3[12] = 1.0, sing_s4[12] = 0.0, sing_s5[12] = 0.0, sing_s6[12] = 1.0;
    sing_s1[13] = 1.0, sing_s2[13] = 0.0, sing_s3[13] = 0.0, sing_s4[13] = 1.0, sing_s5[13] = 0.0, sing_s6[13] = 1.0;
    sing_s1[14] = 0.0, sing_s2[14] = 1.0, sing_s3[14] = 0.0, sing_s4[14] = 1.0, sing_s5[14] = 0.0, sing_s6[14] = 1.0;
    sing_s1[15] = 0.0, sing_s2[15] = 0.0, sing_s3[15] = 1.0, sing_s4[15] = 1.0, sing_s5[15] = 0.0, sing_s6[15] = 1.0;
    sing_s1[16] = 1.0, sing_s2[16] = 0.0, sing_s3[16] = 0.0, sing_s4[16] = 0.0, sing_s5[16] = 1.0, sing_s6[16] = 1.0;
    sing_s1[17] = 0.0, sing_s2[17] = 1.0, sing_s3[17] = 0.0, sing_s4[17] = 0.0, sing_s5[17] = 1.0, sing_s6[17] = 1.0;
    sing_s1[18] = 0.0, sing_s2[18] = 0.0, sing_s3[18] = 1.0, sing_s4[18] = 0.0, sing_s5[18] = 1.0, sing_s6[18] = 1.0;
    sing_s1[19] = 0.0, sing_s2[19] = 0.0, sing_s3[19] = 0.0, sing_s4[19] = 1.0, sing_s5[19] = 1.0, sing_s6[19] = 1.0;



// 


    for ( totn = 0; totn < m; totn++ ) {
        H_KA[totn] = 0.0;
    }


    for ( dian = 0; dian < size2; dian++ ) {
        H_KA[ ( size2 + 1 )*( dian ) ] = u_rep;
    }


    for ( sn = 0; sn < NCn; sn++ ) {
        H_KA[ ( size2 + 1 )*( NCn + 1 )*( sn ) ] = 3.0 * u_rep;
    }


    for ( sn = 0; sn < NCn; sn++ ) {
        H_KA[ ( size2 + 1 )*( NCn - 1 )*( sn + 1 ) ] = 0.0;
    }
// 


    for ( sn = 0; sn < 1; sn++ ) {

        H_KA[ ( size2 + 1 )*( 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( 12 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*1 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*1 + 14 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*2 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*2 + 15 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*3 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*3 + 15 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*4 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*4 + 17 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*5 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*5 + 18 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*6 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*6 + 18 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*7 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*7 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*8 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*8 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*9 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 18 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*9 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*10 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*10 + 17 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*11 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*11 + 18 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*12 + 0 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*12 + 18 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*13 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*13 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*14 + 1 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*14 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*15 + 2 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 3 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 18 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*15 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*16 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 18 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*16 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*17 + 4 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 10 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 18 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*17 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*18 + 5 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 6 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 11 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 12 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*18 + 19 ) ] = 2.0 * u_rep;


        H_KA[ ( size2 + 1 )*( NCn*19 + 7 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 8 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 9 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 13 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 14 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 15 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 16 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 17 ) ] = 2.0 * u_rep;
        H_KA[ ( size2 + 1 )*( NCn*19 + 18 ) ] = 2.0 * u_rep;
    }

//

    for ( hopi = 0; hopi < NCn; hopi++ ) {

        H_KA[ size2*(NCn*hopi) + 1 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi) + 12 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+1) + 0 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+1) + 2 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+1) + 4 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+1) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+2) + 1 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+2) + 3 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+2) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+2) + 15 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+3) + 2 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+3) + 6 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*(NCn*hopi+4) + 1 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+4) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+4) + 10 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+4) + 17 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+5) + 2 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+5) + 4 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+5) + 6 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+5) + 7 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+5) + 11 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+5) + 18 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+6) + 3 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+6) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+6) + 8 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+6) + 12 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );

        H_KA[ size2*(NCn*hopi+7) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+7) + 8 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+7) + 13 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+7) + 19 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );

        H_KA[ size2*(NCn*hopi+8) + 6 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+8) + 7 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+8) + 9 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+8) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );

        H_KA[ size2*(NCn*hopi+9) + 8 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+9) + 15 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );

        H_KA[ size2*(NCn*hopi+10) + 4 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+10) + 11 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );

        H_KA[ size2*(NCn*hopi+11) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+11) + 10 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+11) + 12 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+11) + 13 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*(NCn*hopi+12) + 0 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+12) + 6 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+12) + 11 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+12) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*(NCn*hopi+13) + 7 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+13) + 11 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+13) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+13) + 16 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*(NCn*hopi+14) + 1 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+14) + 8 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+14) + 12 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*(NCn*hopi+14) + 13 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+14) + 15 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+14) + 17 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*(NCn*hopi+15) + 2 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+15) + 9 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s6[hopi] );
        H_KA[ size2*(NCn*hopi+15) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+15) + 18 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*(NCn*hopi+16) + 13 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+16) + 17 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );

        H_KA[ size2*(NCn*hopi+17) + 4 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+17) + 14 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+17) + 16 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*(NCn*hopi+17) + 18 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );

        H_KA[ size2*(NCn*hopi+18) + 5 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+18) + 15 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*(NCn*hopi+18) + 17 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*(NCn*hopi+18) + 19 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*(NCn*hopi+19) + 7 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] + sing_s6[hopi] ) );
        H_KA[ size2*(NCn*hopi+19) + 18 + NCn*hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

    }

// 

    for ( hopi = 0; hopi < NCn; hopi++ ) {

        H_KA[ size2*hopi + 20 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*hopi + 240 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[0] + sing_s3[0] + sing_s4[0] + sing_s5[0] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*1 ) + 0 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*1 ) + 40 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*1 ) + 80 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*1 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[1] + sing_s3[1] + sing_s4[1] + sing_s5[1] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*2 ) + 20 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*2 ) + 60 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*2 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*2 ) + 300 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[2] + sing_s3[2] + sing_s4[2] + sing_s5[2] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*3 ) + 40 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*3 ) + 120 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*( hopi + NCn*4 ) + 20 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*4 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*4 ) + 200 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*4 ) + 340 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[4] + sing_s3[4] + sing_s4[4] + sing_s5[4] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*5 ) + 40 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*5 ) + 80 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*5 ) + 120 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*5 ) + 140 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*5 ) + 220 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*5 ) + 360 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[5] + sing_s3[5] + sing_s4[5] + sing_s5[5] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*6 ) + 60 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*6 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*6 ) + 160 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*6 ) + 240 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*( hopi + NCn*7 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*7 ) + 160 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*7 ) + 260 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*7 ) + 380 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[7] + sing_s3[7] + sing_s4[7] + sing_s5[7] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );

        H_KA[ size2*( hopi + NCn*8 ) + 120 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*8 ) + 140 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*8 ) + 180 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*8 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*( hopi + NCn*9 ) + 160 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*9 ) + 300 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );

        H_KA[ size2*( hopi + NCn*10 ) + 80 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*10 ) + 220 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );

        H_KA[ size2*( hopi + NCn*11 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*11 ) + 200 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*11 ) + 240 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*11 ) + 260 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );

        H_KA[ size2*( hopi + NCn*12 ) + 0 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[12] + sing_s3[12] + sing_s4[12] + sing_s5[12] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*12 ) + 120 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*12 ) + 220 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*12 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );

        H_KA[ size2*( hopi + NCn*13 ) + 140 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*13 ) + 220 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*13 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*13 ) + 320 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*( hopi + NCn*14 ) + 20 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[14] + sing_s3[14] + sing_s4[14] + sing_s5[14] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*14 ) + 160 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*14 ) + 240 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );
        H_KA[ size2*( hopi + NCn*14 ) + 260 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*14 ) + 300 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*14 ) + 340 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*( hopi + NCn*15 ) + 40 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[15] + sing_s3[15] + sing_s4[15] + sing_s5[15] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*15 ) + 180 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s5[hopi] );
        H_KA[ size2*( hopi + NCn*15 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*15 ) + 360 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );

        H_KA[ size2*( hopi + NCn*16 ) + 260 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*16 ) + 340 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );

        H_KA[ size2*( hopi + NCn*17 ) + 80 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[17] + sing_s3[17] + sing_s4[17] + sing_s5[17] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*17 ) + 280 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*17 ) + 320 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s1[hopi] );
        H_KA[ size2*( hopi + NCn*17 ) + 360 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );

        H_KA[ size2*( hopi + NCn*18 ) + 100 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[18] + sing_s3[18] + sing_s4[18] + sing_s5[18] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*18 ) + 300 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s4[hopi] );
        H_KA[ size2*( hopi + NCn*18 ) + 340 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s2[hopi] );
        H_KA[ size2*( hopi + NCn*18 ) + 380 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );

        H_KA[ size2*( hopi + NCn*19 ) + 140 + hopi ] = - t_hub * ( cos( latt_const * Ax ) + 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , ( sing_s2[19] + sing_s3[19] + sing_s4[19] + sing_s5[19] + sing_s1[hopi] + sing_s2[hopi] + sing_s3[hopi] + sing_s4[hopi] + sing_s5[hopi] ) );
        H_KA[ size2*( hopi + NCn*19 ) + 360 + hopi ] = - t_hub * ( cos( latt_const * Ax ) - 1.0 * I * sin ( latt_const * Ax ) ) * pow( -1.0 , sing_s3[hopi] );


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

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, RHO_T, H_T,
                    GSL_COMPLEX_ZERO, COMM2);

// 

    for ( rhi = 0 ; rhi < NCn ; rhi++ ) {

        for ( rhj = 0 ; rhj < NCn ; rhj++ ) {

            site1[ rhi * NCn + rhj ] = sing_s1[ rhi ] + sing_s1[ rhj ];
            site2[ rhi * NCn + rhj ] = sing_s2[ rhi ] + sing_s2[ rhj ];
            site3[ rhi * NCn + rhj ] = sing_s3[ rhi ] + sing_s3[ rhj ];
            site4[ rhi * NCn + rhj ] = sing_s4[ rhi ] + sing_s4[ rhj ];
            site5[ rhi * NCn + rhj ] = sing_s5[ rhi ] + sing_s5[ rhj ];
            site6[ rhi * NCn + rhj ] = sing_s6[ rhi ] + sing_s6[ rhj ];

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

        }
    }

// 


    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {

            gsl_matrix_complex_set(D, rhi, rhj, gsl_complex_rect( sub_D[ rhi * size2 + rhj ] * GSL_REAL(gsl_matrix_complex_get(RHO_T,rhi,rhj)) , sub_D[ rhi * size2 + rhj ] * GSL_IMAG(gsl_matrix_complex_get(RHO_T,rhi,rhj)) ));

        }
    }


// 


    for ( rhi = 0 ; rhi < size2 ; rhi++ ) {
        for ( rhj = 0 ; rhj < size2 ; rhj++ ) {
            arr[rhi * size2 + rhj] = -1.0*I * ( GSL_REAL(gsl_matrix_complex_get(COMM1,rhi,rhj)) + GSL_IMAG(gsl_matrix_complex_get(COMM1,rhi,rhj))*I - GSL_REAL(gsl_matrix_complex_get(COMM2,rhi,rhj)) - GSL_IMAG(gsl_matrix_complex_get(COMM2,rhi,rhj))*I )
                                                + gamma * ( GSL_REAL(gsl_matrix_complex_get(D,rhi,rhj)) + GSL_IMAG(gsl_matrix_complex_get(D,rhi,rhj))*I );
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
