#include <stdlib.h>

/******************************************************************************/

double _Complex *rk4veckxky ( double t0, int m, double _Complex u0[], double dt, double tau0, double kx, double ky, double alpha,
                          double _Complex *f ( double t, int m, double _Complex u[], double tau, double kx, double ky, double alpha ) )

/******************************************************************************/

{
    double _Complex *f0;
    double _Complex *f1;
    double _Complex *f2;
    double _Complex *f3;
    int i;
    double t1;
    double t2;
    double t3;
    double _Complex *u;
    double _Complex *u1;
    double _Complex *u2;
    double _Complex *u3;

    f0 = f ( t0, m, u0, tau0, kx, ky, alpha );

    t1 = t0 + dt / 2.0;
    u1 = ( double _Complex * ) malloc ( m * sizeof ( double _Complex ) );
    for ( i = 0; i < m; i++ )
    {
        u1[i] = u0[i] + dt * f0[i] / 2.0;
    }
    f1 = f ( t1, m, u1, tau0, kx, ky, alpha );

    t2 = t0 + dt / 2.0;
    u2 = ( double _Complex * ) malloc ( m * sizeof ( double _Complex ) );
    for ( i = 0; i < m; i++ )
    {
        u2[i] = u0[i] + dt * f1[i] / 2.0;
    }
    f2 = f ( t2, m, u2, tau0, kx, ky, alpha );

    t3 = t0 + dt;
    u3 = ( double _Complex * ) malloc ( m * sizeof ( double _Complex ) );
    for ( i = 0; i < m; i++ )
    {
        u3[i] = u0[i] + dt * f2[i];
    }
    f3 = f ( t3, m, u3, tau0, kx, ky, alpha );

    u = ( double _Complex * ) malloc ( m * sizeof ( double _Complex ) );
    for ( i = 0; i < m; i++ )
    {
        u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
    }

    free ( f0 );
    free ( f1 );
    free ( f2 );
    free ( f3 );
    free ( u1 );
    free ( u2 );
    free ( u3 );

    return u;
}

