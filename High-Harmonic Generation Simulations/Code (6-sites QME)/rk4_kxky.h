#ifndef UNTITLED_RK4_KXKY_H
#define UNTITLED_RK4_KXKY_H

double _Complex *rk4veckxky ( double t0, int n, double _Complex u0[], double dt, double tau0, double kx, double ky, double alpha,
                          double _Complex *f ( double t, int n, double _Complex u[], double tau, double kx, double ky, double alpha ) );

#endif //UNTITLED_RK4_KXKY_H
