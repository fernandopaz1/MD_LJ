#ifndef INTERACCION_H
#define INTERACCION_H

#include "math.h"	

// Declaracion de funciones
double delta(double *x, int i, int j, double L,double *delta_x);
double Lenard_Jones(double *fuerzas, double *potencial, double e, double sigma, double rc, double *x, int i, int j, double V0, double L,double *p);
double V0_LJ(double rc,double sigma, double epsilon);
int tabla(double epsilon, double sigma, double paso);
int Lenard_Jones_interp(double *fuerzas, double *potencial, double e, double sigma, double rc, double *x, int i, int j, double V0, double *finterp, double paso, double L);

#endif
