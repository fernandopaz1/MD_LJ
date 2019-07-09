#ifndef AVANZAR_H
#define AVANZAR_H

#include "math.h"

// Declaracion de funciones
int velocity_verlet(double *x, double *velocidad, double *fuerza, double *fuerza_vieja, double h, int N, double L);
int paso(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial, int N, double rc, double V0, double L, double h);
int descorrelacionar(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial,int N, double rc, double V0, double L, double h);
int simulacion(int N, int N_frames,double rc, double L, double h, double T);
double coeficiente_verlet(double L, int N, double *x);
#endif
