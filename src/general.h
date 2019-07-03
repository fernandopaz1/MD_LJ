#ifndef GENERAL_H
#define GENERAL_H

#include "math.h"

// Declaracion de funciones
double aleatorio(void);
double gaussiana(double mu, double sigma);
double norma(double *v, int longitud);
void histograma(double *y,double *x,int n,double a,double b,int m);
int minmax(double *vector,int N, double *min, double *max);
int pruebageneral(int N, double mu, double sigma,double radio);
double energia_cinetica(double *velocidad,double *E_cinetica, int N);
double energia_potencial(double *potencial,double *E_potencial, int N, double masa);
double energia(double *velocidad, double *potencial,double *E_cinetica,double *E_potencial, int N, double masa);
int normalizacion_velocidad(double *velocidad, double *E_cinetica,double T, int N);
#endif
