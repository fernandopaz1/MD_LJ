#ifndef INICIALIZAR_H
#define INICIALIZAR_H

#include "math.h"

// Declaracion de funciones
int set_box(double *x,int N,double L);
int set_v(double *velocidad,int N,double T);
int save_checkpoint(char *nombre, double *x, double *velocidad,int N);
int load_checkpoint(char *nombre, double *x, double *velocidad,int N);
#endif
