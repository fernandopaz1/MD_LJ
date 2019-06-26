#include "general.h"
#include "interaccion.h"
#include "avanzar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Definicion de funciones


int velocity_verlet(double *x, double *velocidad, double *fuerza, double *fuerza_vieja, double h, int N){
	int i;

	for(i=0; i<3*N; i++){
		*(x+i)+=*(velocidad+i)*h+0.5*h*h*(*(fuerza+i));
		*(velocidad+i)+=0.5*h*((*(fuerza+i))+(*(fuerza_vieja+i)));
	}
return 0;
}


