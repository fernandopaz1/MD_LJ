#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
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


int paso(double *x, double *velocidad,double *fuerza, double *fuerza_vieja, double *potencial,int N, double rc, double V0, double L, double h){
	int i, j;

	

	for(i=0;i<3*N;i++){	
		*(fuerza_vieja+i)=*(fuerza+i);
		*(fuerza+i)=0.0;
	}																											

	
	for(i=0;i<3*N;i++){
		for(j=i+1;j<3*N;j++){
			Lenard_Jones(fuerza, potencial, 1.0, 1.0, rc, x, i, j, V0,L);
		}	
	}


	velocity_verlet(x,velocidad,fuerza,fuerza_vieja,h,N);


return 0;
}

int simulacion(int N, int N_frames,double rc, double L, double h, double T){
	int l;
	double *x, *velocidad, *fuerza, *fuerza_vieja, *potencial,V0;
	

	x=(double*)malloc(3*N*sizeof(double));
	velocidad=(double*)malloc(3*N*sizeof(double));
	fuerza=(double*)malloc(3*N*sizeof(double));
	fuerza_vieja=(double*)malloc(3*N*sizeof(double));
	potencial=(double*)malloc(N*sizeof(double));
	

	
	char filename[255];
  	sprintf(filename, "prueba.lammpstrj");

	set_box(x, N, L);
	set_v(velocidad,N,T);


	V0=V0_LJ(rc,1.0, 1.0);

	for(l=0;l<N_frames;l++){
		

		paso(x, velocidad, fuerza, fuerza_vieja,potencial, N, rc, V0, L, h);

		save_lammpstrj(filename, x, velocidad, N, L, l);	

	}
	

free(x);
free(velocidad);
free(fuerza);
free(fuerza_vieja);
free(potencial);
return 0;
}


